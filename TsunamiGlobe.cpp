// Copyright (c) 2015-2017 John M. Wilson, Kasey W. Schultz
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include "TsunamiGlobe.h"
#include <cassert>

#define assertThrow(COND, ERR_MSG) assert(COND);

// ----------------------------------------------------------------------
// -------------------- Main Functions --------------------------------
// ----------------------------------------------------------------------

// Set the height for all elements equal to the depth of the bathymetry below the center of the square.
// Result is squares with just enough water so that the water sits at sea level.
void tsunamisquares::World::fillToSeaLevel(void) {
    std::map<UIndex, Square>::iterator     sit;

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        // Add water if the  altitude of the Square center is below sea level
        if (sit->second.xyz()[2] < 0.0) {
            sit->second.set_height(fabs(sit->second.xyz()[2]));
        } else {
            sit->second.set_height(0.0);
        }
        // Also initialize the velocities and accel to zero
        sit->second.set_velocity(Vec<2>(0.0,0.0));
        sit->second.set_accel(Vec<2>(0.0,0.0));
    }
}


// Diffusion: Remove a volume of water from each square and distribute it to the neighbors.
// Model: area_change = diff_const*time_step
void tsunamisquares::World::diffuseSquares(const double dt, const double D) {
    std::map<UIndex, Square>::iterator  it;
    SquareIDSet                         neighborIDs;
    std::map<UIndex, Square>::iterator  nit;
    double                              volume_change, new_level, add_height, height_change;
    Vec<2>                              momentum_change;
    SquareIDSet::iterator               id_it;
    bool debug = false;
    
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //TODO: Check that I do not diffuse into dry squares (wetting them artificially)
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // Initialize updated_heights and momenta, will use this to store the net height and momentum changes
    for (it=_squares.begin(); it!=_squares.end(); ++it) {
        it->second.set_updated_height( it->second.height() );
        it->second.set_updated_momentum( it->second.momentum() );
    }

    // Compute the height changes due to diffusion of water to neighbors
    for (it=_squares.begin(); it!=_squares.end(); ++it) {
        if (it->second.height() > 0) { //}) && squareLevel(it->first) != 0.0) {
            // Compute the new height after diffusing the water by 1 time step
            new_level = squareLevel(it->first)/(1 + D*dt/it->second.area());
            volume_change = (it->second.area())*(squareLevel(it->first) - new_level);
            //assertThrow(volume_change >= 0, "Volume change should be positive");
            height_change = new_level - squareLevel(it->first);
            // Transfer the proportional amount of momentum
            momentum_change = (it->second.momentum())*volume_change/(it->second.volume());

            if (debug) {
                std::cout << "----> Diffusing Square " << it->second.id() << std::endl;
                std::cout << "volume change: " << volume_change << std::endl;
                std::cout << "old level: " << squareLevel(it->first) << std::endl;
                std::cout << "new level: " << new_level << std::endl;
                std::cout << "-> neighbors " << std::endl;
            }
            
            // For continuity, must self-add 1/4 of the volume change to edges and 1/2 to corners.
            // This also balances the momentum distribution.
            int minLat = squareLatLon(it->first)[0] == min_lat();
            int maxLat = squareLatLon(it->first)[0] == max_lat();
            int minLon = squareLatLon(it->first)[1] == min_lon();
            int maxLon = squareLatLon(it->first)[1] == max_lon();
            int cornerSum = minLat + minLon + maxLon + maxLat;    
            if (cornerSum == 1) {
                height_change += volume_change/( it->second.area()*4.0);
            } else if (cornerSum == 2) {
                height_change += volume_change/( it->second.area()*2.0);
            }
            
            // Add the height change to the updated height
            it->second.set_updated_height(it->second.updated_height() + height_change);

            neighborIDs = it->second.get_valid_nearest_neighbors();
            
            for (id_it=neighborIDs.begin(); id_it!=neighborIDs.end(); ++id_it) {
                nit = _squares.find(*id_it);
                // Divide up the diffused volume equally amongst neighbors
                add_height = volume_change/( nit->second.area()*4.0);
                nit->second.set_updated_height( nit->second.updated_height() + add_height);
                nit->second.set_updated_momentum( nit->second.updated_momentum() + momentum_change/4.0);
            }
        }
    }
    
    // Reset the heights and velocities based on the changes
    for (it=_squares.begin(); it!=_squares.end(); ++it) {
        it->second.set_height( it->second.updated_height() );
        it->second.set_velocity( it->second.updated_momentum() / it->second.mass());
    }
}


// Move the water from a Square given its current velocity and acceleration.
// Partition the volume and momentum into the neighboring Squares.
void tsunamisquares::World::moveSquares(const double dt, const bool accel_bool) {
    std::map<UIndex, Square>::iterator sit;
    Geodesic geod(EARTH_MEAN_RADIUS, 0);         /* <(^_^<) Happy Coder says: We're good up through this line!*/

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
		if(isnan(sit->second.height())) std::cout << "ID " << sit->first << ": Height NaN"   << std::endl;
		if(isnan(sit->second.velocity()[0])) std::cout << "ID " << sit->first << ": Velocity NaN" << std::endl;
		if(isnan(sit->second.momentum()[0])) std::cout << "ID " << sit->first << ": Momentum NaN" << std::endl;
		if(isnan(sit->second.accel()[0]))    std::cout << "ID " << sit->first << ": Accel NaN"    << std::endl;
		if(isnan(sit->second.xyz()[0]))      std::cout << "ID " << sit->first << ": Position NaN" << std::endl;
    }

    // Initialize the updated height and velocity to zero. These are the containers
    // used to keep track of the distributed height/velocity from moving squares.
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        sit->second.set_updated_height(0.0);
        Vec<2> m; m[0] = m[1] = 0.0;
        sit->second.set_updated_momentum(m);
        // Set acceleration based on the current slope of the water surface
        updateAcceleration(sit->first);
    }

    // Now go through each square and move the water, distribute to neighbors
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        Vec<2> current_velo, current_accel, current_pos, new_velo, average_velo;
        double local_azimuth, distance_traveled, av_velo_mag, lat2, lon2, a12;
        std::map<double, UIndex> distsNneighbors;
		std::map<double, UIndex>::const_iterator dnit;
        point_spheq new_bottom_left, new_bottom_right, new_top_left, new_top_right;

        current_pos = squareCenter(sit->first);
        current_velo = sit->second.velocity();
        if(accel_bool){
        	current_accel = sit->second.accel();
        }else{
        	current_accel = Vec<2>(0,0);
        }
        

        // Move the square: calculate average velocity during motion, find azimuth of that vector, and move
        //  each vertex of the square to it's new location on the sphere.  This forms a Ring, which intersects some number of
        //  boxes in our rtree.
        new_velo = current_velo + current_accel*dt;
        average_velo = current_velo + current_accel*0.5*dt;
        distance_traveled = average_velo.mag()*dt;

        //If the square didn't move, or has run out of water, immediately distribute it's water back to itself
        //   (along with any contributions from other squares)
        if(average_velo == Vec<2>(0.0, 0.0) || sit->second.height() == 0.0){
        	sit->second.set_updated_height(sit->second.updated_height() + sit->second.height());
        	sit->second.set_updated_momentum(sit->second.updated_momentum() + sit->second.momentum());
        }else{
        	// Calculate azimuth for geodesic calculation
			if(atan2(average_velo[1], average_velo[0]) >= -(M_PI/2)){
				local_azimuth = 90-atan2(average_velo[1], average_velo[0])*(180/M_PI);
			}else{
				local_azimuth = -270-atan2(average_velo[1], average_velo[0])*(180/M_PI);
			}

			//bottom left
			geod.Direct(current_pos[1]-_dlat/2, current_pos[0]-_dlon/2, local_azimuth, distance_traveled, lat2, lon2);
			new_bottom_left = point_spheq(lon2, lat2);

			//bottom right
			geod.Direct(current_pos[1]-_dlat/2, current_pos[0]+_dlon/2, local_azimuth, distance_traveled, lat2, lon2);
			new_bottom_right = point_spheq(lon2, lat2);

			//top left
			geod.Direct(current_pos[1]+_dlat/2, current_pos[0]-_dlon/2, local_azimuth, distance_traveled, lat2, lon2);
			new_top_left = point_spheq(lon2, lat2);

			//top right
			geod.Direct(current_pos[1]+_dlat/2, current_pos[0]+_dlon/2, local_azimuth, distance_traveled, lat2, lon2);
			new_top_right = point_spheq(lon2, lat2);

			// Make Ring from new vertices for accurate overlap calc later
			point_spheq ring_verts[5] = {new_bottom_left, new_top_left, new_top_right, new_bottom_right, new_bottom_left};
			ring_spheq new_ring;
			bg::assign_points(new_ring, ring_verts);

			// Do a quick grab of nearest squares to new location, then do the accurate intersection check
			distsNneighbors = getNearest_rtree(Vec<2>(new_bottom_left.get<0>()+_dlon/2, new_bottom_left.get<1>()+_dlat/2), 5);


			// Init these for renormalizing the fractions
			double this_fraction;
			double fraction_sum = 0.0;
			std::map<UIndex, double> originalFractions, renormFractions;
			std::map<UIndex, double>::iterator frac_it;


			// Iterate through the neighbors and compute overlap area between ring and neighbor boxes
			int is_any_overlap = 0;
			for (dnit=distsNneighbors.begin(); dnit!=distsNneighbors.end(); ++dnit) {
				std::map<UIndex, Square>::iterator neighbor_it = _squares.find(dnit->second);
				std::vector<poly_spheq> output;
				double overlap_area=0.0;

				bg::intersection(new_ring, neighbor_it->second.box(), output);

				if(output.size() != 0){
					BOOST_FOREACH(poly_spheq const& p, output)
					{
						overlap_area = bg::area(p)*EARTH_MEAN_RADIUS*EARTH_MEAN_RADIUS;
					}

					this_fraction = overlap_area/neighbor_it->second.area();
					fraction_sum += this_fraction;
					originalFractions.insert(std::make_pair(neighbor_it->first, this_fraction));
					is_any_overlap++;
				}
			}

			// If no overlap, we want to distribute water to nearest square
			// TODO: We should implement a line (geodesicLine? bg::Linestring?) from origin to destination,
			//       and accurately stop (reflect) the square when it hits invalid areas.
			if(!is_any_overlap){
				originalFractions.insert(std::make_pair(distsNneighbors.begin()->second, 1.0));
				fraction_sum = 1.0;
			}

			// Then normalize these fractions to enforce conservation.

			for (frac_it=originalFractions.begin(); frac_it!=originalFractions.end(); ++frac_it) {
				//assertThrow((frac_it->second)/fraction_sum <= 1, "Area fraction must be less than 1.");
				renormFractions.insert(std::make_pair(frac_it->first, (frac_it->second)/fraction_sum));
			}

			// Compute height and momentum imparted to neighbors
			for (frac_it=renormFractions.begin(); frac_it!=renormFractions.end(); ++frac_it) {
				// This iterator will give us the neighbor square
				std::map<UIndex, Square>::iterator neighbor_it = _squares.find(frac_it->first);
				//// This iterates through the renormalized fractions
				//frac_it = renormFractions.find(nit->second);
				double areaFraction = frac_it->second;

				// Update the amount of water in the neighboring square (conserve volume)
				double dV = sit->second.volume()*areaFraction;
				double H = neighbor_it->second.updated_height();
				double A_n = neighbor_it->second.area();
				neighbor_it->second.set_updated_height(H + dV/A_n);

				// Conserve momentum, update the velocity accordingly (at the end)
				Vec<2> dM = new_velo*sit->second.mass()*areaFraction;
				Vec<2> M  = neighbor_it->second.updated_momentum();
				neighbor_it->second.set_updated_momentum(M+dM);


			}
        }
    }

    // Loop again over squares to set new velocity and height from accumulated height and momentum
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
    	Vec<2> velo_to_set;

        sit->second.set_height(sit->second.updated_height());

        velo_to_set = sit->second.updated_momentum()/(sit->second.mass());

        // Check for invalid directions of motion (eg at edges or near land)
        // TODO: add parameter to toggle btwn elastic and inelastic collisions with edges (currently inelastic).
        //if(sit->second.invalid_directions()[0] && velo_to_set[0]<0.0) velo_to_set[0] = 0.0;
        //if(sit->second.invalid_directions()[1] && velo_to_set[0]>0.0) velo_to_set[0] = 0.0;
        //if(sit->second.invalid_directions()[2] && velo_to_set[1]<0.0) velo_to_set[1] = 0.0;
        //if(sit->second.invalid_directions()[3] && velo_to_set[1]>0.0) velo_to_set[1] = 0.0;

        sit->second.set_velocity(velo_to_set);


    }
    
}

void tsunamisquares::World::updateAcceleration(const UIndex &square_id) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    Vec<2> grav_accel, friction_accel, new_accel, gradient;
    double G = 9.80665; //mean gravitational acceleration at Earth's surface [NIST]
    
    // Only accelerate the water in this square IF there is water in this square
    if (square_it->second.height() != 0.0) {
        // gravitational acceleration due to the slope of the water surface
        gradient = getGradient_planeFit(square_id);
        
        grav_accel = gradient*G*(-1.0);

        // frictional acceleration from fluid particle interaction
        friction_accel = square_it->second.velocity()*(square_it->second.velocity().mag())*(square_it->second.friction())/(-1.0*(square_it->second.height()));
        
        new_accel = grav_accel + friction_accel;

        // Check for invalid directions of motion (eg at edges or near land)
        if(square_it->second.invalid_directions()[0] && new_accel[0]<0.0) new_accel[0] = 0.0;
        if(square_it->second.invalid_directions()[1] && new_accel[0]>0.0) new_accel[0] = 0.0;
        if(square_it->second.invalid_directions()[2] && new_accel[1]<0.0) new_accel[1] = 0.0;
        if(square_it->second.invalid_directions()[3] && new_accel[1]>0.0) new_accel[1] = 0.0;

        // Set the acceleration
        square_it->second.set_accel(new_accel);
    } else {
        square_it->second.set_accel( Vec<2>(0.0, 0.0) );
    }
}


tsunamisquares::SquareIDSet tsunamisquares::World::get_neighbors_for_accel(const UIndex &square_id) const {
	double															  thisLevel;
    SquareIDSet                 	      valid_squares, all_neighbors_and_self;
    SquareIDSet::iterator       	                                      id_it;

    thisLevel = squareLevel(square_id);
    // Grab all valid neighbors
    all_neighbors_and_self = _squares.find(square_id)->second.get_neighbors_and_self();
    
    // Only include the next nearest neighbors if they are not "hi and dry".
    // A wave incident on the beach is not pushed backwards by the tall beach in front of it.
    // The wave only falls back into the ocean after it has creeped up the beach and has water
    // above and below it that define a slope for the water surface.
    
    for (id_it=all_neighbors_and_self.begin(); id_it!=all_neighbors_and_self.end(); ++id_it) {
        /*if (!( (squareLevel(*id_it) == 0) && (squareDepth(*id_it) >= 0))) {
            valid_squares.insert(*id_it);
        }*/
    	// TODO: Do we want to use lower dry cells in acceleration calc?
        if (_squares.find(*id_it)->second.height()>0.0 || squareDepth(*id_it) < thisLevel) {
			valid_squares.insert(*id_it);
		}
    }
    
    return valid_squares;
}

tsunamisquares::Vec<2> tsunamisquares::World::fitPointsToPlane(const UIndex &this_id, const SquareIDSet &square_ids) {
    // --------------------------------------------------------------------
    // Based on StackOverflow article:
    // http://stackoverflow.com/questions/1400213/3d-least-squares-plane
    // --------------------------------------------------------------------
    std::vector<double>             x_vals, y_vals, z_vals;
    SquareIDSet::const_iterator                      id_it;
    Vec<2>                                        gradient;
    SquareIDSet::const_iterator iit;
    std::map<UIndex, Vec<2> > neighborsAndCoords, neighbors_for_fitting;
    int                             i, N = square_ids.size();
    Vec<9> A;
    Vec<3> b, x;
    
	neighborsAndCoords = square(this_id).local_neighbor_coords();

    for (id_it=square_ids.begin(); id_it!=square_ids.end(); ++id_it) {
        x_vals.push_back(neighborsAndCoords[*id_it][0]);
        y_vals.push_back(neighborsAndCoords[*id_it][1]);
        z_vals.push_back(squareLevel(*id_it));
    }

    // Build the b vector and the A matrix.
    // Single index for matrix, array style. A[i][j] = A_vec[i*3 + j],  N_cols=3
    for (i=0; i<N; ++i) {
            
        b[0] += x_vals[i]*z_vals[i];
        b[1] += y_vals[i]*z_vals[i];
        b[2] += z_vals[i];
        
        A[0] += x_vals[i]*x_vals[i];
        A[1] += x_vals[i]*y_vals[i];
        A[2] += x_vals[i];
        A[3] += x_vals[i]*y_vals[i];
        A[4] += y_vals[i]*y_vals[i];
        A[5] += y_vals[i];
        A[6] += x_vals[i];
        A[7] += y_vals[i];
        A[8] += 1.0;
    }

    //Cramer's rule for 3x3 system
    float det_A = (A[0]*A[4]*A[8]+A[1]*A[5]*A[6]+A[3]*A[7]*A[2]-A[2]*A[4]*A[6]-A[1]*A[3]*A[8]-A[0]*A[5]*A[7]);

    x[0] = (b[0]*A[4]*A[8]+A[1]*A[5]*b[2]+b[1]*A[7]*A[2]-A[2]*A[4]*b[2]-A[1]*b[1]*A[8]-b[0]*A[5]*A[7])/det_A;
    x[1] = (A[0]*b[1]*A[8]+b[0]*A[5]*A[6]+A[3]*b[2]*A[2]-A[2]*b[1]*A[6]-b[0]*A[3]*A[8]-A[0]*A[5]*b[2])/det_A;

    gradient[0] = x[0];
    gradient[1] = x[1];

    return gradient;

    // Matrix determinant debugging
    /*if(det_A == 0.0){
     	std::cout<<" 0 determinant! " << std::endl;
        //std::cout << "\npre x: " << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
        //std::cout << "\n\nb: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
        std::cout << std::fixed << "A: " << A[0] << ",\t " << A[1] << ",\t " << A[2] << std::endl;
        std::cout << std::fixed << "   " << A[3] << ",\t " << A[4] << ",\t " << A[5] << std::endl;
        std::cout << std::fixed << "   " << A[6] << ",\t " << A[7] << ",\t " << A[8] << std::endl;
        std::cout << std::endl;
        for(std::vector<double>::iterator i=x_vals.begin(); i!=x_vals.end(); i++){
        	std::cout << std::fixed << "\txs: " << *i << std::endl;
        }
        for(std::vector<double>::iterator i=y_vals.begin(); i!=y_vals.end(); i++){
        	std::cout << std::fixed << "\tys: " << *i << std::endl;
        }
    }*/
    
    // Matrix solver below is adapted from Virtual Quake
    /*
    int     j, k;
    double  v, f, sum;
    int     n = 3;
    for (i=0; i<n; ++i) {
        v = A[i+n*i];
        for (j=i+1; j<n; ++j) {
            f = A[i+n*j]/v;
            for (k=0; k<n; ++k) {
                A[k+n*j] -= f*A[k+n*i];
            }
            b[j] -= f*b[i];
        }
    }
    for (i=n-1; i>=0; --i) {
        sum = b[i];
        for (j=i+1; j<n; ++j) {
            sum -= A[j+n*i]*x[j];
        }
        x[i] = sum/A[i+n*i];
    }
    */
}

tsunamisquares::Vec<2> tsunamisquares::World::getGradient_planeFit(const UIndex &square_id) {
    std::map<UIndex, Square>::const_iterator square_it = _squares.find(square_id);
    Vec<2> gradient;
    bool debug = false;
    SquareIDSet square_ids_to_fit;
    
	square_ids_to_fit = get_neighbors_for_accel(square_id);

	gradient = fitPointsToPlane(square_id, square_ids_to_fit);

    return gradient;
}


tsunamisquares::Vec<2> tsunamisquares::World::getGradient(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator square_it = _squares.find(square_id);
    Vec<2> gradient;
    bool debug = false;
    
    // Initialize the 4 points that will be used to approximate the slopes d/dx and d/dy
    // for this square. These are the centers of the neighbor squares.
    Vec<2> center = squareCenter(square_id);   
    
    UIndex leftID   = square_it->second.left();
    UIndex rightID  = square_it->second.right();
    UIndex topID    = square_it->second.top();
    UIndex bottomID = square_it->second.bottom();
    
    // TODO: Better boundary conditions. For now, just set no acceleration along boundary.
    // ALSO: Not accelerating squares whose neighbor is on a boundary.
    if (square_it->second.xy()[1] == min_lat() || square_it->second.xy()[1] == max_lat() || square_it->second.xy()[0] == min_lon() || square_it->second.xy()[0] == max_lon() || leftID == INVALID_INDEX || rightID == INVALID_INDEX || topID == INVALID_INDEX || bottomID == INVALID_INDEX) {
        gradient = Vec<2>(0.0,0.0);
    } else {
        // Altitude of water level of neighbor squares
        double z_left   = squareLevel(leftID);
        double z_right  = squareLevel(rightID);
        double z_top    = squareLevel(topID);
        double z_bottom = squareLevel(bottomID);
        double z_mid    = squareLevel(square_id);

        // Thickness of water in neighbor squares
        double h_left   = _squares.find(leftID)->second.height();
        double h_right  = _squares.find(rightID)->second.height();
        double h_top    = _squares.find(topID  )->second.height();
        double h_bottom = _squares.find(bottomID)->second.height();
        double h_mid    = square_it->second.height();
        
        // X,Y of neighbor squares
        Vec<2> center_L = squareCenter(leftID);
        Vec<2> center_R = squareCenter(rightID);
        Vec<2> center_T = squareCenter(topID);
        Vec<2> center_B = squareCenter(bottomID);
        
        // ================================================================
        // Gradient = (dz/dx, dz/dy)
        // Handle the cases with dry cells on either left/right/top/bottom.
        // IGNORE cells that are hi and dry
        // ================================================================
        if (h_left == 0.0 && h_right == 0.0 && h_top == 0.0 && h_bottom == 0.0) {
        // Case: No water on any side
        // TODO: Is this check needed?
            gradient[0] = 0.0;
            gradient[1] = 0.0;
        } else  {
            if (h_left > 0.0 && h_right > 0.0 && h_top > 0.0 && h_bottom > 0.0) {
            // Case: No dry neighbors, then do normal gradient
                gradient[0] = (z_right-z_left)/( center_L.dist(center_R) );
                gradient[1] = (z_top-z_bottom)/( center_T.dist(center_B) );
            }
            
            // Case: Hi and dry on the right, water to the left
            if (h_right == 0.0 && z_right >= 0.0 && h_left != 0.0) {
                gradient[0] = (z_mid-z_left)/( center_L.dist(center) );
            } else if (h_left == 0.0 && z_left >= 0.0 && h_right != 0.0) {
            // Case: Hi and dry on the left, water to the right
                gradient[0] = (z_right-z_mid)/( center_R.dist(center) );
            }

            
            // Case: Hi and dry on the top, water on bottom
            if (h_top == 0.0 && z_top >= 0.0 && h_bottom != 0.0) {
                gradient[1] = (z_mid-z_bottom)/( center.dist(center_B) );
            } else if (h_left == 0.0 && z_left >= 0.0 && h_right != 0.0) {
            // Case: Hi and dry on the bottom, water on top
                gradient[1] = (z_top-z_mid)/( center_T.dist(center) );
            }
            
        }
        
        if (debug) {
            std::cout << "square  " << square_id << std::endl;
            std::cout << "d/dx " << gradient[0] << std::endl; 
            std::cout << "d/dy " << gradient[1] << std::endl;
        }
    }
    
    return gradient;
    
}


// Raise/lower the sea floor depth at the square's vertex by an amount "height_change"
void tsunamisquares::World::deformBottom(const UIndex &square_id, const double &height_change) {
    std::map<UIndex, Square>::iterator sit = _squares.find(square_id);
    LatLonDepth new_lld;
    double old_altitude;
    
    new_lld = sit->second.lld();
    old_altitude = new_lld.altitude();
    new_lld.set_altitude(old_altitude + height_change);
    sit->second.set_lld(new_lld);
}


void tsunamisquares::World::bumpCenter(const double bump_height) {
	LatLonDepth min_bound, max_bound, centerLLD;
	Vec<3> min_xyz, max_xyz;
	double dx, dy;
	get_bounds(min_bound, max_bound);

	Conversion c(min_bound);
    min_xyz = Vec<3>(min_bound.lon(), min_bound.lat(), min_bound.altitude());
    max_xyz = Vec<3>(max_bound.lon(), max_bound.lat(), max_bound.altitude());
    dx = max_xyz[0]-min_xyz[0];
    dy = max_xyz[1]-min_xyz[1];

	Vec<2> centerLoc = Vec<2>(min_xyz[0]+dx/2.0, min_xyz[1]+dy/2.0);

	UIndex centralID = getNearest_rtree(centerLoc, 1).begin()->second;

	for (int i = 0; i < 2; i ++ ){
		deformBottom(centralID + num_lons()*i,   bump_height/2.0);
		deformBottom(centralID + num_lons()*i+1, bump_height/2.0);
	}
	for (int i = 0; i < 4; i ++ ){
		//if(i==1 || i==2){deformBottom(centralID + num_lons()*(i-1)-1 , bump_height/2.0);}
		deformBottom(centralID + num_lons()*(i-1)-1 , bump_height/2.0);
		deformBottom(centralID + num_lons()*(i-1)   , bump_height/2.0);
		deformBottom(centralID + num_lons()*(i-1)+1 , bump_height/2.0);
		deformBottom(centralID + num_lons()*(i-1)+2 , bump_height/2.0);
		//if(i==1 || i==2){deformBottom(centralID + num_lons()*(i-1)+2 , bump_height/2.0);}
	}
}


// Flatten the bottom to be the specified depth
void tsunamisquares::World::flattenBottom(const double &depth) {
    std::map<UIndex, Square>::iterator sit;
    LatLonDepth new_lld;
    double newDepth = -fabs(depth);
    
    // Assign the depth for all vertices to be newDepth
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        new_lld = sit->second.lld();
        new_lld.set_altitude(newDepth);
        sit->second.set_lld(new_lld);
    }
}
//
//// ----------------------------------------------------------------------
//// -------------------- Utility Functions -------------------------------
//// ----------------------------------------------------------------------

// Find maximum depth in simulation
double tsunamisquares::World::getMaxDepth() const {
    std::map<UIndex, Square>::const_iterator sit;
    double maxDepth, thisDepth;
    maxDepth = DBL_MAX;

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
    	thisDepth = sit->second.lld().altitude();
    	if(thisDepth < maxDepth){
        	maxDepth = thisDepth;
        }
    }
    return maxDepth;
}

//Find minimum square side length in simulation
double tsunamisquares::World::getMinSize() const {
	std::map<UIndex, Square>::const_iterator sit;
	double minSize, thisX, thisY;
	minSize = DBL_MAX;
	Vec<2> center;
	double vert_dist, top_dist, bottom_dist;

	for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
		center = sit->second.xy();

		vert_dist   = bg::distance(point_spheq(center[0], center[1]-_dlat/2), point_spheq(center[0], center[1]+_dlat/2));
		top_dist    = bg::distance(point_spheq(center[0]-_dlon/2, center[1]+_dlat/2), point_spheq(center[0]+_dlon/2, center[1]+_dlat/2));
		bottom_dist = bg::distance(point_spheq(center[0]-_dlon/2, center[1]-_dlat/2), point_spheq(center[0]+_dlon/2, center[1]-_dlat/2));

		if(vert_dist < minSize){
			minSize = vert_dist;
		}
		if(top_dist < minSize){
			minSize = top_dist;
		}
		if(bottom_dist < minSize){
			minSize = bottom_dist;
		}
	}
	//bg::distance returns angular distance in radians
	return minSize*EARTH_MEAN_RADIUS;
}


// Much faster n-nearest neighbor search using an RTree
std::map<double, tsunamisquares::UIndex> tsunamisquares::World::getNearest_rtree(const Vec<2> &location, const int &numNear)const {
    std::vector<UIndex>						  nIDs;
    std::map<double, UIndex>                  neighbors;
    double 									  square_dist;

    //Use RTree query to get numNear nearest neighbors
    nIDs = _square_rtree.getNearest(location, numNear);

    // Compute distance from "location" to the center of each neighbor.
    for (int i=0; i<nIDs.size(); ++i) {
		square_dist = squareCenter(nIDs[i]).dist(location);
		neighbors.insert(std::make_pair(square_dist, nIDs[i]));
	}
    
    return neighbors;
}

std::vector<tsunamisquares::UIndex> tsunamisquares::World::getRingIntersects_rtree(const ring_spheq &ring)const {
    std::vector<UIndex>						  intersects;
    //Use RTree query to get all boxes that the provided ring intersects (edges and contained within)
    intersects = _square_rtree.getRingIntersects(ring);

    return intersects;
}

std::vector<tsunamisquares::UIndex> tsunamisquares::World::getBoxIntersects_rtree(const box_spheq &box)const {
    std::vector<UIndex>						  intersects;
    //Use RTree query to get all boxes that the provided box intersects (edges and contained within)
    intersects = _square_rtree.getBoxIntersects(box);

    return intersects;
}



// Get the square_id for each closest square to some location = (x,y)
tsunamisquares::UIndex tsunamisquares::World::whichSquare(const Vec<2> &location) const {
    std::map<double, UIndex>                  square_dists;
    std::map<UIndex, Square>::const_iterator  sit;
    UIndex                               neighbor;

    // Compute distance from "location" to the center of each square.
    // Since we use a map, the distances will be ordered since they are the keys
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        double square_dist = squareCenter(sit->first).dist(location);
        square_dists.insert(std::make_pair(square_dist, sit->second.id()));
    }
    
    // Return the ID of the nearest square
    return square_dists.begin()->second;
}


void tsunamisquares::World::indexNeighbors(void) {
	computeNeighbors();
	computeNeighborCoords();
	computeInvalidDirections();
}

void tsunamisquares::World::computeNeighbors(void) {
    std::map<UIndex, Square>::iterator                                  sit;
    double                                              this_lat, this_lon; 
    bool                                    isMinLat, isMinLon, isMaxLat, isMaxLon;
    UIndex                       this_id, left, right, top_right, top_left;
    UIndex                          bottom_left, bottom_right, top, bottom;
    
    // Use the in-place element numbering to find the IDs of the neighboring squares.
    // Must handle the border and corner cases and not include off-model neighbors.
    
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        this_id      = sit->first;
        left         = this_id-1;
        right        = this_id+1;
        top          = this_id-num_lons();
        bottom       = this_id+num_lons();
        top_left     = top-1;
        top_right    = top+1;
        bottom_left  = bottom-1;
        bottom_right = bottom+1;
        
        this_lat      = sit->second.xy()[1];
        this_lon      = sit->second.xy()[0];
        isMinLat      = (this_lat == min_lat());
        isMaxLat      = (this_lat == max_lat());
        isMinLon      = (this_lon == min_lon());
        isMaxLon      = (this_lon == max_lon());
        
        // Handle the corner and edge cases
        if (! (isMaxLat || isMaxLon || isMinLon || isMinLat)) {
            // Interior squares
            sit->second.set_right(right);
            sit->second.set_left(left);
            sit->second.set_top(top);
            sit->second.set_bottom(bottom);
            sit->second.set_top_left(top_left);
            sit->second.set_top_right(top_right);
            sit->second.set_bottom_left(bottom_left);
            sit->second.set_bottom_right(bottom_right);
        } else if (isMaxLat && isMinLon) {
            // Top left (North West) corner
            sit->second.set_right(right);
            sit->second.set_bottom(bottom);
            sit->second.set_bottom_right(bottom_right);
        } else if (isMaxLat && isMaxLon) {
            // Top right (North East) corner
            sit->second.set_left(left);
            sit->second.set_bottom(bottom);
            sit->second.set_bottom_left(bottom_left);
        } else if (isMinLat && isMaxLon) {
            // Bottom right (South East) corner
            sit->second.set_left(left);
            sit->second.set_top(top);
            sit->second.set_top_left(top_left);
        } else if (isMinLat && isMinLon) {
            // Bottom left (South West) corner
            sit->second.set_right(right);
            sit->second.set_top(top);
            sit->second.set_top_right(top_right);
        } else if (isMinLon) {
            // Left (West) border
            sit->second.set_right(right);
            sit->second.set_top(top);
            sit->second.set_bottom(bottom);
            sit->second.set_top_right(top_right);
            sit->second.set_bottom_right(bottom_right);
        } else if (isMaxLat) {
            // Top (North) border
            sit->second.set_right(right);
            sit->second.set_left(left);
            sit->second.set_bottom(bottom);
            sit->second.set_bottom_left(bottom_left);
            sit->second.set_bottom_right(bottom_right);
        } else if (isMaxLon) {
            // right (East) border
            sit->second.set_left(left);
            sit->second.set_top(top);
            sit->second.set_bottom(bottom);
            sit->second.set_top_left(top_left);
            sit->second.set_bottom_left(bottom_left);
        } else if (isMinLat) {
            // Bottom (South) border
            sit->second.set_right(right);
            sit->second.set_left(left);
            sit->second.set_top(top);
            sit->second.set_top_left(top_left);
            sit->second.set_top_right(top_right);
        } else {
            std::cout << "Error, no match to any case! (square " << this_id << ")" << std::endl;
        }

    }
    
    
}

void tsunamisquares::World::computeNeighborCoords(void) {
    std::map<UIndex, Square>::iterator     sit;

    // precompute the local coordinates of each neighbor, used in plane fitting
    for(sit=_squares.begin(); sit!=_squares.end(); ++sit){
    	std::map<UIndex, Vec<2> >   thisNeighborsAndCoords;
        SquareIDSet				    neighborIDs;
        SquareIDSet::const_iterator nit;

        neighborIDs = sit->second.get_neighbors_and_self();
    	for(nit=neighborIDs.begin(); nit!=neighborIDs.end(); nit++){
    		double xcoord = bg::distance(point_spheq(sit->second.xy()[0], sit->second.xy()[1]),
    									 point_spheq(squareCenter(*nit)[0], sit->second.xy()[1])) * EARTH_MEAN_RADIUS;
			double ycoord = bg::distance(point_spheq(sit->second.xy()[0], sit->second.xy()[1]),
										 point_spheq(sit->second.xy()[0], squareCenter(*nit)[1])) * EARTH_MEAN_RADIUS;
    	    thisNeighborsAndCoords.insert(std::make_pair(*nit, Vec<2>(xcoord, ycoord)));
    	}

    	sit->second.set_local_neighbor_coords(thisNeighborsAndCoords);
    }

}

void tsunamisquares::World::computeInvalidDirections(void) {
    std::map<UIndex, Square>::iterator sit;

    // Find directions that water in each square is not allowed to flow.

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        std::vector<bool>                  new_invalid_directions(4, false);

    	if(sit->second.left() == INVALID_INDEX || square(sit->second.left()).xyz()[2]>=0.0){
        	new_invalid_directions[0] = true;
        }
        if(sit->second.right() == INVALID_INDEX || square(sit->second.right()).xyz()[2]>=0.0){
        	new_invalid_directions[1] = true;
        }
        if(sit->second.bottom() == INVALID_INDEX || square(sit->second.bottom()).xyz()[2]>=0.0){
        	new_invalid_directions[2] = true;
        }
        if(sit->second.top() == INVALID_INDEX || square(sit->second.top()).xyz()[2]>=0.0){
        	new_invalid_directions[3] = true;
        }

        sit->second.set_invalid_directions(new_invalid_directions);
    }


}


// ----------------------------------------------------------------------
// -------------------- Single Square Functions -------------------------
// ----------------------------------------------------------------------
tsunamisquares::Vec<2> tsunamisquares::World::squareCenter(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    // (x,y) of square center
    return sit->second.xy();
}

double tsunamisquares::World::squareDepth(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    // altitude of the sea floor below this square (negative below sea level)
    return sit->second.xyz()[2];
}

double tsunamisquares::World::squareLevel(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    // altitude of the water surface for this square
    // = altitude of sea floor + height of water
    return (sit->second.xyz()[2])+(sit->second.height());
}
tsunamisquares::Vec<2> tsunamisquares::World::squareLatLon(const UIndex &square_id) const {
	std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
	return Vec<2>(sit->second.xy()[1], sit->second.xy()[0]);
}


// ----------------------------------------------------------------------
// -------------------- Functions to set initial conditions  ------------
// ----------------------------------------------------------------------
void tsunamisquares::World::setSquareVelocity(const UIndex &square_id, const Vec<2> &new_velo) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    square_it->second.set_velocity(new_velo);
}

void tsunamisquares::World::setSquareAccel(const UIndex &square_id, const Vec<2> &new_accel) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    square_it->second.set_accel(new_accel);
}

void tsunamisquares::World::setSquareHeight(const UIndex &square_id, const double &new_height) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    square_it->second.set_height(new_height);
}

// ----------------------------------------------------------------------
// -------------------- Model Building/Editing --------------------------
// ----------------------------------------------------------------------
tsunamisquares::SquareIDSet tsunamisquares::World::getSquareIDs(void) const {
    SquareIDSet square_id_set;
    std::map<UIndex, Square>::const_iterator  sit;

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        square_id_set.insert(sit->second.id());
    }

    return square_id_set;
}

tsunamisquares::Square &tsunamisquares::World::square(const UIndex &ind) throw(std::domain_error) {
    std::map<UIndex, Square>::iterator it = _squares.find(ind);

    if (it == _squares.end()) throw std::domain_error("tsunamisquares::World::square");
    else return it->second;
}

tsunamisquares::Square &tsunamisquares::World::new_square(void) {
    UIndex  max_ind = next_square_index();
    _squares.insert(std::make_pair(max_ind, Square()));
    _squares.find(max_ind)->second.set_id(max_ind);
    return _squares.find(max_ind)->second;
}

void tsunamisquares::World::clear(void) {
    _squares.clear();
}

void tsunamisquares::World::insert(tsunamisquares::Square &new_square) {
    _squares.insert(std::make_pair(new_square.id(), new_square));
}

size_t tsunamisquares::World::num_squares(void) const {
    return _squares.size();
}

void tsunamisquares::World::printSquare(const UIndex square_id) {
    Square this_square = square(square_id);

    std::cout << "~~~ Square " << this_square.id() << "~~~" << std::endl;
    std::cout << "center: " << squareCenter(square_id) << std::endl;
    std::cout << "density: " << this_square.density() << std::endl;
    std::cout << "area: " << this_square.area() << std::endl;
    if (!isnan(this_square.height())) {
        std::cout << "height: " << this_square.height() << std::endl;
        std::cout << "level: " << squareLevel(square_id) << std::endl;
        std::cout << "volume: " << this_square.volume() << std::endl;
        std::cout << "mass: " << this_square.mass() << std::endl;
        std::cout << "velocity: " << this_square.velocity() << std::endl; 
        std::cout << "accel: " << this_square.accel() << std::endl;    
        std::cout << "momentum: " << this_square.momentum() << std::endl; 
    }
}

void tsunamisquares::World::info(void) const{
    std::cout << "World: " << this->num_squares() << " squares. " << std::endl;
}

void tsunamisquares::World::get_bounds(LatLonDepth &minimum, LatLonDepth &maximum) const {
    std::map<UIndex, Square>::const_iterator    sit;
    double      min_lat, min_lon, min_alt;
    double      max_lat, max_lon, max_alt;

    min_lat = min_lon = min_alt = DBL_MAX;
    max_lat = max_lon = max_alt = -DBL_MAX;

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        min_lat = fmin(min_lat, sit->second.lld().lat());
        max_lat = fmax(max_lat, sit->second.lld().lat());
        min_lon = fmin(min_lon, sit->second.lld().lon());
        max_lon = fmax(max_lon, sit->second.lld().lon());
        min_alt = fmin(min_alt, sit->second.lld().altitude());
        max_alt = fmax(max_alt, sit->second.lld().altitude());
    }

    if (min_lat == DBL_MAX || min_lon == DBL_MAX || min_alt == DBL_MAX) {
        minimum = LatLonDepth();
    } else {
        minimum = LatLonDepth(min_lat, min_lon, min_alt);
    }

    if (max_lat == -DBL_MAX || max_lon == -DBL_MAX || max_alt == -DBL_MAX) {
        maximum = LatLonDepth();
    } else {
        maximum = LatLonDepth(max_lat, max_lon, max_alt);
    }
}


// ----------------------------------------------------------------------
// ----------------------------- Model File I/O -------------------------
// ----------------------------------------------------------------------
std::string tsunamisquares::ModelIO::next_line(std::istream &in_stream) {
    std::string line = "";
    size_t      pos;

    do {
        std::getline(in_stream, line);
        _comment = "";
        // Cut off any initial whitespace
        pos = line.find_first_not_of(" \t");

        if (pos != std::string::npos) line = line.substr(pos, std::string::npos);

        // Comment consists of hash mark until the end of the line
        pos = line.find("#");

        if (pos != std::string::npos) _comment = line.substr(pos, std::string::npos);

        // Extract the non-comment part of the line
        line = line.substr(0, line.find("#"));

        // If the line is empty, we keep going
        if (line.length() > 0) break;
    } while (in_stream && !in_stream.eof());

    return line;
}

void tsunamisquares::ModelIO::next_line(std::ostream &out_stream) const {
    if (!_comment.empty()) out_stream << " # " << _comment;

    out_stream << "\n";
}

void tsunamisquares::World::write_square_ascii(std::ostream &out_stream, const double &time, const UIndex &square_id) const {
    unsigned int        i;
    std::map<UIndex, Square>::const_iterator square_it = _squares.find(square_id);
    double waterLevel, waterHeight;

    out_stream << time << "\t";

    //
    out_stream << square_it->second.xy()[1] << "\t\t" << square_it->second.xy()[0] << "\t\t";

    // Don't write water level for the hi and dry squares until they take on water
    waterLevel  = square_it->second.height() + square_it->second.xyz()[2];
    waterHeight = square_it->second.height();
    if (waterHeight == 0.0 && waterLevel >= 0.0) {
        out_stream << waterHeight << "\t\t";
    } else {
        out_stream << waterLevel << "\t\t";
    }
    
    // Write the altitude of the bottom too
    out_stream << square_it->second.xyz()[2] << "\t\t";

    next_line(out_stream);
}

int tsunamisquares::World::read_bathymetry(const std::string &file_name) {
    std::ifstream   in_file;
    UIndex          i, j, num_squares, num_vertices, num_lats, num_lons;
    LatLonDepth     min_latlon, max_latlon;
    float			Lx_tot = 0.0;
	float			Ly_tot = 0.0;
	float			dlon, dlat;

    // Clear the world first to avoid incorrectly mixing indices
    clear();

    in_file.open(file_name.c_str());

    if (!in_file.is_open()) return -1;

    // Read the first line of metadata
    std::stringstream desc_line(next_line(in_file));
    desc_line >> num_lats;
    desc_line >> num_lons;
    desc_line >> dlat;
	desc_line >> dlon;
    _num_latitudes = num_lats;
    _num_longitudes = num_lons;
    _dlat = dlat;
	_dlon = dlon;

    // Set the number of squares
    num_squares = num_lats*num_lons;

    // Read squares, populate world rtree and square map.
    for (i=0; i<num_squares; ++i) {
        Square     new_square;
        new_square.set_id(i);
        new_square.read_bathymetry(in_file);
        new_square.set_box(dlon, dlat);
        _square_rtree.addBox(new_square.box(), i);
        _squares.insert(std::make_pair(new_square.id(), new_square));
    }

    in_file.close();
        
    // Get world lld bounds
    get_bounds(min_latlon, max_latlon);
    min_latlon.set_altitude(0); //Why is this here?

    // Keep track of Lat/Lon bounds in the World
    _min_lat = min_latlon.lat();
    _min_lon = min_latlon.lon();
    _max_lat = max_latlon.lat();
    _max_lon = max_latlon.lon();

    return 0;
}


int tsunamisquares::World::write_file_kml(const std::string &file_name) {
    std::ofstream                             out_file;
    std::map<UIndex, Square>::const_iterator  sit;
    LatLonDepth                               min_bound, max_bound, center;
    Vec<3>                                    min_xyz, max_xyz;
    double                                    dx, dy, range, Lx, Ly;
    unsigned int                              i;
    double                                    depth = 100; //So the squares are off the surface a bit

    out_file.open(file_name.c_str());

    get_bounds(min_bound, max_bound);
    center = LatLonDepth(max_bound.lat() - (max_bound.lat()-min_bound.lat())/2,
                         max_bound.lon() - (max_bound.lon()-min_bound.lon())/2);
    Conversion c(center);
    min_xyz = c.convert2xyz(min_bound);
    max_xyz = c.convert2xyz(max_bound);
    dx = max_xyz[0]-min_xyz[0];
    dy = max_xyz[1]-min_xyz[1];
    range = fmax(dx, dy) * (1.0/tan(c.deg2rad(30)));

    out_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out_file << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    out_file << "<Document>\n";
    out_file << "<LookAt>\n";
    out_file << "\t<latitude>" << center.lat() << "</latitude>\n";
    out_file << "\t<longitude>" << center.lon() << "</longitude>\n";
    out_file << "\t<altitude>0</altitude>\n";
    out_file << "\t<range>" << range << "</range>\n";
    out_file << "\t<tilt>0</tilt>\n";
    out_file << "\t<heading>0</heading>\n";
    out_file << "\t<altitudeMode>absolute</altitudeMode>\n";
    out_file << "</LookAt>\n";
    out_file << "<Style id=\"sectionLabel\">\n";
    out_file << "\t<IconStyle>\n";
    out_file << "\t\t<Icon>\n";
    out_file << "\t\t\t<href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href>\n";
    out_file << "\t\t</Icon>\n";
    out_file << "\t</IconStyle>\n";
    out_file << "</Style>\n";

    out_file << "<Folder id=\"squares\">\n";

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        // Compute the lat/lon/depth of the 4 corners of the square
        LatLonDepth         lld[4];
        Vec<2>              v2, centerXY;
        Vec<3>              v3;
        LatLonDepth         centerLLD, base;

        base        = getBase();
        centerLLD   = sit->second.lld();
        centerXY    = squareCenter(sit->first);    
        Conversion  c(base);
        Lx          = _dlon;
        Ly          = _dlat;
        // Locate the corners in XYZ, then convert to LLD
        v3      = Vec<3>(centerXY[0]-Lx/2.0, centerXY[1]+Ly/2, 0.0); // top left
        lld[0]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]-Lx/2.0, centerXY[1]-Ly/2, 0.0); // bottom left
        lld[1]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]+Lx/2.0, centerXY[1]-Ly/2, 0.0); // bottom right
        lld[2]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]+Lx/2.0, centerXY[1]+Ly/2, 0.0); // top left
        lld[3]  = c.convert2LatLon(v3);
        
        // Output the KML format polygon for this section
        out_file << "\t\t<Placemark>\n";
        out_file << "\t\t<description>\n";
        out_file << "Square: " << sit->first << "\n";
        out_file << "LLD: " << squareLatLon(sit->first)[0] << "," << squareLatLon(sit->first)[1] << "," << squareDepth(sit->first) << " [m]\n";
        out_file << "XYZ: " << squareCenter(sit->first)[0] << "," << squareCenter(sit->first)[1] << ","   << squareDepth(sit->first) << " [m]\n";
        out_file << "Area: " << sit->second.area()*pow(10,-6) << "[km^2]\n";
        out_file << "Density: " << sit->second.density() << "[kg/m^3]\n";
        out_file << "\t\t</description>\n";
        out_file << "\t\t\t<styleUrl>#baseStyle</styleUrl>\n";
        out_file << "\t\t\t<Polygon>\n";
        out_file << "\t\t\t\t<extrude>0</extrude>\n";
        out_file << "\t\t\t\t<altitudeMode>relativeToGround</altitudeMode>\n";
        out_file << "\t\t\t\t<outerBoundaryIs>\n";
        out_file << "\t\t\t\t\t<LinearRing>\n";
        out_file << "\t\t\t\t\t\t<coordinates>\n";

        for (unsigned int i=0; i<4; ++i) out_file << "\t\t\t\t\t\t\t" << lld[i].lon() << "," << lld[i].lat() << "," << depth << "\n";

        out_file << "\t\t\t\t\t\t</coordinates>\n";
        out_file << "\t\t\t\t\t</LinearRing>\n";
        out_file << "\t\t\t\t</outerBoundaryIs>\n";
        out_file << "\t\t\t</Polygon>\n";
        out_file << "\t\t</Placemark>\n";
    }

    out_file << "</Folder>\n";
    out_file << "</Document>\n";
    out_file << "</kml>\n";

    out_file.close();

    return 0;
}

int tsunamisquares::World::deformFromFile(const std::string &file_name) {
    std::ifstream   in_file;
    UIndex          i, num_points, mappedID;
    double          dz;
    Vec<2>          location;
    std::map<UIndex, Square>::iterator sit;
    std::map<double, UIndex>		   nearestMap;
    LatLonDepth     square_lld;

    in_file.open(file_name.c_str());

    if (!in_file.is_open()) return -1;

    // Read the first line describing the number of sections, etc
    std::stringstream desc_line(next_line(in_file));
    desc_line >> num_points;

    // Read the points, find nearest square, deform the bottom
    for (i=0; i<num_points; ++i) {
        Square     new_square;
        new_square.read_bathymetry(in_file);
        
        // Get location (x,y) for the lat/lon point and get the altitude change
        location = new_square.xy();
        dz  = new_square.xyz()[2];
        
        // Find the closest square, grab its vertex
        // Use the RTree here as well for speed
        nearestMap = getNearest_rtree(new_square.xy(), 1);
        mappedID   = nearestMap.begin()->second;

        sit = _squares.find( mappedID );
        
        // Get the current LLD data for this closest square
        square_lld = sit->second.lld();
        
        // Update the altitude of the vertex by the amount dz
        square_lld.set_altitude(square_lld.altitude() + dz);
        
        // Set the new position
        sit->second.set_lld(square_lld);
        
    }

    in_file.close();

    return 0;
}
void tsunamisquares::Square::read_bathymetry(std::istream &in_stream) {
    std::stringstream   ss(next_line(in_stream));

    ss >> _data._lat;
    ss >> _data._lon;
    ss >> _data._alt;
    _pos = Vec<3>(_data._lon, _data._lat, _data._alt);

}




