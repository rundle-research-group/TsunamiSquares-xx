// Copyright (c) 2017 John M. Wilson, Kasey W. Schultz
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
#include <time.h>

#define assertThrow(COND, ERR_MSG) assert(COND);

int main (int argc, char **argv) {

	// Initialize the world (where the squares live), squares and vertices
    tsunamisquares::World                       this_world;
    tsunamisquares::SquareIDSet::const_iterator it;
    tsunamisquares::SquareIDSet                 ids;
    std::ifstream								param_file;
    std::ofstream                               out_file;
    clock_t                                     start,end;


    // -------------------------------------------------------------------------------- //
    ///////////          CONSTANTS         ////////////
    // -------------------------------------------------------------------------------- //
    /*  Wilson: Trying a simple input file for these parameters.  Can be improved in the future for readability and
     *  suseptability to errors*/

    // Ensure we are given the parameter file name
	assertThrow(argc == 2, "usage: param_file");

	param_file.open(argv[argc-1]);

	std::string 				param_name;
	std::string					value;
	std::map<std::string, std::string>	param_values;

	while ( param_file >> param_name >> value )
	{
		param_values[param_name]= value;

	}

    const std::string   out_file_name    	= param_values["out_file_name"];
    const std::string   bathy_file       	= param_values["bathy_file"];
    const std::string   kml_file         	= param_values["kml_file"];
    const std::string   deformation_file 	= param_values["deformation_file"];

    // Turn on or off movement of squares, calculating of accelerations, or diffusion
    bool	move_bool						= atof(param_values["move_bool"].c_str());
	bool	accel_bool						= atof(param_values["accel_bool"].c_str());
	bool	diffuse_bool					= atof(param_values["diffuse_bool"].c_str());

	//number of Ward diffusion sweeps per time step
	int		ndiffusions						= atof(param_values["ndiffusions"].c_str());
    // Diffusion constant for Schultz diffusion
    double 	D 								= atof(param_values["D"].c_str()); //140616.45;

    // Time step in seconds
    double  dt_param						= atof(param_values["dt"].c_str());
    // Number of times to move squares
    int 	N_steps 						= atof(param_values["N_steps"].c_str());
    // Updating intervals, etc.
    int 	current_step 					= atof(param_values["current_step"].c_str());
    int 	update_step 					= atof(param_values["update_step"].c_str());
    int 	save_step 						= atof(param_values["save_step"].c_str());
    double 	time 							= atof(param_values["time"].c_str());
    int 	output_num_digits_for_percent 	= atof(param_values["output_num_digits_for_percent"].c_str());


    std::string initial_conditions			= param_values["initial_conditions"];
    //Boolean to decide whether to flatten the seafloor before running, for testing purposes
    bool	flatten_bool					= atof(param_values["flatten_bool"].c_str());
    // Flattening the bathymetry to a constant depth (negative for below sea level)
	double 	flat_depth 						= atof(param_values["flat_depth"].c_str());
    //Boolean to decide whether to use bump instead deformation file, for testing purposes
	bool	bump_bool						= atof(param_values["bump_bool"].c_str());
    // How high the central bump should be
	double 	bump_height 					= atof(param_values["bump_height"].c_str());
	// Should we use a gaussian pile as initial conditions? (plus pile height and width
	bool 	gauss_bool	 					= atof(param_values["gauss_bool"].c_str());
	double  gauss_height 					= atof(param_values["gauss_height"].c_str());
	double 	gauss_std	 					= atof(param_values["gauss_std"].c_str());

	// Use plane fitting method for calculating accelerations?  If not falls back on simpler seperate x & y linear accel calc
	bool    doPlaneFit                      = atof(param_values["doPlaneFit"].c_str());

	// Number of multiprocessing threads to run the sim on
	int     num_threads                     = atof(param_values["num_threads"].c_str());

	// Write out final simulation state out at end?
	bool  write_sim_state                   = atof(param_values["write_sim_state"].c_str());
	const std::string   finalstate_file_name= param_values["finalstate_file_name"];

	// Read in initial sim state at beginning?
	bool  read_sim_state                    = atof(param_values["read_sim_state"].c_str());
	const std::string initialstate_file_name= param_values["initialstate_file_name"];



	omp_set_num_threads(num_threads);
    
    // -------------------------------------------------------------------------------- //
    ///////                Simulation Initialization and Loading                   ///////
    // -------------------------------------------------------------------------------- //
    // Read in the bathymetry data TODO: remove max depth calc from bathy reading and move to after the bottom flattening step
    this_world.clear();
    std::cout << std::endl << "Reading..."   << bathy_file.c_str() << std::endl;
    this_world.read_bathymetry(bathy_file.c_str());

    // Flatten the bottom for simple simulation test cases, do not do this for tsunami simulations
	if(flatten_bool){
		std::cout << "Flattening the bottom..."<< std::endl;
		this_world.flattenBottom(flat_depth);
	}

    // Index the neighbors by left/right/top etc.
    std::cout << "Indexing neighbors......" << std::endl;
    this_world.indexNeighbors();

    // --------------------------------------------------------------------------------//
    //            Sea Floor Deformation and Initial Conditions                         //
    // --------------------------------------------------------------------------------//

    // Put water into squares to bring water level up to sealevel.
    std::cout << "Filling with water..." << std::endl;
	this_world.fillToSeaLevel();

	// TODO: make a switch between these options
	// Either gaussian pile, central bump, or deform from file
	std::map<std::string, int> initial_conds;
	initial_conds["eq"] = 1;
	initial_conds["bump"] = 2;
	initial_conds["gauss"]= 3;
	initial_conds["saved"]= 4;
	switch(initial_conds[initial_conditions]){
		case 1: std::cout << "Deforming from file" << std::endl;
			       this_world.deformFromFile(deformation_file);
			       break;
		case 2: std::cout << "Deforming central bump " << std::endl;
    			this_world.bumpCenter(bump_height);
    			break;
		case 3: std::cout << "Accumulating gaussian pile... " << std::endl;
				this_world.gaussianPile(gauss_height, gauss_std);
				break;
		case 4: std::cout << "Reading initial conditions from saved state file" << std::endl;
				this_world.read_sim_state_netCDF(initialstate_file_name, flatten_bool);
				break;
		default: std::cout << "Initial conditions didn't match any known.  Exiting" << std::endl;
				 return 0;

	}



    // --------------------------------------------------------------------------------//
    //                         Time and Diffusion Behavior                             //
    // --------------------------------------------------------------------------------//

    // Compute the time step given the diffusion constant D
    //double dt = (double) (int) this_world.square(0).Lx()*this_world.square(0).Ly()/(2*D); //seconds
    // Use file-provided time step, Check time step for minimum against shallow water wave speed.
	double dt;
	double G = 9.80665;
    double wave_speed = sqrt(abs(G*this_world.max_depth()));
    double maxDt = this_world.min_spacing() / wave_speed;

    if(dt_param > maxDt){
    	dt = maxDt;
    	std::cout << "Using maximum time step of " << dt <<" seconds..." << std::endl;
    }else{
    	dt = dt_param;
    	std::cout << "Using provided time step of " << dt <<" seconds..." << std::endl;
    }

    // Gather model information
    this_world.info();
    int num_lats = this_world.num_lats();
    int num_lons = this_world.num_lons();
    std::cout << "Lons by Lats = (" << num_lons << ", " << num_lats << ")...";
    ids = this_world.getSquareIDs();
    double max_time = N_steps*dt;

    // Write KML model
    //std::cout << "Writing KML..."   << kml_file.c_str() << "  ...";
    //this_world.write_file_kml(kml_file.c_str());
    

    // --------------------------------------------------------------------------------//
    // --==                         File I/O Preparation                          --== //
    // --------------------------------------------------------------------------------//
    // Header for the simulation output
    //const std::string   header = "# time \t lon \t\t lat \t\t water height \t altitude \n";
    //out_file.open(out_file_name.c_str());
    //out_file << header.c_str();
	std::cout << "Initilizing netCDF output...";
    this_world.initilize_netCDF_file(out_file_name);
    std::cout.precision(output_num_digits_for_percent);



    // --------------------------------------------------------------------------------//
    // --========-           Begin the Simulation; Move the Squares          ----====- //
    // --------------------------------------------------------------------------------//
    start = omp_get_wtime();
    bool isHealthy = true;
    std::cout << "Moving squares....time_step=" <<dt << "...";
    while (time < max_time) {
        // If this is a writing step, print status
        if (current_step%update_step == 0) {
            std::cout << ".." << current_step << "/"<<N_steps << "..";
            std::cout << std::flush;
        }
        // Check sim health, exit if there are NaNs or Infs floating around
        //isHealthy = this_world.checkSimHealth();
        //if(!isHealthy) break;

        // Write the current state to file
        if (current_step%save_step == 0) {
            //for (it=ids.begin(); it!=ids.end(); ++it) {
            //    this_world.write_square_ascii(out_file, time, *it);
            //}
        	this_world.append_netCDF_file(out_file_name, current_step, time);
        }

        // Move the squares
        if(move_bool) {
        	this_world.moveSquares(dt, accel_bool, doPlaneFit);
        }

        // Diffuse (smooth) the squares
        if(diffuse_bool) {
			//this_world.diffuseSquaresSchultz(dt, D);
			this_world.diffuseSquaresWard(ndiffusions);
        }



        // Increment time
        time += dt;
        current_step += 1;
    }
    // Write the final state to the file
    for (it=ids.begin(); it!=ids.end(); ++it) {
		this_world.write_square_ascii(out_file, time, *it);
	}
    this_world.append_netCDF_file(out_file_name, current_step, time);
    out_file.close();


    if(write_sim_state){
    	this_world.write_sim_state_netCDF(finalstate_file_name);
    }


    // --------------------------------------------------------------------------------//
    // --========---                    Wrap up and Reporting            ---=======--- //
    // --------------------------------------------------------------------------------//
    std::cout << std::endl << "Results written to " << out_file_name << std::endl;
    end = omp_get_wtime();
    std::cout.precision(2+output_num_digits_for_percent);
    std::cout << "Total time: " << (float(end)-float(start)) << " secs." << std::endl << std::endl;
    return 0;
}           /* <(^_^<) Happy Coder says: We're good up through this line!*/

/*
    //   == DEFORM A BOWL =======

    std::cout << "\nmaking a bowl shape.";

    int xDimensions; 
    int yDimensions;
    std::cout << "\ninput the x dimension: "; 
    std::cin >> xDimensions;
    std::cout << "\nnow the y dimension: ";
    std::cin >> yDimensions;

    int landRowLeft = 0; 
    int landRowRight = 0;
    for (int gradient = 0; gradient < 5; gradient++){
      for (tsunamisquares::UIndex centralDIFF = (100 + landRowLeft); centralDIFF < (110 + landRowRight); centralDIFF++){
        this_world.deformBottom(centralDIFF,bump_height - gradient*5);
      }
      landRowLeft  = landRowLeft + 30 + 1;
      landRowRight = landRowRight + 30 - 1;
    }

    landRowLeft = 120; 
    landRowRight = 120;
    for (int gradient = 0; gradient < 5; gradient++){
      for (tsunamisquares::UIndex centralDIFF = (250 + landRowLeft); centralDIFF < (260 + landRowRight); centralDIFF++){
        this_world.deformBottom(centralDIFF,bump_height - gradient*5);
      }
      landRowLeft  = landRowLeft - 30 + 1;
      landRowRight = landRowRight - 30 - 1;
    }

    tsunamisquares::UIndex centralDIFF = 130;
    int j = 8;
    for (int k = 0; k < 4; k++){
      for (int i = 0; i < j; i ++ ){
        this_world.deformBottom(centralDIFF + 30*i ,bump_height - 5*k);
      }
      centralDIFF =  centralDIFF + 31;
      j = j-2;
    }

    centralDIFF =  centralDIFF -31*4 + 9;
    j = 8;
    for (int k = 0; k < 4; k++){
      for (int i = 0; i < j; i ++ ){
        this_world.deformBottom(centralDIFF + 30*i ,bump_height - 5*k);
      }
      centralDIFF =  centralDIFF + 29;
      j = j-2;
    }
    
    //----==  DEFORM A STAIRCASE in the middle. Testing the plane fitting ======-------
    tsunamisquares::UIndex central = (int) (0.5*num_lons*(num_lats + 1));
    std::cout << " about the central square " << central << "...";
    
    double mid_bump = this_world.square(central).Lx();
    double hi_bump = 2.0*mid_bump;
    
    
    tsunamisquares::UIndex left        = this_world.square(central).left();
    tsunamisquares::UIndex right       = this_world.square(central).right();
    tsunamisquares::UIndex top         = this_world.square(central).top();
    tsunamisquares::UIndex top_left    = this_world.square(central).top_left();
    tsunamisquares::UIndex top_right   = this_world.square(central).top_right();
    tsunamisquares::UIndex bottom      = this_world.square(central).bottom();
    tsunamisquares::UIndex bottom_left = this_world.square(central).bottom_left();
    tsunamisquares::UIndex bottom_right= this_world.square(central).bottom_right();

    // Stair case is hi to the left, drops to zero on the right
    //this_world.deformBottom(right,        hi_bump);
    //this_world.deformBottom(top_right,    hi_bump);
    //this_world.deformBottom(bottom_right, hi_bump);
    
    // this_world.deformBottom(central,     mid_bump);
    // this_world.deformBottom(top,         mid_bump);
    // this_world.deformBottom(bottom,      mid_bump);
    
    //this_world.getGradient_planeFit(central);
    
    //double x_result = this_world.fitPointsToPlane(this_world.square(central).get_nearest_neighbors_and_self());
    //std::cout << "Best fit plane to " << this_world.square(central).get_nearest_neighbors_and_self().size() << " squares." << std::endl;
    //std::cout << "grabbed x = (" << x_result[0] << ", " << x_result[1] << ", " << x_result[2] << std::endl;

    //mid_bump = 10.0;
    //this_world.deformBottom(right,        mid_bump);
    //this_world.deformBottom(top_right,    mid_bump);
    this_world.deformBottom(bottom_right, mid_bump);
    this_world.deformBottom(central,      mid_bump);
    this_world.deformBottom(top,          mid_bump);
    this_world.deformBottom(bottom,       mid_bump);
    this_world.deformBottom(left,         mid_bump);
    this_world.deformBottom(top_left,     mid_bump);
    this_world.deformBottom(bottom_left,  mid_bump);
    */
//////////// Chalkboard ////////////////

    // Creating up sloping beach over the bottom 5 rows
//    assertThrow(num_lats == num_lons, "lats and lons mismatch");
//    
//    for (unsigned int i=0; i< (int) this_world.num_squares(); ++i) {
//        int row = (int)(i/num_lons);
//        if (row == num_lats-5) this_world.deformBottom(i, 50);
//        if (row == num_lats-4) this_world.deformBottom(i, 75);
//        if (row == num_lats-3) this_world.deformBottom(i, 90);
//        if (row == num_lats-2) this_world.deformBottom(i, 101);
//        if (row == num_lats-1) this_world.deformBottom(i, 110);
//        
//    }

    
//    for (unsigned int i=0; i< (int) this_world.num_squares(); ++i) {
//        // Create a wave coming from the top down, first row
//        int row = (int)(i/num_lons);
//        int col = (int)(i%num_lons);
//        if (row== num_lats-8 && col > 9 && col < num_lons-10) {
//            // Wave is 1m hi
//            this_world.deformBottom(i, 10);
//            //this_world.setSquareVelocity(i, tsunamisquares::Vec<2>(0.0, -this_world.square(0).Ly()/100));
//        }
//    }

