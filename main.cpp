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
     *  suseptability to errors
    const std::string   out_file_name    = "tsunami_output.txt"; // name this relevant
    const std::string   bathy_file       = "bathymetry/Pacific_900.txt"; // name this relevant
    //const std::string   kml_file         = "local/Pacific_36.kml";
    // const std::string   deformation_file = "local/Channel_Islands_test_bump.txt"; // kasey says we dont need this now
    // until later
    
    // Diffusion constant (fit to a reasonable looking sim)
    double D = 140616.45; //140616.45;
    // Flattening the bathymetry to a constant depth (negative for below sea level)
    double new_depth = -100.0;
    // Bumping up the bottom
    double bump_height = 50.0;
    // Number of times to move squares
    int N_steps = 20; //number of time steps 10 is fine, to see a bit of movement
    // because boundaries aren't defined very well, we limit the time steps whenever the water hits the walls of things
    // Updating intervals, etc.
    int current_step = 0;
    int update_step = 1;
    int save_step = 1;
    double time = 0.0;
    int output_num_digits_for_percent = 3;
    */

    // Ensure we are given the parameter file name
	assertThrow(argc == 2, "usage: param_file");

	param_file.open(argv[argc-1]);


	std::string 				param_name;
	std::string					value;
	std::vector<std::string>	param_values;

	while ( param_file >> param_name >> value )
	{
	  param_values.push_back(value);
	}

    const std::string   out_file_name    	= param_values[0];
    const std::string   bathy_file       	= param_values[1];
    const std::string   kml_file         	= param_values[2];
    const std::string   deformation_file 	= param_values[3];

    bool	move_bool						= atof(param_values[4].c_str());
	bool	accel_bool						= atof(param_values[5].c_str());
	bool	diffuse_bool					= atof(param_values[6].c_str());
    // Diffusion constant (fit to a reasonable looking sim)
    double 	D 								= atof(param_values[7].c_str()); //140616.45;
    // Time step in seconds
    double  dt_param						= atof(param_values[8].c_str());
    // Number of times to move squares
    int 	N_steps 						= atof(param_values[9].c_str()); //number of time steps 10 is fine, to see a bit of movement
    // because boundaries aren't defined very well, we limit the time steps whenever the water hits the walls of things
    // Updating intervals, etc.
    int 	current_step 					= atof(param_values[10].c_str());
    int 	update_step 					= atof(param_values[11].c_str());
    int 	save_step 						= atof(param_values[12].c_str());
    double 	time 							= atof(param_values[13].c_str());
    int 	output_num_digits_for_percent 	= atof(param_values[14].c_str());
    //Boolean to decide whether to flatten the seafloor before running, for testing purposes
    bool	flatten_bool					= atof(param_values[15].c_str());
    // Flattening the bathymetry to a constant depth (negative for below sea level)
	double 	flat_depth 						= atof(param_values[16].c_str());
    //Boolean to decide whether to use bump instead deformation file, for testing purposes
	bool	bump_bool						= atof(param_values[17].c_str());
    // How high the central bump should be
	double 	bump_height 					= atof(param_values[18].c_str());


    // Header for the simulation output
    const std::string   header = "# time \t lon \t\t lat \t\t water height \t altitude \n";
    
    // -------------------------------------------------------------------------------- //
    ///////                Simulation Initialization and Loading                   ///////
    // -------------------------------------------------------------------------------- //
    start = clock();
    // Read in the bathymetry data
    this_world.clear();
    std::cout << std::endl << "Reading..."   << bathy_file.c_str() << std::endl;
    this_world.read_bathymetry(bathy_file.c_str());
    // Index the neighbors by left/right/top etc.
    std::cout << "Indexing neighbors......" << std::endl;
    this_world.indexNeighbors();

    // Flatten the bottom for simple simulation test cases, do not do this for tsunami simulations
	if(flatten_bool){
		std::cout << "Flattening the bottom..."<< std::endl;
		this_world.flattenBottom(flat_depth);
	}

    // Compute the time step given the diffusion constant D
    //double dt = (double) (int) this_world.square(0).Lx()*this_world.square(0).Ly()/(2*D); //seconds
    // Use file-provided time step, Check time step for minimum against shallow water wave speed.
	double dt;
	double G = 9.80665;
    double wave_speed = sqrt(abs(G*this_world.getMaxDepth()));
    double minDt = this_world.getMinSize() / wave_speed;

    if(0){//dt_param < minDt){
    	dt = minDt;
    	std::cout << "Using minimum time step of " << dt <<" seconds..." << std::endl;
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
    
    
    // Put water into squares to bring water level up to sealevel.
    std::cout << "Filling with water..." << std::endl;
	this_world.fillToSeaLevel();

    // --------------------------------------------------------------------------------//
    //            Sea Floor Deformation and Initial Conditions                         //
    // --------------------------------------------------------------------------------//
    std::cout << "Deforming the bottom... " << std::endl;

    if(bump_bool){
    	//   == DEFORM A LAND BUMP =======
    	std::cout << "\t Deforming central bump " << std::endl;
    	this_world.bumpCenter(bump_height);
    }else{
    	//   == DEFORM FROM FILE ==
    	std::cout << "\t Deforming from file" << std::endl;

		this_world.deformFromFile(deformation_file);

    }


    // --------------------------------------------------------------------------------//
    // --==                         File I/O Preparation                          --== //
    // --------------------------------------------------------------------------------//
    out_file.open(out_file_name.c_str());
    out_file << header.c_str();
    std::cout.precision(output_num_digits_for_percent);



    // --------------------------------------------------------------------------------//
    // --========-           Begin the Simulation; Move the Squares          ----====- //
    // --------------------------------------------------------------------------------//
    std::cout << "Moving squares....time_step=" <<dt << "...";
    while (time < max_time) {
        // If this is a writing step, print status
        if (current_step%update_step == 0) {
            std::cout << ".." << (100.0*current_step)/N_steps << "%..";
            std::cout << std::flush;
        }

        // Write the current state to file
        if (current_step%save_step == 0) {
            for (it=ids.begin(); it!=ids.end(); ++it) {
                this_world.write_square_ascii(out_file, time, *it);
            }
        }
        // Move the squares
        if(move_bool) {
        	this_world.moveSquares(dt, accel_bool);
        }
        if(diffuse_bool) {
        	this_world.diffuseSquares(dt, D); //smoothing operation
        }
        time += dt;
        current_step += 1;
    }
    // Write the final state to the file
    for (it=ids.begin(); it!=ids.end(); ++it) {
		this_world.write_square_ascii(out_file, time, *it);
	}
    out_file.close();


    // --------------------------------------------------------------------------------//
    // --========---                    Wrap up and Reporting            ---=======--- //
    // --------------------------------------------------------------------------------//
    std::cout << std::endl << "Results written to " << out_file_name << std::endl;
    end = clock();
    std::cout.precision(2+output_num_digits_for_percent);
    std::cout << "Total time: " << (float(end)-float(start))/CLOCKS_PER_SEC << " secs." << std::endl << std::endl;
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

// Grab all squares and print
//    ids = this_world.getSquareIDs();
//    
//    for (it=ids.begin(); it!=ids.end(); ++it) {
//        this_world.printSquare(*it);
//    }


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

