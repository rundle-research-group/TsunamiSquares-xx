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

#define UNDEFINED_ELEMENT_ID    UINT_MAX

#include <string.h>
#include <fstream>
#include <map>
#include <stdlib.h>

#include "TsunamiGlobeUtil.h"

namespace tsunamisquares {
    typedef unsigned int UIndex;
    static const UIndex INVALID_INDEX = std::numeric_limits<unsigned int>::max();
    //// INVALID_INDEX = 4294967295;
    
    class ModelIO {
        private:
            std::string         _comment;

        protected:
            std::string comment(void) const {
                return _comment;
            };
            void set_comment(std::string new_comment) {
                _comment = new_comment;
            };
            std::string next_line(std::istream &in_stream);
            void next_line(std::ostream &out_stream) const;
    };

    struct FieldDesc {
        std::string name;
        std::string details;
    };

    typedef struct FieldDesc FieldDesc;
    typedef std::set<UIndex> SquareIDSet;

    // Squares, the functional members of Tsunami Square
    struct SquareData {
        UIndex              _id;
        double              _lat, _lon, _alt, _area;
        Vec<2>              _velocity;
        Vec<2>              _accel;
        double              _height;
        double              _friction;
        double              _density;
        bool                _invalid_directions[4];

        // updated data are used to keep track of volume/momentum redistribution from other squares
        double               _updated_height;
        Vec<2>               _updated_momentum;
    };

    class Square : public ModelIO {
        private:
            SquareData     _data;
            //SquareIDSet    _neighbors;
            UIndex         _top, _bottom, _right, _left, _top_right, _top_left, _bottom_left, _bottom_right;

            Vec<3> 	       _pos;

            box_spheq      _box;


        public:
            Square(void) {
                _data._id = INVALID_INDEX;
                _data._lat = _data._lon = _data._alt = _data._area = std::numeric_limits<double>::quiet_NaN();
                _pos = Vec<3>();
                _data._velocity = _data._accel = _data._updated_momentum = Vec<2>(0.0,0.0);
                _data._height = _data._updated_height = std::numeric_limits<double>::quiet_NaN();
                _data._density = 1025.0; // sea water by default
                _data._friction = 0.02;
                
                for(int i=0;i<4;i++){
                	_data._invalid_directions[i] = 0;
                }

                _left = _right = _top = _bottom = INVALID_INDEX;
                _top_right = _top_left = _bottom_left = _bottom_right = INVALID_INDEX;
                
            };
            
            void clear(void);

            SquareData data(void) const {
                return _data;
            };

            UIndex id(void) const {
                return _data._id;
            };
            void set_id(const UIndex &id) {
                _data._id = id;
            };

            double height(void) const {
                return _data._height;
            };
            void set_height(const double &new_height) {
                _data._height = new_height;
            };
            
            double updated_height(void) const {
                return _data._updated_height;
            };
            void set_updated_height(const double &new_height) {
                _data._updated_height = new_height;
            };
            
            double density(void) const {
                return _data._density;
            };
            void set_density(const double &new_density) {
                _data._density = new_density;
            };
            
            Vec<2> velocity(void) const {
                return _data._velocity;
            };
            void set_velocity(const Vec<2> &new_velocity) {
                _data._velocity = new_velocity;
            };
            
            Vec<2> updated_momentum(void) const {
                return _data._updated_momentum;
            };
            void set_updated_momentum(const Vec<2> &new_momentum) {
                _data._updated_momentum = new_momentum;
            };
            
            Vec<2> accel(void) const {
                return _data._accel;
            };
            void set_accel(const Vec<2> &new_accel) {
                _data._accel = new_accel;
            }
            double friction(void) const {
                return _data._friction;
            }
            void set_friction(const double &new_friction) {
                _data._friction = new_friction;
            }
            void set_box(const double &dlon, const double &dlat) {
            	_box = box_spheq(point_spheq(_pos[0]-dlon, _pos[1]-dlat), point_spheq(_pos[0]+dlon, _pos[1]+dlat));
            	// set area when we set box
            	_data._area = bg::area(ring())*EARTH_MEAN_RADIUS*EARTH_MEAN_RADIUS;
            }
            box_spheq box(void) const {
            	return _box;
            }
            ring_spheq ring(void) const {
            	ring_spheq box_ring;
				point_spheq boxverts[5];
				boxverts[0] =_box.min_corner();
				boxverts[1] = point_spheq(_box.min_corner().get<0>(), _box.max_corner().get<1>());
				boxverts[2] =_box.max_corner();
				boxverts[3] = point_spheq(_box.max_corner().get<0>(), _box.min_corner().get<1>());
				boxverts[4] =_box.min_corner();
				bg::assign_points(box_ring, boxverts);
            	return box_ring;
            }
            double area(void) const {
            	return _data._area;
			}

            double volume(void) const {
                return area()*_data._height;
            };
            
            double mass(void) const {
                return _data._density*volume();
            };
            
            Vec<2> momentum(void) const {
                return _data._velocity*mass();
            };
            
            bool* invalid_directions(void) {
            	return _data._invalid_directions;
            }

            void set_invalid_directions(bool* new_invalid_directions) {
            	for(int i=0;i<4;i++){
            		_data._invalid_directions[i] = new_invalid_directions[i];
            	}
            }

            //  All functions with top/bottom/left/right are setting the IDs of the corresponding
            //     neighboring squares. Left/Right/Top/Bottom are nearest neighbors, the others are
            //     next-nearest neighbors.
            
            void set_top(const UIndex &id) {
                _top = id;
            };
            UIndex top(void) const {
                return _top;
            };
            
            void set_right(const UIndex &id) {
                _right = id;
            };
            UIndex right(void) const {
                return _right;
            };
            
            void set_left(const UIndex &id) {
                _left = id;
            };
            UIndex left(void) const {
                return _left;
            };
            
            void set_bottom(const UIndex &id) {
                _bottom = id;
            };
            UIndex bottom(void) const {
                return _bottom;
            };
            
            void set_top_left(const UIndex &id) {
                _top_left = id;
            };
            UIndex top_left(void) const {
                return _top_left;
            };
            
            void set_top_right(const UIndex &id) {
                _top_right = id;
            };
            UIndex top_right(void) const {
                return _top_right;
            };
            
            void set_bottom_left(const UIndex &id) {
                _bottom_left = id;
            };
            UIndex bottom_left(void) const {
                return _bottom_left;
            };
            
            void set_bottom_right(const UIndex &id) {
                _bottom_right = id;
            };
            UIndex bottom_right(void) const {
                return _bottom_right;
            };
            
            SquareIDSet get_valid_neighbors(void) const {
                SquareIDSet valid_neighbors;
                if (left()        != INVALID_INDEX) valid_neighbors.insert(left());
                if (right()       != INVALID_INDEX) valid_neighbors.insert(right());
                if (top()         != INVALID_INDEX) valid_neighbors.insert(top());
                if (bottom()      != INVALID_INDEX) valid_neighbors.insert(bottom());
                if (top_left()    != INVALID_INDEX) valid_neighbors.insert(top_left());
                if (bottom_left() != INVALID_INDEX) valid_neighbors.insert(bottom_left());
                if (top_right()   != INVALID_INDEX) valid_neighbors.insert(top_right());
                if (bottom_right()!= INVALID_INDEX) valid_neighbors.insert(bottom_right());
                return valid_neighbors;
            };
            
            SquareIDSet get_valid_nearest_neighbors(void) const {
                SquareIDSet valid_neighbors;
                if (left()        != INVALID_INDEX) valid_neighbors.insert(left());
                if (right()       != INVALID_INDEX) valid_neighbors.insert(right());
                if (top()         != INVALID_INDEX) valid_neighbors.insert(top());
                if (bottom()      != INVALID_INDEX) valid_neighbors.insert(bottom());
                return valid_neighbors;
            };
            
            SquareIDSet get_valid_nextnearest_neighbors(void) const {
				SquareIDSet valid_neighbors;
				if (top_left()    != INVALID_INDEX) valid_neighbors.insert(top_left());
				if (bottom_left() != INVALID_INDEX) valid_neighbors.insert(bottom_left());
				if (top_right()   != INVALID_INDEX) valid_neighbors.insert(top_right());
				if (bottom_right()!= INVALID_INDEX) valid_neighbors.insert(bottom_right());
				return valid_neighbors;
			};

            SquareIDSet get_neighbors_and_self(void) const {
                SquareIDSet valid_neighbors;
                if (left()        != INVALID_INDEX) valid_neighbors.insert(left());
                if (right()       != INVALID_INDEX) valid_neighbors.insert(right());
                if (top()         != INVALID_INDEX) valid_neighbors.insert(top());
                if (bottom()      != INVALID_INDEX) valid_neighbors.insert(bottom());
                if (top_left()    != INVALID_INDEX) valid_neighbors.insert(top_left());
                if (bottom_left() != INVALID_INDEX) valid_neighbors.insert(bottom_left());
                if (top_right()   != INVALID_INDEX) valid_neighbors.insert(top_right());
                if (bottom_right()!= INVALID_INDEX) valid_neighbors.insert(bottom_right());
                valid_neighbors.insert(id());
                return valid_neighbors;
            };
            
            void print_neighbors(void) {
                std::cout << "--- #" << id() << " ------------" << std::endl;
                std::cout << " left:\t\t" << left() << std::endl;
                std::cout << " right:\t\t" << right() << std::endl;
                std::cout << " top:\t\t" << top() << std::endl;
                std::cout << " bottom:\t" << bottom() << std::endl;
                std::cout << " t left:\t" << top_left() << std::endl;
                std::cout << " t right:\t" << top_right() << std::endl;
                std::cout << " b left:\t" << bottom_left() << std::endl;
                std::cout << " b right:\t" << bottom_right() << std::endl;
            }

            LatLonDepth lld(void) const {
				return LatLonDepth(_data._lat, _data._lon, _data._alt);
			};
			void set_lld(const LatLonDepth &lld) {
				_data._lat = lld.lat();
				_data._lon = lld.lon();
				_data._alt = lld.altitude();
				_pos = Vec<3>(_data._lon, _data._lat, _data._alt);
				// TEMPORARY FIX TO LET (X,Y,Z) = (LON,LAT,ALT)
				//_pos = Vec<3>(lld.lon(), lld.lat(), lld.altitude());
			};

			Vec<3> xyz(void) const {
				return _pos;
			};
			Vec<2> xy(void) const {
				Vec<2> pos_xy;
				pos_xy[0] = _pos[0];
				pos_xy[1] = _pos[1];
				return pos_xy;
			};

            void read_bathymetry(std::istream &in_stream);

//            void add_neighbor(const UIndex &neighbor_id) {
//                _neighbors.insert(neighbor_id);
//            };
//            
//            SquareIDSet get_neighbors(void) {
//                return _neighbors;
//            };
    
    };
    
    
    // Class to contain all Squares and Bathymetry 
    class World : public ModelIO {
        private:
            std::map<UIndex, Square>  _squares;
            LatLonDepth _base;
            double _min_lat, _max_lat, _min_lon, _max_lon, _dlat, _dlon;
            int _num_latitudes, _num_longitudes;
            RTree_spheq _square_rtree;
            // Diffusion constant
            static const double _D = 140616.45;//100000;
        
        public:
            Square &new_square(void);
            
            Square &square(const UIndex &ind) throw(std::domain_error);
            
            UIndex next_square_index(void) const {
                if (_squares.size()) return _squares.rbegin()->first+1;
                else return 0;
            };
            
            size_t num_squares(void) const;
            size_t num_vertices(void) const;
            
            void insert(Square &new_square);
            
            void clear(void);
            
            LatLonDepth getBase(void) const {
                return _base;
            }
            
            int num_lats(void) const {
                return _num_latitudes;
            }
            int num_lons(void) const {
                return _num_longitudes;
            }
            
            double D(void) const {
                return _D;
            }
            double min_lat(void) const {
                return _min_lat;
            }
            double max_lat(void) const {
                return _max_lat;
            }
            double min_lon(void) const {
                return _min_lon;
            }
            double max_lon(void) const {
                return _max_lon;
            }
            
            void printSquare(const UIndex square_id);
            void printVertex(const UIndex vertex_id);
            void info(void) const;
            
            int read_bathymetry(const std::string &file_name);
            
            void get_bounds(LatLonDepth &minimum, LatLonDepth &maximum) const;
            void reset_base_coord(const LatLonDepth &new_base);
            
            SquareIDSet getSquareIDs(void) const;
            SquareIDSet getVertexIDs(void) const;


            std::map<double, UIndex> getNearest_rtree(const Vec<2> &location, const int &numNear) const;
            std::vector<UIndex> getRingIntersects_rtree(const ring_spheq &ring) const;
            std::vector<UIndex> getBoxIntersects_rtree(const box_spheq &box) const;
            SquareIDSet getNeighborIDs(const UIndex &square_id) const;
            
            /// Test ///
            

            // ======= Main functions =========
            void computeNeighbors(void);
            void computeInvalidDirections(void);
            void fillToSeaLevel(void);
            void moveSquares(const double dt, const double frac_threshold, const bool accel_bool);
            void diffuseSquares(const double dt, const double D);
            Vec<2> getGradient(const UIndex &square_id) const;
            Vec<2> fitPointsToPlane(const SquareIDSet &square_ids);
            tsunamisquares::Vec<2> getGradient_planeFit(const UIndex &square_id);
            void updateAcceleration(const UIndex &square_id);
            void deformBottom(const UIndex &square_id, const double &height_change);
            UIndex whichSquare(const Vec<2> &location) const;
            void flattenBottom(const double &depth);
            void bumpCenter(const double bump_height);
            double getMaxDepth(void) const;
            double getMinSize(void) const;
            // ======= Square functions =========
            Vec<2> squareCenter(const UIndex &square_id) const;
            Vec<2> squareLatLon(const UIndex &square_id) const;
            double squareDepth(const UIndex &square_id) const;
            double squareLevel(const UIndex &square_id) const;
            SquareIDSet get_neighbors_for_accel(const UIndex &square_id) const;
            // ======= Initial condition setup functions ======
            void setSquareVelocity(const UIndex &square_id, const Vec<2> &new_velo);
            void setSquareAccel(const UIndex &square_id, const Vec<2> &new_accel);
            void setSquareHeight(const UIndex &square_id, const double &new_height);
            // ======= File I/O ====================
            int write_file_kml(const std::string &file_name);
            void write_square_ascii(std::ostream &out_stream, const double &time, const UIndex &square_id) const;
            int deformFromFile(const std::string &file_name);
    };
}
