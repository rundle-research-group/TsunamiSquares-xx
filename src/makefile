CXX = g++
CXXFLAGS = -I/usr/local/include/boost_1_64_0 -fopenmp
ODIR = obj
LIBS = -lGeographic -lnetcdf_c++4 -lalglib

DEPS = TsunamiUtil.h TsunamiObjects.h

_OBJ = main.o TsunamiUtil.o TsunamiObjects.o
OBJ = $(addprefix $(ODIR)/,$(_OBJ))
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(LIBS) $(CXXFLAGS) 

TsunamiSquares: $(OBJ)
	$(CXX) -o $@ $^ $(LIBS) $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
