# here we access the root configuration, include files, and libraries
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMathMore
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  = $(ROOTCFLAGS) -I$(ODELIB) -g #-Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS)
GXX	   = g++ $(CXXFLAGS)


a11: sobol


sobol: sobol.cpp 
	$(GXX) $(CXXFLAGS) -o sobol sobol.cpp $(LDFLAGS)


clean:
	rm -f sobol
	rm -f *~ *.d *.so *.pcm 

cleanall: clean
	rm -f *png *pdf 
