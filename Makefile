# here we access the root configuration, include files, and libraries
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMathMore
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  = $(ROOTCFLAGS) -I$(ODELIB) -g #-Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS)
GXX	   = g++ $(CXXFLAGS)


a11: sobel


sobel: sobel.cpp 
	$(GXX) $(CXXFLAGS) -o sobel sobel.cpp $(LDFLAGS)


clean:
	rm -f sobel
	rm -f *~ *.d *.so *.pcm 

cleanall: clean
	rm -f *png *pdf 
