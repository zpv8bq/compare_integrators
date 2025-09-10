CPPFLAGS=-O3 -Wall
#CPPFLAGS=-g -Wall

all: gqconstants

integrators.o: integrators.hpp integrators.cpp
	g++ $(CPPFLAGS) -c integrators.cpp

gqconstants: gqconstants.cpp integrators.hpp integrators.o
	g++ $(CPPFLAGS) -o gqconstants gqconstants.cpp integrators.o

clean: 
	rm -f *.o gqconstants *.png *.ps *.pdf *~ *.d *.pcm

