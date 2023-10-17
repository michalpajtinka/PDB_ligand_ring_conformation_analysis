PROGRAM=ConfAnalyser
SOURCES=main.cpp application.cpp point_3D.cpp vector_3D.cpp plane_3D.cpp atom.cpp \
		angle.cpp molecule.cpp ring.cpp six_atom_ring.cpp five_atom_ring.cpp \
		benzene.cpp cyclohexane.cpp cyclopentane.cpp pyrane.cpp helper_functions.cpp

CXX=g++
CXXFLAGS=-Wall -Wextra -ansi -pedantic -O3 -std=c++20
OBJS=$(SOURCES:.cpp=.o)
RM=rm -f

all:$(PROGRAM)

$(PROGRAM):$(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	$(RM) $(OBJS)

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	$(RM) $(PROGRAM) *.o

debug: CXXFLAGS+=-O0 -g
debug: all
