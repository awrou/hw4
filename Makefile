OBJ_1 = heat_serial.o
OBJ_2 = heat_omp.o
OBJ_3 = heat_mpi.o
CC=g++
CXXFLAGS = -g -Wall -fopenmp

all: heat_serial heat_omp heat_mpi


heat_serial : $(OBJ_1)
	$(CXX) -o $@ $^

heat_omp : $(OBJ_2)
	$(CXX) -o $@ $^ -fopenmp

heat_mpi : $(OBJ_3)
	mpicxx -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
