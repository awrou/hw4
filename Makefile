OBJ_1 = heat_serial.o

CXXFLAGS = -g -Wall

all: heat_serial

heat_serial : $(OBJ_1)
	$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
