CXXFLAGS=-Wall -pg
OBJS = poissonPairs_timeAboveThreshold.o

main: $(OBJS)
	g++ $(OBJS) -o poissonPairs_timeAboveThreshold

.PHONY: clean
clean:
	rm -f *.o
	rm -f poissonPairs_timeAboveThreshold

.PHONY: data
data: 
	rm -f *.dat
