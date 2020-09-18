CXXFLAGS=-Wall -pg
OBJS = parameter.o compose.o nrutils.o gasdev.o gamdev.o ran1.o gammln.o bnldev.o motif.o stdp_noisy.o

main: $(OBJS)
	g++ -pg -o stdp_noisy $(OBJS)
	cp stdp_noisy output/

.PHONY: clean
clean:
	rm -f *.o
	rm -f stdp_noisy	

.PHONY: data
data: 
	rm -f *.dat
