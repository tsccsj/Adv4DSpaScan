GCC= g++
OMPFLAGS = -g -fopenmp

all: Adv4DScanBer

Adv4DScanBer.o: Adv4DScanBer.c
	$(GCC) -o $@ -c $<
io.o: io.c
	$(GCC) -o $@ -c $<
scan.o: scan.c
	$(GCC) $(OMPFLAGS) -o $@ -c $<
mc.o: mc.c
	$(GCC) $(OMPFLAGS) -o $@ -c $< -std=c++11
Adv4DScanBer: Adv4DScanBer.o io.o scan.o mc.o
	$(GCC) $(OMPFLAGS) -o ../$@ $+
clean:
	rm -f *.o ../Adv4DScanBer
