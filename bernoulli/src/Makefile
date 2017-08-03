GCC= g++
OMPFLAGS = -g -fopenmp

all: Adv4DScanBer

Adv4DScanBer.o: Adv4DScanBer.c
	$(GCC) -o $@ -c $<
io.o: io.c
	$(GCC) -o $@ -c $<
scan.o: scan.c
	$(GCC) $(OMPFLAGS) -o $@ -c $<
Adv4DScanBer: Adv4DScanBer.o io.o scan.o
	$(GCC) $(OMPFLAGS) -o ../$@ $+
clean:
	rm -f *.o ../Adv4DScanBer
