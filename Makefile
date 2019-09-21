TARGET=sam
CC=g++ -std=c++0x

INCLUDE=/home/haoz/tools/htslib/include

LIBS=/home/haoz/tools/htslib/lib
#LIBS=/home/haoz/tools/htslib/lib/libhts.a

FLAGS= -lhts

parseCigar:parseCigar.o
	$(CC) parseCigar.o $(FLAGS) -I$(INCLUDE) -L$(LIBS)   -o parseCigar 

sam.o:readsam.c
	$(CC) parseCigar.cpp -c $(FLAGS) -I$(INCLUDE) -L$(LIBS)   -o parseCigar.o 

.PHONY:clean
clean:
	rm *.o


#g++ parseCigar.cpp -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib



