TARGET=sam
CC=g++ -std=c++11

INCLUDE=/home/haoz/tools/htslib/include

LIBS=/home/haoz/tools/htslib/lib
#LIBS=/home/haoz/tools/htslib/lib/libhts.a

FLAGS= -std=c++11 parseCigar.cpp -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib -O3

parseCigar:parseCigar.o preprocess.o
	#$(CC) parseCigar.o $(FLAGS) -I$(INCLUDE) -L$(LIBS)   -o parseCigar 
	$(CC) parseCigar.o preprocess.o $(FLAGS)  -o parseCigar

preprocess.o:recordPreprocessor.cpp
	$(CC)  recordPreprocessor.cpp $(FLAGS)  -o preprocess.o

parseCigar.o:parseCigar.cpp
	$(CC)  parseCigar.cpp $(FLAGS)  -o parseCigar.o


.PHONY:clean
clean:
	rm *.o


#g++ parseCigar.cpp -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib



