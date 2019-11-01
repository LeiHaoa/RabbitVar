TARGET=sam
CC=g++ -std=c++11

INCLUDE=/home/haoz/tools/htslib/include

LIBS=/home/haoz/tools/htslib/lib
#LIBS=/home/haoz/tools/htslib/lib/libhts.a

#FLAGS= -std=c++11 -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib -O3
FLAGS= -std=c++11 -lhts -O3 -g

#$(CC) parseCigar.o $(FLAGS) -I$(INCLUDE) -L$(LIBS)   -o parseCigar 

launcher:Launcher.o RegionBuilder.o parseCigar.o recordPreprocessor.o VariationRealigner.o
	$(CC) -o launcher Launcher.o RegionBuilder.o parseCigar.o recordPreprocessor.o VariationRealigner.o $(FLAGS) -I$(INCLUDE) -L$(LIBS) 

recordPreprocessor.o:recordPreprocessor.cpp
	$(CC) -c recordPreprocessor.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

parseCigar.o:parseCigar.cpp
	$(CC) -c parseCigar.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

##launcher:Launcher.o RegionBuilder.o
##	$(CC) -o launcher Launcher.o RegionBuilder.o  $(FLAGS) -I$(INCLUDE) -L$(LIBS) 

RegionBuilder.o: RegionBuilder.cpp
	$(CC) -c RegionBuilder.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

VariationRealigner.o: VariationRealigner.cpp
	$(CC) -c VariationRealigner.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

Launcher.o: Launcher.cpp
	$(CC) -c Launcher.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)


.PHONY:clean
clean:
	rm *.o


#g++ -std=c++11 recordPreprocessor.cpp parseCigar.cpp -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib -O3



