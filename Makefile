CC=icpc -std=c++11

HOME=/home/old_home/haoz

INCLUDE=${HOME}/tools/htslib/include

LIBS=${HOME}/tools/htslib/lib
#LIBS=/home/haoz/tools/htslib/lib/libhts.a

#FLAGS= -std=c++11 -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib -O3
FLAGS= -std=c++11 -lhts -O3 -g -qopenmp -ffast-math

#$(CC) parseCigar.o $(FLAGS) -I$(INCLUDE) -L$(LIBS)   -o parseCigar 

OBJS= Launcher.o RegionBuilder.o cigarModifier.o parseCigar.o recordPreprocessor.o VariationRealigner.o ToVarsBuilder.o  simpleMode.o somaticMode.o ampliconMode.o

launcher: $(OBJS)
	$(CC) -o launcher $(OBJS) $(FLAGS) -I$(INCLUDE) -L$(LIBS) 

recordPreprocessor.o:recordPreprocessor.cpp
	$(CC) -c recordPreprocessor.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

cigarModifier.o:cigarModifier.cpp
	$(CC) -c cigarModifier.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

parseCigar.o:parseCigar.cpp
	$(CC) -c parseCigar.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

RegionBuilder.o: RegionBuilder.cpp
	$(CC) -c RegionBuilder.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

VariationRealigner.o: VariationRealigner.cpp
	$(CC) -c VariationRealigner.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

ToVarsBuilder.o: ToVarsBuilder.cpp
	$(CC) -c ToVarsBuilder.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

Launcher.o: Launcher.cpp
	$(CC) -c Launcher.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

simpleMode.o: ./modes/simpleMode.cpp
	$(CC) -c ./modes/simpleMode.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

somaticMode.o: ./modes/somaticMode.cpp
	$(CC) -c ./modes/somaticMode.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

ampliconMode.o: ./modes/ampliconMode.cpp
	$(CC) -c ./modes/ampliconMode.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS)

.PHONY:clean
clean:
	rm *.o


#g++ -std=c++11 recordPreprocessor.cpp parseCigar.cpp -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib -O3



