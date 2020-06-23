CC=icpc -std=c++11

HOME=/home/old_home/haoz

INCLUDE=${HOME}/tools/htslib/include

LIBS=${HOME}/tools/htslib/lib
#LIBS=/home/haoz/tools/htslib/lib/libhts.a

#FLAGS= -std=c++11 -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib -O3
FLAGS= -std=c++11 -lhts -O3 -g -qopenmp -ffast-math

#$(CC) parseCigar.o $(FLAGS) -I$(INCLUDE) -L$(LIBS)   -o parseCigar 

OBJS= ./objs/Launcher.o ./objs/RegionBuilder.o ./objs/cigarModifier.o ./objs/parseCigar.o ./objs/recordPreprocessor.o ./objs/VariationRealigner.o ./objs/ToVarsBuilder.o  ./objs/simpleMode.o ./objs/somaticMode.o ./objs/ampliconMode.o

FastVC: $(OBJS)
	$(CC) -o FastVC $(OBJS) $(FLAGS) -I$(INCLUDE) -L$(LIBS) 

./objs/recordPreprocessor.o:./src/recordPreprocessor.cpp
	$(CC) -c ./src/recordPreprocessor.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/recordPreprocessor.o

./objs/cigarModifier.o:./src/cigarModifier.cpp
	$(CC) -c ./src/cigarModifier.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/cigarModifier.o

./objs/parseCigar.o:./src/parseCigar.cpp
	$(CC) -c ./src/parseCigar.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/parseCigar.o

./objs/RegionBuilder.o: ./src/RegionBuilder.cpp
	$(CC) -c ./src/RegionBuilder.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/RegionBuilder.o

./objs/VariationRealigner.o: ./src/VariationRealigner.cpp
	$(CC) -c ./src/VariationRealigner.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/VariationRealigner.o

./objs/ToVarsBuilder.o: ./src/ToVarsBuilder.cpp
	$(CC) -c ./src/ToVarsBuilder.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/ToVarsBuilder.o

./objs/Launcher.o: ./src/Launcher.cpp
	$(CC) -c ./src/Launcher.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/Launcher.o

./objs/simpleMode.o: ./src/modes/simpleMode.cpp
	$(CC) -c ./src/modes/simpleMode.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/simpleMode.o

./objs/somaticMode.o: ./src/modes/somaticMode.cpp
	$(CC) -c ./src/modes/somaticMode.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/somaticMode.o

./objs/ampliconMode.o: ./src/modes/ampliconMode.cpp
	$(CC) -c ./src/modes/ampliconMode.cpp $(FLAGS) -I$(INCLUDE) -L$(LIBS) -o ./objs/ampliconMode.o

.PHONY:clean
clean:
	rm ./objs/*.o


#g++ -std=c++11 recordPreprocessor.cpp parseCigar.cpp -lhts -I/home/haoz/tools/htslib/include -L/home/haoz/tools/htslib/lib -O3



