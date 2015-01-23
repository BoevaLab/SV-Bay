gccount: main.o ConfigFile.o Chameleon.o myFunc.o Genome.o Chromosome.o GenomeElement.o
	g++ -O3 -Wall -m64 -o tgsim main.o ConfigFile.o Chameleon.o myFunc.o Genome.o Chromosome.o GenomeElement.o
	rm main.o
	rm ConfigFile.o
	rm Chameleon.o	
	rm myFunc.o	
	rm Chromosome.o 
	rm GenomeElement.o
	rm Genome.o		
main.o: main.cpp tgsim.h
	g++ -g -O3 -c -Wall -m64 main.cpp
ConfigFile.o: ConfigFile.cpp ConfigFile.h
	g++ -g -O3 -c -Wall -m64  ConfigFile.cpp
Chameleon.o: Chameleon.cpp Chameleon.h
	g++ -g -O3 -c -Wall -m64  Chameleon.cpp
myFunc.o: myFunc.cpp myFunc.h
	g++ -g -c -O3 -Wall -m64  myFunc.cpp	
Genome.o: Genome.cpp Genome.h
	g++ -g -c -O3 -Wall -m64  Genome.cpp	
ConnectedElements.o: Chromosome.cpp Chromosome.h
	g++ -g -c -O3 -Wall -m64  Chromosome.cpp	
GenomeElement.o: GenomeElement.cpp GenomeElement.h
	g++ -g -c -O3 -Wall -m64  GenomeElement.cpp	
clean:
	rm -f gccount *.o *~ *#          
