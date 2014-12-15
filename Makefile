#FLAGS = -I./ -O3 -pg
#FLAGS = -I./ -O3 -g
FLAGS = -I ./ -O3 -g -DHAVE_CONFIG_H=1 
SRC = gclust.cpp paraSA.cpp fasta.cpp divsufsort.c sssort.c trsort.c

all: gclust 

gclust: gclust.o paraSA.o fasta.o divsufsort.o sssort.o trsort.o
	g++ -lpthread $(FLAGS) $^ -o $@

.cpp.o:
	g++ $(FLAGS) -Wall -c $<

.c.o:
	gcc $(FLAGS) -Wall -c $<

clean: 
	rm -f *.o gclust

