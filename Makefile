galsim: quadtree.o simulations.o main.o particles.o
	gcc -Wall main.o simulations.o quadtree.o particles.o  -g -o galsim  -lm

quadtree.o: quadtree.c quadtree.h particles.h
	gcc -Wall -c quadtree.c

simulations.o: simulations.c simulations.h quadtree.h particles.h
	gcc -Wall -c simulations.c

particles.o: particles.h quadtree.h
	gcc -Wall -c particles.c

main.o:
	gcc -Wall -c main.c

clean:
	rm -f galsim simulations.o main.o file_handler.o quadtree.o particles.o
