CFLAGS = -O3 -g 

CC = gcc
all : page mi_go_motif_calculator

page : page.c statistics.o dataio.o information.o mi_library.o hashtable.o
	$(CC) $(CFLAGS) -Wall -o page page.c statistics.o dataio.o information.o mi_library.o hashtable.o -lm $(LFLAGS)

mi_go_motif_calculator : mi_go_motif_calculator.c statistics.o dataio.o information.o mi_library.o hashtable.o
	$(CC) $(CFLAGS) -Wall -o mi_go_motif_calculator mi_go_motif_calculator.c statistics.o dataio.o information.o mi_library.o hashtable.o -lm $(LFLAGS)

mi_library.o : mi_library.c mi_library.h
	$(CC) $(CFLAGS) -Wall -c mi_library.c -static

statistics.o : statistics.c statistics.h
	$(CC) $(CFLAGS) -Wall -c statistics.c -static

dataio.o : dataio.c dataio.h
	$(CC) $(CFLAGS) -Wall -c dataio.c -static

information.o : information.c information.h
	$(CC) $(CFLAGS) -Wall -c information.c -static

hashtable.o : hashtable.c hashtable.h
	$(CC) $(CFLAGS) -Wall -c hashtable.c -static

clean: 
	rm *.o page mi_go_motif_calculator mi_go_motif_calculator.exe page.exe
