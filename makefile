CC=g++
CFLAGS=-std=c++11 -g -Wall 

all: 

dsp.o: dsp.cpp dsp.h
	$(CC) $(CFLAGS) -c dsp.cpp

audio.o: audio.cpp audio.h
	$(CC) $(CFLAGS) -c audio.cpp

hashable.o: hashable.cpp hashable.h cloneable.h
	$(CC) $(CFLAGS) -c hashable.cpp 
	
hashtable.o: hashtable.h hashtable.cpp hashable.o map.h
	$(CC) $(CFLAGS) -c hashtable.cpp 

clean:
	rm *.o *.exe *.stackdump