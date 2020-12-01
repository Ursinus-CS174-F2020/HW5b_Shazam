CC=g++
CFLAGS=-std=c++11 -g -Wall 

all: 

hashable.o: hashable.cpp hashable.h cloneable.h
	$(CC) $(CFLAGS) -c hashable.cpp 
	
hashtable.o: hashtable.h hashtable.cpp hashable.o map.h
	$(CC) $(CFLAGS) -c hashtable.cpp 

clean:
	rm *.o *.exe *.stackdump