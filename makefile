CC=g++
CFLAGS=-std=c++11 -g -Wall 

all: testspecgram

dsp.o: dsp.cpp dsp.h
	$(CC) $(CFLAGS) -c dsp.cpp

dspviz.o: dspviz.cpp dspviz.h BMP.h
	$(CC) $(CFLAGS) -c dspviz.cpp

audio.o: audio.cpp audio.h
	$(CC) $(CFLAGS) -c audio.cpp

hashable.o: hashable.cpp hashable.h cloneable.h
	$(CC) $(CFLAGS) -c hashable.cpp 
	
hashtable.o: hashtable.h hashtable.cpp hashable.o map.h
	$(CC) $(CFLAGS) -c hashtable.cpp 

testspecgram: testspecgram.cpp audio.o dsp.o dspviz.o
	$(CC) $(CFLAGS) -o testspecgram testspecgram.cpp audio.o dsp.o dspviz.o

clean:
	rm *.o *.exe *.stackdump