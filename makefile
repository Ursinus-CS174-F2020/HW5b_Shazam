CC=g++
CFLAGS=-std=c++11 -g -Wall 

all: loadFingerprints query

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

fingerprint.o: fingerprint.cpp fingerprint.h hashable.o hashtable.o audio.o dsp.o dspviz.o
	$(CC) $(CFLAGS) -c fingerprint.cpp 

loadFingerprints: loadFingerprints.cpp audio.o dsp.o dspviz.o hashable.o hashtable.o fingerprint.o
	$(CC) $(CFLAGS) -o loadFingerprints loadFingerprints.cpp audio.o dsp.o dspviz.o hashable.o hashtable.o fingerprint.o

query: query.cpp audio.o dsp.o dspviz.o hashable.o hashtable.o fingerprint.o
	$(CC) $(CFLAGS) -o query query.cpp audio.o dsp.o dspviz.o hashable.o hashtable.o fingerprint.o

clean:
	rm *.o *.exe *.stackdump loadFingerprints query