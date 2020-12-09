#include "dsp.h"
#include "dspviz.h"
#include "audio.h"
#include "fingerprint.h"
#include "time.h"


int main(int argc, char** argv) {
    time_t tic, toc;

    // Step 1: Load in and process query song
    time(&tic);
    SongData s(argc, argv);
    s.processAudio();
    time(&toc);
    double seconds = difftime(toc, tic);
    printf("%.g seconds elapsed processing query audio\n", seconds);

    // Step 2: Load in database
    int NSongs = 6;
    songInfo* songs;
    time(&tic);
    HashTable* db = loadDatabase(NSongs, 1000, &songs);
    time(&toc);
    seconds = difftime(toc, tic);
    printf("%.g seconds elapsed loading database\n", seconds);

    // Step 3: Load in query song
    time(&tic);
    s.queryDatabase(db, songs, NSongs).printInfo();
    time(&toc);
    seconds = difftime(toc, tic);
    printf("%.g seconds elapsed performing query\n", seconds);
    delete db;
    delete[] songs;
    return 0;
}