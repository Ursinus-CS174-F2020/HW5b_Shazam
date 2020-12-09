#include "dsp.h"
#include "dspviz.h"
#include "audio.h"
#include "fingerprint.h"


int main(int argc, char** argv) {
    SongData s(argc, argv);
    s.processAudio();
    return 0;
}