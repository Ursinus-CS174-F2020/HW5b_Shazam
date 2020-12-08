#include "dsp.h"
#include "dspviz.h"
#include "audio.h"

int main() {
    size_t N;
    short* audio = getAudio("femalecountdown.wav", &N);
    int win = 2048;
    int swin = win/2+1;
    int hop = 512;
    DSP dsp(win);
    int NWin = 0;
    double** specgram = dsp.specgram(audio, N, win, hop, false, &NWin);

    BMP canvas = plotSpectrogram(specgram, swin, NWin, 1, 1);
    canvas.write("S.bmp");

    deleteSpecgram(specgram, win);
    delete[] audio;
    return 0;
}