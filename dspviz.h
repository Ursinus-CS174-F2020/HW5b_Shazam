#ifndef DSPVIZ_H
#define DSPVIZ_H

#include "dsp.h"
#include "BMP.h"

/**
 * @param S The spectrogram
 * @param fftlen Number of samples in fft
 * @param NWin Number of windows
 * @param hScale: Factor by which to shrink height
 * @param wScale: Factor by which to shrink width
 * @return Canvas with fft image
 */
BMP plotSpectrogram(double** S, int fftlen, int NWin, int hScale, int wScale);

#endif