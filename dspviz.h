#ifndef DSPVIZ_H
#define DSPVIZ_H

#include "dsp.h"
#include "BMP.h"

/**
 * @param S The spectrogram
 * @param maxBin: Max frequency up to which to plot results
 * @param fftlen Number of samples in fft
 * @param NWin Number of windows
 * @param hScale: Factor by which to scale height
 * @param wScale: Factor by which to scale width
 * @return Canvas with fft image
 */
BMP plotSpectrogram(double** S, int maxBin, int NWin, int hScale, int wScale);

#endif
