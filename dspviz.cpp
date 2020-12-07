#include "dspviz.h"

/**
 * @param S The spectrogram
 * @param fftlen Number of samples in fft
 * @param NWin Number of windows
 * @param hScale: Factor by which to shrink height
 * @param wScale: Factor by which to shrink width
 * @return Canvas with fft image
 */
BMP plotSpectrogram(double** S, int fftlen, int NWin, int hScale, int wScale) {
    int H = fftlen/hScale;
    int W = NWin/wScale;
    BMP canvas(W, H);
    double min = S[0][0];
    double max = S[0][0];
    for (int i = 0; i < fftlen; i++) {
        for (int j = 0; j < NWin; j++) {
            double s = S[i][j];
            if (s < min) {
                min = s;
            }
            if (s > max) {
                max = s;
            }
        }
    }
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            double avg = 0.0;
            for (int di = 0; di < hScale; di++) {
                for (int dj = 0; dj < wScale; dj++) {
                    avg += S[i*hScale+di][j*wScale+dj];
                }
            }
            avg /= hScale*wScale;
            uint8_t val = (uint8_t)(255*avg/(max-min));
            canvas.set_pixel(j, i, val, val, val, 0xFF);
        }
    }
    return canvas;
}