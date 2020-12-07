#include "dspviz.h"

/**
 * @param S The spectrogram
 * @param maxBin: Max frequency up to which to plot results
 * @param NWin Number of windows
 * @param hScale: Factor by which to scale height
 * @param wScale: Factor by which to scale width
 * @return Canvas with fft image
 */
BMP plotSpectrogram(double** S, int maxBin, int NWin, int hScale, int wScale) {
    int H = maxBin*hScale;
    int W = NWin*wScale;
    BMP canvas(W, H);
    double min = S[0][0];
    double max = S[0][0];
    for (int i = 0; i < maxBin; i++) {
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
    for (int i = 0; i < maxBin; i++) {
        for (int j = 0; j < NWin; j++) {
            uint8_t val = (uint8_t)(255*S[i][j]/(max-min));
            for (int di = 0; di < hScale; di++) {
                for (int dj = 0; dj < wScale; dj++) {
                    canvas.set_pixel(j*wScale+dj, i*hScale+di, val, val, val, 0xFF);
                }
            }
        }
    }
    return canvas;
}