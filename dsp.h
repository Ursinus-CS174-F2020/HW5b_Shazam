/**
 * Programmer: Chris Tralie
 * Purpose: Vanilla C++ implementation of the Fast-Fourier Transform
 * and Short-Time Fourier Transform
 */
#ifndef DSP_H
#define DSP_H

#include <complex>
#include <iostream>
#include <math.h>

using namespace std;

typedef complex<double> cdouble;
#define FFT_FORWARD 0
#define FFT_INVERSE 1

class DSP {
	private:
		cdouble*** W;//Cache the complex coefficients for the FFT
		double* hannwindow;
		int fftsize;

		/**
		 * Initialize the complex coefficients for the FFT
		 * @param fftsize The length of the fft
		 */
		void initCoeffs(int fftsize);

		/**
		 * Initialize the coefficients in Hann window
		 * @param fftsize Lenght of fft
		 */
		void initWindow(int fftsize);

		/**
		 * Perform an in-place Cooley-Tukey FFT
		 * @param toReturn Array that holds FFT coefficients
		 * @param N Length of array (assumed to be power of 2)
		 * @param inverse Whether this is a forward or inverse FFT
		 */
		cdouble* performfft(cdouble* toReturn, int N, int inverse);
	
	public:
		cdouble fftres;
		DSP(int fftsize);
		~DSP();

		/**
		 * Implement the dft directly from the definition (used for speed comparison)
		 * @param sig Complex signal on which to compute dft
		 * @param N Length of signal
		 * @return Complex DFT coefficients
		 */
		cdouble* dft(cdouble* sig, int N);

		/**
		 * Perform the FFT on a complex signal
		 * @param sig The signal
		 * @param N Length of the signal (assumed to be power of 2)
		 * @return An N-length array with FFT coefficients
		 */
		cdouble* fft(cdouble* sig, int N);
	
		/**
		 * Perform the inverse FFT on an array of complex FFT coefficients
		 * @param sig The FFT coefficients
		 * @param N Length of the FFT coefficients (assumed to be power of 2)
		 * @return An N-length array with FFT coefficients
		 */
		cdouble* ifft(cdouble* sig, int N);
		
		/**
		 * Helper function to create a complex array out of an array of 
		 * real amplitude samples
		 * @param data An array of shorts for the audio data
		 * @param N Total number of samples in data
		 * @param start Index to start in the array
		 * @param len Length to go in the array
		 * @param useWindow Whether to use the window
		 */
		cdouble* toWindowedComplexArray(short* data, int N, int start, int len, bool useWindow);

		/**
		 * Perform a short-time fourier transform on a bunch of samples
		 * @param sig Samples in the signal
		 * @param N Length of signal
		 * @param win Window length
		 * @param hop Hop length
		 * @param useWindow Whether to use the window
		 * @param NWin Number of windows (returned by reference)
		 */
		cdouble** stft(short* sig, int N, int win, int hop, bool useWindow, int* NWin);

		/**
		 * Perform a magnitude short-time fourier transform on a bunch of samples
		 * @param sig Samples in the signal
		 * @param N Length of signal
		 * @param win Window length
		 * @param hop Hop length
		 * @param useWindow Whether to use the window
		 * @param NWin Number of windows (returned by reference)
		 */
		double** specgram(short* sig, int N, int win, int hop, bool useWindow, int* NWin);
};

/**
 * Free the memory associated to an STFT
 * @param S STFT
 * @param win Window length
 */
void deleteSTFT(cdouble** S, int win);

/**
 * Free the memory associated to a spectrogram
 * @param S Spectrogram
 * @param win Window length
 */
void deleteSpecgram(double** S, int win);

#endif