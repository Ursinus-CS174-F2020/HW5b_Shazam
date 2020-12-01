/**
 * Programmer: Chris Tralie
 * Purpose: Vanilla C++ implementation of the Fast-Fourier Transform
 * and Short-Time Fourier Transform
 */

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
		 * @param start Index to start in the array
		 * @param len Length to go in the array
		 * @param useWindow Whether to use the window
		 */
		cdouble* toWindowedComplexArray(short* data, int start, int len, bool useWindow);
};

