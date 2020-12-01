/**
 * Programmer: Chris Tralie
 * Purpose: Vanilla C++ implementation of the Fast-Fourier Transform
 * and Short-Time Fourier Transform
 */

#include <complex>
#include <iostream>
#include <math.h>

using namespace std;


#define FFT_FORWARD 0
#define FFT_INVERSE 1
typedef complex<double> cdouble;

/**
 * Compute the closest power of 2 greater than or equal
 * to some number
 * @param a The number
 * @return The closest power of 2 >= a
 */
int getClosestPowerOf2(int a) {
	double lg = log((double)a) / log(2.0);
	int power = (int)lg;
	if ((double)((int)lg) < lg) {
		power++;
	}
	return power;
}

/**
 * Compute a version of x which is bit-reversed
 * @param x A 32-bit int to reverse
 * @param length Length of bits
 * @return Bit-reversed version of x
 */
int bitReverse(int x, int length) {
	int toReturn = 0;
	int mirror = length / 2;
	for (int mask = 1; mask <= length; mask <<= 1, mirror >>= 1) {
		if ((mask & x) > 0)
			toReturn |= mirror;
	}
	return toReturn;
}


/**
 * Rearrange the terms in-place so that they're sorted by the least
 * significant bit (this is the order in which the terms are accessed
 * in the FFT)
 * @param a An array of complex numbers
 * @param N Number of elements in the array
 */
void rearrange(cdouble* a, int N) {
	for (int i = 0; i < N; i++) {
		int j = bitReverse(i, N);
		if (j > i) { //Don't waste effort swapping two mirrored
		//elements that have already been swapped
			cdouble temp = a[j];
			a[j] = a[i];
			a[i] = temp;
		}
	}
}

class DSP {
	private:
		cdouble*** W;//Cache the complex coefficients for the FFT
		double* hannwindow;
		int fftsize;

		/**
		 * Initialize the complex coefficients for the FFT
		 * @param fftsize The length of the fft
		 */
		void initCoeffs(int fftsize) {
			int maxlevel = getClosestPowerOf2(fftsize) + 1;
			W = new cdouble**[maxlevel+1];
			for (int level = 1; level <= maxlevel; level++) {
				int FFTSize = 1 << level;
				W[level] = new cdouble*[2];
				W[level][0] = new cdouble[FFTSize >> 1];
				W[level][1] = new cdouble[FFTSize >> 1];
				for (int i = 0; i < FFTSize >> 1; i++) {
					double iangle = (double)i * 2.0 * M_PI / (double)FFTSize;
					double fangle = (-1.0) * iangle;
					W[level][FFT_FORWARD][i] = cdouble(cos(fangle), sin(fangle));
					W[level][FFT_INVERSE][i] = cdouble(cos(iangle), sin(iangle)); 
				}
			}
		}

		/**
		 * Initialize the coefficients in Hann window
		 * @param fftsize Lenght of fft
		 */
		void initWindow(int fftsize) {
			int N = fftsize;
			hannwindow = new double[N];
			for (int n = 0; n < N; n++) {
				double angle = 2.0*M_PI * n / (double)(N - 1);
				//Do a Hamming hannwindow for now
				hannwindow[n] = 0.54 - 0.46*cos(angle);
			}
		}

		/**
		 * Perform an in-place Cooley-Tukey FFT
		 * @param toReturn Array that holds FFT coefficients
		 * @param N Length of array (assumed to be power of 2)
		 * @param inverse Whether this is a forward or inverse FFT
		 */
		cdouble* performfft(cdouble* toReturn, int N, int inverse) {
			rearrange(toReturn, N);
			//Do the trivial FFT size of 2 first
			for (int i = 0; i < N; i += 2) {
				cdouble temp = toReturn[i];
				toReturn[i] = temp + toReturn[i + 1];
				toReturn[i + 1] = temp - toReturn[i + 1];
			}
			int Mindex = 2;//Index used to access the cached complex
			//coefficients
			for (int level = 2; level < N; level <<= 1) {
				double angle = 0.0;
				int FFTSize = level << 1;
				for (int start = 0; start < N; start += FFTSize) {
					//This is a little chunk of an FFT of size "FFTSize"
					//to do in-place with the merging algorithm
					//NOTE: "level" gives the length between mirrored terms
					for (int i = 0; i < level; i++) {
						cdouble coeff = W[Mindex][inverse][i];
						cdouble first = toReturn[start + i];
						cdouble second = coeff*toReturn[start + i + level];
						toReturn[start + i] = first + second;
						toReturn[start + i + level] = first - second;
					}
				}
				Mindex++;
			}
			return toReturn;
		}
	
	public:
		cdouble fftres;
		DSP(int fftsize) {
			this->fftsize = fftsize;
			this->initCoeffs(fftsize);
			this->initWindow(fftsize);
		}
		~DSP() {
			int maxlevel = getClosestPowerOf2(fftsize) + 1;
			// Clean up FFT coefficients
			for (int level = 1; level <= maxlevel; level++) {
				for (int type = 0; type < 2; type++) {
					delete[] W[level][type];
				}
				delete[] W[level];
			}
			delete[] W;
			// Clean up window coefficients
			delete[] hannwindow;
		}

		/**
		 * Implement the dft directly from the definition (used for speed comparison)
		 * @param sig Complex signal on which to compute dft
		 * @param N Length of signal
		 * @return Complex DFT coefficients
		 */
		cdouble* dft(cdouble* sig, int N) {
			cdouble* toReturn = new cdouble[N];
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					double angle = -2.0 * M_PI * (double)i * (double)j / (double)N;
					cdouble coeff(cos(angle), sin(angle));
					toReturn[i] = coeff*sig[i];
				}
			}
			return toReturn;
		}

		/**
		 * Perform the FFT on a complex signal
		 * @param sig The signal
		 * @param N Length of the signal (assumed to be power of 2)
		 * @return An N-length array with FFT coefficients
		 */
		cdouble* fft(cdouble* sig, int N) {
			cdouble* toReturn = new cdouble[N];
			for (int i = 0; i < N; i++) {
				toReturn[i] = sig[i];
			}
			return performfft(toReturn, N, FFT_FORWARD);	
		}
	
		/**
		 * Perform the inverse FFT on an array of complex FFT coefficients
		 * @param sig The FFT coefficients
		 * @param N Length of the FFT coefficients (assumed to be power of 2)
		 * @return An N-length array with FFT coefficients
		 */
		cdouble* ifft(cdouble* sig, int N) {
			cdouble* toReturn = new cdouble[N];
			for (int i = 0; i < N; i++) {
				toReturn[i] = sig[i];
				//Scale by 1/N for inverse FFT
				toReturn[i] *= cdouble(1.0/(double)N, 0);
			}
			return performfft(toReturn, N, FFT_INVERSE);
		}
		
		/**
		 * Helper function to create a complex array out of an array of 
		 * real amplitude samples
		 * @param data An array of shorts for the audio data
		 * @param start Index to start in the array
		 * @param len Length to go in the array
		 * @param useWindow Whether to use the window
		 */
		cdouble* toWindowedComplexArray(short* data, int start, int len, bool useWindow) {
			int N = 1 << getClosestPowerOf2(len);
			cdouble* toReturn = new cdouble[N];
			//Make a complex array out of the real array
			for (int i = 0; i < N; i++) {
				if (i < len) {
					short value = data[start + i];
					toReturn[i] = cdouble((double)value, 0.0);
					if (useWindow) {
						toReturn[i] *= hannwindow[i];
					}
				}
				else {
					//Zero pad if not a power of 2
					toReturn[i] = cdouble(0.0, 0.0);
				}
			}
			return toReturn;
		}
};
