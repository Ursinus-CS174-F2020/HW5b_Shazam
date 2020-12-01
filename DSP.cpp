/**
 * Programmer: Chris Tralie
 * Purpose: Vanilla C++ implementation of the Fast-Fourier Transform
 * and Short-Time Fourier Transform
 */

#include <complex>
#include <math.h>

using namespace std;

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
		const int FORWARD = 0;
		const int INVERSE = 1;
		cdouble*** W;//Cache the complex coefficients for 
		cdouble* hannwindow;
		int fftsize;

		/**
		 * Initialize the complex coefficients for the FFT
		 * @param fftsize The length of the fft
		 */
		void initCoeffs(int fftsize) {
			int maxlevel = getClosestPowerOf2(fftsize) + 1;
			W = new Complex[maxlevel + 1][2][];
			W = new cdouble**[maxlevel+1];
			for (int level = 1; level <= maxlevel; level++) {
				int FFTSize = 1 << level;
				W[level] = new cdouble*[2];
				W[level][0] = new cdouble[FFTSize >> 1];
				W[level][1] = new cdouble[FFTSize >> 1];
				for (int i = 0; i < FFTSize >> 1; i++) {
					double iangle = (double)i * 2.0 * M_PI / (double)FFTSize;
					double fangle = (-1.0) * iangle;
					W[level][FORWARD][i] = cdouble(cos(fangle), sin(fangle));
					W[level][INVERSE][i] = cdouble(cos(iangle), sin(iangle)); 
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
				hannwindow[n] = 0.54 - 0.46*Math.cos(angle);
			}
		}

		/**
		 * Perform an in-place Cooley-Tukey FFT
		 * @param toReturn Array that holds FFT coefficients
		 * @param N Length of array
		 * @param inverse Whether this is a forward or inverse FFT
		 */
		cdouble* performfft(cdouble* toReturn, int N, int inverse) {
			rearrange(toReturn, N);
			//Do the trivial FFT size of 2 first
			for (int i = 0; i < N; i += 2) {
				Complex temp = toReturn[i];
				toReturn[i] = temp.add(toReturn[i + 1]);
				toReturn[i + 1] = temp.subtract(toReturn[i + 1]);
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
						Complex coeff = W[Mindex][inverse][i];
						Complex first = toReturn[start + i];
						Complex second = toReturn[start + i + level].mult(coeff);
						toReturn[start + i] = first.add(second);
						toReturn[start + i + level] = first.subtract(second);
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
			this->initComplex(fftsize);
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
					Complex coeff = new cdouble(cos(angle), sin(angle));
					toReturn[i] = sig[i].mult(coeff);
				}
			}
			return toReturn;
		}
	
	public Complex[] fft(Complex[] sig) {
		Complex[] toReturn = new Complex[sig.length];
		for (int i = 0; i < toReturn.length; i++)
			toReturn[i] = sig[i].clone();
		return performfft(toReturn, FORWARD);		
	}
	
	public Complex[] ifft(Complex[] sig) {
		Complex[] toReturn = new Complex[sig.length];
		for (int i = 0; i < toReturn.length; i++) {
			if (sig[i] != null)
				toReturn[i] = sig[i].clone();
			else
				toReturn[i] = new Complex(0, 0);
			//Scale by 1/N for inverse FFT
			toReturn[i].scale(1.0 / (double)sig.length);
		}
		return performfft(toReturn, INVERSE);
	}

	
	//Helper function to create a complex array out of an array of 
	//real amplitude samples
	public Complex[] toWindowedComplexArray(short[] data, int start, int len) {
		int N = 1 << getClosestPowerOf2(len);
		Complex[] toReturn = new Complex[N];
		//Make a complex array out of the real array
		for (int i = 0; i < N; i++) {
			if (i < len) {
				short value = data[start + i];
				toReturn[i] = new Complex((double)value, 0.0);
				//toReturn[i].scale(hannwindow[i]);
			}
			else
				//Zero pad if not a power of 2 (this shouldn't happen)
				toReturn[i] = new Complex(0.0, 0.0);
		}
		return toReturn;
	}
	
	//Return the strongest bin
	public int getStrongestBin(Complex[] spec, int minbin) {
		double maxPower = 0.0;
		int maxbin = 0;
		//Only go up to the Nyquist Rate (spec.length / 2)
		for (int i = minbin; i < spec.length / 2; i++) {
			double power = spec[i].magSquared();
			if (power > maxPower) {
				maxPower = power;
				maxbin = i;
			}
		}
		return maxbin;
	}
	
	public double getMean(byte[] data) {
		double mean = 0.0;
		for (int i = 0; i < data.length; i++) {
			mean += (short)data[i];
		}
		mean /= (double)data.length;
		return mean;
	}
	
	public double getVariance(byte[] data) {
		double mean = getMean(data);
		double var = 0.0;
		for (int i = 0; i < data.length; i++) {
			double diff = (short)data[i] - mean;
			var += diff*diff;
		}
		return Math.sqrt(var);
	}
	
	public double getStrongestFreq(byte[] data, double sampleRate) {
		short[] sdata = new short[data.length];
		for (int i = 0; i < sdata.length; i++)
			sdata[i] = (short)data[i];
		Complex[] carray = toWindowedComplexArray(sdata, 0, sdata.length);
		Complex[] spec = fft(carray);
		int minbin = (int)((220.0 / 8000.0) * (double)data.length);
		int strongestBin = getStrongestBin(spec, minbin);
		fftres = spec;
		return ((double)strongestBin / (double)data.length)*sampleRate;
	}
	
	public int getMaxBin(byte[] data, double sampleRate, double centerFreq) {
		short[] sdata = new short[data.length];
		for (int i = 0; i < sdata.length; i++)
			sdata[i] = (short)data[i];
		Complex[] carray = toWindowedComplexArray(sdata, 0, sdata.length);
		Complex[] spec = fft(carray);
		fftres = spec;
		
		int bin = -2;
		double maxSpec = 0.0;
		int maxbin = -2;
		//-12 to -7
		//-7 to -2
		//-2 to 3
		//3 to 8
		//8 to 13
		for (int halfstep = -12; halfstep < 12; halfstep += 5) {
			double freqLeft = centerFreq * Math.pow(2.0, (double)halfstep / 12.0);
			double freqRight = centerFreq * Math.pow(2.0, (double)(halfstep+5) / 12.0);
			int binLeft = (int)((freqLeft / sampleRate)*data.length);
			int binRight = (int)((freqRight / sampleRate)*data.length);
			//Now integrate the spectrum over that interval and average
			double specPow = 0.0;
			for (int i = binLeft; i < binRight; i++) {
				specPow += spec[i].magSquared();
			}
			specPow /= (binRight - binLeft);
			if (specPow > maxSpec) {
				maxSpec = specPow;
				maxbin = bin;
			}
			bin++;
		}
		return maxbin;
	}
	
	public static void main(String[] args) {
		double centerFreq = 785.0;
		double sampleRate = 8000;
		short[] data = new short[256];
		for (int halfstep = -12; halfstep < 12; halfstep += 5) {
			double freqLeft = centerFreq * Math.pow(2.0, (double)halfstep / 12.0);
			double freqRight = centerFreq * Math.pow(2.0, (double)(halfstep+5) / 12.0);
			int binLeft = (int)((freqLeft / sampleRate)*data.length);
			int binRight = (int)((freqRight / sampleRate)*data.length);
			System.out.println(freqLeft + " - " + freqRight + " : " + binLeft + " - " + binRight);
		}
	}
	
}