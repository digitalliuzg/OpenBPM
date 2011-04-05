/** beatdetect.c
  * a lightweight beat detector
  *
  * see: http://developer.apple.com/library/mac/#documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html%23//apple_ref/doc/uid/TP40005147-CH202-SW1
  * see: http://developer.apple.com/library/mac/#documentation/Performance/Conceptual/vDSP_Programming_Guide/SampleCode/SampleCode.html%23//apple_ref/doc/uid/TP40005147-CH205-CIAEJIGF
  */

#include <stdio.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h>

#include "beatdetect.h"	// contains filters & example signals
//#include "song.h"	// song sample
//#include "filters.h"	// filters & signals

FFTSetup fft_weights;
DSPSplitComplex spectrum, zeros;

/*
int main()
{
	//float data[] = {1.001,0.744,1.23,4.23,0.99};
	float *data;
	//vDSP_Length indicies[] = {0,1,2,3,4};	// indicies corresponding in decimated-diff-rect'd vector
	vDSP_Length dataLength;
	
	data = sinusoid;
	dataLength = num_sinusoid;

//	float peaks[] = {0.0,1.0,2.0};	// malloc this instead?
	int beats[10];	// malloc this?
	int numBeats, i;
	
	//BeatDetect(data, dataLength, beats, &numBeats);
	BeatDetect(song, songLength, beats, &numBeats);
	
	for(i=0; i<numBeats; i++)
		printf("%d\n",beats[i]);
	
	return 0;
}
*/

/**
  * @brief - locates beats in signal
  */
void DoBeatDetect(float signal[], uint32_t signalLength, int *beats, int *num_beats)
{
	// 1: decimate
	// 2: decompose into subbands
	// 3: smooth (hanning window)
	// 4: differentiate / rectify
	// 5: comb
	
//	signal = noise;
//	signalLength = noiseLength;
	float* subband;


	int numbands = 6;
	uint32_t maxfreq = 4096;
	uint32_t bandlimits[] = {0,200,400,800,1600,3200};
	
	
	uint32_t log2n = (uint32_t)floorf(log2((double)signalLength));
	vDSP_Length n = (vDSP_Length) 1 << log2n;

	
	// allocate mem for the spectrum
	spectrum.realp = (float*)malloc(n/2*sizeof(float));
	spectrum.imagp = (float*)malloc(n/2*sizeof(float));

	zeros.realp = (float*)malloc(n/2*sizeof(float));
	zeros.imagp = (float*)malloc(n/2*sizeof(float));
	
	// hanning window
//	DSPSplitComplex fhann;
//	fhann.realp = hann_real;
//	fhann.imagp = hann_imag;
	
	
	float *result = (float*)malloc(n*sizeof(float));
	float *subbands = (float*)malloc(numbands*n*sizeof(float));
	
	float *temp = (float*)malloc(n*sizeof(float));
	
	// get log2 of the dimension for Radix2 FFT
	int dim = (int)log2((float)n);
	
	// setup the twiddle factors
	fft_weights = vDSP_create_fftsetup(dim,kFFTRadix2);
	
	Filterbank(signal, signalLength, subbands);

	// get the hanning window in the frequency domain
	vDSP_ctoz((DSPComplex *)hann, 2, &zeros, 1, n/2);
	vDSP_fft_zrip(fft_weights, &zeros, 1, dim, kFFTDirection_Forward);
	
	for(int i=0; i<numbands; i++)
	{
		subband = subbands+i*n;
		
		Smooth(subband, signalLength);
//		{
//			FILE* dump = fopen("/Users/wilperkins/Documents/Projects/openbpm/OpenBPM/log/debug.log", "w");
//			for (int i=0; i<signalLength/2; i++) fprintf(dump,"%.20f\n",subband[i]);
//			fclose(dump);
//		}
		
		DiffRect(subband, signalLength);
	}
	int bpm;
	CombFilterbank(subbands, signalLength, bpm);

	FindPeaks(signal, signalLength, beats, num_beats);
	
	free(spectrum.realp);
	free(spectrum.imagp);
	free(zeros.realp);
	free(zeros.imagp);
	free(temp);
	free(result);
	free(subbands);
	return;
}

/**
  * @brief - finds peaks (potential beats) in a signal
  */
void FindPeaks(const float signal[], uint32_t signalLength, int *beats, int *num_beats)
{
	//
	// signal must be > 2048 points!
	//
	// 1: decimate signal
	// 2: differentiate signal
	// 3: half-wave rectify signal
	// 4: locate peaks
	// 5: smooth result?
	// 6: get beat markers
	
	float 			*result, *temp, *peaks, mean;
	vDSP_Length 	*indicies;
	uint32_t		resultLength, tempLength, i, j, k;
	vDSP_Stride		signalStride, resultStride;
	int numPeaks, numBeats;
	
	
	signalStride = resultStride = 1;
	resultLength = tempLength = signalLength / DECIMATION_FACTOR;
	numPeaks = 30;
	
	
	result = (float *) malloc(resultLength * sizeof(float));
	memset(result, 0 , resultLength * sizeof(float));
	
	temp = (float *) malloc(tempLength * sizeof(float));
	memset(temp, 0 , tempLength * sizeof(float));
	
	indicies = (vDSP_Length *)malloc(resultLength * sizeof(vDSP_Length));
	
	peaks = (float *) malloc(numPeaks * sizeof(float));
	memset(peaks,0,numPeaks * sizeof(float));
	
	if (signal == NULL || temp == NULL || result == NULL || indicies == NULL || peaks == NULL) {
        printf("\nmalloc failed to allocate memory for the result.\n");
        exit(0);
    }
	
	for(i=0; i<resultLength; i++)
	{
		indicies[i] = i;
	}
	
#if DECIMATION_FACTOR > 1
	//
	// 1: decimate signal
	//
	vDSP_desamp(signal, (vDSP_Stride)DECIMATION_FACTOR, decimFilter,
				temp, (vDSP_Length)tempLength, decimFilterLength);
	
#else
	// bypass decimation
	temp = signal;
#endif
	
	//
	// 2: differentiate signal
	//
	vDSP_vsub(temp+1,signalStride,temp,signalStride,result,resultStride,tempLength-1);
	
	//
	// 3: half-wave rectify signal
	//
	
	//--vDSP_vthr(result,resultStride,zeros,result,resultStride,resultLength);
	float zero = 0.0f;
	
	//vDSP_vthr(signal,signalStride,&zero,result,resultStride,resultLength);
	
//	for(i=0; i<resultLength; i++)
//		printf("%f\n",result[i]);
//	exit(0);
	
	
	vDSP_vthr(result,resultStride,&zero,result,resultStride,resultLength);
	
	//	
	// 4: locate peaks
	//
	
	// get the X largest values
	vDSP_vsorti(result, indicies, NULL, resultLength, DESCENDING);
	for(i=0; i<numPeaks; i++)
	{
		peaks[i] = (float)indicies[i];	// do some kind of vector copy?
	}
	
	vDSP_vsort(peaks, numPeaks, DESCENDING);
	
	//
	// 5: smooth result?
	// 6: get beat markers
	//
	
	/*
	 // cluster the peaks to get beats
	 j = numBeats = 0;
	 for(i=0; i<numPeaks; i++)
	 {
	 if( (peaks[i] - peaks[j]) > SMOOTHING_LAG || i==numPeaks-1)
	 {
	 GetMoments(peaks+j, i-j+1, &mean, NULL);
	 beats[numBeats++] = mean;
	 j = i + 1;
	 }
	 }
	 */
	
	numBeats = 0;
	
	for(i=0; i<numPeaks; i++)
	{
		beats[numBeats++] = peaks[i] * DECIMATION_FACTOR;
	}
	
	*num_beats = numBeats;
	
	free(peaks);
	free(result);
	free(temp);
	free(indicies);
	
	return;	
 }

/**
  * @brief - comb filterbank
  */
void CombFilterbank(const float signal[], uint32_t signalLength, int *bpm)
{
	static int Fs = 44100;
	static int taps = 3;
	float energy, maxEnergy, bandEnergy;
	int N = (int)signalLength;
	
	*bpm = 0.0;
	float currbpm;
	float minbpm = 80.0;	// min test BPM
	float maxbpm = 160.0;	// max test BPM
	
	int M = floorf((float)(60.0/minbpm*(float)Fs));
	float *temp;	// buffer for goodies
	int tempLength = 2*M + signalLength - 1;
	
	temp = (float *) malloc(tempLength * sizeof(float));
	
//	FILE* dump = fopen("/Users/wilperkins/Documents/Projects/openbpm/OpenBPM/log/debug.log", "w");
//	for (int i=0; i<signalLength; i++) fprintf(dump,"%f\n",signal[i]);
	
	// for each bpm
	//   for each band
	//     sum the energy from all bands
	//   pick bpm w/ highest corresponding energy
	for(currbpm = minbpm; currbpm < maxbpm; currbpm = currbpm+2)
	{

		M = floorf((float)(60.0/currbpm*(float)Fs));
		tempLength = 2*M + signalLength - 1;
		energy = 0;
		
		for(int i=0; i<6; i++)
		{
			//y_n = x_n , 0 < n < N
			memcpy(temp, signal, N * sizeof(float));
	
			//y_n = x_n-2M , N < n < 2M+N
			memcpy(temp+N, signal+N-2*M, 2*M * sizeof(float));
	
			// y_n = y_n + x_n-M , M < n < M+N
			vDSP_vadd(temp+M, (vDSP_Stride)1, signal, (vDSP_Stride)1, temp+M, (vDSP_Stride)1, N);
	
			// y_n = y_n + x_n-2M , 2M < n < N
			vDSP_vadd(temp+2*M, (vDSP_Stride)1, signal, (vDSP_Stride)1, temp+2*M, (vDSP_Stride)1, N-2*M);
	
			vDSP_dotpr(temp, (vDSP_Stride)1, temp, (vDSP_Stride)1, &bandEnergy, tempLength);

			energy += bandEnergy;
		}
		//fprintf(dump, "%.20f\n",energy);
		
		if(energy > maxEnergy)
		{
			maxEnergy = energy;
			*bpm = currbpm;
		}

	}
	free(temp);
//	fclose(dump);
	return;
}

/**
  * @brief - splits a signal into several subbands
  */
void Filterbank(const float signal[], uint32_t signalLength, float subbands[])
{
	uint32_t lbound, ubound, n;

	n = 65536;
	int numbands = 6;
	//uint32_t maxfreq = 4096;
	//uint32_t bandlimits[] = {0,200,400,800,1600,3200};
	
	uint32_t maxfreq = 22050;
	uint32_t bandlimits[] = {0,1076,2153,4306,8613,17226};

	
	// get log2 of the dimension for Radix2 FFT
	int dim = (int)log2((float)n);
	
	lbound = ubound = 0;
	n = signalLength;

	//DSPSplitComplex _spectrum;
	// pack real data into split complex vector
	vDSP_ctoz((DSPComplex *)signal, 2, &spectrum, 1, n/2);
	
	vDSP_fft_zrip(fft_weights, &spectrum, 1, dim, kFFTDirection_Forward);


	// loop through the first n-1 subbands
	for(int i=1; i<numbands; i++)
	{
		// shift the boundaries
		lbound = ubound;
		ubound = (uint32_t)floorf((float)bandlimits[i]/(float)maxfreq*(float)n);
		
		// copy sig[lbound].. sig[ubound] to zeros[lbound].. zeros[ubound]
		memcpy(zeros.realp+lbound/2, spectrum.realp+lbound/2, (ubound-lbound)/2 * sizeof(float));
		memcpy(zeros.imagp+lbound/2, spectrum.imagp+lbound/2, (ubound-lbound)/2 * sizeof(float));
	
		// ifft
		vDSP_fft_zrip(fft_weights, &zeros, 1, dim, kFFTDirection_Inverse);
		
		// unpack complex into real vector
		vDSP_ztoc(&zeros, 1, (DSPComplex *)(subbands+(i-1)*n), 2, n/2);
				
		//reset the buffer
		memset(zeros.realp, 0 , n/2 * sizeof(float));
		memset(zeros.imagp, 0 , n/2 * sizeof(float));		
	}
	
	// get the last subband
	ubound = n;
	
	memcpy(zeros.realp+lbound/2, spectrum.realp+lbound/2, (ubound-lbound)/2 * sizeof(float));
	memcpy(zeros.imagp+lbound/2, spectrum.imagp+lbound/2, (ubound-lbound)/2 * sizeof(float));

	// ifft & unpack
	vDSP_fft_zrip(fft_weights, &zeros, 1, dim, kFFTDirection_Inverse);
	vDSP_ztoc(&zeros, 1, (DSPComplex *)(subbands+5*n), 2, n/2);

}

/**
  * @brief - smooths a signal using a hanning window (in frequency domain)
  */
void Smooth(float signal[], uint32_t signalLength)
{		
	// full-wave rectify the signal
	vDSP_vabs(signal, 1, signal, 1, signalLength);
	
	// get log2 of the dimension for Radix2 FFT
	int dim = (int)log2((float)signalLength);
		
	// pack real data into split complex vector
	vDSP_ctoz((DSPComplex *)signal, 2, &spectrum, 1, signalLength/2);
	
	// convert signal to frequency domain
	vDSP_fft_zrip(fft_weights, &spectrum, 1, dim, kFFTDirection_Forward);
	
	// apply hanning window
	vDSP_zvmul(&spectrum, 1, &zeros, 1, &spectrum, 1, signalLength/2, 1);
	
	// convert back to time domain
	vDSP_fft_zrip(fft_weights, &spectrum, 1, dim, kFFTDirection_Inverse);
	
	// unpack complex into real vector
	vDSP_ztoc(&spectrum, 1, (DSPComplex *)signal, 2, signalLength/2);
}

/**
 * @brief - differentiates then half-wave rectifies a signal
 */
void DiffRect(float signal[], uint32_t signalLength)
{
	float zero = 0;
	vDSP_Stride stride = 1;
	
	
	// differentiate signal
	vDSP_vsub(signal+1,stride,signal,stride,signal,stride,signalLength-1);

	// half-wave rectify signal
	vDSP_vthr(signal,stride,&zero,signal,stride,signalLength);	
}

/**
  * @brief - calculates mean & variance (1st & 2nd moments) of data vector 
  */
void GetMoments(float data[], int length, float *mean, float *var)
{	
	int i, acc;
	i = acc = 0;
	for(i=0; i<length; i++)
	{
		acc += data[i];
	}
	*mean = (float)acc/length;
	// *var = ...
	
	return;
}
