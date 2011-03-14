/** beatdetect.c
  * a lightweight beat detector
  *
  * see: http://developer.apple.com/library/mac/#documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html%23//apple_ref/doc/uid/TP40005147-CH202-SW1
  * see: http://developer.apple.com/library/mac/#documentation/Performance/Conceptual/vDSP_Programming_Guide/SampleCode/SampleCode.html%23//apple_ref/doc/uid/TP40005147-CH205-CIAEJIGF
  */

#include <stdio.h>
#include <Accelerate/Accelerate.h>

#include "beatdetect.h"	// contains filters & example signals
//#include "song.h"	// song sample

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