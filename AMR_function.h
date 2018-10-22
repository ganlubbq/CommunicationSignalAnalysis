#ifndef _AMR_FUNCTION
#define _AMR_FUNCTION
//
// Combine several header files to one header file 
// Last modified by AWH 10-Sep-2015
//
//
#if AMR_Cplus
	#include "AMR_class.h"
#endif

#if AMR_Cplus
	inline void AMRalgorithm::initializeVariable()
#else
	void initializeVariable()
#endif
{
	// Caution : check array length
	int i=0;

	for(i=0; i<AMR_MAXArray; i++){
		AMR_filteredInputData[i]	= 0.0;
	}

	for(i=0; i<AMR_MAX_NUMBINS+2; i++){
		AMR_binCenter[i] = 0.0;
		AMR_bins[i]		 = 0.0;	// 구역 분할 범위, numAMR_bins의 중심값
		AMR_binCount[i]	 = 0;		
	}

	for(i=0; i<AMR_symbolRate_SamplesPerSegment>>1; i++){
		AMR_halfFreqVector[i]    = 0.0;
	}

	for(i=0; i<AMR_symbolRate_SamplesPerSegment; i++){
		AMR_iFrequency[i]		 = 0.0;
		AMR_iPhase[i]      		 = 0.0;
		AMR_ck[i] 				 = 0.0;
		AMR_filteredSamples[i]   = 0.0;	
		AMR_spectrumMagnitude[i] = 0.0;
		AMR_spectrumBuffer[i]    = 0.0;
		AMR_iAmplitude[i]      	 = 0.0;		

		// for moving average filter used for symbol rate estimation
		//AMR_MovInB[i]			= 0.0;
		//AMR_cbegin1[i]			= 0.0;
		//AMR_cbegin2[i]			= 0.0;
		//AMR_cend1[i]			= 0.0;
		//AMR_cend2[i]			= 0.0;
		//AMR_Movfiltout[i]		= 0.0;
	}
	memset(AMR_MovInB,0x0,sizeof(AMR_MovInB));
	memset(AMR_cbegin1,0x0,sizeof(AMR_cbegin1));
	memset(AMR_cbegin2,0x0,sizeof(AMR_cbegin2));
	memset(AMR_cend1,0x0,sizeof(AMR_cend1));
	memset(AMR_cend2,0x0,sizeof(AMR_cend2));
	memset(AMR_Movfiltout,0x0,sizeof(AMR_Movfiltout));

	for(i=0; i<((AMR_MAX_welchPSDsize<<1) + 1); i++) AMR_IQSegment_PSD[i]   = 0.0;
		
	for(i=0; i<AMR_MAX_welchPSDsize; i++);
	{
		AMR_powerSpectrum[i]   			= 0.0;
		AMR_powerSpectrum_dB[i] 		= 0.0;
		AMR_Sxx_PSD[i] 			  		= 0.0;
		AMR_hanningWindow_PSD[i]    	= 0.0;

		AMR_freqVector[i] 	  			= 0.0;			
		AMR_TSE_PSDbuffer[i]			= 0.0;
		AMR_TSE_indexBuffer[i]         	= 0;		
		AMR_IndexBuffer_PSD[i]          = 0;
	}

	//memset(AMR_powerSpectrum,0,sizeof(AMR_powerSpectrum));
	//memset(AMR_powerSpectrum_dB,0,sizeof(AMR_powerSpectrum_dB));
	//memset(AMR_Sxx_PSD,0,sizeof(AMR_Sxx_PSD));
	//memset(AMR_hanningWindow_PSD,0,sizeof(AMR_hanningWindow_PSD));
	//memset(AMR_freqVector,0,sizeof(AMR_freqVector));

	for(i=0; i<(AMR_MAX_welchPSDsize>>1); i++){
		AMR_peakValues[i]				= 0.0;
		AMR_peakIndices[i] 				= -1;
	}
	//for bandwidth estimation
	for(i=0; i<2; i++){
		AMR_BWSampleIdx[i] = 0;
		AMR_BWSampleIdx_TMP[i] = 0;
	}

	//for LPF
	for(i=0; i<AMR_LPFLEN; i++) {
		AMR_hammingWindow[i] = 0.0;
		AMR_LPFCoeff[i] = 0.0;
		AMR_h[i] = 0.0;
	}
	for(i=0; i<AMR_LPFLEN>>1; i++) {
		AMR_k[i] = 0.0;					//for LPF
		AMR_b[i] = 0.0;					//for LPF
		AMR_a[i] = 0.0;					//for LPF
	}
	for(i=0; i<AMR_LPFLEN>>3; i++) {
		AMR_F[i] = 0.0;				    //for LPF
		AMR_M[i] = 0.0;				    //for LPF
	}

	// for features
	//for(i=0 ; i<AMR_MAXNumOfFeature; i++){
	//	AMR_sigma_a[i]				= 0.0f;
	//	AMR_absCC20[i] 				= 0.0f;
	//	AMR_absCC40[i]				= 0.0f;
	//	AMR_lineLMSerror[i]         = 0.0f;
	//	AMR_dataConditioningTmp[i]	= 0.0f;
	//	AMR_mu42f[i]                = 0.0f;
	//}
	memset(AMR_sigma_a,0.0f,sizeof(AMR_sigma_a));
	memset(AMR_absCC20,0.0f,sizeof(AMR_absCC20));
	memset(AMR_absCC40,0.0f,sizeof(AMR_absCC40));
	memset(AMR_lineLMSerror,0.0f,sizeof(AMR_lineLMSerror));
	memset(AMR_dataConditioningTmp,0.0f,sizeof(AMR_dataConditioningTmp));
	memset(AMR_mu42f,0.0f,sizeof(AMR_mu42f));

	//  for cyclic cumulant
	for (i=0; i<((AMR_CC_samplesPerSegment<<1) + 1); i++){
		AMR_IQsegmentCC[i]   = 0.0;
		AMR_IQbufferCC[i]	 = 0.0;
		AMR_IQDataTmp2[i]	 = 0.0;  // nonlinear transform
		AMR_IQDataTmp1[i]	 = 0.0;  // nonlinear transform
		AMR_IQDataResQ[i]	 = 0.0;  // nonlinear transform
	}
	for (i=0; i<AMR_CC_samplesPerSegment; i++){
		AMR_spectrumMagnitudeCC[i] = 0.0;
	}
	// IQ segmentation
	for (i=0; i<AMR_MAX_NUM_SEGMENT; i++){
		beginTOAidx[i] = 0;
		endTOAidx[i] = 0;
	}

#if VISUALSTUDIO
	for (i=0; i<2; i++){
		AMR_interpFiltOut[i] = 0.0;
		AMR_fineCompOut[i]   = 0.0;
	    AMR_pDelayBuffer3[i] = 0.0;
	    AMR_pDelayBuffer2[i] = 0.0;
	    AMR_pDelayBuffer1[i] = 0.0;
		AMR_pTEDDelay1[i]    = 0.0;  // interpolant buffer
		AMR_pTEDDelay2[i]    = 0.0;
	}
	for(i=0; i<AMR_MAXArray; i++){
		AMR_timingLoopOut[i]		= 0.0;
	}

	//Demodulation
	for(i=0; i<IQ_NUM_MAX; i++){
		AMR_FMdemodiPhase[i] = 0.0f;
		AMR_FMdemodiFreq[i] = 0.0f;
		AMR_FMdemodCk[i] = 0.0f;
		AMR_wavBuffer[i] = 0;
		AMR_FMdemodiFreqBuff[i] = 0.0f;
	}

	for(i=0; i<IQ_NUM_MAX; i++){
		AMR_FMdemodiPhase[i] = 0.0f;
		AMR_FMdemodiFreq[i] = 0.0f;
		AMR_FMdemodCk[i] = 0.0f;
	//	AMR_wavBuffer[i] = 0.0f;
	}

	for(i=0; i<(IQ_NUM_MAX<<1)+1; i++){
		AMR_demodIQ[i] = 0.0f;
	}

#endif

	AMR_coarseBandWidth	  = 0.0f;
	AMR_coarseFreqOffset  = 0.0f;
	AMR_coarseSymbolRate  = 0.0f;
	AMR_coarseToneSpacing = 0.0f;
	AMR_coareseSNR_dB     = 0.0f;
	AMR_requireIQsamples  = 0;
	AMR_decimationIdx     = 0;

	AMR_modNum = 0;
	AMR_modOrder= 0;
	AMR_LEN = 0;

	AMR_symratePeak        = 0.0f;   
    AMR_symratePeakIdx     = 0; 
}


#ifndef activateCDW  // for visual studio
	float log2(float lbX){
		return log(lbX) / log(2.0f);
	}
#endif

#if AMR_Cplus	
	inline signed int AMRalgorithm::sign(float inSIGN)
#else
	signed int sign(float inSIGN)
#endif
{

	if (inSIGN > 0.0f){
		return 1;
	}
	else if(inSIGN <0.0f){
		return -1;
	}
	else return 0;
}


#if AMR_Cplus	
	inline float AMRalgorithm::sinc(float inX)
#else
	float sinc(float inX)
#endif

{
/********************************************************************************/
// \ DESCRIPTION
//  > Calculate value of sinc function. It is used to calculate the coefficients of LPF
/********************************************************************************/
	float pi 	 = AMR_PI;   
	float den    = 0.0;
	float tmp    = 0.0;
	float sinVal = 0.0;

	if (inX == 0)
		return 1;
	else
	{
		den = pi*inX;
		sinVal = sin(pi*inX);
		tmp = sinVal / den;
		return tmp;
	}
}

#if AMR_Cplus
	inline float AMRalgorithm::myCabs(float re, float im)
#else
	float myCabs(float re, float im)
#endif
{
/******************************************************************************/
// \ DESCRIPTION
//  > Calculate absolute value of complex value
//
// \ INPUT ARGUMENTS
//	> re : real part 
//  > im : imaginary part
//
// \ OUTPUT ARGUMENTS
//  > ans  : absolute value 
// 
// \ Reference
//  > http://www.nr.com/public-domain.html, complex.c, complex.h
/******************************************************************************/
	// 

	float x,y,ans,temp;
	x=fabs(re);
	y=fabs(im);
	if (x == 0.0f)
		ans=y;
	else if (y == 0.0f)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0f+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0f+temp*temp);
	}
	return ans;
}
#if AMR_Cplus
	inline int AMRalgorithm::myCompare(const void *apple, const void *banna)
#else
	int myCompare(const void *apple, const void *banna)
#endif
// for qsort
{
	//descending, if you want to ascending order, swap banna and apple
	if ( *(float*)banna <  *(float*)apple ) return -1;
	else if	( *(float*)banna >  *(float*)apple ) return 1;
	else return 0;
	
	//return ( *(float*)apple - *(float*)banna );
}

#if AMR_Cplus
	inline int AMRalgorithm::myCompareInt(const void *apple, const void *banna)
#else
	int myCompareInt(const void *apple, const void *banna)
#endif
{
	//descending, if you want to ascending order, swap banna and apple
	if ( *(int*)banna <  *(int*)apple ) return -1;
	else if ( *(int*)banna >  *(int*)apple ) return 1;
	else return 0;
	
}

#if AMR_Cplus
	inline float AMRalgorithm::medianF(float* inMedian,int startIdxMed, int endIdxMed )
#else
	float medianF(float* inMedian,int startIdxMed, int endIdxMed )
#endif
{
	int medianIdx = 0;
	int L = endIdxMed-startIdxMed+1;
	float medianVal = 0.0f;
	medianIdx = startIdxMed + (L >>1);		// implict type casting : float -> int

	if (L % 2 == 1){  											// 요소 개수가 홀수면
		medianVal = inMedian[medianIdx]; 				// 홀수 개수인 배열에서는 중간 요소를 그대로 반환
	}else{
		// Ignore floating data
		medianVal = inMedian[medianIdx - 1] + inMedian[medianIdx]; // 짝수 개 요소는, 중간 두 수의 평균 반환
		medianVal /= (2.0f);
	}

	 return medianVal;
}
#if AMR_Cplus
	inline int AMRalgorithm::medianI(int* mediInput,int startIdxMed2, int endIdxMed2)
#else
	int medianI(int* mediInput,int startIdxMed2, int endIdxMed2)
#endif

{
	signed int medianIdx = 0;
	signed int L = endIdxMed2-startIdxMed2+1;
	signed int medianVal =0;

	medianVal = L / 2;
	medianIdx = startIdxMed2 + medianVal;

	if (L % 2 == 1){  											// 요소 개수가 홀수면
		medianVal = mediInput[medianIdx]; 				// 홀수 개수인 배열에서는 중간 요소를 그대로 반환
	}else{
		// Ignore floating data
		medianVal = mediInput[medianIdx - 1] + mediInput[medianIdx]; // 짝수 개 요소는, 중간 두 수의 평균 반환
		medianVal /= 2;
	}

	 return medianVal;
}

#if AMR_Cplus
	inline void AMRalgorithm::four1(float* data, unsigned long nn, int isign)
#else
	void four1(float* data, unsigned long nn, int isign)
#endif

/*********************************************************************************
 FFT/IFFT routine. (see pages 507-508 of Numerical Recipes in C)

 Inputs:
	data[] : array of complex* data points of size 2*AMR_NFFT+1.
		data[0] is unused,
		* the n'th complex number x(n), for 0 <= n <= length(x)-1, is stored as:
			data[2*n+1] = real(x(n))
			data[2*n+2] = imag(x(n))
		if length(Nx) < AMR_NFFT, the remainder of the array must be padded with zeros

	nn : FFT order AMR_NFFT. This MUST be a power of 2 and >= length(x).
	isign:  if set to 1,
				computes the forward FFT
			if set to -1,
				computes Inverse FFT - in this case the output values have
				to be manually normalized by multiplying with 1/AMR_NFFT.
 Outputs:
	data[] : The FFT or IFFT results are stored in data, overwriting the input.

  from : http://www-ee.uta.edu/eeweb/ip/Courses/DSP_new/Programs/fft.cpp
************************************************************************************/
// double type is recommended as input

{
    int n, mmax, m, j, istep, i, nd2;
    //float  tmpf1 =0.0, tmpf2=0.0;
	//double wtemp, wr, wpr, wpi, wi, theta;	//Dobule precision for the trigonometric recurrences.
    //double tempr, tempi;
	float wtemp, wr, wpr, wpi, wi, theta;	//Dobule precision for the trigonometric recurrences.
    float tempr, tempi;
	float twoPI = AMR_TWOPI;
    n = nn << 1; nd2 = nn >> 1;
    j = 1;
    for (i = 1; i < n; i += 2) {			// This is bit-reversal section of the routine
		if (j > i) {						// Exchange the two complex numbers.
			tempr = data[j];     data[j] = data[i];     data[i] = tempr;
			tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
    }

    mmax = 2;
    while (n > mmax) {					//Outer loop executed log2 n times
		istep = 2*mmax;
		theta = twoPI/(isign*mmax);		//Initialized the trigonometric recurrence
		wtemp = sin(0.5f*theta);
		wpr = -2.0f*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0f;
		wi = 0.0f;
		for (m = 1; m < mmax; m += 2) {					// Here are the two neted inner loops
			for (i = m; i <= n; i += istep) {			// This is the Danielson-Lanczos formula
				j =i + mmax;
				tempr = wr*data[j]   - wi*data[j+1];
				tempi = wr*data[j+1] + wi*data[j];
				data[j]   = data[i]   - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi*wpi + wr;		// Trigonometirc recurrence
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}

    // Because fft operation of this method is calculated in reverse order with respect to MATLAB.
    // To relocate sequence order, we add below codes.
    // ex)
    // 	    [input]
    //                  DC											 Center
    // 		real part : 1     3     5     7     9    11    13    15    17    19    21    23    25    27    29    31
    //      imag part : 2     4     6     8    10    12    14    16    18    20    22    24    26    28    30    32
    //
    //      [output]
    // 		real part : 1    31    29    27    25    23    21    19    17    15    13    11     9     7     5     3
    //      imag part : 2    32    30    28    26    24    22    20    18    16    14    12    10     8     6     4
	
	// Put below code outside 
    for(i=1; i<=(nd2-1); i++){   // except DC
    	AMR_SWAP_FLOAT(data[2*i+1], data[n - (2*i-1)]);  // real
    	AMR_SWAP_FLOAT(data[2*i+2], data[n - (2*i-2)]);  // imag
    }
}

#if AMR_Cplus
	inline int AMRalgorithm::complexMultiply(float* inComp, float* inRef, unsigned int compSigLen, signed int compNum)
#else
	int complexMultiply(float* inComp, float* inRef, unsigned int compSigLen, signed int compNum)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  > Multiply inComp and inRef in the manner of recursive function
//
// \ INPUT ARGUMENTS
//	> inComp : Complex array,
//             As start from a second element,
//             odd elements are real number and even elements are imaginary number
//  > inRef      : Similar to inComp
//  > compSigLen : the length of input data
//  > compNum    : compNum-th order multiplication of two complex numbers.
//
//
// \ OUTPUT ARGUMENTS
//  > inComp  : Overwrite results to input array
//
// \ Author(s) : AWH
// 
// \ Reference
/******************************************************************************/

	unsigned int i=0;
	float tmp =0.0;
	if (compNum <=1) return 0;
	else{
		for(i=0; i<compSigLen; i++){
			tmp = inComp[2*i+1] * inRef[2*i+1] - inComp[2*i+2] * inRef[2*i+2];
			inComp[2*i+2] = inComp[2*i+1] * inRef[2*i+2] + inComp[2*i+2] * inRef[2*i+1];
			inComp[2*i+1] = tmp;
		}
		return	complexMultiply(inComp, inRef, compSigLen, compNum-1);
	}
}
#if AMR_Cplus
	inline int AMRalgorithm::nonlinearTransform(float* inNT,signed int ntLen, signed int nTn, signed int nTq)
#else
	int nonlinearTransform(float* inNT,signed int ntLen, signed int nTn, signed int nTq)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  > Implementation on power operations of a complex array.
// 	   Below codes are equal to ( nonlinearTransform = x.^(n-q) .* conj(x).^q ) in MATLAB
//
// \ INPUT ARGUMENTS
//	> inNT  - Complex array,
//            As start from a second element,
//            odd elements are real number and even elements are imaginary number
//  > ntLen - The length of input array
//  > nTn   - N-th power
//  > nTq   - q conjugate
//
// \ OUTPUT ARGUMENTS
//	> inNT  : Overwrite results to input array
//
// \ Author(s) : AWH
//
/*******************************************************************************/

	unsigned int i    = 0;
	unsigned int L	 = ntLen;
	signed int nMinusQ = 0;

	nMinusQ = nTn-nTq;
	if(nMinusQ < 0 ) return 0;

	//initialize
	for (i= 0 ; i<L; i++){
		AMR_IQDataResQ[2*i+1] = inNT[2*i+1];			// real
		AMR_IQDataResQ[2*i+2] = (-1.0f)*inNT[2*i+2];	//imag

		AMR_IQDataTmp1[2*i+1] = inNT[2*i+1];
		AMR_IQDataTmp1[2*i+2] = (-1.0f)*inNT[2*i+2];

		AMR_IQDataTmp2[2*i+1] = inNT[2*i+1];
		AMR_IQDataTmp2[2*i+2] = inNT[2*i+2];
	}

	////////// (n-q)-th order
	if( nMinusQ > 0 ) complexMultiply(inNT, AMR_IQDataTmp2, L, nMinusQ);

	////////// q-th conjugate
	if(nTq >0){
		complexMultiply(AMR_IQDataResQ, AMR_IQDataTmp1, L, nTq);

		if(nMinusQ){
			////////// multiply   (n-q)th order by qth order
			complexMultiply(inNT, AMR_IQDataResQ, L, 2);
		}else{	// in case of zero
			for(i=0; i<L; i++){ // copy AMR_resQ to inNT
				inNT[2*i+1] = AMR_IQDataResQ[2*i+1];
				inNT[2*i+2] = AMR_IQDataResQ[2*i+2];
			}
		}

	}
	return 1;
}
#if AMR_Cplus
	inline float AMRalgorithm:: dataConditining(float* dcInput,int dcInputLen)
#else
	float dataConditining(float* dcInput,int dcInputLen)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  > Remove outliers in feature vector
//
// \ INPUT ARGUMENTS
//  > dcInput : feature array
//  > dcInputLen : the length of feature array 
//
// \ OUTPUT ARGUMENTS
//  > q2 : Median number
//  
// \ Example
//
// \ See also boxplot in MATLAB
//
// \ Author(s) : AWH
// 
// \ Reference
/******************************************************************************/

	int L = dcInputLen;
	int startIdx = 0, endIdx = 0;
	//float tmpf=0.0, tmpf2=0.0;

	//for boxplot
	float q2=0.0f; //, q1=0.0, q3=0.0, IQR=0.0;

	/*for(i=0; i< dcInputLen; i++){
		AMR_dataConditioningTmp[i] = dcInput[i];
	}*/

	#if AMR_Cplus
		heapsort(dcInput, (int *)dcInput, L); //descending order
	#else
		qsort(dcInput, L, sizeof(float),myCompare); // ascending order
	#endif

	startIdx = 0; endIdx = L-1;
	// compute 50th percentile (second quartile)
	q2 = medianF(dcInput,startIdx,endIdx);

	//// compute 25th percentile (first quartile)
	//for(i=0; i<L; i++){
	//	if(dcInput[i] <= q2){
	//		startIdx = i; endIdx = L-1;
	//		q1 =  medianF(dcInput,startIdx,endIdx);
	//		break;
	//	}
	//}
	//// compute 75th percentile (third quartile)
	//for(i=0; i<L; i++){
	//	if(dcInput[L-i-1] >= q2){
	//		startIdx = 0; endIdx = L-i-1;
	//		q3 =  medianF(dcInput,startIdx,endIdx);
	//		break;
	//	}
	//}
	//// compute Interquartile Range (IQR)
	//IQR = q3-q1;

	//tmpf = q1-(1.5f*IQR);
	//tmpf2 = q3+(1.5f*IQR);
	//// Replace upper outlier and lower outlier with median value
	//for(i=0; i<L; i++){
	//	if((dcInput[i] > tmpf2) || (dcInput[i] < tmpf)){
	//		dcInput[i] = q2;//
	//	}
	//}

	return q2;
}

#if AMR_Cplus
	inline void AMRalgorithm::histogram(float* inHist, unsigned int inHistLength, unsigned int histNumBin, float maxy, float miny)
#else
	void histogram(float* inHist, unsigned int inHistLength, unsigned int histNumBin, float maxy, float miny)
#endif

{

/********************************************************************************/
// \ DESCRIPTION
//  > Calculate the histogram
//
// \ INPUT ARGUMENTS
//	> inHist : input data
//  > inHistLength : the length of input data
//  > histNumBin : the number of bins 
//  > maxy : Maximum value of input data
//  > miny : Minimum value of input data
//
// \ OUTPUT ARGUMENTS
//
// \ Example
//
// \ See also hist.m in MATLAB 
//
// \ Author(s) : Won yu-jun and edited by AWH
// 
// \ Reference
//		
// \ Note : AMR_binCount의 마지막 두 비트는 의미가 없는 데이터이다.
//  center  : the number of samples within one bin interval
/******************************************************************************/

	unsigned int i=0,idx=0;
//	float	maxy=0.0, miny=0.0
	float 	binwidth = 0.0;

	// Find the largest elements and smallest elements in array --> put this code out of function - 140108
//	maxy = inHist[0];
//	for (idx=1 ; idx < inHistLength ; idx+= 1){
//		if(inHist[idx] > maxy ){
//			maxy = inHist[idx];
//		}
//	}
//	miny = maxy;
//	for (idx=0 ; idx < inHistLength ; idx+= 1){
//		if(inHist[idx] < miny ){
//			miny = inHist[idx];
//		}
//	}
	for (idx=0; idx<histNumBin+2; idx++)	AMR_binCount[idx] = 0;

	//tmp =  ; // +- 0.3 error
	binwidth =  (maxy - miny) / (float) histNumBin;

	for (idx=0 ; idx < histNumBin+1 ; idx+= 1)	AMR_bins[idx] = miny + (binwidth*idx);

	AMR_bins[histNumBin] = maxy; // ok

	for (idx=0 ; idx < histNumBin ; idx+= 1) AMR_binCenter[idx] = AMR_bins[idx] + (binwidth / 2.0f); //ok

	// Shift AMR_bins so the interval is ( ] instead of [ ).

	// handles with input samples lower than first bin number
	// -inf <=   < AMR_bins[0]
	for (idx = 0; idx<inHistLength;idx++){
		if(inHist[idx] < AMR_bins[0]){
				AMR_binCount[0] += 1;
		}
	}

	// <=   <
	for (i=1; i<=histNumBin; i++){
		for (idx = 0; idx<inHistLength;idx++){
			if(( AMR_bins[i-1] <= inHist[idx]) && (inHist[idx] < AMR_bins[i]))
				AMR_binCount[i] += 1;
		}
	}

	// handles with input samples greater than last bin
	// AMR_bins[histNumBin] <=
	for (idx = 0; idx<inHistLength;idx++){
		if(AMR_bins[histNumBin] <= inHist[idx]){
				AMR_binCount[histNumBin+1] += 1;
		}
	}

	// Combine first bin with 2nd bin and last bin with next to last bin
	AMR_binCount[1] = AMR_binCount[0] + AMR_binCount[1];
	AMR_binCount[histNumBin] = AMR_binCount[histNumBin] + AMR_binCount[histNumBin+1];

	// shift array by one samples in left
	for (i=0; i<histNumBin; i++){
		AMR_binCount[i] = AMR_binCount[i+1];
	}
}

#if AMR_Cplus
	inline void AMRalgorithm::siftDown( float *a, int *idx, int start, int end) // Call by Reference
#else
	void siftDown( float *a, int *idx, int start, int end) // Call by Reference
#endif
{
//void siftDown( ValType a[], int idx[], int start, int end){
	int root = start;

    while ( root*2+1 < end ) {

	    int child = root*2 + 1;

        //if ((child + 1 < end) && AMR_IS_LESS(a[child+1],a[child])) {
		if ((child + 1 < end) && (a[child+1] < a[child])) {
            child += 1;
        }

		//if (AMR_IS_LESS(a[child], a[root])) {
		if (a[child] < a[root]) {
            AMR_SWAP_FLOAT(a[root], a[child]  );
			AMR_SWAP_INT( idx[root], idx[child] );
            root = child;
        }
        else
            return;
    }
}

#if AMR_Cplus
	inline void AMRalgorithm::siftDownInt( int *a, int *idx, int start, int end) // Call by Reference
#else
	void siftDownInt( int *a, int *idx, int start, int end) // Call by Reference
#endif
{
	int root = start;

    while ( root*2+1 < end ) {

	    int child = root*2 + 1;

        //if ((child + 1 < end) && AMR_IS_LESS(a[child+1],a[child])) {
		if ((child + 1 < end) && (a[child+1] < a[child])) {
            child += 1;
        }

		//if (AMR_IS_LESS(a[child], a[root])) {
		if (a[child] < a[root]) {
            AMR_SWAP_INT( a[root], a[child]  );
			AMR_SWAP_INT( idx[root], idx[child] );
            root = child;
        }
        else
            return;
    }
}

#if AMR_Cplus
	inline void AMRalgorithm::heapsort( float *a, int *idx, int count) // descend order
#else
	void heapsort( float *a, int *idx, int count) // descend order
#endif
{
/******************************************************************************/
// \ DESCRIPTION
//  > heap sort algorithm
//
// \ INPUT ARGUMENTS
//	> a   : input samples
//  > idx : index of input samples
//  > count : The number of input samples
//
// \ OUTPUT ARGUMENTS
//  > a  : overwrite results to the input  
// 
// \ Reference
//  > http://en.wikipedia.org/wiki/Heapsort
/******************************************************************************/
	
    int start, end;

	for (start = (count-2)/2; start >=0; start--) {
        siftDown( a,idx, start, count);
    }

    for (end=count-1; end > 0; end--) {
        AMR_SWAP_FLOAT(a[end],a[0]);
		AMR_SWAP_INT(idx[end],idx[0]);

        siftDown(a, idx, 0, end);
    }
}

#if AMR_Cplus
	inline void AMRalgorithm::heapsortInt( int *a, int *idx, int count) // descend order
#else
	void heapsortInt( int *a, int *idx, int count) // descend order
#endif
{
	// http://en.wikipedia.org/wiki/Heapsort
    int start, end;

	for (start = (count-2)/2; start >=0; start--) {
        siftDownInt( a,idx, start, count);
    }

    for (end=count-1; end > 0; end--) {
        AMR_SWAP_INT(a[end],a[0]);
		AMR_SWAP_INT(idx[end],idx[0]);

        siftDownInt(a, idx, 0, end);
    }
}

#if AMR_Cplus
	inline void AMRalgorithm::fftshift(float* infftshift, float* outfftshift, unsigned int lenfftsh)
#else
	void fftshift(float* infftshift, float* outfftshift, unsigned int lenfftsh)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - Shift zero-frequency component to center of spectrum
//
// \ INPUT ARGUMENTS
//	> infftshift	- input spectrum
//	> outfftshift	- output spectrum 
//	> lenfftsh		- The number of FFT
//  
// \ OUTPUT ARGUMENTS
//
// \ Example
//	> input  = [ 1 2 3 4 5 6 7 8 9 10]
//	> output = [ 6 7 8 9 10 1 2 3 4 5]
//
// \ See also fftshift in MATLAB
//
// \ Author(s) : AWH
// 
// \ Reference
/******************************************************************************/


	unsigned int idx=0, p=0;
	p = lenfftsh >> 1; // lenfftsh is even, so p is odd
	//DC value
	//infftshift[1]
	for(idx=0 ; idx< lenfftsh; idx++){
		// 0 1 2 3 4
		if(idx < p)		outfftshift[idx] = infftshift[p+idx];
		// 5 6 7 8 9
		else 			outfftshift[idx] = infftshift[idx-p];	
	}
}

#if AMR_Cplus
	inline void AMRalgorithm::fit(float x[], float y[], unsigned int ndata, float *a,
	float *b, float *siga, float*sigb, float *lmsErr)
#else
	void fit(float x[], float y[], unsigned int ndata, float *a,
	float *b, float *siga, float*sigb, float *lmsErr)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  > Given a set of data points x[1..ndata],y[1..ndata] fit them to 
//    a straight line y=a+bx by minimizing chi square. 
//
// \ INPUT ARGUMENTS
//	> x : 
//  > y : 
//  > ndata : length of input array
//
// \ OUTPUT ARGUMENTS
//  > siga : probable uncertainity
//	> sigb : probable uncertainity
//	> lmsErr : LMS error
//
// \ Example
//	// y belong to data2.h
//	int i=0;
//	float x[16385];
//	for(i=0; i<16385; i++)		x[i] = i;
//	int ndata = 16384;
//	float a, b, siga, sigb, chi2, q;
//	fit(x, y, ndata, &a, &b, &siga, &sigb, &lmsErr);

// \ See also fit.m in MATLAB
//
// \ Reference
//
//   [1] Numerical Recipes in C, 2ed, pp.661 - 666
//
// \ Additional information
//	 >In order to fitting the curve of data, matlab use QR decomposition 
//	  based on householder matrix. Plz. refer to Householder.m, qr.m
//
/******************************************************************************/

	unsigned int i;
	float t, sxoss, sx=0.0, sy=0.0, st2=0.0, ss=0.0f, chi2=0.0f;

	*b = 0.0f;

	for(i=0; i<ndata; i++){		
		sx += x[i];
		sy += y[i];	
	}
	ss = (float)ndata;	
	sxoss = sx/ss;

	for(i=0; i<ndata; i++){
		t = x[i]-sxoss;
		st2 += t*t;
		*b += t*y[i];
	}

	*b /= st2;			// slove for a,b, sigm_a, and sigma_b
	*a = (sy-sx*(*b))/ss;
	*siga = sqrt((1.0f+sx*sx/(ss*st2))/ss);
	*sigb = sqrt( 1.0f / st2);

	for(i=0; i<ndata; i++)	
		//chi2 += AMR_SQR(y[i]-(*a)-(*b)*x[i]);		// Calculate chi square
		chi2 += (y[i]-(*a)-(*b)*x[i])*(y[i]-(*a)-(*b)*x[i]);		// Calculate chi square
	*lmsErr = sqrt((chi2)/ (ndata-2));	//For unwighted data evaluate typical sig using chi2,
	*siga *= (*lmsErr);					//and adjust the standard deviations.
	*sigb *= (*lmsErr);

}


#if AMR_Cplus	
	inline void AMRalgorithm::unitPower(float* inIQ, unsigned int L)
#else
	void unitPower(float* inIQ, unsigned int L)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - Normalize signal power to 1
//
// \ INPUT ARGUMENTS
//	> inIQ	     - input IQ samples 
//	> L			 - the number of input samples
//
// \ OUTPUT ARGUMENTS
//	> inIQ		- overwrite results to input
//
// \ Author(s) : AWH
// 
// \ Reference
//
/******************************************************************************/
	unsigned int idx = 0;
	float sigPower = 0.0, sigPowerNormFactor = 0.0;

	#if AMR_DEBUG
		printf("\r\n #unitPower \r\n");
	#endif
	//// out = sqrt( length(bandLimitedSig))/norm(bandLimitedSig, 2) * bandLimitedSig; in MATLAB
	for(idx=0; idx<L; idx++){
		//tmpNorm = AMR_SQR(inIQ[2*idx+1]) + AMR_SQR(inIQ[2*idx+2]);
		//tmpNorm = (inIQ[2*idx+1]*inIQ[2*idx+1]) + (inIQ[2*idx+2]*inIQ[2*idx+2]);
		//sigPower += tmpNorm;
		sigPower += (inIQ[2*idx+1]*inIQ[2*idx+1]) + (inIQ[2*idx+2]*inIQ[2*idx+2]);
	}
	sigPowerNormFactor =  sqrt((float)L / sigPower);	
	for(idx=0; idx<L; idx++)	{
	// 0715 edited
		inIQ[2*idx+1] *=  sigPowerNormFactor;
		inIQ[2*idx+2] *=  sigPowerNormFactor;
	}
}

#if AMR_Cplus
	inline void AMRalgorithm::welchPSD(float* inPSD, unsigned int inPSDlength, unsigned int FFTlength, signed int inPower, float* outPSD)
#else
	void welchPSD(float* inPSD, unsigned int inPSDlength, unsigned int FFTlength, signed int inPower, float* outPSD)
#endif

{
/********************************************************************************/
// \ DESCRIPTION
//  > Welch's power spectrum, estimates the power spectrum density using Welch's method
//
// \ INPUT ARGUMENTS
//  > inPDF 	  - input IQ data
//  > inPSDlength - the number of IQ samples
//	> inPower     - power 
//
// \ OUTPUT ARGUMENTS
//  > outPSD   	  - average PSD
//
// \ Used global variable
//  > AMR_hanningWindow_PSD, AMR_IQSegment_PSD, AMR_Sxx_PSD
//
// \ Used sub-function
//  > nonlinearTransform, four1, fftshift
//
// \ Example
//
// \ See also PWELCH in MATLAB
//
// \ Author(s) : AWH
// 
// \ Reference
/******************************************************************************/
	
	unsigned int tmp  = 0;
	unsigned int L 	  = FFTlength;
	unsigned int noverlap 	  = (L >> 1);					// length of the segments
	unsigned int LminusOverlap = L-noverlap;
	unsigned int tmp2 =0; //		  = inPSDlength-noverlap;

	// check this with global variables
	unsigned int numSeg =0; //		  = tmp2 / LminusOverlap; // the number of segments
	unsigned int i=0,j=0,k=0,segStart=0,segEnd=0;

	float xx=0.0;        // constant for hanning window
	float uu=0.0;        // constant for hanning window
	float pi = AMR_PI;   //

    float sumWin = 0.0;    // sum of window
	float Sxx = 0.0f;

	if(inPSDlength < noverlap){
	}else{
		tmp2   = inPSDlength-noverlap;
		numSeg = tmp2 / LminusOverlap; // the number of segments
		// noverlap개의 표본수가 되지 않는 마지막 segment는 그냥 버린다. 복잡해서
	}

    #if AMR_DEBUG
		printf("\r\n #PSD \r\n");
	#endif

	//AMR_DEBUG_numseg  = tmp2 / LminusOverlap;
	
	//tmp = 2*L+1;
	//for(i=0; i<tmp;  i++)	AMR_IQDataSegment[i] = 0.0;

	tmp = L+1;
	uu =  (float) 1.0f / tmp;
	//AMR_DEBUG_floatTmp = uu;

	for(j=0; j<L; j++)	AMR_Sxx_PSD[j] = 0.0f;

	for(j=1; j<=noverlap; j++) {
		//sumWin += (WINDOW(j,facm,facp));
		xx=(float)j * uu;	
		AMR_hanningWindow_PSD[j-1] = ( 0.5f - 0.5f*cos(2*pi*xx));		// symmetric hanning window
		AMR_hanningWindow_PSD[L-j] = AMR_hanningWindow_PSD[j-1];
		sumWin += AMR_hanningWindow_PSD[j-1] * AMR_hanningWindow_PSD[j-1] ; // sum of hanning window
	}
	sumWin = 2.0f*sumWin; // e+4 order

	if(numSeg > 19) numSeg = 19; // to increase the AMR speed, recommend odd number at 150610

	if(numSeg){
		for(i=0; i<numSeg; i++) {
			segStart = i*LminusOverlap;
			segEnd = segStart + L;

			//Apply the window to the data.
			k=1;
			for(j=segStart+1 ; j<= segEnd ; j++){
				AMR_IQSegment_PSD[2*k-1] = inPSD[2*j-1];
				AMR_IQSegment_PSD[2*k]	= inPSD[2*j];
				k++;
			}
			
			//Take nonLinear transform
			if(inPower >1) nonlinearTransform(AMR_IQSegment_PSD,L,inPower,0);
			
			//Windowing
			for(j=0 ; j<L ; j++){

				// ex) FFTlength(L) = 8
				//	2*j+1 : 	  1      3     5     7     9    11    13    15    
				//	2*j+2 :       2      4     6     8    10    12    14    16   
		
				//real
				AMR_IQSegment_PSD[2*j+1]  *= AMR_hanningWindow_PSD[j]; //  e-7 order

				//imag
				AMR_IQSegment_PSD[2*j+2]  *= AMR_hanningWindow_PSD[j];
			}

			four1(AMR_IQSegment_PSD,L,1); // e-1 ~ e+1 order

			// MaxHold trace 
			// evaluate magnitude of spectrum
	//		if (i > 0){
	//			for(j=0; j<L; j++){
	//				currentVal =  AMR_segment[2*j+1]*AMR_segment[2*j+1]; tmp3 =  AMR_segment[2*j+2]*AMR_segment[2*j+2] / sumWin;
	//				if (AMR_Sxx_PSD[j] < currentVal)	AMR_Sxx_PSD[j] = currentVal;
	//			}
	//		}else{
	//			for(j=0; j<L; j++){
	//				AMR_Sxx_PSD[j] =  AMR_segment[2*j+1]*AMR_segment[2*j+1]; tmp3 =  AMR_segment[2*j+2]*AMR_segment[2*j+2] / sumWin;
	//			}
	//		}
			// Spectrum trace is ClearWhite
			// evaluate magnitude of spectrum

			for(j=0; j<L; j++){
				//tmpf = AMR_SQR(AMR_IQSegment_PSD[2*j+1]) + AMR_SQR(AMR_IQSegment_PSD[2*j+2]); // slow
				Sxx  =  (AMR_IQSegment_PSD[2*j+1]*AMR_IQSegment_PSD[2*j+1] +  // e+1 order
						AMR_IQSegment_PSD[2*j+2]*AMR_IQSegment_PSD[2*j+2]);
				//AMR_Sxx_PSD[j] +=  tmpf; // /  sumWin;    // e-3 ~ e+4 order 
				AMR_Sxx_PSD[j] += Sxx;
			}

		} // end for(i=0; i<numSeg; i++) 

		//for(i=1; i<=(noverlap-1); i++){   // except DC
  //  		AMR_SWAP_FLOAT(AMR_Sxx_PSD[2*i+1], AMR_Sxx_PSD[L - (2*i-1)]);  // real
  //  		AMR_SWAP_FLOAT(AMR_Sxx_PSD[2*i+2], AMR_Sxx_PSD[L - (2*i-2)]);  // imag
		//}

		// get average spectrum
		for(i=0; i<L; i++) {
			AMR_Sxx_PSD[i] /= (numSeg*sumWin) ;	// e-5 order
			// convert to dB scale
			//AMR_Sxx_PSD[i] = (10.0)*log10(AMR_Sxx_PSD[i]);
		}

		fftshift(AMR_Sxx_PSD,outPSD,L);
	}else{

	}
	// for debugging

	#if VISUALSTUDIO
		fa=fopen("filename.txt","w");
		fprintf(fa,"data =[");
		for(i=0; i<L; i++) fprintf(fa,"%0.4f,...\n",(outPSD[i]));
		fprintf(fa,"]; \n");
		fprintf(fa,"figure; plot(data)");
		fclose(fa);		    
	#endif	

}

#if AMR_Cplus
	inline float AMRalgorithm::cm(float* inCm, int inCMLength, unsigned int nOrder, unsigned int qConj)
#else
	float cm(float* inCm, int inCMLength, unsigned int nOrder, unsigned int qConj)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  > Calculate the maximum absolute value of cyclic moments 
//
// \ INPUT ARGUMENTS
//	> inCm  : IQ data in a segment
//  > inCMLength : the length of segment
//  > nOrder : N-th power
//  > qConj  : the number of conjugate terms with n-th power
//
// \ OUTPUT ARGUMENTS
//  > maximum : maximum absolute value of cyclic moments 
//
// \ Example
//
// \ See also my_cyclic_cumulant_v3.m in MATLAB
//
// \ Author(s) : AWH
// 
// \ Reference
//   [1] Octavia A. Dobre, Ali Abdi, Yeheskel Bar-Ness, Wei Su,
//       "Cyclostationarity-Based Modulation Classification of Linear
//       Digital Modulations in Flat Fading Channels", Wireless Personal
//       Communications, 2009
//
//   [2] O.A. Dobre, Y. Bar-Ness and Wei Su, "Robust QAM modulation
//       classification algorithm using cyclic cumulants," in Wireless
//       Communications and Networking Conference, 2004. WCNC. 2004 IEEE,
//       pp. 745-748 Vol.2, 2004.
//
//   [3] "CLASSIFICATION OF DIGITAL MODULATION TYPES IN MULTIPATH ENVIRONMENTS", 2008
/******************************************************************************/

	int idx=0, L = inCMLength, tmpi=0;
	float constant = 0.0;
	float maximum=0.0;
	
	for(idx=0; idx<L; idx++){
		AMR_IQbufferCC[2*idx+1] = inCm[2*idx+1];
		AMR_IQbufferCC[2*idx+2] = inCm[2*idx+2];
	}

	nonlinearTransform(AMR_IQbufferCC, L, nOrder, qConj);

	four1(AMR_IQbufferCC, L, 1); // calculate FFT

    //for(idx=1; idx<=(Ld2-1); idx++){   // except DC
    //	AMR_SWAP_FLOAT(AMR_IQbufferCC[2*idx+1], AMR_IQbufferCC[L - (2*idx-1)]);  // real
    //	AMR_SWAP_FLOAT(AMR_IQbufferCC[2*idx+2], AMR_IQbufferCC[L - (2*idx-2)]);  // imag
    //}

	constant = (1.0f)/(float)L;
	for(idx=0 ; idx<L ; idx++){
		AMR_spectrumMagnitudeCC[idx] = constant*myCabs( AMR_IQbufferCC[2*idx+1],AMR_IQbufferCC[2*idx+2]);
		// sqrt, pow 등 미리 선언된 함수를 사용할때 데이터 타입에 맞는 함수를 쓰기 바람
		// DSP에서는 sqrt, Visual studio 에서는 sqrtl
	}

	// Find the largest elements in array
	// By assuming the frequency offset is order of Hz,
	// we can reduce the peak searching range
	// ex) fs = 140e6/256, L = 32768, resoultion = fs/L = 16.69 Hz,
	//     search range = fs/22 = 24.85 kHz  
	//     Therefore, peak searching range is 
	//     0 ~ 24.85kHz and fs-24.85kHz ~ fs  --> -24.85kHz ~ 24.85kHz
	tmpi = L/22;
	maximum = AMR_spectrumMagnitudeCC[0];
	for (idx = 1 ; idx<tmpi ; idx++){
		if(AMR_spectrumMagnitudeCC[idx] > maximum ){
			maximum = AMR_spectrumMagnitudeCC[idx];
		}
	}
	tmpi = 21*L/22;
	for (idx = tmpi ; idx<L ; idx++){
		if(AMR_spectrumMagnitudeCC[idx] > maximum ){
			maximum = AMR_spectrumMagnitudeCC[idx];
		}
	}

	maximum = fabs(maximum);

	return maximum;

	//#if VISUALSTUDIO
	//	fa=fopen("filename.txt","w");
	//	fprintf(fa,"data =[");
	//	for(idx=0; idx<L>>1; idx++) fprintf(fa,"%0.4f,...\n",(AMR_spectrumMagnitudeCC[idx]));
	//	fprintf(fa,"]; \n");
	//	fprintf(fa,"figure; plot(data)");
	//	fclose(fa);		    
	//#endif
}

#if AMR_Cplus
	inline void AMRalgorithm::filter(float* inB, float* inX,unsigned int LenB, unsigned int LenX, float* filterout, unsigned int realFlag)
#else
	void filter(float* inB, float* inX,unsigned int LenB, unsigned int LenX, float* filterout, unsigned int realFlag)
#endif
{

////////////////////////////////////////////////////////////////////
//
// Implementation smooth function in MATALB
// written by YuJun Won
// editted by WooHyun Ahn
//
// Usage
//
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
//   Y = FILTER(B,A,X) filters the data in vector X with the
//   filter described by vectors A and B to create the filtered
//   data Y.  The filter is a "Direct Form II Transposed"
//   implementation of the standard difference equation:
//
//   a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
//                         - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
//
//   We assume that coefficients a(1) = 1 and a(2),.. a(na+1) are all zero
//
//	 y(0) =  b(0)*x(0)
//   y(1) =  b(0)*x(1) + b(1)*x(0)
//   y(2) =  b(0)*x(2) + b(1)*x(1) + b(2)*x(0)
//   ....
//  
//	 See also filter.m in MATLAB
////////////////////////////////////////////////////////////////////

	float realSiftSum = 0.0f;
	float imagSiftSum = 0.0f;
	unsigned int i=0,j=0,supposLen = 1;      // length of superposition

	if (realFlag){  // real data
		for (i=0; i<LenX; i++)
		{
			for(j=0; j<supposLen; j++)	realSiftSum += inB[j]*inX[i-j];

			if(supposLen < LenB)	supposLen++;

			filterout[i] = realSiftSum;
			realSiftSum = 0.0f;
		}
	}
	else{  // complex data
		for (i=0; i<LenX; i++)
		{
			for(j=0; j<supposLen; j++){
				realSiftSum += inB[j]*inX[2*(i-j)+1];  //=(2*i+1)-2*j,  real
				imagSiftSum += inB[j]*inX[2*(i-j)+2];  //=(2*i+2)-2*j, imag
			}

			if(supposLen < LenB)	supposLen++;

			filterout[2*i+1] = realSiftSum;				// e-3 order
			filterout[2*i+2] = imagSiftSum;				// e-4 order
			 
			realSiftSum = 0.0f;
			imagSiftSum = 0.0f;
		}
	}
}

#if AMR_Cplus
	inline void AMRalgorithm::movingAverageFilter(float* MovInSeq, signed int MovInSeqLen, signed int MovWidth, float* MovOutSeq)
#else
	void movingAverageFilter(float* MovInSeq, signed int MovInSeqLen, signed int MovWidth, float* MovOutSeq)
#endif
{
/*********************************************************************************
 
 Inputs:
	MovInSeq[] : Array of real data points, size of AMR_NFFT/2-1.
	MovInSeqLen : The length of MovInSeq
	MovWidth : Filter length. It should be odd value

 Outputs:
	MovOutSeq[] : 

  See smooth.m in MATLAB

  Mostly used for coarse symbol rate estimation and spectrum estimation
************************************************************************************/
	signed int i;
	signed coutidx1=0, count =0;
	signed int loopLimit=0;
	signed int width = 0;

	if (MovWidth == 1){ // not processing
		for(i=0; i<MovInSeqLen ; i++)
		{
			MovOutSeq[i] = MovInSeq[i];
		}
	}else{

		width = MovWidth-1 + MovWidth%2; // force it to be odd.
		//int tmp= MovWidth;
		//MovWidth =  (MovWidth -1 + tmp%2);

		for(i=0;i<width;i++)	AMR_MovInB[i] = 1.0f/(float)width;

		filter(AMR_MovInB,MovInSeq,width,MovInSeqLen,AMR_Movfiltout,1); // real

	//  The default length for the moving average is 5.
	//  The first few elements of yy are given by
	//	yy(1) = y(1)
	//	yy(2) = (y(1) + y(2) + y(3))/3
	//	yy(3) = (y(1) + y(2) + y(3) + y(4) + y(5))/5
	//	yy(4) = (y(2) + y(3) + y(4) + y(5) + y(6))/5

		// cumulative sum
		AMR_cbegin1[0] = MovInSeq[0];
		for(i=1;i<width-2;i++)	AMR_cbegin1[i] = AMR_cbegin1[i-1]+MovInSeq[i];

		// because AMR_width is odd (AMR_width-1)/2 is even

		loopLimit = (width-1)/2;
		for(i=0;i<loopLimit;i++){
			//tmpi2 = i*2+1;
			AMR_cbegin2[i] = AMR_cbegin1[i*2]/(float)(i*2+1) ;
		}

	//  The last few elements of yy are given by
	//	yy(end) = y(end)
	//	yy(end-1) = (y(end) + y(end-1) + y(end-2))/3
	//	yy(end-2) = (y(end) + y(end-1) + y(end-2) + y(end-3) + y(end-4))/5
	//	yy(end-3) = (y(end-1) + y(end-2) + y(end-3) + y(end-4) + y(end-5))/5

		AMR_cend1[0] = MovInSeq[MovInSeqLen-1]; //sample of input(end)
		//for(i=1;i<AMR_spanOfMAfilter-3;i++)
		//MovInSeqLen

		/* cend = cumsum(y(n:-1:n-width+3));*/
		// cumulative sum

		for(i=MovInSeqLen-1;i>MovInSeqLen-width+1;i--)	{
			if(count == 0 )	AMR_cend1[count] = MovInSeq[i];
			else 			AMR_cend1[count] = AMR_cend1[count-1] + MovInSeq[i];
			count += 1;
		}

		/* cend = cend(end:-2:1)./(width-2:-2:1)';*/
		count = 0;
		for(i=width-2; i>0; i=i-2) {
			AMR_cend2[count] = AMR_cend1[i-1]/(float)i;
			count++;
		}

		// below codes are equal to "c = [cbegin;c(width:end);cend];" in MATLAB
		for(i=0;i<loopLimit;i++)	MovOutSeq[i] = AMR_cbegin2[i];

		coutidx1 = width-1;
		//for(i=(AMR_width-1)/2;i<(AMR_width-1)/2+(MovInSeqLen-AMR_width+1);i++)
		for(i=loopLimit ; i<(MovInSeqLen-loopLimit) ; i++)	{
			MovOutSeq[i] = AMR_Movfiltout[coutidx1];
			coutidx1++;
		}

		coutidx1=0;
		//for(i=(AMR_width-1)/2+(MovInSeqLen-AMR_width+1);i<MovInSeqLen ;i++)
		for(i=(MovInSeqLen-loopLimit); i<MovInSeqLen ; i++)
		{
			MovOutSeq[i] = AMR_cend2[coutidx1];
			coutidx1++;
		}
	}
}

#if AMR_Cplus	
	inline float AMRalgorithm::coarseSNRestimation(signed int snrFFTLength)
#else
	float coarseSNRestimation(signed int snrFFTLength)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - Estimate SNR
//
// \ INPUT ARGUMENTS
//	> snrFFTLength  - The length of FFT
//  
// \ OUTPUT ARGUMENTS
//	> AMR_coareseSNR_dB - Estimated SNR
//  > AMR_coareseSNRIdx - Index of estimated SNR
//
// \ Author(s) : AWH
// 
/******************************************************************************/
	signed int i     = 0;
	signed int L     = snrFFTLength;
	signed int denominator = 0;
	float sigPower = 0.0, noisePower =0.0;
	float coareseSNR_dB = 0.0;
#if AMR_DEBUG
	printf(" \r\n #SNRE     ");
#endif

	// average signal power
	for(i=AMR_BWSampleIdx[0]; i<=AMR_BWSampleIdx[1]; i++)	sigPower += AMR_powerSpectrum[i];	
	denominator = AMR_BWSampleIdx[1] - AMR_BWSampleIdx[0] + 1;
	sigPower /= denominator;  // Watt

	// average noise power
	for(i=0 ; i<AMR_BWSampleIdx[0]; i++)		noisePower += AMR_powerSpectrum[i];
	for(i=AMR_BWSampleIdx[1]+1 ; i< L ; i++)	noisePower += AMR_powerSpectrum[i];	
	denominator = AMR_BWSampleIdx[0] + L - AMR_BWSampleIdx[1] - 1;
	if (denominator == 0) denominator = 1;
	noisePower /= (float) denominator;  // Watt

	// 24.77 in DSP	   14.12.08
	// 24.67 in MATLAB
	coareseSNR_dB = 10.0f*log10(sigPower / noisePower);
    
	// Calculate noise power with user setting SNR
	sigPower = 10.0f*log10(sigPower); // dB
	AMR_noisePower = sigPower - AMR_SNR_dB; // dB
	AMR_noisePower = pow(10.0f,AMR_noisePower/10.0f);
	//tmpi2 = sizeof(AMR_SNRTable) / sizeof(AMR_SNRTable[0]);
	//	
	//for(i=0 ; i < tmpi2-1 ; i++){
	//	tmpf2 = ( AMR_SNRTable[i] + AMR_SNRTable[i+1] ) / 2.0f;
	//	if (AMR_coareseSNR_dB < tmpf2){
	//		AMR_coareseSNRIdx = i; break;
	//	}
	//}
	//tmpi = (int) AMR_coareseSNR_dB;   // (float)->(uint)
 //   // force i to even number
 //   if (tmpi%2 == 1) tmpi = tmpi + 1;
 //   // exception handling
 //   if (tmpi > AMR_SNRTable[tmpi2-1]) tmpi = AMR_SNRTable[tmpi2-1];
 //   else if (tmpi < AMR_SNRTable[0] )tmpi =  AMR_SNRTable[0];

 //   for(i=0; i<tmpi2; i++){
 //   	if(AMR_SNRTable[i] == tmpi) AMR_coareseSNRIdx = i;
 //   }
	

#if AMR_DEBUG
			printf(" SNR : %d[dB] \r\n ",(int)coareseSNR_dB);
#endif
	return coareseSNR_dB;
}

#if AMR_Cplus	
	inline float AMRalgorithm::coarseBandWidthEstimation(float* BWEdata,unsigned int BWE_FFTLength)
#else
	float coarseBandWidthEstimation(float* BWEdata,unsigned int BWE_FFTLength)
#endif

{
/********************************************************************************/
// \ DESCRIPTION
//  > Estimate coarse bandwidth based on histogram
//
// \ INPUT ARGUMENTS
//  > BWEdata : Welch's power spectrum in dB scale
//  > BWE_FFTLength : should be power of 2
//
// \ OUTPUT ARGUMENTS
//	> AMR_BWSampleIdx[1]  : upper frequency of bandwidth
//  > AMR_BWSampleIdx[0]  : lower frequency of bandwidth
//  > AMR_samPeriod       : sampling period
//  > AMR_coarseBandWidth
//
// \ Example
//
// \ See also histogram
//
// \ Author(s) : AWH
// 
// \ Reference
//  > ITU-R SM.443-4 recommended
/******************************************************************************/

	float yMax = 0.0, yMin = 0.0;
	float threshold = 0.0f;
	float minMaxGap = 0.0f;//, offset = 0.0f;
//	float minCount = 0.0;
//	float leftSigma3 = 0.0, rightSigma3 = 0.0;
	float noiseFloor      = 0.0f, signalFloor=0.0f;
//	float levelB2NoiseNSig   		   = 0.0f;
	float lineSpectrumHeight		   = 0.0f;
	float lineSpecMagThreshold		   = 0.0f;
	float subLineSpecMagThreshold	   = 0.0f;
//	float normFactor				   = 0.0f;
	float minSpect					   = 0.0f;
  //	float leftBinCenter				   = 0.0f;
//	float rightBinCenter 			   = 0.0f;
	float coarseBW 					   = 0.0f;
  //	float binWidth 					   = 0.0f;
	float firstMaxCnt = 0.0f, firstMinCnt = 0.0f, secondMaxCnt = 0.0f;
	float binWidth = 0.0f;
	float bandwidthConstant = 0.0f;

	unsigned int numOfBins            = 0;
	unsigned int isSingleLineSpectrum = 0;
	//unsigned int isCW = 0;;
	unsigned int axisLowerIdx = 0, axisUpperIdx = 0;

	unsigned int i=0,k=0;
	unsigned int freqCount = 0;
//	int thresholdIdx = 0;

	//freq vector
	unsigned int L = BWE_FFTLength;
	//unsigned int count = 0;	

	int noiseFloorIdx = 0;
	int noiseFloorCount = 0;
	int signalFloorCount = 0; //, signalFloorIdx=0;
	int yMaxIdx = 0;
	int firstMaxCntIdx = 0;

	signed int leftIdx =0, rightIdx=0;
	//unsigned int bwFlag=0;
	//unsigned int centerIdx = 0;
	unsigned int lowerFlag = 0, upperFlag = 0;
	unsigned int isLocalMaximum =0;
	signed int   bwIndex = 0; 
	signed int   idxOffset = 0; // 07-Sep-2015 for CW
	int idxVal = 0;
	int ampLeftIdx = 0, ampRightIdx = 0;
	int iAmpLen = 0;
	int offset = 0;
	// Frequency smoothing
	//movingAverageFilter(BWEdata, L, 1, AMR_powerSpectrum_dB);
#if AMR_DEBUG
	printf(" \r\n #BWE  \r\n");
#endif


	bandwidthConstant = 0.7f;

	// Convert linear scale to log scale
	for(i=0; i<L; i++)	AMR_powerSpectrum_dB[i] = (10.0f)*log10(BWEdata[i]);
#if VISUALSTUDIO
		fa=fopen("filename.txt","w");
		fprintf(fa,"data =[");
		for(i=0; i<L; i++) fprintf(fa,"%0.4f,...\n",(AMR_powerSpectrum_dB[i]));
		fprintf(fa,"]; \n");
		fprintf(fa,"figure; plot(data)");
		fclose(fa);			
#endif
	// Find min. and max. Duplicated code in histogram
	yMax = AMR_powerSpectrum_dB[0];
	for (i=1 ; i < L ; i+= 1){
		if(AMR_powerSpectrum_dB[i] > yMax ){
			yMaxIdx = i;
			yMax = AMR_powerSpectrum_dB[i];
		}
	}

	yMin = yMax;
	for (i=0 ; i < L ; i+= 1){
		if(AMR_powerSpectrum_dB[i] < yMin ){
			yMin = AMR_powerSpectrum_dB[i];
		}
	}

	// Min-max normalization     ///check--> ok
	//normFactor = 100.0f / (yMax-yMin);
	//for (i=0 ; i < L ; i+= 1)	AMR_powerSpectrum_dB[i] = (AMR_powerSpectrum_dB[i]-yMin) * normFactor;
	//yMax = 100.0f; yMin = 0;

	numOfBins = AMR_MAX_NUMBINS;
	minMaxGap = yMax-yMin;
	binWidth = minMaxGap/numOfBins; // dB
	//numOfBins = minMaxGap / binWidth;

	if (numOfBins > AMR_MAX_NUMBINS) numOfBins = AMR_MAX_NUMBINS;
#if AMR_DEBUG
	printf(" #OfBins : %u \r\n",numOfBins);
#endif	
	
	//binWidth = minMaxGap / numOfBins;
	// yMax      yMin
	// 31.3     -30.03	in MATLAB
	// 36.11	-26.71	in DSP
	histogram(AMR_powerSpectrum_dB, L, numOfBins, yMax, yMin);
	
	// Find noise floor
	noiseFloorCount = AMR_binCount[0];
	for (i=1; i<numOfBins; i++){
		if (AMR_binCount[i] > noiseFloorCount){
			noiseFloorCount = AMR_binCount[i];
			noiseFloor      = AMR_binCenter[i];
			noiseFloorIdx   = i;
		}
	}

	signalFloorCount = 0.0f; isLocalMaximum = 0;
	threshold = yMax*0.8f + yMin*0.2f;
	offset = (int)(6.0f/binWidth);
	if (noiseFloor > threshold){
		// In case of that power of signal is larger than power of noise 
	
		for(i=(noiseFloorIdx-4); i>0; i--){
			leftIdx = i-offset;  
			if(leftIdx<0) leftIdx=0;

            rightIdx = i+offset; 
			if(rightIdx>(numOfBins-1)) rightIdx = numOfBins-1;

			// Check current i-th samples is local maximum in the range of i-3 to i+3
			// See histogram of spectrum
			if(AMR_binCount[i] > AMR_MINIMUM_HIST_COUNT){			
				isLocalMaximum = 1;
				for(k=leftIdx; k<=rightIdx; k++){
					if(AMR_binCount[k] >  AMR_binCount[i]){
						isLocalMaximum = 0;
					}else{

					}
				}

				if(isLocalMaximum){
					if(AMR_binCount[i] > signalFloorCount){
						signalFloorCount = AMR_binCount[i];
						signalFloor      = AMR_binCenter[i];
						break;
					}
				}
			}
		}	
	}else{
		// In case of that power of signal is smaller than power of noise 
		for(i=numOfBins; i>=(noiseFloorIdx+4); i--){
			leftIdx = i-offset;  
			if(leftIdx<0) leftIdx=0;

            rightIdx = i+offset; 
			if(rightIdx>(numOfBins-1)) rightIdx = numOfBins-1;

			// Check current i-th samples is local maximum in the range of i-3 to i+3
			// See histogram of spectrum
			if(AMR_binCount[i] > AMR_MINIMUM_HIST_COUNT){		// 13 is determined by experience
				isLocalMaximum = 1;
				for(k=leftIdx; k<=rightIdx; k++){
					if(AMR_binCount[k] >  AMR_binCount[i]){
						isLocalMaximum = 0;
					}else{

					}
				}
				// If it is local maximum, assign new value to variable
				if(isLocalMaximum){
					if(AMR_binCount[i] > signalFloorCount){
						signalFloorCount = AMR_binCount[i];
						signalFloor      = AMR_binCenter[i];
						break;
					}
				}
			}else{

			}
		}
	}

	if(isLocalMaximum == 0){//in case of there is no signal floor
		// is FM?
		signalFloor		  = yMax;
		signalFloorCount  = AMR_binCount[numOfBins-1];
		bandwidthConstant = 0.3f;
	}

#if AMR_DEBUG
	printf("\r Noise Floor : %d, Signal Floor : %d \r\n",(int)noiseFloor,(int)signalFloor);
#endif

	// Exception handling - When BW is similar to fs/2
	//if (levelB2NoiseNSig > (0.9f*minMaxGap) ){ // edit 15.05.21
	//	#if AMR_DEBUG
	//	//	printf("\r\n  (level for Noise and Signal) > %d   \r\n",int(levelB2NoiseNSig));
	//	#endif
	//	// When count of signalFloor is greater than that of noiseFloor,
	//	//  2 times of bandwidth is closed to the sampling frequency
	//	signalFloor = noiseFloor;
	//	levelB2NoiseNSig = signalFloor - (3.0f*binWidth);

	//	for (i=0; i<numOfBins; i++){
	//		if (AMR_binCenter[i] < levelB2NoiseNSig){
	//			if (count == 0){    // to save the first sample index higher than threshold
	//				thresholdIdx = i;
	//				count = thresholdIdx;
	//			}
	//			signalFloorCount = AMR_binCount[count];
	//			signalFloorIdx = count;
	//			if(AMR_binCount[count] > signalFloorCount ){  // Find local maximum
	//				signalFloorCount = AMR_binCount[count];
	//				signalFloorIdx = count;
	//			}
	//			count += 1;
	//		}
	//	}
	//	noiseFloor = AMR_binCenter[signalFloorIdx];
	//}else{
	//
	//	for (i=0; i<numOfBins; i++){
	//		if (AMR_binCenter[i] > levelB2NoiseNSig){
	//			if (count == 0){    // to save the first sample index higher than threshold
	//				thresholdIdx = i;
	//				count = thresholdIdx;
	//			}
	//
	//			if(AMR_binCount[count] > signalFloorCount ){  // Find local maximum
	//				signalFloorCount = AMR_binCount[count];
	//				signalFloorIdx = count;
	//			}
	//			count += 1;
	//		}
	//	}
	//	signalFloor = AMR_binCenter[signalFloorIdx];
	//}
	

	// exception handling
	idxOffset = 5;
	if((yMaxIdx-idxOffset)<0){  leftIdx = 0; } 
	else {leftIdx = (unsigned int) (yMaxIdx-idxOffset);}
	if((yMaxIdx+idxOffset)>L-1) { rightIdx = L-1;}
	else {rightIdx = (unsigned int) (yMaxIdx+idxOffset);}

	// find a minimum spectrum in the range of +-idxOffset peak
	minSpect = AMR_powerSpectrum_dB[leftIdx];
	for(i = leftIdx; i <= rightIdx; i++){
		if(AMR_powerSpectrum_dB[i] < minSpect)
			minSpect = AMR_powerSpectrum_dB[i];
	}
		
	lineSpectrumHeight = yMax - minSpect;

	// check existence of another carrier components
    isSingleLineSpectrum = 1;
	// Threshold is modified at 09-Sep-2015
	subLineSpecMagThreshold = yMax*0.7; // 18-Sep-2015,  - (lineSpectrumHeight/2.0f);
    for(i=1; i<=leftIdx; i++){
        if( AMR_powerSpectrum_dB[i] >  subLineSpecMagThreshold) isSingleLineSpectrum = 0;
	}
	for(i=L; i>=rightIdx; i--){
		if( AMR_powerSpectrum_dB[i] >  subLineSpecMagThreshold) isSingleLineSpectrum = 0;
	}

	// For FM vs FSK classification
	//leftSigma3  = signalFloor - 20.0;	// 15 -> 20 at 150226	
	//rightSigma3 = signalFloor + 20.0;
	//minCount = signalFloorCount;
	//// find a minimum count in the valley
	//for (i=0; i<numOfBins; i++){
	//	if ( (leftSigma3 <= AMR_binCenter[i]) & (AMR_binCenter[i] < rightSigma3) ){
	//		printf("%d",i);
	//		if ( AMR_binCount[i] < minCount ){
	//			minCount = AMR_binCount[i];
	//		}
	//	}
	//}
	//if (minCount ==0) minCount =1.0f;
	//// For FM vs M-FSK classification
	//AMR_sigFloorRatio = (float) signalFloorCount / (float) minCount;

	// Edit decision tree to prevent the mis-classification of AM under SNR > 25, 150410
	// AM signal count on the signal floor is greater than AMR_SIGNALFLOOR_THRESHOLD under SNR > 25, 150410

	// Depth of AM increase, Height of line spectrum decrease
	lineSpecMagThreshold = AMR_CARRIER_LINE_SPECTRUM_THRESHOLD;
#if AMR_DEBUG
	printf("\r (Spectrum Height (%d) > %d ) ? AM,CW,2ASK : others \r\n",(int)lineSpectrumHeight,(int)lineSpecMagThreshold);
	printf("\r (Is line spectrum single? (%u) == 1) ? AM,CW,2ASK : others \r\n",isSingleLineSpectrum);
#endif
	
	if( (lineSpectrumHeight > lineSpecMagThreshold) && isSingleLineSpectrum)
	{
#if AMR_DEBUG
	printf("\r Single line \r\n");
#endif
		 // Carrier exists, AM, 2ASK, CW
		AMR_existsCarrierFlag = 1;	
		AMR_linearFlag  = 1;

		// find maximum elements in the range
		//offset = 0.07f*minMaxGap;			//0.14f ---> 0.07f at 07-Sep-2015 by AWH
		//leftBinCenter  = signalFloor - offset; //14.0f;
		//rightBinCenter = signalFloor + offset; //14.0f;

		//for (i=0; i<numOfBins; i++){
		//	if ( (leftBinCenter <= AMR_binCenter[i]) & (AMR_binCenter[i] <= rightBinCenter) ){
		//		if (AMR_binCount[i] > signalFloorCount) notLocalMaximum = 1;
		//	}
		//}
		/*
		#if AMR_DEBUG
			printf("\r\n (signalFloorCount(%u) < %u) ? AM : 2ASK \r\n",
					signalFloorCount,(unsigned int)AMR_SIGNALFLOOR_THRESHOLD);
		#endif
			
		#if AMR_DEBUG
			printf("\r\n (notLocalMaximum(%d) == 1) ? AM : 2ASK \r\n",notLocalMaximum);
		#endif
		*/

		// Set the length of iAmp, 
		if (AMR_LEN > AMR_symbolRate_SamplesPerSegment){
			iAmpLen = AMR_symbolRate_SamplesPerSegment;
		}else{
			iAmpLen = AMR_LEN;
		}

		//for(j=segStart; j<segEnd; j++){	
		for(i=0; i< iAmpLen; i++){	// 5 is offset to prevent invalid samples
			 // instantaneous amplitude
			AMR_iAmplitude[i] = myCabs(AMR_inputData[2*i+1],AMR_inputData[2*i+2]);    
		}
		for(i=0; i<iAmpLen-10; i++){ // NFFT is always larger than 70
			//AMR_IQsegmentSym[2*j+1] = AMR_IQsegmentSym[2*(j+70)+1];
			AMR_iAmplitude[i] = AMR_iAmplitude[i+10];
		}
		for(i=0; i<10; i++){    
			//AMR_IQsegmentSym[2*(L-70+j)+1] = AMR_IQsegmentSym[2*(L-1)+1];
			AMR_iAmplitude[iAmpLen-10+i] = AMR_iAmplitude[iAmpLen-1];   
		}
			

		movingAverageFilter(AMR_iAmplitude, iAmpLen, 10, AMR_filteredSamples);

		// Find min. and max. Duplicated code in histogram
		firstMaxCnt = AMR_filteredSamples[0];
		firstMaxCntIdx = 0;
		for (i=1 ; i < iAmpLen ; i++){
			if(AMR_filteredSamples[i] > firstMaxCnt ){
				firstMaxCntIdx = i;
				firstMaxCnt = AMR_filteredSamples[i];
			}
		}

		firstMinCnt = firstMaxCnt;
		for (i=0 ; i < iAmpLen ; i++){
			if(AMR_filteredSamples[i] < firstMinCnt ){
				firstMinCnt = AMR_filteredSamples[i];
			}
		}

		//binWidth = minMaxGap / numOfBins;
		// yMax      yMin
		// 31.3     -30.03	in MATLAB
		// 36.11	-26.71	in DSP

		histogram(AMR_filteredSamples, iAmpLen, numOfBins, firstMaxCnt, firstMinCnt);

		firstMaxCnt = AMR_binCount[0];
		for (i=1 ; i < numOfBins ; i++){
			if(AMR_binCount[i] > firstMaxCnt ){
				firstMaxCntIdx = i;
				firstMaxCnt = AMR_binCount[i];
			}
		}

		if ((firstMaxCntIdx-1) > (numOfBins-firstMaxCntIdx)){
           // Right normal dist.
          // distWidth = numOfBins-firstMaxCntIdx;
		   // idxVal = (firstMaxCntIdx-distWidth);

           idxVal = (firstMaxCntIdx<<1)-numOfBins;
           
		   if (idxVal < 0){ idxVal = 0; }
		   else{}
		   ampLeftIdx = 0;  ampRightIdx = idxVal;

		}else{
           // Left normal dist.
           //distWidth = firstMaxCntIdx-1;
//           idxVal = firstMaxCntIdx+distWidth;

           idxVal = firstMaxCntIdx<<1;
           
		   if (idxVal > numOfBins){ idxVal = numOfBins; }
		   else{}
		   ampLeftIdx = idxVal;  ampRightIdx = numOfBins;
		}
		//AMR_binCount
		//AMR_binCenter
		secondMaxCnt = AMR_binCount[ampLeftIdx];
		//secondMaxCntIdx = ampLeftIdx;
		for (i=ampLeftIdx+1; i<ampRightIdx; i++){
			if(AMR_binCount[i] > secondMaxCnt ){
		  //		secondMaxCntIdx = i;
				secondMaxCnt = AMR_binCount[i];
			}
		}

		//if (firstMaxCntIdx > secondMaxCntIdx){
		//	peakDist = (firstMaxCntIdx-secondMaxCntIdx);
		//}else{	
		//	peakDist = (secondMaxCntIdx-firstMaxCntIdx);
		//}
		//distWidth = distWidth-2; // apply offset

		#if AMR_DEBUG
			printf("\r firstMaxCnt : %d, secondMaxCnt : %d \r\n", (int)firstMaxCnt, (int)secondMaxCnt);
			//printf("\r peakDist : %d, distWidth : %d \r\n",peakDist, distWidth);
		#endif

        if((0.35f*firstMaxCnt < secondMaxCnt )){ // 0.5f --> 0.35f at 21-Sep-2015
			// ASK
			AMR_digitalFlag = 1;
			AMR_modNum	  	= 2;
			AMR_modOrder 	= 2;
			

			AMR_bwThreshold_dB = signalFloor*bandwidthConstant + 
				noiseFloor*(1.0f-bandwidthConstant); //noiseFloor + (( signalFloor - noiseFloor ) * 0.45f); 
		}else{
			// AM
				AMR_digitalFlag = 0;  // analog
				AMR_modNum   	= 7;
				AMR_modOrder 	= 0;
//			// Check mis-classification
//			AMR_LEN = 5000;
//			// Calculate the instantaneous phase
//			k=0;
//			for(i=0 ; i<AMR_LEN ; i++){
//				//AMR_iPhase[k] = 0.0f;
//				AMR_iPhase[k] = atan2(AMR_inputData[2*i+2],AMR_inputData[2*i+1]); // wraped instantaneous phase
//				AMR_ck[k] = 0.0f;
//				k++;
//
//			}
//
//			// make range from -PI to +PI (=unwrap)
//			AMR_ck[0] = 0.0f;
//			for(i=1; i<AMR_LEN ; i++){
//				iFreq = AMR_iPhase[i] - AMR_iPhase[i-1];
//				if( iFreq > AMR_PI)
//					AMR_ck[i] = AMR_ck[i-1] - AMR_TWOPI;
//				else if ( iFreq < (-1.0f)*AMR_PI)
//					AMR_ck[i] = AMR_ck[i-1] + AMR_TWOPI;
//				else	AMR_ck[i] = AMR_ck[i-1];
//
//				AMR_iPhase[i-1] += AMR_ck[i-1]; // unwraped instantaneous phase
//			}
//			AMR_iPhase[AMR_LEN-1] += AMR_ck[AMR_LEN-1]; // unwraped instantaneous phase
//
//				// Reuse array for linear curve fitting
//			for(i=0; i<AMR_LEN; i++)	AMR_spectrumBuffer[i] = (float)i;
//			chi2 = 0.0f;
//			fit(AMR_spectrumBuffer, AMR_iPhase, AMR_LEN, &a, &b, &siga, &sigb, &chi2);
//#if VISUALSTUDIO
//		fa=fopen("filename.txt","w");
//		fprintf(fa,"data =[");
//		for(i=0; i<L; i++) fprintf(fa,"%0.4f,...\n",(AMR_iPhase[i]));
//		fprintf(fa,"]; \n");
//		fprintf(fa,"figure; plot(data)");
//		fclose(fa);			
//#endif
//			if( chi2 >= 0.5f ){	// refer to AMvsFM.m in MATLAB
//				// FM
//#if AMR_DEBUG
//				printf("\r instant Freq. Check AM-->FM, chi %d \r\n",(int)(chi2*10.0f));
//#endif				
//
//
//				//AMR_digitalFlag = 0;  // analog
//				//AMR_modNum   	= 0;
//				//AMR_modOrder 	= 0;
//			}

		}

		// Use instantaneous amplitude to classify AM vs CW

		// Check whether it is CW written at 20-August-2015
		//isCW = 0;
		//subLineSpecMagThreshold = 0.45f * offset; //.0f; //

		//for(i=1; i<=leftIdx; i++){
		//	if( AMR_powerSpectrum_dB[i] >  subLineSpecMagThreshold) isCW = 1;
		//}
		//for(i=L; i>=rightIdx; i--){
		//	if( AMR_powerSpectrum_dB[i] >  subLineSpecMagThreshold) isCW = 1;
		//}
		//#if AMR_DEBUG
		//	printf("\r\n (isCW(%d) == 1) ? CW : AM \r\n",isCW);
		//#endif
		//if(isCW){
		//	AMR_digitalFlag = 0; // analog
		//	AMR_modNum   = 7;    ////////////// CW, but assign as unknown
		//	AMR_modOrder = 0;
		//	coarseBW     = 0;
		//}
		AMR_maxSpecMagIdx = yMaxIdx; // for carrier offset estimation
	}
	else
	{
		// non-Carrier
		AMR_existsCarrierFlag = 0;
		//AMR_digitalFlag = 1;
		AMR_modNum	  	= 0;		//unknown
		AMR_bwThreshold_dB = signalFloor*bandwidthConstant + 
				noiseFloor*(1.0f-bandwidthConstant);
#if AMR_DEBUG
	printf("\r Multiple line \r\n");
#endif
	}

	// Calculate bandwidth
	if(AMR_modNum == 7){ // is AM?
		idxVal = (int)( 7500.0f/ AMR_freqInterval);
	
		AMR_BWSampleIdx[0] = (L>>1)-idxVal;
		AMR_BWSampleIdx[1] = (L>>1)+idxVal;

		coarseBW = 15000.0f; // fixed
	}
	else 
	{
		// Store the indices of samples in bandwidth
		for(i= 0; i< L; i++){
			if ( AMR_powerSpectrum_dB[i] > AMR_bwThreshold_dB ){
				AMR_IndexBuffer_PSD[freqCount] = i;
				freqCount++;
				if(freqCount == 1) yMax = AMR_powerSpectrum_dB[i];
				else{
					if(AMR_powerSpectrum_dB[i] > yMax){
						yMax = AMR_powerSpectrum_dB[i];
						yMaxIdx = freqCount;
					}
				}
			}
		}

		// toward left edge, N-(N-1)
		for(i=0; i<yMaxIdx; i++) {
			if (( AMR_IndexBuffer_PSD[yMaxIdx-i]- AMR_IndexBuffer_PSD[yMaxIdx-i-1]) > AMR_minBinInterval){
				axisLowerIdx = AMR_IndexBuffer_PSD[yMaxIdx-i];  
				lowerFlag=1; break;
			}
		}

		// toward right edge, (N+1)-N
		for(i=yMaxIdx; i<freqCount-1; i++) {
			if ( (AMR_IndexBuffer_PSD[i+1] - AMR_IndexBuffer_PSD[i]) > AMR_minBinInterval){
				axisUpperIdx = AMR_IndexBuffer_PSD[i];     
				upperFlag=1; break;
			}
		}
		if (lowerFlag == 0) axisLowerIdx = AMR_IndexBuffer_PSD[0];
		if (upperFlag == 0) axisUpperIdx = AMR_IndexBuffer_PSD[freqCount-1];

		// save global variable
		// 1948, 2103 in DSP   at 141208
		// 1949, 2104 in MATLAB

		// 1620, 2016 in DSP   at 150109
		// 1609, 2017 in MATLAB
		AMR_BWSampleIdx[0] = (signed int)axisLowerIdx;
		AMR_BWSampleIdx[1] = (signed int)axisUpperIdx;

		// 27237 in MATLAB
		// 26436 in DSP
		bwIndex = (AMR_freqVector[AMR_BWSampleIdx[1]] - AMR_freqVector[AMR_BWSampleIdx[0]]);
		if (bwIndex < 0)
		{
			bwIndex = 0;
			#if AMR_DEBUG
				printf("\r Error in bandwidth estimation : negative value \r\n ");
			#endif
		}
	
		coarseBW = (float)bwIndex*AMR_freqInterval;
	}

	AMR_coareseSNR_dB = coarseSNRestimation(AMR_NFFT);
	AMR_noiseFloor_dB = noiseFloor;
	// SNR estimation
#if AMR_DEBUG
	
	printf("\r Bandwidth %d[Hz] \r\n ",(int)coarseBW);
#endif

	return coarseBW;

// method 1, ITU-R SM.443-4 recommended, select all samples higher than threshold
//	for (i = 0 ; i < BWEsegmentLength ; i++){
//		if (AMR_powerSpectrum_dB[i] >= AMR_bwThreshold){
//			axisLowerIdx = i; break;
//		}
//	}
//
//	for (i = 1 ; i <= BWEsegmentLength ; i++){
//		if (AMR_powerSpectrum_dB[BWEsegmentLength-i] >= AMR_bwThreshold){
//			axisUpperIdx = BWEsegmentLength-i; break;
//		}
//	}
}

#if AMR_Cplus	
	inline signed int AMRalgorithm::coarsefreqOffsetEstimation(signed int cfoDataLen,signed int cfoFFTlength)
#else
	signed int coarsefreqOffsetEstimation(signed int cfoDataLen,signed int cfoFFTlength)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  > Estimate coarse frequency offset estimation 
//
// \ INPUT ARGUMENTS
//  > cfoDataLen : The length of IQ data
//  > cfoFFTlength : The length of FFT. It should be power of 2
//
// \ OUTPUT ARGUMENTS
//  > shiftIdx : Index of frequency axis corresponding to offset
//	> AMR_coarseFreqOffset : frequency offset
//
// \ Example
//
// \ See also histogram
//
// \ Author(s) : AWH
// 
// \ Reference
/******************************************************************************/


	float estimatedFreqOffset=0.0;
	float cosArg = 0.0;
	float zc  = 0.0;
	float zs  = 0.0;
	float realPart = 0.0;
	float twoPI = AMR_TWOPI;

	signed int offsetIdx = 0;
	signed int i=0;
	signed int L;
	signed int Ld2;
	signed int shiftIdx = 0;

#if AMR_DEBUG
	printf("\r\n #CFO Est. \r\n");
#endif

	L 	= cfoFFTlength;
	Ld2 = L >> 1;

	if (AMR_existsCarrierFlag){ // AM, 2ASK, CW
		// Based on maximum point
		// Find min. and max. Duplicated code in histogram

		estimatedFreqOffset =  (float)AMR_freqVector[AMR_maxSpecMagIdx]*AMR_freqInterval; // AMR_maxSpecMagIdx from bandwidth estimation
		offsetIdx = AMR_maxSpecMagIdx;
	}else{
		// Based on middle point of bandwidth
		offsetIdx = ( AMR_BWSampleIdx[0] + AMR_BWSampleIdx[1] ) / 2;  // float -> int
		estimatedFreqOffset = (float)AMR_freqVector[offsetIdx]*AMR_freqInterval; // -15354.16
	}
	// sample index corresponding to frequency offset
	shiftIdx = offsetIdx - Ld2;

	// shift banwidth sample index
	// 1850 _ 2246 in DSP    140109
	// 1844 _ 2252 in MATLAB
	// for(i=0; i<2; i++) AMR_BWSampleIdx[i] -= shiftIdx; //  not neccessary in C, but need for MATLAB 

	// freqOffsetCompensation
	for(i=0; i<cfoDataLen; i++){
		cosArg = -i*twoPI*shiftIdx / L;   // e-3order,   +- 0.002 error compared to MATLAB
		zc	   = cos(cosArg);		zs = sin(cosArg);
		realPart = (AMR_inputData[2*i+1] * zc) - (AMR_inputData[2*i+2] * zs);		// real
		AMR_inputData[2*i+2] = (AMR_inputData[2*i+1] * zs) + (AMR_inputData[2*i+2] * zc);      // imag
		AMR_inputData[2*i+1] = realPart;
	}

	// -15354.16 in DSP
	// -15754.699 in MATLAB
	AMR_coarseFreqOffset = estimatedFreqOffset;
	AMR_BWSampleIdx[0] -= shiftIdx;
	AMR_BWSampleIdx[1] -= shiftIdx;

	// estimate coarse freq. offset
	// -1485 in DSP    at 14.12.08
	// -1502.037 in MATLAB

#if AMR_DEBUG
	printf("\r CFO %d[Hz] \r\n ",(int)AMR_coarseFreqOffset);
#endif
	return shiftIdx;
}

#if AMR_Cplus	
	inline void AMRalgorithm::spectrumFeatureExt(float* ifeIQdata, unsigned int ifeIQDataLen, unsigned int ifeNFFT) // only for digital signal
#else
	void spectrumFeatureExt(float* ifeIQdata, unsigned int ifeIQDataLen, unsigned int ifeNFFT) // only for digital signal
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - Calculate instantaneous amplitude, phase, and frequency
//
// \ INPUT ARGUMENTS
//	> ifeIQdata			 - IQ data 
//	> ifeIQDataLen		 - The length of IQ data 
//	> ifeNFFT			 - The length of FFT
//  
// \ OUTPUT ARGUMENTS
//	> 
//
// \ See also 
//
// \ Author(s) : AWH
// 
// \ Reference
//
/******************************************************************************/
	unsigned int i=0, j=0, k=0, m=0;
	unsigned int segStart=0, segEnd=0, outLoop = 0;
	unsigned int T	 = ifeIQDataLen;

	unsigned int numSeg	= 0; // T / L; 			// the number of segments, float -> int
	unsigned int NFFT   = ifeNFFT;
	unsigned int iLen   = AMR_instantPhaseObservationInterval;
	int leniFreq = 0;

	float a=0.0f, b=0.0f, siga=0.0f, sigb=0.0f, chi2=0.0f; // for linear curve fitting
	float iFreq = 0.0f;
	//for instantaneous phase
	float PI		  = AMR_PI;
	float TWOPI 	  = AMR_TWOPI;
	// Kurtosis
	float meaniFreq   = 0.0f;
	float p=0.0f, ep=0.0, var=0.0, curt = 0.0f;


#if AMR_DEBUG
	printf("\r\n #InstFeatures \r\n");
#endif
	if (AMR_modNum != 7) // not AM 
	{
		// The size of input array is smaller than default frame size
		if (T < NFFT){
			numSeg = 1;
			//zero-padding for FFT at 15.05.26
			for(j=T ; j<=NFFT ; j++){
				// maximum +- 0.6 error compare to MATLAB
				AMR_IQsegmentSym[2*j+1] = 0.0f; 
				AMR_IQsegmentSym[2*j+2] = 0.0f; 
			}
		}else{
			numSeg = T/NFFT;
			if(numSeg > 19) numSeg = 19; // to increase the AMR speed, recommend odd number at 150610
		}

		if(iLen > AMR_LEN){iLen=AMR_LEN;}
	

		for(j=0; j<NFFT; j++)	AMR_spectrumMagnitude[j] = 0.0f;  // for FFT
		for(j=0; j<iLen; j++)	AMR_spectrumBuffer[j] = (float)j; // Reuse array for linear curve fitting,  ASK

		for(outLoop=0; outLoop<numSeg; outLoop++) {
			// Divide signal with size of frame size
			segStart = outLoop* NFFT;
			segEnd = segStart + NFFT;

			//////////////////////////////////////////////////
			// Calculate the instantaneous phase
			//////////////////////////////////////////////////

			k=0; AMR_IQsegmentSym[0] = 0.0f; m=0; //tmpf = 0.0; tmpf2=0.0;
			for(j=segStart ; j<segEnd ; j++){
					// maximum +- 0.6 error compared with MATLAB
				AMR_iPhase[k] = atan2(ifeIQdata[2*j+2],ifeIQdata[2*j+1]); // wraped instantaneous phase
				k++;
				if (AMR_linearFlag){
					// remove mean
					// maximum +- 0.6 error compared with MATLAB
					AMR_IQsegmentSym[2*m+1] = ifeIQdata[2*j+1]; 
					AMR_IQsegmentSym[2*m+2] = ifeIQdata[2*j+2]; 
					//tmpf  += AMR_IQsegmentSym[2*m+1];
					//tmpf2 += AMR_IQsegmentSym[2*m+2];
					m++;			
				}
			}

			// 처음 몇 표본의 값이 비정상적으로 크다. 70표본만큼 왼쪽 쉬프트, 끝 70표본은 마지막 값으로 채움
			// 50에서 70으로 변경 at 15.06.16
			if (outLoop == 0){
				for(j=0; j<NFFT-70; j++){ // NFFT is always larger than 70
					//AMR_IQsegmentSym[2*j+1] = AMR_IQsegmentSym[2*(j+70)+1];
					AMR_iPhase[j] = AMR_iPhase[j+70];
				}
				for(j=0; j<70; j++){    
					//AMR_IQsegmentSym[2*(L-70+j)+1] = AMR_IQsegmentSym[2*(L-1)+1];
					AMR_iPhase[NFFT-70+j] = AMR_iPhase[NFFT-1];   
				}
			}	

	
			// make range from -PI to +PI (=unwrap)
			//AMR_ck[0] = 0.0f; meaniFreq = 0.0f;
			//for(j=1; j<NFFT ; j++){
			//	iFreq = AMR_iPhase[j] - AMR_iPhase[j-1];
			//	if( iFreq > PI)
			//		AMR_ck[j] = AMR_ck[j-1] - TWOPI;
			//	else if ( iFreq < (-1.0f)*PI)
			//		AMR_ck[j] = AMR_ck[j-1] + TWOPI;
			//	else	AMR_ck[j] = AMR_ck[j-1];

			//	AMR_iPhase[j-1] += AMR_ck[j-1]; // unwraped instantaneous phase
			//	// instantaneous frequency = differentiation of instantaneous phase
			//	// length is reduce to L-1
			//	if (j>1){ 
			//		AMR_iFrequency[j-2] = AMR_iPhase[j-1]-AMR_iPhase[j-2];
			//		meaniFreq += AMR_iFrequency[j-2];
			//	}
			//}
			AMR_ck[0] = 0.0f; meaniFreq = 0.0f;
			for(j=1; j<NFFT ; j++){
				iFreq = AMR_iPhase[j] - AMR_iPhase[j-1];
				if( iFreq > PI)
					AMR_ck[j] = AMR_ck[j-1] - TWOPI;
				else if ( iFreq < (-1.0f)*PI)
					AMR_ck[j] = AMR_ck[j-1] + TWOPI;
				else	AMR_ck[j] = AMR_ck[j-1];

				AMR_iPhase[j-1] += AMR_ck[j-1]; // unwraped instantaneous phase
				// instantaneous frequency = differentiation of instantaneous phase
				// length is reduce to L-1
				if (j>1){ 
					AMR_iFrequency[j-2] = AMR_iPhase[j-1]-AMR_iPhase[j-2];
					meaniFreq += AMR_iFrequency[j-2];
				}
			}
			AMR_iPhase[NFFT-1] += AMR_ck[NFFT-1]; //

			// consider last two elements of instant. freq.
			AMR_iFrequency[NFFT-2] = AMR_iFrequency[NFFT-3];
			AMR_iFrequency[NFFT-1] = AMR_iFrequency[NFFT-2];
	
			leniFreq = (NFFT-2);
			meaniFreq /= leniFreq;  // float type divided by integer type
		
			// ASK에도 위상이 점프하는 구간이 생긴다. 왜그럴까.. 그래서 긴 구간 보다는 짧은 구간내에서 관찰
			// 하도록 함. 길이 iLen -- > 2015.05.27, 
			// 마지막 배열 원소값을 처리하지 않아서 생긴 문제. 15.06.23
			fit(AMR_spectrumBuffer, AMR_iPhase, iLen, &a, &b, &siga, &sigb, &chi2);
			AMR_lineLMSerror[outLoop] = chi2;	// refer to testLMScurve.m in MATLAB
			
			////////////////////////////////////////////////////////////////////////
			// Compute spectrum
			////////////////////////////////////////////////////////////////////////
			if (AMR_linearFlag){				
				// remove mean  --> unnecessary 15.05.26
				//tmpf  /= L;	 tmpf2 /= L;
				//for(j = 0; j < L ; j++){
				//	AMR_IQsegmentSym[2*j+1] 	-= tmpf;
				//	AMR_IQsegmentSym[2*j+2] 	-= tmpf2;
				//}				
			}else{ // if (AMR_linearFlag)
				// Filter length is heuristic value, you can change the value you want
				movingAverageFilter(AMR_iFrequency, NFFT, 25, AMR_filteredSamples);  // 20 -->35 for DMR at 08-Sep-2015
#if VISUALSTUDIO
		fa=fopen("filename.txt","w");
		fprintf(fa,"data =[");
		for(i=0; i<NFFT; i++) fprintf(fa,"%0.4f,...\n",(AMR_filteredSamples[i]));
		fprintf(fa,"]; \n");
		fprintf(fa,"figure; plot(data)");
		fclose(fa);			
#endif
				// Kurtosis of the normalized-centered instantaneous frequency
				// Source from "Numerical receipts"
				ep = 0.0f; p=0.0f; var = 0.0; curt = 0.0f;
				for(j=0; j<leniFreq ; j++){
					AMR_iFrequency[j] -= meaniFreq;		
					ep += AMR_iFrequency[j];
					p = AMR_iFrequency[j]*AMR_iFrequency[j];
					var += p;
					curt += p*p;
				}
				var = (var-ep*ep/(leniFreq))/(leniFreq-1);		//corrected two-pass formula
				AMR_mu42f[outLoop]=(curt)/(leniFreq*(var)*(var)); //-3.0f;
				
				// remove mean --> unnecessary 15.05.26
				//mean = 0.0;
				//for(j = 0; j < L-1 ; j++){
				//	mean += AMR_filteredSamples[j];
				//}
				//tmpi = L-1;
				//mean /= tmpi;
				//for(j=0; j<NFFT-1; j++){
				for(j=0; j<NFFT; j++){
				//	AMR_filteredSamples[j]		-= mean;
					AMR_IQsegmentSym[2*j+1] 	= AMR_filteredSamples[j];
					AMR_IQsegmentSym[2*j+2] 	= 0.0f;
				}

				/////////////////////////// Under construction. 15.06.25
				// //start = 0, end = 0;
				//if(!isFM){// FM 아닌 경우에만 실행
				//	countFlag = 0;
				//	for(j=1; j<NFFT; j++){
				//	// Time domain analysis for FM vs FSK, refer to testFindZeroCrossing.m in MATLAB
				//		if(sign(AMR_filteredSamples[j]) != sign(AMR_filteredSamples[j-1])){
				//			if(countFlag == 0){
				//				nSamp++;
				//				countFlag = 1;
				//				//start = j;
				//			}else{
				//				//end = j;
				//				if (nSamp < minSPS ) isFM = 1;
				//				nSamp = 0;
				//			} // if(countFlag == 0)
				//		} // if(sign(AMR_filteredSamples[j]) != sign(AMR_filteredSamples[j-1]))
				//		if( countFlag == 1) nSamp++;
				//	}// for(j=1; j<NFFT; j++)
				//} // if(!isFM)
				////////////////////////////////////////////////////
				//zero-padding for FFT
				//AMR_IQsegmentSym[2*NFFT-1] = AMR_IQsegmentSym[2*L-3]; //0.0;
				//AMR_IQsegmentSym[2*NFFT]   = AMR_IQsegmentSym[2*L-2];
			} // linearFlag
	
			// Oerder and H.Meyr algorithm = Squaring method in MATLAB
			nonlinearTransform(AMR_IQsegmentSym, NFFT, 2, 1);
					
			four1(AMR_IQsegmentSym, NFFT, 1); // calculate FFT
 
			// get magnitude
			for (j = 0; j < NFFT; j++){
				AMR_spectrumMagnitude[j] += myCabs(AMR_IQsegmentSym[2*j+1],AMR_IQsegmentSym[2*j+2]);
			}
		}// end for(i=0; i<numSeg; i++)

/////////////////////////// Under construction. 15.06.25
//		if(isFM) {
//			AMR_modNum = 1;		// FM
//			AMR_modOrder = 0;
//			AMR_digitalFlag = 0;
//			AMR_coarseSymbolRate = 0.0;
//#if AMR_DEBUG
//			printf("\r\n The modulation type is FM \r\n");
//#endif
//		}else{
//			AMR_modNum = 4;		// FSK
//			AMR_digitalFlag = 1;
//#if AMR_DEBUG
//			printf("\r\n The modulation type is FSK \r\n");
//#endif
//		}
/////////////////////////////////////////////////
	}

	// Exclude outliers		
	AMR_lineLMSerror[0] = dataConditining(AMR_lineLMSerror, numSeg); // edited at 150319
	AMR_mu42f[0] = dataConditining(AMR_mu42f, numSeg); // edited at 150319

#if AMR_DEBUG
	//printf(" lineLMSerror(x10) : %d \r\n", (int)(AMR_mu42f[0]);			
	printf(" mu42f(x10) : %d \r\n", (int)(AMR_mu42f[0]*10.0f));
#endif

	for (j = 0; j < NFFT; j++){
		AMR_spectrumMagnitude[j] /= numSeg;								   // take average
		//tmpf =  (10.0f)*log10(AMR_spectrumMagnitude[j]);
		AMR_spectrumMagnitude[j] = (10.0f)*log10(AMR_spectrumMagnitude[j]); // convert to dB scale
	}

#if VISUALSTUDIO
		fa=fopen("filename.txt","w");
		fprintf(fa,"data =[");
		for(i=0; i<NFFT; i++) fprintf(fa,"%0.4f,...\n",(AMR_spectrumMagnitude[i]));
		fprintf(fa,"]; \n");
		fprintf(fa,"figure; plot(data)");
		fclose(fa);			
#endif

}


#if AMR_Cplus	
	inline unsigned int AMRalgorithm::peakDetection(float* pdInput, signed int pdInputlength, signed int pdMinPeakDist,float threshold)
#else
	unsigned int peakDetection(float* pdInput, signed int pdInputlength, signed int pdMinPeakDist,float threshold)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - Detect the peaks in the magnitude spectrum
//
// \ INPUT ARGUMENTS
//	> pdInput			- spectrum
//	> pdInputlength	    - 
//	> pdMinPeakDist		-
//  
// \ OUTPUT ARGUMENTS
//	>  (global) - 
//
// \ Example
//
// \ See also 
//
// \ Author(s) : AWH
// 
// \ Reference
//  [1] "Numerical Recipes in C", 1992 pp. 613
//
/******************************************************************************/
	signed int i=0, outloop=0, k=0, count=0;
	signed int L = pdInputlength;
	signed int THRnumPeaks = 0, LMRnumPeaks=0;

	int lowerLocs =0, upperLocs = 0;
	int leftIdx =0, rightIdx = 0;
	short isNotLocalMax = 0;
	//signed int dcIdx = (L >> 1) - 1;
//	signed int distLimit = 0;
	//signed int dist=0;
	float minVal = 0.0f;
	float minHeight = 0.0f;

	// Threshold rule
	for (i=0; i<L; i++){		//140724
		if (pdInput[i] > threshold ){
			AMR_peakIndices[THRnumPeaks] = i;
			AMR_peakValues[THRnumPeaks]  = pdInput[i];
			THRnumPeaks++;
		}
	}
	
	///////////////////// Local maximum rule (LMR) /////////////////////////////
	// Sort descend order ( 3, 2, 1 )
	if(THRnumPeaks){
		// Start withthe larger peaks to make sure we don't accidentally keep a small peak
		// and remove a large peak in tis neighborhood.

		heapsort(AMR_peakValues, AMR_peakIndices, THRnumPeaks);

		for(outloop=0; outloop<THRnumPeaks; outloop++){
			if(AMR_peakIndices[outloop] >= 0){
		        // If the peak is not in the neighborhood of a larger peak, find
		        // secondary peaks to eliminate.
				lowerLocs = AMR_peakIndices[outloop] - pdMinPeakDist;
				if(lowerLocs < 0) lowerLocs = 0;

				upperLocs = AMR_peakIndices[outloop] + pdMinPeakDist;
				if(upperLocs > L) upperLocs = L;

				//Check if the current sample higher than left and right samples
				isNotLocalMax =0;
				for (i=lowerLocs; i<AMR_peakIndices[outloop]; i++){
					if(AMR_peakValues[outloop] < pdInput[i]){
						isNotLocalMax = 1;
					}else{}
				}

				for (i=AMR_peakIndices[outloop]+1; i<=upperLocs; i++){
					if(AMR_peakValues[outloop] < pdInput[i]){
						isNotLocalMax = 1;
					}else{}
				}

				for (i = outloop+1; i < THRnumPeaks; i++){
					if( ( AMR_peakIndices[i] >= AMR_peakIndices[outloop] - pdMinPeakDist )
						&&  (AMR_peakIndices[i] <= AMR_peakIndices[outloop] + pdMinPeakDist) ){
							AMR_peakValues[i]=0.0f;
							// set peakIdx. smaller than local maximum peak as negative
							AMR_peakIndices[i]=-1;
					}
				}
				if (isNotLocalMax){
					AMR_peakIndices[outloop] = -1;	// Delete current samples
				}else{
					// Check current peak about 1.6 times larger than surrounding values
					//minHeight = 1.3f;  // heuristic value
			
					//if (AMR_peakValues[outloop] < minHeight){
					//	AMR_peakIndices[outloop] = -1; // Delete current samples
					//	AMR_peakValues[outloop]  = 0.0f;
					//}else{
						// Keep current samples
					//}
				} // end for (k = 13; k < 18; k++)
			}
		} // end for for(j = 0; j < THRnumPeaks; j++)
	} // end for if(THRnumPeaks)

	count = 0; LMRnumPeaks = 0;
	//Calculate the number of peaks
	for(i=0 ; i<THRnumPeaks ; i++){
		if(AMR_peakIndices[i] >= 0){
			AMR_peakIndices[count] = AMR_peakIndices[i];
			AMR_peakValues[count]  = AMR_peakValues[i];
			LMRnumPeaks++;
			count++;
		}
	}

	if(count == 2){
		// Check tone magnitude
		for(outloop=0; outloop<=1; outloop++){
			// To divide pdMinPeakDist by 4 is no reason
			leftIdx = AMR_peakIndices[outloop] - pdMinPeakDist>>2;
			rightIdx = AMR_peakIndices[outloop] + pdMinPeakDist>>2;

			minVal = pdInput[leftIdx];
			for(i=leftIdx+1; i<=rightIdx; i++){
				if(pdInput[i] < minVal) minVal = pdInput[i];
			}

			if( (AMR_peakValues[outloop] - minVal) < 20.0f)	{count--; LMRnumPeaks--;}
		}
	}
	if(count >0){
		#if AMR_Cplus
			heapsortInt(AMR_peakIndices,AMR_peakIndices, count); //descending order
		#else
			qsort(AMR_peakIndices, count, sizeof(int), myCompareInt);
		#endif
	}else{}
//	// Calculate distance between peaks
//	// 피크가 2개일때 중심을 기준으로 대칭인지 확인한다. --> 뭔가 생각과는 다른 결과
//	// 생각보다 오차가 큼, 일단 보류
//	if(LMRnumPeaks == 2){
//		//distLimit = (signed int)(500.0f / AMR_freqInterval) ;
//		//9 표본이면, 512D, NFFT=1024 에서 267*9 = 2.4kHz에 해당한다.
//		//			, 256D, NFFT=8192 에서는 66*9 = 600Hz에 해당한다.
//		//어떻게 설정하는것이 좋을까
//		distLimit = 9;
//#if AMR_DEBUG
//				printf("not symmetric peaks");
//#endif
//		dist = AMR_peakIndices[0] + AMR_peakIndices[1] - 2*dcIdx;
//		if (dist >= 0) // positive
//		{
//			if (dist > distLimit){
//#if AMR_DEBUG
//				printf("not symmetric peaks, distance : %d ",dist);
//#endif
//				//invailid peaks
//				LMRnumPeaks = 0;
//				for(i=0; i<count; i++){
//					AMR_peakIndices[i] = -1.0f;
//					AMR_peakValues[i] = 0.0f;
//				}
//			}
//		}
//		else if (dist < 0 ){  //negative
//			dist = -dist;
//			if ( dist > distLimit){
//#if AMR_DEBUG
//					printf("not symmetric peaks, distance : %d ",dist);
//#endif
//				//invailid peaks
//				LMRnumPeaks = 0;
//				for(i=0; i<count; i++){
//					AMR_peakIndices[i] = -1.0f;
//					AMR_peakValues[i] = 0.0f;
//				}
//			}
//		}
//
//		//if ( rightSpacing - leftSpacing > distLimit) // 
//		//	LMRnumPeaks = 0;
//	}
		
	
	AMR_numPeaks = LMRnumPeaks;

	return LMRnumPeaks;
}

#if AMR_Cplus	
	inline float AMRalgorithm::coarseSymRateEstimation(float* cseSpectrum, int cseNFFT) // only for digital signal
#else
	float coarseSymRateEstimation(float* cseSpectrum, int cseNFFT) // only for digital signal
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - Estimate symbol rate based on squaring method 
//
// \ INPUT ARGUMENTS
//	> cseSpectrum	   - 
//	> cseNFFT		   - 
//  
// \ OUTPUT ARGUMENTS
//	> AMR_coarseSymbolRate (global) - The estimated symbol rate
//
// \ See also 
//
// \ Author(s) : AWH
// 
// \ Reference
//  [1] Oerder, M. and H. Myer, "Digital Filter and Square Timing Recovery,"  
//		IEEE Transactions on Communications, Vol. COM-36, No. 5, May 1988, pp.  605-612.  
//  [2] J.R. Barry, E. A. Lee, and D. G. Messerschmitt, "Digital communication, 3rd",
//	    Kluwer Academic, pp. 743-748
/******************************************************************************/
	int i=0, j=0, winLen =0;
	int L 		= cseNFFT;
	int Ld2 	= L >> 1;
	int segLen = 0;
	int leftIdx=0,rightIdx=0;
	int numSymbolRateLine = 0;
	int locs		      = 0;
	int offset		      = 0;
	int rectPulseCount    = 0;
	int minDist           = 0;

	//float maxVal = 0.0f, minVal = 0.0f;
	float coarseSR		= 0.0f;
	float lowerBound    = 0.0f, upperBound=0.0f;
	float amp_threshold = 0.0f;


#if AMR_DEBUG
	printf("\r\n #SRE \r\n");
#endif

	// Exclude AM,CW
	if (AMR_modNum != 7) // 
	{
		// set half frequency vector
		for(i=0; i<Ld2; i++) AMR_halfFreqVector[i] = (i * AMR_samplingFrequency) / (float)L;

		if (AMR_linearFlag){
			segLen=0; 
			lowerBound = 0.4f*AMR_coarseBandWidth;  // 0.5f --> 0.4f 2nd board
			upperBound = 1.5f*AMR_coarseBandWidth; // 1.5->1.1f->1.5f at 01-09-2015 by AWH
			for (i=0; i < Ld2 ; i++){
				if ( (lowerBound<= AMR_halfFreqVector[i]) && 
						(AMR_halfFreqVector[i] <= upperBound) ){
					cseSpectrum[segLen] = cseSpectrum[i];
					AMR_halfFreqVector[segLen] = AMR_halfFreqVector[i];
					segLen++;
				}
			}
		}
		else
		{
	        // Move spectrum near to zero and search peaks in the range from
	        // minimum symbol rate and maximum symbol rate
			segLen=0;
			lowerBound = AMR_minimumSymbolRate; //AMR_coarseBandWidth/5.0f; // invailid value due to POCSAG protocol
			upperBound = AMR_maximumSymbolRate;
            for (i=0; i < Ld2 ; i++){
            	if ( (lowerBound <= AMR_halfFreqVector[i]) && 
						(AMR_halfFreqVector[i] <= upperBound)){
            		cseSpectrum[segLen] = cseSpectrum[i];
            		AMR_halfFreqVector[segLen] = AMR_halfFreqVector[i];
				    segLen++;  // 397 in DSP, 397 in MATLAB
            	}
            }
		}
		
	    // use a moving average filter with a tmpi span to smooth all the data at once
		winLen = 13; // at 14-Sep-2015 ,  //at 15.06.24, ,,  29; // segLen / 2;  at 11-Sep-2015     //filter span, float->int, 
		movingAverageFilter(cseSpectrum, segLen, winLen, AMR_filteredSamples);

		// Subtract continuous part from AMR_spectrumMagnitude
        for (i=0; i < segLen ; i++){
        	cseSpectrum[i] -= AMR_filteredSamples[i];
			//tmpf =  pow(10.0f,cseSpectrum[i]/10.0f); // convert to linear scale
			cseSpectrum[i] = pow(10.0f,cseSpectrum[i]/10.0f);
		}
#if VISUALSTUDIO
		fa=fopen("filename.txt","w");
		fprintf(fa,"data =[");
		for(i=0; i<segLen; i++) fprintf(fa,"%0.4f,...\n",(cseSpectrum[i]));
		fprintf(fa,"]; \n");
		fprintf(fa,"figure; plot(data)");
		fclose(fa);			
#endif
        // Search the global maximum peak within range
        AMR_symratePeak = cseSpectrum[0];
        for (i=1; i < segLen ; i++){
        	if(cseSpectrum[i] >AMR_symratePeak){
        		AMR_symratePeak  = cseSpectrum[i];
        		coarseSR = AMR_halfFreqVector[i];	
				AMR_symratePeakIdx = (signed int)i; //type conversion for FM vs FSK classification
        	}
  		}

		if (AMR_linearFlag==0){
			// Check harmonic symbol rate peaks	150616 for RECT. pulse
			
			/////////////////////////// Under construction /////////////////////////////
			// 현재 검출된 심볼율 피크 좌측에 또 다른 심볼율 피크가 있는지 확인한다.
			// 다 수의 심볼율 피크가 있다. 구형 펄스를 사용했을 경우
	
				// peak detection
			minDist       = (signed int) (AMR_minimumSymbolRate / AMR_freqInterval);
			amp_threshold = 2.0f; //heuristic value, for MSK, 18-Sep-2015. AMR_symratePeak / 2.0f;
								  // IQ data num = 40,000
			#if AMR_DEBUG
					printf("\r Len:%d, MINPEAKDIST:%d \r\n",segLen,minDist);
			#endif
			numSymbolRateLine = peakDetection(cseSpectrum, segLen, minDist,amp_threshold);
						
			if(numSymbolRateLine>1){
				#if AMR_DEBUG
					printf("\r Harmonic peaks search \r\n");
				#endif
				offset = 3;
				// AMR_peakIndices[] is dscending order
				rectPulseCount = 0;
				locs = AMR_peakIndices[numSymbolRateLine-1] + minDist;
				for(i=numSymbolRateLine-2; i>=0; i--){
					AMR_peakIndices[i] += minDist;

					leftIdx  = (numSymbolRateLine-i)*locs - offset;
					if (leftIdx<0) leftIdx = 0;
					
					rightIdx = (numSymbolRateLine-i)*locs + offset;
					if (rightIdx > segLen-1 ) rightIdx = segLen-1;
					
					if((AMR_peakIndices[i] >= leftIdx) &&
						(AMR_peakIndices[i] <= rightIdx)){
						rectPulseCount = rectPulseCount + 1;   
					}
				}
				
			     if (rectPulseCount >= (numSymbolRateLine>>1)){
					#if AMR_DEBUG
						printf("\r New SR \r\n");
					#endif	
                   //isRectPulse = 1;
                   AMR_symratePeakIdx = AMR_peakIndices[numSymbolRateLine-1] ;
				   coarseSR = AMR_halfFreqVector[AMR_symratePeakIdx];
				 }
			}

			// Is local maximum?
	/*		leftIdx = (unsigned int)AMR_symratePeakIdx;
			rightIdx = (unsigned int)(AMR_symratePeakIdx+5);
			for(i=leftIdx; i<=rightIdx; i++){
				if( cseSpectrum[i] > AMR_symratePeak){   			
					AMR_symratePeak  = cseSpectrum[i];
        			coarseSR = AMR_halfFreqVector[i];
					AMR_symratePeakIdx = (signed int)i; 
				}
			}*/
		}
		#if AMR_DEBUG
			printf(" Symbol rate %d[sps] \r\n ",(int)coarseSR);
		#endif
		AMR_symRateEstSegLen = segLen;

		return coarseSR;
	
	}
	else // if (AMR_digitalFlag) , case of analog signal
	{
		#if AMR_DEBUG
			printf(" Symbol rate %d[symbol/s] \r\n ",0);
		#endif
		return 0.0f;
	}
}

#if AMR_Cplus	
	inline void AMRalgorithm::lowPassFiltering(float* lpfIQdata)
#else
	void lowPassFiltering(float* lpfIQdata)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - filter the IQ samples using a digital low-pass filter
//
// \ INPUT ARGUMENTS
//	> lpfIQdata		- IQ data 
//	> lpfFFTLength  - the length of FFT 
//  
// \ OUTPUT ARGUMENTS
//	> AMR_filteredInputData (global) - filtered IQ samples
//
// \ Example
//
// \ See also FILTER //  FIR1, FILTER in MATLAB
//
// \ Author(s) : AWH
// 
// \ Reference
//
/******************************************************************************/
	

	float xx	 = 0.0;       // constant for hamming window
	float uu	 = 0.0;       // constant for hamming window


	float Wn     = 0.0;	     // cut-off frequency, must be between 0< Wn < 1.0 (fs/2)
	float normalizedCutoffFreq  = 0.0;
	float cutOffExpansionRatio  = 1.3f;	// edited at 15.14.09
	float sumC					= 0.0;
	float pi 					= AMR_PI;   
	float m = 0.0, b1=0.0, cosArg = 0.0f, sincArg = 0.0f;

	unsigned int i=0,j=0;
	unsigned int N 	 = AMR_LPFORDER;   // LPF order
	unsigned int FL    = AMR_LPFLEN;   // LPF length
	unsigned int FLd2  = FL >> 1;
	unsigned int loopLimit = 0;
	unsigned int s=0;

	signed int   L      = AMR_NFFT;
	signed int   Ld2    = L >> 1;
	signed int   nominator;
	
#if AMR_DEBUG
	printf("\r\n #LPF \r\n");
#endif

	if( AMR_modNum !=7){ // not AM
		// init
		for (i=0; i<FLd2; i++) AMR_b[i] = 0.0f;

		// Distance from center index (lpfFFTLength/2)
		nominator = (AMR_BWSampleIdx[1]) - Ld2;
		// Normalized cutoff frequency
		normalizedCutoffFreq = (float) nominator / (float) Ld2  ; // e-3 order
		Wn = cutOffExpansionRatio*normalizedCutoffFreq;		 // 0.2416, 0.0219,   e-2 order
			
		// set coefficients
		AMR_F[0]=  0.0f; AMR_F[1]= Wn/2.0f; AMR_F[2] = Wn/2.0f; AMR_F[3]= 0.5f;				    
		AMR_M[0] = 1.0f; AMR_M[1] = 1.0f; AMR_M[2] = 0.0f; AMR_M[3] = 0.0f;	

		// get hamming window coefficients
		uu = 1.0f / (float)N;
		if (FL%2 == 1)			// Odd length window
		{
			loopLimit = (FL+1)>>1;
			for(j=0; j<loopLimit; j++) {
				//sumWin += (WINDOW(j,facm,facp));
				xx= (float)j * uu;
				AMR_hammingWindow[j] =  0.54f - 0.46f*cos(2.0f*pi*xx); // symmetric hamming window
				if (j>0) AMR_hammingWindow[FL-j] = AMR_hammingWindow[j-1];
			}
		}
		else{					// Even length window
			loopLimit = FLd2;
			for(j=0; j<loopLimit; j++) {
				//sumWin += (WINDOW(j,facm,facp));
				xx= (float)j * uu;
				//cosVal = cos(2.0f*pi*xx);
				AMR_hammingWindow[j] =  0.54f - 0.46f*cos(2.0f*pi*xx); // symmetric hamming window
				AMR_hammingWindow[FL-j-1] = AMR_hammingWindow[j];    // e-1 order
			}
		}

		// refer to firls.m in MATLAB
		for (i=0; i<FLd2; i++)  AMR_k[i] = (float)i+ 0.5f;	//e+1 order, Type II linear phase FIR
	
		for (s=0; s<(FL>>3); s+=2){
			m = (AMR_M[s+1]-AMR_M[s]) / (AMR_F[s+1]-AMR_F[s]);    //  slope
			b1 = AMR_M[s]-m*AMR_F[s];                  	//  y-intercept
		
			for (i=0; i<FLd2; i++){
				cosArg = 2.0f*pi*AMR_k[i];
				AMR_b[i] +=  m/(4.0f*pi*pi) * ( cos(cosArg*AMR_F[s+1])
					-cos(cosArg*AMR_F[s])) / (AMR_k[i]*AMR_k[i]);
				sincArg = 2.0f*AMR_k[i];
				AMR_b[i] += AMR_F[s+1]*(m*AMR_F[s+1]+b1)*sinc(sincArg*AMR_F[s+1])
					-AMR_F[s]*(m*AMR_F[s]+b1)*sinc(sincArg*AMR_F[s]);
			}
		}

		for (i=0; i<FLd2; i++)	AMR_a[i] = 4.0f*AMR_b[i];	//e-2 order
		for (i=0; i<FLd2; i++){
			AMR_h[i]	  = 0.5f * AMR_a[FLd2-i-1];   // e-3 order
			AMR_h[FLd2+i] = 0.5f * AMR_a[i];
		}
		// Window impulse response to get the filter
		for (i=0; i<FL; i++){
			AMR_LPFCoeff[i]  = AMR_h[i]*AMR_hammingWindow[i]; // e-3 order
			sumC		    += AMR_LPFCoeff[i];
		}
		// unity gain at DC
		for (i=0; i<FL; i++)	AMR_LPFCoeff[i] /= sumC;
			
		filter(AMR_LPFCoeff, lpfIQdata, FL, AMR_LEN, AMR_filteredInputData,0);		// included in movingAverageFilter.h
	}else{

	}
}

#if AMR_Cplus	
	inline void AMRalgorithm::linearVsNonLinearClassification(float* lvnIQ, unsigned int lvnIQLen, unsigned int lvnSamplesPerFramge)
#else
	void linearVsNonLinearClassification(float* lvnIQ, unsigned int lvnIQLen, unsigned int lvnSamplesPerFramge)
#endif

{
/********************************************************************************/
// \ DESCRIPTION
//  - Classify between linear and non-linear modulation type
//
// \ INPUT ARGUMENTS
//	> lvnIQ			       - IQ samples
//	> csrIQDataLen		   - The length of IQ data 
//	> csrSamplesPerFramge  - The number of samples in a segment
//  
// \ OUTPUT ARGUMENTS
//	> AMR_linearFlag (global) - a flag bit,
//						1 --> linear modulation signal
//						0 --> nonlinear modulation signal
// \ Example
//
// \ See also 
//
// \ Author(s) : AWH
// 
// \ Reference
//  [1] Azzouz, E.E., and Nandi, A.K.‘Automatic modulation recognition of
//		communication signals’, Kluwer Academic, 1996
/******************************************************************************/
	unsigned int i=0, j=0, k=0; //, segStart=0, segEnd=0;
	unsigned int T	 		  = lvnIQLen;
	unsigned int L 			  = lvnSamplesPerFramge;

	unsigned int numSeg		  = T / L; 			// the number of segments, float -> int
	unsigned int state1 = 0;
	unsigned int state2 = 0;


	float m_a  = 0.0, meanSquareAcn = 0.0, meanAcn = 0.0 ;
	float median_sigma_a = 0.0;
	float sigma_a_th = AMR_SIGMA_SQUARE_A_THRESHOLD;

	#if AMR_DEBUG
	printf("\r\n #linearVsNonLinearClassification\r\n");
	#endif

	//if (AMR_modNum != 7) // not AM,CW,2ASK
	if (AMR_existsCarrierFlag == 0)
	{
		///////////////////////////////////////////////////////////////////////////////
		//  Instantaneous information
		///////////////////////////////////////////////////////////////////////////////
		if (numSeg >10) numSeg = 10; // 15.06.09
		// calculate sigma_a feature value to classify the between linear and non-linear
		for(i=0; i<numSeg; i++) {
			// Divide signal with size of samplesPerFrame
			//segStart = i*L;
			//segEnd = segStart + L;

			k=0; m_a=0.0;
			//for(j=segStart; j<segEnd; j++){	
			for(j=(i*L); j<((i+1)*L); j++){	
				AMR_iAmplitude[k] = myCabs(lvnIQ[2*j+1],lvnIQ[2*j+2]);     // instantaneous amplitude
				m_a += AMR_iAmplitude[k];        // mean of instant. amp.  ,  0.89 in MATLAB, 0.95 in DSP
				k++;
			}
			m_a /= (float)L;

			meanSquareAcn = 0.0f;  meanAcn = 0.0f;
			for (j=0; j<L; j++){
				// tmpf3 = AMR_iAmplitude[j]/m_a;
				//AMR_iAmplitude[j] = tmpf3 - 1.0f;		// a_cn  Eq. (2.2) in [1]
				AMR_iAmplitude[j] = AMR_iAmplitude[j]/m_a - 1.0f;

				meanSquareAcn  += (AMR_iAmplitude[j] * AMR_iAmplitude[j]); // mean(a_cn.^2)
				meanAcn += AMR_iAmplitude[j] ; 							   // mean(a_cn)
			}
			meanSquareAcn /= L;
			meanAcn /= L;

		    // 6. sigma_a
		    // Standard deviation of the normalized-centered instantaneous --> variance 15.05.27
		    // in the non-weak intervals of a signal segment
			//tmpf3 = meanSquareAcn - (meanAcn*meanAcn); // mean(a_cn.^2) - mean(a_cn)^2
			AMR_sigma_a[i] = meanSquareAcn - (meanAcn*meanAcn); //sqrt(tmpf3); // 0.3110 in MATLAB , 0.3954 in DSP  +- 0.07 error
		} // end for(i=0; i<numSeg; i++)
	
		// data conditioning
		if (numSeg >1) median_sigma_a = dataConditining(AMR_sigma_a,numSeg);
		else median_sigma_a = AMR_sigma_a[0];
		AMR_sigma_a[0] = median_sigma_a; // used for classification between QAM and PSK
		
		#if AMR_DEBUG
			printf("\r (sigma_a(%u) > %u) : linear ? non-linear \r\n ",
				(unsigned int)(median_sigma_a * (10000.0f)),
				(unsigned int)(sigma_a_th * (10000.0f)));
		#endif

		state1 = (unsigned int)(median_sigma_a*10000.0f);
		state2 = (unsigned int)(sigma_a_th*10000.0f);
		// Classify between linear and non-linear modulation signal
		if (state1 > state2){
			AMR_linearFlag = 1;
			state1 = (unsigned int)AMR_coarseBandWidth;

			// exception handling
			if (state1 > AMR_LINEAR_MOD_MAX_BW){
				#if AMR_DEBUG
					printf("\r Linear (maybe non-Linear) \r\n ");
				#endif
				AMR_linearFlag = 0;
			}
		}else{
			AMR_linearFlag = 0;
		}
	}
}


#if AMR_Cplus	
	inline void AMRalgorithm::reSampling(int resIQdataLen, float resEstSymRate)
#else
	void reSampling(int resIQdataLen, float resEstSymRate)
#endif

{
/********************************************************************************/
// \ DESCRIPTION
//  > Decide proper decimation rate of input IQ data  
//
// \ INPUT ARGUMENTS
//	> resIQdataLen	- IQ data 
//	> resEstSymRate	 - Estimated symbol rate
//  
// \ OUTPUT ARGUMENTS
//	> Print Message
//
// \ Example
//
// \ See also 
//
// \ Author(s) : AWH
// 
// \ Reference
//
/******************************************************************************/
	unsigned int powerConst=0, i=0;
	unsigned int L 			  = resIQdataLen;
	unsigned int spsLowerBound = AMR_spsLowerBound;
	unsigned int spsUpperBound = AMR_spsUpperBound;
	unsigned int nowDecimation = 0;
	float Fs 		  = AMR_samplingFrequency;
	unsigned int sps=0;

#ifdef activateCDW
	nowDecimation = uiAMRDecimation;		// power of 2
#else
	nowDecimation = AMR_Decimation;
#endif


#if AMR_DEBUG
	printf("\r\n #ReSampling \r\n ");
#endif

	if (AMR_digitalFlag){		//Digital
		AMR_samplesPerSym = (unsigned int)(Fs / resEstSymRate);
		AMR_numSym 		  = L / AMR_samplesPerSym;

		if (AMR_samplesPerSym <= spsLowerBound){
			powerConst = 1;
			for(i = 0; i< 22; i++){
				powerConst *=  2;
				sps = AMR_samplesPerSym * powerConst;
				if( sps > spsLowerBound) break;
			}
			//AMR_decimationTable[8]={32,64,128,256,512,1024,2048,4096};
			AMR_decimationIdx = (unsigned int) (log2((float)nowDecimation / (float)powerConst) - 5.0f);
			AMR_requireIQsamples = 160000;

			//axisLowerIdx = zeroIdx - AMR_BWSampleIdx[0];
			//axisUpperIdx = AMR_BWSampleIdx[1] + zeroIdx;

			//axisLowerIdx = zeroIdx - (axisLowerIdx / tmpi);
			//axisUpperIdx = zeroIdx + (axisUpperIdx / tmpi);
			//
			//AMR_BWSampleIdx[0] = axisLowerIdx;
			//AMR_BWSampleIdx[1] = axisUpperIdx;

#if AMR_DEBUG
	printf(" Samples per symbol %d \r\n ",AMR_samplesPerSym);
	printf(" # of symbols %d \r\n ",AMR_numSym);
	printf(" Request FPGA to control the decimation rate by factor of %d \r\n ",(nowDecimation / powerConst));
#endif
		}
		else if (AMR_samplesPerSym >= spsUpperBound){
			powerConst = 1;
			for(i = 0; i< 22; i++){				
				powerConst *= 2;
				sps = AMR_samplesPerSym / powerConst;
				if( sps < spsUpperBound) break;
			}
			AMR_decimationIdx = (unsigned int)(log2((float)nowDecimation * (float)powerConst) - 5.0f);
			AMR_requireIQsamples = 160000;
			//Request FPGA to change decimation rate by factor of q/p

			//axisLowerIdx = zeroIdx - AMR_BWSampleIdx[0];
			//axisUpperIdx = AMR_BWSampleIdx[1] - zeroIdx;

			//AMR_BWSampleIdx[0] = zeroIdx - (axisLowerIdx * tmpi);
			//AMR_BWSampleIdx[1] = zeroIdx + (axisUpperIdx * tmpi);
			
#if AMR_DEBUG
	printf(" Samples per symbol %d \r\n ",AMR_samplesPerSym);
	printf(" The number of symbols %d \r\n ",AMR_numSym);
	printf(" Request FPGA to control the decimation rate by factor of %d \r\n ",(nowDecimation * powerConst));
#endif
		}else{
			AMR_decimationIdx = (unsigned int)(log2((float)nowDecimation) - 5.0f);
			AMR_requireIQsamples = 160000;
#if AMR_DEBUG	
	printf(" Samples per symbol %d \r\n ",AMR_samplesPerSym);
	printf(" The number of symbols %d \r\n ",AMR_numSym);
	printf(" Don't change the decimation rate. Keep going \r\n ");
#endif
		}
	}else{					  //Analog
		// To acquire IQ data longer than 2 seconds, we need 136,720 IQ samples with decimation rate = 2048.
		// If possible, we increase the decimation ratio to acquire the more IQ samples
		powerConst = 1024 / nowDecimation;
		//AMR_decimationTable[8]={32,64,128,256,512,1024,2048,4096};
		AMR_decimationIdx = 5;
		AMR_requireIQsamples = 136720; // approx. 2 seconds
		if (powerConst == 1){
#if AMR_DEBUG
	printf(" Don't change the decimation rate. Keep going \r\n ");
#endif
		}else if (powerConst > 1 ) {
			
#if AMR_DEBUG
	printf(" Request FPGA to control the decimation rate by factor of %d \r\n",1024);
#endif
		}else{
#if AMR_DEBUG
	printf("\r Request FPGA to control the decimation rate by factor of %d \r\n",1024);	
#endif
		}
	}
} // EOF
	
#if AMR_Cplus	
	inline signed int AMRalgorithm::setFrequencyVector(int iqLen)
#else
	signed int setFrequencyVector(int iqLen)
#endif
{
/**********************************************************************************/
// \ DESCRIPTION
//  - Set frequency bins and the length of FFT
//
//	fs : sampling frequency 
//	-fs/2 + freqResolution  ~ fs/2				    in matlab,  
//	-fs/2				    ~ fs/2 - freqResolution in C
//
// \ INPUT ARGUMENTS
//	> iqLen	    - The length of IQ data 

// \ OUTPUT ARGUMENTS
//	> L				   		     - The length of FFT
//  > AMR_minBinInterval(global) - The minimun frequency interval
//  > AMR_freqVector(global)	 - The frequency bins range from -NFFT/2 to NFFT/2-1
/***********************************************************************************/

	signed int i=0;
	signed int L =0;
	signed int Ld2 = 0;  // L divided by 2
	signed int avgNum = 0; 
	signed int comp1 = 0, comp2;
	#if AMR_DEBUG
		printf("\r\n #Freq. bin setting       ");
	#endif

    // Determine the length of FFT
	avgNum = iqLen / 8; // the count of averaging
	L = 1;
	for(i=0; i<30 ; i++){
		L *= 2; 
		if( L > avgNum) break;
	}

	// 계산한 NFFT 길이(L)를 현재 값과 가장 가까운 2의 배수로 재설정
	// 가정 : NFFT 길이는 avgNum보다 항상 크다.

	comp1 = L - avgNum;
	if(comp1 < 0 ) comp1 = -comp1;

	comp2 = (L>>1) - avgNum;
	if(comp2 < 0 ) comp2 = -comp2;
	
	if( comp1 > comp2) L = L >>1;
	
	if(L >= 8192) 		L = 8192;
	else if(L <= AMR_MIN_NFFT)   L = AMR_MIN_NFFT;
	Ld2 = L >> 1;

	// Set frequency bins
	#ifdef activateCDW
		AMR_samplingFrequency = (float) AMR_ADCsamplingFrequency /  uiAMRDecimation;		// power of 2
	#else
		AMR_samplingFrequency =  AMR_ADCsamplingFrequency /  (float) AMR_Decimation;
	#endif

	// frequency resolution = the interval between successive frequency bins
	AMR_freqInterval   = (float) AMR_samplingFrequency / L;	
	if( AMR_freqInterval > 900.0f){
		AMR_minBinInterval = (unsigned int)AMR_freqInterval;
	}else{
		AMR_minBinInterval = (unsigned int) (900.0f / AMR_freqInterval);	// 신호 스펙트럼내의 최대 빈간격 = (900 Hz / 33Hz) 150413
	}
	AMR_samPeriod 	   =  1.0f / AMR_freqInterval;
	for(i=0; i<Ld2; i++) AMR_freqVector[i] = -Ld2+i;
	for(i=Ld2; i<L; i++) AMR_freqVector[i] = i-Ld2;

#if AMR_DEBUG
	printf("\r\n NFFT : %d \r\n",L);	
#endif
	return L;
}

#if AMR_Cplus	
	inline void AMRalgorithm::FMvsFSKClassification(float* fcSpectrum, signed int fcSegLen) // only for digital signal
#else
	void FMvsFSKClassification(float* fcSpectrum, signed int fcSegLen) // only for digital signal
#endif
{
	signed int i  = 0;
	signed int L  = fcSegLen;
	//signed int l1 = 0
	signed int l2=0, r1=0; //, r2=0;
	// magnitude of spectrum corresponding to symbol rate 
	unsigned short  isFSK=0;
	float maxHeigthThreshold = 0.0f;


#if AMR_DEBUG
	printf("\r\n # FMvsFSKClassification \r\n");
#endif

	// Classification between FM and FSK
	if (AMR_modNum != 7) // not AM 
	{
		if (!AMR_linearFlag){
			// 심볼율 피크 까지의 진폭 스펙트럼을 AMR_spectrumBuffer 배열로 옮겨 담는다.
			/*for(i=0; i<=AMR_symratePeakIdx-5; i++){
				AMR_spectrumBuffer[count++] = fcSpectrum[i];	
			}*/

			// Threshold rule
			// 심볼율 선 스펙트럼의 진폭(전역최대)과 두번째로 큰 진폭 스펙트럼의 비 = 0.73
			// file:///C:/Users/elint/Documents/MATLAB/Samsung_Thales/html/testFMvsFSKt.html
			// testFMvsFSKt.m
			// 15.06.16
			maxHeigthThreshold = AMR_symratePeak*0.5f;
			//minHeigthThreshold = AMR_symratePeak*0.9f;  // 0.9 --> 0.85f under moving average fitler, lenght of 5
														 // 0.5f --> 0.7f
														 // 0.7f --> 0.68f at 150615
														 // refer to testFMvsFSKt.m in SBC  
			isFSK = 1;

			// harmonicThreshold 보다 큰 선 스펙트럼이 있는 지 검사한다.
			// 원래는 심볼율 선 스펙트럼이 심볼율의 정수배 주파수에서 나타나는 문제를
			// 해결하기 위한 코드였지만... 그런 경우는 드물다?판단하고 사용하지 않도록 한다.
			//
			//for (j=0; j<count; j++){		//140724
			//	if (AMR_spectrumBuffer[j] > harmonicThreshold ){
			//		AMR_peakIndices[THRnumPeaks] = j;
			//		AMR_peakValues[THRnumPeaks]  = AMR_spectrumBuffer[j];
			//	}
			//}

			// center part
		/*	for(i=AMR_symratePeakIdx-1; i<=AMR_symratePeakIdx+1; i++){
				srPower += fcSpectrum[i];
			}	*/

			// left part
			//l1  = AMR_symratePeakIdx-7;	  if(l1 <0)  l1 = 0;
			l2  = AMR_symratePeakIdx-5;   if(l2 < 0) l2 = 0;

			//mean=0; count=0;
			for(i=0; i<=l2; i++){
				/*mean += fcSpectrum[i];
				count++;*/
				//if (fcSpectrum[i] > minHeigthThreshold ) isFM = 1;
				if (fcSpectrum[i] > maxHeigthThreshold ) isFSK = 0;
			}

			// right part
			//if (AMR_coarseSymbolRate < 3000.0f){
			//	r1 = AMR_symratePeakIdx+53;	if(r1 > L-1)  r1 = L-1;
			//	//r2 = AMR_symratePeakIdx+152; if(r2 > Ld2-1) r2 = Ld2-1;
			//}else{
			//	r1 = AMR_symratePeakIdx+3;	if(r1 > L-1)  r1 = L-1;
			//	//r2 = AMR_symratePeakIdx+102; if(r2 > Ld2-1) r2 = Ld2-1;
			//}
			// 302 --> 102 at 150615
			r1 = AMR_symratePeakIdx+5;	if(r1 > L-1)  r1 = L-1;

			for(i=r1; i<L; i++){
	/*			mean += fcSpectrum[i];
				count++;*/
				//if (fcSpectrum[i] > minHeigthThreshold ) isFM = 1;
				if (fcSpectrum[i] > maxHeigthThreshold ) isFSK = 0;
			}			
	

#if AMR_DEBUG
				printf("\r\n (isFSK(%d) == 1 ) ? FSK  \r\n ",(int)isFSK);
				//printf("\r\n (minHeigthCondition(%d) == 1 ) FM   \r\n ",(int)isFM);
#endif


			AMR_modNum = 4;
			if (isFSK){
						// FSK
				AMR_digitalFlag = 1;
				AMR_isFSK = 1;
			}else{
				AMR_isFSK = 0;
			}
		}
	}
}

#if AMR_Cplus	
	void AMRalgorithm::coarseToneSpacing(float* ctsIQdataInput,
		signed int ctsIQdataLen,signed int ctsFFTlength, signed int ctsHighFFTlength, signed int ctsOffset)
#else
	void coarseToneSpacing(float* ctsIQdataInput,signed int ctsIQdataLen,
		signed int ctsFFTlength, signed int ctsHighFFTlength, signed int ctsOffset)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - Estimate tone spacing based on 2P power method (CPM Phase recovery) 
//
// \ INPUT ARGUMENTS
//	> ctsIQdataInput	- IQ data 
//	> ctsIQdataLen	    - The length of IQ data 
//	> ctsFFTlength		- The length of FFT used in bandwidth estimation
//  > ctsHighFFTlength  - The length of FFT larger than ctsFFTlength
//  > ctsOffset         - index corresponding to frequency offset
//
// \ OUTPUT ARGUMENTS
//	> AMR_coarseToneSpacing (global) - The estimated symbol rate
//
// \ Example
//
// \ See also 
//
// \ Author(s) : AWH
// 
// \ Reference
//  [1] Mengali, Umberto, and N. D’Andrea, "Synchronization Techniques for 
//		Digital Receivers, New York", Plenum Press, 1997.
/******************************************************************************/

	signed int numTones		  = 0;
	signed int i			  = 0,j = 0;
	signed int T              = ctsIQdataLen;
	signed int Lplus          = ctsHighFFTlength;
	signed int L			  = ctsFFTlength;
	signed int count		  = 0; //, halfIdx=0;
	signed int toneSpacingIdx = 0;
	signed int minDist        = 0;
	signed int offset         = 0;
	signed int halfIdx        = 0;
	signed int len			  = 0;
  //	signed int bwIndex        = 0;
	unsigned int power		  = 1;
	short isFM_1=0, isFM_2=0, isFM_3=0;
	float yMax = 0.0f; // BW est
//	float ratioTS2SR = 0.0f;	 // ratio tone spacing to symbol rate
	float ratioBW2SR = 0.0f;  // ratio bandwidth to symbol rate
	float ampThreshold = 0.0f;
	float mean=0.0, ss=0.0, ep=0.0, var=0.0, sdev = 0.0;
//	float threshold = 0.0; 
	float minVal = 0.0f;
	float freqInterval = AMR_freqInterval;
	float toneMag = 0.0f;
	float toneDistFreq = 0.0f;
	//int leftIdx = 0, rightIdx = 0;
	// BW est
	int toneDist = 0.0;
    signed int yMaxIdx = 0;
	signed int axisLowerIdx, axisUpperIdx=0;

#if AMR_DEBUG
		printf("\r\n #ToneSpacingEst \r\n");
#endif

	AMR_bwThreshold_dB = AMR_bwThreshold_dB-6.0f; //  bandwidth threshold

	AMR_bwThreshold	 = pow(10.0f,AMR_bwThreshold_dB/10.0f);
	// Think more ??? 
	//if (Lplus > ctsIQdataLen) Lplus = L; // If the length of input IQ is smaller than recommended FFT size
										   // keep current frequency resolution
	Lplus = L;
	 //  FM,unknown,FSK  
	if( AMR_linearFlag==0) 
	{	
		// 반송파 주파수 옵셋 추정하고 옵셋만큼 IQ 표본을 보상하는데,
		// 스펙트럼은 다시 계산하지 않으므로 추정한 대역폭을 옵셋이 있는
		// 값으로 원상복구 한다.

		AMR_BWSampleIdx[0] += ctsOffset;
		AMR_BWSampleIdx[1] += ctsOffset;

		
		while(numTones <= 1){
			if (power > 1){
				 // Decrease fft resolution at 14-Sep-2015
				if(power > 4) {
#if AMR_DEBUG
					printf("\r Tone detection fail \r\n");
#endif
					break;
				}
				if(power == 2) {
					Lplus = Lplus >>1;
					freqInterval = 2.0f*freqInterval;
				}
				// Calculate frequency resolution again
				//AMR_freqInterval = AMR_samplingFrequency / (float)Lplus;
				
				// 신호를 power승 취하고 스펙트럼을 새로 구한다.
				welchPSD(ctsIQdataInput, T, Lplus, power, AMR_powerSpectrum);

  		//	    // convert to decibel scale
				//for (i=0 ; i<Lplus ; i++)	AMR_powerSpectrum_dB[i] = (10.0f)*log10(AMR_powerSpectrum[i]);
		
				// Find min. and max. Duplicated code in histogram
			/*	yMax = AMR_powerSpectrum_dB[0]; //cancel by AWH at 03-Sep-2015
				for (i=1 ; i<Lplus ; i++){
					if(AMR_powerSpectrum_dB[i] > yMax ){
						yMax = AMR_powerSpectrum_dB[i];
					}
				}
				yMin = yMax;
				for (i=0 ; i<Lplus ; i++){
					if(AMR_powerSpectrum_dB[i] < yMin ){
						yMin = AMR_powerSpectrum_dB[i];
					}
				}*/

				// normalization			   //cancel by AWH at 03-Sep-2015
				//tmpf = 100 / (yMax-yMin);
				//for (i=0 ; i<Lplus ; i++){
				///	AMR_powerSpectrum_dB[i] = (AMR_powerSpectrum_dB[i]-yMin) * tmpf;
				//}
	
//
//				// obtain threshold again
//				histogram(AMR_powerSpectrum_dB, Lplus, AMR_MAX_NUMBINS, yMax, yMin);

//				// Find noise floor
//				tmpi = AMR_binCount[0];
//				for (i=1; i<AMR_MAX_NUMBINS; i++){
//					if (AMR_binCount[i] > tmpi){
//						tmpi = AMR_binCount[i];
//						noiseFloor = AMR_binCenter[i];
//					}
//				}
//
//				// Find local maximum
//				signalFloorCount = 0;	// bin count is always positive
//				for (i=0; i<AMR_MAX_NUMBINS; i++){
//					if (AMR_binCenter[i] > levelB2NoiseNSig){
//						if (count == 0){    // to save the first sample index higher than threshold
//							thresholdIdx = i;
//							count = thresholdIdx;
//						}
//						if(AMR_binCount[count] > signalFloorCount ){  // Find local maximum
//							signalFloorCount = AMR_binCount[count];
//							signalFloorIdx = count;
//						}
//						count += 1;
//					}
//				}
//				signalFloor = AMR_binCenter[signalFloorIdx];
//
//				AMR_bwThreshold_dB = ( signalFloor - noiseFloor ) * 0.7; 
				// get index of samples higher than threshold
			/*	count = 0;
				for (i=0; i<Lplus; i++){
					if ( AMR_powerSpectrum_dB[i] >= AMR_bwThreshold_dB ){
						AMR_TSE_indexBuffer[count++] = i; 
					}				
				}	*/		

				// Which one is the better among the either dB or linear scale
				count = 0;
				for (i=0; i<Lplus; i++){
					if ( AMR_powerSpectrum[i] >= AMR_bwThreshold ){
						AMR_TSE_indexBuffer[count++] = i; 
					}				
				}		
				// update index of samples with bandwidth
				AMR_BWSampleIdx[0] = AMR_TSE_indexBuffer[0];
				AMR_BWSampleIdx[1] = AMR_TSE_indexBuffer[count-1];
			}
			
			count=0;
			for(i=AMR_BWSampleIdx[0]; i<=AMR_BWSampleIdx[1]; i++){
				AMR_TSE_PSDbuffer[count++] = AMR_powerSpectrum[i];	
			}
			// Frequency smoothing
			len = count >>2; // window length is value of count/4
			movingAverageFilter(AMR_TSE_PSDbuffer, count, len, AMR_filteredSamples);

			// Remove average components of spectrum
			for(i=0; i<count; i++){
				AMR_TSE_PSDbuffer[i] -= AMR_filteredSamples[i];
			}
			
			// Replace DC component with minimum value
			//min = AMR_TSE_PSDbuffer[0];
			//for(i=1; i<count; i++){
			//	if(AMR_TSE_PSDbuffer[i] < min)	min = AMR_TSE_PSDbuffer[i];
			//}

			//halfIdx = (count >> 1) - 1;
			//offset = (500.0f / AMR_freqInterval) * (power>>1);
			//for(i=halfIdx-offset; i<=halfIdx+offset; i++)	AMR_TSE_PSDbuffer[i] = min;	// increase the range at 150226
#if VISUALSTUDIO
	fa=fopen("filename.txt","w");
	fprintf(fa,"data =[");
	for(i=0; i<count; i++) fprintf(fa,"%0.4f,...\n",(AMR_TSE_PSDbuffer[i]));
	fprintf(fa,"]; \n");
	fprintf(fa,"figure; plot(data)");
	fclose(fa);			
#endif

			// Calculate minimum value
			minVal = AMR_TSE_PSDbuffer[0];
			for(i=1; i<count; i++){
				if(minVal > AMR_TSE_PSDbuffer[i]){
					minVal = AMR_TSE_PSDbuffer[i];
				}
			}
			for(i=0; i<count; i++) mean += AMR_TSE_PSDbuffer[i];
			mean /= count;

			// Get standard variance in [1]
			for(i=0; i<count; i++){
				ss = AMR_TSE_PSDbuffer[i] - mean;
				ep += ss;
				var += ss*ss;
			}
			var = (var-ep*ep/count)/(count-1);		//corrected two-pass formula
			sdev = sqrt(var);

			ratioBW2SR = AMR_coarseBandWidth / AMR_coarseSymbolRate;
			// ratio 
			if (ratioBW2SR > 1.5f)	ampThreshold = mean+(2.6f)*sdev;  // 8-FSK , POCSAG
			else					ampThreshold = 1.0f; // at 14-Sep-2015;   // DMR,MSK
             //Wide tone spacing  
			// 2.7 --> 2.58 at 150226
			// 2.58 --> 2.0 at 150615 for DMR
			// 2.0 --> 1.9 for DMR at 03-Sep-2015
			// 1.9 --> 1.55					
   
			// distance between tones
			minDist = (signed int) (0.8f * AMR_coarseSymbolRate / freqInterval);
		
#if AMR_DEBUG
			minVal = (ampThreshold*10.0f);
			printf("\r Count : %d, MIN. dist : %d, AMP. THRESHOLD(x10) : %d  \r\n",count,minDist,(int)minVal);
#endif
			// Replace DC component with minimum value
			minVal = AMR_TSE_PSDbuffer[0];
			for(i=1; i<count; i++){
				if(AMR_TSE_PSDbuffer[i] < minVal)	minVal = AMR_TSE_PSDbuffer[i];
			}

			halfIdx = (count >> 1) - 1;
			offset = (int)(300.0f/freqInterval);
			for(i=halfIdx-offset; i<halfIdx+offset+1; i++)	AMR_TSE_PSDbuffer[i] = minVal;
			numTones = peakDetection(AMR_TSE_PSDbuffer, count, minDist,ampThreshold);

			//Check tone magnitude 
			if(numTones == 2){
				for(i=0; i<2; i++){
					toneMag = AMR_TSE_PSDbuffer[AMR_peakIndices[i]];
					//leftIdx = AMR_peakIndices[i] - (minDist >> 2);
					//if(leftIdx < 0) leftIdx = 0;

					//rightIdx = AMR_peakIndices[i] + (minDist >> 2);
					//if(rightIdx > Lplus-1) rightIdx = Lplus-1;

					//minVal = AMR_TSE_PSDbuffer[leftIdx];
					//for(j=leftIdx+1; j<=rightIdx; j++){
					//	if(AMR_TSE_PSDbuffer[j] < minVal){	
					//		minVal = AMR_TSE_PSDbuffer[j];
					//	}
					//}

					//toneMag -= minVal;
					if (toneMag < 20.0f){ // How about increase the value ?
						numTones -= 1;
					}
				}
			}

			power *= 2;
		}// End for while

		if (numTones >1){		


			// Check!
			toneSpacingIdx = 0;
			for(i=0; i<numTones-1 ;i++){
				toneDist = AMR_peakIndices[i] - AMR_peakIndices[i+1];
				toneDistFreq = freqInterval*(toneDist) / (float)(power >> 1);
				if( toneDistFreq < 530.0f ) isFM_3 = 1;
				toneSpacingIdx += toneDist;
			}


			AMR_modNum = 4;		// FSK
			AMR_digitalFlag = 1;

			toneSpacingIdx /= (numTones-1);
			AMR_coarseToneSpacing = freqInterval * toneSpacingIdx / (float)(power >> 1);

			// General information
			// 1) 선형변조 신호의 대역폭은 심볼율의 2.5배를 넘지 않는다.
			// 2) 비선형변조 신호인 FSK의 대역폭은 심볼율의 8.5배를 넘지 않는다.
			//    ex) 심볼율과 대역폭의 비가 가장 큰 프로토콜은 POCSAG으로 6.8배 이다.
			// 따라서 심볼율과 대역폭의 비를 계산하면 FSK와 FM을 구분할 수 있다. 

			// Check ratio bandwidth to symbol rate
			
			//if( ratioBW2SR > 7.5f) isFM_1 = 1;

			// Check ratio tone spacing to symbol rate
			//ratioTS2SR = AMR_coarseToneSpacing / AMR_coarseSymbolRate;
			//switch(power>>1){
			//	case 1: 
			//		if(ratioTS2SR < 0.8f) isFM_2 = 1; break;
			//	case 2:
			//		if(ratioTS2SR > 0.8f) isFM_2 = 1; break;
			//	case 4:
			//		if(ratioTS2SR > 0.4f) isFM_2 = 1; break;
			//	default:
			//		// invalid power value
			//		break;
			//}
			//if (AMR_coarseToneSpacing < 530.0f) isFM_3 = 1;
			
			// Check kurtosis for instantaneous frequency,
			if (AMR_mu42f[0] > AMR_FREQ_KURTOSIS_THRESHOLD) isFM_1=1;
			if (isFM_3 || isFM_1 ){
				// FM
#if AMR_DEBUG
			printf("\r\n Invalid parameters (kurtosis:%d, case2:%d, TS:%d) -> FM \r\n",isFM_1,isFM_2,isFM_3);
#endif
				AMR_modNum = 1;		// FM
				AMR_modOrder = 0;
				AMR_digitalFlag = 0;
				AMR_coarseSymbolRate = 0.0;
				AMR_coarseToneSpacing = 0.0;
				AMR_isFSK = 0;
			}else{
				// FSK
				if( numTones == 2 )						 AMR_modOrder = 2;
				else if( 3 <= numTones && numTones <=4)  AMR_modOrder = 4;
				else if( 5 <= numTones)                  AMR_modOrder = 8;

#if AMR_DEBUG
			printf("\r Tone spacing : %d[Hz], power : %d \r\n",(int)AMR_coarseToneSpacing,power>>1);
			printf("\r Mod. order of FSK : %d \r\n",AMR_modOrder);
#endif
			}
		}
		else
		{
			//if (AMR_isFSK){ // from FM vs FSK classification
			//	// FSK이지만 피크 검출에 실패한 경우
			//}else{
				
			#if AMR_DEBUG
				printf("\r Tone detection fail -> FM \r\n");
			#endif
		
			AMR_modNum = 1;		// FM
			AMR_modOrder = 0;
			AMR_digitalFlag = 0;
			AMR_coarseSymbolRate = 0.0f;
			AMR_coarseToneSpacing = 0.0f;

				
			//
			// Estimate bandwidth again.
			//
			unitPower(AMR_inputData, AMR_LEN);
			welchPSD(AMR_inputData, AMR_LEN, AMR_NFFT, 1, AMR_powerSpectrum);
#if VISUALSTUDIO
	fa=fopen("filename.txt","w");
	fprintf(fa,"data =[");
	for(i=0; i<AMR_NFFT; i++) fprintf(fa,"%0.4f,...\n",(AMR_powerSpectrum[i]));
	fprintf(fa,"]; \n");
	fprintf(fa,"figure; plot(data)");
	fclose(fa);			
#endif
			offset = (AMR_NFFT>>2);
			yMax = AMR_powerSpectrum[0]; yMaxIdx = 0;
			for(i= offset; i< (AMR_NFFT-offset); i++){					
				if ( AMR_powerSpectrum[i] > yMax ){
					yMaxIdx = i;
					yMax = AMR_powerSpectrum[i];
				}
			}
			AMR_noiseFloor_dB = AMR_noiseFloor_dB - 2.0f;
			AMR_noiseFloor_dB = pow(10.0f,AMR_noiseFloor_dB / 10.0f);
			axisLowerIdx = 0;
			for(i=yMaxIdx; i>0; i--) {
				if ( AMR_powerSpectrum[i] < AMR_noiseFloor_dB){
					axisLowerIdx = i; break;
				}
			}

			axisUpperIdx = AMR_NFFT-1;
			// toward right edge, (N+1)-N
			for(i=yMaxIdx; i<AMR_NFFT; i++) {
				if ( AMR_powerSpectrum[i] < AMR_noiseFloor_dB){
					axisUpperIdx = i; break;
				}
			}
			AMR_coarseBandWidth = (float)((axisUpperIdx-axisLowerIdx)*AMR_freqInterval);
	
			#if AMR_DEBUG
				printf("\r New estimated FM bandwidth %d \r\n ",(int)AMR_coarseBandWidth);
			#endif
			//}
		}
	}else{ // end for ( !AMR_linearFlag && AMR_modNum != 1 )
		AMR_coarseToneSpacing = 0.0;		
	}
}

#if AMR_Cplus	
	void AMRalgorithm::linearModClassification(float* feIQdata, int feIQdataLen, int feSamplesPerFrame)
#else
	void linearModClassification(float* feIQdata, int feIQdataLen, int feSamplesPerFrame)
#endif
{
/********************************************************************************/
// \ DESCRIPTION
//  - Classify the linear modulation signals 
//
// \ INPUT ARGUMENTS
//	> feIQdata	         - IQ samples
//	> feIQdataLen	     - The length of IQ data 
//  > feSamplesPerFrame  - The number of samples in a frame
//
// \ OUTPUT ARGUMENTS ( global variable )
//	> AMR_modNum	- number corresponding to modulation type
//  > AMR_modOrder  - number corresponding to modulation order
//
// \ Example
//
// \ See also 
//
// \ Author(s) : AWH
// 
// \ Reference
//
/******************************************************************************/

	signed int i=0,j=0,k=0,segStart=0,segEnd=0;
	signed int powerConst = 0;
	signed int T	= feIQdataLen;
	signed int L    = feSamplesPerFrame;
	signed int NFFT = 0;
	signed int numSeg = 0;				// the number of segments, float -> int
	signed int condition1=0, condition2=0;

	float iMean = 0.0, qMean = 0.0;
	float sigPowerNormFactor = 0.0, sigPower = 0.0;

	float ABSCC20_THRESHOLD = AMR_ABSCC20_THRESHOLD;

#if AMR_DEBUG
	printf("\r\n #linearModClassification \r\n");
#endif
	
	// The size of input array is smaller than default frame size
	if (T < L){
		L = T;
		powerConst = (signed int)log2((float)T);
		NFFT = 1;
		for(i=0; i< (powerConst+1) ; i++){
			NFFT = 2 * NFFT ;
		}
		numSeg = 1;
		//zero-padding for FFT at 15.05.26
		for(j=T ; j<=NFFT ; j++){
			// maximum +- 0.6 error compare to MATLAB
			AMR_IQsegmentSym[2*j+1] = 0.0f; 
			AMR_IQsegmentSym[2*j+2] = 0.0f; 
		}
	}else{
		NFFT = L;
		numSeg = T/L;
		if (numSeg>2) numSeg = 2;
	}
	
	if (AMR_linearFlag && AMR_modNum != 7){											   
		
		// calculate feature values
		for(i=0; i<numSeg; i++) {
			// Divide signal with size of samplesPerFrame
			segStart = i*L;
			segEnd = segStart + L;

			//  normalization and remove mean
			//// out = sqrt( length(bandLimitedSig))/norm(bandLimitedSig, 2) * bandLimitedSig; in MATLAB
			sigPower = 0.0;
			for(j=segStart ; j< segEnd ; j++){
				 // = pow(AMR_inputData[2*idx+1],2) + pow(AMR_inputData[2*idx+2],2);
				sigPower += (feIQdata[2*j+1]*feIQdata[2*j+1])  + (feIQdata[2*j+2]*feIQdata[2*j+2]);
			}
			sigPowerNormFactor =  sqrt((float)L / sigPower);
			k=0; 
			iMean = 0.0f; qMean = 0.0f;
			for(j=segStart ; j< segEnd ; j++)	{
			// 0715 edited
				AMR_IQsegmentCC[2*k+1] = feIQdata[2*j+1] * sigPowerNormFactor;  //Inphase components, real
				AMR_IQsegmentCC[2*k+2] = feIQdata[2*j+2] * sigPowerNormFactor;  //Quadrature components, imaginary

				iMean  += AMR_IQsegmentCC[2*k+1]; //Inphase components, real
				qMean += AMR_IQsegmentCC[2*k+2]; //Quadrature components, imaginary

				k++;
			}

			iMean  /= L;
			qMean /= L;
			for(j = 0; j < L ; j++){
				AMR_IQsegmentCC[2*j+1] 	-= iMean;     //Inphase components, real
				AMR_IQsegmentCC[2*j+2] 	-= qMean;    //Quadrature components, imaginary
			}
			
			///////////////////////////////////////////////////////////////////////////////
			//  Cyclostationarity
			///////////////////////////////////////////////////////////////////////////////
			AMR_absCC20[i] = cm(AMR_IQsegmentCC, NFFT, 2, 0);
			AMR_absCC40[i] = cm(AMR_IQsegmentCC, NFFT, 4, 0); 		
			AMR_absCC40[i] -= 3.0f*AMR_absCC20[i]*AMR_absCC20[i]; 	
			// cc40 = CM40 - 3*CC20^2;
		//		cc40		cc61	   cc80	     cc82   in MATLAB
		//		0.2446      2.7188    2.4560   22.7849
		//		0.0907      1.5731    0.3923    9.5382

		//		0.144       1.76      3.17      9.4358  in DSP      at 150113
		// 		0.014       0.91      2.57      6.69

		//		cc40		cc61	   cc80	     cc82
		//		0.0625      0.0814    0.2342    0.0495  in MATLAB with test.h
		//		0.0624      0.081     0.234     0.0494  in DSP		at 150116

		} // for(i=0; i<numSeg; i++)

		// Exclude outliers		
		AMR_absCC20[0]		= dataConditining(AMR_absCC20, numSeg); // edited at 150407
		AMR_absCC40[0]		= dataConditining(AMR_absCC40, numSeg); // edited at 150420
	
		// Classification modulation type: ASK vs PSK, QAM
		condition1 = (signed int)(AMR_lineLMSerror[0]*10000.0f);
		condition2 = (signed int)(AMR_LMS_ERROR_THRESHOLD*10000.0f);
		#if AMR_DEBUG
			printf("\r (AMR_lineLMSerror= %d < %d) ? ASK : PSK,QAM  \r\n",condition1, condition2);		
		#endif

		if (condition1 < condition2){
			// ASK
			AMR_modNum	  = 2;
			AMR_modOrder  = 2;
		}
		else 
		{
			condition1 = (signed int)(AMR_absCC20[0]*10000.0f);
			condition2 = (signed int)(ABSCC20_THRESHOLD*10000.0f);
			#if AMR_DEBUG
				printf("\r (AMR_absCC20 %d < %d) ? 4,8PSK,QAM : BPSK \r\n",condition1,condition2);			
			#endif
			// Classification modulation type: BPSK vs 4,8 PSK, QAM
			if (condition1 > condition2){		
				AMR_modNum   = 3; 	AMR_modOrder = 2;	// BPSK 
			}
			else	
			{	// QPSK,8PSK,QAM			
				// refer to SigmaForLinearMod.fig at 150420    
				condition1 = (signed int)(AMR_sigma_a[0]*10000.0f);
				condition2 = (signed int)(AMR_QAMvsPSK_THRESHOLD*10000.0f);	//threshold
				#if AMR_DEBUG
					printf("\r (AMR_sigma_a %d > %d) ? QAM : 4,8PSK \r\n",condition1,condition2);			
				#endif

				if (condition1 > condition2 ){      // 0.37 --> 0.37^2 = 0.1369 variance of inst. amplitude                
					AMR_modNum = 5; AMR_modOrder = 32;    // QAM 
				}
				else
				{   // refer to CC40 forMPSK.fig at 150420      
					condition1 = (signed int)(AMR_absCC40[0]*10000.0f);
					condition2 = (signed int)(AMR_ABSCC40_THRESHOLD*10000.0f);
					#if AMR_DEBUG
						printf("\r (AMR_absCC40 %d > %d) ? QPSK : 8PSK \r\n",	condition1,condition2);			
					#endif

					if (condition1 > condition2)  AMR_modNum = 3, AMR_modOrder = 4;	// QPSK
					else						  AMR_modNum = 3, AMR_modOrder = 8;	// 8PSK
				}                           
			}// end -- if (AMR_absCC20[0] > ABSCC20_THRESHOLD)
		}// end -- if (AMR_lineLMSerror[0] < AMR_LMS_ERROR_THRESHOLD
	}// end --	if (AMR_digitalFlag && AMR_linearFlag)
}

#if AMR_Cplus	
	void AMRalgorithm::displayResults()
#else
	void displayResults()
#endif
{
/******************************************************************************/
// \ DESCRIPTION
//  - Display the results of AMR
//
// \ INPUT ARGUMENTS (global)
//	> AMR_modOrder	   
//	> AMR_modNum		
//
// \ Author(s) : AWH
// 
/******************************************************************************/
	#if AMR_DEBUG
		printf("\r\n #displayResults \r\n");
		// Description of modNum
		// Unknown -> 0 // FM   → 1 // ASK  → 2     // PSK  → 3
		// FSK  → 4    // QAM  → 5 // MSK,GSMK → 6 // AM   → 7

		// If major class is not clear, assign new class. idea from sim hong suck
		printf("\r Modulation type is ");
		switch(AMR_modNum){
			case 0: printf(" Unknown \r\n");			 break;
			case 1: printf(" FM \r\n");				     break;
			case 2: printf(" %dASK \r\n",AMR_modOrder);  break;
			case 3: printf(" %dPSK \r\n",AMR_modOrder);  break;
			case 4: printf(" %dFSK \r\n",AMR_modOrder);  break;
			case 5: printf(" QAM \r\n",AMR_modOrder);	 break;
			case 6: printf(" (G)MSK \r\n");				 break;
			case 7: printf(" AM \r\n");  	             break;		
		}
	
	#endif
}

#if AMR_Cplus	
	void AMRalgorithm::awgn(float* inAWGN,float reqSNR, float calNoisePower)
#else
	void awgn(float* inAWGN,float reqSNR, float calNoisePower)
#endif

{
/******************************************************************************/
// \ DESCRIPTION
//  - Add white Gaussian noise to input signal
//
// \ INPUT ARGUMENTS
//	> inAWGN	     - input IQ samples 
//	> reqSNR	     - required SNR
//  > AMR_noisePower - noise power of input IQ samples
//                     refer to coarseSNRestimation()
//
// \ OUTPUT ARGUMENTS
//	> inAWGN		- overwrite results to input
//
// \ Author(s) : AWH
// 
/******************************************************************************/
	unsigned int i = 0, j = 0;
	//float sigPower = 0.0, tmpNorm =0.0;
	float tmpf=0.0, sum=0.0, noiseDeviation = 0.0;
	
	static long seed=AMR_SEED;	// for white Gaussian noise

	#if AMR_DEBUG
		printf("\r\n #AWGN \r\n");
	#endif

	//for(i=0; i<AMR_LEN; i++){
	//	tmpNorm = inAWGN[2*i+1]*inAWGN[2*i+1]  + inAWGN[2*i+2]*inAWGN[2*i+2];
	//	sigPower += tmpNorm;
	//}

	//sigPower = (sigPower/ (float)AMR_LEN );
	//reqSNR = pow(10.0f,reqSNR/10.0f);
	//noisePower = sigPower / reqSNR; //e+8
	tmpf = sqrt( calNoisePower / 2.0f);

	noiseDeviation = 1;
	for(i=0; i<AMR_LEN; i++) {	
		// generate noise with zero-mean and unit variance
		sum=0.0;
		for(j=0; j<12; j++)				    //분산이 1/12인 process를 12번 더해서 분산을 1로 만듦
		{
			seed = (seed*AMR_J+1L)%AMR_MM;  //seed 대신0에서1까지의 rand() 함수사용가능
 			sum += (float)seed;
		}
		sum /= (float)AMR_MM;
		sum = sum - 6.0f;
		sum *= noiseDeviation;

		// add noise to real of signal
		inAWGN[2*i+1] += tmpf*sum;

		// generate noise with zero-mean and unit variance
		sum=0;
		for(j=0; j<12; j++)         
		{
			seed=(seed*AMR_J+1L)%AMR_MM;  //seed 대신0에서1까지의 rand() 함수사용가능
 			sum += (float)seed;
		}
		sum /= (float)AMR_MM;
		sum  = sum - 6.0f;
		sum *= noiseDeviation;

		// add noise to imag of signal
		inAWGN[2*i+2] += tmpf*sum;
	}

}

#if AMR_Cplus	
	inline float AMRalgorithm::mod(float* inX,float inY)
#else
	float mod(float* inX,float inY)
#endif
{
	float n = 0;
	float out = 0.0;

	if( inY != 0)
		n = floor(*inX / inY);  // implict type casting
	else{
		//error
	}
	out = *inX - n*inY;
	return out;
}


#if AMR_Cplus	
	inline float AMRalgorithm::my_round(float* inRound)
#else
	float my_round(float* inRound)
#endif
{
	double fractpart, intpart;

	fractpart = modf((double)*inRound, &intpart);
	if (fractpart > 0.0f){
		if (fractpart > 0.5f) intpart += 1.0f;
	}
	else if (fractpart < 0.0f)
	{
		if (fractpart < -0.5f) intpart -= 1.0f;
	}else {

	}
	return (float)intpart;
	//return (*inRound + 0.5f * sign(*inRound));
}

//
// Digital modulation signal demodulation Process
//

#if VISUALSTUDIO	

#if AMR_Cplus
	inline void AMRalgorithm::AGC(float* inAGC,int updatePeriod)
#else
	void AGC(float* inAGC,int updatePeriod)
#endif
{
	int p=0, q=0, k=0, begin=0, end=0;
	// for AGC
	int ref = 1;	// ReferenceLevel
	int L = AMR_LEN;
	int indices = 0;
	int numSubFrames = L/updatePeriod;

	float MaximumGain = 30.0f;
	float MinGain = 0.00001f;
	float mean = 0.0;
	float g = 1.0f;
	float e = 0.0;
	float K = 0.1f;	 //StepSize

	//////////////////////////////////////////////////////////////////
	// Automatic Gain Control
	//////////////////////////////////////////////////////////////////
	#if AMR_DEBUG
		printf("\r\n #AGC \r\n");
	#endif
	for(p=0; p<=numSubFrames-1; p++){
		mean = 0;
		begin = p*updatePeriod;
		end = begin + updatePeriod;
		for(q=begin+1; q<=end; q++){ 
			inAGC[2*q-1] *= g;
			inAGC[2*q]   *= g;
			mean += myCabs(inAGC[2*q-1], inAGC[2*q]);
		}

		mean /= updatePeriod;
		e = (float)ref - mean;
		g += K * e;
		
		if (g < MinGain){
			g = MinGain;
		}
		else if (g > MaximumGain){
			g = MaximumGain;
		}
	}

	for(p=end; p<L; p++){
		inAGC[2*p+1] = 0.0;
		inAGC[2*p+2] = 0.0;
	}
}

#if AMR_Cplus	
	inline void AMRalgorithm::coarseFreqOffsetCompensation(float* inIQ,float Fs,float Rs)
#else
	void coarseFreqOffsetCompensation(float* inIQ,float Fs,float Rs)
#endif
{
	int p=0, q=0, k=0, begin=0, end=0;
	int L = AMR_LEN;
	int tmp = 0;
	int numSubFrames = 0;
	// CFO est.
	int NFFT = 0;
	signed int maxIdx = 0;
	float e=0.0;
	float desireFreqResolution = 0.0;
	float maximum = 0;
	float maxVal = 0.0;
	float val = 0.0;
	float cosArg = 0.0;
	float zc = 0.0, zs = 0.0;
	float realPart = 0.0;
	float tmpf = 0.0, tmpf2 = 0.0;
	float twoPI = AMR_TWOPI;
	//////////////////////////////////////////////////////////////////
	// Coarse frequency offset compensation
	//////////////////////////////////////////////////////////////////
	#if AMR_DEBUG
		printf("\r\n #CFO est. \r\n");
	#endif

	desireFreqResolution = (Rs/100) / 10; // 심볼율 대비 10% 이내의 옵셋
	NFFT = 512;
	for(p=1; p<=500; p++){
	   NFFT= NFFT * 2;
	   if ((Fs/desireFreqResolution) < NFFT){
			break; 
	   }
	}

	//fourth-order nonlinearity method
	numSubFrames = L / NFFT;	// the number of segments
    for(p=0; p<=numSubFrames-1; p++){
		begin = p*NFFT;
		end = begin + NFFT;
		k=1;
		for(q=begin+1; q<=end; q++){ 
			AMR_IQbufferCC[2*k-1] = inIQ[2*q-1];
			AMR_IQbufferCC[2*k]	= inIQ[2*q];
			k++;
		}

		nonlinearTransform(AMR_IQbufferCC, NFFT,4,0);		

		four1(AMR_IQbufferCC,NFFT,1);

		tmpf = 0.0;
		for(q=0; q<NFFT; q++){
			// tmp4 = tmp2 + tmp3; <<-- cause numerical error
			tmpf  =  (AMR_IQbufferCC[2*q+1]*AMR_IQbufferCC[2*q+1] +  // e+1 order, caution overflow
					AMR_IQbufferCC[2*q+2]*AMR_IQbufferCC[2*q+2]);
			AMR_spectrumMagnitudeCC[q] += tmpf;    //reuse, e-3 ~ e+4 order 
			// AMR_spectrumMagnitudeCC[j] += myCabs(AMR_IQbufferCC[2*j+1],AMR_IQbufferCC[2*j+2]);    //reuse, e-3 ~ e+4 order 
		}
	}
	// 15.05.27
	//tmp = NFFT >> 1;
	//for(p=1; p<=(tmp-1); p++){   // except DC
	//	AMR_SWAP_FLOAT(AMR_spectrumMagnitudeCC[2*p+1], AMR_spectrumMagnitudeCC[NFFT - (2*p-1)]);  // real
	//	AMR_SWAP_FLOAT(AMR_spectrumMagnitudeCC[2*p+2], AMR_spectrumMagnitudeCC[NFFT - (2*p-2)]);  // imag
	//}	

	// get average spectrum
	for(p=0; p<NFFT; p++) {
		AMR_Sxx_PSD[p] /= numSubFrames;	// e-5 order
	}
	maxVal = AMR_spectrumMagnitudeCC[0];
	maxIdx = 0;
	for(q=1; q<NFFT; q++){ 
		val = AMR_spectrumMagnitudeCC[q];
		if( val > maxVal){
			maxVal = val;
			maxIdx = q;
		}
	}

	maxIdx /= 4;
	
	if (maxIdx > (NFFT>>1)){ 
		// negative axis
		maxIdx -= NFFT;
	}

	// shift banwidth sample index
	// 1850 _ 2246 in DSP    140109
	// 1844 _ 2252 in MATLAB
	// for(i=0; i<2; i++) AMR_BWSampleIdx[i] -= shiftIdx; //  not neccessary in C, but need for MATLAB 

	// freqOffsetCompensation
	for(p=0; p<L; p++){
		cosArg = -p*twoPI*maxIdx / NFFT;   // e-3order,   +- 0.002 error compared to MATLAB
		zc	   = cos(cosArg);		zs = sin(cosArg);
		realPart = (inIQ[2*p+1] * zc) - (inIQ[2*p+2] * zs);		// real
		inIQ[2*p+2] = (inIQ[2*p+1] * zs) + (inIQ[2*p+2] * zc);      // imag
		inIQ[2*p+1] = realPart;
	}
}

#if AMR_Cplus	
	inline int AMRalgorithm::samplingRateConvertor(float* inSRC,float Fs, float Rs)
#else
	int samplingRateConvertor(float* inSRC,float Fs, float Rs)
#endif
{
	// down-sampling
	int p=0;
	int L = AMR_LEN;
	int dFactor = 0;

	//////////////////////////////////////////////////////////////////
	// sampling-rate conversion process
	//////////////////////////////////////////////////////////////////
	dFactor = (Fs / Rs)/2; // implict type casting
	if (dFactor == 0) { dFactor = 7; } // exception handling 
	L /= dFactor; L++;	// modified the number of samples
	for(p=0; p<=L-1; p++){
		inSRC[2*p+2] = inSRC[2*(p*dFactor)+2];
		inSRC[2*p+1] = inSRC[2*(p*dFactor)+1];
	}	
	return L;
}

#if AMR_Cplus	
	inline int AMRalgorithm::fineCarrierOffsetNSymbolTimingRecovery(float* inIQ,int L)
#else
	int fineCarrierOffsetNSymbolTimingRecovery(float* inIQ,int L)
#endif
{
	//////////////////////////////////////////////////////////////////
	// Fine carrier frequency offset and Symbol timing recovery
	//////////////////////////////////////////////////////////////////

	int p,q;
// carrier freq offset est.
	int DigitalSynthesizerGain = -1;
	// Look into model for details for details of PLL parameter choice. 
	// Refer equation 7.30 of [1] 
	// PhaseErrorDetectorGain(K_p) for Fine Frequency Compensation PLL, determined by 2KA^2 (for binary PAM),
	// QPSK could be treated as two individual binary PAM
	int aK = 1;
	// int PhaseErrorDetectorGain = 1; // 8PSK
	int PhaseRecoveryGain = 1; // K_0 for Fine Frequency Compensation PLL
	//int PhaseRecoveryGain = 2; // equal to the number of samples per symbol
	float cosArg = 0.0;
	float zc = 0.0, zs = 0.0;
	float realPart = 0.0;
	float pPhase = 0.0;
	float e=0.0;
	float DDSOut = 0.0;	
	float phErr = 0.0;
	float loopFiltOut = 0.0;			
	float A = 1/sqrt(2.0f);	//Amplitude
	float PhaseErrorDetectorGain = 4*aK*A*A; //% QPS	
	float PhaseRecoveryLoopBandwidth = 0.005f; // Normalized loop bandwidth for fine frequency compensation
	float PhaseRecoveryDampingFactor = 1.0f;     // Damping Factor for fine frequency compensation
	int PostFilterOversampling = 2;
	// Refer C.57 to C.61 in [1] for K1 and K2
	float theta = PhaseRecoveryLoopBandwidth / (PhaseRecoveryDampingFactor +
		0.25f/PhaseRecoveryDampingFactor)/PostFilterOversampling;
	float d = 1 + 2*PhaseRecoveryDampingFactor*theta + theta*theta;
	float ProportionalGain = (4*PhaseRecoveryDampingFactor*theta/d) / 
		(PhaseErrorDetectorGain*PhaseRecoveryGain);
	float IntegratorGain = (4*theta*theta/d) / (PhaseErrorDetectorGain*PhaseRecoveryGain); // e-6 order
	float loopFiltOutFreq=0;

	//%% Parameter-Timing recovery
	
	int N = 2;  //% samples per symbol
	int pDelayStrobe = 0;
	int pStrobe = 0;
	int underflow = 0;
	int pCount = 1; // % timing loop output 
	signed int k0 = -1; //% the controller is a decrementing modulo-1 counter

	float pRegTemp = 0;
	float pNCODelay = 0;
	float counter = 0;		// NCO
	float pAlpha = 0.5;
	float interpFiltOut = 0;
	float noiseBandwidth = 0.01f; //% BnTs, normalized to the bit rate
	float dampingFactor = 1; //%1/sqrt(2); % Zeta
	float phaseDetectorGain = 2.7f; //% kp % time to approach the steady state
	float cc = noiseBandwidth/(dampingFactor + (1/(4*dampingFactor)));
	//float den = (cc*4*dampingFactor/N);
	float nom = (1.0f+ cc*2*dampingFactor/N + (cc/N)*(cc/N));
	float ProportionalGainTiming = ((cc*4*dampingFactor/N)/nom) /(phaseDetectorGain*k0);
	float den2 = (cc*cc*4)/(N*N);
	float IntegratorGainTiming = (den2/nom) /(phaseDetectorGain*k0);
	float pDelayBuffer1 = 0.0, pDelayBuffer2 = 0.0, pDelayBuffer3 = 0.0;    
	float pTEDDelay1 = 0.0, pTEDDelay2 = 0.0 ; 
	float pMU = 0;
	float timingLoopOut = 0;
	float Delta = 0.0;
	float loopFiltOutTiming = 0.0;
	float floatReminder = 0.0;
	// Carrier Phase recovery
	float re = 0.0;
	float im = 0.0;
	float phEst = 0.0;

	loopFiltOutFreq = 0.0; DDSOut = 0.0; loopFiltOutTiming = 0.0; pCount = 0;
	for(p=0; p<=L-1; p++){
		// Carrier offset recovery  - 2 samples per symbol   
		// Complex phase shift
		zc	   = cos(pPhase);		zs = sin(pPhase);   // e-3order,   +- 0.002 error compared to MATLAB
		AMR_fineCompOut[0] = (inIQ[2*p+1] * zc) - (inIQ[2*p+2] * zs);		// real
		AMR_fineCompOut[1] = (inIQ[2*p+1] * zs) + (inIQ[2*p+2] * zc);      // imag
	
		//// Find phase error
		//if isDQPSK        
		//	%for 8PSK
		//	% http://kr.mathworks.com/help/comm/ref/comm.carriersynchronizer-class.html#bulwrmb
		//	a = sqrt(2)-1;
		//	if re >= im
		//		phErr =  a*sign(re) *im - sign(im)*re; 
		//	else
		//		phErr = sign(re)*im - a*sign(im)*re;
		//	end
		//else
		//	% QAM and QPSK
		//	phErr = sign(re)*im - sign(im)*re;
		phErr = (sign(AMR_fineCompOut[0])*AMR_fineCompOut[1]) - (sign(AMR_fineCompOut[1])*AMR_fineCompOut[0]);
		//end
		// caution initialization  
		loopFiltOutFreq += (phErr*IntegratorGain);        // loop filter output

		//% Direct Digital Synthesizer
		DDSOut += (phErr*ProportionalGain + loopFiltOutFreq); //% Integrator
		pPhase =  DigitalSynthesizerGain * DDSOut;

		// Timing recovery  - ZCTED, 2 samples per symbols       
		// Parabolic interpolator - Farrow structure, p 471 in [1]
		for(q=0; q<=1; q++){
			AMR_interpFiltOut[q] = AMR_pDelayBuffer2[q] +        
				pMU*(-pAlpha*(AMR_fineCompOut[q]+AMR_pDelayBuffer2[q]+AMR_pDelayBuffer3[q])+
				(1+pAlpha)*AMR_pDelayBuffer1[q])- pAlpha *
				(AMR_pDelayBuffer1[q]+AMR_pDelayBuffer2[q]-AMR_fineCompOut[q]-AMR_pDelayBuffer3[q])
				*pMU*pMU;

			//% Update delay buffers
			AMR_pDelayBuffer3[q] = AMR_pDelayBuffer2[q];
			AMR_pDelayBuffer2[q] = AMR_pDelayBuffer1[q];
			AMR_pDelayBuffer1[q] = AMR_fineCompOut[q];   
		}
		
		// Timing Error Detector
		//% e = TED(obj,interpFiltOut);
		if (pStrobe && (pDelayStrobe != pStrobe)){
			e = AMR_pTEDDelay1[0] * (sign(AMR_pTEDDelay2[0]) - 
				sign(AMR_interpFiltOut[0])) + AMR_pTEDDelay1[1] * 
				(sign(AMR_pTEDDelay2[1]) - sign(AMR_interpFiltOut[1]));
		}else{
			e = 0.0f;
		}     

	
		if (pDelayStrobe != pStrobe){
			//% Shift contents in delay register
			for(q=0; q<=1; q++){
				AMR_pTEDDelay2[q] = AMR_pTEDDelay1[q];
				AMR_pTEDDelay1[q] = AMR_interpFiltOut[q];
			}
		}else if (pStrobe){
			//% Two consecutive high strobes
			for(q=0; q<=1; q++){
				AMR_pTEDDelay2[q] = 0.0f; // Stuff missing sample
				AMR_pTEDDelay1[q] = AMR_interpFiltOut[q];
			}	
		}
		//% If both current and previous enable signals are 0, skip current sample
		//% and keep the delayed signals unchanged. (Bit stripping)
		pDelayStrobe = pStrobe;   

		//% Loop filter
		loopFiltOutTiming = loopFiltOutTiming + IntegratorGainTiming*e;    // integrator componenet of loop filter

		//% Updates timing difference for Interpolation Filter
		//% [underflow,pMU] = NCO_control(obj,e,loopFiltOut);...
		Delta = e*ProportionalGainTiming + loopFiltOutTiming; //% If loop is in lock, Delta would be small
		floatReminder = mod(&pNCODelay,1); // mod(pNCODelay,1)
		
		//tmpf = Delta+(1.0f/(float)PostFilterOversampling);
		
		counter = floatReminder-Delta-(1.0/(double)PostFilterOversampling); // decrementing counter
		
		//if( fabs(counter) < 0.0001f) counter = 0.0f;
		if (counter < 0){
			underflow = 1;
			pRegTemp =	floatReminder/(Delta+(1.0f/(double)PostFilterOversampling));
		}
		else{
			underflow = 0;
		} 
		//tmpf = (1.0f/(float)PostFilterOversampling) + Delta;  // NCO control word
		//counter = pNCODelay - tmpf;             //% update counter value for next cycle
		//if (counter<0){               //   % test to see if underflow has occurred
		//	counter = 1 + counter;    //% reduce counter value modulo-1 if underflow
		//	pStrobe = underflow;  
		//	underflow = 1;              //% set underflow flag
		//	pRegTemp = pNCODelay/tmpf;            //% update mu
		//}else{
		//	pStrobe = underflow;
		//	underflow = 0;
		//	pRegTemp = pMU;
	 //  }   

		pMU = pRegTemp;
		//  % Update delay buffer
		pNCODelay = counter;  
		if (pStrobe>0){		
		   AMR_timingLoopOut[2*pCount+1] = AMR_interpFiltOut[0];	// real
		   AMR_timingLoopOut[2*pCount+2] = AMR_interpFiltOut[1];	// imag
		   pCount++;
		}

		pStrobe = underflow;   

	}
	L = pCount;
	
	return L;
/*#if VISUALSTUDIO
	fa=fopen("real.txt","w");
	fprintf(fa,"data =[");
	for(p=0; p<L; p++) fprintf(fa,"%0.4f,...\n",(AMR_timingLoopOut[2*p+1]));
	fprintf(fa,"]; \n");
	fprintf(fa,"figure; plot(data)");
	fclose(fa);			
#endif	
#if VISUALSTUDIO
	fa=fopen("imag.txt","w");
	fprintf(fa,"data2 =[");
	for(p=0; p<L; p++) fprintf(fa,"%0.4f,...\n",(AMR_timingLoopOut[2*p+2]));
	fprintf(fa,"]; \n");
	fprintf(fa,"figure; cc=data+1i*data2; scatterplot(cc)");
	fclose(fa);			
#endif*/	
}

#if AMR_Cplus	
	inline void AMRalgorithm::carrierPhaseRecovery(float* inCPR,int L,int updatePeriod)
#else
	void carrierPhaseRecovery(float* inCPR,int L,int updatePeriod)
#endif
{
	int p=0, q=0, k=0;
    int numSubFrames = 0;
	int begin=0, end =0;

	float re, im, tmpf, phEst, zc, zs;

	re=im=tmpf=phEst=zc=zs=0.0;

	//////////////////////////////////////////////////////////////////
	// Carrier phase recovery
	//////////////////////////////////////////////////////////////////
    numSubFrames = L/updatePeriod; // implict type casting


	AMR_modOrder = 4; ///////////////////////////////////// tmp

 	for(p=0; p<=numSubFrames-1; p++){

		begin = p*updatePeriod;
		end = begin + updatePeriod;
		k=1;
		for(q=begin+1; q<=end; q++){ 
			AMR_IQbufferCC[2*k-1] = inCPR[2*q-1];
			AMR_IQbufferCC[2*k]	= inCPR[2*q];
			k++;
		}
		
		nonlinearTransform(AMR_IQbufferCC,updatePeriod,AMR_modOrder,0);
		re = 0.0f; im = 0.0f;
		for(q=1; q<=updatePeriod; q++){ 
			re += AMR_IQbufferCC[2*q-1];
			im += AMR_IQbufferCC[2*q];
		}	
		tmpf = atan2(im,re);
		phEst = (1.0f/(float)AMR_modOrder) * tmpf; // in degree
		zc	   = cos(phEst);		zs = sin(phEst);   // e-3order,   +- 0.002 error compared to MATLAB

		for(q=begin; q<end; q++){			
			re = (inCPR[2*q+1] * zc) + (inCPR[2*q+2] * zs);		// real
			inCPR[2*q+2] = (inCPR[2*q+2] * zc) - (inCPR[2*q+1] * zs);      // imag
			inCPR[2*q+1] = re;
		}
	}


//#if VISUALSTUDIO
//	fa=fopen("scatterPlot.txt","w");
//	fprintf(fa,"data =[");
//	for(p=0; p<L; p++) fprintf(fa,"%0.4f + 1i*%0.4f,...\n",(AMR_timingLoopOut[2*p+1]),(AMR_timingLoopOut[2*p+2]));
//	fprintf(fa,"]; \n");
//	fprintf(fa,"figure; plot(data,'.')");
//	fclose(fa);			
//#endif	
}

#if AMR_Cplus	
	inline void AMRalgorithm::symbolDecision(float* inIQ,int L,int modOrder)
#else
	void symbolDecision(float* inIQ,int L,int modOrder)
#endif
{
	int p;
	float tmpf;

	// normalization factor to convert from PI-domain to linear domain
	float normFactor = (float)modOrder / (float)AMR_TWOPI; 

    // convert input signal angle to linear domain; round the value to get ideal
    //  constellation points 
	for(p=0; p<L; p++){
		tmpf = atan2(inIQ[2*p+2],inIQ[2*p+1]);
		tmpf *= normFactor;
	    inIQ[p] = my_round( &tmpf );
		// move the negative integers by Modulation order
		if (inIQ[p] <0 ) inIQ[p] += modOrder; 
	}
	// how to handle -0.0?

}
#if AMR_Cplus	
	inline void AMRalgorithm::PSKdemodulation(float* inPSKdemod)
#else
	void PSKdemodulation(float* inPSKdemod)
#endif
{
	/********************************************************************************/
	// \ DESCRIPTION
	//  - Normalize signal power to 1
	//
	// \ INPUT ARGUMENTS
	//	> inIQ	     - input IQ samples 
	//	> L			 - the number of input samples
	//
	// \ OUTPUT ARGUMENTS
	//	> inIQ		- overwrite results to input
	//
	// \ Author(s) : AWH
	// 
	// \ Reference
	//  > [1] Michael Rice, "Digital Communications - A Discrete-Time Approach"
	/******************************************************************************/

	int p=0, q=0, k=0, begin=0, end=0;
	int sigLen = 0;		

	AMR_samplingFrequency = 546875.0f;
	AMR_coarseSymbolRate = 35000.0f;
	AMR_modOrder = 4;

	AGC(inPSKdemod,100);
	
	coarseFreqOffsetCompensation(inPSKdemod,AMR_samplingFrequency,AMR_coarseSymbolRate);
	 		
	sigLen = samplingRateConvertor(inPSKdemod,AMR_samplingFrequency,AMR_coarseSymbolRate);

	sigLen = fineCarrierOffsetNSymbolTimingRecovery(inPSKdemod,sigLen);
	 
	carrierPhaseRecovery(AMR_timingLoopOut,sigLen,100);
	
	symbolDecision(AMR_timingLoopOut,sigLen,AMR_modOrder);

}

#if AMR_Cplus
	inline unsigned int AMRalgorithm::FMdemodulation(float* inFMdemod, unsigned int sigLen)
#else
	unsigned int FMdemodulation(float* inFMdemod, unsigned int sigLen)
#endif
{
	unsigned int j,i,k;
	float freqdev = 20000.0f;  //% 20kHz, Default frequency deviation
	float cc = 0.0;
	float Fs;
	float threshold =0.0f;
	float mean = 0.0, ss=0.0, ep=0.0, var=0.0, sdev = 0.0; //statistics
	Fs = 136720.0f;
	cc = (1.0f/((float)AMR_TWOPI*freqdev));

    // %% FM demodulation using frequency discriminator
	for(j=0; j<=sigLen-1; j++){
			// maximum +- 0.6 error compared with MATLAB
		AMR_FMdemodiPhase[j] = atan2(inFMdemod[2*j+2],inFMdemod[2*j+1]); // wraped instantaneous phase
	}
	for(j=0; j<=sigLen-1; j++) AMR_FMdemodCk[j] = 0.0f;

	// make range from -PI to +PI (=unwrap)
	for(j=1; j<=sigLen-1; j++){
		if( AMR_FMdemodiPhase[j] - AMR_FMdemodiPhase[j-1] > (float)AMR_PI)
			AMR_FMdemodCk[j] = AMR_FMdemodCk[j-1] - (float)AMR_TWOPI;
		else if ( AMR_FMdemodiPhase[j] - AMR_FMdemodiPhase[j-1] < (-1.0f)*(float)AMR_PI)
			AMR_FMdemodCk[j] = AMR_FMdemodCk[j-1] + (float)AMR_TWOPI;
		else	AMR_FMdemodCk[j] = AMR_FMdemodCk[j-1];
	}

	for(j=0; j <=sigLen-1 ;j++){
		AMR_FMdemodiPhase[j] += AMR_FMdemodCk[j]; // unwraped instantaneous phase
		// instantaneous frequency = differentiation of instantaneous phase
		// length is reduce to L-1
		if (j>0){
			AMR_FMdemodiFreq[j-1] = AMR_FMdemodiPhase[j]-AMR_FMdemodiPhase[j-1];
			AMR_FMdemodiFreq[j-1] *=  cc*Fs;
		}
	}
	AMR_FMdemodiFreq[0] = AMR_FMdemodiFreq[1];
	AMR_FMdemodiFreq[sigLen-1] = AMR_FMdemodiFreq[sigLen-2]; 

	//if(AMR_numTOA > 10){
		//
		// Square envelop detection for hopping signal
		//
	for(j=0; j <=sigLen-1 ;j++){
		AMR_FMdemodiFreq[j] = AMR_FMdemodiFreq[j] * AMR_FMdemodiFreq[j];
	}
	//}
	//
	// Hard limiter
	//
	mean=0.0f;
	for(i=0; i<=sigLen-1; i++){
		mean += AMR_FMdemodiFreq[i];
	}			
	mean /= sigLen;

	ep = 0.0; var=0.0;
	for(i=0; i<=sigLen-1; i++){
		AMR_FMdemodiFreq[i] -=  mean;
		ss = AMR_FMdemodiFreq[i];
		ep += ss;
		var += ss*ss;
	}
	var = (var-ep*ep/sigLen)/(sigLen-1);		//corrected two-pass formula
	sdev = sqrt(var);
	threshold =  3.0f*sdev;		
	j=0;
	for(i=0; i<=sigLen-1; i++){
		if(AMR_FMdemodiFreq[i] <= threshold && AMR_FMdemodiFreq[i] >= -threshold ){
			AMR_FMdemodiFreqBuff[j] = AMR_FMdemodiFreq[i];
			j++;
		}
	}
	
	sigLen = j;
	return (sigLen);
}

#if AMR_Cplus
	inline void AMRalgorithm::AMdemodulation(float* inAMdemod)
#else
	void AMdemodulation(float* inAMdemod)
#endif
{
	int j;
	int L = AMR_LEN;
	float freqdev = 20e3;  //% Default frequency deviation
	float cc = 0.0;
	float Fs;
	float m_a=0.0, mean = 0.0;
	AMR_samplingFrequency = 68359;
	Fs = AMR_samplingFrequency;

	for(j=0 ; j<=L-1 ; j++){
		AMR_iAmplitude[j] = myCabs(inAMdemod[2*j+1],inAMdemod[2*j+2]);     // instantaneous amplitude
		m_a += AMR_iAmplitude[j] / L;        // mean of instant. amp.  ,  0.89 in MATLAB, 0.95 in DSP
	}
	
	for (j = 0; j<=L-1; j++){
		AMR_iAmplitude[j] -= m_a;
	}

}


//
// Make *.wav file
//

#if AMR_Cplus
	inline void AMRalgorithm::write_little_endian(unsigned int word, int numBytes, FILE *wavFile)
#else
	void write_little_endian(unsigned int word, int numBytes, FILE *wavFile)
#endif
{
	// Reference
	// http://blog.daum.net/sualchi/13720173
	//

	unsigned buf;
	while(numBytes > 0)
	{
		buf = word & 0xff;
		fwrite(&buf, 1, 1, wavFile);
		numBytes--;
		word >>= 8;

	}
}

#if AMR_Cplus
	inline void AMRalgorithm::write_float_little_endian(float word, int numBytes, FILE *wavFile)
#else
	void write_float_little_endian(float word, int numBytes, FILE *wavFile)
#endif
{
	// Reference
	// http://blog.daum.net/sualchi/13720173
	//

	fwrite(&word, sizeof(float), 1, wavFile);
}

#if AMR_Cplus
	inline void AMRalgorithm::write_wav(char *filename, unsigned long num_samples, float *data, int s_rate)
#else
	void write_wav(char *filename, unsigned long num_samples, float *data, int s_rate)
#endif
{

	// Reference
	// http://blog.daum.net/sualchi/13720173
	// http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html
	//

	FILE* wav_file;
	unsigned int sample_rate;
	unsigned int num_channels;
	unsigned int bytes_per_sample;
	unsigned int byte_rate;
	unsigned long i;		//* counter for samples *//
	
	num_channels	 = 1;		// monoaural
	bytes_per_sample = 2;	// bytes/sample

	if (s_rate <=0 ) sample_rate = 132300; //44100;
	else sample_rate = (unsigned int) s_rate;

	byte_rate = sample_rate * num_channels * bytes_per_sample;

	wav_file = fopen(filename, "w");
	assert(wav_file); //* make sure it opened *//

	/* write RIFF header */
	fwrite("RIFF", 1, 4, wav_file);
	write_little_endian(36 + bytes_per_sample* num_samples*num_channels, 4, wav_file);
	fwrite("WAVE", 1, 4, wav_file);

	/* write fmt subchunk */
	fwrite("fmt ",1,4,wav_file);
	write_little_endian(16,4,wav_file);	// subChunk1Size is 16
	write_little_endian(1,2,wav_file);
	write_little_endian(num_channels,2,wav_file);
	write_little_endian(sample_rate,4,wav_file);
	write_little_endian(byte_rate,4,wav_file);
	write_little_endian(num_channels*bytes_per_sample,2,wav_file); //block align
	write_little_endian(8*bytes_per_sample, 2, wav_file); // bits/samples

	//write data subchunk
	fwrite("data", 1,4,wav_file);
	write_little_endian(bytes_per_sample * num_samples*num_channels, 4,wav_file);
	for(i=0; i<num_samples; i++){
		write_little_endian((unsigned int)(data[i]),bytes_per_sample,wav_file);
	}
	fclose(wav_file);
}

#endif // for #define VISUALSTUDIO

	
#endif // for #ifndef _AMR_FUNCTIONS

/* 
FILTER.C 
An ANSI C implementation of MATLAB FILTER.M (built-in)
Written by Chen Yangquan <elecyq@nus.edu.sg>
1998-11-11
*/
//
//#include<stdio.h>
//#define ORDER 3
//#define NP 1001
//
//
////void filter(int,float *,float *,int,float *,float *);
//filter(int ord, float *a, float *b, int np, float *x, float *y)
//{
//        int i,j;
//	y[0]=b[0]*x[0];
//	for (i=1;i<ord+1;i++)
//	{
//        y[i]=0.0;
//        for (j=0;j<i+1;j++)
//        	y[i]=y[i]+b[j]*x[i-j];
//        for (j=0;j<i;j++)
//        	y[i]=y[i]-a[j+1]*y[i-j-1];
//	}
//	// end of initial part
//
//	for (i=ord+1;i<np+1;i++)
//	{
//		y[i]=0.0;
//			for (j=0;j<ord+1;j++)
//			y[i]=y[i]+b[j]*x[i-j];
//			for (j=0;j<ord;j++)
//			y[i]=y[i]-a[j+1]*y[i-j-1];
//	}
//} // end of filter
//
//
//main()
//{
//	FILE *fp;
//	float x[NP],y[NP],a[ORDER+1],b[ORDER+1];
//	int i,j;
//
//	// printf("hello world \n"); 
//
//	if((fp=fopen("acc1.dat","r"))!=NULL)
//	{
//		for (i=0;i<NP;i++)
//		{
//			fscanf(fp,"%f",&x[i]);
//	//         printf("%f\n",x[i]); 
//		}
//	}
//	else
//	{
//		printf("\n file not found! \n");
//		exit(-1);
//	}
//	close(fp);
//
//	//  test coef from
//	// [b,a]=butter(3,30/500);  in MATLAB
//
///*	b[0]=0.0007;
//	b[1]=0.0021;
//	b[2]=0.0021;
//	b[3]=0.0007;
//	a[0]=1.0000;
//	a[1]=-2.6236;
//	a[2]=2.3147;
//	a[3]=-0.6855*/;
//
//	filter(ORDER,a,b,NP,x,y);
//	/* NOW y=filter(b,a,x);*/
//
//	/* reverse the series for FILTFILT */
//	for (i=0;i<NP;i++)
//	{ x[i]=y[NP-i-1];}
//	/* do FILTER again */
//	filter(ORDER,a,b,NP,x,y);
//	/* reverse the series back */
//	for (i=0;i<NP;i++)
//	{ x[i]=y[NP-i-1];}
//	for (i=0;i<NP;i++)
//	{ y[i]=x[i];}
//	/* NOW y=filtfilt(b,a,x); boundary handling not included*/
//
//	if((fp=fopen("acc10.dat","w+"))!=NULL)
//	{
//		for (i=0;i<NP;i++)
//		{
//			fprintf(fp,"%f\n",y[i]);
//		}
//	}
//	else
//	{
//		printf("\n file cannot be created! \n");
//		exit(-1);
//	}
//	close(fp);
//}  
//end of filter.c


// filter.m MATLAB version
/*
// http://www.mathworks.com/matlabcentral/answers/9900-use-filter-constants-to-hard-code-filter#answer_13623
	function [Y, z] = myFilter(b, a, X, z)
	% Author: Jan Simon, Heidelberg, (C) 2011
	n    = length(a);
	z(n) = 0;      % Creates zeros if input z is omitted

	b = b / a(1);  % [Edited, Jan, 26-Oct-2014, normalize parameters]
	a = a / a(1);

	Y    = zeros(size(X));
	for m = 1:length(Y)
	   Y(m) = b(1) * X(m) + z(1);
	   for i = 2:n
		  z(i - 1) = b(i) * X(m) + z(i) - a(i) * Y(m);
	   end
	end
	z = z(1:n - 1);
*/

