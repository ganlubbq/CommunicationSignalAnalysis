#ifndef _AMRCLASS
#define _AMRCLASS

//
// Object-oriented programming
//
// Notice
//  - qsort is replaced by heapsort deu to the error
//  - 
// 15.06.25

#include <stdlib.h>

class AMRalgorithm{
public:
	// default constructor
	AMRalgorithm();

public:   // memeber variable
	FILE *fa;

	int AMR_decimationTable[8];
	int AMR_DEBUG_numseg;			 //for AMR_DEBUGging
	int AMR_ADCsamplingFrequency;	 //140MHz
	int AMR_Decimation;
	signed int AMR_symRateEstSegLen;

	unsigned int AMR_existsCarrierFlag ; //for AM,CW vs others by using whether existance of the carrier component
	unsigned int AMR_digitalFlag ;		//for analog vs digital signal identification bits	
	unsigned int AMR_linearFlag  ;		//for linear vs non-linear signal identification bits
	unsigned int AMR_minBinInterval;
	unsigned int AMR_numIQFPGA;       // total number of IQ samples in *.diq file

	unsigned int AMR_SNRTable[6];
	unsigned int AMR_decimationIdx;
	unsigned int AMR_requireIQsamples;
	unsigned int AMR_maxSpecMagIdx;
	float AMR_samplingFrequency;
	float AMR_freqInterval;
	float AMR_samPeriod;	
	float AMR_DEBUG_floatTmp;	//for AMR_DEBUGging
	float AMR_noisePower;
	float AMR_SNR_dB;
	// SNR control
	float AMR_noiseFloor_dB;       

	unsigned int AMR_binCount[AMR_MAX_NUMBINS+2];	 // for histogram
	float		 AMR_binCenter[AMR_MAX_NUMBINS+2];	 // for histogram	
	float		 AMR_bins[AMR_MAX_NUMBINS+2];	// 구역 분할 범위, numAMR_bins의 중심값

	unsigned int			AMR_numPeaks;						    	// for mode operation

	/// for welch PSD and tone spacing 
	float AMR_IQSegment_PSD[(AMR_MAX_welchPSDsize<<1) + 1]; // for powerSpectrum	
	float AMR_hanningWindow_PSD[AMR_MAX_welchPSDsize];
	float AMR_Sxx_PSD[AMR_MAX_welchPSDsize];
	float AMR_powerSpectrum[AMR_MAX_welchPSDsize];
	float AMR_powerSpectrum_dB[AMR_MAX_welchPSDsize];
	float AMR_peakValues[AMR_MAX_welchPSDsize>>1];		// for PeakDetection
	float AMR_TSE_PSDbuffer[AMR_MAX_welchPSDsize];
	
	signed int  AMR_peakIndices[AMR_MAX_welchPSDsize>>1];	// for PeakDetection
	signed int AMR_freqVector[AMR_MAX_welchPSDsize];
	unsigned int AMR_IndexBuffer_PSD[AMR_MAX_welchPSDsize];
	int   AMR_TSE_indexBuffer[AMR_MAX_welchPSDsize];	
	signed int AMR_BWSampleIdx[2];
	signed int AMR_BWSampleIdx_TMP[2];

	float AMR_hammingWindow[AMR_LPFLEN];				//for LPF
	float AMR_LPFCoeff[AMR_LPFLEN];					//for LPF
	float AMR_h[AMR_LPFLEN];						        //for LPF
	float AMR_F[AMR_LPFLEN>>3];						    //for LPF
	float AMR_M[AMR_LPFLEN>>3];						    //for LPF
	float AMR_k[AMR_LPFLEN>>1];						    //for LPF
	float AMR_b[AMR_LPFLEN>>1];						    //for LPF
	float AMR_a[AMR_LPFLEN>>1];						    //for LPF

	float AMR_MovInB[AMR_symbolRate_SamplesPerSegment];			//for moving average
	float AMR_cbegin1[AMR_symbolRate_SamplesPerSegment];			//for moving average
	float AMR_cbegin2[AMR_symbolRate_SamplesPerSegment];			//for moving average
	float AMR_cend1[AMR_symbolRate_SamplesPerSegment];			//for moving average
	float AMR_cend2[AMR_symbolRate_SamplesPerSegment];			//for moving average
	float AMR_Movfiltout[AMR_symbolRate_SamplesPerSegment];		//for moving average

	////Coarse symbol rate estimation(CSRE)
	float AMR_ck[AMR_symbolRate_SamplesPerSegment];				 // Coarse symbol rate estimation
	float AMR_iPhase[AMR_symbolRate_SamplesPerSegment];			 // Coarse symbol rate estimation
	float AMR_iFrequency[AMR_symbolRate_SamplesPerSegment];	     // instantaneous frequency
	float AMR_iAmplitude[AMR_symbolRate_SamplesPerSegment];		 // linear vs non-linear classification
	float AMR_spectrumBuffer[AMR_symbolRate_SamplesPerSegment];
	float AMR_filteredSamples[AMR_symbolRate_SamplesPerSegment];		 // Coarse symbol rate estimation
	float AMR_spectrumMagnitude[AMR_symbolRate_SamplesPerSegment]; // Coarse symbol rate estimation

	float AMR_halfFreqVector[AMR_symbolRate_SamplesPerSegment>>1];  // for symbol rate estimation
	float AMR_IQsegmentSym[(AMR_symbolRate_SamplesPerSegment<<1) + 1]; // for symbol rate estimation

	 ////features
	float AMR_sigma_a[AMR_MAXNumOfFeature];
	float AMR_absCC20[AMR_MAXNumOfFeature];
	float AMR_absCC40[AMR_MAXNumOfFeature];
	float AMR_lineLMSerror[AMR_MAXNumOfFeature];
	float AMR_mu42f[AMR_MAXNumOfFeature];
	float AMR_dataConditioningTmp[AMR_MAXNumOfFeature];

	//// for cyclic cumulant
	float AMR_IQbufferCC[(AMR_CC_samplesPerSegment<<1) + 1];   // nonlinear transform 
	float AMR_IQsegmentCC[(AMR_CC_samplesPerSegment<<1) + 1];  // for cyclic cumulant
	float AMR_IQDataTmp2[(AMR_CC_samplesPerSegment<<1) + 1];   // nonlinear transform 
	float AMR_IQDataTmp1[(AMR_CC_samplesPerSegment<<1) + 1];   // nonlinear transform 
	float AMR_IQDataResQ[(AMR_CC_samplesPerSegment<<1) + 1];   // nonlinear transform 
	float AMR_spectrumMagnitudeCC[AMR_CC_samplesPerSegment]; // for cyclic cumulant

	float AMR_oneIQData[3];		// for cumulant
	float AMR_oneIQDataSec[3];		// for cumulant
	
	//// parameter estimation
	float AMR_bwThreshold_dB;     // threshold of normalized spectrum  when 0 (min) to 1 (max)
	float AMR_bwThreshold;    

	unsigned int AMR_BW_NFFT;
	unsigned int AMR_NFFT;
	unsigned int AMR_numSeg;   				// the number of segments
	unsigned int AMR_isFSK;

	unsigned int IQ_Buffer[IQ_NUM_MAX];
	unsigned int IQ_BufferConv[IQ_NUM_MAX];
	short I_BufferConv[IQ_NUM_MAX];
	short Q_BufferConv[IQ_NUM_MAX];

	 ////TOA searching
	unsigned int beginTOAidx[AMR_MAX_NUM_SEGMENT];
	unsigned int endTOAidx[AMR_MAX_NUM_SEGMENT];

	//demodulation
	float AMR_interpFiltOut[2];
	float AMR_fineCompOut[2];
	float AMR_pDelayBuffer3[2];
	float AMR_pDelayBuffer2[2];
	float AMR_pDelayBuffer1[2]; 
	float AMR_pTEDDelay1[2];  // interpolant buffer
	float AMR_pTEDDelay2[2];
	float AMR_timingLoopOut[AMR_MAXArray];

public:
	unsigned int AMR_numTOA;
	 ////struct AMR_struct { ... } 구조체를 amrObj 이라는 데이터형으로 정의한다.
	 ////http://www.tipssoft.com/bulletin/board.php?bo_table=FAQ&wr_id=968

	unsigned int IQbeginIdx[AMR_MAX_NUM_SEGMENT];
	unsigned int IQendIdx[AMR_MAX_NUM_SEGMENT];		// the last index in the IQ array
	unsigned int IQnum[AMR_MAX_NUM_SEGMENT];		// the number of IQ samples
	//AMR_struct   amrObj[AMR_MAX_NUM_SEGMENT];

	unsigned int AMR_LEN;
	unsigned int AMR_modNum;
	unsigned int AMR_modOrder;
	unsigned int AMR_samplesPerSym;
	unsigned int AMR_numSym;
	unsigned int AMR_coarseSamplesPerSymbol ;
	unsigned int AMR_coarseNumberOfsymbols  ;
	unsigned int AMR_coareseSNRIdx;
	signed int AMR_symratePeakIdx;


	float AMR_coarseBandWidth; 	
	float AMR_coarseFreqOffset;	
	float AMR_coareseSNR_dB;     
	float AMR_coarseSymbolRate;  
	float AMR_coarseToneSpacing; 
	float AMR_symratePeak;   // 15.05.27

	
	float AMR_inputData[AMR_MAXArray];
	float AMR_filteredInputData[AMR_MAXArray];
	float demodSym[AMR_MAXArray];
	
	float I_BufferConvF[IQ_NUM_MAX];
	float Q_BufferConvF[IQ_NUM_MAX];

	

	//Demodulation
	unsigned int AMR_wavBuffer[IQ_NUM_MAX];
	float AMR_FMdemodiPhase[IQ_NUM_MAX];
	float AMR_FMdemodiFreq[IQ_NUM_MAX];
	float AMR_FMdemodCk[IQ_NUM_MAX];
	float AMR_demodIQ[(IQ_NUM_MAX<<1) +1];
	float AMR_FMdemodiFreqBuff[IQ_NUM_MAX];

public:
	unsigned int LOCAL_load_IQ(char* fname);
	signed int coarsefreqOffsetEstimation(int cfoDataLen,int cfoFFTlength);
	signed int sign(float inSIGN);
	unsigned int peakDetection(float* pdInput, signed int pdInputlength, signed int pdMinPeakDist, float threshold);

	void initializeVariable();
	void awgn(float* inAWGN,float reqSNR, float calNoisePower);
	
	void unitPower(float* inIQ, unsigned int L);
	void welchPSD(float* inPSD, unsigned int inPSDlength, unsigned int FFTlength, signed int inPower, float* outPSD);
	
	void linearModClassification(float* feIQdata, int feIQdataLen, int feSamplesPerFrame);
	void coarseToneSpacing(float* ctsIQdataInput,signed int ctsIQdataLen,signed int ctsFFTlength,signed int ctsHighFFTlength, signed int offset);
	void reSampling(int resIQdataLen, float resEstSymRate);

	void linearVsNonLinearClassification(float* lvnIQ, unsigned int lvnIQLen, unsigned int csrSamplesPerFrame);
	void lowPassFiltering(float* lpfIQdata);
	
	void displayResults();

	
	signed int setFrequencyVector(int iqLen);
	void spectrumFeatureExt(float* ifeIQdata, unsigned int ifeIQDataLen, unsigned int ifeNFFT);
	int myCompareInt(const void *apple, const void *banna);		// for qsort in C
	int myCompare(const void *apple, const void *banna);		// for qsort in C

	float coarseSNRestimation(signed int snrFFTLength);	
	float coarseBandWidthEstimation(float* BWEdata,unsigned int BWE_FFTLength);
	float coarseSymRateEstimation(float* cseSpectrum, int cseNFFT);  // only for digital signal

public:
	void filter(float* inB, float* inX,unsigned int LenB, unsigned int LenX, float* filterout, unsigned int realFlag);
	void four1(float* data, unsigned long nn, int isign);
	void movingAverageFilter(float* MovInSeq, signed int MovInSeqLen, signed int MovWidth, float* MovOutSeq);	
	void fit(float x[], float y[], unsigned int ndata, float *a,float *b, float *siga, float*sigb, float *lmsErr);
	void fftshift(float* infftshift, float* outfftshift, unsigned int lenfftsh);	
	void siftDown( float *a, int *idx, int start, int end); // Call by Reference
	void siftDownInt( int *a, int *idx, int start, int end); // Call by Reference
	void heapsort( float *a, int *idx, int count); // descend order
	void heapsortInt( int *a, int *idx, int count); // descend order
	void histogram(float* inHist, unsigned int inHistLength, unsigned int histNumBin, float maxy, float miny);
	void FMvsFSKClassification(float* fcSpectrum, signed int fcNFFT); // 15.05.27


	// for demodulation
	void PSKdemodulation(float* inPSKdemod);
	unsigned int FMdemodulation(float* inFMdemod, unsigned int sigLen);
	void AMdemodulation(float* inAMdemod);
	void AGC(float* inAGC,int updatePeriod);
	void coarseFreqOffsetCompensation(float* inIQ,float Fs,float Rs);
	int samplingRateConvertor(float* inSRC,float Fs, float Rs);
	int fineCarrierOffsetNSymbolTimingRecovery(float* inIQ,int L);
	void carrierPhaseRecovery(float* inPSKdemod,int L,int updatePeriod);
	void symbolDecision(float* inIQ,int L,int modOrder);
	

	int medianI(int* mediInput,int startIdxMed2, int endIdxMed2);
	int complexMultiply(float* inComp, float* inRef, unsigned int compSigLen, signed int compNum);
	int nonlinearTransform(float* inNT,int ntLen, int nTn, int nTq);

	float dataConditining(float* dcInput,int dcInputLen);
	float medianF(float* inMedian,int startIdxMed, int endIdxMed );
	float myCabs(float re, float im);
	float sinc(float inX);
	float cm(float* inCm, int inCMLength, unsigned int nOrder, unsigned int qConj);
	// for demodulation
	float mod(float* inX, float inY);
	float my_round(float* inRound);

	// wav file
	void write_little_endian(unsigned int word, int numBytes, FILE *wavFile);
	void write_wav(char *filename, unsigned long num_samples, float *data, int s_rate);
	void write_float_little_endian(float word, int numBytes, FILE *wavFile);
};


AMRalgorithm::AMRalgorithm()
{
	int i=0;
	
	//memset(amrObj,0,sizeof(amrObj));

	for(i=5; i<=12; i++)	AMR_decimationTable[i-5] = (int)pow(2.0,i);
	for(i=0; i<5; i++)		AMR_SNRTable[i] = 10+5*i;
	AMR_SNRTable[5] = 9999;
	for(i=0; i<AMR_MAX_NUM_SEGMENT; i++){
		IQbeginIdx[i] = IQendIdx[i] = IQnum[i] =0 ;
	}

    AMR_ADCsamplingFrequency = 140000000;	 //140MHz
    AMR_Decimation 		     = 256;
    AMR_requireIQsamples = AMR_decimationIdx	             = 0;
    AMR_samPeriod = AMR_freqInterval = AMR_samplingFrequency = 0.0;
    AMR_LEN				     = 64;
    AMR_NFFT				 = 2;  
	AMR_numIQFPGA            = 0;
	AMR_symRateEstSegLen     = 0;
	AMR_isFSK                = 0;
	AMR_noisePower			 = 0;
	AMR_SNR_dB				 = 0;
	AMR_noiseFloor_dB        = 0;

}

#endif
