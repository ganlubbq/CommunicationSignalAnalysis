#ifndef _assignVariable
#define _assignVariable

// Last modified by AWH at 20-Aug-2015

#if VISUALSTUDIO
	FILE *fa;
	int AMR_numIQFPGA = 0;
#endif
//assign global variables

#ifdef activateCDW
// for input IQ data
	float AMR_inputData[AMR_MAXArray]; // complex IQ data, as start from index of 1(not 0),
				  						// real and imag is stored in odd and even element, respectively
#endif

// struct AMR_struct { ... } 구조체를 amrObj 이라는 데이터형으로 정의한다.
// http://www.tipssoft.com/bulletin/board.php?bo_table=FAQ&wr_id=968
//

typedef struct AMR_struct{		//150416
	unsigned int IQbeginIdx;
	unsigned int IQendIdx;		// the last index in the IQ array
	unsigned int IQnum;			// the number of IQ samples
	unsigned int modNum;		// modulation type
	unsigned int modOrder;		// modulation order

	float coarseSymbolRate;
	float coarseBandWidth;
	float coarseFreqOffset;
	float coarseToneSpacing;     // float -> int , 14072
} amrData;

int			 AMR_decimationTable[8] = {32,64,128,256,512,1024,2048,4096};
//unsigned int AMR_SNRTable[6]       = {10,15,20,25,30,9999};
short        AMR_enable_awgn = 1;

// global variable
int   AMR_ADCsamplingFrequency = 140000000;	 //140MHz
int   AMR_Decimation 		   = 1024;
unsigned int AMR_decimationIdx = 0;
unsigned int AMR_requireIQsamples = 0;
unsigned int AMR_maxSpecMagIdx  = 0;
float AMR_samplingFrequency		= 0.0f;
float AMR_freqInterval  	    = 0.0f;
float AMR_samPeriod 	        = 0.0f;
// SNR control
float AMR_SNR_dB				= 0.0f;
float AMR_noisePower            = 0.0f;
float AMR_noiseFloor_dB         = 0.0f;

unsigned int   AMR_LEN  = 64;
unsigned int   AMR_samplesPerSym = 0;
unsigned int   AMR_numSym        = 0;
signed int    AMR_NFFT =  2;

// parameter estimation
float AMR_bwThreshold_dB     = 0.0;  // threshold of normalized spectrum  when 0 (min) to 1 (max)
float AMR_bwThreshold        = 0.0;
float AMR_coarseBandWidth 	 = 0.0;
float AMR_coarseFreqOffset	 = 0.0;
float AMR_coareseSNR_dB      = 0.0;
float AMR_coarseSymbolRate   = 0.0;
float AMR_coarseToneSpacing  = 0.0;
float AMR_symratePeak        = 0.0;   // 15.05.27
float AMR_sigFloorRatio      = 0.0;	  // 15.06.16 FM vs FSK classification


int   AMR_coarseSamplesPerSymbol = 0;
int   AMR_coarseNumberOfsymbols  = 0;
signed int   AMR_symratePeakIdx = 0; // 15.05.27
signed int AMR_symRateEstSegLen = 0; // 15.06.16
unsigned int   AMR_numSeg	     = 0;   				// the number of segments
unsigned int   AMR_coareseSNRIdx = 0;

//for AMR_DEBUGging
int   AMR_DEBUG_numseg	 = 0;
float AMR_DEBUG_floatTmp = 0.0;


//for AM,CW vs others by using whether existance of the carrier component
unsigned int AMR_existsCarrierFlag = 0;
//for analog vs digital signal identification bits
unsigned int AMR_digitalFlag = 0;
//for linear vs non-linear signal identification bits
unsigned int AMR_linearFlag  = 0;

unsigned int AMR_minBinInterval = 0;
int AMR_modNum         = 0;
int	AMR_modOrder       = 0;
int AMR_isFSK          = 0;

float AMR_inputData[AMR_MAXArray];
float AMR_filteredInputData[AMR_MAXArray];

// for histogram
unsigned int	AMR_binCount[AMR_MAX_NUMBINS+2];
float			AMR_binCenter[AMR_MAX_NUMBINS+2];
float			AMR_bins[AMR_MAX_NUMBINS+2];	// 구역 분할 범위, numAMR_bins의 중심값

unsigned int	AMR_numPeaks = 0;							// for mode operation
// for welchPSD function 
float AMR_IQSegment_PSD[(AMR_MAX_welchPSDsize<<1) + 1]; // welchPSD inner array
float AMR_hanningWindow_PSD[AMR_MAX_welchPSDsize];  // welchPSD inner array
float AMR_Sxx_PSD[AMR_MAX_welchPSDsize];			// welchPSD inner array
float AMR_powerSpectrum[AMR_MAX_welchPSDsize];		// welchPSD output
float AMR_powerSpectrum_dB[AMR_MAX_welchPSDsize];   // welchPSD output

float AMR_TSE_PSDbuffer[AMR_MAX_welchPSDsize];
int   AMR_TSE_indexBuffer[AMR_MAX_welchPSDsize];

float AMR_peakValues[AMR_MAX_welchPSDsize >>1];			// for tone spacing estimation, PeakDetection
signed int AMR_peakIndices[AMR_MAX_welchPSDsize >>1];	// for tone spacing estimation, PeakDetection
signed int AMR_freqVector[AMR_MAX_welchPSDsize];
unsigned int AMR_IndexBuffer_PSD[AMR_MAX_welchPSDsize];

signed int   AMR_BWSampleIdx[2];
signed int   AMR_BWSampleIdx_TMP[2];

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

//Coarse symbol rate estimation(CSRE)
float AMR_ck[AMR_symbolRate_SamplesPerSegment];				 // Coarse symbol rate estimation
float AMR_iPhase[AMR_symbolRate_SamplesPerSegment];			 // Coarse symbol rate estimation
float AMR_iFrequency[AMR_symbolRate_SamplesPerSegment];	     // instantaneous frequency
float AMR_iAmplitude[AMR_symbolRate_SamplesPerSegment];		 // linear vs non-linear classification
float AMR_spectrumBuffer[AMR_symbolRate_SamplesPerSegment];
float AMR_filteredSamples[AMR_symbolRate_SamplesPerSegment];		 // Coarse symbol rate estimation
float AMR_spectrumMagnitude[AMR_symbolRate_SamplesPerSegment]; // Coarse symbol rate estimation

float AMR_halfFreqVector[AMR_symbolRate_SamplesPerSegment>>1];  // for symbol rate estimation
float AMR_IQsegmentSym[(AMR_symbolRate_SamplesPerSegment<<1) + 1]; // for symbol rate estimation

// features
float AMR_sigma_a[AMR_MAXNumOfFeature];
float AMR_absCC20[AMR_MAXNumOfFeature];
float AMR_absCC40[AMR_MAXNumOfFeature];
float AMR_lineLMSerror[AMR_MAXNumOfFeature];
float AMR_mu42f[AMR_MAXNumOfFeature];

float AMR_dataConditioningTmp[AMR_MAXNumOfFeature];

// for cyclic cumulant
float AMR_IQbufferCC[(AMR_CC_samplesPerSegment<<1) + 1];   // nonlinear transform 
float AMR_IQsegmentCC[(AMR_CC_samplesPerSegment<<1) + 1];  // for cyclic cumulant
float AMR_IQDataTmp2[(AMR_CC_samplesPerSegment<<1) + 1];   // nonlinear transform 
float AMR_IQDataTmp1[(AMR_CC_samplesPerSegment<<1) + 1];   // nonlinear transform 
float AMR_IQDataResQ[(AMR_CC_samplesPerSegment<<1) + 1];   // nonlinear transform 
float AMR_spectrumMagnitudeCC[AMR_CC_samplesPerSegment]; // for cyclic cumulant

//float AMR_oneIQData[3];		// for cumulant
//float AMR_oneIQDataSec[3];		// for cumulant

// for IQ segmentation
unsigned int beginTOAidx[AMR_MAX_NUM_SEGMENT];
unsigned int endTOAidx[AMR_MAX_NUM_SEGMENT];
unsigned int AMR_numTOA = 0;
amrData amrObj[AMR_MAX_NUM_SEGMENT];

#if VISUALSTUDIO
	//// Demodulation
	/*float AMR_FMdemodiPhase[IQ_NUM_MAX];
	float AMR_FMdemodiFreq[IQ_NUM_MAX];
	float AMR_FMdemodCk[IQ_NUM_MAX];
	float AMR_demodIQ[(IQ_NUM_MAX<<1) +1];*/
	////float AMR_inputPhase[AMR_MAXArray];
	unsigned int AMR_wavBuffer[IQ_NUM_MAX];
	float AMR_FMdemodiPhase[IQ_NUM_MAX];
	float AMR_FMdemodiFreq[IQ_NUM_MAX];
	float AMR_FMdemodCk[IQ_NUM_MAX];
	float AMR_demodIQ[(IQ_NUM_MAX<<1) +1];
	float AMR_FMdemodiFreqBuff[IQ_NUM_MAX];
#endif

//float AMR_demodSig[AMR_MAXArray];
//float AMR_filteredInstantFreq[AMR_MAXArray];
float AMR_interpFiltOut[2];
float AMR_fineCompOut[2];
float AMR_pDelayBuffer3[2];
float AMR_pDelayBuffer2[2];
float AMR_pDelayBuffer1[2]; 
float AMR_pTEDDelay1[2]; // interpolant buffer
float AMR_pTEDDelay2[2];
float AMR_timingLoopOut[AMR_MAXArray];


#if VISUALSTUDIO
	unsigned int IQ_Buffer[IQ_NUM_MAX];
	unsigned int IQ_BufferConv[IQ_NUM_MAX];
	short I_BufferConv[IQ_NUM_MAX];
	short Q_BufferConv[IQ_NUM_MAX];
	float I_BufferConvF[IQ_NUM_MAX];
	float Q_BufferConvF[IQ_NUM_MAX];
#endif

#endif
