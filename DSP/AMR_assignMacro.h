#ifndef _assignMacro
#define _assignMacro

// Last modified by AWH at 1-Sep-2015

// 256 Decimation, 1ksps, 200symbols condition --> unused
#define AMR_CW_minimumIQ	    10000	   //15.08.20109375

#if VISUALSTUDIO
	#define IQ_NUM_MAX			8000000	   // only for VISUALSTUDIO
	#define AMR_WAV_S_RATE (44100)
	#define AMR_WAV_BUF_SIZE (AMR_WAV_S_RATE*100) // 100 second buffer
#endif

// AMR_PI to machine precision, defined in math.h
#define AMR_PI					3.141592653589793  //  the ratio of the circumference of a circle to its diameter
#define AMR_TWOPI				6.283185307179586

/////////////////////////// Important factor in AMR

#define AMR_MAXArray 					400001  // 
												// 262145  //  2^18+1, complex array,   --> (real+imag) + reserved 1 bits (0 index)
#define AMR_MAXNumOfFeature				20      // 2000 -> 30   150407
#define AMR_MAX_NUM_SEGMENT				5000    // the maximum number of segment   1500602
#define AMR_MAX_NFFT					8192    // 15.08.20
#define AMR_MIN_NFFT					1024    // AMR_CW_minimumIQ / 8 refer to setFreqVector()
#define AMR_MAX_welchPSDsize			AMR_MAX_NFFT   // array size for welch PSD
#define AMR_BW_SamplesPerSegment		8192	//should be power of 2, for PSD and coarse bandwidth estimation
#define AMR_samplesPerSegment			16384
#define AMR_oneHopMinIQ					10000    // The minimum IQ samples in a hop
#define AMR_MINIMUM_HIST_COUNT          13     // the counts of histogram in bandwidth estimation


#define AMR_instantAmpObservationInterval		4000  // Adjust for 8-FSK under SNR=15dB at 10-Sep-2015	// sigma_a
#define AMR_instantPhaseObservationInterval		3000
#define AMR_symbolRate_SamplesPerSegment		16384
#define AMR_toneSpacing_FFTsize					8192
#define AMR_CC_samplesPerSegment				32768	//which one is better? 65536 and 32768

// Spectrum trace ( reserved )
#define AMR_MaxHold 			771
#define AMR_ClearWhite 			772

// LPF
#define AMR_LPFORDER	31				// LPF order
#define AMR_LPFLEN		(AMR_LPFORDER + 1)

#define AMR_LINEAR_MOD_MAX_BW   61500 // 30kHz --> 61.5kHz at 150331
#define AMR_minimumSymbolRate   900 // 900Hz
#define AMR_maximumSymbolRate  	41000 // 20kHz --> 41kHz at 150331
#define AMR_spsLowerBound	    13 	// depends on cyclic cumulant
#define AMR_spsUpperBound	    32

#define AMR_CARRIER_LINE_SPECTRUM_THRESHOLD	    13.0f  //dB, 01-Oct-2015

// Value for classification between analog and digital modulation signal
#define AMR_SIGMA_SQUARE_A_THRESHOLD			0.0676f // 0.26^2 for 8-FSK under SNR=15dB at 17-Sep-2015

// Value for classification between AM,CW and others
#define AMR_SIGNALFLOOR_THRESHOLD			10

// Value for classification ASK vs PSK, QAM
#define AMR_LMS_ERROR_THRESHOLD		1.5f

// Value for BPSK vs M-PSK, QAM
#define AMR_ABSCC20_THRESHOLD		0.40f   // 045f --> 0.40 at 21-Sep-2015
#define AMR_ABSCC40_THRESHOLD		0.30f	// 0.35f --> 0.30f at 150611
#define AMR_QAMvsPSK_THRESHOLD  	0.1369f

// Value for kurtosis of instantaneous frequency
#define AMR_FREQ_KURTOSIS_THRESHOLD 2.65f  // 2.8 --> 2.65 refer to AMRtestbed_script_FMvsFSK

// knn
#define	AMR_refDataOffset		45		// the number of reference data corresponding to each modulation type
#define	AMR_K_NEIGHBOR			3		// the number of k-neighbor
#define AMR_NUM_LINEAR_MOD		3		// coressponding to QPSK, 8PSK, 32QAM

// histogram
#define	AMR_MAX_NUMBINS				35		// x of matlab (ÃÊ±â °ª)

//#define AMR_IS_LESS(v1, v2)  (v1 < v2)								//heap sort
#define AMR_SWAP_DOUBLE(r,s) do{float t=r; r=s; s=t; } while(0)		//heap sort
#define AMR_SWAP_FLOAT(r,s)  do{float t=r; r=s; s=t; } while(0)		//heap sort
#define AMR_SWAP_INT(r,s)    do{int t=r; r=s; s=t; }   while(0)	    //heap sort
#define AMR_VSIZE(d1vector) (sizeof(d1vector)/sizeof(d1vector[0]))	//heap sort

// Signal-to-noise ratio (SNR)
#define AMR_J		2045L
#define AMR_MM		1048576                /*(0x01L<<SHIFT)= 1048576*/
#define AMR_SEED	78923

// straight line fitting
//static float sqrarg;
//#define AMR_SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg) // slow method in embedded system


#endif
