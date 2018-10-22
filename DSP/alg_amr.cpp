/*
// \ DESCRIPTION
//  > Automatic modulation recognition (AMR) algorithm 
//    for modulation signals such as AM, FM, FSK, ASK, PSK, and QAM
//
// \ Notice
//  > Before DSP test, you should modify this file extension to (*.c) 
//
// \ See also 
//   > revision.txt
// 
// \ Reference
//   [1] John G. Proakis, Dimitris G. Manolakis,
//       "Digital Signal Processing", Pearson Prentice Hall,
//       p.267, p.475, 2007
//   [2] time measure, http://e2e.ti.com/support/dsp/tms320c6000_high_performance_dsps/f/112/p/37848/136807
//
// Last modified by AWH at 07-Sep-2015
*/

//#define activateCDW     // Enable for DSP
#define VISUALSTUDIO		1  // Disable for DSP
#define AMR_DEBUG			1
#define AMR_PROFILE			0
#define AMR_Cplus			1  // Disable for DSP

#include <stdlib.h>  // for qsort
#include <stdio.h>   // UART, smooth
#include <math.h>	 // sin, sqrt, pow, FT
#include <string.h>
#include "amr_assignMacro.h" // It should be located in front of program

#if VISUALSTUDIO
	#include <assert.h>
#endif

#if AMR_Cplus
	#include "AMR_class.h"			 // AMR algorithm object
	//#include "testData.h"
	//#include "FM_IQ.h"
	// http://breadlab.net/36
	//#include <time.h>
	
    #include <Windows.h>
	#include <mmsystem.h>			// play *.wav
	#include <conio.h>
	//#pragma comment(lib,"winmm.lib") 
	
#else
	#include "AMR_assignVariable.h"

#endif


#ifdef activateCDW
	// function profiling
	#include <c6x.h>

	// Include file for this module
	// General type
	#include <tistdtypes.h>
	#include <usertype.h>
	#include "../ns/alg_proc_stc.h"

	#include "app_debug.h"
	#include "app_main.h"
	#include "app_memmap.h"
#endif

#ifdef activateCDW
	// from FPGA
	Uint32 uiAMRSYSCHBw;
	Uint32 uiAMRSamplingFreq;
	Uint32 uiAMRSYSFFTW;
	Uint32 uiAMRSYSOVPSD;
	Uint32 uiAMRSYSFFTL;
	Uint32 uiAMRMVAVGFIL;
	Uint32 uiAMRTHSBW;
	Uint32 uiAMRTHPKD;
	Uint32 uiAMRMinPeakDis;
	Uint32 uiSNRfromGUI; //  uiAMRKneighbor;
	Uint32 uiAMRCDWindex;//CDW index
	Uint32 uiAMRDecimation;
	//float uiAMRDecimation;    //141205
	Uint32 uiAMRrsvd2;

	extern float I_BufferConvF[IQ_NUM_MAX*50];
	extern float Q_BufferConvF[IQ_NUM_MAX*50];

#endif

//#include "refData_150121.h"
//#include "refData_150131.h" 
//#include "refData_150203.h"		// 2ASK 제거
//#include "refData_150331.h"
//#include "refData_150407.h"
#include "AMR_function.h"
//#include "demodulation.h"

#ifdef activateCDW
	#include "AMR_parsingIQdata.h" 	// 150416, read data from *.diq file 
#else
	#include "readIQFile.h"
#endif

#ifdef activateCDW
//CDW_SEND
int LOCAL_send_CDW_AMR(int index)
{
	float tmpf = 0.0f;

	if (AMR_modOrder == 0)
	{	
		tmpf = 0.0f;
	}
	else{
		tmpf = log2(AMR_modOrder);
	}
	// fill the structure for MDD
	//unsigned int
	atMDD[0].iModeType			 = AMR_modNum;
	atMDD[0].iModOrder			 = AMR_modOrder;
	atMDD[0].iSampleNum  		 = AMR_LEN;
	atMDD[0].iChannelCoding 	 = (AMR_requireIQsamples << 1) + AMR_decimationIdx;
	atMDD[0].iPulseShape 		 = 0x0;

	// float
	atMDD[0].fSymbolRate  		 = (float) AMR_coarseSymbolRate;
	atMDD[0].fBitRate			 = AMR_coarseSymbolRate * tmpf;
	atMDD[0].fBandWidth		 	 = (float) AMR_coarseBandWidth;
	atMDD[0].fFreqOffset  		 = (float) AMR_coarseFreqOffset;
	atMDD[0].fRollOffFactor	 	 = 0.0f;
	atMDD[0].fBandWidth_time 	 = 0.0f;
	atMDD[0].fMaxFreqDev 		 = 0.0f;
	atMDD[0].fFreqSep		     = (float) AMR_coarseToneSpacing;     // float -> int , 140722

	/* Printout debug log */
	trace(TRACE_MSG_LOG, "=>LOCAL_send_FDD()");
	/******************/

 	printf("LOCAL_send_CDW_AMR Start....\n");

	//Copy index of ATCDW
	atCDW[0].uiIndex = FPGA_uiCDWindex;
	atCDW[0].uiComintState = 2;
	atCDW[0].rev1 = 0;	 // add 140722
	atCDW[0].rev2 = 0;
	atCDW[0].rev3 = 0;

//	memcpy( &atCDW[0].atSMDD, &atMDD[0], sizeof(SSC_MDD_TBL));
	atCDW[0].atSMDD.fBandWidth     = atMDD[0].fBandWidth;
	atCDW[0].atSMDD.fBandWidth_time= atMDD[0].fBandWidth_time;
	atCDW[0].atSMDD.fBitRate       = atMDD[0].fBitRate;
	atCDW[0].atSMDD.fFreqOffset    = atMDD[0].fFreqOffset;
	atCDW[0].atSMDD.fMaxFreqDev    = atMDD[0].fMaxFreqDev;
	atCDW[0].atSMDD.fRollOffFactor = atMDD[0].fRollOffFactor;
	atCDW[0].atSMDD.fSymbolRate    = atMDD[0].fSymbolRate;

	atCDW[0].atSMDD.iChannelCoding = atMDD[0].iChannelCoding;
	atCDW[0].atSMDD.iModeType	   = atMDD[0].iModeType;
	atCDW[0].atSMDD.iModOrder	   = atMDD[0].iModOrder;
	atCDW[0].atSMDD.iPulseShape    = atMDD[0].iPulseShape;
	atCDW[0].atSMDD.iSampleNum     = atMDD[0].iSampleNum;
	atCDW[0].atSMDD.fFreqSep       = atMDD[0].fFreqSep;

	FPGA_CDW_SIZE = (unsigned short)1;

	memcpy( (void *)FPGA_PCI_DATA_BASE, atCDW, sizeof(SSC_CDW_TBL));

	trace(TRACE_MSG_LOG, "=>LOCAL_send_FDD():Sending PCI data");
	printf("LOCAL_send_CDW End....\n");
	FPGA_PCI_SEND_ON = 0;
	TSK_sleep(1);
	FPGA_PCI_SEND_ON = 1;
	TSK_sleep(1);
	FPGA_PCI_SEND_ON = 0;
	TSK_sleep(100);
	return 0;
}

void AMRparameterSet(void)
{
	uiAMRSYSCHBw = FPGA_uiAMRSYSCHBw;
	uiAMRSamplingFreq = FPGA_uiAMRSamplingFreq;
	uiAMRSYSFFTW = FPGA_uiAMRSYSFFTW;
	uiAMRSYSOVPSD = FPGA_uiAMRSYSOVPSD;
	uiAMRSYSFFTL = FPGA_uiAMRSYSFFTL;
	uiAMRMVAVGFIL = FPGA_uiAMRMVAVGFIL;
	uiAMRTHSBW = FPGA_uiAMRMVAVGFIL;
	uiAMRTHPKD = FPGA_uiAMRTHPKD;
	uiAMRMinPeakDis = FPGA_uiAMRMinPeakDis;
	//uiAMRKneighbor = FPGA_uiAMRKneighbor;
	uiSNRfromGUI  = FPGA_uiAMRKneighbor;
	uiAMRCDWindex = FPGA_uiCDWindex;//CDW index
	uiAMRDecimation = AMR_decimationTable[FPGA_uiDecimation]; 		// 141205
	uiAMRrsvd2 = FPGA_uiAMRrsvd2;
}

#endif


#if AMR_Cplus	

	int main(){
		unsigned int end   = 0, begin = 0;
		signed int tmpi=0;	
		signed int idx=0, i=0, toaLoop;
		signed int offset = 0;
		// time measure
		time_t startTime =0, endTime = 0;
		time_t totStartTime, totEndTime = 0;
		float gap, tmpf;

		float min = 0; // for demodulation
		float max = 0;
		
		//AMRalgorithm amr; // overflow occurs in stack
		AMRalgorithm *amr;  // use dynamic memory in heap 
		amr = new AMRalgorithm();

	#if AMR_DEBUG
			//amr->AMR_LEN = 16384;
			//amr->initializeVariable();
			//amr->PSKdemodulation(AMR_inputData);
			//amr->FMdemodulation(AMR_inputData);
			//amr->AMdemodulation(AMR_inputData);
	#endif

	#if 0 //AMR_PROFILE			
			startTime = clock();
	#endif

		amr->AMR_numTOA = 0;			
		amr->AMR_numTOA = amr->LOCAL_load_IQ("20150608_005215(FM_CW_1024D).diq");
		//amr->AMR_numTOA = amr->LOCAL_load_IQ("20150615_022431(BPSK_2ksps_1024D).diq"); amr->AMR_Decimation = 1024;
		//amr->AMR_numTOA = amr->LOCAL_load_IQ("20150509_065751(FM_hopping_audio_1024D).diq");
		amr->AMR_Decimation = 1024;
	#if  0 //AMR_PROFILE			
			endTime = clock();
			gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("\r\n Exe. time for amr->LOCAL_load_IQ() : %.7f sec\r\n",gap);
			//amr->AMR_numTOA = amr->LOCAL_load_IQ("20150430_145331(4FSK_0.35rc_10ksps_256deci).diq");
	#endif		
		//amr->AMR_numTOA = amr->LOCAL_load_IQ("20150508_063040(FM_hopping_1024D_3000threshold).diq");

		for(toaLoop=0; toaLoop<amr->AMR_numTOA; toaLoop++){

	#if AMR_PROFILE						
				totStartTime = clock();
	#endif	
			amr->initializeVariable();


	#if 0 //AMR_PROFILE		
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->initializeVariable() : %.7f sec\r\n",gap);
	#endif			
			amr->AMR_LEN = amr->IQnum[toaLoop];
			/////////////////////////////////////////////////////////////////////////
			// Determine whether the number of IQ samples is sufficient or not
			/////////////////////////////////////////////////////////////////////////

			// is IQ samples valid?
			tmpi = (AMR_MAXArray-1) >> 1;
			if (amr->AMR_LEN > tmpi){
				amr->AMR_LEN = tmpi;
				#if AMR_DEBUG
					printf("\r\n Too many IQ samples, the num. of IQ samples is restricted to %d!!! \r\n",amr->AMR_LEN);
				#endif
			}
			else if (amr->AMR_LEN < AMR_CW_minimumIQ)
			{
				#if AMR_DEBUG
					printf("\r\n Insufficient IQ samples %d!!! \r\n",amr->AMR_LEN);
				#endif

				#ifdef activateCDW
					LOCAL_send_CDW_AMR(FPGA_uiCDWindex);
				#endif
				return 0;
			}
			else
			{
				#if AMR_DEBUG
				printf("\r\n Sufficient IQ samples!!! \r\n");
				#endif
			}

			#if AMR_DEBUG
				printf("\r\n Run AMR #[%d] with %d samples \r\n",toaLoop, amr->AMR_LEN );
			#endif

			/////////////////////////////////////////////////////////////////////////
			// Load IQ samples
			/////////////////////////////////////////////////////////////////////////	
			
			begin = amr->IQbeginIdx[toaLoop]; 
			end   = begin + amr->AMR_LEN; // amrObj[toaLoop].IQendIdx;

			i=0;
			for(idx=begin; idx<=end; idx++){
				amr->AMR_inputData[2*i+1] = amr->I_BufferConvF[idx];
				amr->AMR_inputData[2*i+2] = amr->Q_BufferConvF[idx]; 
				i++;
			}

			amr->AMR_SNR_dB = 15; // dB

	#if 0 //AMR_PROFILE			
				startTime = clock();
	#endif

			// Set frequency bins
			amr->AMR_NFFT = amr->setFrequencyVector(amr->AMR_LEN);

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->setFrequencyVector() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif

			// Force signal power to unit
			amr->unitPower(amr->AMR_inputData, amr->AMR_LEN);

	#if 0 //AMR_PROFILE		
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->unitPower() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif

			// Calculate welch's PSD for bandwidth estimation
			amr->welchPSD(amr->AMR_inputData, amr->AMR_LEN, amr->AMR_NFFT, 1, amr->AMR_powerSpectrum);

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->welchPSD() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif

			// Bandwidth estimation
			amr->AMR_coarseBandWidth = amr->coarseBandWidthEstimation(amr->AMR_powerSpectrum,amr->AMR_NFFT);

	#if 0 //AMR_PROFILE						
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->coarseBandWidthEstimation() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif

			// Add white Gaussian noise to signal	
			amr->awgn(amr->AMR_inputData, amr->AMR_SNR_dB, amr->AMR_noisePower);

	#if 0 //AMR_PROFILE		
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->awgn() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif
				
//	#if AMR_PARAM_ON
//			// SNR estimation
//			amr->AMR_coareseSNR_dB = amr->coarseSNRestimation(amr->AMR_NFFT);
//	#endif

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->coarseSNRestimation() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif		
			// carrier frequency offset estimation
			offset = amr->coarsefreqOffsetEstimation(amr->AMR_LEN, amr->AMR_NFFT); // signed argument

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->coarsefreqOffsetEstimation() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif	
			// low-pass filtering
			amr->lowPassFiltering(amr->AMR_inputData);

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->lowPassFiltering() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif	
			amr->linearVsNonLinearClassification(amr->AMR_filteredInputData, amr->AMR_LEN, AMR_instantAmpObservationInterval);

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->linearVsNonLinearClassification() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif	
			amr->spectrumFeatureExt(amr->AMR_filteredInputData, amr->AMR_LEN, amr->AMR_NFFT);

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->spectrumFeatureExt() : %.7f sec\r\n",gap);
			
				startTime = clock();
	#endif	
			amr->AMR_coarseSymbolRate = amr->coarseSymRateEstimation(amr->AMR_spectrumMagnitude, amr->AMR_NFFT);

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->coarseSymRateEstimation() : %.7f sec\r\n",gap);
			
				startTime = clock();
	#endif	
			amr->FMvsFSKClassification(amr->AMR_spectrumMagnitude, amr->AMR_symRateEstSegLen);

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->FMvsFSKClassification() : %.7f sec\r\n",gap);
	#endif	

			amr->reSampling(amr->AMR_LEN, amr->AMR_coarseSymbolRate);

	#if 0 //AMR_PROFILE			
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->reSampling() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif	
			// require to edit the program 
			amr->coarseToneSpacing(amr->AMR_filteredInputData, amr->AMR_LEN, amr->AMR_NFFT, AMR_toneSpacing_FFTsize,offset);  	// only for FSK digital signal

	#if 0 //AMR_PROFILE				
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->coarseToneSpacing() : %.7f sec\r\n",gap);

				startTime = clock();
	#endif		
			amr->linearModClassification(amr->AMR_filteredInputData, amr->AMR_LEN, AMR_CC_samplesPerSegment);   // only for linear digital signal

	#if 0 //AMR_PROFILE	
				endTime = clock();
				gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->linearModClassification() : %.7f sec \r\n",gap);
	#endif

	#if AMR_PROFILE	
				totEndTime = clock();
				gap = (float)(totEndTime-totStartTime)/(CLOCKS_PER_SEC);
				printf("\r\n Exe. time for amr->linearModClassification() : %.7f sec \r\n",gap);
	#endif
			amr->displayResults();

		//	///// add code to store the result to amrObj
		//	amr->amrObj[toaLoop].coarseBandWidth = amr->AMR_coarseBandWidth;
		//	amr->amrObj[toaLoop].coarseFreqOffset = amr->AMR_coarseFreqOffset;
		//	amr->amrObj[toaLoop].coarseSymbolRate = amr->AMR_coarseSymbolRate;
		//	amr->amrObj[toaLoop].coarseToneSpacing = amr->AMR_coarseToneSpacing;

		//	amr->amrObj[toaLoop].modNum = amr->AMR_modNum;
		//	amr->amrObj[toaLoop].modOrder = amr->AMR_modOrder;

			//
			// FM Demodulation process
			//

			// acquire IQ samples into AMR_demodIQ
			if(amr->AMR_numTOA > 10){	//hopping 
				i=0;
				for(toaLoop=0; toaLoop<amr->AMR_numTOA; toaLoop++){
					begin = amr->IQbeginIdx[toaLoop]; 
					end   = begin + amr->AMR_LEN; // amrObj[toaLoop].IQendIdx;
					for(idx=begin; idx<=end; idx++){
						amr->AMR_demodIQ[2*i+1] = amr->I_BufferConvF[idx];
						amr->AMR_demodIQ[2*i+2] = amr->Q_BufferConvF[idx]; 
						i++;
					}
				}

				amr->AMR_LEN = i;
			}
			else
			{						//CW
				amr->AMR_LEN = amr->AMR_numIQFPGA;
				if (amr->AMR_numIQFPGA > IQ_NUM_MAX) amr->AMR_numIQFPGA = IQ_NUM_MAX;
				for(i=0; i<amr->AMR_numIQFPGA; i++){
					amr->AMR_demodIQ[2*i+1] = amr->I_BufferConvF[i];
					amr->AMR_demodIQ[2*i+2] = amr->Q_BufferConvF[i]; 
				}
			}

			// Execute FM demodulation
			amr->AMR_LEN = amr->FMdemodulation(amr->AMR_demodIQ,amr->AMR_LEN);
			
	
			// find min & max
			min = amr->AMR_FMdemodiFreqBuff[0];
			for(i=1; i<amr->AMR_LEN; i++){
				if ( min > amr->AMR_FMdemodiFreqBuff[i]){
					min = amr->AMR_FMdemodiFreqBuff[i];					
				}			
			}
			max =  min;
			for(i=1; i<amr->AMR_LEN; i++){
				if ( max < amr->AMR_FMdemodiFreqBuff[i]){
					max = amr->AMR_FMdemodiFreqBuff[i];					
				}			
			}
			for(i=0; i<amr->AMR_LEN; i++){
				tmpf = (amr->AMR_FMdemodiFreqBuff[i] - min ) / (max-min); // scale conversion --> 0 to 1						
				tmpf *= 65536.0f;
				amr->AMR_FMdemodiFreqBuff[i] = tmpf;
			}

		//amr->write_wav("demodSig.wav", amr->AMR_LEN, amr->AMR_FMdemodiFreqBuff, 132300);
		amr->write_wav("demodSig.wav", amr->AMR_LEN, amr->AMR_FMdemodiFreqBuff, 136720);
		
		// play *.wav, from http://kkikkodev.tistory.com/54
		//PlaySound(TEXT("demodSig.wav"),NULL,SND_FILENAME | SND_ASYNC | SND_LOOP | SND_NODEFAULT);
		//printf("아무 키나 입력하시면 소리 재생이 멈춥니다.\n");
		//while (!_kbhit());
		//PlaySound(NULL, 0, 0);

		//
		} // end -- toaLoop
		delete(amr);
		return 0;
	}
#else
	#if VISUALSTUDIO
		int main()
	#else
		int process_amr()
	#endif
	{
		unsigned int end   = 0, begin = 0;
		unsigned int numFPGAiq = 0,toaLoop=0, i=0;
		unsigned int compConst=0;	
		unsigned  int idx=0 ;
		signed int offset = 0;
		//////// for de-hopping signal - temp
		unsigned int  j=0, k=0, fmCnt=0, amCnt=0;
		float a=0.0f, b=0.0f, siga=0.0f, sigb=0.0f, chi2=0.0f; // for linear curve fitting
		float iFreq = 0.0f;
		//for instantaneous phase
		float PI		  = AMR_PI;
		float TWOPI 	  = AMR_TWOPI;
		//////////////////		
		
#ifdef activateCDW
		#if AMR_PROFILE
			unsigned int startLow=0,startHigh=0, stopLow=0, stopHigh=0; // for profile
		#endif
		memset(amrObj,0,sizeof(amrObj));

		AMRparameterSet();
		/////////////////////////////////////////////////////////////////////////
		// Receive data from DDR memory
		/////////////////////////////////////////////////////////////////////////
		// Initialize IQ buffer
		#if AMR_DEBUG
				printf(" While receiving IQ samples, please wait... \r\n ");
		#endif

		numFPGAiq = LOCAL_load_IQ(); // implict type casting float => int

		#if AMR_DEBUG
				printf(" The num. of recevied IQ samples %d \r\n ",numFPGAiq);
		#endif
		AMR_numTOA = 0;
		AMR_numTOA = readTOA(numFPGAiq);
#if AMR_DEBUG
		printf(" # of TOA : %d \r\n",AMR_numTOA);
#endif		
	#else
		AMR_numTOA = 0;
		
		//AMR_numTOA = LOCAL_load_IQ("20150611_093704(4piDQPSK_18ksps_256D_nIQ_15000).diq"); AMR_Decimation = 256;
		//AMR_numTOA = LOCAL_load_IQ("20150610_113601(2FSK_h_0.5_20ksps_256D_nIQ_8000).diq"); AMR_Decimation = 256;
		//AMR_numTOA = LOCAL_load_IQ("[3-5] 4FSK(DMR)_info_PN_sps_273k512_rs_4.8k_sf_0.35_fd_1944_date_140105.diq"); AMR_Decimation = 512;
		//AMR_numTOA = LOCAL_load_IQ("[3-3-2] 2FSK_info_PN_sps_273k512_rs_1.2k_sf_0.35_fd_4.5k_date_140105.diq"); AMR_Decimation = 512; // POCSAG
		//AMR_numTOA = LOCAL_load_IQ("20150615_003925(DMR_1024D).diq"); AMR_Decimation = 1024;
		//AMR_numTOA = LOCAL_load_IQ("20150615_022431(BPSK_2ksps_1024D).diq"); AMR_Decimation = 1024;
		//AMR_numTOA = LOCAL_load_IQ("20150615_022711(BPSK_2ksps_1024D).diq"); AMR_Decimation = 1024;
		
		//AMR_numTOA = LOCAL_load_IQ("20150615_011657(2FSK_h_1_1.2ksps_1024D).diq"); AMR_Decimation = 1024;
		//AMR_numTOA = LOCAL_load_IQ("20150726_171045(FM_problem).diq"); AMR_Decimation = 1024;
		//AMR_numTOA = LOCAL_load_IQ("20150820_184606.diq"); AMR_Decimation = 512;
		
		
		//AMR_numTOA = LOCAL_load_IQ("20150608_005215(FM_CW_1024D).diq"); AMR_Decimation = 1024;
		//AMR_numTOA = LOCAL_load_IQ("20150616_015214.diq"); AMR_Decimation = 1024;
	
		
	    //AMR_numTOA = LOCAL_load_IQ("20150615_011117(2FSK_h_1_1.2ksps_1024D).diq"); AMR_Decimation = 1024;
		//AMR_numTOA = LOCAL_load_IQ("20150615_075905(MSK_1.2ksps_2048D).diq"); AMR_Decimation = 2048;
		
		//AMR_numTOA = LOCAL_load_IQ("20150527_103915.diq");
		//AMR_numTOA = LOCAL_load_IQ("20150623_032630(AM_Hopping_100hop_1024D).diq");
		//AMR_numTOA = LOCAL_load_IQ("20150509_065751(FM_hopping_audio_1024D).diq");
		
		//AMR_numTOA = LOCAL_load_IQ("20150907_054604_ASK_512D.diq");AMR_Decimation = 512;
		//AMR_numTOA = LOCAL_load_IQ("20150306_090656(4FSK_h_1_5ksps_rect_info_PN_1024deci).diq"); AMR_Decimation = 1024;
		//AMR_numTOA = LOCAL_load_IQ("20150306_090925(8FSK_h_1_5ksps_rect_info_PN_1024deci).diq"); AMR_Decimation = 1024;
		//AMR_numTOA = LOCAL_load_IQ("20150915_075212_AM_0dBm_90depth_1024D.diq"); AMR_Decimation = 1024;
		AMR_numTOA = LOCAL_load_IQ("20150918_064313_2FSK_5ksps_h_1.diq"); AMR_Decimation = 512;
		
	#endif
	if( AMR_numTOA == 0){
		initializeVariable();

		#ifdef activateCDW
			LOCAL_send_CDW_AMR(FPGA_uiCDWindex);
		#endif
	}else if((0 < AMR_numTOA) && (AMR_numTOA <= 5)){
		// CW signal	
		for(toaLoop=0; toaLoop<AMR_numTOA; toaLoop++){
		
				initializeVariable();

				AMR_LEN = amrObj[toaLoop].IQnum;
				/////////////////////////////////////////////////////////////////////////
				// Determine whether the number of IQ samples is sufficient or not
				/////////////////////////////////////////////////////////////////////////

				// is IQ samples valid?
				compConst = (AMR_MAXArray-1) >> 1;
				if (AMR_LEN > compConst){
					AMR_LEN = compConst;
					#if AMR_DEBUG
						printf(" Too many IQ samples, the num. of IQ samples is restricted to %d!!! \r\n",AMR_LEN);
					#endif
				}
				else if (AMR_LEN < AMR_CW_minimumIQ)
				{
					#if AMR_DEBUG
						printf(" Insufficient IQ samples %d!!! \r\n",AMR_LEN);
						// Send amr invalid message to GUI
						
					#endif

					#ifdef activateCDW
						LOCAL_send_CDW_AMR(FPGA_uiCDWindex);
					#endif
					return 0;
				}
				else
				{
					#if AMR_DEBUG
					printf(" Sufficient IQ samples!!! \r\n");
					#endif
				}

				#if AMR_DEBUG
					printf(" Run AMR #[%d] with %d samples \r\n",toaLoop, AMR_LEN );
				#endif

				/////////////////////////////////////////////////////////////////////////
				// Load IQ samples
				/////////////////////////////////////////////////////////////////////////	
			
				begin = amrObj[toaLoop].IQbeginIdx;
				end   = begin+AMR_LEN-1; // amrObj[toaLoop].IQendIdx; //begin + AMR_LEN; // 
				#if AMR_DEBUG
					printf("\r begin [%d], end[%d] \r\n",begin,end);
				#endif

				i=0;
				for(idx=begin; idx<=end; idx++){
					AMR_inputData[2*i+1] = I_BufferConvF[idx];
					AMR_inputData[2*i+2] = Q_BufferConvF[idx]; 
					i++;
				}

	end=0;
#ifdef activateCDW
				// Set SNR
				AMR_enable_awgn = 1;
				if(uiSNRfromGUI >= 40){
					AMR_enable_awgn = 0;
				}
				else if (uiSNRfromGUI <= 5){
					AMR_SNR_dB = (float)uiSNRfromGUI - 45.0f; 
				}else{
					AMR_SNR_dB = (float)uiSNRfromGUI;
				}

#else
				AMR_SNR_dB = 16; // dB
#endif
								 // values for SNR=15dB condition
								 // ASK  : 23 
							     // BPSK : 15
							     // FM : 13
								 // MSK : 16
				// Set frequency bins
				AMR_NFFT=setFrequencyVector(AMR_LEN);

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n unitPower begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
				// Force signal power to unit
				unitPower(AMR_inputData, AMR_LEN);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n unitPower end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n welchPSD begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
				// Calculate welch's PSD for bandwidth estimation
				welchPSD(AMR_inputData, AMR_LEN, AMR_NFFT, 1, AMR_powerSpectrum);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n welchPSD end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n coarseBandWidthEstimation begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
	
				// Bandwidth estimation
				AMR_coarseBandWidth = coarseBandWidthEstimation(AMR_powerSpectrum, AMR_NFFT);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n coarseBandWidthEstimation end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif
		
	#if AMR_PROFILE
		TSCL = 0; // need to write to it to start counting
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n awgn begin high : %u low : %u \r\n",startHigh, startLow);
	#endif					
				// Add white Gaussian noise to signal	
				awgn(AMR_inputData, AMR_SNR_dB,AMR_noisePower);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n awgn end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif
				// Check SNR
				// Calculate welch's PSD for bandwidth estimation
				welchPSD(AMR_inputData, AMR_LEN, AMR_NFFT, 1, AMR_powerSpectrum);

				// Bandwidth estimation
				AMR_coarseBandWidth = coarseBandWidthEstimation(AMR_powerSpectrum, AMR_NFFT);

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n coarsefreqOffsetEstimation begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
				offset = coarsefreqOffsetEstimation(AMR_LEN, AMR_NFFT);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n coarsefreqOffsetEstimation end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n lowPassFiltering begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
				lowPassFiltering(AMR_inputData);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n lowPassFiltering end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n linearVsNonLinearClassification begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
				linearVsNonLinearClassification(AMR_filteredInputData, AMR_LEN, AMR_instantAmpObservationInterval);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n linearVsNonLinearClassification end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n spectrumFeatureExt begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
				spectrumFeatureExt(AMR_filteredInputData, AMR_LEN, AMR_NFFT);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n spectrumFeatureExt end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n coarseSymRateEstimation begin high : %u low : %u \r\n",startHigh, startLow);
	#endif
				AMR_coarseSymbolRate = coarseSymRateEstimation(AMR_spectrumMagnitude, AMR_NFFT);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n coarseSymRateEstimation end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n FMvsFSKClassification begin high : %u low : %u \r\n",startHigh, startLow);
	#endif				
				//FMvsFSKClassification(AMR_spectrumMagnitude, AMR_symRateEstSegLen);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n FMvsFSKClassification end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n coarseToneSpacing begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
				// require to edit the program 
				coarseToneSpacing(AMR_filteredInputData, AMR_LEN, AMR_NFFT,AMR_toneSpacing_FFTsize,offset);  	// only for FSK digital signal

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n coarseToneSpacing end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

		
	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n reSampling end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif
				reSampling(AMR_LEN, AMR_coarseSymbolRate);

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n reSampling end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif

	#if AMR_PROFILE
		startLow = TSCL;
		startHigh = TSCH;
		printf("\r\n linearModClassification begin high : %u low : %u \r\n",startHigh, startLow);
	#endif	
				linearModClassification(AMR_filteredInputData, AMR_LEN, AMR_CC_samplesPerSegment);   // only for linear digital signal

	#if AMR_PROFILE
		stopLow = TSCL;
		stopHigh = TSCH; // stop will have the total number of CPU cycles
		printf("\r\n linearModClassification end high : %u low : %u \r\n",stopHigh, stopLow);
	#endif
				displayResults();

				///// add code to store the result to amrObj
				amrObj[toaLoop].coarseBandWidth = AMR_coarseBandWidth;
				amrObj[toaLoop].coarseFreqOffset = AMR_coarseFreqOffset;
				amrObj[toaLoop].coarseSymbolRate = AMR_coarseSymbolRate;
				amrObj[toaLoop].coarseToneSpacing = AMR_coarseToneSpacing;

				amrObj[toaLoop].modNum = AMR_modNum;
				amrObj[toaLoop].modOrder = AMR_modOrder;
				if( toaLoop == 0){
					#ifdef activateCDW
						LOCAL_send_CDW_AMR(FPGA_uiCDWindex);
					#endif
				}
			} // end -- toaLoop
		}else if(AMR_numTOA > 5){  // Hopping signal - analog

			if(AMR_numTOA>20) AMR_numTOA = 20; //Limit the repetition number

			for(toaLoop=0; toaLoop<AMR_numTOA; toaLoop++){
		
				initializeVariable();

				AMR_LEN = amrObj[toaLoop].IQnum;
				/////////////////////////////////////////////////////////////////////////
				// Determine whether the number of IQ samples is sufficient or not
				/////////////////////////////////////////////////////////////////////////

				// Is IQ samples valid?
				compConst = (AMR_MAXArray-1) >> 1;
				if (AMR_LEN > compConst){
					AMR_LEN = compConst;
					#if AMR_DEBUG
						printf("\r\n Too many IQ samples, the num. of IQ samples is restricted to %d!!! \r\n",AMR_LEN);
					#endif
				}
				else if (AMR_LEN < AMR_oneHopMinIQ)
				{
					#if AMR_DEBUG
						printf("\r\n Insufficient IQ samples %d!!! \r\n",AMR_LEN);
					#endif

					#ifdef activateCDW
						LOCAL_send_CDW_AMR(FPGA_uiCDWindex);
					#endif
					return 0;
				}
				else
				{
					//#if AMR_DEBUG
					//	printf("\r\n Sufficient IQ samples!!! \r\n");
					//#endif
				}

				#if AMR_DEBUG
					printf("\r\n Run AMR #[%d] with %d samples \r\n",toaLoop, AMR_LEN );
				#endif

				/////////////////////////////////////////////////////////////////////////
				// Load IQ samples
				/////////////////////////////////////////////////////////////////////////	
			
				begin = amrObj[toaLoop].IQbeginIdx;
				end   = begin + AMR_LEN; // amrObj[toaLoop].IQendIdx;
				//#if AMR_DEBUG
				//	printf("\r\n begin [%d], end[%d] \r\n",begin,end);
				//#endif

				i=0;
				for(idx=begin; idx<end; idx++){
					AMR_inputData[2*i+1] = I_BufferConvF[idx];
					AMR_inputData[2*i+2] = Q_BufferConvF[idx]; 
					i++;
				}
				// Add white Gaussian noise to signal	
				awgn(AMR_inputData, AMR_SNR_dB, AMR_noisePower);

				// Force signal power to unit
				unitPower(AMR_inputData, AMR_LEN);

				//#if AMR_DEBUG
				//	printf("\r\n #AMvsFM \r\n");
				//#endif


				//////////////////////////////////////////////////
				// Calculate the instantaneous phase
				//////////////////////////////////////////////////

				k=0;
				for(j=0 ; j<AMR_LEN ; j++){
					//AMR_iPhase[k] = 0.0f;
					AMR_iPhase[k] = atan2(AMR_inputData[2*j+2],AMR_inputData[2*j+1]); // wraped instantaneous phase
					AMR_ck[k] = 0.0f;
					k++;

				}

				// make range from -PI to +PI (=unwrap)
				AMR_ck[0] = 0.0f;
				for(j=1; j<AMR_LEN ; j++){
					iFreq = AMR_iPhase[j] - AMR_iPhase[j-1];
					if( iFreq > PI)
						AMR_ck[j] = AMR_ck[j-1] - TWOPI;
					else if ( iFreq < (-1.0f)*PI)
						AMR_ck[j] = AMR_ck[j-1] + TWOPI;
					else	AMR_ck[j] = AMR_ck[j-1];

					AMR_iPhase[j-1] += AMR_ck[j-1]; // unwraped instantaneous phase
				}
				AMR_iPhase[AMR_LEN-1] += AMR_ck[AMR_LEN-1]; // unwraped instantaneous phase

				 // Reuse array for linear curve fitting
				for(j=0; j<AMR_LEN; j++)	AMR_spectrumBuffer[j] = (float)j;
				chi2 = 0.0f;
				fit(AMR_spectrumBuffer, AMR_iPhase, AMR_LEN, &a, &b, &siga, &sigb, &chi2);

				if( chi2 >= 0.5f ) fmCnt++;	// refer to AMvsFM.m in MATLAB
				else amCnt++;
				
			} // for(toaLoop=0; toaLoop<AMR_numTOA; toaLoop++
#if AMR_DEBUG
				printf("\r\n fmCnt %u \r\n",fmCnt);
				printf("\r\n amCnt %u \r\n",amCnt);
#endif

			if (fmCnt > amCnt){ 
				AMR_modNum = 1;		// FM
				AMR_coarseBandWidth = 25000;
			}
			else{ 
				AMR_modNum   = 7;    ////////////// AM
				AMR_coarseBandWidth = 15000;
			}

			amrObj[toaLoop].coarseBandWidth = AMR_coarseBandWidth;
			amrObj[toaLoop].coarseFreqOffset = 0.0f;
			amrObj[toaLoop].coarseSymbolRate = 0.0f;
			amrObj[toaLoop].coarseToneSpacing = 0.0f;

			amrObj[toaLoop].modNum = AMR_modNum;
			amrObj[toaLoop].modOrder = 0;

			#ifdef activateCDW
				LOCAL_send_CDW_AMR(FPGA_uiCDWindex);
			#endif
			
		} // else if(AMR_numTOA > 5)
		return 0;
	}
#endif


//#if VISUALSTUDIO
//	fa=fopen("filename.txt","w");
//	fprintf(fa,"data =[");
//	for(i=0; i<AMR_BW_SamplesPerSegment; i++) fprintf(fa,"%0.4f,...\n",(AMR_powerSpectrum[i]));
//	fprintf(fa,"]; \n");
//	fprintf(fa,"figure; plot(data)");
//	fclose(fa);			
//#endif
	
	
