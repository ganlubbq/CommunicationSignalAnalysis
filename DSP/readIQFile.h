#ifndef _READDIQFILE
#define _READDIQFILE

#if AMR_Cplus
	inline unsigned int AMRalgorithm::LOCAL_load_IQ(char *fname)
#else
	int LOCAL_load_IQ(char *fname)
#endif
{
	/********************************************************************************/
	// \ DESCRIPTION
	//  > Read IQ data from *.diq file
	//    You should write the filename below code
	//    a=fopen("filename.diq","rb");
	//
	// \ INPUT ARGUMENTS
	//	> none
	//
	// \ OUTPUT ARGUMENTS
	//
	// \ Example
	//
	// \ See also 
	//
	// \ Author(s) : AWH
	// 
	// \ Reference

	/******************************************************************************/

	unsigned int i,j,k, ddr3addr=0;
	unsigned int index = 0;
	unsigned int iq_num=0;
	unsigned int toa_num=0;
	unsigned int remain_num=0;
	unsigned int mod_num=0;
	unsigned int numIQFPGA= 0;
	unsigned int numIQ=0;

	unsigned int nowTOAidx=0;
	unsigned int oldTOAidx=0;
	unsigned int numTOA = 0;
	unsigned int result = 0;
	unsigned int skip   = 0;
	unsigned int iqIdx = 0;
	unsigned int offset = 0;

	unsigned short q_tmp=0, i_tmp=0;
	unsigned short i_sig=0, q_sig=0;
	
	
	// To read the *.diq file
	memset(IQ_Buffer, 0.0f, sizeof(IQ_Buffer));
	memset(I_BufferConvF, 0.0f, sizeof(I_BufferConvF));
	memset(Q_BufferConvF, 0.0f, sizeof(Q_BufferConvF));

	//fa=fopen("test.diq","rb");
	//fa=fopen("20150407_105734(BPSK_PN20_30ksps_256deci).diq","rb");
	fa=fopen(fname,"rb");
	//fa = fopen("20150407_105734(BPSK_PN20_30ksps_256deci).diq","rb");
	if (fa == NULL) {printf("\r\n File error \r\n"); exit(1);}
	// obtain file size;
	fseek(fa, 0, SEEK_END);
	numIQFPGA = ftell(fa);
	rewind(fa);
	numIQFPGA = numIQFPGA / 4;
	
	// read 4 byte by numIQFPGA times
	//numIQFPGA = 240000;
	result = fread(IQ_Buffer,4,numIQFPGA,fa) ;
	if( result != numIQFPGA) {printf("\r\n Reading error \r\n"); exit(2);}

	fclose(fa);

	/* Printout debug log */	
	printf("=>LOCAL_load_IQ()");		
	/******************/
	// Clear DDR2 Data Buffer

/*%%%%%%%%%%%%%%%%%%%%%%%% IQ data structure %%%%%%%%%%%%%%%%%%%%%%%%
% Type      |  header  |        Inphase      |       Quadrature   |     
% Bit index |  31   30 |     29     ~    15  |   14      ~      0 |           
%           | TOA : 11 | (sign bit)          | (sign bit)         |
%           | IQ  : 01 |                     |                    |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/ 

// Flip data by 16byte 

	numIQ =0; j=0; k=0; 
	//Discard the front 16byte
	for(i=4; i<numIQFPGA; i++){
		j = k / 4;
		IQ_BufferConv[i-4] = IQ_Buffer[11 + 8*j - i];	
		k++;
	}
	numIQFPGA = numIQFPGA - 4; AMR_numIQFPGA = numIQFPGA;
	for(i=0; i<numIQFPGA; i++)
	{
		if ((IQ_BufferConv[i]&0xc0000000)==0xc0000000) // TOA
		{
			nowTOAidx = i;
			if(numTOA > 0)
			{
                if ((nowTOAidx - oldTOAidx) < AMR_oneHopMinIQ)// # of IQ samples in a segemnt
                {   
                    oldTOAidx = nowTOAidx;  // skip
				}
                else
				{
                    beginTOAidx[numTOA-1] = oldTOAidx;   
                    endTOAidx[numTOA-1] = nowTOAidx;
                    oldTOAidx = nowTOAidx;
                    numTOA++;
				}
			}
			else
			{
				oldTOAidx = nowTOAidx;
                numTOA++;
			}
		} // end      for if ((IQ_BufferConv[i]&0xc0000000)==0xc0000000) // TOA
	} // end for      for(i=0; i<numIQFPGA; i++)

	// When the next TOA is not available
    if (numTOA > 0){ 
		if (nowTOAidx == oldTOAidx){                      
			if ((numIQFPGA - nowTOAidx) > AMR_oneHopMinIQ)
			{
				beginTOAidx[numTOA-1] = nowTOAidx;   
				endTOAidx[numTOA-1] = numIQFPGA;            
			}
			else 
			{
				numTOA--;
			}
		}

		offset = 8;
		for(j=0; j<numTOA; j++){
			numIQ = 0;
			for(i=beginTOAidx[j]+offset; i<endTOAidx[j]-offset+1; i++) // skip front 
			{ 
				//check IQ header
				if((IQ_BufferConv[i]&0xc0000000)==0x40000000) // IQ Data
				{
					i_tmp = (unsigned short)((IQ_BufferConv[i]&0x3fff8000)>>15); //uppper 15bits
					q_tmp = (unsigned short)((IQ_BufferConv[i]&0x7FFF));			// lowwer 15bits

					//2의 보수 처리
					if ((i_tmp & 0x4000) == 0x4000)
					{
							I_BufferConvF[iqIdx] = i_tmp&0x3fff;
							I_BufferConvF[iqIdx] = (-1.0f)*(((float)0x3fff-I_BufferConvF[iqIdx])+1);
					}
					else
					{
							I_BufferConvF[iqIdx] = i_tmp&0x3fff;
					}

					if ((q_tmp & 0x4000) == 0x4000)
					{
							Q_BufferConvF[iqIdx] = q_tmp&0x3fff;
							Q_BufferConvF[iqIdx] = (-1.0f)*(((float)0x3fff-Q_BufferConvF[iqIdx])+1.0f);
					}
					else
					{
							Q_BufferConvF[iqIdx] = q_tmp&0x3fff;
					}				
					iqIdx++; numIQ++;
				}
				else // tresh samples, skip this
				{		 		

				}
			} // end for	for(i=endTOAidx[i]+5; i<beginTOAidx[j]-1; i++)

			// assign data
			#if AMR_Cplus
				IQbeginIdx[j] = iqIdx - numIQ; 
				IQendIdx[j]   = iqIdx;		
				IQnum[j]      = numIQ;
			#else
				amrObj[j].IQbeginIdx = iqIdx - numIQ; 
				amrObj[j].IQendIdx   = iqIdx;			 
				amrObj[j].IQnum		 = numIQ;
			#endif


		} // end for   for(i=0; i<numTOA; i++)
		return numTOA;
	}
	else
	{
		return 0;
	}
}

#endif