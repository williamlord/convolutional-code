/*! 
  \file
  \brief Testing of Convolutional codes (CONV) class
  \author Michael Ng

  1.02

  2009/09/02
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "itpp/itbase.h"
#include "itpp/itcomm.h"
#include "conv.h"
#include "fileser.h"

/*! 
  \brief This function test the CONV module
*/

cvec my_channel(int channel_type, int N, double path_factor, TDL_Channel &ray_channel);

int main(int argc, char *argv[])
{
	int i,j,k;			//Counter variable
	
	bvec b_source_bits, b_decoded_bits, b_encoded_bits;
	ivec i_decoded_bits, i_decoded_symbols;
	cvec c_received_signals;
	vec  d_received_llr;
	ivec No_of_Errors, No_of_Bits, No_of_BlockErrors, No_of_Blocks;
	
	//Channel setup
	cvec channel_gains;
	TDL_Channel ray_channel; // default: uncorrelated Rayleigh fading channel
        AWGN_Channel awgn_channel; // AWGN channel
	
	//-------------------------------------------------------CONV
	//Convolutional codes module, CONV
	CONV conv;
	char *ctrlfile;
	if(argc==1) ctrlfile = "control.conv.par";
	else        ctrlfile = argv[1];
	
	conv.set_parameters_from_file(ctrlfile);
	conv.initialise();
	
	//result file
	FILE *f;  char *resultfile= conv.get_result_filename();
	f=fopen(resultfile, "w");
	
	//print out parameters
	conv.print_parameters(stdout);    
	conv.print_parameters(f);
	//-------------------------------------------------------CONV	
	
	// read simulation parameters
	FILE *sim_file;
	FILESERVICE fileser;
	sim_file = fopen("control.sim.par", "r");
	if(sim_file==NULL) it_error("control.sim.par not found");
	int max_nrof_frame = fileser.scan_integer(sim_file);
	double SNR_start  = fileser.scan_double(sim_file); 
	double SNR_step   = fileser.scan_double(sim_file);
	double SNR_end    = fileser.scan_double(sim_file);
	double G_sd        = fileser.scan_double(sim_file);
	double G_sr        = fileser.scan_double(sim_file);
	double G_rd        = fileser.scan_double(sim_file);
	int channel_type   = fileser.scan_integer(sim_file);
	int mapper_bps     = fileser.scan_integer(sim_file);
	fclose(sim_file);
	char text[100];
	sprintf( text, "%f:%f:%f", SNR_start, SNR_step, SNR_end );
	
	//SNR setup
	vec SNRdB = text;	//Setting the simulation Eb/N0 range in dB


	int M = 1<<mapper_bps;
        PSK psk(M);
	cvec source_frame;
	
	double rate = (double)mapper_bps * conv.get_rate(); 

	printf (  "! %d-PSK (mapper_bps=%d) :: Overall_Rate=%.4f\n!\n", M, mapper_bps, rate);
	fprintf(f,"! %d-PSK (mapper_bps=%d) :: Overall_Rate=%.4f\n!\n", M, mapper_bps, rate);

	vec SNR    = pow(10.0, SNRdB/10.0);
	vec N0     = 1.0/SNR;
	vec sigma2 = N0/2;

	
	BERC berc;			//BER counter
	BLERC blerc;			//FER counter
	Real_Timer tt;			//Time counter

	vec ber(SNRdB.length());	//allocate memory for vector to store BER
	vec bler(SNRdB.length());	//allocate memory for vector to store FER
	
	ber.clear();			//Clear up buffer of BER counter
	bler.clear();			//Clear up buffer of FER counter	
	
	blerc.set_blocksize((long)conv.get_info_bit_length());	//set blocksize of the FER counter
	
	tt.tic();					//Start timer
	
	//RNG_randomize();				//construct random source 
	RNG_reset(0);                                   //reset random seed 
	
	b_source_bits.set_size(conv.get_info_bit_length(), false);
	b_encoded_bits.set_size(conv.get_coded_bit_length(), false);
	c_received_signals.set_size(conv.get_sym_length(), false);
	d_received_llr.set_size(conv.get_coded_bit_length(), false);
	b_decoded_bits.set_size(conv.get_info_bit_length(), false);
	
	No_of_Errors.set_size(SNRdB.length(), false);		//Set the length
	No_of_Bits.set_size(SNRdB.length(), false);		//for ivectors storing the no	
	No_of_BlockErrors.set_size(SNRdB.length(),false);	//of errors bits or error frames
	No_of_Blocks.set_size(SNRdB.length(),false);
			
	No_of_Errors.clear();
	No_of_Bits.clear();
	No_of_BlockErrors.clear();
	No_of_Blocks.clear();
	
	printf (  "!SNR(dB)\tEbN0(dB)\tBER\t\tFER\t\tnrof Frames\n");
	fprintf(f,"!SNR(dB)\tEbN0(dB)\tBER\t\tFER\t\tnrof Frames\n");
	
	for(i=0; i< SNRdB.length(); i++)
	{
		//Set channel noise level
		awgn_channel.set_noise(N0(i));
		
		for(j=0;j<max_nrof_frame;j++)
		{			
		        //Generate random source bits
			b_source_bits = randb(conv.get_info_bit_length());
			
			//CONV encode
			conv.encode_bits(b_source_bits, b_encoded_bits);
			
			source_frame = psk.modulate_bits(b_encoded_bits);
			
			// Fast Rayleigh channel + AWGN transmission 
			channel_gains = my_channel(channel_type, source_frame.length(), G_sd, ray_channel); 
			c_received_signals = elem_mult(source_frame, channel_gains);
			c_received_signals = awgn_channel( c_received_signals );
			
			//Demodulation to get llr
			psk.demodulate_soft_bits(c_received_signals, channel_gains, N0(i), d_received_llr);
			//Convert to Log(Pr=1/Pr=0)
                        d_received_llr = -1 * d_received_llr;

			//CONV decode
			conv.assign_apr_codeword(d_received_llr);
			conv.decode(i_decoded_bits, i_decoded_symbols);

			b_decoded_bits = to_bvec(i_decoded_bits.left(conv.get_info_bit_length())); 
			// remove dummy bits if trellis/code termination is used
			
			berc.clear();
			blerc.clear();
			
			berc.count (b_source_bits, b_decoded_bits);	//Count error bits in a word
			blerc.count (b_source_bits, b_decoded_bits);	//Count frame errors
			
			No_of_Errors(i) += berc.get_errors();		//Updating counters	
			No_of_Bits(i) += berc.get_errors()+berc.get_corrects();
			No_of_BlockErrors(i) +=blerc.get_errors();
			No_of_Blocks(i)++;
			
			if(No_of_Errors(i)>100000)
				break;
		}
		
		ber(i) = (double)No_of_Errors(i)/No_of_Bits(i);
		bler(i) = (double)No_of_BlockErrors(i)/No_of_Blocks(i);
		
		double EbN0dB = SNRdB(i) - 10*log10(rate);
		
		printf("%f\t%f\t%e\t%e\t%d\n",    SNRdB(i), EbN0dB, ber(i), bler(i), No_of_Blocks(i));
		fprintf(f,"%f\t%f\t%e\t%e\t%d\n", SNRdB(i), EbN0dB, ber(i), bler(i), No_of_Blocks(i));
		
		if(ber(i)<1e-5)  break;			
	}

	fprintf(f,"!Elapsed time = %d s\n", tt.get_time());		//output simulation time
	tt.toc();							//Stop timer and output simulation time
	
	fclose(f);							//close output file	
	return 0 ;							//exit program
}


cvec my_channel(int channel_type, int N, double path_factor, TDL_Channel &ray_channel){
    Array<cvec> channel_gains;
    int k;

    if(channel_type==0){// fast fading
        ray_channel.generate(N, channel_gains);
    }
    else if(channel_type==1){//quasi-static fading
        ray_channel.generate(1, channel_gains);
        channel_gains(0).set_size(N,true);
        for(k=0;k<N;k++) channel_gains(0)[k] = channel_gains(0)[0];
    }
    channel_gains(0) = sqrt(path_factor)*channel_gains(0);

    return channel_gains(0);
}














