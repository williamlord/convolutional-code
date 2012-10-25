/*! 
  \file 
  \brief Implementation of Convolutional codes (CONV) class
  \author Michael Ng

  1.02

  2009/09/02
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <itpp/itcomm.h>

#include "conv.h"
#include "fileser.h"

// this is just to satisfiy the new compiler
double pow(int a, int b){
  return pow((double) a, (double) b);
}



CONV::~CONV()
{
}


/*-----------CONV parameter setup----------------------------------*/
void CONV::set_parameters_from_file(char *fname)
{
  FILE *f;
  FILESERVICE fileser;

  f = fopen(fname, "r");
  if(f==NULL) { sprintf(fname,"%s not found",fname); it_error(fname);}
  
  int mode_in = fileser.scan_integer(f);
  int k_in    = fileser.scan_integer(f);
  int n_in    = fileser.scan_integer(f);
  int L_in    = fileser.scan_integer(f);
  int Terminated_in = fileser.scan_integer(f); 
  int decoder_type_in     = fileser.scan_integer(f);
  int no_of_symbols_in    = fileser.scan_integer(f);
  fileser.scan_text(f, result_filename);
  
  set_parameters(mode_in, k_in, n_in, L_in, Terminated_in, 
		 decoder_type_in, no_of_symbols_in);
		 //modulation_type_in, interleaver_mode_in);
  
  fclose(f);
}

void CONV::set_parameters(const int mode_in, const int k_in, const int n_in, const int L_in, const int Terminated_in, 
			  const int decoder_type_in, const int no_of_symbols_in)
{
  mode=mode_in;

  if(mode!=0 && mode!=NSC && mode!=SC && mode!=RSC && mode!= RNSC) s_error("mode: Uncoded, NSC, SC, RSC or RNSC");
  
  if(mode==0){ //uncoded 
      k=n=1;
      L=0;
  }
  else{
      k=k_in;
      n=n_in;
      L=L_in;
  }

  Terminated=Terminated_in;
  if(Terminated!=ON && Terminated!=OFF) s_error("coder termination: ON or OFF");
  
  decoder_type=decoder_type_in;
  if(decoder_type!=C_EXACT_LOG_MAP && decoder_type!=C_APPROX_LOG_MAP && decoder_type!=C_MAX_LOG_MAP )
    s_error("decoder type: C_EXACT_LOG_MAP, C_APPROX_LOG_MAP or C_MAX_LOG_MAP");
  
  max_no_of_symbols = no_of_symbols_in;
  no_of_symbols     = no_of_symbols_in;
  
  if(Terminated==1){ 
    no_of_info_symbols=no_of_symbols_in-L;
  }
  else no_of_info_symbols=no_of_symbols_in;
  
  no_of_input_bits = no_of_symbols * k;
  no_of_info_bits = no_of_info_symbols * k;
  no_of_coded_bits = no_of_symbols * n;
  no_of_tail_bits = (no_of_symbols - no_of_info_symbols) * k;
  
}

void CONV::print_parameters(FILE *fp)
{
  int i,j;
  
  if(mode==0){
      fprintf(fp, "! Uncoded \n"); 
      return;
  }

  fprintf(fp, "!Convolutional codes:\n");
  if(mode==NSC)       fprintf(fp, "! NSC: Non Systematic Convolutional code (k=%d, n=%d), code memory=%d\n", k,n,L);
  else if(mode==SC)   fprintf(fp, "! SC: Systematic Convolutional code (k=%d, n=%d), code memory=%d\n", k,n,L);
  else if(mode==RSC)  fprintf(fp, "! RSC: Recursive Systematic Convolutional code (TCM type) (k=%d, n=%d), code memory=%d\n", k,n,L);
  else if(mode==RNSC) fprintf(fp, "! RNSC: Recursive Systematic Convolutional code (from NSC) (k=%d, n=%d), code memory=%d\n", k,n,L);
  
  if(decoder_type==C_EXACT_LOG_MAP)       fprintf(fp, "!  Decoder type: exact Log MAP\n");
  else if(decoder_type==C_APPROX_LOG_MAP) fprintf(fp, "!  Decoder type: approximated Log MAP\n");
  else if(decoder_type==C_MAX_LOG_MAP)    fprintf(fp, "!  Decoder type: maximum Log MAP\n");
  
  if(mode==SC || mode==RSC){
    int i;
    fprintf(fp, "!  Gen_poly = ");
    for(i=0;i<n;i++) fprintf(fp, "%d ", trellis->GenPoly[i]);
    fprintf(fp,"\n");
  }
  else if(mode==NSC || mode==RNSC){
    for(i=0;i<k;i++){
	fprintf(fp, "!  Gen_poly[%d] = ",i);
	for(j=0;j<n;j++) fprintf(fp,"%d ", trellis->GenPolyMat(i,j));
	fprintf(fp,"\n");
    }
  }
  
  fprintf(fp, "!  Code termination: %s\n", Terminated==ON ? "ON":"OFF");
  fprintf(fp, "!  Block length: output coded (%d-bit) symbol (%d), input data (%d-bit) symbol (%d)\n", n, no_of_symbols, k, no_of_info_symbols);
  fprintf(fp, "!  Bit Block length: output coded bits (%d), input data bits (%d)\n", n*no_of_symbols, k*no_of_info_symbols);
  

  /*
  if(interleaver_mode==OFF)        fprintf(fp, "!  Channel interleaver: OFF\n");
  else if(interleaver_mode==BIT)   fprintf(fp, "!  Channel interleaver: one BIT-based random interleaver\n");
  else if(interleaver_mode==PBIT)  fprintf(fp, "!  Channel interleaver: %d Parallel BIT-based random interleaver\n",n);
  else if(interleaver_mode==SYM)   fprintf(fp, "!  Channel interleaver: SYMbol-based random interleaver\n");
  else if(interleaver_mode==IQ_SYM)fprintf(fp, "!  Channel interleaver: IQ SYMbol-based random interleaver\n");
  */

  fprintf(fp, "!\n");
}

void CONV::print_coding_tables(FILE *fp)
{
  int i, j;
  
  fprintf(fp, "CurrentState Branch : Codeword NextState PreviousState\n"); 
  for(i=0;i<trellis->nrof_states;i++){
      for(j=0;j<trellis->nrof_codewords;j++){
	
	if(trellis->MaxPs(i,j)==NO_LINK)
	  fprintf(fp, "%d %d : %d %d NO_LINK", 
		  i, j, trellis->Lb(i,j), trellis->Ns(i,j));
	else	  
	  fprintf(fp, "%d %d : %d %d %d", 
		  i, j, trellis->Lb(i,j), trellis->Ns(i,j), trellis->MaxPs(i,j));

	fprintf(fp, "\n");
      }
  }
  //cout << trellis->Lb << endl;
  //cout << trellis->Ns << endl;
  //cout << trellis->MaxPs << endl;
}


int CONV::encode( bvec b_Input_bits, ivec &i_Output_symbols )
{
  ivec i_Input_bits, i_Input_symbols;

  if(mode==0) it_error("Uncoded system use encode_bits(bvec b_Input_bits, bvec &b_Output_bits)");
  
  if(b_Input_bits.length()!=no_of_info_bits){ 
    printf("b_Input_bits.length()=%d no_of_info_bits=%d\n",b_Input_bits.length(), no_of_info_bits);
    s_error("CONV::encode: make sure that b_Input_bits.length()==no_of_info_bits ");
  }
  
  i_Input_bits    = to_ivec(b_Input_bits);
  
  if(Terminated==ON){ 
    /* Insert tail bits */
    ivec tail_bits;
    tail_bits.set_length(no_of_tail_bits, false);
    tail_bits.zeros();
    i_Input_bits = concat(i_Input_bits,tail_bits); 
    //tail_bits.~ivec();
  }

  i_Input_symbols = bits_seq_to_symbol(k, no_of_symbols, i_Input_bits); 

#ifdef debug_conv  
  cout << "i_Input_bits   (" << i_Input_bits.length() <<")= " << i_Input_bits << endl;
  cout << "i_Input_symbols(" << i_Input_symbols.length() << ")= " << i_Input_symbols << endl;
#endif

  i_Output_symbols.set_size(no_of_symbols,false);


  //----------------------
  int i,m,s;
  for(i=s=0; i<no_of_symbols; i++){
    m = i_Input_symbols[i];
    i_Output_symbols[i] = trellis->Lb(s,m); // Lb = codeword table
    s = trellis->Ns(s,m);                   // Ns = next state table
  }
  //----------------------
  
  return 0;
}

int CONV::encode_bits( bvec b_Input_bits, bvec &b_Output_bits )
{
  ivec i_Output_symbols(no_of_symbols);

  if(mode==0){ // uncoded
      b_Output_bits = b_Input_bits;
      return 0;
  }
  
  encode(b_Input_bits, i_Output_symbols);
  b_Output_bits = to_bvec( symbol_to_bits_seq(n, no_of_symbols, i_Output_symbols) );
  
  return 0;
}

ivec CONV::encode_bits( bvec b_Input_bits)
{
  ivec i_Output_symbols(no_of_symbols);
  
  // uncoded
  if(mode==0) return to_ivec(b_Input_bits);
    
  encode(b_Input_bits, i_Output_symbols);
  return symbol_to_bits_seq(n, no_of_symbols, i_Output_symbols);
}

void CONV::decode_symbol(mat Apo, ivec &symbols, int N, int M)
{
    int k,i,m;
    double max;

    for(k=0; k<N; k++)
    {
        i=0; max=Apo(k,0);
        for(m=1; m<M; m++)
        {    
            if(Apo(k,m) > max)
            {
                max = Apo(k,m);
                i = m;
            }
        }
        symbols[k] = i;
    }
}


int CONV::decode(ivec &i_decoded_bits, ivec &i_decoded_symbols)
{
        int N, M ,S;
        //int j, m, i;

        N  = no_of_symbols;  
        M  = 1<<k;
        S  = 1<<L;

	//b_decoded_bits.set_size(no_of_info_bits,false);
	i_decoded_bits.set_size(no_of_info_bits,false);
	if(mode==0) //uncoded
	    map->Apo_dataword = map->Apr_codeword;
	if(mode==NSC)
	  SISO_dec_log(no_of_symbols, 
		       trellis->nrof_datawords, trellis->nrof_codewords,
		       trellis->nrof_states, trellis->MaxPs, trellis->Ns, trellis->Lb, 
		       map->Apr_codeword, map->Apr_dataword,
		       map->Apo_codeword, map->Apo_dataword, 
		       Terminated, 0);	
	else if(mode==RNSC)
	  SISO_dec_log_RNSC(no_of_symbols, 
		       trellis->nrof_datawords, trellis->nrof_codewords,
		       trellis->nrof_states, trellis->ptrPs, trellis->Ns, trellis->Lb, 
		       map->Apr_codeword, map->Apr_dataword,
		       map->Apo_codeword, map->Apo_dataword, 
		       Terminated, 0);	
	else if(mode==RSC)
	  SISO_dec_log_RSC(no_of_symbols, 
		       trellis->nrof_datawords, trellis->nrof_codewords,
		       trellis->nrof_states, trellis->Ps, trellis->Ns, trellis->Lb, 
		       map->Apr_dataword, map->Apr_codeword, 
		       map->Apo_dataword, map->Apo_codeword);
	
	// Final decoding
	i_decoded_symbols.set_length(no_of_symbols,false);
	decode_symbol(map->Apo_dataword, i_decoded_symbols, no_of_symbols, M);
	i_decoded_bits = symbol_to_bits_seq(k, i_decoded_symbols.length(), i_decoded_symbols);		
	//b_decoded_bits = to_bvec(i_decoded_bits.left(no_of_info_bits));
	
	return 1; 
}


void CONV::SISO_dec_log(int N, int D, int C,
			int S, imat Ps, imat Ns, imat Lb, 
			mat Apr_codeword, mat Apr_dataword,
			mat &Apo_codeword, mat &Apo_dataword, 
			int Terminated, int frame_index)
{
        int i, k, m, j, c;
        mat alpha, beta;
        double abc, max, codeword_max;
	
        alpha = mat(N+1, S);
        beta  = mat(N+1, S);
	
        /* compute and normalise alpha */
        for(i=0; i<S; i++) alpha(0,i)=MINF;
        alpha(0,0)=0.;
        for(k=1; k<=N; k++) 
        {   max = MINF;
            for(i=0; i<S; i++)
            {
                  alpha(k,i)=MINF;
                  for( c=0; c<C; c++){
		      if(Ps(i,c) != NO_LINK){

			  m=NO_LINK;
			  for(j=0;j<D;j++){ if(c==Lb( Ps(i,c) ,j)) { m=j;break;} } 
			  if(m==NO_LINK) it_error("Check code table!");
			  
                          alpha(k,i) = jacolog( alpha(k,i), 
						alpha(k-1,Ps(i,c)) + Apr_codeword(k-1,c) + Apr_dataword(k-1,m));
		      }
                  }
                  if(max < alpha(k,i)) max = alpha(k,i);
            }
            for(i=0; i<S; i++) alpha(k,i) -= max;
        }

        /* compute and normalise beta */
        for(i=0; i<S; i++) beta(N,i)=0.;
        if(Terminated==ON&&mode==NSC) for(i=1; i<S; i++) beta(N,i)=MINF;

        for(k=N-1; k>=0; k--) 
        {   max = MINF;
            for(i=0; i<S; i++)
            { 
                 beta(k,i)=MINF;
                 for( m=0; m<D; m++)
                     beta(k,i) = jacolog( beta(k,i),
					  beta(k+1,Ns(i,m)) + Apr_codeword(k, Lb(i,m) ) +  Apr_dataword(k,m) );
		 
                 if(max < beta(k,i)) max = beta(k,i);
            }
            for(i=0; i<S; i++) beta(k,i) -= max;
        }

        /* compute and normalise Apo_dataword & Apo_codeword */
        for(k=0; k<N; k++) 
        {
            max = MINF; for(m=0;m<D;m++) Apo_dataword(k,m) = MINF;
	    codeword_max = MINF; for(c=0;c<C;c++) Apo_codeword(k,c) = MINF;
	    
            for(i=0; i<S; i++){
                for(c=0; c<C; c++){
                        
                    if( Ps(i,c)!=NO_LINK ){
		      
		        m=NO_LINK;
                        for(j=0;j<D;j++){ if(c==Lb( Ps(i,c) ,j)) { m=j;break;} }
                        if(m==NO_LINK) it_error("Check code table!");

			abc = alpha(k,Ps(i,c)) + beta(k+1,i) + Apr_codeword(k, c);

                        Apo_dataword(k,m) =
                            jacolog(Apo_dataword(k,m),
				    abc + Apr_dataword(k,m));				
                        if(max < Apo_dataword(k,m)) max = Apo_dataword(k,m);
			
			Apo_codeword(k,c) =
			    jacolog(Apo_codeword(k,c),
				    abc + Apr_dataword(k,m));			
			if(codeword_max < Apo_codeword(k,c)) codeword_max = Apo_codeword(k,c);
		    }
                }
            }
            for(m=0;m<D;m++)  Apo_dataword(k,m) -= max;
	    for(c=0;c<C;c++)  Apo_codeword(k,c) -= codeword_max;
        }
	
	return;
}

void CONV::SISO_dec_log_RNSC(int N, int D, int C,
			int S, struct trellis_state **ptrPs, imat Ns, imat Lb, 
			mat Apr_codeword, mat Apr_dataword,
			mat &Apo_codeword, mat &Apo_dataword, 
			int Terminated, int frame_index)
{
        int i, k, m, j, c;
        mat alpha, beta;
        double abc, max, codeword_max;

	int iw, w;
        alpha = mat(N+1, S);
        beta  = mat(N+1, S);
	
        /* compute and normalise alpha */
        for(i=0; i<S; i++) alpha(0,i)=MINF;
        alpha(0,0)=0.;
        for(k=1; k<=N; k++) 
        {   max = MINF;
            for(i=0; i<S; i++)
            {
                  alpha(k,i)=MINF;
		  
		  for(m=0;m<C;m++){
		    iw = ptrPs[i][m].infoword;     
		    w  = ptrPs[i][m].word;
		    
		    if(w == NO_LINK) break;
		    
		    alpha(k,i) = jacolog( alpha(k,i),
					  alpha(k-1, ptrPs[i][m].index) +
					  Apr_codeword(k-1,w) + Apr_dataword(k-1,iw));
		  }  
                  if(max < alpha(k,i)) max = alpha(k,i);
            }
            for(i=0; i<S; i++) alpha(k,i) -= max;
        }

        /* compute and normalise beta */
        for(i=0; i<S; i++) beta(N,i)=0.;

        for(k=N-1; k>=0; k--) 
        {   max = MINF;
            for(i=0; i<S; i++)
            { 
                 beta(k,i)=MINF;
                 for( m=0; m<D; m++)
                     beta(k,i) = jacolog( beta(k,i),
					  beta(k+1,Ns(i,m)) + Apr_codeword(k, Lb(i,m) ) +  Apr_dataword(k,m) );
		 
                 if(max < beta(k,i)) max = beta(k,i);
            }
            for(i=0; i<S; i++) beta(k,i) -= max;
        }

        /* compute and normalise Apo_dataword & Apo_codeword */
        for(k=0; k<N; k++) 
        {
            max = MINF; for(m=0;m<D;m++) Apo_dataword(k,m) = MINF;
	    codeword_max = MINF; for(c=0;c<C;c++) Apo_codeword(k,c) = MINF;
	    
            for(i=0; i<S; i++){
	      for(m=0;m<C;m++){
		iw = ptrPs[i][m].infoword;
		w  = ptrPs[i][m].word;
		
		if(w == NO_LINK) break;

		abc = alpha(k,ptrPs[i][m].index) + beta(k+1,i) + Apr_codeword(k, w);
		
		Apo_dataword(k,iw) =
		  jacolog(Apo_dataword(k,iw), 
			  abc + Apr_dataword(k,iw));	       
		if(max < Apo_dataword(k,iw)) max = Apo_dataword(k,iw);

		Apo_codeword(k,w) =
		  jacolog(Apo_codeword(k,w),
			  abc + Apr_dataword(k,iw));
		if(codeword_max < Apo_codeword(k,w)) codeword_max = Apo_codeword(k,w);
	      }
	    }
	    
            for(m=0;m<D;m++)  Apo_dataword(k,m) -= max;
	    for(c=0;c<C;c++)  Apo_codeword(k,c) -= codeword_max;
        }
	
	return;
}

void CONV::SISO_dec_log_RSC(int N, int M, int C, int S, imat Ps, imat Ns, imat Lb,
			    mat Apr_dataword, mat Apr_codeword, mat &Apo, mat &Apo_codeword)
{
        int i, k, m, c;
        mat alpha, beta;  
        double max=MINF, abc, max_Apo_codeword=MINF;
        
        alpha = mat(N+1,S+1);
        beta  = mat(N+1,S+1);
        
        /* compute alpha */
        for(i=0; i<S; i++) alpha(0,i)=MINF;
        alpha(0,0)=0.;
        for(k=1; k<=N; k++) 
        {       
                max = MINF;
                for(i=0; i<S; i++)
                {
                        /* init: m=0 and 1 */
		        alpha(k,i) = jacolog( alpha(k-1, Ps(i,0)) + Apr_codeword(k-1, Lb(Ps(i,0),0)) + Apr_dataword(k-1,0),
					      alpha(k-1,Ps(i,1)) + Apr_codeword(k-1, Lb(Ps(i,1),1)) + Apr_dataword(k-1,1));
                        for( m=2; m<M; m++)
                            alpha(k,i) = jacolog( alpha(k,i), 
						  alpha(k-1,Ps(i,m)) + Apr_codeword(k-1, Lb(Ps(i,m),m)) + Apr_dataword(k-1,m));
                              
                        if ( max < alpha(k,i) ) max = alpha(k,i);
                }
                for(i=0; i<S; i++) alpha(k,i) -= max;
        }

        /* compute beta */
        for(i=0; i<S; i++) beta(N,i)=0.;
        for(k=N-1; k>=0; k--) 
        {  
                max = MINF;
                for(i=0; i<S; i++)
                {
                        beta(k,i)= jacolog( beta(k+1,Ns(i,0)) + Apr_codeword(k, Lb(i,0)) + Apr_dataword(k,0),
					    beta(k+1,Ns(i,1)) + Apr_codeword(k, Lb(i,1)) + Apr_dataword(k,1));
                        for( m=2; m<M; m++) 
                            beta(k,i) = jacolog( beta(k,i), 
						 beta(k+1,Ns(i,m)) + Apr_codeword(k, Lb(i,m)) + Apr_dataword(k,m));
                        
                        if ( max < beta(k,i) ) max = beta(k,i);
                }
                for(i=0; i<S; i++) beta(k,i) -= max;
        }

        /* compute apo */
        for(k=0; k<N; k++) 
        {
                max = MINF;
                max_Apo_codeword = MINF;
                
                for(c=0; c<C; c++) Apo_codeword(k,c)=MINF; /* some links are not possible */
                
                for(m=0; m<M; m++)
                {
                        Apo(k,m) = MINF; 
                        for(i=0; i<S; i++) 
                        {
                            abc = alpha(k,Ps(i,m)) + beta(k+1,i) + Apr_codeword(k, Lb(Ps(i,m),m));
                            Apo(k,m) = jacolog( Apo(k,m), abc );
                            Apo_codeword(k, Lb( Ps(i,m), m) ) = jacolog( Apo_codeword(k, Lb(Ps(i,m), m) ), abc + Apr_dataword(k,m));
                        } 
                        			
                        Apo(k,m) += Apr_dataword(k,m);
                        if ( max < Apo(k,m) ) max = Apo(k,m);

                }
		for(m=0; m<M; m++) Apo(k,m) -= max;
                
                for(c=0; c<C; c++) if(max_Apo_codeword < Apo_codeword(k,c)) max_Apo_codeword = Apo_codeword(k,c);
                for(c=0; c<C; c++) Apo_codeword(k,c) -= max_Apo_codeword;
        }

  return;
}


int CONV::assign_apr_dataword(vec llr_data_bit)
{  
  int j,m;
  
  LLR_to_SymProb(no_of_info_symbols, k, llr_data_bit, map->Apr_dataword);

  for(j=no_of_info_symbols;j<no_of_symbols;j++){
    map->Apr_dataword(j,0) = 0.0;
    for(m=1;m<trellis->nrof_datawords;m++)
      map->Apr_dataword(j,m) = MINF;
  }
  
  return 1;
}
int CONV::assign_apr_dataword(mat apr_dataword)
{ 
  int j, m;

  for(j=0;j<no_of_info_symbols;j++)
    for(m=0;m<trellis->nrof_datawords;m++)
      map->Apr_dataword(j,m) = apr_dataword(j,m);
   
  for(j=no_of_info_symbols;j<no_of_symbols;j++){
    map->Apr_dataword(j,0) = 0.0;
    for(m=1;m<trellis->nrof_datawords;m++)
      map->Apr_dataword(j,m) = MINF;
  }
  
  return 1;
}

int CONV::assign_apr_codeword(vec llr_coded_bit)
{  
  LLR_to_SymProb(no_of_symbols, n, llr_coded_bit, map->Apr_codeword);
  return 1;
}
int CONV::assign_apr_codeword(mat apr_codeword)
{  
  map->Apr_codeword = apr_codeword;
  return 1;
}

int CONV::get_apo_codeword_llr(vec &llr_coded_bit)
{  
  llr_coded_bit.set_size(no_of_symbols*n);
  
  SymProb_to_LLR(no_of_symbols, n, llr_coded_bit, map->Apo_codeword);
  return 1;
}
int CONV::get_apo_codeword(mat &apo_codeword)
{  
  apo_codeword = map->Apo_codeword;
  return 1;
}

int CONV::get_apo_dataword_llr(vec &llr_data_bit)
{  
  llr_data_bit.set_size(no_of_symbols*k);
  
  SymProb_to_LLR(no_of_symbols, k, llr_data_bit, map->Apo_dataword);
  return 1;
}
int CONV::get_apo_dataword(mat &apo_dataword)
{  
  apo_dataword = map->Apo_dataword;
  return 1;
}

void CONV::initialise(){
  srand(0);
  
  if(mode==0){
      map     = new struct MAP_decoder;
      map->Apr_codeword = mat(no_of_symbols, 2);
      map->Apo_dataword = mat(no_of_symbols, 2);
      
      return; // uncoded
  }

  switch(mode)
  {
    case SC:    
      allocate_memory();
      //initSC();
      break;
    case RSC:    
      allocate_memory();
      initRSC();
      break;
    case NSC:
      allocate_memory();
      initNSC();
      break;
    case RNSC:
      allocate_memory();
      initRNSC();
      break;
    default:
      s_error("no such convolutional code");
      break;   
  }
  
  printf("init done");
  //init_channel_interleaver(no_of_symbols);
}

void CONV::allocate_memory()
{
    int S, M, C, N;
    
    S = 1<<L;
    M = 1<<k;
    C = 1<<n;    
    N = no_of_symbols;
    
    //-------trellises
    trellis = new (struct TRELLIS_parameters);
    
    trellis->L = L;
    trellis->k = k;
    trellis->n = n;

    trellis->nrof_states   = S;
    trellis->nrof_branches = M;
    trellis->nrof_datawords= M;
    trellis->nrof_codewords= C;
    
    trellis->Lb = imat(S, M);
    trellis->Ns = imat(S, M);
    trellis->Ps = imat(S, M);
    trellis->MaxPs = imat(S, C);

    if(mode == RNSC){
      trellis->ptrPs = new struct trellis_state *[S]; 
      for(int i=0; i<S; i++){ 
	trellis->ptrPs[i] = new struct trellis_state[C];
	for(int j=0; j<C; j++){ 
	  trellis->ptrPs[i][j].word = NO_LINK;
	  trellis->ptrPs[i][j].infoword = NO_LINK;
	}  
      }
    }
    
    
    trellis->GenPoly.set_size(n, false);
    trellis->GenPolyMat = imat(k, n);
    
    //--------MAP decoder
    map     = new struct MAP_decoder;
    
    map->Apr_dataword = mat(N, M);
    map->Apo_dataword = mat(N, M);
    map->Extr_dataword= mat(N, M);
    
    map->Apr_codeword = mat(N, C);
    map->Apo_codeword = mat(N, C);
    map->Extr_codeword= mat(N, C);
    
    int j,m,c;
    for(j=0; j<N; j++){
      for(m=0; m<M; m++) map->Apr_dataword(j,m) = 0.0;
      for(c=0; c<C; c++) map->Apr_codeword(j,c) = 0.0; 
    }
    
}

void CONV::initRSC()
{
        int i, j, S, h, s, m, M, c, L;
	imat H;
        ivec D, K, cpoly;

	if(trellis->n!=trellis->k+1) s_error("Codes having n > k+1 is not implemented yet. Pls use n=k+1.");
	
        /*---------Initialise and memory allocation--------------------------*/   
	L = trellis->L; 
        H = imat(trellis->n, trellis->L+1);

        D.set_size(trellis->L+1,false);
        K.set_size(trellis->n, false);
        cpoly.set_size(trellis->n, false);
        
        M = 1<<trellis->k;
        S = 1<<trellis->L;
        /*---------END_Initialise and memory allocation----------------------*/
        
        // get generator polynomial depending on L and k
        get_RSC_poly();
	
        // make a copy of poly
        for(i=0; i<trellis->n; i++) cpoly[i] = trellis->GenPoly[i];

        /* compute explicitely poli coeff from octal rep */ 
        h=0; 
        for (i=0; i <= L; i++)
        {
                /* extract last dec once every three times */
                h=h%3;
                if ( h == 0 )
                {
                        for( m=0; m<trellis->n; m++ )
                        {       K[m] = cpoly[m]%10;     /* extract */
                            cpoly[m] = cpoly[m]/10;     /* shift */
                        }
                }
                h++;

                /* compute poli coeff */
                for( m=0; m<trellis->n; m++ )
                {       
		    H(m,i) = K[m]%2;       /* compute */
                    K[m] = K[m]/2;              /* shift */
                }
        }
        /* a check */
        if ( H(0,0) != 1 ) s_error( "RSC: the feedback poly is inacceptable" );

        /* compute tables */
        for ( s = 0; s < S; s++ )
        {
            /* fill vector D (state) */
                h = s;
                for(i = 0; i<trellis->L; i++)
                {
                        D[i] = h%2;     
                        h = h / 2;
                }

                /* compute next state and output */
                for( m=0; m<M; m++ )
                {
                        /* fill vector K (input bits) */
                        h=m;
                        for(i=1; i<trellis->n; i++)
                        {
                                K[i] = h%2;
                                h=h/2;
                        }

                        /* compute output bit (store in K[0] ) */
                        h = D[0];
                        for( i=1; i<trellis->k; i++)   h = ( h + K[i]*H(i,0) ) %2;
                        trellis->Lb(s,m) = 2*m + h;
                        K[0] = h;

                        /* compute new state */                 
                        c = 1; trellis->Ns(s,m) = 0;
                        for( j=0; j<trellis->L; j++)
                        {
                                /* bit from previous reg. */
                                if ( j < trellis->L-1 ) h=D[j+1];
                                else h=0;

                                /* input and feedback bits */
                                for( i=0; i<trellis->n; i++)   h = ( h + K[i]*H(i,j+1) ) %2;
                        
                                /* add to state */
                                trellis->Ns(s,m) += h*c;
                                c = c*2;
                        }
                }
        }
	
        /* compute previous state table */
        for(i=0;i<S;i++)
          for(j=0;j<M;j++)
            trellis->Ps( trellis->Ns(i,j) ,j)=i;


	/* compute previous state */
	int d;
	for(i=0; i<trellis->nrof_states; i++) for(d=0; d<trellis->nrof_codewords;d++)
	  trellis->MaxPs(i,d)=NO_LINK;                        
	for(i=0; i<trellis->nrof_states; i++) for(d=0; d<trellis->nrof_datawords;d++)
	{
	    trellis->MaxPs( trellis->Ns(i,d) , trellis->Lb(i,d) )
	      = i;
	}

}

void CONV::get_RSC_poly()
{
    int k=trellis->k;
    int L=trellis->L;

    /* generator polynimials */
    switch(k)
    {
        case 1:
            switch(L)
            {
                case 1:
                    trellis->GenPoly[0] = 1;  trellis->GenPoly[1] = 2;
                    break;            
                case 3:
                    trellis->GenPoly[0] = 13;  trellis->GenPoly[1] = 6;
                    break;
                case 4:
                    trellis->GenPoly[0] = 23;  trellis->GenPoly[1] = 6;
                    break;
                case 6:
                    trellis->GenPoly[0] = 117; trellis->GenPoly[1] = 26; 
                    break;
                case 7:
                    trellis->GenPoly[0] = 217; trellis->GenPoly[1] = 110;
                    break;
                case 8:
                    trellis->GenPoly[0] = 427; trellis->GenPoly[1] = 230;
                    break;
                case 9:
                    trellis->GenPoly[0] = 1017; trellis->GenPoly[1] = 120;
                    break;
                default:
                    s_error("no generator for such code yet: k=1, n=2");
                    break;                
            }
            break;
        case 2:
            switch(L)
            {
	        //case 1: //dummy
                //    trellis->GenPoly[0] = 1;  trellis->GenPoly[1] = 2;  trellis->GenPoly[2] = 0;
                //    break;
                case 3:
                    trellis->GenPoly[0] = 11;  trellis->GenPoly[1] = 2;  trellis->GenPoly[2] = 4;
                    break;
                case 4:
                    trellis->GenPoly[0] = 23;  trellis->GenPoly[1] = 2;  trellis->GenPoly[2] = 10;
                    break;
                case 6:
                    trellis->GenPoly[0] = 103; trellis->GenPoly[1] = 30; trellis->GenPoly[2] = 66;
                    break;
                case 7:
                    trellis->GenPoly[0] = 277; trellis->GenPoly[1] = 54; trellis->GenPoly[2] = 122;
                    break;
                case 8:
                    trellis->GenPoly[0] = 435; trellis->GenPoly[1] = 72; trellis->GenPoly[2] = 130;
                    break;
                default:
                    s_error("no generator for such code yet: k=2, n=3");
                    break;                
            }
            break;
        case 3:
            switch(L)
            {
                case 3:
                    trellis->GenPoly[0] = 11;  trellis->GenPoly[1] = 2;  trellis->GenPoly[2] = 4;  trellis->GenPoly[3] = 10;
                    break;
                case 4:
                    trellis->GenPoly[0] = 23;  trellis->GenPoly[1] = 2;  trellis->GenPoly[2] = 4;  trellis->GenPoly[3] = 10;
                    break;
                case 6:
                    trellis->GenPoly[0] = 101; trellis->GenPoly[1] = 16; trellis->GenPoly[2] = 64; trellis->GenPoly[3] = 0;
                    break;
                case 7:
                    trellis->GenPoly[0] = 203; trellis->GenPoly[1] = 14; trellis->GenPoly[2] = 42; trellis->GenPoly[3] = 10; 
                    break;
                default:
                    s_error("no generator for such code yet: k=3, n=4");
                    break;                
            }
            break;
        case 5:
            switch(L)
            {
                case 3:
                    trellis->GenPoly[0] = 11; trellis->GenPoly[1] = 2;  trellis->GenPoly[2] = 4;
                    trellis->GenPoly[3] = 0;  trellis->GenPoly[4] = 0;  trellis->GenPoly[5] = 0;
                    break;
                case 6:  /*this is derived from k=3,L=6 case above 16.3.00*/
                    trellis->GenPoly[0] = 101; trellis->GenPoly[1] = 16;  trellis->GenPoly[2] = 64; 
                    trellis->GenPoly[3] = 0;   trellis->GenPoly[4] = 0;   trellis->GenPoly[5] = 0;
                    break;    
                default:
                    s_error("no generator for such code yet: k=5, n=6");
                    break;
            }
            break;
        case 7:
            switch(L)
            {
                case 3:
                    trellis->GenPoly[0] = 11; trellis->GenPoly[1] = 2;  trellis->GenPoly[2] = 4;
                    trellis->GenPoly[3] = 0;  trellis->GenPoly[4] = 0;  trellis->GenPoly[5] = 0;
                    trellis->GenPoly[6] = 0;  trellis->GenPoly[7] = 0;
                    break;
                default:
                    s_error("no generator for such code yet: k=7, n=8");
                    break;                
            }
            break;
        default:
            s_error("no generator for such code yet: k>7");
            break;
    }
    
}


void CONV::initNSC()
{
  int i, s, m, k, test, M;
  int state;     /* number of state */
  int data;      /* number of data word */
  int d, c, sum;
  ivec lM;       /* actual code memory per data bit */
  imat cwb;     /* codeword bits of individual state */
  imat *cw;     /* codeword of individual state */
  ivec com_s;    /* state corresponding to a component */
  ivec com_d;    /* data corresponding to a component */

  int lD, lC, K;
  imat gen;  

  /*---------Initialise and memory allocation--------------------------*/
  K = trellis->L;    lD = trellis->k;    lC = trellis->n;

  gen = imat(lD, lC);
  
  lM.set_size(lD, false); 
  cwb= imat(lD, lC);

  com_s.set_size(lD);
  com_d.set_size(lD);

  /* get generator polynomial depending on L and k */
  get_NSC_poly(trellis->GenPolyMat);

      /* make a copy of poly */
  for(d=0; d<lD; d++) for(c=0; c<lC; c++) gen(d,c) = trellis->GenPolyMat(d,c);
  M = trellis->M;
  
  //cw = (imat *) calloc(sizeof(imat), lD);
  cw = new imat[lD];
  for(i=0; i<lD; i++) cw[i].set_size((int)pow(2,M), 2, false);

  /*---------END_Initialise and memory allocation----------------------*/

  
  /*get the equivalent of decimal generator, for shift function at a later state*/
  for(i=0;i<lD;i++) for(c=0;c<lC;c++)
    gen(i,c) = octal_to_decimal(gen(i,c));
  
  /*get the exact shift register state for each input bit and the exact generator polynomial*/
  for(i=0; i<lD; i++) for(m=0,test=0,lM[i]=M; m<M; m++){
    for(c=0; c<lC; c++)
      test += gen(i,c)&1; 
      
    if(test==0){ lM[i]--; for(c=0;c<lC;c++) gen(i,c)=gen(i,c)>>1; }
    else m=M;
  }

  for(m=0,sum=0; m<lD; m++) {
#ifdef debug_conv
    printf("\nlM[%d] = %d", m, lM[m]); 
#endif
    sum+=lM[m];
  } 
  if(K!=sum) s_error("(initCC) Conflict state number K with lM.");

#ifdef debug_conv
  printf("\ngen in decimal:");
  for(i=0;i<lD;i++) {printf("\ngen[%d] = ",i); for(c=0;c<lC;c++) printf("%d ",gen(i,c));}printf("\n");
#endif
  
  /* generate label table for individual section */
  for(i=0;i<lD;i++){
    state = (int)pow(2,lM[i]);
    for(d=0;d<2;d++){
      for(s=0;s<state;s++){
	
	for(c=0;c<lC;c++){
	  cwb(i,c) = (d&1)*((gen(i,c)>>lM[i])&1);

	  for(k=lM[i]-1; k>=0; k--){
	    cwb(i,c) += ((s>>k)&1)*((gen(i,c)>>k)&1);
	  }

	  if(cwb(i,c)%2) cwb(i,c)=1;  /*simple modulo two adder*/
	  else cwb(i,c)=0;

	  //printf("\nd=%d, s=%d,\tcwb(%d,%d)=%d",d,s,i,c,cwb(i,c));
	}

	cw[i](s,d) = cwb(i,0);
	for(c=1;c<lC;c++)
	  cw[i](s,d) += cwb(i,c)*(int)pow(2,c) ;
      }    
    }
  }

  /* print sub table */
  /*
  for(i=0;i<lD;i++){
    state = pow(2,lM[i]);
    printf("\n-------------------");
    for(d=0;d<2;d++){
      printf("\n");
      for(s=0;s<state;s++){
	printf("cw[%d](%d,%d)=%d ",i,s,d,cw[i](s,d));
      }
    }
  }printf("\n");
  */

  /* generate label table */
  state = (int) pow(2,K);   data = (int) pow(2,lD);
  for(s=0;s<state;s++){
    for(d=0;d<data;d++){

      for(i=0,sum=0,k=0;i<lD;i++){
	com_d[i] = (d>>i) & 1;

	if(i==0) com_s[i] = s & ( (int)pow(2,lM[i]) - 1 );
	else{
	  k += lM[i-1];
	  com_s[i] = ( s>>k ) & ( (int)pow(2,lM[i]) - 1 );
	}
	//printf("\nd=%d -> com_d[%d]=%d\t\ts=%d -> com_s[%d]=%d",d,i,com_d[i],s,i,com_s[i]);
	sum ^= cw[i]( com_s[i] , com_d[i] );
      }
      
      trellis->Lb(s,d) = sum;
      //printf("\nd=%d s=%d, trellis->Lb=%d",d,s,trellis->Lb(s,d));
    }
  }
  
  /* compute next state */
  for(s=0;s<state;s++){
    for(d=0;d<data;d++){
      
      for(i=0,sum=0,k=0;i<lD;i++){

	if(i==0) com_s[i] = s & ( (int)pow(2,lM[i]) - 1 );
	else{
	  k += lM[i-1];
	  com_s[i] = ( s>>k ) & ( (int)pow(2,lM[i]) - 1 );
	}

	com_s[i]>>= 1;
	com_s[i] ^= (int)pow(2,lM[i]-1)*( (d>>i)&1 );

	//printf("\nd=%d s=%d -> com_s[%d]=%d",d,s,i,com_s[i]);	
	sum ^= com_s[i]<<k;
      }
      
      trellis->Ns(s,d)=sum;
      //printf("\nd=%d s=%d, Ns=%d",d,s,trellis->Ns(s,d));
    }
  }

  
  /* compute previous state */
  for(i=0; i<trellis->nrof_states; i++) for(d=0; d<trellis->nrof_codewords;d++)
      trellis->MaxPs(i,d)=NO_LINK;                        
  for(i=0; i<trellis->nrof_states; i++) for(d=0; d<trellis->nrof_datawords;d++)
  {
      trellis->MaxPs( trellis->Ns(i,d) , trellis->Lb(i,d) )
          = i;
  }

  delete [] cw;
}


int calc_state_transition(const int instate, const int input, ivec &parity, 
			  ivec &gen_pol_rev, int K, int n)
{
     int i, j, parity_temp, parity_bit;
     int m = K-1;
     int in = 0;
     int temp = (gen_pol_rev(0) & (instate<<1));
 
     for (i=0; i<K; i++) {
         in = (temp & 1) ^ in;
         temp = temp >> 1;
     }
     in = in ^ input;
 
     parity.set_size(n-1,false);
     for (j=0; j<(n-1); j++) {
       parity_temp = ((instate<<1) + in) & gen_pol_rev(j+1);
       parity_bit = 0;
       for (i=0; i<K; i++) {
         parity_bit = (parity_temp & 1) ^ parity_bit;
         parity_temp = parity_temp >> 1;
       }
       parity(j) = parity_bit;
     }
     return in + ((instate << 1) & ((1<<m)-1));
}




void CONV::initRNSC()
{
  int i,j;
  int K = trellis->L+1;
  
  if(trellis->k != 1) it_error("RNSC for k==1 only"); 

  // get generator polynomial depending on L and k
  get_NSC_poly(trellis->GenPolyMat);
  
  ivec gen(trellis->n);
  for(i=0;i<trellis->n;i++) gen(i) = trellis->GenPolyMat(0,i);

  Rec_Syst_Conv_Code rnsc;
  rnsc.set_generator_polynomials (gen, K);
  
  ///////////////
  int Nstates = (int) pow(2,trellis->L);  
  int s0, s1, s_prim, temp;
  ivec p0, p1;

  for (s_prim=0; s_prim<Nstates; s_prim++) {
     s0 = calc_state_transition(s_prim,0,p0, gen,K,trellis->n);
     trellis->Ns(s_prim,0) = s0;
     //trellis->Ps(s0,0) = s_prim;
     for (j=temp=0; j<(trellis->n-1); j++) temp += (1<<j) * p0(j);     
     trellis->Lb(s_prim,0) = temp;
     trellis->Lb(s0,0)     = temp;
     
     s1 = calc_state_transition(s_prim,1,p1, gen,K,trellis->n);
     trellis->Ns(s_prim,1) = s1;
     //trellis->Ps(s1,1) = s_prim;
     for (j=temp=0; j<(trellis->n-1); j++) temp += (1<<j) * p1(j);
     trellis->Lb(s_prim,1) = temp + (1<<(trellis->n-1));
     trellis->Lb(s1,1)     = temp + (1<<(trellis->n-1));
  }

  // compute previous state 
  int Ns, Cw, m;
  for(i=0; i<Nstates; i++) for(j=0; j<trellis->nrof_datawords;j++){

    Ns = trellis->Ns(i,j);
    Cw = trellis->Lb(i,j);
    
    for(m=0;m<trellis->nrof_codewords;m++) if(trellis->ptrPs[Ns][m].word == NO_LINK) break;   
    
    if(trellis->ptrPs[Ns][m].word != NO_LINK){
        printf("i%d j%d m%d Ns%d Cw%d\n", i,j,m, Ns, Cw);
        printf("Not enough nrof_ps_branches (%d). Pls edit the code to make it bigger!\n", trellis->nrof_codewords);
        exit(-1);
    }
    
    trellis->ptrPs[Ns][m].infoword = j;
    trellis->ptrPs[Ns][m].word = Cw;
    trellis->ptrPs[Ns][m].index = i;
  }

  /*
  for (i=0; i<Nstates; i++){
    printf("Lb(%d,0) -> %d  ", i, trellis->Lb(i,0));
    printf("Lb(%d,1) -> %d\n", i, trellis->Lb(i,1));
  }
  for (i=0; i<Nstates; i++){
    printf("Ns(%d,0) -> %d  ", i, trellis->Ns(i,0));
    printf("Ns(%d,1) -> %d\n", i, trellis->Ns(i,1));
  }
  for (i=0; i<Nstates; i++){
    for(j=0; j<trellis->nrof_codewords;j++)
      printf("Ps(%d,%d)->%d ", i, j, trellis->ptrPs[i][j].index);
    printf("\n");
  }
  */

}

void CONV::get_NSC_poly(imat &GenPolyMat)
{ 
    int k = trellis->k;
    int n = trellis->n;
    int L = trellis->L;

    /* generator polynimials */
    switch(k)
    {
        case 1:
	// from page 492-493 Digital Communications 4th edition: John G. Proakis
	    if(n==2){ // R=1/2
	      switch(L)
	      {
		case 2:
		    trellis->M = 2;
		    GenPolyMat(0,0) = 7;   GenPolyMat(0,1) = 5;
		    break;
		case 3:
		    trellis->M = 3;		    
		    if(mode==RNSC){ // used for Turbo code
		      GenPolyMat(0,0) = 13;   GenPolyMat(0,1) = 15; 
		    }
		    else{ // highest free distance
		      GenPolyMat(0,0) = 15;   GenPolyMat(0,1) = 17;
		    }
		    break;
		case 4:
		    trellis->M = 4;
		    GenPolyMat(0,0) = 23;   GenPolyMat(0,1) = 35;
		    break;
		case 5:
		    trellis->M = 5;
		    GenPolyMat(0,0) = 53;   GenPolyMat(0,1) = 75;
		    break;
		case 6:
		    trellis->M = 6;
		    GenPolyMat(0,0) = 133;  GenPolyMat(0,1) = 171;
		    break;
		case 7:
		    trellis->M = 7;
		    GenPolyMat(0,0) = 247;  GenPolyMat(0,1) = 371;
		    break;
		case 8:
		    trellis->M = 8;
		    GenPolyMat(0,0) = 561;  GenPolyMat(0,1) = 753;
		    break;
		case 9:
		    trellis->M = 9;
		    GenPolyMat(0,0) = 1167;  GenPolyMat(0,1) = 1545;
		    break;
		case 10:
		    trellis->M = 10;
		    GenPolyMat(0,0) = 2335;  GenPolyMat(0,1) = 3661;
		    break;
		case 11:
		    trellis->M = 11;
		    GenPolyMat(0,0) = 4335;  GenPolyMat(0,1) = 5723;
		    break;

		default:
		    s_error("no generator for such code yet: k=1, n=2");
		break;
	      }
	    }
	    else if(n==3){ // R=1/3
	      switch(L)
	      {
		case 3:
		    trellis->M = 3;
		    GenPolyMat(0,0) = 13;   GenPolyMat(0,1) = 15;   GenPolyMat(0,2) = 17;
		    break;
		case 4:
		    trellis->M = 4;
		    GenPolyMat(0,0) = 25;   GenPolyMat(0,1) = 33;   GenPolyMat(0,2) = 37;
		    break;
		case 5:
		    trellis->M=5;
		    GenPolyMat(0,0) = 47;   GenPolyMat(0,1) = 53;   GenPolyMat(0,2) = 75;
		    break;
		case 6:
		    trellis->M = 6;
		    GenPolyMat(0,0) = 133;  GenPolyMat(0,1) = 145;  GenPolyMat(0,2) = 175;
		    break;
		default:
		    s_error("no generator for such code yet: k=1, n=3");
		break;
	      }
	    }
	    else if(n==4){ // R=1/4
	      switch(L)
	      {
		case 3:
		    trellis->M = 3;
		    GenPolyMat(0,0) = 13;   GenPolyMat(0,1) = 15;   GenPolyMat(0,2) = 15;   GenPolyMat(0,3) = 17;
		    break;
		case 4:
		    trellis->M = 4;
		    GenPolyMat(0,0) = 25;   GenPolyMat(0,1) = 27;   GenPolyMat(0,2) = 33;   GenPolyMat(0,3) = 37;
		    break;
		case 5:
		    trellis->M=5;
		    GenPolyMat(0,0) = 53;   GenPolyMat(0,1) = 67;   GenPolyMat(0,2) = 71;   GenPolyMat(0,3) = 75;
		    break;
		case 6:
		    trellis->M = 6;
		    GenPolyMat(0,0) = 135;  GenPolyMat(0,1) = 135;  GenPolyMat(0,2) = 147;  GenPolyMat(0,3) = 163;
		    break;
		default:
		    s_error("no generator for such code yet: k=1, n=4");
		break;
	      }
	    }
	    else{
	      printf("code for k=%d n=%d\n", k, n);
	      s_error("no generator for such code yet");
	    }

            break;
        case 2:
	// from page 331 Error Control Coding: Shu Lin and Daniel J. Costello
	    if(n!=k+1) s_error("no generator for such code yet: k=2, n!=3");
            switch(L)
            {
                case 3:
                    trellis->M = 2;
                    GenPolyMat(0,0) = 4;  GenPolyMat(0,1) = 2;  GenPolyMat(0,2) = 6;
                    GenPolyMat(1,0) = 1;  GenPolyMat(1,1) = 4;  GenPolyMat(1,2) = 7;
                    break;
                case 4:
                    trellis->M = 2;
                    GenPolyMat(0,0) = 7;  GenPolyMat(0,1) = 1;  GenPolyMat(0,2) = 4;
                    GenPolyMat(1,0) = 2;  GenPolyMat(1,1) = 5;  GenPolyMat(1,2) = 7;
                    break;
                case 5:
                    trellis->M = 3;
                    GenPolyMat(0,0) = 14;  GenPolyMat(0,1) = 06;  GenPolyMat(0,2) = 16;
                    GenPolyMat(1,0) = 03;  GenPolyMat(1,1) = 10;  GenPolyMat(1,2) = 17;
                    break;
	        case 6: // see my note book for these values
                    trellis->M = 3;
                    GenPolyMat(0,0) = 15;  GenPolyMat(0,1) = 06;  GenPolyMat(0,2) = 15;
                    GenPolyMat(1,0) = 06;  GenPolyMat(1,1) = 15;  GenPolyMat(1,2) = 17;
                    break;
                default:
                    s_error("no generator for such code yet: k=2, n=3");
                    break;                
            }
            break;
        case 3:
	// from page 331 Error Control Coding: Shu Lin and Daniel J. Costello
	    if(n!=k+1) s_error("no generator for such code yet: k=3, n!=4");
            switch(L)
            {
                case 3:
                    trellis->M = 2;
                    GenPolyMat(0,0) = 4;  GenPolyMat(0,1) = 4;  GenPolyMat(0,2) = 4;  GenPolyMat(0,3) = 4;
                    GenPolyMat(1,0) = 0;  GenPolyMat(1,1) = 6;  GenPolyMat(1,2) = 2;  GenPolyMat(1,3) = 4;
                    GenPolyMat(2,0) = 0;  GenPolyMat(2,1) = 2;  GenPolyMat(2,2) = 5;  GenPolyMat(2,3) = 5;
                    break;
		case 5:
		    trellis->M = 2;
		    GenPolyMat(0,0) = 6;  GenPolyMat(0,1) = 2;  GenPolyMat(0,2) = 2;  GenPolyMat(0,3) = 6;
		    GenPolyMat(1,0) = 1;  GenPolyMat(1,1) = 6;  GenPolyMat(1,2) = 0;  GenPolyMat(1,3) = 7;
		    GenPolyMat(2,0) = 0;  GenPolyMat(2,1) = 2;  GenPolyMat(2,2) = 5;  GenPolyMat(2,3) = 5;
		    break;
                case 6:
                    trellis->M = 2;
                    GenPolyMat(0,0) = 6;  GenPolyMat(0,1) = 1;  GenPolyMat(0,2) = 0;  GenPolyMat(0,3) = 7;
                    GenPolyMat(1,0) = 3;  GenPolyMat(1,1) = 4;  GenPolyMat(1,2) = 1;  GenPolyMat(1,3) = 6;
                    GenPolyMat(2,0) = 2;  GenPolyMat(2,1) = 3;  GenPolyMat(2,2) = 7;  GenPolyMat(2,3) = 4;
                    break;
                default:
                    s_error("no generator for such code yet: k=3, n=4");
                    break;                
            }
            break;
        default:
            s_error("no generator for such code yet: k>3");
            break;
    }
}




































































/*-----------symbol and bits conversion------------------*/
void CONV::symbol_to_bits(int word_length, int block_length, ivec symbols, imat &bits_block){
  int i, j;
  
  for(j=0; j<block_length; j++)
      for(i=0; i<word_length; i++)
          bits_block(i,j) = (symbols[j]>>i) & 1 ;
}

void CONV::bits_to_symbol(int word_length, int block_length, imat bits_block, ivec &symbols){
  int i, j;
  
  for(j=0; j<block_length; j++)
      for(i=0,symbols[j]=0; i<word_length; i++)
          symbols[j] += bits_block(i,j) * (1 << i);
}

ivec CONV::symbol_to_bits_seq(int word_length, int block_length, ivec symbols){
  int i, j, k;
  int N=block_length*word_length;
  ivec bits_seq;
  bits_seq.set_length(N, false);

  if(N!=bits_seq.length()){ 
    printf("%d %d ",N, bits_seq.length());
    s_error("symbol_to_bits_seq: check bits_seq.length()");
  }
  
  for(k=j=0; j<block_length; j++)
      for(i=0; i<word_length; i++)
          bits_seq[k++] = (symbols[j]>>i) & 1 ;

  return bits_seq;
}

ivec CONV::bits_seq_to_symbol(int word_length, int block_length, ivec bits_seq){
  int i, j, k;
  int N=block_length*word_length;
  ivec symbols;
  symbols.set_length(block_length, false);

  if(N!=bits_seq.length()){
    printf("%d %d ",N, bits_seq.length());
    s_error("bits_seq_to_symbol: check bits_seq.length()");
  }

  for(j=k=0; j<block_length; j++)
    for(i=0,symbols[j]=0; i<word_length; i++){
          symbols[j] += bits_seq[k++] * (1 << i);
	  //printf("%d(%d) ",bits_seq[k-1], symbols[j]);
    }
  
  return symbols;
}




/*---------utilities-------------*/
void CONV::s_error ( const char * message )
{
    fprintf (stderr, "Fatal: %s\n", message);
    fflush ( stderr );
    exit (-1);
}

double CONV::correct_rounding_error(double distance){
  // to the precision of 6 decimal places
  double dist;    
  dist  = floor(distance * 1.0e8);    
  dist /= 1.0e8;
  
  dist  = ceil(dist * 1.0e6);    
  dist /= 1.0e6;
  
  return (dist);
}

int CONV::octal_to_decimal(int gen){
  
  if(gen<10) return gen;
  else if(gen<100){
    return(8*(gen/10) + gen%10);
  }
  else if(gen<1000){    
    return(64*(gen/100) + 8*( (gen%100)/10 ) + gen%10);
  }
  else if(gen<10000){
    return(512*(gen/1000) + 64*( (gen%1000)/100 ) + 8*( (gen%100)/10 ) + gen%10);
  }
  else s_error("Octal to decimal not complete for octal number>10000.");

  return 0;
}


/*-----LLR <-> Prob----------------*/
void CONV::SymProb_to_LLR ( int N, int bps, vec &LLR, mat SymPr )
{
  mat BitPr;
  
  BitPr = mat(N*bps, 2);
  
  Pr_to_BitPr_log(N, bps, SymPr, BitPr, 0);
  Prob_to_LLR (N, bps, LLR, BitPr);
}

void CONV::LLR_to_SymProb ( int N, int bps, vec LLR, mat & SymPr )
{
  mat BitPr;
  
  BitPr = mat(N*bps, 2);
  
  LLR_to_Prob(N, bps, LLR, BitPr);
  BitPr_to_Pr_log(N, bps, SymPr, BitPr, 0);
}

void CONV::Prob_to_LLR ( int N, int bps, vec &LLR, mat BitPr )
{
  int k, i;
  for ( k = 0; k < N; k++ )
    for ( i = 0; i < bps; i++ )
        LLR[k*bps+i] = BitPr(k*bps+i, 1) - BitPr(k*bps+i, 0); 
}

void CONV::LLR_to_Prob ( int N, int bps, vec LLR, mat &BitPr )
{
  int k, i;
  
  if((N*bps)!=LLR.length()) it_error("N*bps!=LLR.length() in LLR_to_Prob");

  for ( k = 0; k < N; k++ )
  {
    for ( i = 0; i < bps; i++ )
    {
      BitPr(k*bps+i, 0) = - jacolog_1 ( 0.0, LLR[k * bps + i] );
      BitPr(k*bps+i, 1) = - jacolog_1 ( 0.0, - LLR[k * bps + i] );
    }
  }

}

void CONV::Pr_to_BitPr_log(int N, int bps, mat Pr, mat &BitPr, int frame_index){
    int k, m, i, M;
    double max;

    M = (int)pow(2, bps);
    
    for(k=0;k<N;k++){
        for(i=0;i<bps;i++){
	  
            max=MINF; BitPr(i+k*bps+frame_index*N, 0) = BitPr(i+k*bps+frame_index*N, 1) = MINF;
            for(m=0;m<M;m++)
	      BitPr(i+k*bps+frame_index*N, (m>>i)&1 ) = jacolog( BitPr(i+k*bps+frame_index*N, (m>>i)&1 ) , Pr(k,m));
            
            if( BitPr(i+k*bps+frame_index*N, 0) < BitPr(i+k*bps+frame_index*N, 1) ) max = BitPr(i+k*bps+frame_index*N,1);
            else max = BitPr(i+k*bps+frame_index*N,0);
            
            for(m=0;m<2;m++) BitPr(i+k*bps+frame_index*N,m) -= max;
        }
    }
}

void CONV::BitPr_to_Pr_log(int N, int bps, mat &Pr, mat BitPr, int frame_index){
    int k, m, i, M;
    double max;

    M = (int)pow(2, bps);
  
    for(k=0;k<N;k++){
        max = MINF;
        for(m=0;m<M;m++){
            Pr(k,m) = 0; /*log 1 = 0*/
            for(i=0;i<bps;i++){
                Pr(k,m) += BitPr(i+k*bps+frame_index*N, (m>>i)&1 );
            }

            if(max < Pr(k,m)) max = Pr(k,m);
        }
        for(m=0;m<M;m++) Pr(k,m) -= max;
        for(m=0;m<M;m++){ if(Pr(k,m) < MINF) Pr(k,m)=MINF; }
    }
}

/*------ jacobian logarithm ------------------ */
double CONV::jacolog_1( double x, double y)
{
    double r;
    
    if(x>y) r = x + log ( 1 + exp( y - x ) );
    else    r = y + log ( 1 + exp( x - y ) );

    return r;
}

double CONV::jacolog_2( double x, double y)
{       
        double r;
        double diff;

        if(x>y){ r = x; diff=x-y; }
        else   { r = y; diff=y-x; }  
              
        if(diff > 3.7 )      r += 0.00;
        else if(diff > 2.25) r += 0.05;
        else if(diff > 1.5 ) r += 0.15;
        else if(diff > 1.05) r += 0.25;
        else if(diff > 0.7 ) r += 0.35;
        else if(diff > 0.43) r += 0.45;
        else if(diff > 0.2 ) r += 0.55;
        else                 r += 0.65;
        
        return r;
}

double CONV::jacolog_3( double x, double y)
{
    double r;
    
    if(x>y) r = x;
    else    r = y;
    
    return r;
}

double CONV::jacolog(double x, double y)
{
  if(decoder_type==C_EXACT_LOG_MAP) return jacolog_1(x, y);
  else if(decoder_type==C_APPROX_LOG_MAP) return jacolog_2(x, y);
  else if(decoder_type==C_MAX_LOG_MAP) return jacolog_3(x, y);
  else s_error("unknown decoder type");
}
