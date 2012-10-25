/*! 
  \file
  \brief Definition of Convolutional codes (CONV) class
  \author Michael Ng

  1.01

  2004/02/24
*/

#ifndef _CONV_H_
#define _CONV_H_

#include <itpp/itbase.h>
#include "conv.h"

/*! \defgroup FEC Forward Error Correction
 */

using namespace itpp;
using namespace std;

//#define debug_conv

//! define OFF
#define OFF      0
//! define ON
#define ON       1

//! define Non-Systematic Convolutional code
#define NSC    1
//! define Systematic Convolutional code
#define SC     2
//! define Recursive Systematic Convolutional code (from TCM)
#define RSC    3 
//! define Recursive Systematic Convolutional code from NSC
#define RNSC    4

//! define Exact Logarithmic (Log) Maximum A-posteriori Probability (MAP) decoding
#define C_EXACT_LOG_MAP  1
//! define Approximate Log MAP decoding
#define C_APPROX_LOG_MAP 2
//! define Maximum Log MAP decoding
#define C_MAX_LOG_MAP    3

//! define Minimum Log probability (-infinity)
#define MINF   -100000
//! define No Linkage for a particular trellis transition
#define NO_LINK 100000

//! struct of TRELLIS
struct TRELLIS_parameters
{ 
  int k;          //!< length of dataword
  int n;          //!< length of codeword
  int L;          //!< code memory length
  int M;          //!< code memory length for individual databit for NSC
  
  int nrof_states;   //!< number of possible trellis states 
  int nrof_branches; //!< number of possible output trellis branches : 2^k
  int nrof_datawords;//!< number of possible datawords : 2^k
  int nrof_codewords;//!< number of possible codewords : 2^n
  
  ivec GenPoly;    //!< generator polynomials in vec form
  imat GenPolyMat; //!< generator polynomiala in matrix form

  imat Lb;        //!< codeword (trellis transition Label) table
  imat Ns;        //!< Table of next state 
  imat Ps;        //!< Table of previous state
  imat MaxPs;     //!< Table of all possible previous state
  struct trellis_state **ptrPs; 
};

//! struct of MAP decoder
struct MAP_decoder
{
  mat Apr_dataword;    //!< a-priori probability of dataword
  mat Apo_dataword;    //!< a-posteriori probability of dataword
  mat Extr_dataword;   //!< extrinsic probability of dataword

  mat Apr_codeword;    //!< a-priori probability of codeword
  mat Apo_codeword;    //!< a-posteriori probability of codeword
  mat Extr_codeword;   //!< extrinsic probability of codeword
};

/*! 
  \ingroup FEC
  \brief CONV Class
  
  This class provides four convolutional schemes: 
  - \ref NSC, \ref SC, \ref RSC and \ref RNSC
    
  Decoding in logarithmic domain is performed: \a decoder_type
  available are \ref C_EXACT_LOG_MAP, \ref C_APPROX_LOG_MAP, \ref C_MAX_LOG_MAP.
  
*/
class CONV
{
 public:
  //! Class constructor
  CONV(){}
  //! Destructor
  ~CONV();

  /*! 
    \brief Setup parameters for the CONV encoder/decoder and also calculate the bit and symbol block lengths based on a control file
    
    \param fname the control file name.

    \note after reading the parameters from a control file, set_parameters() will be called.
  */  
  void set_parameters_from_file(char *fname);
  
  /*! 
    \brief Setup parameters for the CONV encoder/decoder and also calculate the bit and symbol block lengths
    
    \param mode_in CM.mode: \ref NSC, \ref SC, \ref RSC, \ref RNSC
    \param k_in    dataword length
    \param n_in    codeword length
    \param L_in    code memory length
    \param Terminated_in        code termination: \ref OFF, \ref ON
    \param decoder_type_in      decoding type: \ref C_EXACT_LOG_MAP, \ref C_APPROX_LOG_MAP, \ref C_MAX_LOG_MAP
    \param no_of_symbols_in     number of coded symbols in a block 
    
    \note
    - \a k=k_in, \a n=n_in, \a L=L_in 
    - if \a Terminated is OFF: 
      -# \a no_of_info_symbols=no_of_symbols_in
    - if \a Terminated is ON: 
      -# \a no_of_info_symbols=no_of_symbols_in-L
    - \a no_of_info_bits=no_of_info_symbols*k 
    - \a no_of_coded_bits=no_of_symbols*n
    - \a no_of_tail_bits=(no_of_symbols-no_of_info_symbols)*k
  */  
  void set_parameters(const int mode_in, const int k_in, const int n_in, const int L_in, const int Terminated_in, 
		      const int decoder_type_in, const int no_of_symbols_in);
  
  /*!
    \brief Initialise CONV: allocate memory and initialise coder
  */
  void initialise();
  
  /*! 
    \brief CONV encoder: encode
    
    \param b_Input_bits     input data bits (binary vector) 
    \param i_Output_symbols output coded symbol (integer vector)
    
    \note
    - length of \a b_Input_bits=no_of_info_symbols*k
    - length of \a i_Output_symbols=no_of_symbols
  */ 
  int encode ( bvec b_Input_bits, ivec &i_Output_symbols);
  
  /*! 
    \brief CONV encoder: encode
    
    \param b_Input_bits     input data bits (binary vector) 
    \param b_Output_bits    output coded bits (binary vector)
    
    \note
    - length of \a b_Input_bits=no_of_info_symbols*k
    - length of \a b_Output_bits=no_of_symbols*n
  */ 
  int encode_bits ( bvec b_Input_bits, bvec &b_Output_bits);
  
  /*! 
    \brief CONV encoder: encode
    
    \param b_Input_bits     input data bits (binary vector) 
    
    \note
    - length of \a b_Input_bits=no_of_info_symbols*k
    - return i_output_coded_bits (integer vector)
    - length of \a i_output_coded_bits=no_of_symbols*n
  */ 
  ivec encode_bits ( bvec b_Input_bits);
  
  /*! 
    \brief CONV decoder: decode
    
    \param i_decoded_bits decoded bits
    \param i_decoded_symbols decoded symbols
    
    \note 
    - data symbol a-posteriori probability updated and stored in 
      -# map->Apo_dataword
    - coded symbol a-posteriori probability updated and stored in 
      -# map->Apo_codeword
  */   
  int decode( ivec &i_decoded_bits, ivec &decoded_symbols);

  /*! 
    \brief Logarithmic SISO(MAP) decoder for any coding scheme
    
    \param N             the block length
    \param D             the number of datawords (trellis branches)
    \param C             the number of codewords
    \param S             the number of trellis states
    \param MaxPs         Previous state table, consider ALL codewords
    \param Ns            Next state table
    \param Lb            Codeword (branch label)
    \param Apr_codeword  a-priori  probability of codeword
    \param Apr_dataword  a-priori  probability of dataword
    \param Apo_codedword a-posteriori probability of codeword
    \param Apo_dataword  a-posteriori probability of dataword
    \param Terminated    Code termination \a Terminated
    \param frame_index   if channel interleaver is \e J times longer than the coded block, \f$ frame\_index \in\{0\ldots J-1\} \f$
  */
  void SISO_dec_log(int N, int D, int C, int S, 
		    imat MaxPs, imat Ns, imat Lb, 
		    mat Apr_codeword, mat Apr_dataword,
		    mat &Apo_codeword, mat &Apo_dataword, 
		    int Terminated, int frame_index);
  
  void SISO_dec_log_RNSC(int N, int D, int C, int S, 
		    struct trellis_state **ptrPs, imat Ns, imat Lb, 
		    mat Apr_codeword, mat Apr_dataword,
		    mat &Apo_codeword, mat &Apo_dataword, 
		    int Terminated, int frame_index);
  
  /*! 
    \brief Logarithmic MAP decoder for RSC(TCM) faster
    
    \param N   the block length
    \param M   the number of trellis branches (dataword)
    \param C   the number of codewords
    \param S   the number of trellis states
    \param Ps  Previous state table, consider ONLY datawords
    \param Ns  Next state table
    \param Lb  Codeword (branch label)
    \param Apr_dataword a-priori probability of dataword
    \param Apr_codeword a-priori probability of codeword
    \param Apo_dataword a-posteriori probability of dataword
    \param Apo_codeword a-posteriori probability of codeword
  */
  void SISO_dec_log_RSC(int N, int M, int C, int S, 
			imat Ps, imat Ns, imat Lb, 
			mat Apr_dataword, mat Apr_codeword, 
			mat &Apo_dataword, mat &Apo_codeword);
  
  int assign_apr_dataword (vec apr_dataword);
  int assign_apr_dataword (mat apr_dataword);

  int assign_apr_codeword (vec llr_coded_bit);
  int assign_apr_codeword (mat apr_codeword);

  int get_apo_dataword_llr(vec &llr_data_bit);
  int get_apo_dataword(mat &apo_dataword);

  int get_apo_codeword_llr(vec &llr_coded_bit);
  int get_apo_codeword(mat &apo_codeword);
  
  /*! 
    \brief print parameters to a file pointer.

    printout the parameters of coder
  */
  void print_parameters(FILE *fp);  

  /*! 
    \brief print coding tables to a file pointer.
    
    printout the tables of Codeword NextState & PreviousState.
  */
  void print_coding_tables(FILE *fp);  
  
  //! return the result filename.
  char* get_result_filename(){return result_filename;}

  //! return the output coded symbol vector length.
  int get_sym_length(){return no_of_symbols;}
  //! return the input data symbol vector length
  int get_info_sym_length(){return no_of_info_symbols;}
  //! return the output coded bit vector length
  int get_coded_bit_length(){return no_of_coded_bits;}
  //! return the input data bit vector length
  int get_info_bit_length(){return no_of_info_bits;}
  //! return the input data bit vector length
  int get_termination_bit_length(){return no_of_tail_bits;}
  
  //! return the coder mode 
  int get_coder_mode(){return mode;}
  
  //! return the actual coding rate 
  double get_rate(){return (double)no_of_info_bits/(double)no_of_coded_bits; }

  //! return the number of bits per dataword
  int get_k(){return k; }

  //! return the number of bits per codeword
  int get_n(){return n; }
  
  //! return the number of code memory
  int get_memory_length(){return L; }

  //! Jacobian logarithmic summation: based on the \a decoder_type obtions.
  double jacolog(double x, double y);
  //! Jacobian logarithmic summation: \a decoder_type= \ref C_EXACT_LOG_MAP
  double jacolog_1( double x, double y);
  //! Jacobian logarithmic summation: \a decoder_type= \ref C_APPROX_LOG_MAP
  double jacolog_2( double x, double y);
  //! Jacobian logarithmic summation: \a decoder_type= \ref C_MAX_LOG_MAP
  double jacolog_3( double x, double y);

  void symbol_to_bits(int word_length, int block_length, ivec symbols, imat &bits_block);
  void bits_to_symbol(int word_length, int block_length, imat bits_block, ivec &symbols);
  ivec symbol_to_bits_seq(int word_length, int block_length, ivec symbols);
  ivec bits_seq_to_symbol(int word_length, int block_length, ivec bits_seq);

  void SymProb_to_LLR ( int N, int bps, vec &LLR, mat SymPr );
  void LLR_to_SymProb ( int N, int bps, vec LLR, mat &SymPr );
  void Prob_to_LLR ( int N, int bps, vec &LLR, mat BitPr );
  void LLR_to_Prob ( int N, int bps, vec LLR, mat &BitPr );
  void Pr_to_BitPr_log(int N, int bps, mat Pr, mat &BitPr, int frame_index);
  void BitPr_to_Pr_log(int N, int bps, mat &Pr, mat BitPr, int frame_index);
 
 protected:
  struct TRELLIS_parameters  *trellis;     //!< Code Trellis 
  struct MAP_decoder *map;                 //!< MAP decoder 
  
  //! allocate memory
  void allocate_memory();
  //! initialise SC coder
  void initSC();
  //! initialise RSC coder
  void initRSC();
  //! initialise NSC coder
  void initNSC();
  //! initialise RNSC coder
  void initRNSC();


  
  /*! 
    \brief get generator polynomial from a predefined table for RSC
    
    Based on \a k and \a L , the appropriate generator polynomials are
    loaded into GenPoly of \a TRELLIS_parameters.
  */
  void get_RSC_poly();

  /*! 
    \brief get generator polynomial and puncturing pattern from a predefined table for BICM/BICMID
    
    Based on \a k and \a L , the appropriate generator polynomials are
    loaded into GenPolyMat of \a TRELLIS_parameters.
  */
  void get_NSC_poly(imat &GenPolyMat);

  /*! 
    \brief decode symbol based on symbol's a-poteriori probability
    
    Based on data symbol's a-poteriori probability, symbol block
    length and number of possible different input symbols, the most
    likely input data symbol sequence is given.

    \param Apo      data symbol's a-poteriori probability (\a Apo of \a MAP_dec or \a Apo_dataword \a SISO_dec)
    \param symbols  the most likely data symbol sequence
    \param N        symbol block length (\a no_of_symbols)
    \param M        number of possible different input data symbols (\a 2^k)
  */
  void decode_symbol(mat Apo, ivec &symbols, int N, int M);

 private:
  char result_filename[500];  //!< name of result file
  int mode;               //!< CONV mode: \ref NSC, \ref SC, \ref RSC, \ref RNSC
  int Terminated;         //!< code termination: \ref OFF, \ref ON
  int decoder_type;       //!< decoding type: \ref C_EXACT_LOG_MAP, \ref C_APPROX_LOG_MAP, \ref C_MAX_LOG_MAP

  int k;                  //!< length of dataword
  int n;                  //!< length of codeword
  int L;                  //!< code memory length
  
  int max_no_of_symbols;  //!< maximum number of coded symbols in a block, for variable length CM

  int no_of_symbols;      //!< number of coded symbols in a block 
  int no_of_info_symbols; //!< number of information (input data) symbols in one block 
  int no_of_coded_bits;   //!< number of output coded bits in one block
  int no_of_info_bits;    //!< number of information (input data) bits in one block 
  int no_of_tail_bits;    //!< number of termination (tail) bits in one block 
  int no_of_input_bits;   //!< number of input bits (\a no_of_info_bits + \a no_of_tail_bits ) in one block 

  //! printout error message and exit
  void s_error ( const char * message );
  //! correct rounding error to the precision of 6 decimal places
  double correct_rounding_error(double distance);
  //! compute decimal number from octal number
  int octal_to_decimal(int gen);
  


};

//! struct of the trellis states
typedef struct trellis_state{
  int index;    //!< state index
  int word;     //!< codeword 
  int infoword; //!< information word 
}Trellis_State;



#endif
