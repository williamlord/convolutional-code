/*!
 \file
 \brief Provide control file content access
 \author Tong Hooi LIEW

 1.01

 2003/03/26
*/

/*! \defgroup file_service File Service
*/

#ifndef FILESER
#define FILESER

#include <stdio.h>
#include <stdlib.h>
#include "itpp/itbase.h"

/*!
 \ingroup file_service
 \brief FILESERVICE class


 \warning
 While applying \a scan_text and \a scan_word, extra care needs to be taken due to the default value of LENGTH of 500. This is sufficient for 
 conventional simulation control parameter handling. However, while it is desired to read more than 500 characters in a single line, the value 
 of LENGTH should be adjusted.

*/
 
using namespace itpp;


class FILESERVICE
{
	public:

		//! Class Constructor, sets the \a LENGTH to 500;
		FILESERVICE(){LENGTH = 500;}
		
		//!Scans an integer from the file \b InputFile
		int scan_integer(FILE *InputFile);

		//!Scans a string from the file \b InputFile and stores in the string buffer \b text. The length of the string is defined as LENGTH
		void scan_text(FILE *InputFile, char *text);	

		//!Scans a double from the file \b InputFile
		double scan_double(FILE *InputFile);

		//!Scans a unsigned long number from the file \b InputFile
		unsigned long scan_unsigned_long(FILE *InputFile);

		//!Scans a word from the string buffer \b FileLine, the starting position will be \b StartPosition, and store in \b text
		int scan_word(int StartPosition, char *FileLine, char *text);


	private:
		//! The length of the input string length
		int LENGTH ;
};

//void output_parameters(FILE*, struct source_information*);

unsigned long bin_to_decimal(ivec& ptrShiftRegister, int NMinusK);

void decimal_to_bin(ivec& ptrShiftRegister, unsigned long DecimalNotation, int NMinusK);

#endif
