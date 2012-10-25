#include <stdio.h>
#include <stdlib.h>
#include "fileser.h"

int
FILESERVICE::scan_integer(FILE *InputFile)
{
    int integer;
    char FileLine[LENGTH];
    char *ptr;
    
    ptr=fgets(FileLine, LENGTH, InputFile);

    if(ptr==NULL)
    {
	fprintf(stderr, "Error in reading the control file.\n");
	exit(-1);
    }
    
    sscanf(FileLine,"%i", &integer);
    return(integer);
}

void
FILESERVICE::scan_text(FILE *InputFile, char *text)
{
    int i, length; 
    char FileLine[LENGTH];
    char *ptr;

    ptr=fgets(FileLine, LENGTH, InputFile);

    if(ptr==NULL)
    {
	fprintf(stderr, "Error in reading the control file.\n");
	exit(-1);
    }
    length=strlen(FileLine);
    for(i=0;i<length;i++)
    {
        if(FileLine[i]=='\0' || FileLine[i] == ' ' || FileLine[i] == '\n')
            break;
        
        text[i]=FileLine[i];
    }
    text[i]='\0';
/*    sscanf(FileLine,"%s", text);*/
}

/*---------------------------------------------------------------------------*/
/* scan in a double from input file. */

double
FILESERVICE::scan_double(FILE *InputFile)
{
    char FileLine[LENGTH];
    char *ptr;
    double FloatNumber;
    
    ptr=fgets(FileLine, LENGTH, InputFile);

    if(ptr==NULL)
    {
	fprintf(stderr, "Error in reading the control file.\n");
	exit(-1);
    }
    
    sscanf(FileLine,"%lf", &FloatNumber);

    return(FloatNumber);
}
/* 22/7/98 */
/*---------------------------------------------------------------------------*/
/* scan in an unsigned long from input file. */

unsigned long
FILESERVICE::scan_unsigned_long(FILE *InputFile)
{
    char FileLine[LENGTH];
    char *ptr;
    unsigned long BigInteger;
    
    ptr=fgets(FileLine, LENGTH, InputFile);

    if(ptr==NULL)
    {
	fprintf(stderr, "Error in reading the control file.\n");
	exit(-1);
    }
    
    sscanf(FileLine,"%lu", &BigInteger);

    return(BigInteger);
}
/* 22/7/98 */

/*---------------------------------------------------------------------------*/
/* get a word from a string */

int
FILESERVICE::scan_word(int StartPosition, char *FileLine, char *text)
{
    int i, EndPosition, length;

    length=strlen(FileLine);

    for(EndPosition=StartPosition,i=0;EndPosition<length;EndPosition++,i++)
    {
        if(FileLine[EndPosition]=='\0' || FileLine[EndPosition] == ' ' ||
           FileLine[EndPosition] == '\n')
            break;
        
        text[i]=FileLine[EndPosition];
    }
    text[i]='\0';
    
    return(EndPosition+1);
}

/*
void output_parameters(FILE *pf, struct source_information *ptrSI)
{
	int i,j;

	fprintf(pf, "!Simulation of Variable Length Code MAP decoding.......\n");
	
	for(i = 0; i < ptrSI->nrof_possible_symbols; i++)
	{
		fprintf(pf, "!Symbol No %d\tSymbol Sequence is ", i);
		
		j = 0;
		while(ptrSI->symbol[i][j]!=-1)
		{
			fprintf(pf, "%d ",ptrSI->symbol[i][j]);
			j++;
		}

		fprintf(pf, "\tProbability is %f\n", ptrSI->symbol_probability[i]);

	}
	
	fflush(pf);
}
*/

unsigned long bin_to_decimal(ivec& ptrShiftRegister, int NMinusK)
{
	int i;
	unsigned long DecimalNotation,  multiplier;
	
	multiplier=1;
	DecimalNotation=0;
	
	for(i=NMinusK-1;i>=0;i--)
	{
		DecimalNotation+=ptrShiftRegister[i]*multiplier;
		multiplier=multiplier*2;
	}

	
	return(DecimalNotation);
}



void decimal_to_bin(ivec& ptrShiftRegister, unsigned long DecimalNotation, int NMinusK)
{
	    int i;
	    int quot;
	    int rem;
	    
	    quot=DecimalNotation;
	    
	    for(i=NMinusK-1;i>=0;i--)
	    {
		    rem=quot%2;
		    quot=quot/2;
		    
		    ptrShiftRegister[i]=rem;
	    }
}
