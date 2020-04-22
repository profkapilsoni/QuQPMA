/*
*************************************************************
*	About: The code implements Algorithm 2					*
*	Usage: Run in command prompt "algo2.exe <input_file>"	*
*	<input_file> contains biological sequence in which 		*
*	searching must be performed for the fixed pattern P		*
*	The quantum circuit is designed for P = "GCCTC"			*
*************************************************************
*/

//Instead of applying unitaries "U_Loc" and "U_Sub" using their corresponding Unitary matrix,
// a corresponding quantum circuit is applied
//NOTE: This is done to make the simulation memory efficient

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

//The function gen_Psym is not optimized and runs in brute-force 
int gen_Psym(char P[], int size, char **Psym)
{
	*Psym = (char *)malloc(sizeof(char) * size);
	int count = 0;
	for(int i = 0; i < size; i++)
	{
		int found = 0;
		for(int j = 0; j < count; j++)
		{
			if((*Psym)[j] == P[i]) found = 1;
		}
		if(!found) (*Psym)[count++] = P[i];
	}
	return count;
}

//Used in building Algebraic Normal Form
int check_term(int t, int a, int width)
{
	for (int i = 0; i < width; i++)
	{
		if ((t & (1 << i)) == 0 && (a & (1 << i)) != 0) return 0;
	}
	return 1;
}

//Builds ANF
void anf_calc(int f[], int n, int ***anf, int **anf_size)
{
	*anf = (int **)malloc(sizeof(int *) * n);
	for(int i = 0; i < n; i++)
		(*anf)[i] = (int *)malloc(sizeof(int) * (int)pow(2, n));
	*anf_size = (int *)malloc(sizeof(int) * n);
	for(int i = 0; i < n; i++)
		(*anf_size)[i] = 0;

	for (int i = n - 1; i >= 0; i--)
	{
		for (int term = 0; term < (int)pow(2, n); term++)
		{
			int count = 0;
			for (int fi = 0; fi < (int)pow(2, n); fi++)
			{
				if (check_term(term, fi, n) && (f[fi] & (1 << i)) != 0) count++;
			}
			if (count % 2 == 1) (*anf)[i][(*anf_size)[i]++] = term;
		}
	}
}

//Return 2's complement of an integer of "n" bits
int _2s_complement(int num, int bits)
{
	for(int i = 0; i < bits; i++)
		num ^= (1 << i);
	return (num + 1);
}

//Binary Adder
int bin_add(int a, int b, int width)
{
	int carry = 0, sum = 0;
	for(int i = 0; i < width; i++)
	{
		int _a = ((a & (1 << i)) >> i);
		int _b = ((b & (1 << i)) >> i);
		sum ^= ((_a ^ _b ^ carry) << i);
		if((_a + _b + carry) >= 2)
			carry = 1;
		else carry = 0;
	}
	return sum;
}


//STARTS - Quantum Implementation of Filtering Part

int filtering(int n, int **anf, int *anf_size, int Psym_size, int **Loc)
{
	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(2 * n, env);
    initZeroState(qubits);
    printf("\nQuantum Parameters for Filtering Part:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    //Single qubit unitary for XOR gate
    ComplexMatrix2 ux = {
    	.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
    };

    for(int i = n; i < 2 * n; i++)
    	hadamard(qubits, i);

    //Evaluation of U_Loc and U_Sub by below quantum circuit
    for(int idx = 0; idx < n; idx++)
    {
    	for(int i = 0; i < anf_size[idx]; i++)
    	{
    		if(anf[idx][i] == 0)
    		{
    			pauliX(qubits, idx);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n);
    		int ctrl_size = 0;
    		int term = anf[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = (qb + n);
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits, ctrls, ctrl_size, idx, ux);
    		free(ctrls);
    	}
    }

    for(int i = n; i < 2 * n; i++)
    	hadamard(qubits, i);

    //Reading all the highest probability indexes corresponding to the threshold
    int d = Psym_size - 1; //Threshold of the approximate search
    double req_prob = (double)d / (double)pow(2, n);
    req_prob *= req_prob;
    //printf("Threshold probability is %f\n", req_prob);
    //Counting and storing likely indexes
    int res_count = 0;
    int realloc_size = 1;
    (*Loc) = (int *)malloc(sizeof(int) * realloc_size);
    qreal prob;
    for(int i = 0; i < (int)pow(2, n); i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(prob >= req_prob)
    	{
    		if(realloc_size == res_count)
    		{
    			realloc_size *= 2;
    			(*Loc) = (int *)realloc((*Loc), sizeof(int) * realloc_size);
    		}
    		(*Loc)[res_count++] = i;
    	}
    }

	destroyQureg(qubits, env);
    destroyQuESTEnv(env);

    return res_count;
}

//ENDS - Quantum Implementation of Filtering Part

//STARTS - Quantum Implementation of Verification Part

int is_in(int elem, int *a, int n)
{
	for(int i = 0; i < n; i++)
		if(elem == a[i]) return 1;
	return 0;
}

int is_at(char *T, int T_size, int idx, char *P, int P_size)
{
	int found = 1;
	int temp_i = 0;
	if((idx + P_size) > T_size) return 0;
	for(int i = idx; i < idx + P_size; i++)
		if(T[i] != P[temp_i++]) found = 0;
	return found;
}

int verification(char *T, int N, char *P, int M, int *Loc, int Loc_size, int d)
{
	QuESTEnv env = createQuESTEnv();

	int _1 = (int)(log(N)/log(2));

    Qureg qubits = createQureg(_1, env);
    initZeroState(qubits);
    printf("\nQuantum Parameters for Verification Part:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);
    
    //Creating |delta_{init}> state
    //Finding the correct "i"
    int *targs = (int *)malloc(sizeof(int) * _1);
	for(int i = 0; i < _1; i++)
	{
		hadamard(qubits, i);
		targs[i] = i;
	}
	int temp = 0;
	ComplexMatrixN e = createComplexMatrixN(_1);
    for(int i = 0; i < (int)pow(2, _1); i++)
    {
    	if(is_in(i, Loc, Loc_size) && is_at(T, N, i, P, M))
    	{
    		e.real[i][i] = -1;
    		temp++;
    	}
    	else e.real[i][i] = 1;
    }
    int times = (int)(3.14 / 4 * sqrt(N) / sqrt(temp));
    for(int gi = 0; gi < times; gi++)
    {
    	multiQubitUnitary(qubits, targs, _1, e);
        for(int i = 0; i < _1; i++)
            hadamard(qubits, i);
        for(int i = 0; i < _1; i++)
            pauliX(qubits, i);
        multiControlledPhaseFlip(qubits, targs, _1);
        for(int i = 0; i < _1; i++)
            pauliX(qubits, i);
        for(int i = 0; i < _1; i++)
            hadamard(qubits, i);
    }

    qreal prob;
    int index = 0;
    for(int i = 0; i < _1; i++)
    {
    	int outcome = measureWithStats(qubits, i, &prob);
		if(outcome) index ^= (outcome << i);
    }

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

	return index;
}

//ENDS - Quantum Implementation of Verification Part


int main (int narg, char *varg[])
{
	if(narg != 2)
	{
		printf("Usage: algo2.exe <input_file>");
		return 0;
	}
	
	//STARTS - Reading and Storing Input File in the array "T"
	//The code assumes that the number of elements in the input file is in power of 2
	//for further easier processing

	FILE *fp = fopen(varg[1], "r");
	int count = 0;
	char ch;
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			count++;
	}
	fclose(fp);
	
	char *T = (char *)malloc(sizeof(char) * count);
	int idx = 0;
	fp = fopen(varg[1], "r");
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			T[idx++] = ch;
	}
	fclose(fp);

	//ENDS - Reading and Storing Input File in the array "T"

	//STARTS - Finds "Psym" and "i_j" corresponding to the pattern "P"
	
	//P = "GCCTC" must be searched in "T"
	const int P_size = 5;
	char P[] =  { 'G', 'C', 'C', 'T', 'C' };
	char *Psym;
	int Psym_size = gen_Psym(P, P_size, &Psym); //Generates Psym corresponding to P

	int *i_j = (int *)malloc(sizeof(int) * Psym_size);
	for(int i = 0; i < Psym_size; i++)
	{
		for(int j = 0; j < P_size; j++)
		{
			if(Psym[i] == P[j])
			{
				i_j[i] = j;
				break;
			}
		}
	}

	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= count; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}
	printf("Pattern to search: ");
	for(int i = 0; i < P_size; i++)
		printf("%c", P[i]);
	printf("\nPsym is: ");
	for(int i = 0; i < Psym_size; i++)
		printf("%c", Psym[i]);
	printf("\ni_j is: ");
	for(int i = 0; i < Psym_size; i++)
		printf("%d ", i_j[i]);
	printf("\n");

	//ENDS - Finds "Psym" and "i_j" corresponding to the pattern "P"

	//STARTS - Preliminaries for quantum circuit corresponding to U_Loc and U_Sub
	//NOTE: Not the efficient approach, can be improved
	
	int n = (int)(log(count)/log(2));

	int **anf, *anf_size;
    int *f = (int *)malloc(sizeof(int) * (int)pow(2, n));
    for(int i = 0; i < count; i++) //count = 2^n
    {
    	idx = 0;
    	for(int j = 0; j < Psym_size; j++)
    	{	
    		if(tolower(T[i]) == tolower(Psym[j])) 
    		{
    			idx = i_j[j];
    			break;
    		}
    	}
    	f[i] = bin_add(i, _2s_complement(idx, n), n);
    }
    anf_calc(f, n, &anf, &anf_size);

	//ENDS - Preliminaries for quantum circuit corresponding to U_Loc

	//STARTS - Quantum filtering call

    int *Loc;
	int res_count = filtering(n, anf, anf_size, Psym_size, &Loc);
	printf("\nLikely indexes are:\n");
	for(int i = 0; i < res_count; i++)
		printf("Index %d\n", Loc[i]);

	//ENDS - Quantum filtering call

	//STARTS - Quantum Verification call

	char choice = 'y';
	while(choice == 'y')
	{
		int found_index = verification(T, count, P, P_size, Loc, res_count, Psym_size - 1);
		printf("\nFound a highly likely correct index: %d\n", found_index);
		printf("Do you want to find another correct index!\nRepeating Verification doesn't guarantee new correct index.\nEnter choice (y/n): ");
		choice = getch();
		choice = tolower(choice);
		printf("\n");
	}

	//ENDS - Quantum Verification call

    return 0;
}