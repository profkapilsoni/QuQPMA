/*
*************************************************************
*	About: The code implements Algorithm 3					*
*	Usage: Run in command prompt "algo3.exe <input_file>"	*
*	<input_file> contains biological sequence in which 		*
*	searching must be performed for the fixed pattern P		*
*	The quantum circuit is designed for P = "GCCTC"			*
*************************************************************
*/

//The following implementation does not use ancilla qubits
// instead phase flip approach is employed 

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "QuEST.h"

char *T;
int T_size;

int is_at(int idx, char *P, int P_size)
{
	int found = 1;
	int temp_i = 0;
	if((idx + P_size) > T_size) return 0;
	for(int i = idx; i < idx + P_size; i++)
		if(T[i] != P[temp_i++]) found = 0;
	return found;
}

int main (int narg, char *varg[])
{
	if(narg != 2)
	{
		printf("Usage: algo3.exe <input_file>");
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
	
	T = (char *)malloc(sizeof(char) * count);
	T_size = count;
	count = 0;
	fp = fopen(varg[1], "r");
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			T[count++] = ch;
	}
	fclose(fp);

	//ENDS - Reading and Storing Input File in the array "T"

	//P = "GCCTC" must be searched in "T"
	const int P_size = 5;
	char P[] =  { 'G', 'C', 'C', 'T', 'C' };

	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= count; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}
	printf("Pattern to search: ");
	for(int i = 0; i < P_size; i++)
		printf("%c", P[i]);
	printf("\n");

	//STARTS - Implementation of Quantum Exact Pattern Search

	//Number of qubits in the system without using Ancilla qubits
	int n = (int)(log(T_size)/log(2));

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(n, env);
    initZeroState(qubits);
    printf("\nQuantum Simulation Parameters for Quantum Exact Pattern Match Algorithm:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    count = 0;
    //Unitary matrix creation corresponding to the pattern match
	ComplexMatrixN e = createComplexMatrixN(n);
    for(int i = 0; i < (int)pow(2, n); i++)
    {
    	if(is_at(i, P, P_size))
    	{
    		e.real[i][i] = -1;
    		count++;
    	}
    	else e.real[i][i] = 1;
    }

    //Apply Hadamard on all the n-qubits to create superposition
    //Here each element in superposition is considered to represent "x" 
    for(int i = 0; i < n; i++)
        hadamard(qubits, i);

    int *targs = (int *)malloc(sizeof(int) * n);
    for(int i = 0; i < n; i++)
    	targs[i] = i;

    int times = (int)(3.14 * (pow(2, n / 2) / sqrt(count)) / 4);
    for(int gi = 0; gi < times; gi++)
    {
    	//Marking
    	multiQubitUnitary(qubits, targs, n, e);
    	
        //Diffusion
        for(int i = 0; i < n; i++)
            hadamard(qubits, i);
        for(int i = 0; i < n; i++)
            pauliX(qubits, i);
        multiControlledPhaseFlip(qubits, targs, n);
        for(int i = 0; i < n; i++)
            pauliX(qubits, i);
        for(int i = 0; i < n; i++)
            hadamard(qubits, i);
    }

    //Measuring the result
    qreal prob;
    int index = 0;
    for(int i = 0; i < n; i++)
    {
    	int outcome = measureWithStats(qubits, i, &prob);
		if(outcome) index ^= (outcome << i);
    }
    printf("\nPattern Found at Index: %d\n", index);

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    //ENDS - Implementation of Quantum Exact Pattern Search

    return 0;
}