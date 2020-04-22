#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

//Used in building Algebraic Normal Form
int check_term(int t, int a, int width)
{
	for (int i = 0; i < width; i++)
	{
		if ((t & (1 << i)) == 0 && (a & (1 << i)) != 0) return 0;
	}
	return 1;
}

//Builds ANF for a function which takes m-bit input and produces n-bit output
void anf_calc(int f[], int m, int n, int ***anf, int **anf_size)
{
	*anf = (int **)malloc(sizeof(int *) * n);
	for(int i = 0; i < n; i++)
		(*anf)[i] = (int *)malloc(sizeof(int) * (int)pow(2, m));
	*anf_size = (int *)malloc(sizeof(int) * n);
	for(int i = 0; i < n; i++)
		(*anf_size)[i] = 0;

	for (int i = n - 1; i >= 0; i--)
	{
		for (int term = 0; term < (int)pow(2, m); term++)
		{
			int count = 0;
			for (int fi = 0; fi < (int)pow(2, m); fi++)
			{
				if (check_term(term, fi, m) && (f[fi] & (1 << i)) != 0) count++;
			}
			if (count % 2 == 1) (*anf)[n - 1 - i][(*anf_size)[n - 1 - i]++] = term;
		}
	}
}

int main (int narg, char *varg[])
{
	if(narg != 2)
	{
		printf("Usage: %s <input_file>", varg[0]);
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
	int T_size = count;
	count = 0;
	fp = fopen(varg[1], "r");
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			T[count++] = ch;
	}
	fclose(fp);

	//ENDS - Reading and Storing Input File in the array "T"
	
	//"P" must be searched in "T"
	char P[] =  { 'A' };
	int P_size = 1;

	//Printing all the inputs
	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= T_size; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}
	printf("\nPattern to search are:\n");
	for(int i = 0; i < P_size; i++)
		printf("Pattern No. %d: %c\n", i + 1, P[i]);

	//Display Warning if necessary
	
	printf("\nPatterns will searched in parallel with the following parameters information:\n");

	int n = (int)(log(T_size)/log(2));

	ch = 'y';
	int req_qubits = n + 2 + 2 + 1;
	if(req_qubits > 25)
	{
		printf("\nWARNING: We do not recommend to run the \"%s\" Algorithm since it requires %d qubits!", varg[0], req_qubits);
		printf("\nEnter your choice if you want to continue (y/n): ");
		ch = getch();
		ch = tolower(ch);
	}
	else printf("\nThe \"%s\" Algorithm requires %d qubits.\n", varg[0], req_qubits);

	if(ch != 'y')
	{
		printf("\nExiting!\n");
		return 0;
	}
	else if(ch == 'y' && req_qubits > 25)
	{
		printf("\n\nYou chose to ignore our warning and continue!\nGOOD LUCK!\n");
	}

	//STARTS - Preliminaries for quantum circuit
	//NOTE: Not the efficient approach, can be improved

	int **anf, *anf_size;
    int *f = (int *)malloc(sizeof(int) * (int)pow(2, n));
    for(int i = 0; i < T_size; i++) //T_size = 2^n
    {
    	if(tolower(T[i]) == 'a') f[i] = 0;
    	else if(tolower(T[i]) == 't') f[i] = 1;
    	else if(tolower(T[i]) == 'c') f[i] = 2;
    	else if(tolower(T[i]) == 'g') f[i] = 3;
    }
    anf_calc(f, n, 2, &anf, &anf_size);

	//ENDS - Preliminaries for quantum circuit

	//STARTS

    QuESTEnv env = createQuESTEnv();
    Qureg qubits = createQureg(req_qubits, env);
    printf("\nQuantum Parameters for \"%s\":\n", varg[0]);
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    //Single qubit unitary for XOR gate
    ComplexMatrix2 ux = {
    	.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
    };

    for(int np = 0; np < P_size; np++)
    {
    	initZeroState(qubits);

	    //Print status of quantum state construction
	    printf("\nConstructing required quantum state to search \"%c\"..\n", P[np]);

	    //Create superposition of all indexes from qubits indexed from 0 to (n - 1)
	    for(int i = 0; i < n; i++)
	    	hadamard(qubits, i);

	    //Store "d_i" corresponding to each "i"
	    //"d_i" are stored in qubits indexed from (n + 1) to (n)
	    int set_qb = (n + 1);
	    for(int idx = 0; idx < 2; idx++)
	    {
	    	for(int i = 0; i < anf_size[idx]; i++)
	    	{
	    		if(anf[idx][i] == 0)
	    		{
	    			pauliX(qubits, set_qb);
	    			continue;
	    		}
	    		int *ctrls = (int *)malloc(sizeof(int) * n);
	    		int ctrl_size = 0;
	    		int term = anf[idx][i];
	    		int qb = 0;
	    		while(term)
	    		{
	    			if(term & 1) ctrls[ctrl_size++] = qb;
	    			term >>= 1;
	    			qb++;
	    		}
	    		multiControlledUnitary(qubits, ctrls, ctrl_size, set_qb, ux);
	    		free(ctrls);
	    	}
	    	set_qb--;
	    }

	    //Store element "e" of "E" set in qubits indexed from (n + 3) to (n + 2)
	    set_qb = (n + 2 + 2);
    	int set = 0;
    	if(tolower(P[np]) == 'a') set = 0;
    	else if(tolower(P[np]) == 't') set = 1;
    	else if(tolower(P[np]) == 'c') set = 2;
    	else if(tolower(P[np]) == 'g') set = 3;
    	if(set & 2) pauliX(qubits, set_qb - 1);
    	if(set & 1) pauliX(qubits, set_qb - 2);

    	//Evaluate "d_i" == "e"
    	int eval_form[][4] = { { 0, 0, 0, 0 }, { 0, 0, 1, 1 }, { 1, 1, 0, 0}, {1, 1, 1, 1} };
		int anc_ctrls[] = { n + 3, n + 1, n + 2, n };
		multiStateControlledUnitary(qubits, anc_ctrls, eval_form[0], 4, n + 4, ux);
		multiStateControlledUnitary(qubits, anc_ctrls, eval_form[1], 4, n + 4, ux);
		multiStateControlledUnitary(qubits, anc_ctrls, eval_form[2], 4, n + 4, ux);
		multiStateControlledUnitary(qubits, anc_ctrls, eval_form[3], 4, n + 4, ux);
		
		//Prints the quantum state before applying Grover's algorithm
		printf("Required Quantum State is constructed.\n");
		printf("Do you want to view the constructed quantum state?(y/n):");
		char ch = getch();
		if(ch == 'y')
		{
			printf("\nConstructed quantum state before applying Grover's search is:");
			qreal prob;
		    int mask = (int)pow(2, n) - 1;
		    for(int i = 0; i < (int)pow(2, n + 2 + 2 + 1); i++)
		    {
		    	prob = getProbAmp(qubits, i);
		    	if(prob != 0.0)
		    		printf("\ni = %d d = %X e = %X (d_i == e) = %X has prob %f", 
		    			(i & mask),
		    			(i & ((1 << ((n + 1)) ^ (1 << n)))) >> n,
		    			(i >> (n + 2)) & 3,
		    			(i >> (n + 2 + 2)),
		    			prob);
		    }
		}

		//Removing d_i and e
		set_qb = (n + 2 + 2);
    	set = 0;
    	if(tolower(P[np]) == 'a') set = 0;
    	else if(tolower(P[np]) == 't') set = 1;
    	else if(tolower(P[np]) == 'c') set = 2;
    	else if(tolower(P[np]) == 'g') set = 3;
    	if(set & 2) pauliX(qubits, set_qb - 1);
    	if(set & 1) pauliX(qubits, set_qb - 2);
    	set_qb = (n + 1);
	    for(int idx = 0; idx < 2; idx++)
	    {
	    	for(int i = 0; i < anf_size[idx]; i++)
	    	{
	    		if(anf[idx][i] == 0)
	    		{
	    			pauliX(qubits, set_qb);
	    			continue;
	    		}
	    		int *ctrls = (int *)malloc(sizeof(int) * n);
	    		int ctrl_size = 0;
	    		int term = anf[idx][i];
	    		int qb = 0;
	    		while(term)
	    		{
	    			if(term & 1) ctrls[ctrl_size++] = qb;
	    			term >>= 1;
	    			qb++;
	    		}
	    		multiControlledUnitary(qubits, ctrls, ctrl_size, set_qb, ux);
	    		free(ctrls);
	    	}
	    	set_qb--;
	    }

	    count = 0;
	    ComplexMatrixN e = createComplexMatrixN(n);
		for(int i = 0; i < (int)pow(2, n); i++)
		{
		  	if(T[i] == P[np])
		   	{
		   		e.real[i][i] = -1;
		   		count++;
		   	}
		   	else e.real[i][i] = 1;
		}

		//Running Grover's search
		int *targs = (int *)malloc(sizeof(int) * n);
		for(int i = 0; i < n; i++)
			targs[i] = i;
	    printf("\nRunning Grover's Algorithm..\n");
	    if(count != 0)
	    {
		    int times =  (int)(3.14 * (pow(2, n / 2) / sqrt(count)) / 4);
		    if(times == 0) times = 1;
		    for(int gi = 0; gi < times; gi++)
		    {
		    	//Marking
		    	controlledMultiQubitUnitary(qubits, n + 4, targs, n, e);
		    	//multiControlledPhaseFlip(qubits, ctrls, n + 1);

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

		    qreal prob, max = 0.0;
		    for(int i = 0; i < (int)pow(2, n + 2 + 2 + 1); i++)
		    {
		    	prob = getProbAmp(qubits, i);
		    	if(max <= prob) max = prob;
		    }
		    int mask = (int)pow(2, n) - 1;
		    count = 0;
		    for(int i = 0; i < (int)pow(2, n + 2 + 2 + 1); i++)
		    {
		    	prob = getProbAmp(qubits, i);
		    	if(fabs(max - prob) <= 0.0000001)
		    	{
		    		printf("Correct Index %d\n", (i & mask));
		    		break;
		    	}
		    }
		}
	    else printf("No index exist containing the pattern.\n");

	    destroyComplexMatrixN(e);
	}

	destroyQureg(qubits, env);
    destroyQuESTEnv(env);

    //ENDS

    return 0;
}