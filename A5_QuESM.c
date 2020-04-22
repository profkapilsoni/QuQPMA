#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

//Generates all binary string of length "n" and store in "_bin_str"
int **_bin_str;
int bin_idx = 0;
void gen_all_bin_str(int n, int a[], int i)
{
	if(i == n)
	{
		for(int i = 0; i < n; i++)
			_bin_str[bin_idx][i] = a[i];
		bin_idx++;
		return;
	}
	a[i] = 0;
	gen_all_bin_str(n, a, i + 1);
	a[i] = 1;
	gen_all_bin_str(n, a, i + 1);
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

int is_at(char *T, int T_size, int idx, char *P, int P_size)
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
	int M = 1;
	const int P_size = 3;
	char P[][4] =  { {'G', 'C', 'T', '\0' } };

	//Printing all the inputs
	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= T_size; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}
	printf("\nPattern to search are:\n");
	for(int i = 0; i < M; i++)
		printf("Pattern No. %d: %s\n", i + 1, P[i]);

	//Display Warning if necessary
	
	printf("\nPatterns will searched in parallel with the following parameters information:\n");

	int n = (int)(log(T_size)/log(2));

	ch = 'y';
	int req_qubits = n + 4 * P_size + 1;
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

	printf("Started Classical Preprocessing..\n");
	
	int ***anf, **anf_size;
	anf = (int ***)malloc(sizeof(int **) * P_size);
	anf_size = (int **)malloc(sizeof(int *) * P_size);
	int **f = (int **)malloc(sizeof(int *) * P_size);
	for(int i = 0; i < P_size; i++)
		f[i] = (int *)malloc(sizeof(int) * (int)pow(2, n));
	for(int i = 0; i < T_size; i++)
	{
	 	for(int j = 0; j < P_size; j++)
	   	{
	    	if(tolower(T[i + j]) == 'a') f[j][i] = 0;
	   		else if(tolower(T[i + j]) == 't') f[j][i] = 1;
	   		else if(tolower(T[i + j]) == 'c') f[j][i] = 2;
	   		else if(tolower(T[i + j]) == 'g') f[j][i] = 3;
	   	}
	}
	for(int i = 0; i < P_size; i++)
		anf_calc(f[i], n, 2, &anf[i], &anf_size[i]);
	
	_bin_str = (int **)malloc(sizeof(int *) * (int)pow(2, 2 * P_size));
	for(int i = 0; i < (int)pow(2, 2 * P_size); i++)
		_bin_str[i] = (int *)malloc(sizeof(int) * 2 * P_size);
	int temp_a[6]; //Set accordingly = (2 * P_size)
	gen_all_bin_str(2 * P_size, temp_a, 0);
	
	printf("Completed Classical Preprocessing!\n");

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

    for(int np = 0; np < M; np++)
    {
    	initZeroState(qubits);

	    //Print status of quantum state construction
	    printf("\nConstructing required quantum state to search \"%s\"..\n", P[np]);

	    //Create superposition of all indexes from qubits indexed from 0 to (n - 1)
	    for(int i = 0; i < n; i++)
	    	hadamard(qubits, i);

	    //Store T_i to T_{i + P_size - 1}
	    int set_qb = n + 2 * P_size - 1;
	    for(int ip = 0; ip < P_size; ip++)
	    {
		    for(int idx = 0; idx < 2; idx++)
		    {
		    	for(int i = 0; i < anf_size[ip][idx]; i++)
		    	{
		    		if(anf[ip][idx][i] == 0)
		    		{
		    			pauliX(qubits, set_qb);
		    			continue;
		    		}
		    		int *ctrls = (int *)malloc(sizeof(int) * n);
		    		int ctrl_size = 0;
		    		int term = anf[ip][idx][i];
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
		}

	    //Store "P"
	    set_qb = (n + 4 * P_size - 1);
    	for(int i = 0; i < P_size; i++)
    	{
    		int set = 0;
	    	if(tolower(P[np][i]) == 'a') set = 0;
	    	else if(tolower(P[np][i]) == 't') set = 1;
	    	else if(tolower(P[np][i]) == 'c') set = 2;
	    	else if(tolower(P[np][i]) == 'g') set = 3;
	    	if(set & 2) pauliX(qubits, set_qb);
	    	if(set & 1) pauliX(qubits, set_qb - 1);
	    	set_qb -= 2;
    	}

    	//Evaluate "T[i, i + P_size - 1] == P"
		int *anc_ctrls = (int *)malloc(sizeof(int) * 4 * P_size);
		set_qb = (n + 4 * P_size - 1);
		int idx = 0;
		for(int i = 0; i < 2 * P_size; i++)
		{
			anc_ctrls[idx++] = set_qb - i;
			anc_ctrls[idx++] = set_qb - i - (2 * P_size);
		}
		for(int i = 0; i < bin_idx; i++)
		{
			int *eval_form = (int *)malloc(sizeof(int) * 4 * P_size);
			int eval_idx = 0;
			for(int j = 0; j < 2 * P_size; j++)
			{
				eval_form[eval_idx++] = _bin_str[i][j];
				eval_form[eval_idx++] = _bin_str[i][j];
			}
			multiStateControlledUnitary(qubits, anc_ctrls, eval_form, 4 * P_size, n + 4 * P_size, ux);
			free(eval_form);
		}
		
		//Prints the quantum state before applying Grover's algorithm
		printf("Required Quantum State is constructed.\n");
		printf("Do you want to view the constructed quantum state?(y/n):");
		char ch = getch();
		if(ch == 'y')
		{
			printf("\nConstructed quantum state before applying Grover's search is:");
		    qreal prob;
		    int mask = (int)pow(2, n) - 1;
		    char Alpha[] = { 'A', 'T', 'C', 'G' };
		    for(int i = 0; i < (int)pow(2, n + 4 * P_size + 1); i++)
		    {
		    	prob = getProbAmp(qubits, i);
		    	if(prob != 0.0)
		    	{
		    		printf("\ni = %d, T[i, i + P_size - 1] = ", (i & mask));
		    		for(int j = P_size - 1; j >= 0; j--)
		    			printf("%c", Alpha[(i >> (j * 2 + n)) & 0x3]);
		    		printf(", P = ");
		    		for(int j = 2 * P_size - 1; j >= P_size; j--)
		    			printf("%c", Alpha[(i >> (j * 2 + n)) & 0x3]);
		    		printf(", (T == P) = %d has prob %f", (i >> (n + 4 * P_size)), prob);
		    	}
		   	}
		}

		//Removing "T" and "P"
		set_qb = (n + 4 * P_size - 1);
    	for(int i = 0; i < P_size; i++)
    	{
    		int set = 0;
	    	if(tolower(P[np][i]) == 'a') set = 0;
	    	else if(tolower(P[np][i]) == 't') set = 1;
	    	else if(tolower(P[np][i]) == 'c') set = 2;
	    	else if(tolower(P[np][i]) == 'g') set = 3;
	    	if(set & 2) pauliX(qubits, set_qb);
	    	if(set & 1) pauliX(qubits, set_qb - 1);
	    	set_qb -= 2;
    	}
    	set_qb = n + 2 * P_size - 1;
	    for(int ip = 0; ip < P_size; ip++)
	    {
		    for(int idx = 0; idx < 2; idx++)
		    {
		    	for(int i = 0; i < anf_size[ip][idx]; i++)
		    	{
		    		if(anf[ip][idx][i] == 0)
		    		{
		    			pauliX(qubits, set_qb);
		    			continue;
		    		}
		    		int *ctrls = (int *)malloc(sizeof(int) * n);
		    		int ctrl_size = 0;
		    		int term = anf[ip][idx][i];
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
		}

	    count = 0;
	    ComplexMatrixN e = createComplexMatrixN(n);
		for(int i = 0; i < (int)pow(2, n); i++)
		{
		  	if(is_at(T, T_size, i, P[np], P_size))
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
		    if(times % 2 == 0) times++;
		    for(int gi = 0; gi < times; gi++)
		    {
		    	//Marking
		    	controlledMultiQubitUnitary(qubits, n + 4 * P_size, targs, n, e);

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

		free(anc_ctrls);
	    destroyComplexMatrixN(e);
	}

	for(int i = 0; i < (int)pow(2, 2 * P_size); i++)
		free(_bin_str[i]);
	free(_bin_str);
	destroyQureg(qubits, env);
    destroyQuESTEnv(env);

    //ENDS

    return 0;
}