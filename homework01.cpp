#include <stdio.h>
#include <stdlib.h>

//attribute number
#define n 18 
//gene number 
#define y 6178 
//array for gene names 
char name[y][100];
char name_new[y][100];


//function to write gene names and their functional expression data in a new file 
void print_gene_expression(char name[y][100], double **data, int rows) {
	FILE *out;
	int i, j;

	//new file for insertion 
	if ((out = fopen("gene_expression_no1.txt", "w")) == NULL) {
		printf("can not open out_file\n");
		exit(1);
	}

	//write gene names and their functional expression data in a new file 
	for (i = 0; i<rows; i++) {
		fprintf(out, "%s\t", name[i]);
		for (j = 0; j<n; j++) {
			fprintf(out, "%lf\t", data[i][j]);
		}
		fprintf(out, "\n");
	}
}



int main(void) {
	FILE *in;
	int i, j, k, x, idx[y][n], ptr, ptr1, sum, sum_count;
	double **data, **data_new, init, fin, avg;

	printf("\n %s\n %s\n", "Make sure that file yeast_microarray_data.txt exist", "Press enter to continue");
	getchar();

	if ((in = fopen("yeast_microarray_data.txt", "r")) == NULL) {
		printf("can not open in_file");
		exit(1);
	}

	//dynamic storage allocation of array "data" (for one dimension of array "data")
	data = (double **)malloc(sizeof(double *)*y);
	data_new = (double **)malloc(sizeof(double *)*y);
	for (i = 0; i<y; i++) {

		//dynamic storage allocation of array "data" (for two dimension of array "data")
		data[i] = (double *)malloc(sizeof(double)*n);
		data_new[i] = (double *)malloc(sizeof(double)*n);

		//store gene names in array "name" (e.g., YDL185W is gene name)
		fscanf(in, "%s", &name[i]);
		for (j = 0; j<n; j++) {

			//load gene functionnal expression data on array "data"
			fscanf(in, "%lf", &data[i][j]);
		}
	}
	//print_gene_expression(name, data);	

	// make index
	for (i = 0; i < y; i++) {
		for (j = 0; j < n; j++) {
			if (data[i][j] == 1000) {
				idx[i][j] = 1;
			}
			else {
				idx[i][j] = 0;
			}
		}
	}

	// replacing 1000
	init = 0; fin = 0; avg = 0; sum_count = 0;
	x = 0; // counter for data_new[][]
	
	for (i = 0; i<y; i++) { // counter for data[][]
		// bypass all missing data > 20% ('1000' > 3)
		sum = 0; // recount sum
		for (j = 0; j < n; j++) { sum += idx[i][j]; }
		if (sum > 3) { printf("%s %d\n","missing data > 20 persent on ",i); sum_count++; continue; }
		
		// check 1000 in begin
		ptr = 0;
		if (idx[i][0] == 1) {
			// cout << "Begining miss data on " <<i<< endl;
			// find until neighbor !0
			ptr = 0;
			while (ptr < n) {
				if (idx[i][ptr] == 0) {
					// copy idx[]=0 to idx[]=1
					for (k = 0; k < ptr; k++) {
						idx[i][k] = 0;
						data[i][k] = data[i][ptr];
					}
					break;
				} ptr++; // move until find idx[]=1
			}
		}

		// check 1000 in last
		ptr = 0;
		if (idx[i][n - 1] == 1) {
			// cout <<"Ending miss data on " << i << endl;
			// find until neighbor !0
			int ptr = n - 1;
			while (ptr > -1) {
				if (idx[i][ptr] == 0) {
					// copy idx[]=0 to idx[]=1
					for (k = n - 1; k > ptr; k--) {
						idx[i][k] = 0;
						data[i][k] = data[i][ptr];
					}
					break;
				} ptr--; // move until find  idx[]=1
			}
		}
		
		// check 1000 in middle	
		ptr = 0; ptr1 = 0;
		while (ptr < n) {
			if (idx[i][ptr] == 0) {
				init = data[i][ptr];
			}
			else {
				ptr1 = ptr;
				while (1) {
					if (idx[i][ptr1] == 0) {
						fin = data[i][ptr1];
						avg = (init + fin) / 2;
						// cout << endl << "init = " << init << " : final = " << fin << " : avg = " << avg << endl;
						break;
					}
					ptr1++;
				}
				ptr1 = ptr;
				while (1) {
					if (idx[i][ptr1] == 1) {
						idx[i][ptr1] = 0;
						data[i][ptr1] = avg;
					}
					else { break; }
					ptr1++;
				}
			}
			++ptr;
		}

		for (j = 0; j < n; j++) { data_new[x][j]=data[i][j]; }
		for (j = 0; j < 100; j++) { name_new[x][j] = name[i][j]; }
		x++;
	}
	print_gene_expression(name_new, data_new, x);

	printf("\n %s %d %s %d\n", "saving data = ",x,"; removing data (missing data) = ",sum_count);
	printf(" %s \n %s \n","File gene_expression_no1.txt already saved", "Press enter to exit");
	getchar(); 

	free(data);
	free(data_new);
	return 0;
}