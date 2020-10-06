#include<stdio.h>
#include <stdlib.h>
#include <math.h> 

//attribute number 
#define x 18
//gene number 
#define y 6019 
//threshold of correlation coefficient between gene pair 
#define corr_threshold 0.9 

char name[y][100]; 

//prototype declaration 
double corr_single(double a[],double b[],int n);
void insert_gene_corr_structure_array(double **data); 

//structure type to store line numbers of genes and correlation coefficient 
typedef struct GENE_CORRELATION{ 
	int gene_num1; 
	int gene_num2; 
	double corr; 
}gene_corr; 

//function to calculate correlation coefficient 
double corr_single(double a[],double b[],int n){ 
	int i; 
	double d=0.0; 
	double mean_a,mean_b; 
	double var_a,var_b; 
	double corr; 

	mean_a=0;mean_b=0;var_a=0;var_b=0;corr=0;
	for(i=0;i<n;i++) mean_a+=a[i]/n; 
	for(i=0;i<n;i++) mean_b+=b[i]/n; 
	for(i=0;i<n;i++) var_a+=(a[i]-mean_a)*(a[i]-mean_a)/n; 
	for(i=0;i<n;i++) var_b+=(b[i]-mean_b)*(b[i]-mean_b)/n; 
	for(i=0;i<n;i++) corr+=(a[i]-mean_a)/sqrt(var_a)*(b[i]-mean_b)/sqrt(var_b)/n; 

	return corr; 
} 

//function to calculate correlation coefficient between genes and store it in array of structure type 
void insert_gene_corr_structure_array(double **data){
	FILE *out; 
	int i, j; 
	gene_corr *gene_corr_array; 
	double correlation; 

	//file to insert results 
	if((out=fopen("cluster_genes_no2.txt", "w"))==NULL){ 
		printf("can not open out_file\n"); 
		exit(1); 
	}

	//dynamic storage allocation for array of structure type 
	gene_corr_array = (gene_corr *)malloc(sizeof(gene_corr)*(y*y)); 

	for(i=0; i<y; i++){ 
		for(j=0; j<y; j++){ 
			if(i<j){
				//calculate correlation coefficient between genes and store it in variable correlation 
				correlation = corr_single(data[i], data[j], x); 

				/*When correlation coefficient was over threshold, Line numbers of two genes and their 
				correlation coefficient are stored in their corresponding structure type.*/ 
				if(corr_threshold < correlation){ 
					gene_corr_array[i*y+j].corr=correlation; 
					gene_corr_array[i*y+j].gene_num1=i; 
					gene_corr_array[i*y+j].gene_num2=j; 
				} 

				/*When correlation coefficient was under threshold, the value of -5 is stored in variable corr of their 
				corresponding structure type.*/ 
				else{ 
					gene_corr_array[i*y+j].corr=-5; 
				} 
			} 

			//the value of -5 is stored in variable corr of over-prepared structure type 
			else{ 
				gene_corr_array[i*y+j].corr = -5; 
			}
		}
	} 

	//insert line numbers of gene pairs and their correlation coefficient into the new file 
	for(i=0; i<y; i++){ 
		for(j=0; j<y; j++){ 
			if(gene_corr_array[i*y+j].corr != -5){ 
				fprintf(out, "%d\t%d\t%lf\n", gene_corr_array[i*y+j].gene_num1, gene_corr_array[i*y+j].gene_num2, gene_corr_array [i*y+j].corr); 
			}
		} 
	} 

	//release dynamic storage allocation 
	free(gene_corr_array); 
}

int main(void) {
	FILE *in;
	int i, j; //, k, x, idx[y][n], ptr, ptr1, sum, sum_count;
	double **data; //, **data_new, init, fin, avg;

	printf("\n %s\n %s\n", "Make sure that file gene_expression_no1.txt exist", "Press enter to continue");
	getchar();

	if ((in = fopen("gene_expression_no1.txt", "r")) == NULL) {
		printf("can not open in_file");
		exit(1);
	}

	//dynamic storage allocation of array "data" (for one dimension of array "data")
	data = (double **)malloc(sizeof(double *)*y);
	for (i = 0; i<y; i++) {

		//dynamic storage allocation of array "data" (for two dimension of array "data")
		data[i] = (double *)malloc(sizeof(double)*x);

		//store gene names in array "name" (e.g., YDL185W is gene name)
		fscanf(in, "%s", &name[i]);
		for (j = 0; j<x; j++) {

			//load gene functionnal expression data on array "data"
			fscanf(in, "%lf", &data[i][j]);
		}
	}

	insert_gene_corr_structure_array(data);
	printf(" %s \n %s \n", "File cluster_genes_no2.txt already saved", "Press enter to exit");
	getchar();

	free(data);
	return 0;
}
