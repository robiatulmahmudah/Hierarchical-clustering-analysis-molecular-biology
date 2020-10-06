#include <stdio.h>
#include <stdlib.h>
#include <math.h> 

//attribute number 
#define x 18
//gene number 
#define y 6019 
//threshold of correlation coefficient between gene pair 
#define corr_threshold 0.9 
//number gene_expression
#define n_gene 2732

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

//quick sort algorithm
void quick(gene_corr *gene_corr_array,int start,int end) { 
	int pl=start; 
	int pr=end; 
	double center=gene_corr_array[(pl+pr)/2].corr; 
	gene_corr t;

	do{ 
		while(gene_corr_array[pl].corr > center) pl++; 
		while(gene_corr_array[pr].corr < center) pr--; 
		if(pl<=pr){ 
			t=gene_corr_array[pl]; 
			gene_corr_array[pl]=gene_corr_array[pr]; 
			gene_corr_array[pr]=t;

			pl++; 
			pr--; 
		}
	}while(pl<=pr); 

	if(start < pr) quick(gene_corr_array, start, pr); 
	if(pl < end) quick(gene_corr_array, pl, end); 
} 


int main(void) {
	FILE *in, *out;
	int i; 
	gene_corr *gene_corr_array; 
	double correlation; 

	printf("\n %s\n %s\n", "Make sure that file cluster_genes_no2.txt exist", "Press enter to continue");
	getchar();

	if ((in = fopen("cluster_genes_no2.txt", "r")) == NULL) {
		printf("can not open in_file");
		exit(1);
	}

	//dynamic storage allocation for array of structure type 
	gene_corr_array = (gene_corr *)malloc(sizeof(gene_corr)*(n_gene)); 

	//store info  line number gen a, gen b, corr_coef into gene_corr_array
	for (i = 0; i<n_gene; i++) {
		fscanf(in, "%d %d %lf\n", &gene_corr_array[i].gene_num1, &gene_corr_array[i].gene_num2, &gene_corr_array [i].corr); 
	}


	//print initial condition
	printf("%s\n", "data success put into gene_corr_array");
	printf("%s %d %d %lf\n","this is 1st data : ", gene_corr_array[0].gene_num1, gene_corr_array[0].gene_num2, gene_corr_array[0].corr);
	printf("%s %d %d %lf\n\n","this is last data : ", gene_corr_array[n_gene-1].gene_num1, gene_corr_array[n_gene-1].gene_num2, gene_corr_array[n_gene-1].corr);

	quick(gene_corr_array,0,n_gene-1);

	//print final condition
	printf("%s\n", "gene_corr_array already sorted");
	printf("%s %d %d %lf\n","this is 1st data : ", gene_corr_array[0].gene_num1, gene_corr_array[0].gene_num2, gene_corr_array[0].corr);
	printf("%s %d %d %lf\n\n","this is last data : ", gene_corr_array[n_gene-1].gene_num1, gene_corr_array[n_gene-1].gene_num2, gene_corr_array[n_gene-1].corr);

	//saved sorted gene_corr_array into file
	if ((out = fopen("cluster_genes_no3.txt", "w")) == NULL) {
		printf("can not open out_file\n");
		exit(1);
	}
	//insert line numbers of gene pairs and their correlation coefficient into the new file 
	for (i = 0; i<n_gene; i++) {
		fprintf(out, "%d\t%d\t%lf\n", gene_corr_array[i].gene_num1, gene_corr_array[i].gene_num2, gene_corr_array[i].corr);
	}


	printf(" %s\n %s\n", "File cluster_genes_no3.txt already saved", "Press enter to exit");
	getchar();

	free(gene_corr_array);
	return 0;
}

