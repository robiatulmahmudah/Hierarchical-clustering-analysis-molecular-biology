#include<stdio.h>
#include <stdlib.h>
#include <math.h> 

//attribute number 
#define x 18
//gene number 
#define y 6019  // 6019 to 5
//threshold of correlation coefficient between gene pair 
#define corr_threshold 0.9 // 0.9 to 0.7
//number gene_expression
#define n_gene 2732 // 2732 to 3

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


int main(void) {
	FILE *in,*out, *in1;
	int i, j, k[y];
	gene_corr *gene_corr_array; 
	double correlation, **data; 

	printf("\n %s\n %s\n", "Make sure that file gene_expression_no1.txt and cluster_genes_no3 exist", "Press enter to continue");
	getchar();
	printf("%s", ".....wait....");

	// cluster_genes_no3.txt --> dummy_no3.txt consist of 5 gene
	if ((in = fopen("cluster_genes_no3.txt", "r")) == NULL) {
		printf("can not open in_file");
		exit(1);
	}

	// gene_expression_no1.txt --> name_dummy_no3.txt"
	if ((in1 = fopen("gene_expression_no1.txt", "r")) == NULL) {
		printf("can not open in_file");
		exit(1);
	}

	//dynamic storage allocation for array of structure type 
	gene_corr_array = (gene_corr *)malloc(sizeof(gene_corr)*(n_gene)); 

	//store info  line number gen a, gen b, corr_coef into gene_corr_array
	for (i = 0; i<n_gene; i++) {
		fscanf(in, "%d %d %lf\n", &gene_corr_array[i].gene_num1, &gene_corr_array[i].gene_num2, &gene_corr_array [i].corr); 
	}

	//dynamic storage allocation of array "data" (for one dimension of array "data")
	data = (double **)malloc(sizeof(double *)*y);
	for (i = 0; i<y; i++) {

		//dynamic storage allocation of array "data" (for two dimension of array "data")
		data[i] = (double *)malloc(sizeof(double)*x);

		//store gene names in array "name" (e.g., YDL185W is gene name)
		fscanf(in1, "%s", &name[i]);
		for (j = 0; j<x; j++) {

			//load gene functionnal expression data on array "data"
			fscanf(in1, "%lf", &data[i][j]);
		}
	}


	//implement sample program 4
	int c1,c2,**cluster, num_genes;
	cluster = ( int **)malloc(sizeof(int *)*y); 
	for( i =0 ; i<y; i++){ 
		cluster[i] = ( int *)malloc(sizeof(int)*y); 
		for(j=0; j<y; j++){ 
			if(i==j){ 
				cluster[i][j] = 1;
			} 
			else{ 
				cluster[i][j] = 0; 
			}
		}
	} 

	/** //interupt for checking
	int rows = sizeof(cluster)/sizeof(cluster[0]);
    int cols = sizeof(cluster[0])/sizeof(cluster[0][0]);
    printf("%s %d \t %d \n", "size of CLUSTER in rows and cols are ",rows,cols);
	printf( "sizeof(**cluster): %lu\n", sizeof(**cluster) );
	**/

	//clustering gene into matrix 0,1	
	for(i=0;i<n_gene;i++){
		//printf("%s %d %s %d\n","awal i = ",i," dari n_gene = ", n_gene );
		//printf("%s %f %s %f\n","  gene_corr_array[i].corr = ",gene_corr_array[i].corr,"corr_threshold = ",corr_threshold);

		if(gene_corr_array[i].corr < corr_threshold){ 
			break;
		}
		else{
			c1=c2=0;
			//printf("%s %d %s %d\n","    sblm  while    c1 = ",c1,"; c2 = ",c2); 
			while(cluster[c1][gene_corr_array[i].gene_num1] != 1) c1++; 
			while(cluster[c2][gene_corr_array[i].gene_num2] != 1) c2++; 
			//printf("%s %d %s %d\n","    ssdah while    c1 = ",c1,"; c2 = ",c2); 
			if (c1 != c2){ 
				for(j=0; j<y; j++){ 
					//printf("%s %d %s %d\n","      dalam loop c1 = ",c1,"; c2 = ",c2); 
					cluster[c1][j]=cluster[c2][j] + cluster[c1][j]; 
					cluster[c2][j] = 0; 
				}
			}
		//printf("%s %d %s %d\n","akhir i = ",i," dari n_gene-1 = ", n_gene-1 );
		}
	} 
	printf("\n %s\n", "Cluster already built");

	//Interprate matrix 0,1 into name
	for (i = 0; i<y; i++) {
		//printf("%s \n", name[i]);
	}
	for (i = 0; i<y; i++) {
		k[i]=0;
		for(j=0;j<y;j++){
			if(cluster[i][j]==1) k[i]++;
		}
		//printf("%s %d %s %d\n", "jumlah baris ke - ",i, " adalah ",k[i]);
	}

	int longest = 0;
    for(i=0;i<y;i++)
    {
        if(k[i]>longest)
        longest=k[i];
    }
    printf("%s %d\n","The great number of genes in a cluster is : ",longest);



	//print final condition
	//printf("\n %s\n", "gene_corr_array already sorted");
	//printf("%s %d %d %lf\n","this is 1st data : ", gene_corr_array[0].gene_numl, gene_corr_array[0].gene_num2, gene_corr_array[0].corr);
	//printf("%s %d %d %lf\n","this is last data : ", gene_corr_array[n_gene-1].gene_numl, gene_corr_array[n_gene-1].gene_num2, gene_corr_array[n_gene-1].corr);


	// saving and exit
	//cluster code in 0,1 saved [optional]
	/**
	if((out=fopen("cluster_genes_no4.txt", "w"))==NULL){ 
		printf("can not open out_file\n"); 
		exit(1); 
	}
	//write content of cluster in a new file  
	for( i =0 ; i<y; i++){ 
		for(j=0; j<y; j++){ 
			fprintf(out, "%d\t", cluster[i][j]);
		}
		fprintf(out, "\n");
	} 
	**/
	// blm ktm logic unt saving nama gene harusnya berkorelasi dengan cluster_genes_no3 dikaitkan dengan cluster[?][y] dan name[y][100]
	// 1. cluster_genes_no3 -> gene_corr_array[n_gene].gene_num1, &gene_corr_array[n_gene].gene_num2, &gene_corr_array [n_gene].corr
	// 2. gene_corr_array[n_gene].gene_num1, &gene_corr_array[n_gene].gene_num2 BERKORELASI dengan name[y][100], n_gene subset dari y
	// 3. gene_corr_array[n_gene].gene_num1, &gene_corr_array[n_gene].gene_num2 --> cluster[?][y]
	// goal : file[y][name[y][100] pada peta cluster[?][y] ]
	//cluster name of gene saved
	if((out=fopen("name_cluster_genes_no4.txt", "w"))==NULL){ 
		printf("can not open out_file\n"); 
		exit(1); 
	}

	//write content of cluster in a new file 
	for( i =0 ; i<y; i++){ 
		c1=0;
		//printf("%s %d\n","i = ",i);
		while(cluster[c1][i] != 1) c1++; 
		fprintf(out, "%s\t", name[i]);
	    for(j=0;j<y;j++){
	    	if((cluster[c1][j] == 1)&&(i!=j))fprintf(out, "%s\t", name[j]);
	    	//printf("%s %d %s %d\n", "i = ",i,"; j = ",j);
	   	}
	    fprintf(out, "\n");
  	} 

	printf("\n %s\n %s\n", "File cluster_genes_no4.txt and name_cluster_genes_no4.txt already saved", "Press enter to exit");
	getchar();

	free(data);
	free(gene_corr_array);
	free(cluster);
	return 0;



}





