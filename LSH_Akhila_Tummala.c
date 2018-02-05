#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<limits.h>
#include<string.h>

#define ROWS 25000 //No of points retrived from the dataset
#define DIM 74 //No of dimensions for each point
#define file_loc "bio_train_bin.bin" //Location of the bin file of the dataset
#define file_query "query_text.txt" //Location of the query file

void read_binary(char *file_path, double *data);
double* read_query(char *file_path);
void lsh_search(int dim, int ndata, double *data, int m, double **r, double *b, double w, int *cluster_size, int *cluster_start, int **h, int *cluster_assign, double *query);
void lsh(int dim, int ndata, double *data, int m, double **r, double *b, double w, int **h, int *cluster_assign);

int num_clusters = 0; //No of clusters
void main(){
	int i = 0, j = 0;
	int dim = DIM, ndata = ROWS, m = 37;
	double w = 100000;
	double *data; //Stores the data from the dataset
	data = (double*) malloc((dim*ndata) * sizeof(double));
	read_binary(file_loc, data); //reads the binary file
	double **r;
	r = (double**) malloc(m * sizeof(double*));
	for(i=0;i<m;i++)
		r[i] = (double*) malloc(dim * sizeof(double));
	srand(time(NULL));
	for(i=0;i<m;i++)
		for(j=0;j<dim;j++)
			r[i][j] = rand() % 101; //randomly generates r values
	double *b;
	b = (double*) malloc(m * sizeof(double));
	srand(time(NULL));
	for(i=0;i<m;i++)
		b[i] = 0;
	int **h;
	h = (int**) malloc(ndata * sizeof(int*));
	for(i=0;i<ndata;i++)
		h[i] = (int*) malloc(m * sizeof(int));
	for(i=0;i<ndata;i++)
		for(j=0;j<m;j++)
			h[i][j] = 0;
	int *cluster_assign;
	cluster_assign = (int*) malloc(ndata * sizeof(int));
	for(i=0;i<ndata;i++)
		cluster_assign[i]=0; //used to store the cluster assign values
	lsh(dim, ndata, data, m, r, b, w, h, cluster_assign); //call the lsh function
	int *cluster_size; //stores the sizes of the clusters
	cluster_size = (int*) malloc(num_clusters * sizeof(int));
	int *cluster_start; //stores the starting indexes of the clusters
	cluster_start = (int*) malloc(num_clusters * sizeof(int));
	for(i=0;i<num_clusters;i++){
		cluster_size[i]=0;
		cluster_start[i]=0;
	}
	
	for(i=1;i<=num_clusters;i++){
		int temp_size = 0;
		for(j=0;j<ndata;j++){
			if(cluster_assign[j]==i)
				temp_size++;
		}
		cluster_size[i-1] = temp_size;
	}
	
	cluster_start[0] = 0;
	for(i=1;i<ndata;i++)
		cluster_start[i] = cluster_start[i-1] + cluster_size[i-1]*dim;
	printf("Total Clusters = %d\n", num_clusters);
	double *query_total;
	query_total = (double*) malloc((dim*10) * sizeof(double)); //stores 10 query points
	query_total = read_query(file_query); //read the query file
	double *query = (double*) malloc(dim * sizeof(double));
	for(i=0;i<10;i++){
		for(j=0;j<dim;j++)
			query[j] = query_total[i*dim+j]; //pass each query point to the search algorithm
		lsh_search(dim, ndata, data, m, r, b, w, cluster_size, cluster_start, h, cluster_assign, query);
	}
}
//<-------- Reads the binary file ---------->
void read_binary(char *file_path, double *data){
	FILE *file = fopen(file_path, "rb");
	int i;
	double temp_data = -1.0;
	for(i=0;i<ROWS*DIM;i++){
		fread(&temp_data, sizeof(double), 1, file);
		data[i] = temp_data;
	}
	fclose(file);
}
//<---------- Reads the query file -------------->
double* read_query(char* file_path){
	FILE *file = fopen(file_path, "r");
	int i;
	double *data = (double*) malloc(10*DIM * sizeof(double));
	char *split_data = malloc(256 * sizeof(char));
	for(i=0;i<10*DIM;i++){
		fscanf(file, "%s", split_data);
		data[i] = strtod(split_data, NULL);
	}
	fclose(file);
	return data;
}

void lsh(int dim, int ndata, double *data, int m, double **r, double *b, double w, int **h, int *cluster_assign){
	int i=0, j=0, k=0;
	int h_compare[ndata][m]; //stores the h values assigned until now for dividing clusters
	for(i=0;i<ndata;i++)
		for(j=0;j<m;j++)
			h_compare[i][j]=0;
	
	//<--------- Unit vectors or r ------->
	for(i=0;i<m;i++){
		double temp_sum = 0.0;
		for(j=0;j<dim;j++)
			temp_sum = temp_sum + pow(r[i][j],2);
		temp_sum = sqrt(temp_sum);
		for(j=0;j<dim;j++)
			r[i][j] = r[i][j] / temp_sum;
	}
	
	//<------------- Calculate the h values -------->
	for(i=0;i<ndata;i++){
		for(j=0;j<m;j++){
			double temp_h = 0.0;
			for(k=0;k<dim;k++)
				temp_h = temp_h + data[i*dim+k] * r[j][k];
			temp_h = floor((temp_h - b[j]) / w);
			h[i][j] = (int) temp_h;
		}
	}
	
	//<------- Fill the cluster assigns ------>
	for(i=0;i<m;i++)
		h_compare[0][i] = h[0][i];
	cluster_assign[0] = 1;
	num_clusters = 1;
	
	for(i=1;i<ndata;i++){
		int m_matches = 0, data_point=-1;
		for(j=0;j<i;j++){
			m_matches = 0;
			for(k=0;k<m;k++){
				if(h[i][k]==h_compare[j][k])
					m_matches++;
			}
			if(m_matches==m){
				data_point = j;
				break;
			}
		}
		if(data_point==-1){
			num_clusters++;
			cluster_assign[i] = num_clusters;
		}
		else
			cluster_assign[i] = cluster_assign[data_point];
	}
	
	//<--------- Sort the cluster assign, h values and the data array ------>
	int cluster_counter=1;
	int sorted_index = 0;
	int temp_assign = 0;
	double temp_sort_data = 0.0;
	int temp_sort_h = 0;
	while(cluster_counter!=num_clusters){
		for(i=sorted_index;i<ndata;i++){
			if(cluster_assign[i]==cluster_counter){
				temp_assign = cluster_assign[sorted_index];
				cluster_assign[sorted_index] = cluster_assign[i];
				cluster_assign[i] = temp_assign;
				for(j=0;j<dim;j++){
					temp_sort_data = data[i*dim+j];
					data[i*dim+j] = data[sorted_index*dim+j];
					data[sorted_index*dim+j] = data[i*dim+j];
				}
				for(j=0;j<m;j++){
					temp_sort_h = h[sorted_index][j];
					h[sorted_index][j] = h[i][j];
					h[i][j] = temp_sort_h;
				}
				sorted_index++;
			}
		}
		cluster_counter++;
	}
	//<---------- End of sorting --------------->
}

void lsh_search(int dim, int ndata, double *data, int m, double **r, double *b, double w, int *cluster_size, int *cluster_start, int **h, int *cluster_assign, double *query){
	int i=0, j=0;
	int h_query[m]; //stores the h values for query
	for(i=0;i<m;i++)
		h_query[i] = 0;
	
	//<------- calculate the h values for query point ------->
	for(i=0;i<m;i++){
		double temp_h = 0.0;
		for(j=0;j<dim;j++)
			temp_h = temp_h + query[j] * r[i][j];
		temp_h = floor((temp_h - b[i]) / w);
		h_query[i] = (int) temp_h;
	}
	
	//<---------- Compare the query h values with actual h values ------->
	int search_cluster = -1;
	for(i=0;i<num_clusters;i++){
		int m_matches = 0;
		for(j=0;j<m;j++)
			if(h[cluster_start[i]/dim][j]==h_query[j])
				m_matches++;
		if(m_matches==m){
			search_cluster = i;
			break;
		}
	}
	
	if(search_cluster==-1)
		printf("No cluster found for query point\n");
	else{
		//<---------- Exhaustive search in the cluster ------->
		double min_dist = LONG_MAX;
		for(i=0;i<cluster_size[search_cluster];i++){
			double temp_dist = 0.0;
			for(j=0;j<dim;j++)
				temp_dist = temp_dist + pow((query[j] - data[cluster_start[search_cluster]+i*dim+j]),2);
			temp_dist = sqrt(temp_dist);
			if(temp_dist<min_dist)
				min_dist = temp_dist;
		}	
		printf("Total visits = %d\n",cluster_size[search_cluster]);
		printf("Minimum distance = %f\n", min_dist);
	}
}
