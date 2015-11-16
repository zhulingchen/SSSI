
//Schlumberger 2012 - written by Andriy Gelman (AGelman@slb.com)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mex.h>
#include "matrix_toolbox_ag.h"

FILE *f1;

int L2_regression(double *output, double *input, double *h, double *kernel_size_mat, int rows, int cols, double *input_map, double *C00, double *C01, double *C11, int ext, double *output_label);


int find_minimum(double *input, double *input_minimum_value, int r, int c);
double abs_value(double input);

//int identity(double *input, int r, int c);
int create_exponential_kernel(double *output, double *C00, double *C01, double *C11, double h, double *input_map, int rows, int cols, int i, int j, int half_ksize);
double square_val(double input);
int find_maximum_vector(double *minimum, double *input, int samples);

int create_reweighted_kernel(double *W_reweighted, double *z, int number_of_samples, double *X, double *s);

int copy_data(double *output, double *input, int r, int c);

//int absolute_value(int *input);
int absolute_value(double *input);
//int polar_ribiere_optimization(double *output, double *Xw, double *X, double *z, double *s_init, int number_of_samples );


//Initialize mex functions
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 
    mxArray *output_mxArray, *output_label_mxArray;

    const mwSize *dims;
    
    double *input, *C00, *C01, *C11;
    double *output, *output_label;
    double *h;
    double *input_map;
    
    //optimization parameter defines whether L1 or L2 denoising is performed. To select L1 input 1, L2 = 2. L1 minimization requires an initialization value. 
    //int kernel_size;
    double *kernel_size_mat;
   
    int rows, cols, ext;
    
    dims = mxGetDimensions(prhs[0]);
    rows = (int)dims[0];
    cols = (int)dims[1];
    
    //Initialize input into pointer
    input = mxGetPr(prhs[0]);
    
    //say 7 is the map which specifies which are pixels are noisy and which are not. 
    input_map = mxGetPr(prhs[1]);
    
    //get h
    h = mxGetPr(prhs[2]);
    
    //kernel_size = (int)(*mxGetPr(prhs[9]));
    kernel_size_mat = mxGetPr(prhs[3]);
    
    C00 = mxGetPr(prhs[4]);
	C01 = mxGetPr(prhs[5]);
	C11 = mxGetPr(prhs[6]);
    
	ext = (int)(*mxGetPr(prhs[7]));
    

    //Initialize output_mxArray and output;  
    output_mxArray = plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    output_label_mxArray = plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    
    output = mxGetPr(output_mxArray);
    output_label = mxGetPr(output_label_mxArray);
    
    //copy the initialised data into output
    copy_data(output, input, rows, cols);
    
	L2_regression(output, input, h, kernel_size_mat, rows, cols, input_map, C00, C01, C11, ext, output_label); 

}


int L2_regression(double *output, double *input, double *h, double *kernel_size_mat, int rows, int cols, double *input_map, double *C00, double *C01, double *C11, int ext, double *output_label)
{
    
    
    //loop variables
    int i,j, k=0;
    int ii;
    int k_W = 0;
    
    //loop variables in the window
    int i_win, j_win;
    
	//int ind_i, ind_j;

    //int number_of_samples = kernel_size*kernel_size;
    int number_of_samples;
    
    int half_kernel_size;
    int original_half_kernel_size;
    
    double *W;
    double *W_kernel;
    double *z;
    
    //define covariance matrix array
    double cov_mat[4] = {1, 0, 0, 1};
    double *x_vec;
    double *y_vec;
    double h_val;
    double threshold = 1.0e-5;
    double weight_threshold;

    int y_length;
    int x_length;
    int stop_here;
    int max_number_of_samples = (int)square_val(61);
    int min_number_of_samples;
    double sum_weight;

    W = ( double* )malloc( max_number_of_samples*sizeof( double) );
    W_kernel = ( double* )malloc( max_number_of_samples*sizeof( double) );
    z = ( double* )malloc( max_number_of_samples*sizeof( double) );
    
    
    x_vec = ( double* )malloc( max_number_of_samples*sizeof( double) );
    y_vec = ( double* )malloc( max_number_of_samples*sizeof( double) );
    
    
    //matrix_init_ones(W, number_of_samples, 1);
    //evaluate the weighting matrix based on Gaussian Kernel
    //loop through each pixel
    //matrix_initialize_vector(x_vec, -half_kernel_size, half_kernel_size);
    
    //scan ordering of pixels
    //matrix_initialize_vector(y_vec, half_kernel_size, -half_kernel_size);
    
	//f1 = fopen("z.txt", "w");
    
	//value ext specifies the extensions which are applied to the input image. no filtering is applied in the extended 
	
    //initialize the output_label
    //this matrix contains values which are deemed to be unreliable to the selection of h parameter
    matrix_init_ones(output_label, rows, cols);
    
    for(i=ext; i< (rows-ext); i++){
        for(j=ext; j<(cols-ext); j++){
    
            
            
            
	//		mexPrintf("%d %d.\n", 
    //             i, j);

			//fprintf(f1, "%d, %d \n", i, j);
	
			
			original_half_kernel_size = (int)(kernel_size_mat[i + j*rows]-1)/2;
    
			if ((i ==  249) & (j == 50))
                 stop_here = 1;

            //if ( (h[i+j*rows] > 0.8) && (input[i+j*rows] >= 0) ){
            if ( (h[i+j*rows] > 0.2)){    						
				//mexPrintf("%d %d.\n", i, j);

                h_val = h[i+j*rows];
				
                half_kernel_size = original_half_kernel_size;
                
				min_number_of_samples = 5;
				number_of_samples = 0;

                                          
				create_exponential_kernel(W_kernel, C00, C01, C11, h_val, input_map, rows, cols, i, j, half_kernel_size);

                //create z vector
                k_W = 0;
                k = 0;
                
                //evaluate number_of_samples
                //also find minimum distance to the nearest point

                sum_weight = 0;
                for(i_win= i-half_kernel_size; i_win <= i+half_kernel_size; i_win++){
					if((i_win >= 0) && (i_win < rows)){
						for(j_win= j-half_kernel_size; j_win <= j+half_kernel_size; j_win++){
							if((j_win >= 0) && (j_win < cols)){
								if (input_map[i_win+j_win*rows] == 1){
									if (W_kernel[k_W] > threshold){
										z[k] = input[i_win+j_win*rows];
										W[k] = W_kernel[k_W];
										sum_weight = sum_weight + W[k];
                                        k++;
									}
									k_W++;
								}
							}
						}
					}
                }
                
                number_of_samples = k;
                
                if (number_of_samples < min_number_of_samples)
                    output_label[i+j*rows] = 0;
                
                output[i+j*rows] = 0;
                for(ii = 0; ii < number_of_samples; ii++){
                    output[i+j*rows] = output[i+j*rows] + W[ii]*z[ii];
                }
                output[i+j*rows] = output[i+j*rows]/(sum_weight);
                
            }
            else{
                output[i+j*rows] = input[i+j*rows];
                
            }
        }
    }
    
	free(x_vec);
    free(y_vec);
    free(W);
    free(W_kernel);
    free(z);
    
    return 0;
    
}


int find_minimum(double *input, double *input_minimum_value, int r, int c)
{
 
int i, j;

double minimum_value = 1e10;

for(i=0;i<r;i++){
	for(j=0;j<c;j++){
		if (input[i+j*r] < minimum_value){
			minimum_value = input[i+j*r]; 
		}
	}
}

*input_minimum_value = minimum_value;

return 0;

}


double abs_value(double input){

	if (input > 0)
		return input;
	if (input < 0)
		input = -input;
	return input;
}

double square_val(double input){
    return input*input;
}


int create_exponential_kernel(double *output, double *C00, double *C01, double *C11, double h, double *input_map, int rows, int cols, int i, int j, int half_ksize)
{
	int ii,jj,k = 0;
    
    double cov_mat[4];
	double cov_det;


	for(ii= -half_ksize; ii <= half_ksize; ii++){          
		for(jj = -half_ksize; jj <= half_ksize; jj++){
			//check the label to see, whether this pixel is located in the label
			if (((i + ii) >= 0) && ((i + ii) < rows) && ((j+jj) >= 0) && ((j+jj) < cols)){
                
                cov_mat[0] = C00[(i+ii) +(j+jj)*rows];
                cov_mat[1] = C01[(i+ii) +(j+jj)*rows];
                cov_mat[2] = C01[(i+ii) +(j+jj)*rows];
                cov_mat[3] = C11[(i+ii) +(j+jj)*rows];
                
                cov_det = cov_mat[0]*cov_mat[3] - cov_mat[1]*cov_mat[2];
                
                if (input_map[ (i+ii) +(j+jj)*rows] == 1){
                    output[k] = 100*(1/(2*3.142*square_val(h)))* pow(cov_det, -0.5) * exp(-(1/cov_det)*(square_val(jj)*cov_mat[3] - 2*(ii)*(jj)*cov_mat[1] + square_val(ii)*cov_mat[0])/(2*square_val(h)));
                    //output[k] = exp(-(1/cov_det)*(square_val(jj)*cov_mat[3] - 2*(ii)*(-jj)*cov_mat[1] + square_val(ii)*cov_mat[0])/(2*square_val(h)));
                    k++;
                }
			}
		}        
	}
}

int create_reweighted_kernel(double *W_reweighted, double *z, int number_of_samples, double *X, double *s){

    double epsilon = 1.0e-2;
    int i;
    //first evaluate z - Xs
    //use W_reweighted as a buffer
    matrix_multiplication(W_reweighted, X, s, number_of_samples, 6, 6, 1);
    matrix_matrix_multiply_by_scalar(W_reweighted, W_reweighted, number_of_samples, 1, -1);          
    matrix_add_matrix(W_reweighted, z, W_reweighted, number_of_samples, 1);
    matrix_matrix_abs_output(W_reweighted, W_reweighted, number_of_samples, 1);
        
    for(i = 0; i<number_of_samples; i++){
        //add epsilon to each point in the matrix - this ensures that the inverse is stable
        W_reweighted[i] = W_reweighted[i] + epsilon; 
    }
    
    //invert each value
    for(i = 0; i<number_of_samples; i++){
        //add epsilon to each point in the matrix - this ensures that the inverse is stable
        W_reweighted[i] = 1/W_reweighted[i]; 
    }
    return 0;
    
}

int find_maximum_vector(double *maximum, double *input, int samples){

	int i;
	*maximum = 0;
	for (i = 1; i<samples; i++){
		if(*maximum < (input[i]-input[i-1]))
			*maximum = (input[i]-input[i-1]);
	}
	return 0;
}

int copy_data(double *output, double *input, int r, int c){
    
    int i,j;
    for(i=0;i<r;i++){
        for(j=0;j<c;j++){
            output[i+j*r] = input[i+j*r];
        }
    }
}


int absolute_value(double *input){

	if (*input < 0)
		*input = -*input;

}