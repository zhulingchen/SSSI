
//Schlumberger 2012 - written by Andriy Gelman (AGelman@slb.com)



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mex.h>
#include "matrix_toolbox_ag.h"

FILE *f1;


int	estimate_steering_parameters(double *C00, double *C01, double *C11, double *sigma_1_mat, double *sigma_2_mat, double *theta_mat, double *input_y, double *input_x, double *input_map, double *pixels_to_be_estimated, double min_sigma, int win_size, int rows, int cols);
    
int find_minimum(double *input, double *input_minimum_value, int r, int c);
double abs_value(double input);

//int identity(double *input, int r, int c);
//int create_exponential_kernel(double *output, double *cov_mat, double h, double *input_map, int rows, int cols, int i, int j, int half_ksize);
double square_val(double input);
int find_maximum_vector(double *minimum, double *input, int samples);

//int polar_ribiere_optimization(double *output, double *Xw, double *X, double *z, double *s_init, int number_of_samples );


//Initialize mex functions
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 
    mxArray *C00_mxArray, *C01_mxArray, *C11_mxArray, *sigma_1_mxArray, *sigma_2_mxArray, *theta_mxArray; 

    const mwSize *dims;
    
    double *input_y, *input_x, *input_map, *pixels_to_be_estimated, *win_size, *min_sigma;
    double *C00, *C01, *C11, *sigma_1_mat, *sigma_2_mat, *theta_mat;
    
    //optimization parameter defines whether L1 or L2 denoising is performed. To select L1 input 1, L2 = 2. L1 minimization requires an initialization value. 
    int rows, cols;
    
    dims = mxGetDimensions(prhs[0]);
    rows = (int)dims[0];
    cols = (int)dims[1];
    
    //Initialize input into pointer
	input_y = mxGetPr(prhs[0]);          
    input_x = mxGetPr(prhs[1]);
    input_map = mxGetPr(prhs[2]);
    pixels_to_be_estimated = mxGetPr(prhs[3]);
    win_size = mxGetPr(prhs[4]);
    min_sigma = mxGetPr(prhs[5]);
    
    //Initialize output_mxArray and output;  
    C00_mxArray = plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    C01_mxArray = plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    C11_mxArray = plhs[2] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    sigma_1_mxArray = plhs[3] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    sigma_2_mxArray = plhs[4] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    theta_mxArray = plhs[5] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    
    
    C00 = mxGetPr(C00_mxArray);
    C01 = mxGetPr(C01_mxArray); 
    C11 = mxGetPr(C11_mxArray);
    sigma_1_mat = mxGetPr(sigma_1_mxArray);
    sigma_2_mat = mxGetPr(sigma_2_mxArray);
    theta_mat = mxGetPr(theta_mxArray);
            
	estimate_steering_parameters(C00,  C01, C11, sigma_1_mat, sigma_2_mat, theta_mat, input_y, input_x, input_map, pixels_to_be_estimated, *min_sigma, (int)*win_size, rows, cols); 

}

int	estimate_steering_parameters(double *C00, double *C01, double *C11, double *sigma_1_mat, double *sigma_2_mat, double *theta_mat, double *input_y, double *input_x, double *input_map, double *pixels_to_be_estimated, double min_sigma, int win_size, int rows, int cols){
    
    
    //loop variables
    int i,j, k=0;
    
    //coefficients used to evaluate the eigenvalues of the matrix
    double a, b,c,d;
 
    //loop variables in the window
    int i_win, j_win;
    
    //int number_of_samples = kernel_size*kernel_size;
    int number_of_samples;
    
    int half_kernel_size;
    
    double *S;
    double *S_transpose;
    double A[4];
    
    double eig_1[2];
    double eig_2[2];
    double lam_1;
    double lam_2;
    double norm;
    double theta;
    double epsilon = 1.0;
    double epsilon_2 = 0.05;
    double sigma_1;
    double sigma_2;
    //double min_sigma = 0.7;
    //double min_sigma;
    double max_sigma;
    double R_theta[4];  
	double R_theta_transpose[4];
    double cov_mat[4];
    double lambda_mat[4];
    double buffer[4];
    double buffer_2[4];
	double mean_value;
    double gamma;
    
    int stop_here;
    int max_number_of_samples = (int)square_val(51);
    
    //evaluate scaling factor
    double scaling_factor;
    
    max_sigma = 1/min_sigma;
    
    S = ( double* )malloc( max_number_of_samples*2*sizeof( double) );
    S_transpose = ( double* )malloc( 2*max_number_of_samples*sizeof( double) );
    
    half_kernel_size = (win_size-1)/2;
    
    for(i=0; i< rows; i++){
        for(j=0; j<cols; j++){
            C00[i + j*rows] = 1;
            C01[i + j*rows] = 0;
            C11[i + j*rows] = 1;
            sigma_1_mat[i + j*rows] = 1;
            sigma_2_mat[i + j*rows] = 1; 
            theta_mat[i + j*rows] = 0;
        }
    }
    
    for(i=half_kernel_size; i< (rows-half_kernel_size); i++){
        for(j=half_kernel_size; j<(cols-half_kernel_size); j++){
            
            if (pixels_to_be_estimated[i+j*rows] == 1){
                
                if ((i ==  135) & (j == 177))
                    stop_here = 1;
                
                //create S matrix
                k = 0;
                
                for(i_win= i-half_kernel_size; i_win <= i+half_kernel_size; i_win++){
                    for(j_win= j-half_kernel_size; j_win <= j+half_kernel_size; j_win++){
                        if (input_map[i_win+j_win*rows] == 1){
                            if (abs_value(input_x[i_win+j_win*rows]) + abs_value(input_y[i_win+j_win*rows]) > 2)
                            {
                                S[0 + k*2] = input_x[i_win+j_win*rows];
                                S[1 + k*2] = input_y[i_win+j_win*rows];
                                k++;
                            }
                        }
                    }
                }
                
                number_of_samples = k;
                
                //evaluate mean of S_x and S_y and then subtract this value from each row
                /*
                mean_value = 0;
                for(k = 0; k <number_of_samples; k++){
                    mean_value = mean_value + S[0 + k*2];
                }
                mean_value = mean_value/number_of_samples;
                for(k = 0; k <number_of_samples; k++){
                    S[0 + k*2] = S[0 + k*2] - mean_value;
                }
                
                mean_value = 0;
                for(k = 0; k <number_of_samples; k++){
                    mean_value = mean_value + S[1 + k*2];
                }
                mean_value = mean_value/number_of_samples;
                for(k = 0; k <number_of_samples; k++){
                    S[1 + k*2] = S[1 + k*2] - mean_value;
                }
                */
                
                matrix_transpose(S_transpose, S, 2, number_of_samples);
                matrix_multiplication(A, S, S_transpose, 2, number_of_samples, number_of_samples, 2);
                
                //evaluate the eigenvalues and eigenvectors of the matrix A
                a = A[0];
                b = A[1];
                c = A[3];
                
                
                //if the determinant is greater than zero
                if ((a*c - square_val(b)) > 0) {
                    d = sqrt(square_val(a+c) - 4*(a*c - square_val(b)));
                    
                    lam_1 = (a+c-d)/2;
                    lam_2 = (a+c+d)/2;
                    
                    if (abs(b) > 0){
                        eig_1[0] = -b;
                        eig_1[1] = (a-lam_1);
                    }
                    else
                    {
                        eig_1[0] = 1;
                        eig_1[1] = 0;
                    }
                        
                    norm = square_val(eig_1[0])+square_val(eig_1[1]);
                    eig_1[0] = eig_1[0]/sqrt(norm);
                    eig_1[1] = eig_1[1]/sqrt(norm);
                    
                    if (abs(b)>0){
                        eig_2[0] = -b;
                        eig_2[1] = (a-lam_2);
                    }
                    else
                    {
                        eig_2[0] = 0;
                        eig_2[0] = 1;
                    }
                        
                    
                    norm = square_val(eig_2[0])+square_val(eig_2[1]);
                    eig_2[0] = eig_2[0]/sqrt(norm);
                    eig_2[1] = eig_2[1]/sqrt(norm);
                    
                    theta = atan(eig_1[1]/eig_1[0]);
                    
                    sigma_2 = (lam_1 + epsilon)/(lam_2 + epsilon);
                    if (sigma_2 < min_sigma)
                        sigma_2 = min_sigma;
                    if (sigma_2 > max_sigma)
                        sigma_2 = max_sigma;
                    
                    sigma_1 = 1/sigma_2;
                    
                    //sigma_1_mat[i + j*rows] = sqrt(lam_1/number_of_samples);
                    //sigma_2_mat[i + j*rows] =  sqrt(lam_2/number_of_samples);

                    //scaling_factor = sqrt(((lam_1*lam_2 + 1)));
                    //scaling_factor = (1/scaling_factor)*100;
                   
                    //scaling_factor = (sigma_1_mat[i + j*rows] - sigma_2_mat[i + j*rows] + epsilon_2)/(sigma_1_mat[i + j*rows] + sigma_2_mat[i + j*rows] + epsilon_2);
                    
                    scaling_factor = 1;
                    //scaling_factor = sqrt(((lam_1*lam_2 + epsilon_2)/number_of_samples));

                    
                    lambda_mat[0] = sigma_1*scaling_factor;
                    lambda_mat[1] = 0;
                    lambda_mat[2] = 0;
                    lambda_mat[3] = sigma_2*scaling_factor;
                    
                    R_theta[0] = cos(theta);
                    R_theta[1] = sin(theta);
                    R_theta[2] = -sin(theta);
                    R_theta[3] = cos(theta);
                    
                    sigma_1_mat[i + j*rows] = sigma_1;
                    sigma_2_mat[i + j*rows] =  sigma_2;
                    theta_mat[i + j*rows] =  theta;
                    
                    //evaluate cov_mat = R_theta*lambda_mat*lambda_mat*R_theta';
                    matrix_multiplication(buffer, R_theta, lambda_mat, 2, 2, 2, 2);
                    buffer[0];
                    buffer[1];
                    buffer[2];
                    buffer[3];
                    
                    matrix_multiplication(buffer_2, buffer, lambda_mat, 2, 2, 2, 2);
                    buffer_2[0];
                    buffer_2[1];
                    buffer_2[2];
                    buffer_2[3];
                    
                    matrix_transpose(R_theta_transpose, R_theta, 2, 2);
                    R_theta_transpose[0];
                    R_theta_transpose[1];
                    R_theta_transpose[2];
                    R_theta_transpose[3];
                    
      
                    
                    
                    matrix_multiplication(cov_mat, buffer_2, R_theta_transpose, 2, 2, 2, 2);
                    cov_mat[0];
                    cov_mat[1];
                    cov_mat[2];
                    cov_mat[3];
                    
                    C00[i + j*rows] = cov_mat[0];
                    C01[i + j*rows] = cov_mat[1];
                    C11[i + j*rows] = cov_mat[3];
                    

                    
                }
                else{
                    C00[i + j*rows] = 1;
                    C01[i + j*rows] = 0;
                    C11[i + j*rows] = 1;
                }
                
                /*
                 * f1 = fopen("z.txt", "w");
                 * for(ii=0;ii<number_of_samples;ii++)
                 * {
                 * for(jj=0;jj<1;jj++)
                 * {
                 * fprintf(f1, "%lf ", z[ii+jj*number_of_samples]);
                 * }
                 * fprintf(f1, "\n");
                 * }
                 * fclose(f1);
                 */
            }
        }
    }
    
    free(S);
    free(S_transpose);

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


int create_exponential_kernel(double *output, double *cov_mat, double h, double *input_map, int rows, int cols, int i, int j, int half_ksize)
{
	int ii,jj,k = 0;

	double cov_det;

	cov_det = cov_mat[0]*cov_mat[3] - cov_mat[1]*cov_mat[2];

	for(ii= -half_ksize; ii <= half_ksize; ii++){          
		for(jj = -half_ksize; jj <= half_ksize; jj++){
			//check the label to see, whether this pixel is located in the label
			if (input_map[ (i+ii) +(j+jj)*rows] == 1){

				output[k] = exp(-(1/cov_det)*(square_val(jj)*cov_mat[3] - 2*(ii)*(jj)*cov_mat[1] + square_val(ii)*cov_mat[0])/(2*square_val(h)));
				//output[k] = exp(-(1/cov_det)*(square_val(jj)*cov_mat[3] - 2*(ii)*(-jj)*cov_mat[1] + square_val(ii)*cov_mat[0])/(2*square_val(h)));
				k++;
			}
		}        
	}
}

int create_reweighted_kernel(double *W_reweighted, double *z, int number_of_samples, double *X, double *s){

    double epsilon = 1.0e-8;
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


