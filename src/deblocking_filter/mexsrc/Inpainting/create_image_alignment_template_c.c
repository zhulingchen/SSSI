
//Schlumberger 2012 - written by Andriy Gelman (AGelman@slb.com)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mex.h>
#include "matrix_toolbox_ag.h"

//FILE *f1;

int create_image_allignment_template(double *output_image, double *output_image_label, double *input_image, double *input_label, double *target_template_label, double *h, double *kernel_size_mat, double *C00, double *C01, double *C11, double *sigma_1, double *sigma_2, double *theta, int ext_length, int rows, int cols, double *A, double *phase_shift, double *dip_compensate_kernel_flag);
 
int find_minimum(double *input, double *input_minimum_value, int r, int c);
double abs_value(double input);

//int identity(double *input, int r, int c);
int create_exponential_kernel(double *output, double *C00, double *C01, double *C11, double *sigma_1, double *sigma_2, double *theta, double h, double *input_map, int rows, int cols, int i, int j, int half_ksize, double *A, double *phase_shift, double *dip_compensate_kernel_flag, int ext_length);

int copy_data_into_vec(double *z, double *input_data, double *input_map, int *number_of_elements, int rows, int cols, int i, int j, int half_ksize);

int initialize_to_zero(double *in, int rows, int cols);

double square_val(double input);

//int find_maximum_vector(double *minimum, double *input, int samples);
//int create_reweighted_kernel(double *W_reweighted, double *z, int number_of_samples, double *X, double *s);
//int copy_data(double *output, double *input, int r, int c);
//int absolute_value(int *input);
//int absolute_value(double *input);
//int polar_ribiere_optimization(double *output, double *Xw, double *X, double *z, double *s_init, int number_of_samples );


//Initialize mex functions
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 
    mxArray *output_image_mxArray, *output_image_label_mxArray; // *feature_location_x_mxArray, *feature_location_y_mxArray;

    const mwSize *dims;
    
    double *input_image, *input_label;
    double *output_image, *output_image_label, *target_template_label; 
    //covariance matrices
    
    double *C00, *C01, *C11;
    //kernel 
    double *h, *kernel_size_mat;
   
    double *sigma_1, *sigma_2, *theta;
    
    double *A, *phase_shift, *dip_compensate_kernel_flag;
    
    int rows, cols, ext_length;
   
    dims = mxGetDimensions(prhs[0]);
    rows = (int)dims[0];
    cols = (int)dims[1];
    
    //Initialize input into pointer
    input_image = mxGetPr(prhs[0]);
    input_label = mxGetPr(prhs[1]);
 
    target_template_label = mxGetPr(prhs[2]);
    
    //get h
    h = mxGetPr(prhs[3]);
    
    //kernel_size = (int)(*mxGetPr(prhs[9]));
    kernel_size_mat = mxGetPr(prhs[4]);
    
    C00 = mxGetPr(prhs[5]);
	C01 = mxGetPr(prhs[6]);
	C11 = mxGetPr(prhs[7]);
    
    sigma_1 = mxGetPr(prhs[8]);
    sigma_2 = mxGetPr(prhs[9]);
    theta = mxGetPr(prhs[10]);
    
	ext_length = (int)(*mxGetPr(prhs[11]));
    
    //dip parameters
    A = mxGetPr(prhs[12]);
    phase_shift = mxGetPr(prhs[13]);
    dip_compensate_kernel_flag = mxGetPr(prhs[14]);
    
    //Initialize output_mxArray and output;  
    output_image_mxArray = plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    output_image_label_mxArray = plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    
    output_image = mxGetPr(output_image_mxArray);
    output_image_label = mxGetPr(output_image_label_mxArray);
    
    
    create_image_allignment_template(output_image, output_image_label, input_image, input_label, target_template_label, h, kernel_size_mat, C00, C01, C11, sigma_1, sigma_2, theta, ext_length, rows, cols, A, phase_shift, dip_compensate_kernel_flag);
    
}


//int calculate_prediction_cost(double *output, double *input, double *h, double *kernel_size_mat, int rows, int cols, double *input_map, double *C00, double *C01, double *C11, int ext, double *output_label)
int create_image_allignment_template(double *output_image, double *output_image_label, double *input_image, double *input_label, double *target_template_label, double *h, double *kernel_size_mat, double *C00, double *C01, double *C11, double *sigma_1, double *sigma_2, double *theta, int ext_length, int rows, int cols, double *A, double *phase_shift, double *dip_compensate_kernel_flag){
    
    int max_number_of_samples = (int)square_val(61);

    double *W_kernel;
	double *z;
	double h_val;
	
    int kernel_size;
    int half_kernel_size;
    int y_ind, x_ind, l;
    
    int number_of_elements;
	double prediction, sum_weight;
    
    int stop_here;
    
    //int pixel_ind;
	//int stop_here;
	//define loop variable
	//int l;
	//double sum_weight;
	//double prediction;
	//double *error_sig;
	
	W_kernel = (double*)malloc( max_number_of_samples*sizeof( double) );
    z = (double*)malloc(max_number_of_samples*sizeof( double) );
	//error_sig =  (double*)malloc( max_number_of_samples*sizeof( double) );
    
    //pixel_ind = 0;
    initialize_to_zero(output_image, rows, cols);
    initialize_to_zero(output_image_label, rows, cols);
    
    for(y_ind=ext_length; y_ind<(rows-ext_length); y_ind++){
        for(x_ind=ext_length; x_ind<(cols-ext_length); x_ind++){
            
            
            if (y_ind == 814 && x_ind == 35)
                stop_here = 1;
            
            //check whether the current pixel needs to be predicted
            if (target_template_label[y_ind + x_ind*rows] == 1){
                
                h_val = h[y_ind + x_ind*rows];
                kernel_size = (int)kernel_size_mat[y_ind + x_ind*rows];
                half_kernel_size = (kernel_size - 1)/2;
                
                create_exponential_kernel(W_kernel, C00, C01, C11, sigma_1, sigma_2, theta, h_val, input_label, rows, cols, y_ind, x_ind, half_kernel_size, A, phase_shift, dip_compensate_kernel_flag, ext_length);
                copy_data_into_vec(z, input_image, input_label, &number_of_elements, rows, cols, y_ind, x_ind, half_kernel_size);
                
                //loop through the two dimensions
                prediction = 0;
                sum_weight = 0;
                      
                //set the minimum number of elements
                if (number_of_elements > 5){
                   
                    for(l=0;l<number_of_elements;l++){
                        prediction = prediction + W_kernel[l]*z[l];
                        sum_weight = sum_weight + W_kernel[l];
                    }
                    prediction = prediction/sum_weight;
                    
                    //pixel_ind = pixel_ind + 1;
                    output_image_label[y_ind + x_ind*rows] = 1;
                    output_image[y_ind + x_ind*rows] = prediction;                    
                }            
            }                                              
        }
    }
    
	free(W_kernel);
	free(z);

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


int create_exponential_kernel(double *output, double *C00, double *C01, double *C11, double *sigma_1, double *sigma_2, double *theta, double h, double *input_map, int rows, int cols, int i, int j, int half_ksize, double *A, double *phase_shift, double *dip_compensate_kernel_flag, int ext_length)
{
	int ii,jj,k = 0;
    
    double cov_mat[4];
	double cov_det;
    
    double pi = 3.142;
    
    double shift;
    double theta_val;
    
    double tangent_gradient, f_x_x_i, f_x_i;
    double image_length;
    
    double R_theta[4];
    double R_theta_transpose[4];
    double buffer[4];
    double buffer_2[4];

    double lambda_mat[4];
    
    image_length = (double)(cols -2*ext_length); 
    
	for(ii= -half_ksize; ii <= half_ksize; ii++){          
		for(jj = -half_ksize; jj <= half_ksize; jj++){
			//check the label to see, whether this pixel is located in the label
			if (((i + ii) >= 0) && ((i + ii) < rows) && ((j+jj) >= 0) && ((j+jj) < cols)){
                
				//check wether the pixel belongs to the central pixel
				if (!((ii == 0) && (jj == 0))){

					cov_mat[0] = C00[(i+ii) +(j+jj)*rows];               
					cov_mat[1] = C01[(i+ii) +(j+jj)*rows];                
					cov_mat[2] = C01[(i+ii) +(j+jj)*rows];                
					cov_mat[3] = C11[(i+ii) +(j+jj)*rows];
                
					
                
					if (input_map[ (i+ii) +(j+jj)*rows] == 1){
                        
                        if (*dip_compensate_kernel_flag == 0){
                            shift = 0;
                            cov_mat[0] = C00[(i+ii) +(j+jj)*rows];
                            cov_mat[1] = C01[(i+ii) +(j+jj)*rows];
                            cov_mat[2] = C01[(i+ii) +(j+jj)*rows];
                            cov_mat[3] = C11[(i+ii) +(j+jj)*rows];
                            cov_det = cov_mat[0]*cov_mat[3] - cov_mat[1]*cov_mat[2];
                        }
                        else{
                            
                            
                            tangent_gradient = (*A)*(2*pi/image_length)*cos((2*pi/image_length)*(j+1+jj-(*phase_shift)-(double)ext_length));
                            
                            theta_val = atan(tangent_gradient);
                            //theta_val = theta[(i+ii) +(j+jj)*rows];         
                     
                            lambda_mat[0] = sigma_1[(i+ii) +(j+jj)*rows];
                            lambda_mat[1] = 0;
                            lambda_mat[2] = 0;
                            lambda_mat[3] = sigma_2[(i+ii) +(j+jj)*rows];
                            
                            
                         
                            R_theta[0] = cos(theta_val);
                            R_theta[1] = sin(theta_val);
                            R_theta[2] = -sin(theta_val);
                            R_theta[3] = cos(theta_val);
                            
                            
                      
                            matrix_multiplication(buffer, R_theta, lambda_mat, 2, 2, 2, 2);
                       
                            matrix_multiplication(buffer_2, buffer, lambda_mat, 2, 2, 2, 2);
                    
                            matrix_transpose(R_theta_transpose, R_theta, 2, 2);
                  
                            matrix_multiplication(cov_mat, buffer_2, R_theta_transpose, 2, 2, 2, 2);
                            
                            
                            
                            //tangent_gradient = A*(2*pi/image_length)*cos((2*pi/image_length)*(x_i-phase_shift-ext_length));
        
                            
                            
                            cov_det = cov_mat[0]*cov_mat[3] - cov_mat[1]*cov_mat[2];
                            
                            f_x_i = (*A)*sin((2*pi/image_length)*(j+1-(*phase_shift)-(double)ext_length));
                            f_x_x_i = (*A)*sin((2*pi/image_length)*((j+1+jj-(*phase_shift)-(double)ext_length)));
                            shift = tangent_gradient*(-jj) - ( f_x_i - f_x_x_i );
                            shift = -shift;
                            
                            if (((ii+i == 89) && (jj+j == 235))){
                                    
                                if (((i == 108) && (j == 211))){
                                    
                                    mexPrintf("%d", tangent_gradient);
                                }
                            }
                            
                            
                        }
                            
                        
						output[k] = 100*(1/(2*pi*square_val(h)))* pow(cov_det, -0.5) * exp(-(1/cov_det)*(square_val(jj)*cov_mat[3] - 2*(ii+shift)*(jj)*cov_mat[1] + square_val(ii+shift)*cov_mat[0])/(2*square_val(h)));
						
						//output[k] = 100*(1/(2*3.142*square_val(h)))* pow(cov_det, -0.5) * exp(-(1/cov_det)*(square_val(ii)*cov_mat[3] - 2*(ii)*(jj)*cov_mat[1] + square_val(jj)*cov_mat[0])/(2*square_val(h)));
						
						//output[k] = exp(-(1/cov_det)*(square_val(jj)*cov_mat[3] - 2*(ii)*(-jj)*cov_mat[1] + square_val(ii)*cov_mat[0])/(2*square_val(h)));
						k++;
					}
				
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


int copy_data_into_vec(double *z, double *input_data, double *input_map, int *number_of_elements, int rows, int cols, int i, int j, int half_ksize)
{
	int ii,jj,k = 0;
    int stop_here;

	*number_of_elements = 0;

	for(ii= -half_ksize; ii <= half_ksize; ii++){          
		for(jj = -half_ksize; jj <= half_ksize; jj++){

			
			//check the label to see, whether this pixel is located in the label
			if (((i + ii) >= 0) && ((i + ii) < rows) && ((j+jj) >= 0) && ((j+jj) < cols)){

				if (((ii == 0) && (jj == 0))){
					stop_here = 1;
				}

				if (!((ii == 0) && (jj == 0))){
					if (input_map[ (i+ii) +(j+jj)*rows] == 1){
						z[k] = input_data[(i+ii) +(j+jj)*rows];
						//output[k] = exp(-(1/cov_det)*(square_val(jj)*cov_mat[3] - 2*(ii)*(-jj)*cov_mat[1] + square_val(ii)*cov_mat[0])/(2*square_val(h)));
						k++;
					}
				}
			}
		
		}        
	}

	*number_of_elements = k;
}


int initialize_to_zero(double *in, int rows, int cols){
    
    int x_ind, y_ind;
    for(x_ind=0; x_ind<cols; x_ind++){
        for(y_ind=0; y_ind<rows; y_ind++){
            in[y_ind + x_ind*rows] = 0;
        }
    }

}


