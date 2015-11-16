
//Schlumberger 2012 - written by Andriy Gelman (AGelman@slb.com)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mex.h>
#include "matrix_toolbox_ag.h"

FILE *f1;

int template_matching(double *cost, double *number_of_predictions_per_feature_point, double *unreliable_template, double *input, double *input_x, double *input_y, double *input_label, double *template, double *template_x, double *template_y, double *template_label, int ext_length, double *feature_location_x, double *feature_location_y, int number_of_feature_points_per_row, int shift_y, int window_size_x, int window_size_y, int rows, int cols);

int find_minimum(double *input, double *input_minimum_value, int r, int c);
double abs_value(double input);

//int identity(double *input, int r, int c);
int create_exponential_kernel(double *output, double *C00, double *C01, double *C11, double h, double *input_map, int rows, int cols, int i, int j, int half_ksize);

int copy_data_into_vec(double *z, double *input_data, double *input_map, int *number_of_elements, int rows, int cols, int i, int j, int half_ksize);

double square_val(double input);

int find_maximum_vector(double *minimum, double *input, int samples);

int create_reweighted_kernel(double *W_reweighted, double *z, int number_of_samples, double *X, double *s);

int copy_data(double *output, double *input, int r, int c);

//int absolute_value(int *input);
int absolute_value(double *input);
//int polar_ribiere_optimization(double *output, double *Xw, double *X, double *z, double *s_init, int number_of_samples );


//Initialize mex functions
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 
    mxArray *cost_mxArray, *number_of_predictions_per_feature_point_mxArray, *unreliable_template_mxArray; // *feature_location_x_mxArray, *feature_location_y_mxArray;

    const mwSize *dims;
    
    double *input, *input_x, *input_y, *input_label;
    double *template, *template_x, *template_y, *template_label;
    
    double *alpha;
    
	//double *output, *output_label;
    int shift_y;
	int number_of_feature_points_per_row;
    double *feature_location_x, *feature_location_y;
    
	double *cost, *number_of_predictions_per_feature_point, *unreliable_template;

	//calculate_cost(input, input_label_lower, input_label_upper, h, kernel_size_mat, C00, C01, C11, ext_length, feature_location_x, feature_location_y, shift_y);

    //optimization parameter defines whether L1 or L2 denoising is performed. To select L1 input 1, L2 = 2. L1 minimization requires an initialization value. 
    //int kernel_size;
  
    int rows, cols, ext_length;
	int rows_cost, cols_cost;
    int window_size_x, window_size_y, i;

    dims = mxGetDimensions(prhs[0]);
    rows = (int)dims[0];
    cols = (int)dims[1];
    
    //Initialize input into pointer
    input = mxGetPr(prhs[0]); 
    input_x = mxGetPr(prhs[1]);
    input_y = mxGetPr(prhs[2]);
	input_label = mxGetPr(prhs[3]);
    
    template = mxGetPr(prhs[4]);
    template_x = mxGetPr(prhs[5]);
    template_y = mxGetPr(prhs[6]);
    
    template_label = mxGetPr(prhs[7]);
    
	ext_length = (int)(*mxGetPr(prhs[8]));
    
	feature_location_x = mxGetPr(prhs[9]);
	feature_location_y = mxGetPr(prhs[10]);
    
    dims = mxGetDimensions(prhs[9]);
    rows_cost = (int)dims[0];
    cols_cost = (int)dims[1];
    number_of_feature_points_per_row = rows_cost*cols_cost;
    
	shift_y = (int)(*mxGetPr(prhs[11]));

	window_size_x = (int)(*mxGetPr(prhs[12]));
    window_size_y = (int)(*mxGetPr(prhs[13]));
    
    alpha = mxGetPr(prhs[14]);
    
    //Initialize output_mxArray and output;  
    cost_mxArray = plhs[0] = mxCreateDoubleMatrix(1, number_of_feature_points_per_row, mxREAL);
    number_of_predictions_per_feature_point_mxArray = plhs[1] = mxCreateDoubleMatrix(1, number_of_feature_points_per_row, mxREAL);
    unreliable_template_mxArray = plhs[2] = mxCreateDoubleMatrix(1, number_of_feature_points_per_row, mxREAL);
 
    cost = mxGetPr(cost_mxArray);
	number_of_predictions_per_feature_point = mxGetPr(number_of_predictions_per_feature_point_mxArray);
    unreliable_template = mxGetPr(unreliable_template_mxArray);
    
    //for(i=0;i<number_of_feature_points_per_row; i++){
    //    feature_location_x[i] = feature_location_x[i] - 1;
    //    feature_location_y[i] = feature_location_y[i] - 1;
    //}
    
    template_matching(cost, number_of_predictions_per_feature_point, unreliable_template, input, input_x, input_y, input_label, template, template_x, template_y, template_label, ext_length, feature_location_x, feature_location_y, number_of_feature_points_per_row, shift_y, window_size_x, window_size_y, rows, cols);
}



 
int template_matching(double *cost, double *number_of_predictions_per_feature_point, double *unreliable_template, double *input, double *input_x, double *input_y, double *input_label, double *template, double *template_x, double *template_y, double *template_label, int ext_length, double *feature_location_x, double *feature_location_y, int number_of_feature_points_per_row, int shift_y, int window_size_x, int window_size_y, int rows, int cols)
{
    
 	int x_ind, y_ind;
    int iter;
	int number_of_elements;
	int half_kernel_size;
	int pixel_ind;
	int stop_here;
	//define loop variable
	int l;
	double sum_weight;
	double prediction;
	double *error_sig;
    double *error_sig_x;
    double *error_sig_y;
    double *error_sig_gradient; 
    double cost_intensity; 
    double cost_gradient;
    double tau = 30;
	double beta = 1;
    
    double buffer_x;
    double buffer_y;
    
    int max_number_of_samples = (int)square_val(61);

    //this specifies that 20% of the pixels in the original input have to be non-zero
    double min_number_of_non_zero_pixels = (window_size_x*2) * (window_size_y*2) * 0.4;
    int non_zero_pixels;
    
    //allocate memory to the error signal
	error_sig =  (double*)malloc( max_number_of_samples*sizeof( double) );
    error_sig_gradient =  (double*)malloc( max_number_of_samples*sizeof( double) );
    error_sig_x =  (double*)malloc( max_number_of_samples*sizeof( double) );
    error_sig_y =  (double*)malloc( max_number_of_samples*sizeof( double) );

    
    for(iter = 0; iter<number_of_feature_points_per_row; iter++){
        
        /*
        pixel_ind = 0;     
        for(x_ind = feature_location_x[iter] - window_size_x; x_ind <= feature_location_x[iter] + window_size_x; x_ind++){
            for(y_ind = feature_location_y[iter] - window_size_y - shift_y; y_ind <= feature_location_y[iter] + window_size_y - shift_y; y_ind++){
                
                //y_ind = y_ind - shift_y;
                
                if (input_label[y_ind + shift_y + x_ind*rows] == 1){                    
                    if (template_label[y_ind + x_ind*rows] == 1){                        
                        error_sig[pixel_ind] = input[y_ind + shift_y + x_ind*rows] - template[y_ind + x_ind*rows];
                        pixel_ind = pixel_ind + 1;                        
                        if (pixel_ind == 42){
                            stop_here = 1;
                        }
                    }
                    
                }
            }
        }
        */
        
        
        //calculate number of non-zero pixels
        non_zero_pixels = 0;
        for(x_ind = feature_location_x[iter] - window_size_x; x_ind <= feature_location_x[iter] + window_size_x; x_ind++){
            for(y_ind = feature_location_y[iter] - window_size_y; y_ind <= feature_location_y[iter] + window_size_y; y_ind++){
                
                //y_ind = y_ind - shift_y;              
                if (input_label[y_ind + x_ind*rows] == 1){                    
                    if (input[y_ind + x_ind*rows] > 0){
                        non_zero_pixels++;
                    }             
                }
            }
        }
        //mexPrintf("%d \n", non_zero_pixels);
        
        unreliable_template[iter] = 1;
        if (non_zero_pixels > (int)min_number_of_non_zero_pixels){
            unreliable_template[iter] = 0;
        }
        
        
        pixel_ind = 0;
        for(x_ind = feature_location_x[iter] - window_size_x; x_ind <= feature_location_x[iter] + window_size_x; x_ind++){
            for(y_ind = feature_location_y[iter] - window_size_y; y_ind <= feature_location_y[iter] + window_size_y; y_ind++){
                
                //y_ind = y_ind - shift_y;
                //mexPrintf("%d \n", y_ind + shift_y + x_ind*rows);
                
                if (input_label[y_ind + x_ind*rows] == 1){
                    
                    if ((y_ind + shift_y >= 0) && (y_ind + shift_y <= rows)){
                        
                        if (template_label[y_ind + shift_y + x_ind*rows] == 1){
                            error_sig[pixel_ind] = input[y_ind + x_ind*rows] - template[y_ind + shift_y + x_ind*rows];
                            pixel_ind = pixel_ind + 1;
                            if (pixel_ind == 42){
                                stop_here = 1;
                            }
                        }
                    }
                }
            }        
        }
        //unreliable_template[iter] = 1;
        if (pixel_ind < (int)min_number_of_non_zero_pixels){
            unreliable_template[iter] = 1;
        }

        
        pixel_ind = 0;     
        for(x_ind = feature_location_x[iter] - window_size_x; x_ind <= feature_location_x[iter] + window_size_x; x_ind++){
            for(y_ind = feature_location_y[iter] - window_size_y; y_ind <= feature_location_y[iter] + window_size_y; y_ind++){
                
                //y_ind = y_ind - shift_y;
                
                if (input_label[y_ind + x_ind*rows] == 1){
                    if ((y_ind + shift_y >= 0) && (y_ind + shift_y <= rows)){
                        if (template_label[y_ind + shift_y + x_ind*rows] == 1){
                            error_sig_x[pixel_ind] = input_x[y_ind + x_ind*rows] - template_x[y_ind + shift_y + x_ind*rows];
                            pixel_ind = pixel_ind + 1;
                            if (pixel_ind == 42){
                                stop_here = 1;
                            }
                        }
                    }
                    
                }
            }
        }
        if (pixel_ind < (int)min_number_of_non_zero_pixels){
            unreliable_template[iter] = 1;
        }
        
        pixel_ind = 0;     
        for(x_ind = feature_location_x[iter] - window_size_x; x_ind <= feature_location_x[iter] + window_size_x; x_ind++){
            for(y_ind = feature_location_y[iter] - window_size_y; y_ind <= feature_location_y[iter] + window_size_y; y_ind++){
                
                //y_ind = y_ind - shift_y;
                
                if (input_label[y_ind + x_ind*rows] == 1){
                    if ((y_ind + shift_y >= 0) && (y_ind + shift_y <= rows)){
                        if (template_label[y_ind + shift_y + x_ind*rows] == 1){
                            error_sig_y[pixel_ind] = input_y[y_ind + x_ind*rows] - template_y[y_ind + shift_y + x_ind*rows];
                            pixel_ind = pixel_ind + 1;
                            if (pixel_ind == 42){
                                stop_here = 1;
                            }
                        }
                    }

                    
                }
            }
        }
        if (pixel_ind < (int)min_number_of_non_zero_pixels){
            unreliable_template[iter] = 1;
        }
        
        //cost[iter] = 0;
        cost_intensity = 0;
        //calculate total cost
        for(l=0; l <pixel_ind; l++){
            
            if (error_sig[l] > 0){
                cost_intensity = cost_intensity + error_sig[l];
            }
            else{
                cost_intensity = cost_intensity - error_sig[l];
            }
            
            //cost_intensity = cost_intensity + error_sig[l]*error_sig[l];
            if (l==35){
                stop_here = 1;
            }
            //mexPrintf("%d \n", l);
            //mexPrintf("%d \n", cost[0]);
        }
        
        
        //cost[iter] = 0;
        cost_gradient = 0;
        //calculate total cost
        for(l=0; l <pixel_ind; l++){
            cost_gradient = cost_gradient + sqrt(error_sig_x[l]*error_sig_x[l] + error_sig_y[l]*error_sig_y[l]);
            if (l==35){
                stop_here = 1;
            }
            //mexPrintf("%d \n", l);
            //mexPrintf("%d \n", cost[0]);
        }
        
        
        //cost[iter] = 0;
        //cost_gradient = 0;
        //calculate total cost
        error_sig_gradient[l] = 0;
        for(l=0; l <pixel_ind; l++){
            
            error_sig_gradient[l] = exp(-beta*tau*sqrt(error_sig_x[l]*error_sig_x[l] + error_sig_y[l]*error_sig_y[l]));
            
            //mexPrintf("%d \n", l);
            //mexPrintf("%d \n", cost[0]);
        }
        
        for (l=0; l <pixel_ind; l++){
            if (error_sig[l] < 0){
                error_sig[l] = -error_sig[l];
            }
            error_sig[l] = exp(-beta*error_sig[l]);
        }
        
        cost[iter] = 0;
        for (l=0; l <pixel_ind; l++){
          
              if (error_sig[l] == 0 && error_sig_gradient[l] == 0){
                mexPrintf("Zero \n");
            }
            cost[iter] = cost[iter] - log(error_sig[l] + error_sig_gradient[l])/beta;
        }
        
        
        //cost[iter] = 0;
        //tau = 20;
        //cost[iter] = cost_intensity + tau*cost_gradient;
        cost[iter] = cost_intensity;
        
        //normalize by the number of pixels
        //cost[iter] = cost[iter]/(pixel_ind-1);
		number_of_predictions_per_feature_point[iter] = pixel_ind-1;
		    
    }
    
	free(error_sig);
    free(error_sig_gradient);
    free(error_sig_x);
    free(error_sig_y);
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
                
				//check wether the pixel belongs to the central pixel
				if (!((ii == 0) && (jj == 0))){

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
