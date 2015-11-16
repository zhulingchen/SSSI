

//Schlumberger 2012 - written by Andriy Gelman (AGelman@slb.com)

int matrix_multiplication(double *output, double *A, double *B, int row_a, int col_a, int row_b, int col_b);
int matrix_copy_matrix(double *output, double *input, int rows, int cols);
int matrix_multiplication(double *output, double *A, double *B, int row_a, int col_a, int row_b, int col_b);
int matrix_add_two_rows(double *input, int rows, int cols, int row_location_1, int row_location_2, double scalar);
int matrix_row_permutation(double *input, int rows, int cols, int row_1, int row_2, double *row_buffer);
int matrix_scale_row(double *input, int rows, int cols, int i, double scalar);
int matrix_inversion(double *output, double *input, int rows, int cols, double *buffer_input, double *row_buffer);
int matrix_init_ones(double *output, int rows, int cols);
int matrix_scale_cols_by_vector_component(double *output, double *input, int rows_mat, int cols_mat, double *vector, int cols_vec);
int matrix_initialize_vector(double *output, double min_val, double max_val);

int matrix_add_matrix(double *output, double *in_1, double *in_2, int rows, int cols);
int matrix_matrix_multiply_by_scalar(double *output, double *input, int rows, int cols, double scalar);
int matrix_matrix_multiply_element_by_element(double *output, double *in_1, double *in_2, int rows, int cols);

int identity(double *input, int r, int c);
FILE *f1;

int identity(double *input, int r, int c)
{
    int i ,j;
    
    for (i=0;i<r;i++){
        for (j=0;j<r;j++){
            if (i != j)
                input[i+j*r] = 0;
            else
                input[i+j*r] = 1;       
    }
    }
    
    return 0;
    
}


int matrix_multiplication(double *output, double *A, double *B, int row_a, int col_a, int row_b, int col_b)
{
    int i,j,k;
    
    for(i=0;i<row_a;i++)
			{
			for(j=0;j<col_b;j++)
				{
				output[i + j*row_a]=0;
				for(k=0;k<col_a;k++)
					{
                        output[i + j*row_a]=output[i + j*row_a]+ (A[i + k*row_a]*B[k + j*row_b]);
					}
				}
			}
    return 0;
}


int matrix_inversion(double *output, double *input, int rows, int cols, double *buffer_input, double *row_buffer)
{
    int i, j;
   
    double pivot_value, scale;
    double lambda_inv = 1.0e-8;
    
	pivot_value = 0;

    //initialize the output to be the identity matrix
    matrix_copy_matrix(buffer_input, input, rows, cols);
    
    //add a scalar to each row to ensure that the matrix is invertible
    for (i = 0; i < rows; i++){
        buffer_input[i + i*rows] = buffer_input[i + i*rows] + lambda_inv; 
        input[i + i*rows] = input[i + i*rows] + lambda_inv; 
               
    }

    //init output as identity
    identity(output, rows, cols);

    //work with the diagonal elements as pivot points
    for(j = 0; j<cols; j++){
        
        //find the first non-zero row in the pivot 
        for (i = j; i < rows; i++){
            if (buffer_input[i + j*rows] > 1e-6){
                
                //permute the rows and then break

                matrix_row_permutation(buffer_input, rows, cols, j, i, row_buffer);
                matrix_row_permutation(output, rows, cols, j, i, row_buffer);
                break;
            }         
        }
        
        pivot_value = buffer_input[j + j*rows];
        //pivot_value = pivot_value + 1.0e-5;

        matrix_scale_row(buffer_input, rows, cols, j, 1/pivot_value);
        matrix_scale_row(output, rows, cols, j, 1/pivot_value);
       
        pivot_value = 1;
        
        for (i = 0; i < rows; i++){
            if (i != j){
                
                //i is not equal to j
                
                //apply a row operation by using the pivot value to set the value in 
                //part of the matrix to zero#
                //j is the matrix value we are trying to remove
                //i is the pivot location
				scale =-buffer_input[i+j*cols];
                matrix_add_two_rows(buffer_input, rows, cols, i, j, scale);
                matrix_add_two_rows(output, rows, cols, i, j, scale);
         
            }    
        }
    }


}


int matrix_scale_row(double *input, int rows, int cols, int i, double scalar){

    //function scales the i_th row in the matrix by the scalar 
    
    int j;
    
    for (j = 0; j < cols; j++){
        input[i + j*cols] = input[i + j*cols]*scalar; 
    }
}


int matrix_add_two_rows(double *input, int rows, int cols, int row_location_1, int row_location_2, double scalar){

    //function adds row_location_2 to row_location_1 in the matrix
    //row_location_1 = row_location_1 + scalar*row_location_2
    int j;
    
    for (j = 0; j < cols; j++){
        input[row_location_1 + j*cols] = input[row_location_1 + j*cols] + scalar*input[row_location_2 + j*cols]; 
    }
    return 0;
}

int matrix_row_permutation(double *input, int rows, int cols, int row_1, int row_2, double *row_buffer){

    int j;
    //function switches two rows in the matrix
    //create a buffer matrix
    
	for(j=0; j<cols; j++){
         row_buffer[j] = input[row_1 + j*cols];
    }

    for(j=0; j<cols; j++){
        input[row_1 + j*cols] = input[row_2 + j*cols];
    }
    
    for(j=0; j<cols; j++){
        input[row_2 + j*cols] = row_buffer[j];
    }

}

int matrix_copy_matrix(double *output, double *input, int rows, int cols){

    int i,j;
    
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            output[i+j*rows] = input[i+j*rows];
        }
    }
}

int matrix_transpose(double *output, double *input, int rows_in, int cols_in){

    int i,j;

    for(i=0;i<rows_in;i++){
        for(j=0;j<cols_in;j++){
            output[j+i*cols_in] = input[i+j*rows_in];
        }
    }
}

int matrix_init_ones(double *output, int rows, int cols)
{
     
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            output[i+j*rows] = 1;
        }
    }
}

int matrix_scale_cols_by_vector_component(double *output, double *input, int rows_mat, int cols_mat, double *vector, int cols_vec)
{

    int i,j;
     
    for(j=0;j<cols_mat;j++){
        for(i=0;i<rows_mat;i++){
            output[i+j*rows_mat] = input[i+j*rows_mat]*vector[j];
        }
    }

}

int matrix_initialize_vector(double *output, double min_val, double max_val)
{
	//note the memory must already be allocated to output, which must have the same dimensions as (max_val - min_val + 1)
    double i;
    int k = 0; 
    for(i=min_val;i<=max_val;i++){
            output[k] = (double)i;
			k++;

    }
}

int matrix_add_matrix(double *output, double *in_1, double *in_2, int rows, int cols){

    int i,j;
    
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            output[i+j*rows] = in_1[i+j*rows] + in_2[i+j*rows];
        }
    }
}


int matrix_matrix_multiply_by_scalar(double *output, double *input, int rows, int cols, double scalar){

    int i,j;
    
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            output[i+j*rows] = input[i+j*rows]*scalar;
        }
    }
}

int matrix_matrix_sign_output(double *output, double *input, int rows, int cols)
{

	int i,j;

    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
			if (input[i+j*rows] > 0){
				output[i+j*rows] = 1;
			}
			if (input[i+j*rows] < 0)
			{
				output[i+j*rows] = -1;
			}
			
			if (input[i+j*rows] == 0)
			{
				output[i+j*rows] = 0;
			}
        }
    }
}

int matrix_matrix_abs_output(double *output, double *input, int rows, int cols)
{

	int i,j;

    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
			if (input[i+j*rows] > 0){
				output[i+j*rows] = input[i+j*rows];

			}
			if (input[i+j*rows] < 0)
			{
				output[i+j*rows] = - input[i+j*rows];
			}
			
			if (input[i+j*rows] == 0)
			{
				output[i+j*rows] = 0;
			}
        }
    }
}

int matrix_matrix_multiply_element_by_element(double *output, double *in_1, double *in_2, int rows, int cols){

    int i,j;
    
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            output[i+j*rows] = in_1[i+j*rows]*in_2[i+j*rows];
        }
    }
}