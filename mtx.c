#include <stdio.h>
#include <math.h>
#include "mtx.h"

/*
 * Function: mtx_create
 * ------------------------
 *  Construct mtx_matrix from data input
 *
 *  Result: mtx_matrix struct
 */

void mtx_create(int rows, int cols, float *data, struct mtx_matrix* mtx_out){
    int row;
    int col;
   
    if (rows < MAX_ROWS)  mtx_out->rows = rows;
    else mtx_out->rows = MAX_ROWS;

    if (cols < MAX_COLS)  mtx_out->cols = cols; 
    else mtx_out->cols = MAX_COLS;
    
    for (row=1;row <= rows;row++){
        for (col=1;col <= cols;col++){
            mtx_out->data[(row-1)*mtx_out->cols + (col-1)] = data[(row-1)*mtx_out->cols + (col-1)];
        }
    }

    return;
}

/*
 * Function: mtx_create_ones
 * ------------------------
 *  Construct mtx_matrix of ones
 *
 *  Result: mtx_matrix struct
 */

void mtx_create_ones(int rows, int cols, struct mtx_matrix* mtx_out){
    int row;
    int col;
   
    if (rows < MAX_ROWS)  mtx_out->rows = rows;
    else mtx_out->rows = MAX_ROWS;

    if (cols < MAX_COLS)  mtx_out->cols = cols; 
    else mtx_out->cols = MAX_COLS;
    
    for (row=1;row <= rows;row++){
        for (col=1;col <= cols;col++){
            mtx_out->data[(row-1)*mtx_out->cols + (col-1)] = 1.0;
        }
    }

    return;
}

/*
 * Function: mtx_get
 * ------------------------
 * Get value from mtx_matrix
 *
 * Result: float from i,j
 *
 */

float mtx_get(int row, int col, struct mtx_matrix* mtx_a){
    if (row<=mtx_a->rows && col<=mtx_a->cols){
        return mtx_a->data[(row-1)*mtx_a->cols + (col-1)];    
    } else {
        return -1.0;
    }
}

/*
 * Function: mtx_set
 * ------------------------
 * Set value in mtx_matrix
 *
 */

void mtx_set(int row, int col, struct mtx_matrix* mtx_a, float data){
    if (row<=mtx_a->rows && col<=mtx_a->cols){
        mtx_a->data[(row-1)*mtx_a->cols + (col-1)] = data;
    }
}

/*
 * Function: mtx_print
 * ------------------------
 *  Print mtx_matrix to stdout
 *
 */

void mtx_print(struct mtx_matrix* mtx_a){
    int row;
    int col;

    if (mtx_a->rows < MAX_ROWS && mtx_a->rows > 0 &&
        mtx_a->cols < MAX_COLS && mtx_a->cols >0) {
        printf("\n\r");
        for (row=1;row <= mtx_a->rows;row++){
            printf("|  ");
            for (col=1;col <= mtx_a->cols;col++){
                printf("%f  ",mtx_get(row,col,mtx_a));    
            }
            printf("|\n\r");
        }
        printf("\n\r");
   } else {
       printf("EMPTY MATRIX\n\r");
   }
   return;
}

/*
 * Function: mtx_sum
 * ------------------------
 *  Add two matrices
 *
 *  Result: Matrix of sum of inputs
 */

int mtx_sum(struct mtx_matrix* mtx_a, struct mtx_matrix* mtx_b, 
                struct mtx_matrix* mtx_out){
    int row;
    int col;
    float data;

    if ((mtx_a->rows != mtx_b->rows) || (mtx_a->cols != mtx_b->cols) || 
        (mtx_out->rows != mtx_b->rows) || (mtx_out->cols != mtx_b->cols)) {
        return -1;
    }

    for (row=1;row <= mtx_a->rows;row++){
        for (col=1;col <= mtx_a->cols;col++){
            data = mtx_get(row,col,mtx_a) + mtx_get(row,col,mtx_b);
            mtx_set(row,col,mtx_out,data);    
        }
    }

    return 0;
}

/*
 * Function: mtx_scaler
 * ------------------------
 *  Multiply matrix elements by scaler
 *
 *  Result: Matrix of scaled values
 */

int mtx_scale(struct mtx_matrix* mtx_a, float scaler, 
                struct mtx_matrix* mtx_out){
    int row;
    int col;
    float data;

    if ((mtx_out->rows != mtx_a->rows) || (mtx_out->cols != mtx_a->cols)) {
        return -1;
    }

    for (row=1;row <= mtx_a->rows;row++){
        for (col=1;col <= mtx_a->cols;col++){
            data = scaler*mtx_get(row,col,mtx_a);
            mtx_set(row,col,mtx_out,data);    
        }
    }

    return 0;
}
