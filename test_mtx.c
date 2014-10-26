#include <stdio.h>
#include "mtx.h"

int main(void) {
    struct mtx_matrix mtx_1_r; 
    struct mtx_matrix mtx_2_r; 
    struct mtx_matrix mtx_3_r;
    struct mtx_matrix mtx_4_r;
    struct mtx_matrix mtx_5_r;
    struct mtx_matrix mtx_6_r;
    struct mtx_matrix mtx_7_r;

    float data_1[] = {1, 2, 3,
                      -4, 5, 6,
                      7, 8, -9};

    float data_2[] = {10, 20, 30};

    mtx_create(3,3,data_1,&mtx_1_r);
    mtx_create(3,1,data_2,&mtx_2_r);
    mtx_create_ones(3,3,&mtx_3_r);
    mtx_create_ones(3,3,&mtx_4_r);
    mtx_create_ones(3,1,&mtx_5_r);
    mtx_create_ones(3,3,&mtx_6_r);
    mtx_create_ones(3,3,&mtx_7_r);

    printf("Matrix 1\n\r");
    mtx_print(&mtx_1_r);
    printf("Matrix 2\n\r");
    mtx_print(&mtx_2_r);

    if (mtx_sum(&mtx_1_r,&mtx_2_r,&mtx_3_r)==0){
        printf("Matrix 1 + Matrix 2\n\r\n\r");
        mtx_print(&mtx_3_r);
    } else {
        printf("Sum Failed\n\r\n\r"); 
    }

    if (mtx_scale(&mtx_1_r,5.0,&mtx_4_r)==0){
        printf("Matrix 1 * 5.0\n\r");
        mtx_print(&mtx_4_r);
    } else {
        printf("Scale Failed\n\r\n\r");
    }

    if (mtx_mult(&mtx_1_r,&mtx_2_r,&mtx_5_r)==0){
        printf("Matrix 1 * Matrix 2\n\r");
        mtx_print(&mtx_5_r);
    } else {
        printf("Multiply Failed\n\r\n\r");
    }

    if (mtx_trans(&mtx_1_r,&mtx_6_r)==0){
        printf("Transpose Matrix 1\n\r");
        mtx_print(&mtx_6_r);
    } else {
        printf("Transpose Failed\n\r\n\r");
    }
    
    if (mtx_inv(&mtx_1_r,&mtx_7_r)==0){
        printf("Inverse Matrix 1\n\r");
        mtx_print(&mtx_7_r);
    } else {
        printf("Inverse Failed\n\r\n\r");
    }

    float det;
    det = mtx_det(&mtx_1_r);
    printf("Matrix 1 Determinant: %f\n\r",det);

    float norm;
    norm = mtx_norm(&mtx_2_r);
    printf("Matrix 2 Norm: %f\n\r",norm);

    float trace;
    trace = mtx_trace(&mtx_1_r);
    printf("Matrix 1 Trace: %f\n\r",trace);

    float max;
    int lr;
    int lc;
    max = mtx_max(&mtx_1_r,&lr,&lc);
    printf("Matrix 1 Max: %f Row: %d Col: %d\n\r",max,lr,lc);
}
