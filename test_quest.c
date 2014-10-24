#include <stdio.h>
#include "mtx.h"

int main(void) {
   struct mtx_matrix mtx_1_r; 
   struct mtx_matrix mtx_2_r; 
   struct mtx_matrix mtx_3_r;
   struct mtx_matrix mtx_4_r;

   float data_1[] = {1, 2, 3,
                     4, 5, 6,
                     7, 8, 9}; 

   float data_2[] = {10, 20, 30,
                     40, 50, 60,
                     70, 80, 90}; 

   mtx_create(3,3,data_1,&mtx_1_r);
   mtx_create(3,3,data_2,&mtx_2_r);
   mtx_create_ones(3,3,&mtx_3_r);
   mtx_create_ones(3,3,&mtx_4_r);
   mtx_scale(&mtx_1_r,5.0,&mtx_4_r);
   mtx_sum(&mtx_1_r,&mtx_2_r,&mtx_3_r);
   printf("Matrix 1\n\r");
   mtx_print(&mtx_1_r);
   printf("Matrix 2\n\r");
   mtx_print(&mtx_2_r);
   printf("Matrix 1 + 2\n\r");
   mtx_print(&mtx_3_r);
   printf("Matrix 1 * 5.0\n\r");
   mtx_print(&mtx_4_r);
}
