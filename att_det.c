#include <stdio.h>
#include <math.h>
#include "mtx.h"
#include "att_det.h"

/*
 * Function: quat_2_dcm
 * ------------------------
 *  Converts a 4 element quaternion
 *  array to a 3x3 DCM array
 *
 *  input q - 4x1 mtx_matrix
 *  input dcm - 3x3 mtx_matrix
 *
 *  result: DCM
 */

void q_2_dcm(struct mtx_matrix* q, struct mtx_matrix* dcm){
	mtx_set(1,1,dcm, 	
        mtx_get(1,1,q)*mtx_get(1,1,q)
            -mtx_get(2,1,q)*mtx_get(2,1,q)
            -mtx_get(3,1,q)*mtx_get(3,1,q)
            +mtx_get(4,1,q)*mtx_get(4,1,q)); 
    mtx_set(1,2,dcm,                                
        2*(mtx_get(1,1,q)*mtx_get(2,1,q)+mtx_get(3,1,q)*mtx_get(4,1,q)));
    mtx_set(1,3,dcm,                                
        2*(mtx_get(1,1,q)*mtx_get(3,1,q)-mtx_get(2,1,q)*mtx_get(4,1,q)));
    mtx_set(2,1,dcm,                 		    
        2*(mtx_get(1,1,q)*mtx_get(2,1,q)-mtx_get(3,1,q)*mtx_get(4,1,q)) );
    mtx_set(2,2,dcm,                                
        -1*mtx_get(1,1,q)*mtx_get(1,1,q)
            +mtx_get(2,1,q)*mtx_get(2,1,q)
            -mtx_get(3,1,q)*mtx_get(3,1,q)
            +mtx_get(4,1,q)*mtx_get(4,1,q));
    mtx_set(2,3,dcm,                                
        2*(mtx_get(2,1,q)*mtx_get(3,1,q)+mtx_get(1,1,q)*mtx_get(4,1,q)));
    mtx_set(3,1,dcm,                 		    
        2*(mtx_get(1,1,q)*mtx_get(3,1,q)+mtx_get(2,1,q)*mtx_get(4,1,q)));
    mtx_set(3,2,dcm,                                
        2*(mtx_get(2,1,q)*mtx_get(3,1,q)-mtx_get(1,1,q)*mtx_get(4,1,q))); 
    mtx_set(3,3,dcm,                                
        -1*mtx_get(1,1,q)*mtx_get(1,1,q)
            -mtx_get(2,1,q)*mtx_get(2,1,q)
            +mtx_get(3,1,q)*mtx_get(3,1,q)
            +mtx_get(4,1,q)*mtx_get(4,1,q));
}

/*
 * Function: dcm_2_q
 * ------------------------
 *  Converts a 3x3 DCM array
 *  to a 4 element quaternion array
 *
 *  input dcm - 3x3 mtx_matrix
 *  input q - 4x1 mtx_matrix
 *
 *  result: quaternion
 */

void dcm_2_q(struct mtx_matrix* dcm, struct mtx_matrix* q){
	int max_qt_row;
    int max_qt_col;
	float max_qt_val = 0;
    struct mtx_matrix qt;

    max_qt_row = 0;
    max_qt_col = 0;

	/*Find Maximum value element in q temp*/
    mtx_create_ones(4,1,&qt);
    mtx_set(1,1,&qt,  .25*(1+2*mtx_get(1,1,dcm)-mtx_trace(dcm)));
    mtx_set(2,1,&qt,  .25*(1+2*mtx_get(2,2,dcm)-mtx_trace(dcm)));
    mtx_set(3,1,&qt,  .25*(1+2*mtx_get(3,3,dcm)-mtx_trace(dcm)));
    mtx_set(4,1,&qt,  .25*(1+mtx_trace(dcm)));
    max_qt_val = mtx_max(&qt,&max_qt_row,&max_qt_col);
    printf("Max Val = %f Max Row: %d\n\r",max_qt_val, max_qt_row);

	switch(max_qt_row){
		case 1 :
			mtx_set(1,1,q,  
                sqrt(mtx_get(1,1,&qt)));
            mtx_set(2,1,q,  
                (mtx_get(1,2,dcm)+mtx_get(2,1,dcm))/4/mtx_get(1,1,q));   
            mtx_set(3,1,q,  
                (mtx_get(1,3,dcm)+mtx_get(3,1,dcm))/4/mtx_get(1,1,q));   
            mtx_set(4,1,q, 
                (mtx_get(2,3,dcm)-mtx_get(3,2,dcm))/4/mtx_get(1,1,q));   
		case 2 :
			mtx_set(2,1,q,  
                sqrt(mtx_get(2,1,&qt)));
            mtx_set(1,1,q,  
                (mtx_get(1,2,dcm)+mtx_get(2,1,dcm))/4/mtx_get(2,1,q));   
            mtx_set(3,1,q,  
                (mtx_get(2,3,dcm)+mtx_get(3,2,dcm))/4/mtx_get(2,1,q));   
            mtx_set(4,1,q, 
                (mtx_get(3,1,dcm)-mtx_get(1,3,dcm))/4/mtx_get(2,1,q));   
		case 3 :
			mtx_set(3,1,q,  
                sqrt(mtx_get(3,1,&qt)));
            mtx_set(1,1,q,  
                (mtx_get(1,3,dcm)+mtx_get(3,1,dcm))/4/mtx_get(3,1,q));   
            mtx_set(2,1,q,  
                (mtx_get(2,3,dcm)+mtx_get(3,2,dcm))/4/mtx_get(3,1,q));   
            mtx_set(4,1,q, 
                (mtx_get(1,2,dcm)-mtx_get(2,1,dcm))/4/mtx_get(3,1,q));   
		case 4 :
			mtx_set(4,1,q,  
                sqrt(mtx_get(4,1,&qt)));
            mtx_set(1,1,q,  
                (mtx_get(2,3,dcm)-mtx_get(3,2,dcm))/4/mtx_get(4,1,q));   
            mtx_set(2,1,q,  
                (mtx_get(3,1,dcm)-mtx_get(1,3,dcm))/4/mtx_get(4,1,q));   
            mtx_set(3,1,q, 
                (mtx_get(1,2,dcm)-mtx_get(2,1,dcm))/4/mtx_get(4,1,q));   
	}
}

/*
 * Function: body_rate_dcm_rot
 * ------------------------
 *  Rotates a DCM array by a
 *  body rate vector array
 *
 *  result: DCM
 */

void body_rate_dcm_rot(struct mtx_matrix* body_rates, struct mtx_matrix* prior_dcm, 
                        struct mtx_matrix* rot_dcm){
	float timestep = 1.0/UPDATE_RATE;

    struct mtx_matrix w_ss;
    struct mtx_matrix diff_dcm;
    struct mtx_matrix diff_dcm_scaled;

    mtx_create_ones(3,3,&w_ss);
    mtx_create_ones(3,3,&diff_dcm);
    mtx_create_ones(3,3,&diff_dcm_scaled);

    mtx_ss(body_rates, &w_ss);
    mtx_mult(&w_ss,prior_dcm,&diff_dcm);
    mtx_scale(&diff_dcm,timestep,&diff_dcm_scaled);
    mtx_sum(&diff_dcm_scaled,prior_dcm,rot_dcm);
	return;
}

/*
 * Function: est_sun_vec_ls
 * ------------------------
 *  Estimates a Sun pointing vector 
 *  based on sun sensor readings using
 *  a Least Squares algorithm
 *
 *  result: Sun vector in the body frame
 */

int est_sun_vec_ls(struct mtx_matrix* sun_sens_volt, 
                    struct mtx_matrix* sun_sens_norm, 
                    struct mtx_matrix* sun_vec){

	int i;
	int ss_count = 0;
    struct mtx_matrix norm_trans;
    struct mtx_matrix sun_right;
    struct mtx_matrix sun_left;
    struct mtx_matrix sun_left_inverse;

	/* Initialize Matrices */
    mtx_create_ones(3,1,&sun_right);
    mtx_create_ones(3,3,&sun_left);
    mtx_create_ones(3,3,&sun_left_inverse);
    mtx_create_ones(3,SS_COUNT,&norm_trans);

    /* Check Sun Sensors against threshold voltage*/
	for (i=0;i<SS_COUNT;i++){
		if(mtx_get(i,1,sun_sens_volt)>SS_V_CUTOFF) ss_count++;
	}

	if (ss_count < 3) { 
		for (i=0;i<3;i++){
			mtx_set(1,1,sun_vec,0);
			mtx_set(2,1,sun_vec,0);
			mtx_set(3,1,sun_vec,0);
		}
		return 1;
	}

	/*Calculate Least Squares solution - inv(norm'*norm)*norm'*sun_sens*/
    mtx_trans(sun_sens_norm,&norm_trans)==0;
    mtx_print(&norm_trans);
    mtx_mult(&norm_trans, sun_sens_volt, &sun_right);
    mtx_print(&sun_right);
    mtx_mult(&norm_trans, sun_sens_norm, &sun_left);
    mtx_print(&sun_left);
    mtx_inv(&sun_left, &sun_left_inverse);
    mtx_print(&sun_left_inverse);
    mtx_mult(&sun_left_inverse, &sun_right, sun_vec);

	return 0;
}

