#include "att_det.h"
#include "mtx.h"
#include <stdio.h>
#include <math.h>

/*
 * Function: quat_2_dcm
 * ------------------------
 *  Converts a 4 element quaternion
 *  array to a 3x3 DCM array
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
 *  result: quaternion
 */

void dcm_2_q(struct mtx_matrix* dcm, struct mtx_matrix* q){
	int max_qt_row = 0;
    int max_qt_col = 0;
	float max_qt_val = 0;
    struct mtx_matrix qt;

	/*Find Maximum value element in q temp*/
    mtx_create_ones(4,1,&qt);
    mtx_set(1,1,&qt,  .25*(1+2*mtx_get(1,1,dcm)-mtx_trace(dcm)));
    mtx_set(2,1,&qt,  .25*(1+2*mtx_get(2,2,dcm)-mtx_trace(dcm)));
    mtx_set(3,1,&qt,  .25*(1+2*mtx_get(3,3,dcm)-mtx_trace(dcm)));
    mtx_set(4,1,&qt,  .25*(1+mtx_trace(dcm)));
    max_qt_val = mtx_max(&qt,&max_qt_row,&max_qt_col);

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

//void body_rate_dcm_rot(float* body_rates, float* prior_dcm, float* rot_dcm){
//	float timestep = 1.0/UPDATE_RATE;
//    struct mtx_matrix w_ss_mtrx;
//    struct mtx_matrix prior_dcm_mtrx;
//    struct mtx_matrix diff_dcm_mtrx;
//    struct mtx_matrix diff_dcm_mtrx_scaled;
//    struct mtx_matrix sum_dcm_mtrx;
//
//	float w_ss[] = {0            ,-body_rates[2], body_rates[1],
//              	    body_rates[2],0             ,-body_rates[0],
//              	    -body_rates[1], body_rates[0],0             };
//
//    mtx_create(3,3,w_ss,&w_ss_mtrx);
//    mtx_create(3,3,prior_dcm,&prior_dcm_mtrx);
//    mtx_create_ones(3,3,&diff_dcm_mtrx);
//    mtx_create_ones(3,3,&diff_dcm_mtrx_scaled);
//    mtx_create_ones(3,3,&sum_dcm_mtrx);
//
//    mtx_mult(&w_ss_mtrx,&prior_dcm_mtrx,&diff_dcm_mtrx);
//    mtx_scale(&diff_dcm_mtrx,timestep,&diff_dcm_mtrx_scaled);
//    mtx_sum(&diff_dcm_mtrx_scaled,&prior_dcm_mtrx);
//	
//	rot_dcm[0] = mtx_get(0,0,&sum_dcm_mtrx);
//	rot_dcm[1] = mtx_get(0,1,&sum_dcm_mtrx);
//	rot_dcm[2] = mtx_get(0,2,&sum_dcm_mtrx);
//	rot_dcm[3] = mtx_get(1,0,&sum_dcm_mtrx);
//	rot_dcm[4] = mtx_get(1,1,&sum_dcm_mtrx);
//	rot_dcm[5] = mtx_get(1,1,&sum_dcm_mtrx);
//	rot_dcm[6] = mtx_get(2,0,&sum_dcm_mtrx);
//	rot_dcm[7] = mtx_get(2,1,&sum_dcm_mtrx);
//	rot_dcm[8] = mtx_get(2,2,&sum_dcm_mtrx);
//	
//	return;
//}
//
/*
 * Function: est_sun_vec_ls
 * ------------------------
 *  Estimates a Sun pointing vector 
 *  based on sun sensor readings using
 *  a Least Squares algorithm
 *
 *  result: Sun vector in the body frame
 */
//
//int est_sun_vec_ls(gsl_matrix* sun_sens_volt_mtrx_p, gsl_matrix* sun_sens_norm_mtrx_p, gsl_matrix* sun_vec_mtrx_p){
//
//	/* Initialize Matrices */
//	int sun_signum;
//	gsl_permutation* sun_perm = gsl_permutation_alloc(3);
//	gsl_matrix* sun_right_mtrx_p = gsl_matrix_alloc(3,1);
//	gsl_matrix* sun_left_mtrx_p = gsl_matrix_alloc(3,3);
//	gsl_matrix* sun_left_inverse_mtrx_p = gsl_matrix_alloc(3,3);
//
//	int i;
//	int ss_count = 0;
//	for (i=0;i<9;i++){
//		if(sun_sens_volt_mtrx_p->data[i]>SS_V_CUTOFF){
//			ss_count++;
//		}
//	}
//
//	if (ss_count < 3) { 
//		for (i=0;i<3;i++){
//			sun_vec_mtrx_p->data[i] = 0;
//		}
//		return 1;
//	}
//
//	/*Calculate Least Squares solution - inv(norm'*norm)*norm'*sun_sens*/
//	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, sun_sens_norm_mtrx_p, sun_sens_volt_mtrx_p, 0.0, sun_right_mtrx_p);
//	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, sun_sens_norm_mtrx_p, sun_sens_norm_mtrx_p, 0.0, sun_left_mtrx_p);
//	gsl_linalg_LU_decomp(sun_left_mtrx_p,sun_perm,&sun_signum);
//	gsl_linalg_LU_invert(sun_left_mtrx_p,sun_perm,sun_left_inverse_mtrx_p);
//        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, sun_left_inverse_mtrx_p, sun_right_mtrx_p, 0.0, sun_vec_mtrx_p);
//
//	/*Cleanup*/
//	gsl_matrix_free(sun_right_mtrx_p);
//	gsl_matrix_free(sun_left_mtrx_p);
//	gsl_matrix_free(sun_left_inverse_mtrx_p);
//	return 0;
//}
//
