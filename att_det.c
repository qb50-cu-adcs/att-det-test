#include "att_det.h"
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/*
 * Function: est_sun_vec_ls
 * ------------------------
 *  Estimates a Sun pointing vector 
 *  based on sun sensor readings using
 *  a Least Squares algorithm
 *
 *  result: Sun vector in the body frame
 */

int est_sun_vec_ls(gsl_matrix* sun_sens_volt_mtrx_p, gsl_matrix* sun_sens_norm_mtrx_p, gsl_matrix* sun_vec_mtrx_p){

	/* Initialize Matrices */
	int sun_signum;
	gsl_permutation* sun_perm = gsl_permutation_alloc(3);
	gsl_matrix* sun_right_mtrx_p = gsl_matrix_alloc(3,1);
	gsl_matrix* sun_left_mtrx_p = gsl_matrix_alloc(3,3);
	gsl_matrix* sun_left_inverse_mtrx_p = gsl_matrix_alloc(3,3);

	int i;
	int ss_count = 0;
	for (i=0;i<9;i++){
		if(sun_sens_volt_mtrx_p->data[i]>SS_V_CUTOFF){
			ss_count++;
		}
	}

	if (ss_count < 3) { 
		for (i=0;i<3;i++){
			sun_vec_mtrx_p->data[i] = 0;
		}
		return 1;
	}

	/*Calculate Least Squares solution - inv(norm'*norm)*norm'*sun_sens*/
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, sun_sens_norm_mtrx_p, sun_sens_volt_mtrx_p, 0.0, sun_right_mtrx_p);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, sun_sens_norm_mtrx_p, sun_sens_norm_mtrx_p, 0.0, sun_left_mtrx_p);
	gsl_linalg_LU_decomp(sun_left_mtrx_p,sun_perm,&sun_signum);
	gsl_linalg_LU_invert(sun_left_mtrx_p,sun_perm,sun_left_inverse_mtrx_p);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, sun_left_inverse_mtrx_p, sun_right_mtrx_p, 0.0, sun_vec_mtrx_p);

	/*Cleanup*/
	gsl_matrix_free(sun_right_mtrx_p);
	gsl_matrix_free(sun_left_mtrx_p);
	gsl_matrix_free(sun_left_inverse_mtrx_p);
	return 0;
}

/*
 * Function: est_svd_2v
 * -----------------
 *  Calculates an attitude estimate based on
 *  the SVD solution to Whaba's problem using
 *  2 Body Vectors and 2 ECI vectors.
 *
 *  result: updated est_state 
 */
void est_svd_2v(gsl_matrix* body_1_vec_mtrx_p,gsl_matrix* body_2_vec_mtrx_p,gsl_matrix* eci_1_vec_mtrx_p,gsl_matrix* eci_2_vec_mtrx_p, struct est_state* state){

	/*Initialize Matrices*/
	int est_signum;
	gsl_matrix* b_mtrx_p = gsl_matrix_alloc(3,3);
	gsl_matrix_set_all(b_mtrx_p,0);
	gsl_matrix* u_lu_mtrx_p = gsl_matrix_alloc(3,3);
	gsl_matrix* v_mtrx_p = gsl_matrix_alloc(3,3);
	gsl_matrix* v_lu_mtrx_p = gsl_matrix_alloc(3,3);
	gsl_matrix* dcm_right_mtrx_p = gsl_matrix_alloc(3,3);
	gsl_vector* s_vec_p= gsl_vector_alloc(3);
	gsl_vector* work_vec_p = gsl_vector_alloc(3);
	gsl_permutation* est_perm = gsl_permutation_alloc(3);

	/*Start SVD Wahba Solution 
          Calc B = B + eci'*body*/
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, body_1_vec_mtrx_p, eci_1_vec_mtrx_p, 1.0, b_mtrx_p);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, body_2_vec_mtrx_p, eci_2_vec_mtrx_p, 1.0, b_mtrx_p);

	/*Calc SVD*/
	gsl_linalg_SV_decomp(b_mtrx_p, v_mtrx_p, s_vec_p, work_vec_p);

	/*Calc M = diag([1,1,det(U)*det(V)]);*/
	gsl_matrix_memcpy(u_lu_mtrx_p, b_mtrx_p);
	gsl_linalg_LU_decomp(u_lu_mtrx_p,est_perm,&est_signum);
	double u_det = gsl_linalg_LU_det(u_lu_mtrx_p,est_signum);
	gsl_matrix_memcpy(v_lu_mtrx_p, v_mtrx_p);
	gsl_linalg_LU_decomp(v_lu_mtrx_p,est_perm,&est_signum);
	double v_det = gsl_linalg_LU_det(v_lu_mtrx_p,est_signum);
	double m[] = {1,0,0,
                      0,1,0,
                      0,0,u_det*v_det};
	gsl_matrix_view m_mtrx = gsl_matrix_view_array(m,3,3);

	/*Calc DCM = U*M*V'*/
	double dcm[] = {0,0,0,
                      0,0,0,
                      0,0,0};
	gsl_matrix_view dcm_mtrx = gsl_matrix_view_array(dcm,3,3);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &m_mtrx.matrix, v_mtrx_p, 0.0, dcm_right_mtrx_p);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, b_mtrx_p, dcm_right_mtrx_p, 0.0, &dcm_mtrx.matrix);

	printf("DCM = [%.2f, %.2f, %.2f\n     %.2f,%.2f,%.2f\n     %.2f,%.2f,%.2f]\n",
               dcm[0], dcm[1], dcm[2], dcm[3], dcm[4], dcm[5], dcm[6], dcm[7], dcm[8]);

	/*Convert DCM to Quaternion*/

	/*Assign Output - TEMP*/
	state->att_quaternion[0] = 0.1;
	state->att_quaternion[1] = 0.2;
	state->att_quaternion[2] = 0.3;
	state->att_quaternion[3] = 1.4;

	/*Cleanup	*/
	gsl_matrix_free(b_mtrx_p);
	gsl_matrix_free(u_lu_mtrx_p);
	gsl_matrix_free(v_mtrx_p);
	gsl_matrix_free(v_lu_mtrx_p);
	gsl_matrix_free(dcm_right_mtrx_p);
	gsl_vector_free(s_vec_p);
	gsl_vector_free(work_vec_p);
	return;
}

/*
 * Function: est_svd
 * -----------------
 *  Calculates an attitude estimate based on
 *  the SVD solution to Whaba's problem.
 *  
 *  The body frame sun vector is calculated 
 *  independantly by a Least Squares estimator
 *
 *  If the SC is in eclipse, the body rates are 
 *  integrated in place of a sun vector
 *
 *  result: updated est_state 
 */

void est_svd(double* gyro_body_rates, double* mag_vec, double* sun_sens_volt, double* sun_sens_norm, double* sun_eci_vec, double* mag_eci_vec, struct est_state* state){

	/*Initialize Matrices*/
	gsl_matrix* sun_vec_mtrx_p = gsl_matrix_alloc(3,1);

	/*Vector to Matrix Conversion*/
	gsl_matrix_view mag_vec_mtrx = gsl_matrix_view_array(mag_vec,3,1);
	gsl_matrix_view sun_eci_vec_mtrx = gsl_matrix_view_array(sun_eci_vec,3,1);
	gsl_matrix_view mag_eci_vec_mtrx = gsl_matrix_view_array(mag_eci_vec,3,1);
	gsl_matrix_view sun_sens_volt_mtrx = gsl_matrix_view_array(sun_sens_volt,SS_COUNT,1);
	gsl_matrix_view sun_sens_norm_mtrx = gsl_matrix_view_array(sun_sens_norm,SS_COUNT,3);

	/*Calc Sun Vector*/ 
	int sv_ret;
	sv_ret = est_sun_vec_ls(&sun_sens_volt_mtrx.matrix, &sun_sens_norm_mtrx.matrix, sun_vec_mtrx_p);
	printf("Sun Vector:\n");
	gsl_matrix_fprintf(stdout, sun_vec_mtrx_p,"%.2f");

	/*Check if in Eclipse*/
	if (sv_ret) { 
        //Eclipse
	 	printf("ECLIPSE: %d\n",sv_ret);
		//Integrate Rate Gyros for 2nd Vector
		//Calculate Estimate
		est_svd_2v(&mag_eci_vec_mtrx.matrix,&sun_eci_vec_mtrx.matrix,&mag_vec_mtrx.matrix,sun_vec_mtrx_p,state);
	} else {
	//Sun
	 	printf("DAY: %d\n",sv_ret);
		//Calculate Estimate
		est_svd_2v(&mag_eci_vec_mtrx.matrix,&sun_eci_vec_mtrx.matrix,&mag_vec_mtrx.matrix,sun_vec_mtrx_p,state);
	}

	/*Assign output Body Rates -- TEMP*/
	state->sat_body_rates[0] = gyro_body_rates[0];
	state->sat_body_rates[1] = gyro_body_rates[1];
	state->sat_body_rates[2] = gyro_body_rates[2];

	/*Cleanup	*/
	gsl_matrix_free(sun_vec_mtrx_p);
	return;
}

