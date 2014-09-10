#include "att_det.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/*
 * Function: quat_2_dcm
 * ------------------------
 *  Converts a 4 element quaternion
 *  array to a 3x3 DCM array
 *
 *  result: DCM
 */

void q_2_dcm(double* q, double* dcm){
	double dcm_new[] = { 	q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3], 2*(q[0]*q[1]+q[2]*q[3]), 2*(q[0]*q[2]-q[1]*q[3]),
                 		2*(q[0]*q[1]-q[2]*q[3]), -q[0]*q[0]+q[1]*q[1]-q[2]*q[2]+q[3]*q[3], 2*(q[1]*q[2]+q[0]*q[3]),
                 		2*(q[0]*q[2]+q[1]*q[3]), 2*(q[1]*q[2]-q[0]*q[3]), -q[0]*q[0]-q[1]*q[1]+q[2]*q[2]+q[3]*q[3] };
	memcpy(dcm,&dcm_new,9*sizeof(double));
}

/*
 * Function: dcm_2_q
 * ------------------------
 *  Converts a 3x3 DCM array
 *  to a 4 element quaternion array
 *
 *  result: quaternion
 */

void dcm_2_q(double* dcm, double* q){
	int max_q_idx = 0;
	double max_q_val = 0;

	/*Find Maximum element in q*/
	int i;
	for(i=0;i<4;i++){
		if (q[i]>max_q_val || -q[i]>max_q_val){
			max_q_val = q[i];
			max_q_idx = i;
		}	
	};

	switch(max_q_idx){
		case 0 :
			q[0] = sqrt((1.0/4.0)*(1+2*dcm[0]-(dcm[0]+dcm[4]+dcm[8])));
			q[1] = ((dcm[1]+dcm[3])/4.0)/q[0];
			q[2] = ((dcm[2]+dcm[6])/4.0)/q[0];
			q[3] = ((dcm[5]-dcm[7])/4.0)/q[0];
		case 1 :
			q[1] = sqrt((1.0/4.0)*(1+2*dcm[1]-(dcm[0]+dcm[4]+dcm[8])));
			q[0] = ((dcm[1]+dcm[3])/4.0)/q[1];
			q[2] = ((dcm[5]+dcm[7])/4.0)/q[1];
			q[3] = ((dcm[6]-dcm[2])/4.0)/q[1];
		case 2 :
			q[2] = sqrt((1.0/4.0)*(1+2*dcm[2]-(dcm[0]+dcm[4]+dcm[8])));
			q[0] = ((dcm[2]+dcm[6])/4.0)/q[2];
			q[1] = ((dcm[5]+dcm[7])/4.0)/q[2];
			q[3] = ((dcm[1]-dcm[3])/4.0)/q[2];
		case 3 :
			q[3] = sqrt((1.0/4.0)*(1+(dcm[0]+dcm[4]+dcm[8])));
			q[0] = ((dcm[5]-dcm[7])/4.0)/q[3];
			q[1] = ((dcm[6]-dcm[2])/4.0)/q[3];
			q[2] = ((dcm[1]-dcm[3])/4.0)/q[3];
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

void body_rate_dcm_rot(double* body_rates, double* prior_dcm, double* rot_dcm){
	double timestep = 1.0/UPDATE_RATE;
	gsl_matrix* diff_dcm_mtrx_p = gsl_matrix_alloc(3,3);

	double w_ss[] = { 0            ,-body_rates[2], body_rates[1],
              		  body_rates[2],0             ,-body_rates[0],
              	         -body_rates[1], body_rates[0],0             };
        gsl_matrix_view w_ss_mtrx = gsl_matrix_view_array(w_ss,3,3);

        gsl_matrix_view prior_dcm_mtrx = gsl_matrix_view_array(prior_dcm,3,3);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &w_ss_mtrx.matrix, &prior_dcm_mtrx.matrix, 0.0, diff_dcm_mtrx_p);
	gsl_matrix_scale(diff_dcm_mtrx_p,timestep);
	gsl_matrix_add(diff_dcm_mtrx_p,&prior_dcm_mtrx.matrix);
	
	rot_dcm[0] = gsl_matrix_get(diff_dcm_mtrx_p,0,0);
	rot_dcm[1] = gsl_matrix_get(diff_dcm_mtrx_p,0,1);
	rot_dcm[2] = gsl_matrix_get(diff_dcm_mtrx_p,0,2);
	rot_dcm[3] = gsl_matrix_get(diff_dcm_mtrx_p,1,0);
	rot_dcm[4] = gsl_matrix_get(diff_dcm_mtrx_p,1,1);
	rot_dcm[5] = gsl_matrix_get(diff_dcm_mtrx_p,1,1);
	rot_dcm[6] = gsl_matrix_get(diff_dcm_mtrx_p,2,0);
	rot_dcm[7] = gsl_matrix_get(diff_dcm_mtrx_p,2,1);
	rot_dcm[8] = gsl_matrix_get(diff_dcm_mtrx_p,2,2);
	
	/*Cleanup*/
	gsl_matrix_free(diff_dcm_mtrx_p);
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

	/* Start SVD Wahba Solution */
        /* Calc B = B + eci'*body*/
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
	double dcm_plus[] = {0,0,0,
                             0,0,0,
                             0,0,0};
	gsl_matrix_view dcm_mtrx = gsl_matrix_view_array(dcm_plus,3,3);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &m_mtrx.matrix, v_mtrx_p, 0.0, dcm_right_mtrx_p);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, b_mtrx_p, dcm_right_mtrx_p, 0.0, &dcm_mtrx.matrix);

	printf("DCM = [%.4f, %.4f, %.4f\n     %.4f,%.4f,%.4f\n     %.4f,%.4f,%.4f]\n",
               dcm_plus[0], dcm_plus[1], dcm_plus[2], dcm_plus[3], dcm_plus[4], dcm_plus[5], dcm_plus[6], dcm_plus[7], dcm_plus[8]);

	/*Convert DCM to Quaternion*/
	double q_plus[4];
	dcm_2_q(dcm_plus,q_plus);

	/*Assign Output*/
	state->att_quaternion[0] = q_plus[0];
	state->att_quaternion[1] = q_plus[1];
	state->att_quaternion[2] = q_plus[2];
	state->att_quaternion[3] = q_plus[3];

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
	int in_eclipse;
	in_eclipse = est_sun_vec_ls(&sun_sens_volt_mtrx.matrix, &sun_sens_norm_mtrx.matrix, sun_vec_mtrx_p);
	printf("Sun Vector:\n");
	gsl_matrix_fprintf(stdout, sun_vec_mtrx_p,"%.2f");

	/*Check if in Eclipse*/
	if (in_eclipse) { 
		//Integrate Rate Gyros
		double prior_dcm[9];
		double updated_dcm[9];
		q_2_dcm(state->att_quaternion,prior_dcm);
		body_rate_dcm_rot(gyro_body_rates, prior_dcm, updated_dcm);

		//Pull reference vector from DCM
		double x_intg[] = {  updated_dcm[0],
				     updated_dcm[3],
				     updated_dcm[6]};
	        gsl_matrix_view x_intg_vec_mtrx = gsl_matrix_view_array(x_intg,3,1);

		double x_eci[] = {  1,
				    0,
				    0};
	        gsl_matrix_view x_eci_vec_mtrx = gsl_matrix_view_array(x_eci,3,1);

		//Calculate State Estimate - Magnetometer + reference vector from rate gyro integration
		est_svd_2v(&mag_eci_vec_mtrx.matrix,&x_intg_vec_mtrx.matrix,&mag_vec_mtrx.matrix,&x_eci_vec_mtrx.matrix,state);

	} else {
		//Calculate State Estimate - Magnetometer and Sun Vector
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

