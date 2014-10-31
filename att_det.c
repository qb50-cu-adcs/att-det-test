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
 * Function: cal_sun_sens
 * ------------------------
 *  calibrate sun sensors, remove albedo effects
 *
 *  result: mtx_matrix of eci->body frame rotation
 */
int cal_sun_sens( struct mtx_matrix* sun_sens_volt, 
                  struct mtx_matrix* sun_sens_norm, 
                  struct mtx_matrix* sun_sens_cal){

    int i;
    int j;
    float norm_sq_val;
    float min_sum;
    float min_sum_row;
    struct mtx_matrix cur_sens;
    struct mtx_matrix cur_sum;
    struct mtx_matrix tmp_133a;
    struct mtx_matrix tmp_131a;
    mtx_create_ones(1,3,&cur_sens);
    mtx_create_ones(SS_COUNT,3,&tmp_133a);
    mtx_create_ones(SS_COUNT,1,&tmp_131a);
    mtx_create_ones(SS_COUNT,3,&cur_sum);
    
    /*Compare sun sensor against sun sensor on opposite side of craft*/
    for (i=1;i<=SS_COUNT;i++){
        mtx_set(i,1,&cur_sens,   mtx_get(i,1,sun_sens_norm));
        mtx_set(i,2,&cur_sens,   mtx_get(i,2,sun_sens_norm));
        mtx_set(i,3,&cur_sens,   mtx_get(i,3,sun_sens_norm));
        /*Add currect sensor to all norm vectors and find min val*/
        mtx_mult(&tmp_131a,&cur_sens,&tmp_133a); 
        mtx_sum(&tmp_133a,sun_sens_norm,&cur_sum);
        for (j=1;j<=SS_COUNT;j++){
            norm_sq_val = mtx_get(j,1,&cur_sum)*mtx_get(j,1,&cur_sum)
                            + mtx_get(j,2,&cur_sum)*mtx_get(j,2,&cur_sum)
                            + mtx_get(j,3,&cur_sum)*mtx_get(j,3,&cur_sum);
            if(j==1 || norm_sq_val<min_sum) {
                min_sum = norm_sq_val;
                min_sum_row = j;
            }
        }
        /*If voltage is less than voltage at opposite location, set to 0*/
        if(mtx_get(i,1,sun_sens_volt)>mtx_get(min_sum_row,1,sun_sens_volt) &&
            mtx_get(i,1,sun_sens_volt)>SS_V_CUTOFF){
            mtx_set(i,1,sun_sens_cal, mtx_get(i,1,sun_sens_volt));
        } else {
            mtx_set(i,1,sun_sens_cal, 0);
        }
    }
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
    struct mtx_matrix sun_sens_cal;

	/* Initialize Matrices */
    mtx_create_ones(3,1,&sun_right);
    mtx_create_ones(3,3,&sun_left);
    mtx_create_ones(3,3,&sun_left_inverse);
    mtx_create_ones(3,SS_COUNT,&norm_trans);
    mtx_create_ones(SS_COUNT,1,&sun_sens_cal);

    cal_sun_sens(sun_sens_volt,sun_sens_norm,&sun_sens_cal);

    /* Check Sun Sensors against threshold voltage*/
	for (i=0;i<SS_COUNT;i++){
		if(mtx_get(i,1,&sun_sens_cal)>SS_V_CUTOFF) ss_count++;
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
    mtx_mult(&norm_trans, &sun_sens_cal, &sun_right);
    mtx_mult(&norm_trans, sun_sens_norm, &sun_left);
    mtx_inv(&sun_left, &sun_left_inverse);
    mtx_mult(&sun_left_inverse, &sun_right, sun_vec);

	return 0;
}

/*
 * Function: est_quest
 * ------------------------
 *  Estimates a DCM matching the rotation between the observed 
 *  body frame vectors and the reference eci vectors using the 
 *  QUEST method. 
 *  Rodrigues parameters are used for Eigenvector estimate
 *  returns -1 when rotation = 180 deg.
 *
 *  result: mtx_matrix of eci->body frame rotation
 */

int est_quest_rp(struct mtx_matrix* b_k, struct mtx_matrix* eci_k, 
                    struct mtx_matrix* dcm_out){

    struct mtx_matrix q_out;
    struct mtx_matrix b;
    struct mtx_matrix s;
    struct mtx_matrix z;
    struct mtx_matrix tmp_33a;
    struct mtx_matrix tmp_33b;
    struct mtx_matrix eye_33;
    struct mtx_matrix m_inv;
    struct mtx_matrix norm_obs;
    struct mtx_matrix norm_eci;
    struct mtx_matrix tmp_13a;
    struct mtx_matrix p;
    struct mtx_matrix q_scale_sub;
    float sigma;
    float a_i;

    int k;
    float norm_obs_in;
    float norm_eci_in;
    float q_scale;

    /*Eigenvalue of K ~= 1*/
    float eig_val = 1.0;

    mtx_create_ones(3,1,&norm_obs);
    mtx_create_ones(3,1,&norm_eci);
    mtx_create_ones(4,1,&q_out);
    mtx_create_val(3,3,&b,0);
    mtx_create_ones(3,3,&s);
    mtx_create_ones(3,1,&z);
    mtx_create_eye(3,3,&eye_33);
    mtx_create_ones(3,3,&m_inv);
    mtx_create_ones(3,1,&p);
    mtx_create_ones(1,1,&q_scale_sub);
    mtx_create_ones(1,3,&tmp_13a);
    mtx_create_val(3,3,&tmp_33a,0);
    mtx_create_val(3,3,&tmp_33b,0);

    /*B = SUM(a_i*|b_k|*|eci_k|')*/
    a_i = (1.0/b_k->cols);
    for (k=1;k<=b_k->cols;k++){
        norm_obs_in = sqrt(mtx_get(1,k,b_k)*mtx_get(1,k,b_k) 
                            + mtx_get(2,k,b_k)*mtx_get(2,k,b_k)
                            + mtx_get(3,k,b_k)*mtx_get(3,k,b_k)); 
        norm_eci_in = sqrt(mtx_get(1,k,eci_k)*mtx_get(1,k,eci_k) 
                            + mtx_get(2,k,eci_k)*mtx_get(2,k,eci_k)
                            + mtx_get(3,k,eci_k)*mtx_get(3,k,eci_k)); 
        mtx_set(1,1,&norm_obs,   mtx_get(1,k,b_k)/norm_obs_in); 
        mtx_set(2,1,&norm_obs,   mtx_get(2,k,b_k)/norm_obs_in); 
        mtx_set(3,1,&norm_obs,   mtx_get(3,k,b_k)/norm_obs_in); 
        mtx_set(1,1,&norm_eci,   mtx_get(1,k,eci_k)/norm_eci_in); 
        mtx_set(2,1,&norm_eci,   mtx_get(2,k,eci_k)/norm_eci_in); 
        mtx_set(3,1,&norm_eci,   mtx_get(3,k,eci_k)/norm_eci_in); 
        mtx_trans(&norm_eci,&tmp_13a);
        mtx_mult(&norm_obs,&tmp_13a,&tmp_33a);
        mtx_scale(&tmp_33a,a_i,&tmp_33b);
        mtx_sum(&b,&tmp_33b,&tmp_33a);
        mtx_copy(&tmp_33a,&b);
    }

    /*sigma, S, Z*/
    sigma = mtx_trace(&b);
    mtx_trans(&b,&tmp_33a);
    mtx_sum(&b,&tmp_33a,&s);
    mtx_set(1,1,&z, mtx_get(2,3,&b)-mtx_get(3,2,&b));
    mtx_set(2,1,&z, mtx_get(3,1,&b)-mtx_get(1,3,&b));
    mtx_set(3,1,&z, mtx_get(1,2,&b)-mtx_get(2,1,&b));

    /*M_inv = (eig_val + sigma)*I - S*/
    mtx_scale(&eye_33, eig_val + sigma, &tmp_33a);
    mtx_scale(&s, -1, &tmp_33b);
    mtx_sum(&tmp_33a, &tmp_33b, &m_inv);

    /*Check if matrix is non-invertable (rot=180)*/
    if (mtx_det(&m_inv)<MIN_MINV_DET) return -1;

    /* p = inv(M_inv)*Z*
     * q = 1/sqrt(1+p'*p)*[p;1] */
    mtx_inv(&m_inv,&tmp_33a);
    mtx_mult(&tmp_33a,&z,&p);
    mtx_trans(&p,&tmp_13a);
    mtx_mult(&tmp_13a,&p,&q_scale_sub);
    q_scale = 1.0/sqrt(1.0+mtx_get(1,1,&q_scale_sub));
    mtx_set(1,1,&q_out, q_scale*mtx_get(1,1,&p));
    mtx_set(2,1,&q_out, q_scale*mtx_get(2,1,&p));
    mtx_set(3,1,&q_out, q_scale*mtx_get(3,1,&p));
    mtx_set(4,1,&q_out, q_scale);

    q_2_dcm(&q_out, dcm_out);
    return 0;  
}

/*
 * Function: rot_obs_vec
 * ------------------------
 *  Rotate observation vectors by theta (degrees)
 *
 *  result: mtx_matrix of rotated 3 vectors
 */
void rot_obs_vec(struct mtx_matrix* b_k, float theta, int axis,
                    struct mtx_matrix* b_out){
    int i;
    struct mtx_matrix r;
    float theta_rad = M_PI*theta/180;
    mtx_create_eye(3,3,&r); 
    
    /*Create DCM for rotation about x,y,z axis according to input*/
    switch(axis){
        case 1 : 
            mtx_set(2,2,&r,   cos(theta_rad));
            mtx_set(2,3,&r,  -sin(theta_rad));
            mtx_set(3,2,&r,   sin(theta_rad));
            mtx_set(3,3,&r,   cos(theta_rad));
            break;
        case 2 : 
            mtx_set(1,1,&r,   cos(theta_rad));
            mtx_set(1,3,&r,   sin(theta_rad));
            mtx_set(3,1,&r,  -sin(theta_rad));
            mtx_set(3,3,&r,   cos(theta_rad));
            break;
        case 3 : 
            mtx_set(1,1,&r,   cos(theta_rad));
            mtx_set(1,2,&r,  -sin(theta_rad));
            mtx_set(2,1,&r,   sin(theta_rad));
            mtx_set(2,2,&r,   cos(theta_rad));
            break;
    }
    /*Multiply input vectors by R*/
    mtx_mult(&r,b_k,b_out); 
}

/*
 * Function: est_quest
 * ------------------------
 *  A wrapper around the est_quest_rp function to handle the case 
 *  when the rotation from ECI->Body Frame is at 180 deg.
 *
 *  result: mtx_matrix of eci->body frame rotation
 */
int est_quest(struct mtx_matrix* b_k, struct mtx_matrix* eci_k, 
                    struct mtx_matrix* dcm_out){
    int k;
    struct mtx_matrix obs_vec;
    mtx_create_ones(3,2,&obs_vec);

    /*Try rotation around each axis*/
    mtx_copy(b_k,&obs_vec);
    for (k=0;k<=3;k++){ 
        if(est_quest_rp(&obs_vec,eci_k,dcm_out)==0){
            return 0;
        } else {
            rot_obs_vec(b_k,ROT_ANGLE,k+1,&obs_vec);
        }
    }
    return -1;
}


/*
 * Function: est_simp
 * ------------------------
 *  Simplified Estimator
 *  Integrates rate gyros if necessary and calls QUEST estimator
 *
 *  result: mtx_matrix of eci->body frame rotation
 */

int est_simp(struct mtx_matrix* sun_sens_volt, struct mtx_matrix* sun_sens_norm,
                    struct mtx_matrix* mag_sens, 
                    struct mtx_matrix* sun_eci, struct mtx_matrix* mag_eci, 
                    struct est_state* state, struct mtx_matrix* dcm_out){

    struct mtx_matrix sun_vec;
    struct mtx_matrix b_k;
    struct mtx_matrix eci_k;
    struct mtx_matrix prior_dcm;
    struct mtx_matrix post_dcm;
    mtx_create_ones(3,1,&sun_vec);
    mtx_create_ones(3,2,&b_k);
    mtx_create_ones(3,2,&eci_k);
    mtx_create_ones(3,3,&prior_dcm);
    mtx_create_ones(3,3,&post_dcm);

    /*Add magnetic field vector to observation matrix*/
    mtx_set(1,1,&b_k,   mtx_get(1,1,mag_sens));    
    mtx_set(2,1,&b_k,   mtx_get(2,1,mag_sens));    
    mtx_set(3,1,&b_k,   mtx_get(3,1,mag_sens));    
    mtx_set(1,1,&eci_k, mtx_get(1,1,mag_eci));
    mtx_set(2,1,&eci_k, mtx_get(2,1,mag_eci));
    mtx_set(3,1,&eci_k, mtx_get(3,1,mag_eci));

    /*If three sensors are active, add sun sensor to observation matrix
     *otherwise add integrated gyro X vector*/
    if(est_sun_vec_ls(sun_sens_volt, sun_sens_norm, &sun_vec)==0){
        mtx_set(1,2,&b_k,   mtx_get(1,1,&sun_vec));    
        mtx_set(2,2,&b_k,   mtx_get(2,1,&sun_vec));    
        mtx_set(3,2,&b_k,   mtx_get(3,1,&sun_vec));    
        mtx_set(1,2,&eci_k, mtx_get(1,1,sun_eci));
        mtx_set(2,2,&eci_k, mtx_get(2,1,sun_eci));
        mtx_set(3,2,&eci_k, mtx_get(3,1,sun_eci));
    } else {
        q_2_dcm(&state->att_quaternion,&prior_dcm);
        body_rate_dcm_rot(&state->sat_body_rates,&prior_dcm,&post_dcm);
        mtx_set(1,2,&b_k,   mtx_get(1,1,&post_dcm));    
        mtx_set(2,2,&b_k,   mtx_get(1,2,&post_dcm));    
        mtx_set(3,2,&b_k,   mtx_get(1,3,&post_dcm));    
        mtx_set(1,2,&eci_k, 1);
        mtx_set(2,2,&eci_k, 0);
        mtx_set(3,2,&eci_k, 0);
    }
    est_quest(&b_k,&eci_k,dcm_out);
}


/*
 * Function: create_att_state
 * ------------------------
 * Initialize estimator state struct 
 *
 */

int create_att_state(struct est_state* state){
    mtx_create_val(3,3,&state->att_err,0);
    mtx_create_val(3,1,&state->sat_body_rates,0);
    mtx_create_val(3,1,&state->gyro_bias,0);
    mtx_create_val(4,1,&state->att_quaternion,0);
    mtx_set(4,1,&state->att_quaternion,1);
}
