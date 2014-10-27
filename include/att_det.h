#define SS_V_CUTOFF	.01
#define SS_COUNT	13
#define UPDATE_RATE	2.0

struct est_state{
	float att_quaternion[4];
	float att_err[9];
	float sat_body_rates[3];
	float gyro_bias[3];
}; 

void q_2_dcm(struct mtx_matrix* q, struct mtx_matrix* dcm);
void dcm_2_q(struct mtx_matrix* dcm, struct mtx_matrix* q);
void body_rate_dcm_rot(struct mtx_matrix* body_rates, struct mtx_matrix* prior_dcm, struct mtx_matrix* rot_dcm);
int est_sun_vec_ls(struct mtx_matrix* sun_sens_volt, struct mtx_matrix* sun_sens_norm, struct mtx_matrix* sun_vec);
