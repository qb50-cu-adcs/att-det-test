#define SS_V_CUTOFF	.01
#define SS_COUNT	15
#define UPDATE_RATE	2.0

struct est_state{
	float att_quaternion[4];
	float att_err[9];
	float sat_body_rates[3];
	float gyro_bias[3];
}; 

void q_2_dcm(struct mtx_matrix* q, struct mtx_matrix* dcm);
void dcm_2_q(struct mtx_matrix* dcm, struct mtx_matrix* q);
