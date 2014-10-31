#define SS_V_CUTOFF	204
#define SS_COUNT	13
#define UPDATE_RATE	2.0
#define THETA_RP .25
#define MAX_ROTS 3
#define ROT_ANGLE .5
#define MIN_MINV_DET .00001

struct est_state{
	struct mtx_matrix att_quaternion; //4x1
	struct mtx_matrix att_err; //3x3
	struct mtx_matrix sat_body_rates; //3x1
	struct mtx_matrix gyro_bias; //3x1
}; 

void q_2_dcm(struct mtx_matrix* q, struct mtx_matrix* dcm);
void dcm_2_q(struct mtx_matrix* dcm, struct mtx_matrix* q);
void body_rate_dcm_rot(struct mtx_matrix* body_rates, struct mtx_matrix* prior_dcm, struct mtx_matrix* rot_dcm);
int est_sun_vec_ls(struct mtx_matrix* sun_sens_volt, struct mtx_matrix* sun_sens_norm, struct mtx_matrix* sun_vec);
int est_quest_rp(struct mtx_matrix* b_k, struct mtx_matrix* eci_k, struct mtx_matrix* dcm_out);
int est_quest(struct mtx_matrix* b_k, struct mtx_matrix* eci_k, struct mtx_matrix* dcm_out);
//int est_simp(struct mtx_matrix* sun_sens_volt, struct mtx_matrix* sun_sens_norm, struct mtx_matrix* mag_sens, struct mtx_matrix* gyro_rates, struct mtx_matrix* sun_eci, struct mtx_matrix* mag_eci, struct att_state state, struct mtx_matrix* dcm_out);
int create_att_state(struct est_state* state);
