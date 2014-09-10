#define SS_V_CUTOFF	.01
#define SS_COUNT	15
struct est_state{
	double att_quaternion[4];
	double att_err[9];
	double sat_body_rates[3];
	double gyro_bias[3];
}; 
void est_svd(double* gyro_body_rates, double* mag_vec, double* sun_sens_volt, double* sun_sens_norm, double* sun_eci_vec, double* mag_eci_vec, struct est_state* state);
