#include "att_det.h"
#include <stdio.h>

int main(void) {
	//Input
	double gyro_body_rates[] = {.0029,-.0025,.0006}; //DPS
	double mag_vec[] = {-31901, 9583, 1017}; //nT
	double sun_sens_volt[] = {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}; //V
	double mag_eci_vec[] = {-0.2814, -0.0320, -0.9591}; //Normalized Cartesian
	double sun_eci_vec[] = {4, 5, 6}; //Normalized Cartesian
	double sun_sens_norm[] = {1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1};
	//Output
	struct est_state state;
	state.att_quaternion[0] = .4868462445;
	state.att_quaternion[1] = -.5004750141;
	state.att_quaternion[2] = .5088538058;
	state.att_quaternion[3] = -.5036062082;

	est_svd(gyro_body_rates, mag_vec, sun_sens_volt, sun_sens_norm, sun_eci_vec, mag_eci_vec, &state);

	printf("Attitude: [%.4f, %.4f, %.4f, %.4f]   Body Rates: [%.4f, %.4f, %.4f]\n", state.att_quaternion[0], state.att_quaternion[1], state.att_quaternion[2], state.att_quaternion[3], state.sat_body_rates[0], state.sat_body_rates[1], state.sat_body_rates[2]);
	return 0;
}
