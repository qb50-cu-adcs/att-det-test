#include "att_det.h"
#include <stdio.h>

int main(void) {
	//Input
	double gyro_body_rates[] = {.01,.02,.03}; //DPS
	double mag_vec[] = {.1, 1, 10}; //T
	double sun_sens_volt[] = {0.02, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}; //V
	double mag_eci_vec[] = {1, 2, 3}; //Normalized Cartesian
	double sun_eci_vec[] = {4, 5, 6}; //Normalized Cartesian
	double sun_sens_norm[] = {1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1};
	//Output
	struct est_state state;

	est_svd(gyro_body_rates, mag_vec, sun_sens_volt, sun_sens_norm, sun_eci_vec, mag_eci_vec, &state);

	printf("Attitude: [%.2f, %.2f, %.2f, %.2f]   Body Rates: [%.2f, %.2f, %.2f]\n", state.att_quaternion[0], state.att_quaternion[1], state.att_quaternion[2], state.att_quaternion[3], state.sat_body_rates[0], state.sat_body_rates[1], state.sat_body_rates[2]);
	return 0;
}
