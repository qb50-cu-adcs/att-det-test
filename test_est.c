#include "att_det.h"
#include <stdio.h>

int main(void) {
	//Input
	double gyro_body_rates[] = {.01,.02,.03}; //DPS
	double mag_vec[] = {.1, 1, 10}; //T
	double sun_sens_volt[] = {1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; //V
	double mag_eci_vec[] = {1, 2, 3}; //Normalized Cartesian
	double sun_eci_vec[] = {4, 5, 6}; //Normalized Cartesian
	double sun_sens_norm[] = {1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1,
				  1,0,0, 0,1,0, 0,0,1};

	//Output
	double att_quaternion[4]; //Quaternion
	double sat_body_rates[3]; //DPS 	

	est_svd(gyro_body_rates, mag_vec, sun_sens_volt, sun_sens_norm, sun_eci_vec, mag_eci_vec, att_quaternion, sat_body_rates);

	printf("Attitude: [%.2f, %.2f, %.2f, %.2f]   Body Rates: [%.2f, %.2f, %.2f]\n", att_quaternion[0],att_quaternion[1],att_quaternion[2],att_quaternion[3],sat_body_rates[0],sat_body_rates[1],sat_body_rates[2]);
	return 0;
}
