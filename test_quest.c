#include <stdio.h>
#include "mtx.h"
#include "att_det.h"

int main(void) {
    struct mtx_matrix dcm_est;
    struct mtx_matrix sun_sens_norm;
    struct mtx_matrix sun_sens_volt;
    struct mtx_matrix mag_sens;
    struct mtx_matrix mag_eci;
    struct mtx_matrix sun_eci;
    struct est_state state;

    create_att_state(&state);
    mtx_set(1,1,&state.sat_body_rates,   0.006767668354614);
    mtx_set(2,1,&state.sat_body_rates,  -0.001273360011606);  
    mtx_set(3,1,&state.sat_body_rates,  -0.001961282664780);

    mtx_set(1,1,&state.att_quaternion,  -0.052532019728704);
    mtx_set(2,1,&state.att_quaternion,   0.061096801318250);
    mtx_set(3,1,&state.att_quaternion,   0.708975579587270);
    mtx_set(4,1,&state.att_quaternion,  -0.700614869468809);

    float norm_sun_sens_cart[] = {  0.9418,   0.2879,  -0.1736,
                                    0.9397,        0,   0.3420,
                                    0.9418,  -0.2879,   0.1736,
                                   -0.9366,  -0.3043,  -0.1736,
                                   -0.9397,        0,   0.3420,
                                   -0.9366,   0.3043,  -0.1736,
                                         0,   0.8660,  -0.5000,
                                         0,   0.8660,   0.5000,
                                   -0.5000,   0.8660,        0,
                                         0,  -0.8660,  -0.5000,
                                         0,  -0.8660,   0.5000,
                                   -0.5000,  -0.8660,        0,
                                         0,        0,   1.0000};

    float sun_sens_read[] =  { 0,
                               8.000000000000,
                               0005.000000000000,
                               872.999987870389,
                               0050.000000000000,
                               0027.000000000000,
                               0001.000000000000,
                               0017.999987870389,
                               0015.000000000000,
                               2042.000000000000,
                               0760.000012129611,
                               1480.000000000000,
                               48.000000000000};


    float mag_read[] = {      503.00012633273,
                              -308.00000467899,
                              19352.00001871596};
    float mag_eci_vals[] = { .0000007722672871089,
                             .0000053648825749822,
                             .0000284237836229281};
    float sun_eci_vals[] = {  0,
                              0,
                              0};


    mtx_create(3,1,mag_read,&mag_sens);
    mtx_create(3,1,mag_eci_vals,&mag_eci);
    mtx_create(3,1,sun_eci_vals,&sun_eci);
    mtx_create(13,3,norm_sun_sens_cart,&sun_sens_norm);
    mtx_create(13,1,sun_sens_read,&sun_sens_volt);
    mtx_create_ones(3,3,&dcm_est);
    mtx_create_val(13,1,&sun_sens_volt,0);

    printf("State: Attidue Error\n\r");
    mtx_print(&state.att_err);
    printf("State: Gyro Bias\n\r");
    mtx_print(&state.gyro_bias);
    printf("State: Body Rates\n\r");
    mtx_print(&state.sat_body_rates);
    printf("State: Quaternion\n\r");
    mtx_print(&state.att_quaternion);

    est_simp(&sun_sens_volt,&sun_sens_norm,&mag_sens,
                &sun_eci,&mag_eci,&state,&dcm_est);

    printf("Outputs:\n\r");
    mtx_print(&dcm_est);
    
}
