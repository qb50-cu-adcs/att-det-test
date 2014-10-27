#include <stdio.h>
#include "mtx.h"
#include "att_det.h"

int main(void) {
    struct mtx_matrix rates_vec;
    struct mtx_matrix prior_dcm;
    struct mtx_matrix dcm_out;
    struct mtx_matrix sun_sens_norm;
    struct mtx_matrix sun_sens_volt;
    struct mtx_matrix sun_vec;

    float rates[] = { .0006,
                      .0010,
                     -.0004};

    float prior[] = { -.0367, -.9987, -.0367,
                      -.3616, -.0210,  .9321,
                      -.9316,  .0475, -.3604};

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

    float sun_sens_read[] = {   0.0333,
                                     0,
                                     0,
                                2.1565,
                                2.0621,
                                0.5380,
                                0.3883,
                                0.1457,
                                0.4012,
                                1.5100,
                                2.7357,
                                2.8886,
                                1.4108};

    mtx_create(3,1,rates,&rates_vec);
    mtx_create(3,3,prior,&prior_dcm);
    mtx_create_ones(3,3,&dcm_out);
    mtx_create_ones(3,1,&sun_vec);
    mtx_create(13,3,norm_sun_sens_cart,&sun_sens_norm);
    mtx_create(13,1,sun_sens_read,&sun_sens_volt);

    body_rate_dcm_rot(&rates_vec, &prior_dcm, &dcm_out);
    est_sun_vec_ls(&sun_sens_volt, &sun_sens_norm, &sun_vec);

    printf("Inputs:\n\r");
    mtx_print(&rates_vec);
    mtx_print(&prior_dcm);
    printf("Outputs:\n\r");
    mtx_print(&dcm_out);
    mtx_print(&sun_vec);

}
