#include <iostream>
#include "math.h"
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace std;

/* cccc  SK-I relic */
/* Subroutine */ int pdf_rel_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.88993e-5f;
/* Computing 2nd power */
    d_1 = (*e + 76.23f) / 31.29f;
    fit = exp(d_1 * d_1 * -.5f) * 151600;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_rel_sk1_med_ */

/* Subroutine */ int pdf_rel_sk1_med_old_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.91733e-5f;
/* Computing 2nd power */
    d_1 = (*e + 71.56f) / 30.68f;
    fit = exp(d_1 * d_1 * -.5f) * 97450;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_rel_sk1_med_old_ */

/* Subroutine */ int pdf_rel_sk1_low_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
/*      c = 0.0000691733 */
    c_ = 5.88993e-5f;
    fit = exp(6.332f - *e * .2283f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_rel_sk1_low_ */

/* Subroutine */ int pdf_rel_sk1_high_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.88993e-5f;
/*      c = 0.0000691733 */
    fit = exp(6.198f - *e * .1972f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_rel_sk1_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccc */
/* cccc SK-I nu-mu (from decay electron data) */
/* Subroutine */ int pdf_de_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2, d_3;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    c_ = 9.24637e-7f;
/* Computing 2nd power */
    d_1 = (*e - 45.954f) / 6.993f;
/* Computing 2nd power */
    d_2 = (*e - 54.734f) / 5.4754f;
/* Computing 2nd power */
    d_3 = (*e - 30.715f) / 9.6517f;
    fit = exp(d_1 * d_1 * -.5f) * 31049 - exp(d_2 * d_2 * -.5f) * 4196.6f 
	    + exp(d_3 * d_3 * -.5f) * 26559;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 18. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 16. && *e <= 18.) {
	*value = spalllow * c_ * fit;
    }
    return 0;
} /* pdf_de_sk1_med_ */

/* Subroutine */ int pdf_de_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2, d_3;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    c_ = 9.24637e-7f;
/* Computing 2nd power */
    d_1 = (*e - 19.75f) / 5.31f;
/* Computing 2nd power */
    d_2 = (*e - 14.8f) / 14.4f;
/* Computing 2nd power */
    d_3 = (*e + 11.62f) / 21.25f;
    fit = exp(d_1 * d_1 * -.5f) * 49.345f + exp(d_2 * d_2 * -.5f) * 
	    54.285f - exp(d_3 * d_3 * -.5f) * 29.74f;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 18. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 16. && *e <= 18.) {
	*value = spalllow * c_ * fit;
    }
    return 0;
} /* pdf_de_sk1_low_ */

/* Subroutine */ int pdf_de_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    c_ = 9.24637e-7f;
/* Computing 2nd power */
    d_1 = (*e - 30) / 12.2f;
    fit = exp(d_1 * d_1 * -.5f) * 160.2f;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 18. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 16. && *e <= 18.) {
	*value = spalllow * c_ * fit;
    }
    return 0;
} /* pdf_de_sk1_high_ */

/* ccccccccccccccccccccccccccccccccccccccc */
/* cccc SK-I nue (from atmpd MC) */
/* Subroutine */ int pdf_nue_sk1_med_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.22464e-4f;
    fit = *e * .034f - 6.06f + *e * .0578f * *e - *e * 3.51e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk1_med_ */

/* Subroutine */ int pdf_nue_sk1_low_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.22464e-4f;
    fit = .13513513513513514f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk1_low_ */

/* Subroutine */ int pdf_nue_sk1_high_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.22464e-4f;
    fit = 13.776f - *e * .6f + *e * .01086f * *e - *e * 6.151e-5f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk1_high_ */

/* Subroutine */ int pdf_nue_sk1_med_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.30624e-5f;
    fit = *e * .5481f - 8.602f + *e * .1177f * *e - *e * 7.135e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk1_med_old_ */

/* Subroutine */ int pdf_nue_sk1_low_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.30624e-5f;
    fit = .27f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk1_low_old_ */

/* Subroutine */ int pdf_nue_sk1_high_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.30624e-5f;
    fit = 27.351f - *e * 1.1214f + *e * .02026f * *e - *e * 1.15e-4f * *e * *
	    e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk1_high_old_ */

/* ccccccccccccccccccccccccccccccccccccccc */
/* cccccccccc NC elastic ccccccccccccccc */
/* Subroutine */ int pdf_ncel_sk1_med_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.49925e-4f;
    fit = exp(11.348f - *e * .3837f) + exp(1.78f - *e * .02323f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk1_med_ */

/* Subroutine */ int pdf_ncel_sk1_low_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.49925e-4f;
    fit = 11.f - *e * .44f + *e * .01f * *e - *e * 6.83e-5f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk1_low_ */

/* Subroutine */ int pdf_ncel_sk1_high_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.49925e-4f;
    fit = exp(10.37f - *e * .21f) + exp(5.676f - *e * .06384f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk1_high_ */

/* Subroutine */ int pdf_ncel_sk1_med_old_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.16289e-5f;
    fit = exp(12.178f - *e * .39023f) + exp(3.286f - *e * .033467f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk1_med_old_ */

/* Subroutine */ int pdf_ncel_sk1_low_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.16289e-5f;
    fit = 29.316f - *e * 1.151f + *e * .02525f * *e - *e * 1.7717e-4f * *e * *
	    e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk1_low_old_ */

/* Subroutine */ int pdf_ncel_sk1_high_old_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.16289e-5f;
    fit = exp(11.05f - *e * .20783f) + exp(6.226f - *e * .058322f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk1_high_old_ */

/* cccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccc mu/pi cccccccccccccccccccccccc */
/* Subroutine */ int pdf_mupi_sk1_med_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.039e-4f;
    fit = 15.04f - *e * .5956f + *e * .01675f * *e - *e * 1.18e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk1_med_ */

/* Subroutine */ int pdf_mupi_sk1_low_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.039e-4f;
    fit = *e * .0447f + 50.14f - *e * 1.374e-4f * *e - *e * 2.676e-5f * *e * *
	    e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk1_low_ */

/* Subroutine */ int pdf_mupi_sk1_high_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.039e-4f;
    fit = 74.48f - *e * 2.788f + *e * .03636f * *e - *e * 1.605e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk1_high_ */

/* Subroutine */ int pdf_mupi_sk1_med_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.00174e-5f;
    fit = 64.8f - *e * 3.9f + *e * .1263f * *e - *e * .0015358f * *e * *e + *
	    e * 6.3776e-6f * *e * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk1_med_old_ */

/* Subroutine */ int pdf_mupi_sk1_low_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.00174e-5f;
    fit = *e * 4.5457f + 118 - *e * .087851f * *e + *e * 4.3357e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk1_low_old_ */

/* Subroutine */ int pdf_mupi_sk1_high_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.00174e-5f;
    fit = 223.71f - *e * 8.195f + *e * .1053f * *e - *e * 4.6262e-4f * *e * *
	    e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk1_high_old_ */

/* cccc  SK-II relic */
/* Subroutine */ int pdf_rel_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
/*      c = 0.0000927416 */
    c_ = 7.83959e-5f;
/* Computing 2nd power */
    d_1 = (*e + 26.47f) / 14.82f;
    fit = exp(d_1 * d_1 * -.5f) * 4453;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_rel_sk2_low_ */

/* Subroutine */ int pdf_rel_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 7.83959e-5f;
/* Computing 2nd power */
    d_1 = (*e + 30.78f) / 24.74f;
    fit = exp(d_1 * d_1 * -.5f) * 8993;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_rel_sk2_med_ */

/* Subroutine */ int pdf_rel_sk2_med_old_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.27416e-5f;
/* Computing 2nd power */
    d_1 = (*e - 21.84f) / 5.194f;
/* Computing 2nd power */
    d_2 = (*e + 60.82f) / 30.14f;
    fit = exp(d_1 * d_1 * -.5f) * 122.4f + exp(d_2 * d_2 * -.5f) * 30340;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_rel_sk2_med_old_ */

/* Subroutine */ int pdf_rel_sk2_high_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
/*      c = 0.0000927416 */
    c_ = 7.83959e-5f;
    fit = exp(6.286f - *e * .1704f);
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_rel_sk2_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccc Decaye electron PDFs ccccccccccccc */
/* cccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf_de_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 5.39269e-5f;
/* Computing 2nd power */
    d_1 = (*e - 13.68f) / 15.18f;
    fit = exp(d_1 * d_1 * -.5f) * 17.49f;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_de_sk2_low_ */

/* Subroutine */ int pdf_de_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 5.39269e-5f;
/* Computing 2nd power */
    d_1 = (*e - 46.11f) / 7.142f;
/* Computing 2nd power */
    d_2 = (*e - 32.2f) / 10.64f;
    fit = exp(d_1 * d_1 * -.5f) * 361.f + exp(d_2 * d_2 * -.5f) * 498.f;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_de_sk2_med_ */

/* Subroutine */ int pdf_de_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 5.39269e-5f;
/* Computing 2nd power */
    d_1 = (*e - 32.15f) / 13.31f;
    fit = exp(d_1 * d_1 * -.5f) * 4.767f;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_de_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccc */
/* cccccccccccccccc NU-E CC cccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf_nue_sk2_low_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.48658e-4f;
    fit = .17567567567567569f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_nue_sk2_low_ */

/* Subroutine */ int pdf_nue_sk2_med_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.48658e-4f;
    fit = *e * .18f - 3.4226f + *e * .02076f * *e - *e * 1.039e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_nue_sk2_med_ */

/* Subroutine */ int pdf_nue_sk2_high_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.48658e-4f;
    fit = *e * .02841f + .782f + *e * 4.257e-4f * *e - *e * 2.192e-7f * *e * *
	    e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_nue_sk2_high_ */

/* Subroutine */ int pdf_nue_sk2_low_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.99758e-5f;
    fit = .1891891891891892f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_nue_sk2_low_old_ */

/* Subroutine */ int pdf_nue_sk2_med_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.99758e-5f;
    fit = 4.514f - *e * .4219f + *e * .06915f * *e - *e * 3.418e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_nue_sk2_med_old_ */

/* Subroutine */ int pdf_nue_sk2_high_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.99758e-5f;
    fit = *e * .0637f + .069f - *e * 1.364e-4f * *e + *e * 2.68e-6f * *e * *e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_nue_sk2_high_old_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccc NC ELASTIC cccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf_ncel_sk2_low_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 4.47633e-4f;
    fit = exp(12.02f - *e * .586f) + exp(1.048f - *e * .0077f);
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_ncel_sk2_low_ */

/* Subroutine */ int pdf_ncel_sk2_med_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 4.47633e-4f;
    fit = exp(10.377f - *e * .3579f) + exp(.4096f - *e * .004236f);
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_ncel_sk2_med_ */

/* Subroutine */ int pdf_ncel_sk2_high_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 4.47633e-4f;
    fit = exp(8.816f - *e * .1724f) + exp(3.391f - *e * .04f);
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_ncel_sk2_high_ */

/* Subroutine */ int pdf_ncel_sk2_low_old_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.63703e-4f;
    fit = exp(11.31f - *e * .5312f) + exp(1.01f - *e * .006344f);
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_ncel_sk2_low_old_ */

/* Subroutine */ int pdf_ncel_sk2_med_old_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.63703e-4f;
    fit = exp(11.81f - *e * .3875f) + exp(3.837f - *e * .04744f);
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_ncel_sk2_med_old_ */

/* Subroutine */ int pdf_ncel_sk2_high_old_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.63703e-4f;
    fit = exp(8.819f - *e * .1727f) + exp(3.668f - *e * .0429f);
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_ncel_sk2_high_old_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccc */
/* cccccccccccccccc MU / PI ccccccccccccccccccccccc */
/* cccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf_mupi_sk2_low_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.77227e-4f;
    fit = *e * .5983f + 8.05f - *e * .004924f * *e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_mupi_sk2_low_ */

/* Subroutine */ int pdf_mupi_sk2_med_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.77227e-4f;
    fit = *e * .0658f + 7.781f - *e * 8.636e-4f * *e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_mupi_sk2_med_ */

/* Subroutine */ int pdf_mupi_sk2_high_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.77227e-4f;
    fit = 22.481f - *e * .703f + *e * .008161f * *e - *e * 2.9887e-5f * *e * *
	    e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_mupi_sk2_high_ */

/* Subroutine */ int pdf_mupi_sk2_low_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.94671e-4f;
    fit = *e * .8674f + 19.47f - *e * .00749f * *e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_mupi_sk2_low_old_ */

/* Subroutine */ int pdf_mupi_sk2_med_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.94671e-4f;
    fit = 31.52f - *e * .9216f + *e * .02158f * *e - *e * 1.317e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_mupi_sk2_med_old_ */

/* Subroutine */ int pdf_mupi_sk2_high_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.94671e-4f;
    fit = 30.15f - *e * .7076f + *e * .004573f * *e;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_mupi_sk2_high_old_ */

/* cccc  SK-III relic */
/* Subroutine */ int pdf_rel_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.92062e-5f;
/* Computing 2nd power */
    d_1 = (*e + 60.28f) / 29.046f;
    fit = exp(d_1 * d_1 * -.5f) * 51280;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_rel_sk3_med_ */

/* Subroutine */ int pdf_rel_sk3_low_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.92062e-5f;
    fit = exp(5.6884f - *e * .20983f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_rel_sk3_low_ */

/* Subroutine */ int pdf_rel_sk3_high_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.92062e-5f;
    fit = exp(5.0965f - *e * .13258f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_rel_sk3_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccc */
/* cccc SK-III nu-mu (from decay electron data) */
/* Subroutine */ int pdf_de_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 4.39444e-6f;
/* Computing 2nd power */
    d_1 = (*e - 46.81f) / 5.784f;
/* Computing 2nd power */
    d_2 = (*e - 34.4f) / 9.28f;
    fit = exp(d_1 * d_1 * -.5f) * 5095.f + exp(d_2 * d_2 * -.5f) * 5790.f 
	    + 2007.f - *e * 40.47f - *e * .144f * *e + *e * .004445f * *e * *
	    e;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_de_sk3_med_ */

/* Subroutine */ int pdf_de_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    c_ = 4.39444e-6f;
/* Computing 2nd power */
    d_1 = (*e - 6.f) / 17.59f;
    fit = exp(d_1 * d_1 * -.5f) * 28.55f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_de_sk3_low_ */

/* Subroutine */ int pdf_de_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    c_ = 4.39444e-6f;
/* Computing 2nd power */
    d_1 = (*e - 32.02f) / 12.04f;
    fit = exp(d_1 * d_1 * -.5f) * 35.54f;
    *value = c_ * fit;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    *value = c_ * fit;
    if (*e > 70.) {
	*value = 0.;
    }
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_de_sk3_high_ */

/* ccccccccccccccccccccccccccccccccccccccc */
/* cccc SK-III nue (from atmpd MC) */
/* Subroutine */ int pdf_nue_sk3_med_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.13332e-4f;
    fit = 5.8984f - *e * .55787f + *e * .073688f * *e - *e * 4.588e-4f * *e * 
	    *e;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk3_med_ */

/* Subroutine */ int pdf_nue_sk3_low_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.13332e-4f;
    fit = .067567567567567571f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk3_low_ */

/* Subroutine */ int pdf_nue_sk3_high_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.13332e-4f;
    fit = 8.8303f - *e * .41617f + *e * .009215f * *e - *e * 5.697e-5f * *e * 
	    *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk3_high_ */

/* Subroutine */ int pdf_nue_sk3_med_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.13501e-4f;
    fit = 6.6336f - *e * .6248f + *e * .07505f * *e - *e * 4.6687e-4f * *e * *
	    e;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk3_med_old_ */

/* Subroutine */ int pdf_nue_sk3_low_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.13501e-4f;
    fit = .067567567567567571f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk3_low_old_ */

/* Subroutine */ int pdf_nue_sk3_high_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.13501e-4f;
    fit = 8.1214f - *e * .37767f + *e * .0085575f * *e - *e * 5.343e-5f * *e *
	     *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_nue_sk3_high_old_ */

/* ccccccccccccccccccccccccccccccccccccccc */
/* cccccccccc NC elastic ccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf_ncel_sk3_med_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.55626e-4f;
    fit = exp(9.514f - *e * .29444f) + exp(1.161f - *e * .012166f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk3_med_ */

/* Subroutine */ int pdf_ncel_sk3_low_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.55626e-4f;
    fit = 11.337f - *e * .385f + *e * .007695f * *e - *e * 5.096e-5f * *e * *
	    e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk3_low_ */

/* Subroutine */ int pdf_ncel_sk3_high_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.55626e-4f;
    fit = exp(9.835f - *e * .1833f) + exp(5.06f - *e * .05371f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk3_high_ */

/* Subroutine */ int pdf_ncel_sk3_med_old_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.51298e-4f;
    fit = exp(9.6116f - *e * .3086f) + exp(2.22f - *e * .0265f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk3_med_old_ */

/* Subroutine */ int pdf_ncel_sk3_low_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.51298e-4f;
    fit = 9.264f - *e * .2094f + *e * .0044f * *e - *e * 3.202e-5f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk3_low_old_ */

/* Subroutine */ int pdf_ncel_sk3_high_old_(double *e, double *value)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.51298e-4f;
    fit = exp(9.9f - *e * .1874f) + exp(5.367f - *e * .05713f);
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ncel_sk3_high_old_ */

/* cccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccc mu/pi cccccccccccccccccccccccc */
/* cccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf_mupi_sk3_med_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.0449e-4f;
    fit = 12.614f - *e * .34755f + *e * .009822f * *e - *e * 6.8091e-5f * *e *
	     *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk3_med_ */

/* Subroutine */ int pdf_mupi_sk3_low_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.04489e-4f;
    fit = *e * .70352f + 42.f - *e * .013613f * *e + *e * 5.38e-5f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk3_low_ */

/* Subroutine */ int pdf_mupi_sk3_high_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.04489e-4f;
    fit = 59.883f - *e * 1.9f + *e * .02053f * *e - *e * 7.374e-5f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk3_high_ */

/* Subroutine */ int pdf_mupi_sk3_med_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.21094e-4f;
    fit = 24.8f - *e * .9827f + *e * .02329f * *e - *e * 1.571e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk3_med_old_ */

/* Subroutine */ int pdf_mupi_sk3_low_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.21094e-4f;
    fit = *e * 2.894f + 52.38f - *e * .05557f * *e + *e * 2.748e-4f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk3_low_old_ */

/* Subroutine */ int pdf_mupi_sk3_high_old_(double *e, double *value)
{
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.21094e-4f;
    fit = 78.09f - *e * 2.014f + *e * .01471f * *e - *e * 1.846e-5f * *e * *e;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mupi_sk3_high_old_ */

/* cccc  SK-I functions */
/* Subroutine */ int pdf2_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.87816e-4f;
/* Computing 2nd power */
    d_1 = (*e - .165f) / 5.61f;
    fit = exp(d_1 * d_1 * -.5f) * 3894;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_sk1_low_ */

/* Subroutine */ int pdf2_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.87816e-4f;
/* Computing 2nd power */
    d_1 = (*e + 9.83f) / 9.19f;
    fit = exp(d_1 * d_1 * -.5f) * 51340;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_sk1_med_ */

/* Subroutine */ int pdf2_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.87816e-4f;
/* Computing 2nd power */
    d_1 = (*e + 4.95f) / 7.56f;
    fit = exp(d_1 * d_1 * -.5f) * 7356;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_sk1_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf2_half_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.85143e-5f;
/* Computing 2nd power */
    d_1 = (*e + 2.37f) / 6.23f;
    fit = exp(d_1 * d_1 * -.5f) * 8655;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_half_sk1_low_ */

/* Subroutine */ int pdf2_half_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.85143e-5f;
/* Computing 2nd power */
    d_1 = (*e + 17.02f) / 11.68f;
    fit = exp(d_1 * d_1 * -.5f) * 6.42e5f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_half_sk1_med_ */

/* Subroutine */ int pdf2_half_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.85143e-5f;
/* Computing 2nd power */
    d_1 = (*e + 12.6f) / 9.345f;
    fit = exp(d_1 * d_1 * -.5f) * 20350;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_half_sk1_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf3_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.06284e-4f;
/* Computing 2nd power */
    d_1 = (*e - 7.89f) / 10.08f;
    fit = exp(d_1 * d_1 * -.5f) * 4.116f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_sk1_low_ */

/* Subroutine */ int pdf3_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.06284e-4f;
/* Computing 2nd power */
    d_1 = (*e + 22.75f) / 14.05f;
    fit = exp(d_1 * d_1 * -.5f) * 57850;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_sk1_med_ */

/* Subroutine */ int pdf3_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 2.06284e-4f;
/* Computing 2nd power */
    d_1 = (*e + 19.92f) / 11.75f;
    fit = exp(d_1 * d_1 * -.5f) * 2212;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_sk1_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf3_half_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.64213e-5f;
/* Computing 2nd power */
    d_1 = (*e + 5.6081f) / 7.9722f;
    fit = exp(d_1 * d_1 * -.5f) * 4026.9f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_half_sk1_low_ */

/* Subroutine */ int pdf3_half_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.64213e-5f;
/* Computing 2nd power */
    d_1 = (*e + 26.336f) / 16.133f;
    fit = exp(d_1 * d_1 * -.5f) * 414160;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_half_sk1_med_ */

/* Subroutine */ int pdf3_half_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.64213e-5f;
/* Computing 2nd power */
    d_1 = (*e + 14.838f) / 11.792f;
    fit = exp(d_1 * d_1 * -.5f) * 5939;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_half_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf4_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.38783e-4f;
/* Computing 2nd power */
    d_1 = (*e + 3.991f) / 6.434f;
    fit = exp(d_1 * d_1 * -.5f) * 1809;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_sk1_low_ */

/* Subroutine */ int pdf4_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.38783e-4f;
/* Computing 2nd power */
    d_1 = (*e + 20.285f) / 16.345f;
    fit = exp(d_1 * d_1 * -.5f) * 15774;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_sk1_med_ */

/* Subroutine */ int pdf4_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.38783e-4f;
/* Computing 2nd power */
    d_1 = (*e + 48.52f) / 18.71f;
    fit = exp(d_1 * d_1 * -.5f) * 5190;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_sk1_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf4_half_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.19564e-5f;
/* Computing 2nd power */
    d_1 = (*e + 18.016f) / 10.992f;
    fit = exp(d_1 * d_1 * -.5f) * 10520;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_half_sk1_low_ */

/* Subroutine */ int pdf4_half_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.19564e-5f;
/* Computing 2nd power */
    d_1 = (*e + 28.981f) / 19.452f;
    fit = exp(d_1 * d_1 * -.5f) * 190550;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_half_sk1_med_ */

/* Subroutine */ int pdf4_half_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.19564e-5f;
/* Computing 2nd power */
    d_1 = (*e + 35.838f) / 16.777f;
    fit = exp(d_1 * d_1 * -.5f) * 23945;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_half_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf5_half_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 9.65165e-6f;
/* Computing 2nd power */
    d_1 = (*e + 9.1907f) / 9.954f;
/* Computing 2nd power */
    d_2 = (*e - 8.4581f) / -1.9055f;
    fit = exp(d_1 * d_1 * -.5f) * 2331.9f - exp(d_2 * d_2 * -.5f) * 
	    116430;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_half_sk1_low_ */

/* Subroutine */ int pdf5_half_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 9.65165e-6f;
/* Computing 2nd power */
    d_1 = (*e + 22.705f) / 20.976f;
    fit = exp(d_1 * d_1 * -.5f) * 68113;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_half_sk1_med_ */

/* Subroutine */ int pdf5_half_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 9.65165e-6f;
/* Computing 2nd power */
    d_1 = (*e + 48.088f) / 21.192f;
    fit = exp(d_1 * d_1 * -.5f) * 17058;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_half_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf5_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.06955e-4f;
/* Computing 2nd power */
    d_1 = (*e + 54.32f) / 23.06f;
    fit = exp(d_1 * d_1 * -.5f) * 631.6f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_sk1_low_ */

/* Subroutine */ int pdf5_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.06955e-4f;
/* Computing 2nd power */
    d_1 = (*e + 37.17f) / 22.28f;
    fit = exp(d_1 * d_1 * -.5f) * 22580;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_sk1_med_ */

/* Subroutine */ int pdf5_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.06955e-4f;
/* Computing 2nd power */
    d_1 = (*e + 68.46f) / 24.19f;
    fit = exp(d_1 * d_1 * -.5f) * 275;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf6_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.77094e-6f;
/* Computing 2nd power */
    d_1 = (*e + 28.46f) / 13.29f;
    fit = exp(d_1 * d_1 * -.5f) * 23070;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_sk1_low_ */

/* Subroutine */ int pdf6_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.77094e-6f;
/* Computing 2nd power */
    d_1 = (*e + 9.185f) / 19.1f;
/* Computing 2nd power */
    d_2 = (*e + 148.4f) / 42.25f;
    fit = exp(d_1 * d_1 * -.5f) * 22172 + exp(d_2 * d_2 * -.5f) * 4.8e6f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_sk1_med_ */

/* Subroutine */ int pdf6_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.77094e-6f;
/* Computing 2nd power */
    d_1 = (*e + 14.44f) / 16.22f;
    fit = exp(d_1 * d_1 * -.5f) * 2572;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf6_half_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.36936e-6f;
/* Computing 2nd power */
    d_1 = (*e + 32.165f) / 14.473f;
    fit = exp(d_1 * d_1 * -.5f) * 20167;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_half_sk1_low_ */

/* Subroutine */ int pdf6_half_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.36936e-6f;
/* Computing 2nd power */
    d_1 = (*e - 2.3985f) / 17.732f;
/* Computing 2nd power */
    d_2 = (*e + 215.94f) / -48.467f;
    fit = exp(d_1 * d_1 * -.5f) * 8423.2f + exp(d_2 * d_2 * -.5f) * 
	    469640000;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_half_sk1_med_ */

/* Subroutine */ int pdf6_half_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.36936e-6f;
/* Computing 2nd power */
    d_1 = (*e + 32.81f) / 19.977f;
    fit = exp(d_1 * d_1 * -.5f) * 3281.9f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_half_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf7_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.86691e-5f;
/* Computing 2nd power */
    d_1 = (*e + 50.89f) / 32.51f;
    fit = exp(d_1 * d_1 * -.5f) * 17.26f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_sk1_low_ */

/* Subroutine */ int pdf7_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.86691e-5f;
/* Computing 2nd power */
    d_1 = (*e + 52.34f) / 29.32f;
/* Computing 2nd power */
    d_2 = (*e + 25.57f) / 13.22f;
    fit = exp(d_1 * d_1 * -.5f) * 21370 - exp(d_2 * d_2 * -.5f) * 60804;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_sk1_med_ */

/* Subroutine */ int pdf7_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.86691e-5f;
/* Computing 2nd power */
    d_1 = (*e + 95.28f) / 34.97f;
    fit = exp(d_1 * d_1 * -.5f) * 2079;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf7_half_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.49321e-6f;
/* Computing 2nd power */
    d_1 = (*e + 26.763f) / 15.006f;
    fit = exp(d_1 * d_1 * -.5f) * 3804;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_half_sk1_low_ */

/* Subroutine */ int pdf7_half_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.49321e-6f;
/* Computing 2nd power */
    d_1 = (*e + 5.624f) / 22.979f;
/* Computing 2nd power */
    d_2 = (*e - 19.231f) / 7.7718f;
    fit = exp(d_1 * d_1 * -.5f) * 12029 + exp(d_2 * d_2 * -.5f) * 1897;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_half_sk1_med_ */

/* Subroutine */ int pdf7_half_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.49321e-6f;
/* Computing 2nd power */
    d_1 = (*e + 17.556f) / 18.926f;
    fit = exp(d_1 * d_1 * -.5f) * 638.39f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_half_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf8_sk1_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.17383e-5f;
/* Computing 2nd power */
    d_1 = (*e + 30.18f) / 14.08f;
    fit = exp(d_1 * d_1 * -.5f) * 1707;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf8_sk1_low_ */

/* Subroutine */ int pdf8_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.17383e-5f;
/* Computing 2nd power */
    d_1 = (*e + .905f) / 21.73f;
/* Computing 2nd power */
    d_2 = (*e - 22.97f) / 3.889f;
    fit = exp(d_1 * d_1 * -.5f) * 1204 + exp(d_2 * d_2 * -.5f) * 74.63f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf8_sk1_med_ */

/* Subroutine */ int pdf8_sk1_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.17383e-5f;
/* Computing 2nd power */
    d_1 = (*e - 15.34f) / 10.93f;
    fit = exp(d_1 * d_1 * -.5f) * 7.806f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf8_sk1_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cccc SK-II pdfs */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf2_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.43541e-4f;
/* Computing 2nd power */
    d_1 = (*e - .253f) / 6.52f;
    fit = exp(d_1 * d_1 * -.5f) * 5179;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf2_sk2_low_ */

/* Subroutine */ int pdf2_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.43541e-4f;
/* Computing 2nd power */
    d_1 = (*e + 9.9f) / 9.766f;
    fit = exp(d_1 * d_1 * -.5f) * 146300;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf2_sk2_med_ */

/* Subroutine */ int pdf2_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.43541e-4f;
/* Computing 2nd power */
    d_1 = (*e + 9.252f) / 8.36f;
    fit = exp(d_1 * d_1 * -.5f) * 12040;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf2_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf2_half_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 7.3631e-5f;
/* Computing 2nd power */
    d_1 = (*e + 9.76f) / 8.67f;
    fit = exp(d_1 * d_1 * -.5f) * 31900;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf2_half_sk2_low_ */

/* Subroutine */ int pdf2_half_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 7.3631e-5f;
/* Computing 2nd power */
    d_1 = (*e + 17.45f) / 12.5f;
    fit = exp(d_1 * d_1 * -.5f) * 208700;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf2_half_sk2_med_ */

/* Subroutine */ int pdf2_half_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 7.3631e-5f;
/* Computing 2nd power */
    d_1 = (*e + 18.48f) / 11.01f;
    fit = exp(d_1 * d_1 * -.5f) * 20420;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf2_half_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf3_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 4.80922e-4f;
/* Computing 2nd power */
    d_1 = (*e + 21.13f) / 10.097f;
    fit = exp(d_1 * d_1 * -.5f) * 12650;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf3_sk2_low_ */

/* Subroutine */ int pdf3_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 4.80922e-4f;
/* Computing 2nd power */
    d_1 = (*e + 23.77f) / 14.93f;
    fit = exp(d_1 * d_1 * -.5f) * 24100;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf3_sk2_med_ */

/* Subroutine */ int pdf3_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 4.80922e-4f;
/* Computing 2nd power */
    d_1 = (*e + 29.74f) / 12.89f;
    fit = exp(d_1 * d_1 * -.5f) * 12103;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf3_sk2_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf3_half_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.24649e-5f;
/* Computing 2nd power */
    d_1 = (*e + 17.64f) / 11.032f;
    fit = exp(d_1 * d_1 * -.5f) * 48901;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf3_half_sk2_low_ */

/* Subroutine */ int pdf3_half_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.24649e-5f;
/* Computing 2nd power */
    d_1 = (*e + 14.698f) / 14.914f;
    fit = exp(d_1 * d_1 * -.5f) * 63283;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf3_half_sk2_med_ */

/* Subroutine */ int pdf3_half_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 3.24649e-5f;
/* Computing 2nd power */
    d_1 = (*e + 33.876f) / 14.647f;
    fit = exp(d_1 * d_1 * -.5f) * 78701;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf3_half_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf4_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.51304e-4f;
/* Computing 2nd power */
    d_1 = (*e + 6.283f) / 8.116f;
    fit = exp(d_1 * d_1 * -.5f) * 73.32f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf4_sk2_low_ */

/* Subroutine */ int pdf4_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.51304e-4f;
/* Computing 2nd power */
    d_1 = (*e + 9.73f) / 15.22f;
    fit = exp(d_1 * d_1 * -.5f) * 3363;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf4_sk2_med_ */

/* Subroutine */ int pdf4_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.51304e-4f;
/* Computing 2nd power */
    d_1 = (*e + 34.06f) / 17.75f;
    fit = exp(d_1 * d_1 * -.5f) * 1039;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf4_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf4_half_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.00347e-5f;
/* Computing 2nd power */
    d_1 = (*e + 13.988f) / 11.882f;
    fit = exp(d_1 * d_1 * -.5f) * 11300;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf4_half_sk2_low_ */

/* Subroutine */ int pdf4_half_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.00347e-5f;
/* Computing 2nd power */
    d_1 = (*e + 14.225f) / 17.523f;
    fit = exp(d_1 * d_1 * -.5f) * 37051;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf4_half_sk2_med_ */

/* Subroutine */ int pdf4_half_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 2.00347e-5f;
/* Computing 2nd power */
    d_1 = (*e + 43.98f) / 18.473f;
    fit = exp(d_1 * d_1 * -.5f) * 43247;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf4_half_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf5_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.67553e-4f;
/* Computing 2nd power */
    d_1 = (*e + 24.85f) / 13.45f;
    fit = exp(d_1 * d_1 * -.5f) * 4216;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf5_sk2_low_ */

/* Subroutine */ int pdf5_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.67553e-4f;
/* Computing 2nd power */
    d_1 = (*e + 5.288f) / 16.75f;
    fit = exp(d_1 * d_1 * -.5f) * 1860;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf5_sk2_med_ */

/* Subroutine */ int pdf5_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.67553e-4f;
/* Computing 2nd power */
    d_1 = (*e + 39.72f) / 19.08f;
    fit = exp(d_1 * d_1 * -.5f) * 1251.f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf5_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf5_half_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.46195e-5f;
/* Computing 2nd power */
    d_1 = (*e + 16.876f) / 13.094f;
    fit = exp(d_1 * d_1 * -.5f) * 11430;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf5_half_sk2_low_ */

/* Subroutine */ int pdf5_half_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.46195e-5f;
/* Computing 2nd power */
    d_1 = (*e + 5.9453f) / 18.061f;
    fit = exp(d_1 * d_1 * -.5f) * 17350;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf5_half_sk2_med_ */

/* Subroutine */ int pdf5_half_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.46195e-5f;
/* Computing 2nd power */
    d_1 = (*e + 42.442f) / 20.261f;
    fit = exp(d_1 * d_1 * -.5f) * 14434;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf5_half_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf6_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.18783e-5f;
/* Computing 2nd power */
    d_1 = (*e + 30.5f) / 16.48f;
    fit = exp(d_1 * d_1 * -.5f) * 25950;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf6_sk2_low_ */

/* Subroutine */ int pdf6_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.18783e-5f;
/* Computing 2nd power */
    d_1 = (*e + 4.55f) / 21.73f;
/* Computing 2nd power */
    d_2 = (*e - 17.13f) / 9.424f;
    fit = exp(d_1 * d_1 * -.5f) * 7070 + exp(d_2 * d_2 * -.5f) * 2823;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf6_sk2_med_ */

/* Subroutine */ int pdf6_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.18783e-5f;
/* Computing 2nd power */
    d_1 = (*e + 6.227f) / 14.96f;
    fit = exp(d_1 * d_1 * -.5f) * 618.3f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf6_sk2_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf6_half_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.16155e-5f;
/* Computing 2nd power */
    d_1 = (*e + 15.945f) / 14.127f;
    fit = exp(d_1 * d_1 * -.5f) * 5504.5f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf6_half_sk2_low_ */

/* Subroutine */ int pdf6_half_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.16155e-5f;
/* Computing 2nd power */
    d_1 = (*e + 8.3248f) / 22.396f;
/* Computing 2nd power */
    d_2 = (*e - 18.691f) / 8.9412f;
    fit = exp(d_1 * d_1 * -.5f) * 9752 + exp(d_2 * d_2 * -.5f) * 2099.4f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf6_half_sk2_med_ */

/* Subroutine */ int pdf6_half_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.16155e-5f;
/* Computing 2nd power */
    d_1 = (*e + 57.382f) / 25.838f;
    fit = exp(d_1 * d_1 * -.5f) * 10635;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf6_half_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf7_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.05338e-4f;
/* Computing 2nd power */
    d_1 = (*e + 36.6f) / 17.83f;
    fit = exp(d_1 * d_1 * -.5f) * 3121;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf7_sk2_low_ */

/* Subroutine */ int pdf7_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.05338e-4f;
/* Computing 2nd power */
    d_1 = (*e + 38.34f) / 27.4f;
/* Computing 2nd power */
    d_2 = (*e + 32.77f) / 14.2f;
    fit = exp(d_1 * d_1 * -.5f) * 8070 - exp(d_2 * d_2 * -.5f) * 168300;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf7_sk2_med_ */

/* Subroutine */ int pdf7_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.05338e-4f;
/* Computing 2nd power */
    d_1 = (*e + 95.34f) / 36.97f;
    fit = exp(d_1 * d_1 * -.5f) * 1299;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf7_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf7_half_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.84941e-6f;
/* Computing 2nd power */
    d_1 = (*e + 39.254f) / 18.588f;
    fit = exp(d_1 * d_1 * -.5f) * 33224;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf7_half_sk2_low_ */

/* Subroutine */ int pdf7_half_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.84941e-6f;
/* Computing 2nd power */
    d_1 = (*e - 18.414f) / 9.9235f;
/* Computing 2nd power */
    d_2 = (*e - 7.5457f) / 21.343f;
    fit = exp(d_1 * d_1 * -.5f) * 3026 + exp(d_2 * d_2 * -.5f) * 4025;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf7_half_sk2_med_ */

/* Subroutine */ int pdf7_half_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.84941e-6f;
/* Computing 2nd power */
    d_1 = (*e + 43.663f) / 24.034f;
    fit = exp(d_1 * d_1 * -.5f) * 4219.3f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf7_half_sk2_high_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf8_sk2_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.34795e-5f;
/* Computing 2nd power */
    d_1 = (*e - 12.37f) / 7.712f;
    fit = exp(d_1 * d_1 * -.5f) * 36.48f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf8_sk2_low_ */

/* Subroutine */ int pdf8_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.34795e-5f;
/* Computing 2nd power */
    d_1 = (*e + 51.1f) / 32.17f;
/* Computing 2nd power */
    d_2 = (*e + 6.097f) / 11.12f;
    fit = exp(d_1 * d_1 * -.5f) * 9870 - exp(d_2 * d_2 * -.5f) * 3477;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf8_sk2_med_ */

/* Subroutine */ int pdf8_sk2_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.34795e-5f;
/* Computing 2nd power */
    d_1 = (*e + 88.63f) / 35.75f;
    fit = exp(d_1 * d_1 * -.5f) * 1001;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf8_sk2_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cccc SK-III pdfs */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf2_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.01832e-5f;
/* Computing 2nd power */
    d_1 = (*e + 4.111f) / 6.32f;
    fit = exp(d_1 * d_1 * -.5f) * 17350;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_sk3_low_ */

/* Subroutine */ int pdf2_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.01832e-5f;
/* Computing 2nd power */
    d_1 = (*e + 10.7f) / 9.26f;
    fit = exp(d_1 * d_1 * -.5f) * 569000;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_sk3_med_ */

/* Subroutine */ int pdf2_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.01832e-5f;
/* Computing 2nd power */
    d_1 = (*e + 3.f) / 7.f;
    fit = exp(d_1 * d_1 * -.5f) * 5863;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf2_half_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.10378e-5f;
/* Computing 2nd power */
    d_1 = (*e - .992f) / 6.04f;
    fit = exp(d_1 * d_1 * -.5f) * 3104;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_half_sk3_low_ */

/* Subroutine */ int pdf2_half_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.10378e-5f;
/* Computing 2nd power */
    d_1 = (*e + 12.4f) / 10.91f;
    fit = exp(d_1 * d_1 * -.5f) * 320300;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_half_sk3_med_ */

/* Subroutine */ int pdf2_half_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.10378e-5f;
/* Computing 2nd power */
    d_1 = (*e + 13.02f) / 9.58f;
    fit = exp(d_1 * d_1 * -.5f) * 18440;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf2_half_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf3_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.29604e-5f;
/* Computing 2nd power */
    d_1 = (*e + 8.977f) / 7.642f;
    fit = exp(d_1 * d_1 * -.5f) * 24030;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_sk3_low_ */

/* Subroutine */ int pdf3_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.29604e-5f;
/* Computing 2nd power */
    d_1 = (*e + 23.48f) / 13.62f;
    fit = exp(d_1 * d_1 * -.5f) * 585100;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_sk3_med_ */

/* Subroutine */ int pdf3_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 3.29604e-5f;
/* Computing 2nd power */
    d_1 = (*e - 10.85f) / 5.393f;
    fit = exp(d_1 * d_1 * -.5f) * 179.2f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf3_half_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.76459e-5f;
/* Computing 2nd power */
    d_1 = (*e + 12.304f) / 9.8427f;
    fit = exp(d_1 * d_1 * -.5f) * 7402.4f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_half_sk3_low_ */

/* Subroutine */ int pdf3_half_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.76459e-5f;
/* Computing 2nd power */
    d_1 = (*e + 24.266f) / 15.613f;
    fit = exp(d_1 * d_1 * -.5f) * 349140;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_half_sk3_med_ */

/* Subroutine */ int pdf3_half_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.76459e-5f;
/* Computing 2nd power */
    d_1 = (*e + 27.095f) / 13.502f;
    fit = exp(d_1 * d_1 * -.5f) * 33229;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf3_half_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf4_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.45749e-4f;
/* Computing 2nd power */
    d_1 = (*e + 27.f) / 13.68f;
    fit = exp(d_1 * d_1 * -.5f) * 1279;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_sk3_low_ */

/* Subroutine */ int pdf4_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.45749e-4f;
/* Computing 2nd power */
    d_1 = (*e + 24.96f) / 17.09f;
    fit = exp(d_1 * d_1 * -.5f) * 22690;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_sk3_med_ */

/* Subroutine */ int pdf4_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.45749e-4f;
/* Computing 2nd power */
    d_1 = (*e + 25.61f) / 13.92f;
    fit = exp(d_1 * d_1 * -.5f) * 2025;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf4_half_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.24962e-5f;
/* Computing 2nd power */
    d_1 = (*e + 13.154f) / 10.283f;
    fit = exp(d_1 * d_1 * -.5f) * 6980.3f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_half_sk3_low_ */

/* Subroutine */ int pdf4_half_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.24962e-5f;
/* Computing 2nd power */
    d_1 = (*e + 24.942f) / 18.587f;
    fit = exp(d_1 * d_1 * -.5f) * 143260;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_half_sk3_med_ */

/* Subroutine */ int pdf4_half_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.24962e-5f;
/* Computing 2nd power */
    d_1 = (*e + 26.587f) / 16.088f;
    fit = exp(d_1 * d_1 * -.5f) * 6363.8f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf4_half_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf5_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.12112e-4f;
/* Computing 2nd power */
    d_1 = (*e + 57.76f) / 21.73f;
    fit = exp(d_1 * d_1 * -.5f) * 2263;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_sk3_low_ */

/* Subroutine */ int pdf5_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.12112e-4f;
/* Computing 2nd power */
    d_1 = (*e + 23.08f) / 19.59f;
    fit = exp(d_1 * d_1 * -.5f) * 8970;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_sk3_med_ */

/* Subroutine */ int pdf5_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.12112e-4f;
/* Computing 2nd power */
    d_1 = (*e + 43.72f) / 18.86f;
    fit = exp(d_1 * d_1 * -.5f) * 2640;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf5_half_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.00406e-5f;
/* Computing 2nd power */
    d_1 = (*e + 5.4978f) / 9.5758f;
    fit = exp(d_1 * d_1 * -.5f) * 1573.3f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_half_sk3_low_ */

/* Subroutine */ int pdf5_half_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.00406e-5f;
/* Computing 2nd power */
    d_1 = (*e + 19.536f) / 20.216f;
    fit = exp(d_1 * d_1 * -.5f) * 55956;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_half_sk3_med_ */

/* Subroutine */ int pdf5_half_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 1.00406e-5f;
/* Computing 2nd power */
    d_1 = (*e + 29.236f) / 18.327f;
    fit = exp(d_1 * d_1 * -.5f) * 3866.4f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf5_half_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf6_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 9.19493e-6f;
/* Computing 2nd power */
    d_1 = (*e + 18.3f) / 12.88f;
    fit = exp(d_1 * d_1 * -.5f) * 3708;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_sk3_low_ */

/* Subroutine */ int pdf6_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 9.19493e-6f;
/* Computing 2nd power */
    d_1 = (*e - 8.847f) / 11.17f;
/* Computing 2nd power */
    d_2 = (*e - 2.131f) / 18.96f;
    fit = exp(d_1 * d_1 * -.5f) * 7884 + exp(d_2 * d_2 * -.5f) * 5689;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_sk3_med_ */

/* Subroutine */ int pdf6_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 9.19493e-6f;
/* Computing 2nd power */
    d_1 = (*e + 52.64f) / 22.76f;
    fit = exp(d_1 * d_1 * -.5f) * 16910;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf6_half_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.5816e-6f;
/* Computing 2nd power */
    d_1 = (*e + 29.443f) / 13.795f;
    fit = exp(d_1 * d_1 * -.5f) * 26063;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_half_sk3_low_ */

/* Subroutine */ int pdf6_half_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.5816e-6f;
/* Computing 2nd power */
    d_1 = (*e - 12.88f) / 10.445f;
/* Computing 2nd power */
    d_2 = (*e + 4.8189f) / 21.207f;
    fit = exp(d_1 * d_1 * -.5f) * 4547.8f + exp(d_2 * d_2 * -.5f) * 
	    9486.5f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_half_sk3_med_ */

/* Subroutine */ int pdf6_half_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.5816e-6f;
/* Computing 2nd power */
    d_1 = (*e + 6.6456f) / 15.741f;
    fit = exp(d_1 * d_1 * -.5f) * 429.61f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf6_half_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf7_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.05623e-5f;
/* Computing 2nd power */
    d_1 = (*e + 33.3f) / 14.45f;
    fit = exp(d_1 * d_1 * -.5f) * 4927;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_sk3_low_ */

/* Subroutine */ int pdf7_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.05623e-5f;
/* Computing 2nd power */
    d_1 = (*e + 22.27f) / 23.97f;
/* Computing 2nd power */
    d_2 = (*e + 20.14f) / 10.79f;
    fit = exp(d_1 * d_1 * -.5f) * 4316 - exp(d_2 * d_2 * -.5f) * 73930;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_sk3_med_ */

/* Subroutine */ int pdf7_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.05623e-5f;
/* Computing 2nd power */
    d_1 = (*e + 50.42f) / 23.36f;
    fit = exp(d_1 * d_1 * -.5f) * 920.5f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf7_half_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.71323e-6f;
/* Computing 2nd power */
    d_1 = (*e + 32.509f) / 15.64f;
    fit = exp(d_1 * d_1 * -.5f) * 11910;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_half_sk3_low_ */

/* Subroutine */ int pdf7_half_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.71323e-6f;
/* Computing 2nd power */
    d_1 = (*e - 15.378f) / 10.016f;
/* Computing 2nd power */
    d_2 = (*e - 8.0079f) / 20.287f;
    fit = exp(d_1 * d_1 * -.5f) * 4299.5f + exp(d_2 * d_2 * -.5f) * 
	    5028.1f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_half_sk3_med_ */

/* Subroutine */ int pdf7_half_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.71323e-6f;
/* Computing 2nd power */
    d_1 = (*e + 52.159f) / 25.084f;
    fit = exp(d_1 * d_1 * -.5f) * 6443.4f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf7_half_sk3_high_ */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int pdf8_sk3_low_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.35816e-5f;
/* Computing 2nd power */
    d_1 = (*e + 41.05f) / 17.15f;
    fit = exp(d_1 * d_1 * -.5f) * 2110;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf8_sk3_low_ */

/* Subroutine */ int pdf8_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.35816e-5f;
/* Computing 2nd power */
    d_1 = (*e + 54.82f) / 31.85f;
/* Computing 2nd power */
    d_2 = (*e + 25.37f) / 13.72f;
    fit = exp(d_1 * d_1 * -.5f) * 16210 - exp(d_2 * d_2 * -.5f) * 53080;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf8_sk3_med_ */

/* Subroutine */ int pdf8_sk3_high_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.35816e-5f;
/* Computing 2nd power */
    d_1 = (*e + 100.8f) / 34.15f;
    fit = exp(d_1 * d_1 * -.5f) * 4340;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf8_sk3_high_ */

/* Subroutine */ int pdf_mal_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.56887e-5f;
/* Computing 2nd power */
    d_1 = (*e + 26.23f) / 19.73f;
/* Computing 2nd power */
    d_2 = (*e + 13.47f) / .9455f;
    fit = exp(d_1 * d_1 * -.5f) * 25840 - exp(d_2 * d_2 * -.5f) * 34545;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mal_sk1_med_ */

/* Subroutine */ int pdf_ksw_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 4.91498e-5f;
/* Computing 2nd power */
    d_1 = (*e - 1.06f) / 16.5f;
/* Computing 2nd power */
    d_2 = (*e - 24.11f) / 10.25f;
    fit = exp(d_1 * d_1 * -.5f) * 4224.6f - exp(d_2 * d_2 * -.5f) * 
	    437.77f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ksw_sk1_med_ */

/* Subroutine */ int pdf_kawa_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.55934e-5f;
/* Computing 2nd power */
    d_1 = (*e + 73.04f) / 29.63f;
    fit = exp(d_1 * d_1 * -.5f) * 209760;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_kawa_sk1_med_ */

/* Subroutine */ int pdf_woo_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 5.55934e-5f;
/* Computing 2nd power */
    d_1 = (*e + 18.26f) / 20.67f;
/* Computing 2nd power */
    d_2 = (*e + 15.744f) / 7.365f;
    fit = exp(d_1 * d_1 * -.5f) * 7500 + exp(d_2 * d_2 * -.5f) * 404.61f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_woo_sk1_med_ */

/* Subroutine */ int pdf_mal_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 1.12652e-4f;
/* Computing 2nd power */
    d_1 = (*e + 15.25f) / 17.92f;
    fit = exp(d_1 * d_1 * -.5f) * 6715.3f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_mal_sk2_med_ */

/* Subroutine */ int pdf_ksw_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.74586e-5f;
/* Computing 2nd power */
    d_1 = (*e + 22.12f) / 21.85f;
/* Computing 2nd power */
    d_2 = (*e - 17.17f) / 15.26f;
    fit = exp(d_1 * d_1 * -.5f) * 7944.4f - exp(d_2 * d_2 * -.5f) * 
	    187.3f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_ksw_sk2_med_ */

/* Subroutine */ int pdf_mal_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 8.02656e-5f;
/* Computing 2nd power */
    d_1 = (*e + 31.481f) / 20.94f;
    fit = exp(d_1 * d_1 * -.5f) * 23130;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_mal_sk3_med_ */

/* Subroutine */ int pdf_ksw_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 7.045e-5f;
/* Computing 2nd power */
    d_1 = (*e + 34.367f) / 23.072f;
/* Computing 2nd power */
    d_2 = (*e - 20.102f) / 4.508f;
    fit = exp(d_1 * d_1 * -.5f) * 16942 + exp(d_2 * d_2 * -.5f) * 190.13f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_ksw_sk3_med_ */

/* Subroutine */ int pdf_fail_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.34346e-5f;
/* Computing 2nd power */
    d_1 = (*e + 126.68f) / 49.31f;
/* Computing 2nd power */
    d_2 = (*e - .632f) / 12.95f;
    fit = exp(d_1 * d_1 * -.5f) * 40527 + exp(d_2 * d_2 * -.5f) * 2099;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_fail_sk1_med_ */

/* Subroutine */ int pdf_fail_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 8.57219e-5f;
/* Computing 2nd power */
    d_1 = (*e + 119.5f) / 41.22f;
    fit = exp(d_1 * d_1 * -.5f) * 280930;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_fail_sk2_med_ */

/* Subroutine */ int pdf_fail_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.5679e-5f;
/* Computing 2nd power */
    d_1 = (*e + 115.6f) / 40.03f;
    fit = exp(d_1 * d_1 * -.5f) * 331300;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_fail_sk3_med_ */

/* Subroutine */ int pdf_woos_sk1_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.41323e-5f;
/* Computing 2nd power */
    d_1 = (*e + 20.132f) / 20.85f;
/* Computing 2nd power */
    d_2 = (*e - 18.06f) / 6.1f;
    fit = exp(d_1 * d_1 * -.5f) * 7123.1f + exp(d_2 * d_2 * -.5f) * 
	    215.4f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_woos_sk1_med_ */

/* Subroutine */ int pdf_woos_sk2_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1, d_2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s20, s17, s18, s19;
    static double fit;

    spallhigh = .882f;
    spalllow = .762f;
    s17 = .738f;
    s18 = .821f;
    s19 = .878f;
    s20 = .965f;
    c_ = 9.02006e-5f;
/* Computing 2nd power */
    d_1 = (*e + 34.03f) / 22.34f;
/* Computing 2nd power */
    d_2 = (*e - 17.377f) / 2.3f;
    fit = exp(d_1 * d_1 * -.5f) * 22531 - exp(d_2 * d_2 * -.5f) * 269.63f;
    *value = c_ * fit;
    if (*e > 21.2f && *e <= 26.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 20.14f && *e <= 21.2f) {
	*value = spallhigh * c_ * s20 * fit;
    }
    if (*e > 20. && *e <= 20.14f) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 19.08f && *e <= 20.) {
	*value = spalllow * c_ * s19 * fit;
    }
    if (*e > 18.02f && *e <= 19.08f) {
	*value = spalllow * c_ * s18 * fit;
    }
    if (*e > 17.5f && *e <= 18.02f) {
	*value = spalllow * c_ * s17 * fit;
    }
    return 0;
} /* pdf_woos_sk2_med_ */

/* Subroutine */ int pdf_woos_sk3_med_(double *e, double *value)
{
    /* System generated locals */
    double d_1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double spalllow;
    static double c_;
    static double spallhigh, s16, s17, s18, s19;
    static double fit;

    spallhigh = .908f;
    spalllow = .818f;
    s16 = .738f;
    s17 = .821f;
    s18 = .878f;
    s19 = .965f;
    c_ = 6.60453e-5f;
/* Computing 2nd power */
    d_1 = (*e + 17.1f) / 19.78f;
    fit = exp(d_1 * d_1 * -.5f) * 7294.4f;
    *value = c_ * fit;
    if (*e > 20. && *e <= 24.) {
	*value = spallhigh * c_ * fit;
    }
    if (*e > 19. && *e <= 20.) {
	*value = spallhigh * c_ * s19 * fit;
    }
    if (*e > 18. && *e <= 19.) {
	*value = spallhigh * c_ * s18 * fit;
    }
    if (*e > 17. && *e <= 18.) {
	*value = spalllow * c_ * s17 * fit;
    }
    if (*e > 16. && *e <= 17.) {
	*value = spalllow * c_ * s16 * fit;
    }
    return 0;
} /* pdf_woos_sk3_med_ */


double pdf(double en, int sknum, int model, double elow, int ev_type, int region){
    if (model != 0 && ev_type == 4){
        cout << "Model " << model << " not implemented!" << endl;
    }
    double val = 0;
    // Relic
    if (ev_type == 4){
        if (model == 0){
            if (sknum == 1){
                if (region == 0) pdf_rel_sk1_low_(&en, &val);
                if (region == 1) pdf_rel_sk1_med_(&en, &val);
                if (region == 2) pdf_rel_sk1_high_(&en, &val);
            }
            if (sknum == 2){
                if (region == 0) pdf_rel_sk2_low_(&en, &val);
                if (region == 1) pdf_rel_sk2_med_(&en, &val);
                if (region == 2) pdf_rel_sk2_high_(&en, &val);
            }
            if (sknum == 3){
                if (region == 0) pdf_rel_sk3_low_(&en, &val);
                if (region == 1) pdf_rel_sk3_med_(&en, &val);
                if (region == 2) pdf_rel_sk3_high_(&en, &val);
            }
        }
    }
    //nue
    if (ev_type == 0){
        if (sknum == 1){
            if (region == 0) pdf_nue_sk1_low_(&en, &val);
            if (region == 1) pdf_nue_sk1_med_(&en, &val);
            if (region == 2) pdf_nue_sk1_high_(&en, &val);
        }
        if (sknum == 2){
            if (region == 0) pdf_nue_sk2_low_(&en, &val);
            if (region == 1) pdf_nue_sk2_med_(&en, &val);
            if (region == 2) pdf_nue_sk2_high_(&en, &val);
        }
        if (sknum == 3){
            if (region == 0) pdf_nue_sk3_low_(&en, &val);
            if (region == 1) pdf_nue_sk3_med_(&en, &val);
            if (region == 2) pdf_nue_sk3_high_(&en, &val);
        }
    }
    //numu
    if (ev_type == 1){
        if (sknum == 1){
            if (region == 0) pdf_de_sk1_low_(&en, &val);
            if (region == 1) pdf_de_sk1_med_(&en, &val);
            if (region == 2) pdf_de_sk1_high_(&en, &val);
        }
        if (sknum == 2){
            if (region == 0) pdf_de_sk2_low_(&en, &val);
            if (region == 1) pdf_de_sk2_med_(&en, &val);
            if (region == 2) pdf_de_sk2_high_(&en, &val);
        }
        if (sknum == 3){
            if (region == 0) pdf_de_sk3_low_(&en, &val);
            if (region == 1) pdf_de_sk3_med_(&en, &val);
            if (region == 2) pdf_de_sk3_high_(&en, &val);
        }
    }
    //nc
    if (ev_type == 2){
        if (sknum == 1){
            if (region == 0) pdf_ncel_sk1_low_(&en, &val);
            if (region == 1) pdf_ncel_sk1_med_(&en, &val);
            if (region == 2) pdf_ncel_sk1_high_(&en, &val);
        }
        if (sknum == 2){
            if (region == 0) pdf_ncel_sk2_low_(&en, &val);
            if (region == 1) pdf_ncel_sk2_med_(&en, &val);
            if (region == 2) pdf_ncel_sk2_high_(&en, &val);
        }
        if (sknum == 3){
            if (region == 0) pdf_ncel_sk3_low_(&en, &val);
            if (region == 1) pdf_ncel_sk3_med_(&en, &val);
            if (region == 2) pdf_ncel_sk3_high_(&en, &val);
        }
    }
    //mupi
    if (ev_type == 3){
        if (sknum == 1){
            if (region == 0) pdf_mupi_sk1_low_(&en, &val);
            if (region == 1) pdf_mupi_sk1_med_(&en, &val);
            if (region == 2) pdf_mupi_sk1_high_(&en, &val);
        }
        if (sknum == 2){
            if (region == 0) pdf_mupi_sk2_low_(&en, &val);
            if (region == 1) pdf_mupi_sk2_med_(&en, &val);
            if (region == 2) pdf_mupi_sk2_high_(&en, &val);
        }
        if (sknum == 3){
            if (region == 0) pdf_mupi_sk3_low_(&en, &val);
            if (region == 1) pdf_mupi_sk3_med_(&en, &val);
            if (region == 2) pdf_mupi_sk3_high_(&en, &val);
        }
    }
    return val;
}

vector<vector<double>> get_pdfs(vector<double> energies, int sknum, int region, int model, double elow){
    vector<vector<double>> result;
    for(int i = 0;i < energies.size(); i++){
        vector<double> res;
        for(int evtype = 0; evtype < 5; evtype++){
            res.push_back(pdf(energies[i], sknum, model, elow, evtype, region));
        }
        result.push_back(res);
    }
    return result;
}

PYBIND11_MODULE(likes, m){
    m.doc() = "PDFs for SRN spectral analysis";
    m.def("pdf", &pdf);
    m.def("get_pdfs", &get_pdfs);
}
