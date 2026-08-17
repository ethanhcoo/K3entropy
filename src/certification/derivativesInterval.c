/*
 * This program certifies the derivative estimates in Appendix C.2. It reads
 * the pseudo-orbit and chart choices from k3_orbit_data.c and uses the
 * outward-rounded MPFR intervals in certificate_interval.c.
 *
 * It checks four things:
 *   (1) the entries of Table 3 are the correct two-decimal roundings;
 *   (2) K < 115, R < 163, and M < 442;
 *   (3) the discriminant is positive on the axis-aligned squares of radius
 *       10^{-3} around the pseudo-orbit coordinates, and hence on the balls
 *       inside them; and
 *   (4) every entry of L_i - \widetilde L_i has absolute value less than
 *       10^{-2}, where \widetilde L_i is the rational matrix in Table 4.
 *
 * Together with the analytic estimates in the manuscript, the bounds in (2)
 * give |\mathcal D| < 116 and ||D\mathcal D||_{C^1,infinity} < 443 on the
 * 10^{-5}-balls. The manuscript then obtains the bounds 122 and
 * 2.5 * 10^4 for the first and second derivatives of p_+ and p_-; those
 * later estimates are not checked separately by this program.
 *
 * Each finite decimal string is treated as the exact rational number it
 * denotes by parsing it once downward and once upward. An inconclusive
 * interval check stops the program without a certificate.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>

#include "certificate_interval.h"
#include "k3_orbit_data.h"

#define PRECISION 5000

typedef cert_interval interval;

typedef struct {
    interval x, y, z;
} mpfr_point;

/*
 * The exact rational matrices \widetilde L_i displayed in Table 4. A
 * finite decimal string denotes the corresponding exact rational number,
 * not a binary floating-point approximation.
 */
static const char *tilde_L_decimal[K3_ORBIT_LENGTH][2][2] = {
    {{"-0.2159", "-0.3755"}, {"-0.4694", "0.4623"}},
    {{"1.7718", "-2.1227"}, {"-12.3247", "16.3690"}},
    {{"0.3539", "-4.7730"}, {"0.4185", "-4.6570"}},
    {{"-3.1432", "-0.4406"}, {"-0.0916", "-1.1425"}},
    {{"-0.4583", "0.1150"}, {"-0.6163", "0.8316"}},
    {{"1.4772", "-1.9866"}, {"0.3707", "-2.6806"}},
    {{"-0.8852", "0.0258"}, {"-0.1241", "0.3218"}},
    {{"1.0118", "1.1967"}, {"13.6472", "13.3154"}},
    {{"-0.6239", "-4.3400"}, {"0.7475", "5.7641"}},
    {{"-1.3601", "1.6743"}, {"-0.3730", "-2.2037"}}
};

/* Column order for the six quantities displayed in Table 3. */
enum table3_column {
    TABLE3_D,
    TABLE3_DX,
    TABLE3_DY,
    TABLE3_DXX,
    TABLE3_DYY,
    TABLE3_DXY
};

/*
 * Exact decimal entries displayed in Table 3. To certify rounding to two
 * decimal places, the program proves that each exact value is at
 * distance strictly less than 0.005 from the corresponding entry below.
 */
static const char *table3_decimal[K3_ORBIT_LENGTH][6] = {
    {"114.24", "136.59", "-45.98", "-419.46", "-441.48", "-178.78"},
    {"8.71", "-22.00", "15.52", "39.45", "13.84", "-81.14"},
    {"70.23", "86.57", "-124.38", "-58.60", "-2.96", "-171.09"},
    {"8.59", "-40.85", "4.26", "113.57", "-24.69", "-87.82"},
    {"108.31", "39.04", "162.94", "-260.04", "-171.83", "27.57"},
    {"10.43", "23.28", "-23.28", "28.55", "28.55", "-92.99"},
    {"108.31", "-162.94", "-39.04", "-171.83", "-260.04", "27.57"},
    {"8.59", "-4.26", "-40.85", "-24.69", "113.57", "87.82"},
    {"70.23", "124.38", "-86.57", "-2.96", "-58.60", "-171.09"},
    {"8.71", "15.52", "22.00", "13.84", "39.45", "81.14"}
};

_Noreturn static void fail(const char *message)
{
    fprintf(stderr, "CERTIFICATION FAILURE: %s\n", message);
    exit(EXIT_FAILURE);
}

static void require_interval_status(
    cert_interval_status status,
    const char *context
)
{
    if (status != CERT_INTERVAL_OK) {
        fprintf(
            stderr,
            "CERTIFICATION FAILURE: %s (%s)\n",
            context,
            cert_interval_status_string(status)
        );
        exit(EXIT_FAILURE);
    }
}

/*
 * Reject nonfinite or unordered intervals before using them in certificate
 * comparisons.
 */
static void require_valid_interval(const interval *x, const char *context)
{
    if (!cert_interval_valid(x)) {
        fail(context);
    }
}

static void *checked_malloc(size_t size)
{
    void *result = malloc(size);
    if (result == NULL) {
        fail("memory allocation failed");
    }
    return result;
}

static void interval_init(interval *x)
{
    require_interval_status(
        cert_interval_init(x, PRECISION),
        "could not initialize an interval"
    );
}

static void interval_clear(interval *x)
{
    cert_interval_clear(x);
}

static void init_point(mpfr_point *p)
{
    interval_init(&p->x);
    interval_init(&p->y);
    interval_init(&p->z);
}

static void clear_point(mpfr_point *p)
{
    interval_clear(&p->x);
    interval_clear(&p->y);
    interval_clear(&p->z);
}

static void interval_copy(interval *z, const interval *x)
{
    require_interval_status(
        cert_interval_copy(z, x),
        "interval copy failed"
    );
}

static void interval_set_si(interval *x, long value)
{
    require_interval_status(
        cert_interval_set_si(x, value),
        "integer conversion failed"
    );
}

static void interval_set_exact_decimal(interval *x, const char *decimal)
{
    require_interval_status(
        cert_interval_set_exact_decimal(x, decimal),
        "exact-decimal conversion failed"
    );
}

static void interval_add(interval *z, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_add(z, x, y),
        "interval addition failed"
    );
}

static void interval_sub(interval *z, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_sub(z, x, y),
        "interval subtraction failed"
    );
}

static void interval_mul(interval *z, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_mul(z, x, y),
        "interval multiplication failed"
    );
}

static void interval_mul_si(interval *z, const interval *x, long value)
{
    require_interval_status(
        cert_interval_mul_si(z, x, value),
        "signed-integer interval multiplication failed"
    );
}

static void interval_neg(interval *z, const interval *x)
{
    require_interval_status(
        cert_interval_neg(z, x),
        "interval negation failed"
    );
}

static void interval_div(interval *z, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_div(z, x, y),
        "interval division failed"
    );
}

static void interval_sqrt(interval *z, const interval *x)
{
    require_interval_status(
        cert_interval_sqrt(z, x),
        "interval square root failed"
    );
}

static void interval_add_si(interval *z, const interval *x, long value)
{
    require_interval_status(
        cert_interval_add_si(z, x, value),
        "signed-integer interval addition failed"
    );
}

static void interval_ui_sub(interval *z, unsigned long value, const interval *x)
{
    require_interval_status(
        cert_interval_ui_sub(z, value, x),
        "unsigned-integer interval subtraction failed"
    );
}

/*
 * Use reciprocal-then-multiply division here to preserve the established
 * interval enclosures. cert_interval_div_ui uses direct division and can
 * give a different, usually tighter, printed enclosure.
 */
static void interval_div_ui(interval *z, const interval *x, unsigned long value)
{
    interval constant;
    if (value == 0) {
        fail("division by zero");
    }
    interval_init(&constant);
    interval_set_si(&constant, (long)value);
    interval_div(z, x, &constant);
    interval_clear(&constant);
}

/*
 * Use exponentiation by squaring through interval_mul to preserve the
 * established evaluation order and certificate transcript.
 */
static void interval_pow_ui(interval *z, const interval *x, unsigned long power)
{
    interval result;
    interval factor;
    interval temporary;

    require_valid_interval(x, "invalid interval passed to interval power");
    interval_init(&result);
    interval_init(&factor);
    interval_init(&temporary);
    interval_set_si(&result, 1);
    interval_copy(&factor, x);

    while (power > 0) {
        if (power & 1UL) {
            interval_mul(&temporary, &result, &factor);
            interval_copy(&result, &temporary);
        }
        power >>= 1;
        if (power > 0) {
            interval_mul(&temporary, &factor, &factor);
            interval_copy(&factor, &temporary);
        }
    }

    interval_copy(z, &result);
    interval_clear(&result);
    interval_clear(&factor);
    interval_clear(&temporary);
}

static void interval_expand(
    interval *z,
    const interval *center,
    const interval *radius
)
{
    require_interval_status(
        cert_interval_expand(z, center, radius),
        "interval box construction failed"
    );
}

static void interval_abs_upper(mpfr_t result, const interval *x)
{
    require_interval_status(
        cert_interval_abs_upper(result, x),
        "interval absolute-value bound failed"
    );
}

static void update_max_abs(mpfr_t maximum, const interval *x)
{
    mpfr_t candidate;
    mpfr_init2(candidate, PRECISION);
    interval_abs_upper(candidate, x);
    if (mpfr_cmp(candidate, maximum) > 0) {
        mpfr_set(maximum, candidate, MPFR_RNDU);
    }
    mpfr_clear(candidate);
}

static void update_min_lower(mpfr_t minimum, const interval *x)
{
    require_valid_interval(x, "invalid interval in lower-bound accumulator");
    if (mpfr_cmp(x->lo, minimum) < 0) {
        mpfr_set(minimum, x->lo, MPFR_RNDD);
    }
}

static interval **allocate_matrix(int rows, int columns)
{
    interval **matrix = checked_malloc((size_t)rows * sizeof(*matrix));
    for (int i = 0; i < rows; i++) {
        matrix[i] = checked_malloc((size_t)columns * sizeof(*matrix[i]));
        for (int j = 0; j < columns; j++) {
            interval_init(&matrix[i][j]);
            interval_set_si(&matrix[i][j], 0);
        }
    }
    return matrix;
}

static void clear_matrix(interval **matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            interval_clear(&matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}

static void zero_matrix(interval **matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            interval_set_si(&matrix[i][j], 0);
        }
    }
}

/* Interval evaluators for the surface formulas. */
void alphaExact(interval *result, const interval *x, const interval *y);
void betaExact(interval *result, const interval *x, const interval *y);
void DExact(interval *result, const interval *x, const interval *y);
void ix_exact(mpfr_point *p); /* [sigma_1] */
void iy_exact(mpfr_point *p); /* [sigma_2] */
void iz_exact(mpfr_point *p); /* [sigma_3] */
void pPlusExact(interval *result, const interval *x, const interval *y);
void pMinusExact(interval *result, const interval *x, const interval *y);
void lift_exact_charts(
    mpfr_point *p,
    const interval *x,
    const interval *y,
    int i
);

/*
 * Derivative evaluators and matrix operations used in the certifications
 * below.
 */
void DixExact(const mpfr_point *p, interval **matrix); /* [D sigma_1] */
void DiyExact(const mpfr_point *p, interval **matrix); /* [D sigma_2] */
void DizExact(const mpfr_point *p, interval **matrix); /* [D sigma_3] */
void DxExact(interval *result, const interval *x, const interval *y);
void DyExact(interval *result, const interval *x, const interval *y);
void DxxExact(interval *result, const interval *x, const interval *y);
void DyyExact(interval *result, const interval *x, const interval *y);
void DxyExact(interval *result, const interval *x, const interval *y);
void pxExact(interval *result, const interval *s, const interval *t);
void pyExact(interval *result, const interval *s, const interval *t);
void alpha_xExact(
    interval *result,
    const interval *x,
    const interval *y
);
void alpha_yExact(
    interval *result,
    const interval *x,
    const interval *y
);
void multiplyMatricesExact(
    interval **firstMatrix,
    interval **secondMatrix,
    interval **result,
    int ROWS1,
    int COLS1,
    int ROWS2,
    int COLS2
);
void printMatrixExact(interval **matrix, int ROWS1, int COLS1);

static void print_interval(const interval *x)
{
    mpfr_fprintf(stderr, "[%+.17RDe, %+.17RUe]", x->lo, x->hi);
}

/*
 * Increase maximum to an outward-rounded upper bound for
 *
 *     sup_{x in enclosure} |x - q|,
 *
 * where q is the exact rational represented by reference_decimal. Parsing
 * q in both directed modes gives an interval containing q, so subtracting
 * that interval can only enlarge the resulting error bound.
 */
static void update_reference_error(
    mpfr_t maximum,
    const interval *enclosure,
    const char *reference_decimal
)
{
    interval reference;
    interval difference;
    mpfr_t candidate;

    interval_init(&reference);
    interval_init(&difference);
    mpfr_init2(candidate, PRECISION);
    interval_set_exact_decimal(&reference, reference_decimal);
    interval_sub(&difference, enclosure, &reference);
    interval_abs_upper(candidate, &difference);
    if (mpfr_cmp(candidate, maximum) > 0) {
        mpfr_set(maximum, candidate, MPFR_RNDU);
    }

    interval_clear(&reference);
    interval_clear(&difference);
    mpfr_clear(candidate);
}

static void verify_mpfr_status(void)
{
    if (mpfr_nanflag_p()) {
        fail("MPFR reported a NaN");
    }
    if (mpfr_divby0_p()) {
        fail("MPFR reported division by zero");
    }
    if (mpfr_overflow_p()) {
        fail("MPFR reported overflow");
    }
    if (mpfr_underflow_p()) {
        fail("MPFR reported underflow");
    }
    if (mpfr_erangeflag_p()) {
        fail("MPFR reported an exponent-range error");
    }
}

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    /* Every interval-arithmetic operation below supplies directed rounding. */
    mpfr_set_default_rounding_mode(MPFR_RNDN);

    /* Fix the precision and exponent range for this certification. */
    mpfr_set_default_prec(PRECISION);
    if (mpfr_set_emin(-1073) != 0 || mpfr_set_emax(1024) != 0) {
        fail("MPFR rejected the requested exponent range");
    }
    mpfr_clear_flags();

    if (!k3_orbit_data_are_valid()) {
        fail("shared pseudo-orbit/chart data failed its consistency check");
    }

    /* Work matrices for the chain-rule product defining L_i. */
    interval **temp1 = allocate_matrix(3, 3);
    interval **temp2 = allocate_matrix(3, 3);
    interval **temp3 = allocate_matrix(3, 3);
    interval **temp4 = allocate_matrix(3, 3);
    interval **temp5 = allocate_matrix(3, 3);
    interval **x_mat = allocate_matrix(3, 3);
    interval **y_mat = allocate_matrix(3, 3);
    interval **z_mat = allocate_matrix(3, 3);

    /* Table 1 decimals are interpreted as exact rational inputs. */
    interval a[K3_ORBIT_LENGTH];
    interval b[K3_ORBIT_LENGTH];
    for (int i = 0; i < K3_ORBIT_LENGTH; i++) {
        interval_init(&a[i]);
        interval_init(&b[i]);
        interval_set_exact_decimal(&a[i], k3_orbit_a_decimal[i]);
        interval_set_exact_decimal(&b[i], k3_orbit_b_decimal[i]);
    }

    /* Quantities used in the four certifications from Appendix C.2. */
    interval out;
    interval radius;
    interval box_x;
    interval box_y;
    interval box_D;
    interval matrix_error_tolerance;
    interval two_decimal_half_unit;
    mpfr_t K_upper;
    mpfr_t R_upper;
    mpfr_t M_upper;
    mpfr_t minimum_box_D;
    mpfr_t maximum_L_error;
    mpfr_t maximum_table3_error;

    interval_init(&out);
    interval_init(&radius);
    interval_init(&box_x);
    interval_init(&box_y);
    interval_init(&box_D);
    interval_init(&matrix_error_tolerance);
    interval_init(&two_decimal_half_unit);
    interval_set_exact_decimal(&radius, "1e-3");
    interval_set_exact_decimal(&matrix_error_tolerance, "1e-2");
    interval_set_exact_decimal(&two_decimal_half_unit, "5e-3");

    mpfr_init2(K_upper, PRECISION);
    mpfr_init2(R_upper, PRECISION);
    mpfr_init2(M_upper, PRECISION);
    mpfr_init2(minimum_box_D, PRECISION);
    mpfr_init2(maximum_L_error, PRECISION);
    mpfr_init2(maximum_table3_error, PRECISION);
    mpfr_set_zero(K_upper, 1);
    mpfr_set_zero(R_upper, 1);
    mpfr_set_zero(M_upper, 1);
    mpfr_set_inf(minimum_box_D, 1);
    mpfr_set_zero(maximum_L_error, 1);
    mpfr_set_zero(maximum_table3_error, 1);

    fprintf(stderr, "MPFR version: %s\n", mpfr_get_version());
    fprintf(stderr, "precision: %d bits\n", PRECISION);
    fprintf(
        stderr,
        "Every printed interval is outward-rounded and contains the exact value.\n\n"
    );

    for (int i = 0; i < K3_ORBIT_LENGTH; i++) {
        fprintf(
            stderr,
            "(x,y) is (%s, %s)\n",
            k3_orbit_a_decimal[i],
            k3_orbit_b_decimal[i]
        );

        /* Print the encoded Table 1 chart choice beside its enclosures. */
        fprintf(stderr, "chart is %d\n", k3_orbit_chart_code[i]);

        /*
         * Enclose the six exact quantities whose two-decimal roundings are
         * displayed in Table 3. In addition to updating K, R, and M, compare
         * each enclosure directly with the exact rational printed there.
         */
        DExact(&out, &a[i], &b[i]);
        update_max_abs(K_upper, &out);
        update_reference_error(
            maximum_table3_error,
            &out,
            table3_decimal[i][TABLE3_D]
        );
        fprintf(stderr, "DExact(x,y) is ");
        print_interval(&out);
        fprintf(stderr, "\n");

        DxExact(&out, &a[i], &b[i]);
        update_max_abs(R_upper, &out);
        update_reference_error(
            maximum_table3_error,
            &out,
            table3_decimal[i][TABLE3_DX]
        );
        fprintf(stderr, "DxExact(x,y) is ");
        print_interval(&out);
        fprintf(stderr, "\n");

        DyExact(&out, &a[i], &b[i]);
        update_max_abs(R_upper, &out);
        update_reference_error(
            maximum_table3_error,
            &out,
            table3_decimal[i][TABLE3_DY]
        );
        fprintf(stderr, "DyExact(x,y) is ");
        print_interval(&out);
        fprintf(stderr, "\n");

        DxxExact(&out, &a[i], &b[i]);
        update_max_abs(M_upper, &out);
        update_reference_error(
            maximum_table3_error,
            &out,
            table3_decimal[i][TABLE3_DXX]
        );
        fprintf(stderr, "DxxExact(x,y) is ");
        print_interval(&out);
        fprintf(stderr, "\n");

        DyyExact(&out, &a[i], &b[i]);
        update_max_abs(M_upper, &out);
        update_reference_error(
            maximum_table3_error,
            &out,
            table3_decimal[i][TABLE3_DYY]
        );
        fprintf(stderr, "DyyExact(x,y) is ");
        print_interval(&out);
        fprintf(stderr, "\n");

        DxyExact(&out, &a[i], &b[i]);
        update_max_abs(M_upper, &out);
        update_reference_error(
            maximum_table3_error,
            &out,
            table3_decimal[i][TABLE3_DXY]
        );
        fprintf(stderr, "DxyExact(x,y) is ");
        print_interval(&out);
        fprintf(stderr, "\n");

        /*
         * The square box contains B_{10^{-3}}(a_i,b_i), so a positive
         * lower bound here certifies positivity on the Euclidean ball.
         */
        interval_expand(&box_x, &a[i], &radius);
        interval_expand(&box_y, &b[i], &radius);
        DExact(&box_D, &box_x, &box_y);
        update_min_lower(minimum_box_D, &box_D);
        fprintf(stderr, "D on the 1e-3 box is ");
        print_interval(&box_D);
        fprintf(stderr, "\n");

        zero_matrix(temp4, 3, 3);
        zero_matrix(temp5, 3, 3);

        /*
         * temp4 encloses D\phi_i(0) = D\widetilde\phi_i(a_i,b_i).
         * For a negative chart we obtain the derivatives from
         * p_-(s,t) = -p_+(-s,t).
         */
        int k = 0;
        for (int j = 0; j < 3; j++) {
            if (j == 3 - k3_orbit_chart_code[i]) {
                pxExact(&temp4[j][0], &a[i], &b[i]);
                pyExact(&temp4[j][1], &a[i], &b[i]);
            } else if (j == 6 - k3_orbit_chart_code[i]) {
                interval negative_a;
                interval_init(&negative_a);

                /* Differentiate p_-(s,t) = -p_+(-s,t) in s. */
                interval_neg(&negative_a, &a[i]);
                pxExact(&temp4[j][0], &negative_a, &b[i]);

                /* Differentiate p_-(s,t) = -p_+(-s,t) in t. */
                pyExact(&temp4[j][1], &negative_a, &b[i]);
                interval_neg(&temp4[j][1], &temp4[j][1]);
                interval_clear(&negative_a);
            } else {
                if (k == 0) {
                    interval_set_si(&temp4[j][0], 1);
                    interval_set_si(&temp4[j][1], 0);
                }
                if (k == 1) {
                    interval_set_si(&temp4[j][0], 0);
                    interval_set_si(&temp4[j][1], 1);
                }
                if (k == 2) {
                    fail("invalid input-chart derivative");
                }
                k++;
            }
        }

        /*
         * temp5 is D\xi_{i+1}, where \xi_{i+1}: R^3 -> R^2 is the
         * affine coordinate projection used in Section 3.5.1.
         */
        k = 0;
        for (int j = 0; j < 3; j++) {
            if (
                j == 3 - k3_orbit_chart_code[(i + 1) % K3_ORBIT_LENGTH]
            ) {
                interval_set_si(&temp5[0][j], 0);
                interval_set_si(&temp5[1][j], 0);
            } else if (
                j == 6 - k3_orbit_chart_code[(i + 1) % K3_ORBIT_LENGTH]
            ) {
                interval_set_si(&temp5[0][j], 0);
                interval_set_si(&temp5[1][j], 0);
            } else {
                if (k == 0) {
                    interval_set_si(&temp5[0][j], 1);
                    interval_set_si(&temp5[1][j], 0);
                }
                if (k == 1) {
                    interval_set_si(&temp5[0][j], 0);
                    interval_set_si(&temp5[1][j], 1);
                }
                if (k == 2) {
                    fail("invalid output-chart derivative");
                }
                k++;
            }
        }

        /* Enclose x_i = \widetilde\phi_i(a_i,b_i) on X(R). */
        mpfr_point p;
        init_point(&p);
        lift_exact_charts(
            &p,
            &a[i],
            &b[i],
            k3_orbit_chart_code[i]
        );

        /*
         * f = sigma_3 o sigma_2 o sigma_1.  Evaluate the three derivative
         * factors at x_i, sigma_1(x_i), and sigma_2(sigma_1(x_i)),
         * respectively.
         */
        DixExact(&p, x_mat);
        ix_exact(&p);
        DiyExact(&p, y_mat);
        iy_exact(&p);
        DizExact(&p, z_mat);
        iz_exact(&p);

        /* temp3 = Df_{x_i} = Dsigma_3 Dsigma_2 Dsigma_1. */
        multiplyMatricesExact(y_mat, x_mat, temp2, 3, 3, 3, 3);
        multiplyMatricesExact(z_mat, temp2, temp3, 3, 3, 3, 3);

        /*
         * temp2 = L_i = (Df_i^c)_0 = Dxi_{i+1} Df_{x_i} Dphi_i(0).
         */
        multiplyMatricesExact(temp5, temp3, temp1, 2, 3, 3, 3);
        multiplyMatricesExact(temp1, temp4, temp2, 2, 3, 3, 2);
        fprintf(stderr, "L_%d = (Df_%d^c)_0 is enclosed by\n", i, i);
        printMatrixExact(temp2, 2, 2);

        /*
         * Compare each entry of L_i with the corresponding exact
         * four-decimal entry of \widetilde L_i in Table 4.
         */
        for (int row = 0; row < 2; row++) {
            for (int column = 0; column < 2; column++) {
                update_reference_error(
                    maximum_L_error,
                    &temp2[row][column],
                    tilde_L_decimal[i][row][column]
                );
            }
        }
        fprintf(stderr, "\n");
        clear_point(&p);
    }

    /* Fatal MPFR flags must be rejected before any comparison can certify. */
    verify_mpfr_status();

    fprintf(stderr, "Computed enclosure summary\n");
    mpfr_fprintf(stderr, "K upper bound: %.40RUe\n", K_upper);
    mpfr_fprintf(stderr, "R upper bound: %.40RUe\n", R_upper);
    mpfr_fprintf(stderr, "M upper bound: %.40RUe\n", M_upper);
    mpfr_fprintf(
        stderr,
        "minimum 1e-3-box discriminant lower bound: %.40RDe\n",
        minimum_box_D
    );
    mpfr_fprintf(
        stderr,
        "maximum error from a two-decimal Table 3 entry: %.40RUe\n",
        maximum_table3_error
    );
    mpfr_fprintf(
        stderr,
        "maximum entrywise error from Table 4 tilde(L_i): %.40RUe\n",
        maximum_L_error
    );

    /* Table 3 rounding: strict error < 0.005 also rules out a tie. */
    if (mpfr_cmp(maximum_table3_error, two_decimal_half_unit.lo) >= 0) {
        fail("the two-decimal roundings in Table 3 were not certified");
    }
    if (mpfr_cmp_ui(K_upper, 115) >= 0) {
        fail("K < 115 was not certified");
    }
    if (mpfr_cmp_ui(R_upper, 163) >= 0) {
        fail("R < 163 was not certified");
    }
    if (mpfr_cmp_ui(M_upper, 442) >= 0) {
        fail("M < 442 was not certified");
    }
    if (mpfr_sgn(minimum_box_D) <= 0) {
        fail("discriminants were not certified positive on all 1e-3 boxes");
    }
    if (mpfr_cmp(maximum_L_error, matrix_error_tolerance.lo) >= 0) {
        fail("matrix entrywise error < 1e-2 was not certified");
    }

    fprintf(
        stderr,
        "\nCertified conclusions\n"
        "PASS: Table 3 contains the correct two-decimal roundings.\n"
        "PASS: K < 115, R < 163, M < 442.\n"
        "PASS: D > 0 on every enclosing 1e-3 square, hence on each ball.\n"
        "PASS: every entry of L_i - tilde(L_i) has absolute value < 1e-2.\n"
    );

    for (int i = 0; i < K3_ORBIT_LENGTH; i++) {
        interval_clear(&a[i]);
        interval_clear(&b[i]);
    }
    interval_clear(&out);
    interval_clear(&radius);
    interval_clear(&box_x);
    interval_clear(&box_y);
    interval_clear(&box_D);
    interval_clear(&matrix_error_tolerance);
    interval_clear(&two_decimal_half_unit);
    mpfr_clear(K_upper);
    mpfr_clear(R_upper);
    mpfr_clear(M_upper);
    mpfr_clear(minimum_box_D);
    mpfr_clear(maximum_L_error);
    mpfr_clear(maximum_table3_error);

    clear_matrix(temp1, 3, 3);
    clear_matrix(temp2, 3, 3);
    clear_matrix(temp3, 3, 3);
    clear_matrix(temp4, 3, 3);
    clear_matrix(temp5, 3, 3);
    clear_matrix(x_mat, 3, 3);
    clear_matrix(y_mat, 3, 3);
    clear_matrix(z_mat, 3, 3);

    mpfr_free_cache();
    return EXIT_SUCCESS;
}

/* alpha(x,y) = xy/((1+x^2)(1+y^2)). */
void alphaExact(interval *result, const interval *x, const interval *y)
{
    interval temp1, temp2, temp3;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);

    interval_mul(&temp1, x, y);

    interval_mul(&temp2, x, x);
    interval_add_si(&temp2, &temp2, 1);

    interval_mul(&temp3, y, y);
    interval_add_si(&temp3, &temp3, 1);

    interval_mul(&temp2, &temp2, &temp3);
    interval_div(result, &temp1, &temp2);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
}

/* beta(x,y) = 1/((1+x^2)(1+y^2)). */
void betaExact(interval *result, const interval *x, const interval *y)
{
    interval temp2, temp3, one;
    interval_init(&temp2);
    interval_init(&temp3);
    interval_init(&one);

    interval_mul(&temp2, x, x);
    interval_add_si(&temp2, &temp2, 1);

    interval_mul(&temp3, y, y);
    interval_add_si(&temp3, &temp3, 1);

    interval_mul(&temp2, &temp2, &temp3);

    interval_set_si(&one, 1);
    interval_div(result, &one, &temp2);

    interval_clear(&temp2);
    interval_clear(&temp3);
    interval_clear(&one);
}

void ix_exact(mpfr_point *p)
{
    mpfr_point q;
    init_point(&q);

    interval_copy(&q.x, &p->y);
    interval_copy(&q.y, &p->z);
    interval_copy(&q.z, &p->x);

    iz_exact(&q);
    interval_copy(&p->x, &q.z);
    clear_point(&q);
}

void iy_exact(mpfr_point *p)
{
    mpfr_point q;
    init_point(&q);

    interval_copy(&q.x, &p->x);
    interval_copy(&q.y, &p->z);
    interval_copy(&q.z, &p->y);

    iz_exact(&q);
    interval_copy(&p->y, &q.z);
    clear_point(&q);
}

/* sigma_3(x,y,z) = (x,y,-z-10 alpha(x,y)). */
void iz_exact(mpfr_point *p)
{
    interval temp1;
    interval_init(&temp1);

    alphaExact(&temp1, &p->x, &p->y);
    interval_mul_si(&temp1, &temp1, -10);
    interval_sub(&p->z, &temp1, &p->z);
    interval_clear(&temp1);
}

/*
 * D(x,y) = 100x^2y^2 + 8(1+x^2)(1+y^2)
 *          - 4(1+x^2)^2(1+y^2)^2.
 */
void DExact(interval *result, const interval *x, const interval *y)
{
    interval temp1, temp2, temp3, temp4;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);
    interval_init(&temp4);

    interval_mul(&temp1, x, x);
    interval_mul(&temp2, y, y);
    interval_mul(&temp1, &temp1, &temp2);
    interval_mul_si(&temp1, &temp1, 100);

    interval_mul(&temp2, x, x);
    interval_add_si(&temp2, &temp2, 1);
    interval_mul(&temp3, y, y);
    interval_add_si(&temp3, &temp3, 1);
    interval_mul(&temp2, &temp2, &temp3);
    interval_mul_si(&temp2, &temp2, 8);

    interval_mul(&temp3, x, x);
    interval_add_si(&temp3, &temp3, 1);
    interval_mul(&temp3, &temp3, &temp3);

    interval_mul(&temp4, y, y);
    interval_add_si(&temp4, &temp4, 1);
    interval_mul(&temp4, &temp4, &temp4);

    interval_mul(&temp3, &temp3, &temp4);
    interval_mul_si(&temp3, &temp3, 4);

    interval_add(result, &temp1, &temp2);
    interval_sub(result, result, &temp3);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
    interval_clear(&temp4);
}

/* p_+(x,y) = -5 alpha(x,y) + beta(x,y)sqrt(D(x,y))/2. */
void pPlusExact(interval *result, const interval *x, const interval *y)
{
    interval temp1, temp2, temp3;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);

    alphaExact(&temp1, x, y);
    interval_mul_si(&temp1, &temp1, -5);

    betaExact(&temp2, x, y);
    DExact(&temp3, x, y);
    interval_sqrt(&temp3, &temp3);
    interval_mul(&temp2, &temp2, &temp3);
    interval_div_ui(&temp2, &temp2, 2);

    interval_add(result, &temp1, &temp2);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
}

/* p_-(x,y) = -5 alpha(x,y) - beta(x,y)sqrt(D(x,y))/2. */
void pMinusExact(interval *result, const interval *x, const interval *y)
{
    interval temp1, temp2, temp3;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);

    alphaExact(&temp1, x, y);
    interval_mul_si(&temp1, &temp1, -5);

    betaExact(&temp2, x, y);
    DExact(&temp3, x, y);
    interval_sqrt(&temp3, &temp3);
    interval_mul(&temp2, &temp2, &temp3);
    interval_div_ui(&temp2, &temp2, 2);

    interval_sub(result, &temp1, &temp2);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
}

void lift_exact_charts(
    mpfr_point *p,
    const interval *a,
    const interval *b,
    int i
)
{
    interval temp0;
    interval_init(&temp0);

    if (i == 1) {
        pPlusExact(&temp0, a, b);
        interval_copy(&p->x, a);
        interval_copy(&p->y, b);
        interval_copy(&p->z, &temp0);
    }
    if (i == 2) {
        pPlusExact(&temp0, a, b);
        interval_copy(&p->x, a);
        interval_copy(&p->y, &temp0);
        interval_copy(&p->z, b);
    }
    if (i == 3) {
        pPlusExact(&temp0, a, b);
        interval_copy(&p->x, &temp0);
        interval_copy(&p->y, a);
        interval_copy(&p->z, b);
    }
    if (i == 4) {
        pMinusExact(&temp0, a, b);
        interval_copy(&p->x, a);
        interval_copy(&p->y, b);
        interval_copy(&p->z, &temp0);
    }
    if (i == 5) {
        pMinusExact(&temp0, a, b);
        interval_copy(&p->x, a);
        interval_copy(&p->y, &temp0);
        interval_copy(&p->z, b);
    }
    if (i == 6) {
        pMinusExact(&temp0, a, b);
        interval_copy(&p->x, &temp0);
        interval_copy(&p->y, a);
        interval_copy(&p->z, b);
    }
    if (i < 1 || i > 6) {
        fail("invalid chart number");
    }

    interval_clear(&temp0);
}

void multiplyMatricesExact(
    interval **firstMatrix,
    interval **secondMatrix,
    interval **result,
    int ROWS1,
    int COLS1,
    int ROWS2,
    int COLS2
)
{
    if (COLS1 != ROWS2) {
        fail("matrix dimensions do not match");
    }
    for (int i = 0; i < ROWS1; i++) {
        for (int j = 0; j < COLS2; j++) {
            interval_set_si(&result[i][j], 0);
            for (int k = 0; k < COLS1; k++) {
                interval temp;
                interval_init(&temp);
                interval_mul(
                    &temp,
                    &firstMatrix[i][k],
                    &secondMatrix[k][j]
                );
                interval_add(&result[i][j], &result[i][j], &temp);
                interval_clear(&temp);
            }
        }
    }
}

void printMatrixExact(interval **matrix, int ROWS1, int COLS1)
{
    fprintf(stderr, "{");
    for (int i = 0; i < ROWS1; i++) {
        fprintf(stderr, "{");
        for (int j = 0; j < COLS1; j++) {
            print_interval(&matrix[i][j]);
            if (j < COLS1 - 1) {
                fprintf(stderr, ",");
            }
        }
        if (i == ROWS1 - 1) {
            fprintf(stderr, "}}\n");
        } else {
            fprintf(stderr, "},\n");
        }
    }
}

void alpha_xExact(
    interval *result,
    const interval *x,
    const interval *y
)
{
    interval temp1, temp2, temp3;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);

    /* partial_x alpha = y(1-x^2)/((1+x^2)^2(1+y^2)). */
    interval_mul(&temp1, y, y);
    interval_add_si(&temp1, &temp1, 1);

    interval_mul(&temp2, x, x);
    interval_add_si(&temp2, &temp2, 1);
    interval_mul(&temp2, &temp2, &temp2);

    interval_mul(&temp1, &temp1, &temp2);

    interval_mul(&temp2, x, x);
    interval_ui_sub(&temp2, 1, &temp2);

    interval_mul(&temp2, &temp2, y);
    interval_div(result, &temp2, &temp1);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
}

void alpha_yExact(
    interval *result,
    const interval *x,
    const interval *y
)
{
    alpha_xExact(result, y, x);
}

void DixExact(const mpfr_point *p, interval **matrix)
{
    interval_set_si(&matrix[0][0], -1);
    interval temp1, temp2, temp3;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);

    /* matrix[0][1] = -10 * partial_y alpha(p.y,p.z). */
    alpha_xExact(&temp1, &p->y, &p->z);
    interval_mul_si(&matrix[0][1], &temp1, -10);

    /* matrix[0][2] = -10 * partial_z alpha(p.y,p.z). */
    alpha_yExact(&temp2, &p->y, &p->z);
    interval_mul_si(&matrix[0][2], &temp2, -10);

    interval_set_si(&matrix[1][0], 0);
    interval_set_si(&matrix[1][1], 1);
    interval_set_si(&matrix[1][2], 0);
    interval_set_si(&matrix[2][0], 0);
    interval_set_si(&matrix[2][1], 0);
    interval_set_si(&matrix[2][2], 1);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
}

void DiyExact(const mpfr_point *p, interval **matrix)
{
    interval_set_si(&matrix[0][0], 1);
    interval_set_si(&matrix[0][1], 0);
    interval_set_si(&matrix[0][2], 0);

    interval temp1, temp2, temp3;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);

    /* matrix[1][0] = -10 * partial_x alpha(p.x,p.z). */
    alpha_xExact(&temp1, &p->x, &p->z);
    interval_mul_si(&matrix[1][0], &temp1, -10);

    interval_set_si(&matrix[1][1], -1);

    /* matrix[1][2] = -10 * partial_z alpha(p.x,p.z). */
    alpha_yExact(&temp2, &p->x, &p->z);
    interval_mul_si(&matrix[1][2], &temp2, -10);

    interval_set_si(&matrix[2][0], 0);
    interval_set_si(&matrix[2][1], 0);
    interval_set_si(&matrix[2][2], 1);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
}

void DizExact(const mpfr_point *p, interval **matrix)
{
    interval_set_si(&matrix[0][0], 1);
    interval_set_si(&matrix[0][1], 0);
    interval_set_si(&matrix[0][2], 0);
    interval_set_si(&matrix[1][0], 0);
    interval_set_si(&matrix[1][1], 1);
    interval_set_si(&matrix[1][2], 0);

    interval temp1, temp2, temp3;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);

    /* matrix[2][0] = -10 * partial_x alpha(p.x,p.y). */
    alpha_xExact(&temp1, &p->x, &p->y);
    interval_mul_si(&matrix[2][0], &temp1, -10);

    /* matrix[2][1] = -10 * partial_y alpha(p.x,p.y). */
    alpha_yExact(&temp1, &p->x, &p->y);
    interval_mul_si(&matrix[2][1], &temp1, -10);

    interval_set_si(&matrix[2][2], -1);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
}

/*
 * partial_x D = 200xy^2 + 16x(1+y^2)
 *               - 16x(1+x^2)(1+y^2)^2.
 */
void DxExact(interval *result, const interval *x, const interval *y)
{
    interval temp1, temp2, temp3, temp4;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);
    interval_init(&temp4);

    interval_mul(&temp1, x, y);
    interval_mul(&temp1, &temp1, y);
    interval_mul_si(&temp1, &temp1, 200);

    interval_mul(&temp2, y, y);
    interval_add_si(&temp2, &temp2, 1);
    interval_mul(&temp2, &temp2, x);
    interval_mul_si(&temp2, &temp2, 16);

    interval_mul(&temp3, x, x);
    interval_add_si(&temp3, &temp3, 1);
    interval_mul(&temp4, y, y);
    interval_add_si(&temp4, &temp4, 1);
    interval_mul(&temp4, &temp4, &temp4);
    interval_mul(&temp3, &temp3, &temp4);
    interval_mul(&temp3, &temp3, x);
    interval_mul_si(&temp3, &temp3, 16);

    interval_add(result, &temp1, &temp2);
    interval_sub(result, result, &temp3);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
    interval_clear(&temp4);
}

void DyExact(interval *result, const interval *x, const interval *y)
{
    DxExact(result, y, x);
}

/*
 * partial_xx D = 200y^2 + 16(1+y^2) - 32x^2(1+y^2)^2
 *                - 16(1+x^2)(1+y^2)^2.
 */
void DxxExact(interval *result, const interval *x, const interval *y)
{
    interval temp1, temp2, temp3, temp4, temp5;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);
    interval_init(&temp4);
    interval_init(&temp5);

    interval_mul(&temp1, y, y);
    interval_mul_si(&temp1, &temp1, 200);

    interval_mul(&temp2, y, y);
    interval_add_si(&temp2, &temp2, 1);
    interval_mul_si(&temp2, &temp2, 16);

    interval_mul(&temp3, x, x);
    interval_mul(&temp4, y, y);
    interval_add_si(&temp4, &temp4, 1);
    interval_mul(&temp4, &temp4, &temp4);
    interval_mul(&temp3, &temp3, &temp4);
    interval_mul_si(&temp3, &temp3, 32);

    interval_mul(&temp4, x, x);
    interval_add_si(&temp4, &temp4, 1);
    interval_mul(&temp5, y, y);
    interval_add_si(&temp5, &temp5, 1);
    interval_mul(&temp5, &temp5, &temp5);
    interval_mul(&temp4, &temp4, &temp5);
    interval_mul_si(&temp4, &temp4, 16);

    interval_add(result, &temp1, &temp2);
    interval_sub(result, result, &temp3);
    interval_sub(result, result, &temp4);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
    interval_clear(&temp4);
    interval_clear(&temp5);
}

void DyyExact(interval *result, const interval *x, const interval *y)
{
    DxxExact(result, y, x);
}

/* partial_xy D = 432xy - 64x(1+x^2)y(1+y^2). */
void DxyExact(interval *result, const interval *x, const interval *y)
{
    interval temp1, temp2, temp3;
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);

    interval_mul(&temp1, x, y);
    interval_mul_si(&temp1, &temp1, 432);

    interval_mul(&temp2, x, x);
    interval_add_si(&temp2, &temp2, 1);
    interval_mul(&temp2, &temp2, x);
    interval_mul(&temp2, &temp2, y);
    interval_mul(&temp3, y, y);
    interval_add_si(&temp3, &temp3, 1);
    interval_mul(&temp2, &temp2, &temp3);
    interval_mul_si(&temp2, &temp2, 64);

    interval_sub(result, &temp1, &temp2);

    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
}

/*
 * Enclose partial_s p_+(s,t). Set
 *
 *     S(s,t) = D(s,t)/4,
 *
 * Thus p_+ = -5 alpha + beta sqrt(S). After collecting terms,
 * partial_s p_+ is the quotient with numerator
 *
 *     s(-2+23t^2) - s^3(2+27t^2) - 5t sqrt(S) + 5s^2t sqrt(S)
 *
 * and denominator
 *
 *     (1+s^2)^2 (1+t^2) sqrt(S(s,t)).
 */

void pxExact(interval *result, const interval *s, const interval *t)
{
    interval numerator, denominator, temp1, temp2, temp3, temp4;
    interval temp5, temp6, temp7;
    interval_init(&numerator);
    interval_init(&denominator);
    interval_init(&temp1);
    interval_init(&temp2);
    interval_init(&temp3);
    interval_init(&temp4);
    interval_init(&temp5);
    interval_init(&temp6);
    interval_init(&temp7);

    /* numerator = s(-2+23t^2) - s^3(2+27t^2). */
    interval_mul(&temp1, t, t);
    interval_mul_si(&temp1, &temp1, 23);
    interval_add_si(&temp1, &temp1, -2);
    interval_mul(&temp1, &temp1, s);

    interval_mul(&temp2, s, s);
    interval_mul(&temp2, &temp2, s);

    interval_mul_si(&temp3, t, 27);
    interval_mul(&temp3, t, &temp3);
    interval_add_si(&temp3, &temp3, 2);

    interval_mul(&temp2, &temp2, &temp3);
    interval_sub(&temp1, &temp1, &temp2);

    /* temp4 = 1 - t^4 - s^4(1+t^2)^2. */
    interval_pow_ui(&temp4, t, 4);
    interval_pow_ui(&temp5, s, 4);
    interval_mul(&temp6, t, t);
    interval_add_si(&temp6, &temp6, 1);
    interval_mul(&temp6, &temp6, &temp6);
    interval_mul(&temp5, &temp5, &temp6);
    interval_ui_sub(&temp4, 1, &temp4);
    interval_sub(&temp4, &temp4, &temp5);

    /* Add s^2(23t^2-2t^4) to temp4. */
    interval_mul(&temp5, s, s);
    interval_mul_si(&temp6, t, 23);
    interval_mul(&temp6, &temp6, t);
    interval_pow_ui(&temp7, t, 4);
    interval_mul_si(&temp7, &temp7, 2);
    interval_sub(&temp6, &temp6, &temp7);
    interval_mul(&temp6, &temp6, s);
    interval_mul(&temp6, &temp6, s);
    interval_add(&temp4, &temp4, &temp6);

    /* temp4 = sqrt(S). */
    interval_sqrt(&temp4, &temp4);

    /* Add -5t sqrt(S) + 5s^2t sqrt(S) to the numerator. */
    interval_mul(&temp5, t, &temp4);
    interval_mul_si(&temp5, &temp5, -5);
    interval_add(&numerator, &temp1, &temp5);

    interval_mul(&temp5, s, s);
    interval_mul(&temp5, &temp5, t);
    interval_mul(&temp5, &temp5, &temp4);
    interval_mul_si(&temp5, &temp5, 5);
    interval_add(&numerator, &numerator, &temp5);

    /* denominator = (1+s^2)^2(1+t^2)sqrt(S). */
    interval_mul(&temp1, s, s);
    interval_add_si(&temp1, &temp1, 1);
    interval_mul(&temp1, &temp1, &temp1);

    interval_mul(&temp2, t, t);
    interval_add_si(&temp2, &temp2, 1);

    interval_mul(&denominator, &temp1, &temp2);
    interval_mul(&denominator, &denominator, &temp4);

    interval_div(result, &numerator, &denominator);

    interval_clear(&numerator);
    interval_clear(&denominator);
    interval_clear(&temp1);
    interval_clear(&temp2);
    interval_clear(&temp3);
    interval_clear(&temp4);
    interval_clear(&temp5);
    interval_clear(&temp6);
    interval_clear(&temp7);
}

void pyExact(interval *result, const interval *s, const interval *t)
{
    pxExact(result, t, s);
}
