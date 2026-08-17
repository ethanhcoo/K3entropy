/*
 * Certificate for Appendix C.3 of
 * "A Real K3 Automorphism with Most of Its Entropy in the Real Part".
 *
 * This file is the outward-rounded counterpart of periodicExact.c.  It uses
 * the same ten decimal points, the same six-valued chart encoding, and the
 * same chart-lift -> f -> output-projection organization.  The difference is
 * that every formula evaluation is enclosed by an MPFR interval, with lower
 * endpoints rounded toward -infinity and upper endpoints toward +infinity;
 * scalar certificate bounds are rounded in their conservative direction.
 *
 * For each i, the certification program does all of the following:
 *
 *   1. lifts (a_i,b_i) through the chart listed in Table 1;
 *   2. encloses f = sigma_3 o sigma_2 o sigma_1 on that lifted point;
 *   3. proves that the image lies in the required branch of the next chart;
 *   4. encloses f_i^c(0) in two residual intervals;
 *   5. proves the corresponding Appendix C.3 table row and
 *      ||f_i^c(0)||_2 < 10^(-29).
 *
 * Every finite decimal in this file denotes the exact rational represented
 * by that decimal expansion.  Parsing it once downward and once upward gives
 * an interval containing that exact rational.  Any inconclusive domain,
 * branch, or strict-inequality check terminates without a certificate.
 */

#include "certificate_interval.h"
#include "k3_orbit_data.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/* Appendix C.3 uses 500 bits; an override is useful for audit builds. */
#ifndef PRECISION
#define PRECISION 500
#endif
#define EXPONENT_MIN (-1073)
#define EXPONENT_MAX 1024

typedef cert_interval interval;

typedef struct {
    interval x;
    interval y;
    interval z;
} interval_point;

/* The directed-decimal bounds displayed in the Appendix C.3 table. */
static const char *const table_discriminant_lower[K3_ORBIT_LENGTH] = {
    "114", "8.70", "70.2", "8.59", "108",
    "10.4", "108", "8.59", "70.2", "8.70"
};

static const char *const table_first_residual_upper[K3_ORBIT_LENGTH] = {
    "5.67e-32", "3.17e-31", "1.29e-30", "1.09e-30", "1.67e-31",
    "1.90e-31", "1.76e-31", "7.92e-31", "2.00e-30", "5.73e-31"
};

static const char *const table_second_residual_upper[K3_ORBIT_LENGTH] = {
    "3.51e-31", "1.00e-30", "7.40e-31", "1.67e-31", "9.62e-32",
    "4.12e-31", "3.69e-31", "7.05e-30", "2.58e-30", "7.46e-31"
};

static const char *const table_norm_upper[K3_ORBIT_LENGTH] = {
    "3.56e-31", "1.05e-30", "1.49e-30", "1.10e-30", "1.93e-31",
    "4.53e-31", "4.09e-31", "7.09e-30", "3.26e-30", "9.40e-31"
};

static _Noreturn void fail(const char *message)
{
    fprintf(stderr, "CERTIFICATION FAILURE: %s\n", message);
    exit(EXIT_FAILURE);
}

static _Noreturn void failf(const char *format, ...)
{
    va_list arguments;

    fprintf(stderr, "CERTIFICATION FAILURE: ");
    va_start(arguments, format);
    vfprintf(stderr, format, arguments);
    va_end(arguments);
    fputc('\n', stderr);
    exit(EXIT_FAILURE);
}

/* Keep the shared kernel's status API behind certificate-specific failures. */
static void require_interval_status(
    cert_interval_status status,
    const char *context
)
{
    if (status != CERT_INTERVAL_OK) {
        failf("%s: %s", context, cert_interval_status_string(status));
    }
}

/* A malformed interval must never be allowed to make a comparison vacuous. */
static void interval_require_valid(const interval *x, const char *context)
{
    if (!cert_interval_valid(x)) {
        failf("%s: invalid interval", context);
    }
}

static void interval_init(interval *x)
{
    require_interval_status(
        cert_interval_init(x, PRECISION),
        "interval initialization"
    );
}

static void interval_clear(interval *x)
{
    cert_interval_clear(x);
}

static void point_init(interval_point *p)
{
    interval_init(&p->x);
    interval_init(&p->y);
    interval_init(&p->z);
}

static void point_clear(interval_point *p)
{
    interval_clear(&p->x);
    interval_clear(&p->y);
    interval_clear(&p->z);
}

static void interval_copy(interval *z, const interval *x)
{
    require_interval_status(
        cert_interval_copy(z, x),
        "interval copy"
    );
}

static void interval_set_si(interval *x, long value)
{
    require_interval_status(
        cert_interval_set_si(x, value),
        "integer conversion"
    );
}

static void interval_set_exact_decimal(interval *x, const char *decimal)
{
    require_interval_status(
        cert_interval_set_exact_decimal(x, decimal),
        "decimal conversion"
    );
}

static void interval_add(interval *z, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_add(z, x, y),
        "interval addition"
    );
}

static void interval_sub(interval *z, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_sub(z, x, y),
        "interval subtraction"
    );
}

static void interval_mul(interval *z, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_mul(z, x, y),
        "interval multiplication"
    );
}

static void interval_mul_si(interval *z, const interval *x, long value)
{
    require_interval_status(
        cert_interval_mul_si(z, x, value),
        "interval scalar multiplication"
    );
}

static void interval_div(interval *z, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_div(z, x, y),
        "interval division"
    );
}

static void interval_sqrt(interval *z, const interval *x)
{
    require_interval_status(
        cert_interval_sqrt(z, x),
        "interval square root"
    );
}

static void interval_one_plus_square(interval *z, const interval *x)
{
    interval square;
    interval one;

    interval_init(&square);
    interval_init(&one);
    interval_mul(&square, x, x);
    interval_set_si(&one, 1);
    interval_add(z, &one, &square);
    interval_clear(&square);
    interval_clear(&one);
}

/*
 * D(u,v) = 100u^2v^2 + 8(1+u^2)(1+v^2)
 *          - 4(1+u^2)^2(1+v^2)^2.
 */
static void discriminant(interval *result, const interval *u, const interval *v)
{
    interval u2;
    interval v2;
    interval one_plus_u2;
    interval one_plus_v2;
    interval term1;
    interval term2;
    interval term3;
    interval temporary1;
    interval temporary2;

    interval_init(&u2);
    interval_init(&v2);
    interval_init(&one_plus_u2);
    interval_init(&one_plus_v2);
    interval_init(&term1);
    interval_init(&term2);
    interval_init(&term3);
    interval_init(&temporary1);
    interval_init(&temporary2);

    interval_mul(&u2, u, u);
    interval_mul(&v2, v, v);
    interval_one_plus_square(&one_plus_u2, u);
    interval_one_plus_square(&one_plus_v2, v);

    interval_mul(&term1, &u2, &v2);
    interval_mul_si(&term1, &term1, 100);

    interval_mul(&term2, &one_plus_u2, &one_plus_v2);
    interval_mul_si(&term2, &term2, 8);

    interval_mul(&temporary1, &one_plus_u2, &one_plus_u2);
    interval_mul(&temporary2, &one_plus_v2, &one_plus_v2);
    interval_mul(&term3, &temporary1, &temporary2);
    interval_mul_si(&term3, &term3, 4);

    interval_add(&temporary1, &term1, &term2);
    interval_sub(result, &temporary1, &term3);

    interval_clear(&u2);
    interval_clear(&v2);
    interval_clear(&one_plus_u2);
    interval_clear(&one_plus_v2);
    interval_clear(&term1);
    interval_clear(&term2);
    interval_clear(&term3);
    interval_clear(&temporary1);
    interval_clear(&temporary2);
}

static void require_strictly_positive(
    const interval *x,
    const char *failure_message
)
{
    interval_require_valid(x, "positivity-test input");
    if (mpfr_sgn(x->lo) <= 0) {
        fail(failure_message);
    }
}

/* Decode periodicExact.c's six chart values into Psi_j and its branch. */
static int chart_number_from_code(int chart)
{
    int number;

    if (!k3_chart_decode(chart, &number, NULL)) {
        failf("invalid chart code %d", chart);
    }
    return number;
}

static int chart_sign_from_code(int chart)
{
    int sign;

    if (!k3_chart_decode(chart, NULL, &sign)) {
        failf("invalid chart code %d", chart);
    }
    return sign;
}

static void p_value(
    interval *result,
    interval *discriminant_interval,
    const interval *u,
    const interval *v,
    int sign
)
{
    interval square_root;
    interval numerator;
    interval denominator;
    interval uv;
    interval one_plus_u2;
    interval one_plus_v2;
    interval temporary;

    interval_init(&square_root);
    interval_init(&numerator);
    interval_init(&denominator);
    interval_init(&uv);
    interval_init(&one_plus_u2);
    interval_init(&one_plus_v2);
    interval_init(&temporary);

    /* p_+ and p_- are the two roots in the coordinate lost by projection. */
    discriminant(discriminant_interval, u, v);
    require_strictly_positive(
        discriminant_interval,
        "a chart discriminant was not certified strictly positive"
    );
    interval_sqrt(&square_root, discriminant_interval);

    interval_mul(&uv, u, v);
    interval_mul_si(&numerator, &uv, -10);
    if (sign > 0) {
        interval_add(&numerator, &numerator, &square_root);
    } else {
        interval_sub(&numerator, &numerator, &square_root);
    }

    interval_one_plus_square(&one_plus_u2, u);
    interval_one_plus_square(&one_plus_v2, v);
    interval_mul(&denominator, &one_plus_u2, &one_plus_v2);
    interval_mul_si(&denominator, &denominator, 2);
    interval_div(&temporary, &numerator, &denominator);
    interval_copy(result, &temporary);

    interval_clear(&square_root);
    interval_clear(&numerator);
    interval_clear(&denominator);
    interval_clear(&uv);
    interval_clear(&one_plus_u2);
    interval_clear(&one_plus_v2);
    interval_clear(&temporary);
}

static void chart_lift(
    interval_point *point,
    interval *discriminant_interval,
    const interval *u,
    const interval *v,
    int chart
)
{
    interval dependent_coordinate;
    int number = chart_number_from_code(chart);
    int sign = chart_sign_from_code(chart);

    /* Insert p_+ or p_- in the coordinate prescribed by Psi_j. */
    interval_init(&dependent_coordinate);
    p_value(
        &dependent_coordinate,
        discriminant_interval,
        u,
        v,
        sign
    );

    if (number == 1) {
        interval_copy(&point->x, &dependent_coordinate);
        interval_copy(&point->y, u);
        interval_copy(&point->z, v);
    } else if (number == 2) {
        interval_copy(&point->x, u);
        interval_copy(&point->y, &dependent_coordinate);
        interval_copy(&point->z, v);
    } else if (number == 3) {
        interval_copy(&point->x, u);
        interval_copy(&point->y, v);
        interval_copy(&point->z, &dependent_coordinate);
    } else {
        fail("invalid decoded chart number");
    }

    interval_clear(&dependent_coordinate);
}

static void alpha(interval *result, const interval *u, const interval *v)
{
    interval numerator;
    interval denominator;
    interval one_plus_u2;
    interval one_plus_v2;

    interval_init(&numerator);
    interval_init(&denominator);
    interval_init(&one_plus_u2);
    interval_init(&one_plus_v2);

    interval_mul(&numerator, u, v);
    interval_one_plus_square(&one_plus_u2, u);
    interval_one_plus_square(&one_plus_v2, v);
    interval_mul(&denominator, &one_plus_u2, &one_plus_v2);
    interval_div(result, &numerator, &denominator);

    interval_clear(&numerator);
    interval_clear(&denominator);
    interval_clear(&one_plus_u2);
    interval_clear(&one_plus_v2);
}

/* The three coordinate involutions, applied in the manuscript's order. */
static void sigma_1(interval_point *point)
{
    interval alpha_yz;
    interval correction;

    interval_init(&alpha_yz);
    interval_init(&correction);
    alpha(&alpha_yz, &point->y, &point->z);
    interval_mul_si(&correction, &alpha_yz, 10);
    interval_mul_si(&point->x, &point->x, -1);
    interval_sub(&point->x, &point->x, &correction);
    interval_clear(&alpha_yz);
    interval_clear(&correction);
}

static void sigma_2(interval_point *point)
{
    interval alpha_xz;
    interval correction;

    interval_init(&alpha_xz);
    interval_init(&correction);
    alpha(&alpha_xz, &point->x, &point->z);
    interval_mul_si(&correction, &alpha_xz, 10);
    interval_mul_si(&point->y, &point->y, -1);
    interval_sub(&point->y, &point->y, &correction);
    interval_clear(&alpha_xz);
    interval_clear(&correction);
}

static void sigma_3(interval_point *point)
{
    interval alpha_xy;
    interval correction;

    interval_init(&alpha_xy);
    interval_init(&correction);
    alpha(&alpha_xy, &point->x, &point->y);
    interval_mul_si(&correction, &alpha_xy, 10);
    interval_mul_si(&point->z, &point->z, -1);
    interval_sub(&point->z, &point->z, &correction);
    interval_clear(&alpha_xy);
    interval_clear(&correction);
}

static void apply_f(interval_point *point)
{
    /* f = sigma_3 o sigma_2 o sigma_1. */
    sigma_1(point);
    sigma_2(point);
    sigma_3(point);
}

static void chart_projection(
    interval *u,
    interval *v,
    const interval_point *point,
    int chart
)
{
    int number = chart_number_from_code(chart);

    /* Remove the dependent coordinate to recover the unshifted chart pair. */
    if (number == 1) {
        interval_copy(u, &point->y);
        interval_copy(v, &point->z);
    } else if (number == 2) {
        interval_copy(u, &point->x);
        interval_copy(v, &point->z);
    } else if (number == 3) {
        interval_copy(u, &point->x);
        interval_copy(v, &point->y);
    } else {
        fail("invalid decoded chart number");
    }
}

static void verify_output_chart_branch(
    interval *output_discriminant,
    interval *branch_indicator,
    const interval_point *point,
    int chart
)
{
    interval u;
    interval v;
    interval dependent;
    interval one_plus_u2;
    interval one_plus_v2;
    interval denominator_factor;
    interval temporary;
    interval uv;
    int number = chart_number_from_code(chart);
    int sign = chart_sign_from_code(chart);

    interval_init(&u);
    interval_init(&v);
    interval_init(&dependent);
    interval_init(&one_plus_u2);
    interval_init(&one_plus_v2);
    interval_init(&denominator_factor);
    interval_init(&temporary);
    interval_init(&uv);

    /* First recover (u,v) and the coordinate w lost by projection. */
    chart_projection(&u, &v, point, chart);
    if (number == 1) {
        interval_copy(&dependent, &point->x);
    } else if (number == 2) {
        interval_copy(&dependent, &point->y);
    } else {
        interval_copy(&dependent, &point->z);
    }

    discriminant(output_discriminant, &u, &v);
    require_strictly_positive(
        output_discriminant,
        "an output-chart discriminant was not certified strictly positive"
    );

    /*
     * For a root w=p^+ or p^-,
     *
     *   2(1+u^2)(1+v^2)w + 10uv = +/- sqrt(D(u,v)).
     *
     * Its sign therefore certifies that the point is on the required
     * branch of the next chart.
     */
    interval_one_plus_square(&one_plus_u2, &u);
    interval_one_plus_square(&one_plus_v2, &v);
    interval_mul(
        &denominator_factor,
        &one_plus_u2,
        &one_plus_v2
    );
    interval_mul_si(&denominator_factor, &denominator_factor, 2);
    interval_mul(branch_indicator, &denominator_factor, &dependent);
    interval_mul(&uv, &u, &v);
    interval_mul_si(&temporary, &uv, 10);
    interval_add(branch_indicator, branch_indicator, &temporary);

    if (sign > 0 && mpfr_sgn(branch_indicator->lo) <= 0) {
        fail("the output point was not certified to lie on the plus branch");
    }
    if (sign < 0 && mpfr_sgn(branch_indicator->hi) >= 0) {
        fail("the output point was not certified to lie on the minus branch");
    }

    interval_clear(&u);
    interval_clear(&v);
    interval_clear(&dependent);
    interval_clear(&one_plus_u2);
    interval_clear(&one_plus_v2);
    interval_clear(&denominator_factor);
    interval_clear(&temporary);
    interval_clear(&uv);
}

/*
 * Rigorous analogue of periodicExact.c's f_c.
 *
 * On entry, s and t enclose the two input-chart coordinates.  On return,
 * they enclose the two unshifted coordinates in the output chart.  The
 * additional output intervals retain the two checks proving that the image
 * actually belongs to that output chart: D>0 and the correct branch sign.
 */
static void f_c(
    interval *s,
    interval *t,
    interval *input_discriminant,
    interval *output_discriminant,
    interval *branch_indicator,
    int input_chart,
    int output_chart
)
{
    interval_point point;

    point_init(&point);
    chart_lift(&point, input_discriminant, s, t, input_chart);
    apply_f(&point);
    verify_output_chart_branch(
        output_discriminant,
        branch_indicator,
        &point,
        output_chart
    );
    chart_projection(s, t, &point, output_chart);
    point_clear(&point);
}

/* m(R)=max(|lower(R)|,|upper(R)|), rounded upward. */
static void interval_abs_upper(mpfr_t result, const interval *x)
{
    require_interval_status(
        cert_interval_abs_upper(result, x),
        "absolute-value bound"
    );
}

/*
 * Rigorous analogue of periodicExact.c's dist_2plane.  The coordinate
 * differences are retained because Appendix C.3 records separate bounds
 * for them as well as the Euclidean norm bound.
 */
static void dist_2plane(
    mpfr_t result,
    interval *first_residual,
    interval *second_residual,
    const interval *first,
    const interval *second,
    const interval *next_first,
    const interval *next_second
)
{
    mpfr_t first_absolute;
    mpfr_t second_absolute;
    mpfr_t first_square;
    mpfr_t second_square;
    mpfr_t sum;

    mpfr_init2(first_absolute, PRECISION);
    mpfr_init2(second_absolute, PRECISION);
    mpfr_init2(first_square, PRECISION);
    mpfr_init2(second_square, PRECISION);
    mpfr_init2(sum, PRECISION);

    interval_sub(first_residual, first, next_first);
    interval_sub(second_residual, second, next_second);
    interval_abs_upper(first_absolute, first_residual);
    interval_abs_upper(second_absolute, second_residual);

    mpfr_mul(
        first_square,
        first_absolute,
        first_absolute,
        MPFR_RNDU
    );
    mpfr_mul(
        second_square,
        second_absolute,
        second_absolute,
        MPFR_RNDU
    );
    mpfr_add(sum, first_square, second_square, MPFR_RNDU);
    mpfr_sqrt(result, sum, MPFR_RNDU);
    if (!mpfr_number_p(result)) {
        fail("residual norm upper bound is not finite");
    }

    mpfr_clear(first_absolute);
    mpfr_clear(second_absolute);
    mpfr_clear(first_square);
    mpfr_clear(second_square);
    mpfr_clear(sum);
}

/* Prove that the manuscript's displayed D lower bound is conservative. */
static void certify_table_lower_bound(
    const interval *computed,
    const char *displayed_decimal,
    int row
)
{
    interval displayed;

    interval_init(&displayed);
    interval_set_exact_decimal(&displayed, displayed_decimal);
    if (mpfr_cmp(computed->lo, displayed.hi) < 0) {
        interval_clear(&displayed);
        failf(
            "Appendix C.3 row %d has an invalid discriminant lower bound",
            row
        );
    }
    interval_clear(&displayed);
}

/* Prove that a manuscript decimal is an upper bound for a computed value. */
static void certify_table_upper_bound(
    const mpfr_t computed,
    const char *displayed_decimal,
    int row,
    const char *column
)
{
    interval displayed;

    interval_init(&displayed);
    interval_set_exact_decimal(&displayed, displayed_decimal);
    if (!mpfr_number_p(computed) || mpfr_cmp(computed, displayed.lo) > 0) {
        interval_clear(&displayed);
        failf(
            "Appendix C.3 row %d has an invalid %s upper bound",
            row,
            column
        );
    }
    interval_clear(&displayed);
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

static const char *chart_label(int chart)
{
    const char *label = k3_chart_label(chart);

    if (label == NULL) {
        failf("invalid chart code %d", chart);
    }
    return label;
}

int main(void)
{
    interval a[K3_ORBIT_LENGTH];
    interval b[K3_ORBIT_LENGTH];
    interval delta;
    interval manuscript_maximum_norm;
    mpfr_t maximum_norm_upper;

    /* Every operation below supplies its own directed rounding mode. */
    mpfr_set_default_rounding_mode(MPFR_RNDN);

    if (!k3_orbit_data_are_valid()) {
        fail("shared pseudo-orbit/chart data failed its self-check");
    }

    /* Match the precision and exponent range used by periodicExact.c. */
    mpfr_set_default_prec(PRECISION);
    if (
        mpfr_set_emin(EXPONENT_MIN) != 0
        || mpfr_set_emax(EXPONENT_MAX) != 0
    ) {
        fail("MPFR rejected the requested exponent range");
    }
    mpfr_clear_flags();

    interval_init(&delta);
    interval_init(&manuscript_maximum_norm);
    interval_set_exact_decimal(&delta, "1e-29");
    interval_set_exact_decimal(&manuscript_maximum_norm, "7.09e-30");
    mpfr_init2(maximum_norm_upper, PRECISION);
    mpfr_set_zero(maximum_norm_upper, 1);

    for (int i = 0; i < K3_ORBIT_LENGTH; i++) {
        interval_init(&a[i]);
        interval_init(&b[i]);
        interval_set_exact_decimal(&a[i], k3_orbit_a_decimal[i]);
        interval_set_exact_decimal(&b[i], k3_orbit_b_decimal[i]);
    }

    printf("MPFR version: %s\n", mpfr_get_version());
    printf("precision: %d bits\n", PRECISION);
    printf("exponent range: [%d,%d]\n", EXPONENT_MIN, EXPONENT_MAX);
    printf(
        "Each listed decimal is interpreted as the exact finite decimal.\n"
    );
    printf(
        "Every table entry below is checked against the full-precision enclosure.\n\n"
    );
    printf("Appendix C.3 residual-certificate table\n");
    printf(" i  chart   D_i^lb       M_i,1       M_i,2       N_i\n");

    for (int i = 0; i < K3_ORBIT_LENGTH; i++) {
        int next = (i + 1) % K3_ORBIT_LENGTH;
        interval image_first;
        interval image_second;
        interval input_discriminant;
        interval output_discriminant;
        interval branch_indicator;
        interval first_residual;
        interval second_residual;
        mpfr_t first_residual_upper;
        mpfr_t second_residual_upper;
        mpfr_t norm_upper;

        interval_init(&image_first);
        interval_init(&image_second);
        interval_init(&input_discriminant);
        interval_init(&output_discriminant);
        interval_init(&branch_indicator);
        interval_init(&first_residual);
        interval_init(&second_residual);
        mpfr_init2(first_residual_upper, PRECISION);
        mpfr_init2(second_residual_upper, PRECISION);
        mpfr_init2(norm_upper, PRECISION);

        /* As in periodicExact.c, start from copies of (a_i,b_i). */
        interval_copy(&image_first, &a[i]);
        interval_copy(&image_second, &b[i]);

        /*
         * Enclose the absolute output-chart coordinates.  Besides applying
         * f, f_c proves D>0 and the required branch of the next chart; these
         * are precisely the checks establishing 0 in U_i.
         */
        f_c(
            &image_first,
            &image_second,
            &input_discriminant,
            &output_discriminant,
            &branch_indicator,
            k3_orbit_chart_code[i],
            k3_orbit_chart_code[next]
        );

        /* Subtract the next chart center to enclose f_i^c(0). */
        dist_2plane(
            norm_upper,
            &first_residual,
            &second_residual,
            &image_first,
            &image_second,
            &a[next],
            &b[next]
        );
        interval_abs_upper(first_residual_upper, &first_residual);
        interval_abs_upper(second_residual_upper, &second_residual);

        /* Certify, rather than merely reproduce, this manuscript row. */
        certify_table_lower_bound(
            &input_discriminant,
            table_discriminant_lower[i],
            i
        );
        certify_table_upper_bound(
            first_residual_upper,
            table_first_residual_upper[i],
            i,
            "M_i,1"
        );
        certify_table_upper_bound(
            second_residual_upper,
            table_second_residual_upper[i],
            i,
            "M_i,2"
        );
        certify_table_upper_bound(
            norm_upper,
            table_norm_upper[i],
            i,
            "N_i"
        );

        if (mpfr_cmp(norm_upper, maximum_norm_upper) > 0) {
            mpfr_set(maximum_norm_upper, norm_upper, MPFR_RNDU);
        }

        /* Print the exact decimals whose directed inequalities were checked. */
        printf(
            "%2d  %-5s  %-12s %-12s %-12s %s\n",
            i,
            chart_label(k3_orbit_chart_code[i]),
            table_discriminant_lower[i],
            table_first_residual_upper[i],
            table_second_residual_upper[i],
            table_norm_upper[i]
        );

        /* Full enclosures remain visible for independent audit/debugging. */
        mpfr_printf(
            "    raw R_i,1 = [%+.36RDe, %+.36RUe]\n"
            "    raw R_i,2 = [%+.36RDe, %+.36RUe]\n"
            "    raw D_input lower = %.36RDe; raw D_output lower = %.36RDe\n"
            "    raw branch indicator = [%+.36RDe, %+.36RUe]\n"
            "    raw N_i = %.36RUe\n",
            first_residual.lo,
            first_residual.hi,
            second_residual.lo,
            second_residual.hi,
            input_discriminant.lo,
            output_discriminant.lo,
            branch_indicator.lo,
            branch_indicator.hi,
            norm_upper
        );

        interval_clear(&image_first);
        interval_clear(&image_second);
        interval_clear(&input_discriminant);
        interval_clear(&output_discriminant);
        interval_clear(&branch_indicator);
        interval_clear(&first_residual);
        interval_clear(&second_residual);
        mpfr_clear(first_residual_upper);
        mpfr_clear(second_residual_upper);
        mpfr_clear(norm_upper);
    }

    /* Reject every exceptional MPFR condition before trusting comparisons. */
    verify_mpfr_status();

    if (mpfr_cmp(maximum_norm_upper, manuscript_maximum_norm.lo) > 0) {
        fail("max_i N_i <= 7.09e-30 was not certified");
    }
    if (mpfr_cmp(manuscript_maximum_norm.hi, delta.lo) >= 0) {
        fail("7.09e-30 < 1e-29 was not certified");
    }

    printf("\n");
    mpfr_printf(
        "maximum rigorous norm upper bound: %.40RUe\n",
        maximum_norm_upper
    );
    mpfr_printf("downward-rounded delta endpoint: %.40RDe\n", delta.lo);
    printf(
        "\nCertified conclusions\n"
        "PASS: every input and output discriminant is strictly positive.\n"
        "PASS: every image is on the required branch, so 0 is in every U_i.\n"
        "PASS: every directed bound displayed in Appendix C.3 is certified.\n"
        "PASS: max_i ||f_i^c(0)||_2 <= 7.09e-30 < 1e-29.\n"
    );

    for (int i = 0; i < K3_ORBIT_LENGTH; i++) {
        interval_clear(&a[i]);
        interval_clear(&b[i]);
    }
    interval_clear(&delta);
    interval_clear(&manuscript_maximum_norm);
    mpfr_clear(maximum_norm_upper);
    mpfr_free_cache();
    return EXIT_SUCCESS;
}
