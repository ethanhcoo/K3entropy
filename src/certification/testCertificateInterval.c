/* Focused regression tests for certificate_interval.c. */

#include "certificate_interval.h"

#include <stdio.h>
#include <stdlib.h>

#define EXPECT(condition, message)                                           \
    do {                                                                     \
        if (!(condition)) {                                                  \
            fprintf(stderr, "testCertificateInterval: %s\n", (message));    \
            return 0;                                                        \
        }                                                                    \
    } while (0)

#define EXPECT_OK(expression)                                                \
    do {                                                                     \
        cert_interval_status test_status_ = (expression);                    \
        if (test_status_ != CERT_INTERVAL_OK) {                              \
            fprintf(                                                        \
                stderr,                                                      \
                "testCertificateInterval: %s returned %s\n",                \
                #expression,                                                 \
                cert_interval_status_string(test_status_)                    \
            );                                                               \
            return 0;                                                        \
        }                                                                    \
    } while (0)

static void set_bounds_si(cert_interval *x, long lo, long hi)
{
    mpfr_set_si(x->lo, lo, MPFR_RNDD);
    mpfr_set_si(x->hi, hi, MPFR_RNDU);
}

static int has_bounds_si(const cert_interval *x, long lo, long hi)
{
    return cert_interval_valid(x)
        && mpfr_cmp_si(x->lo, lo) == 0
        && mpfr_cmp_si(x->hi, hi) == 0;
}

static int has_bounds_dyadic(
    const cert_interval *x,
    long lo,
    mpfr_exp_t lo_exponent,
    long hi,
    mpfr_exp_t hi_exponent
)
{
    return cert_interval_valid(x)
        && mpfr_cmp_si_2exp(x->lo, lo, lo_exponent) == 0
        && mpfr_cmp_si_2exp(x->hi, hi, hi_exponent) == 0;
}

static int encloses_rational(const cert_interval *x, long numerator, long denominator)
{
    mpq_t rational;
    int result;

    mpq_init(rational);
    mpq_set_si(rational, numerator, denominator);
    mpq_canonicalize(rational);
    result = cert_interval_valid(x)
        && mpfr_cmp_q(x->lo, rational) <= 0
        && mpfr_cmp_q(x->hi, rational) >= 0;
    mpq_clear(rational);
    return result;
}

static int test_decimal_and_validity(mpfr_prec_t precision)
{
    cert_interval x;
    cert_interval invalid;
    mpfr_t old_lo;
    mpfr_t old_hi;

    EXPECT_OK(cert_interval_init(&x, precision));
    EXPECT_OK(cert_interval_init(&invalid, precision));
    EXPECT(!cert_interval_valid(&x), "freshly initialized NaN interval was valid");
    EXPECT(cert_interval_contains_zero(&x), "invalid interval was not conservative");

    EXPECT_OK(cert_interval_set_exact_decimal(&x, "0.1"));
    EXPECT(cert_interval_valid(&x), "parsed decimal interval was invalid");
    EXPECT(encloses_rational(&x, 1, 10), "0.1 did not enclose exact 1/10");
    EXPECT(mpfr_get_prec(x.lo) == precision, "lower precision changed");
    EXPECT(mpfr_get_prec(x.hi) == precision, "upper precision changed");

    mpfr_init2(old_lo, precision);
    mpfr_init2(old_hi, precision);
    mpfr_set(old_lo, x.lo, MPFR_RNDN);
    mpfr_set(old_hi, x.hi, MPFR_RNDN);
    EXPECT(
        cert_interval_set_exact_decimal(&x, "not-a-decimal")
            == CERT_INTERVAL_PARSE_ERROR,
        "malformed decimal had the wrong status"
    );
    EXPECT(
        mpfr_cmp(x.lo, old_lo) == 0 && mpfr_cmp(x.hi, old_hi) == 0,
        "failed decimal parse modified its destination"
    );
    EXPECT(
        cert_interval_set_exact_decimal(&x, "nan") == CERT_INTERVAL_PARSE_ERROR,
        "NaN decimal was accepted"
    );
    EXPECT(
        cert_interval_set_exact_decimal(&x, NULL)
            == CERT_INTERVAL_INVALID_ARGUMENT,
        "NULL decimal had the wrong status"
    );
    EXPECT(
        cert_interval_copy(&x, &invalid) == CERT_INTERVAL_INVALID_INTERVAL,
        "copy accepted an invalid source"
    );

    mpfr_clear(old_lo);
    mpfr_clear(old_hi);
    cert_interval_clear(&invalid);
    cert_interval_clear(&x);
    return 1;
}

static int test_additive_and_aliasing(mpfr_prec_t precision)
{
    cert_interval x;
    cert_interval y;
    cert_interval z;

    EXPECT_OK(cert_interval_init(&x, precision));
    EXPECT_OK(cert_interval_init(&y, precision));
    EXPECT_OK(cert_interval_init(&z, precision));

    set_bounds_si(&x, 2, 3);
    set_bounds_si(&y, 4, 5);
    EXPECT_OK(cert_interval_add(&z, &x, &y));
    EXPECT(has_bounds_si(&z, 6, 8), "addition bounds were wrong");
    EXPECT_OK(cert_interval_add(&x, &x, &y));
    EXPECT(has_bounds_si(&x, 6, 8), "left-aliased addition was wrong");

    set_bounds_si(&x, 2, 3);
    set_bounds_si(&y, 4, 5);
    EXPECT_OK(cert_interval_sub(&y, &x, &y));
    EXPECT(has_bounds_si(&y, -3, -1), "right-aliased subtraction was wrong");
    EXPECT_OK(cert_interval_neg(&y, &y));
    EXPECT(has_bounds_si(&y, 1, 3), "aliased negation was wrong");
    EXPECT_OK(cert_interval_add_si(&y, &y, -2));
    EXPECT(has_bounds_si(&y, -1, 1), "aliased integer addition was wrong");
    EXPECT(cert_interval_contains_zero(&y), "zero containment missed zero");
    EXPECT_OK(cert_interval_ui_sub(&y, 4UL, &y));
    EXPECT(has_bounds_si(&y, 3, 5), "aliased unsigned subtraction was wrong");
    EXPECT(!cert_interval_contains_zero(&y), "zero containment found zero");

    EXPECT_OK(cert_interval_set_si(&z, -17));
    EXPECT(has_bounds_si(&z, -17, -17), "signed integer conversion was wrong");
    EXPECT_OK(cert_interval_copy(&z, &y));
    EXPECT(has_bounds_si(&z, 3, 5), "interval copy was wrong");
    EXPECT_OK(cert_interval_copy(&z, &z));
    EXPECT(has_bounds_si(&z, 3, 5), "self-copy was wrong");

    cert_interval_clear(&z);
    cert_interval_clear(&y);
    cert_interval_clear(&x);
    return 1;
}

static int check_product(
    cert_interval *x,
    cert_interval *y,
    cert_interval *z,
    long x_lo,
    long x_hi,
    long y_lo,
    long y_hi,
    long expected_lo,
    long expected_hi
)
{
    set_bounds_si(x, x_lo, x_hi);
    set_bounds_si(y, y_lo, y_hi);
    if (cert_interval_mul(z, x, y) != CERT_INTERVAL_OK) {
        return 0;
    }
    return has_bounds_si(z, expected_lo, expected_hi);
}

static int test_multiplication_quadrants(mpfr_prec_t precision)
{
    cert_interval x;
    cert_interval y;
    cert_interval z;

    EXPECT_OK(cert_interval_init(&x, precision));
    EXPECT_OK(cert_interval_init(&y, precision));
    EXPECT_OK(cert_interval_init(&z, precision));

    EXPECT(
        check_product(&x, &y, &z, 2, 3, 4, 5, 8, 15),
        "positive-positive product was wrong"
    );
    EXPECT(
        check_product(&x, &y, &z, -3, -2, 4, 5, -15, -8),
        "negative-positive product was wrong"
    );
    EXPECT(
        check_product(&x, &y, &z, 2, 3, -5, -4, -15, -8),
        "positive-negative product was wrong"
    );
    EXPECT(
        check_product(&x, &y, &z, -3, -2, -5, -4, 8, 15),
        "negative-negative product was wrong"
    );
    EXPECT(
        check_product(&x, &y, &z, -2, 3, -4, 5, -12, 15),
        "mixed-sign product was wrong"
    );

    set_bounds_si(&x, 2, 3);
    set_bounds_si(&y, 4, 5);
    EXPECT_OK(cert_interval_mul(&x, &x, &y));
    EXPECT(has_bounds_si(&x, 8, 15), "left-aliased product was wrong");
    set_bounds_si(&x, -3, -2);
    set_bounds_si(&y, -5, -4);
    EXPECT_OK(cert_interval_mul(&y, &x, &y));
    EXPECT(has_bounds_si(&y, 8, 15), "right-aliased product was wrong");

    set_bounds_si(&x, -2, 3);
    EXPECT_OK(cert_interval_mul_si(&x, &x, -3));
    EXPECT(has_bounds_si(&x, -9, 6), "aliased negative scalar product was wrong");
    EXPECT_OK(cert_interval_mul_si(&x, &x, 0));
    EXPECT(has_bounds_si(&x, 0, 0), "zero scalar product was wrong");

    set_bounds_si(&x, 2, 3);
    EXPECT_OK(cert_interval_square(&x, &x));
    EXPECT(has_bounds_si(&x, 4, 9), "positive aliased square was wrong");
    set_bounds_si(&x, -3, -2);
    EXPECT_OK(cert_interval_square(&z, &x));
    EXPECT(has_bounds_si(&z, 4, 9), "negative square was wrong");
    set_bounds_si(&x, -2, 3);
    EXPECT_OK(cert_interval_square(&z, &x));
    EXPECT(has_bounds_si(&z, 0, 9), "mixed-sign square was wrong");

    cert_interval_clear(&z);
    cert_interval_clear(&y);
    cert_interval_clear(&x);
    return 1;
}

static int test_division_and_domains(mpfr_prec_t precision)
{
    cert_interval x;
    cert_interval y;
    cert_interval z;

    EXPECT_OK(cert_interval_init(&x, precision));
    EXPECT_OK(cert_interval_init(&y, precision));
    EXPECT_OK(cert_interval_init(&z, precision));

    set_bounds_si(&x, 2, 4);
    EXPECT_OK(cert_interval_inv(&x, &x));
    EXPECT(
        has_bounds_dyadic(&x, 1, -2, 1, -1),
        "positive aliased reciprocal was wrong"
    );
    set_bounds_si(&x, -4, -2);
    EXPECT_OK(cert_interval_inv(&z, &x));
    EXPECT(
        has_bounds_dyadic(&z, -1, -1, -1, -2),
        "negative reciprocal was wrong"
    );

    set_bounds_si(&x, 2, 4);
    set_bounds_si(&y, 2, 4);
    EXPECT_OK(cert_interval_div(&x, &x, &y));
    EXPECT(
        has_bounds_dyadic(&x, 1, -1, 2, 0),
        "left-aliased interval division was wrong"
    );
    set_bounds_si(&x, 2, 4);
    set_bounds_si(&y, 2, 4);
    EXPECT_OK(cert_interval_div(&y, &x, &y));
    EXPECT(
        has_bounds_dyadic(&y, 1, -1, 2, 0),
        "right-aliased interval division was wrong"
    );
    set_bounds_si(&x, -4, 8);
    EXPECT_OK(cert_interval_div_ui(&x, &x, 4UL));
    EXPECT(has_bounds_si(&x, -1, 2), "aliased unsigned division was wrong");

    EXPECT_OK(cert_interval_set_si(&z, 42));
    set_bounds_si(&x, -1, 1);
    EXPECT(
        cert_interval_inv(&z, &x) == CERT_INTERVAL_DIVISION_BY_ZERO,
        "zero-containing reciprocal had the wrong status"
    );
    EXPECT(has_bounds_si(&z, 42, 42), "failed reciprocal modified destination");
    EXPECT(
        cert_interval_div(&z, &z, &x) == CERT_INTERVAL_DIVISION_BY_ZERO,
        "zero-containing division had the wrong status"
    );
    EXPECT(
        cert_interval_div_ui(&z, &z, 0UL) == CERT_INTERVAL_DIVISION_BY_ZERO,
        "unsigned division by zero had the wrong status"
    );

    set_bounds_si(&x, 4, 9);
    EXPECT_OK(cert_interval_sqrt(&x, &x));
    EXPECT(has_bounds_si(&x, 2, 3), "aliased square root was wrong");
    set_bounds_si(&x, -1, 4);
    EXPECT_OK(cert_interval_set_si(&z, 42));
    EXPECT(
        cert_interval_sqrt(&z, &x) == CERT_INTERVAL_NEGATIVE_RADICAND,
        "partly negative square root had the wrong status"
    );
    EXPECT(has_bounds_si(&z, 42, 42), "failed square root modified destination");
    set_bounds_si(&x, -4, -1);
    EXPECT(
        cert_interval_sqrt(&z, &x) == CERT_INTERVAL_NEGATIVE_RADICAND,
        "negative square root had the wrong status"
    );

    cert_interval_clear(&z);
    cert_interval_clear(&y);
    cert_interval_clear(&x);
    return 1;
}

static int test_power_expand_and_abs(mpfr_prec_t precision)
{
    cert_interval x;
    cert_interval radius;
    cert_interval z;
    mpfr_t bound;

    EXPECT_OK(cert_interval_init(&x, precision));
    EXPECT_OK(cert_interval_init(&radius, precision));
    EXPECT_OK(cert_interval_init(&z, precision));
    mpfr_init2(bound, precision);

    set_bounds_si(&x, -2, 3);
    EXPECT_OK(cert_interval_pow_ui(&z, &x, 0UL));
    EXPECT(has_bounds_si(&z, 1, 1), "zeroth power was wrong");
    EXPECT_OK(cert_interval_pow_ui(&z, &x, 2UL));
    EXPECT(has_bounds_si(&z, 0, 9), "mixed-sign even power was wrong");
    EXPECT_OK(cert_interval_pow_ui(&x, &x, 3UL));
    EXPECT(has_bounds_si(&x, -8, 27), "aliased odd power was wrong");
    set_bounds_si(&x, -3, -2);
    EXPECT_OK(cert_interval_pow_ui(&z, &x, 4UL));
    EXPECT(has_bounds_si(&z, 16, 81), "negative even power was wrong");

    set_bounds_si(&x, 2, 3);
    EXPECT_OK(cert_interval_set_exact_decimal(&radius, "0.5"));
    EXPECT_OK(cert_interval_expand(&x, &x, &radius));
    EXPECT(
        has_bounds_dyadic(&x, 3, -1, 7, -1),
        "aliased interval expansion was wrong"
    );
    set_bounds_si(&radius, -1, 1);
    EXPECT_OK(cert_interval_set_si(&z, 42));
    EXPECT(
        cert_interval_expand(&z, &x, &radius)
            == CERT_INTERVAL_NEGATIVE_RADIUS,
        "negative radius had the wrong status"
    );
    EXPECT(has_bounds_si(&z, 42, 42), "failed expansion modified destination");

    set_bounds_si(&x, -3, 2);
    EXPECT_OK(cert_interval_abs_upper(bound, &x));
    EXPECT(mpfr_cmp_si(bound, 3) == 0, "absolute-value upper bound was wrong");
    EXPECT_OK(cert_interval_abs_upper(x.lo, &x));
    EXPECT(mpfr_cmp_si(x.lo, 3) == 0, "aliased absolute upper bound was wrong");

    EXPECT_OK(cert_interval_set_exact_decimal(&x, "0.1"));
    EXPECT_OK(cert_interval_square(&z, &x));
    EXPECT(encloses_rational(&z, 1, 100), "square of 0.1 missed exact 1/100");

    mpfr_clear(bound);
    cert_interval_clear(&z);
    cert_interval_clear(&radius);
    cert_interval_clear(&x);
    return 1;
}

/* A failed precision-narrowing operation must not partly overwrite output. */
static int test_unchanged_on_range_failure(void)
{
    const mpfr_exp_t old_emax = mpfr_get_emax();
    cert_interval source;
    cert_interval destination;
    mpfr_t bound;
    int success = 1;

    if (mpfr_set_emax(10) != 0) {
        fprintf(stderr, "testCertificateInterval: could not set test emax\n");
        return 0;
    }
    mpfr_clear_flags();
    if (cert_interval_init(&source, 100) != CERT_INTERVAL_OK ||
        cert_interval_init(&destination, 53) != CERT_INTERVAL_OK) {
        fprintf(stderr, "testCertificateInterval: range-test init failed\n");
        (void)mpfr_set_emax(old_emax);
        return 0;
    }
    mpfr_init2(bound, 53);

    /* The 100-bit largest finite value at emax=10 overflows upward at 53 bits. */
    mpfr_set_ui(source.lo, 1, MPFR_RNDN);
    mpfr_nextbelow(source.lo);
    mpfr_mul_2si(source.lo, source.lo, 10, MPFR_RNDN);
    mpfr_set(source.hi, source.lo, MPFR_RNDN);

    if (cert_interval_set_si(&destination, 7) != CERT_INTERVAL_OK ||
        cert_interval_copy(&destination, &source)
            != CERT_INTERVAL_INVALID_INTERVAL ||
        mpfr_cmp_si(destination.lo, 7) != 0 ||
        mpfr_cmp_si(destination.hi, 7) != 0) {
        fprintf(
            stderr,
            "testCertificateInterval: failed copy modified its destination\n"
        );
        success = 0;
    }

    mpfr_set_si(bound, 7, MPFR_RNDN);
    if (cert_interval_abs_upper(bound, &source)
            != CERT_INTERVAL_INVALID_INTERVAL ||
        mpfr_cmp_si(bound, 7) != 0) {
        fprintf(
            stderr,
            "testCertificateInterval: failed abs bound modified its result\n"
        );
        success = 0;
    }

    mpfr_clear(bound);
    cert_interval_clear(&destination);
    cert_interval_clear(&source);
    if (mpfr_set_emax(old_emax) != 0) {
        fprintf(stderr, "testCertificateInterval: could not restore emax\n");
        success = 0;
    }
    mpfr_clear_flags();
    return success;
}

static int run_precision(mpfr_prec_t precision)
{
    if (!test_decimal_and_validity(precision)) {
        return 0;
    }
    if (!test_additive_and_aliasing(precision)) {
        return 0;
    }
    if (!test_multiplication_quadrants(precision)) {
        return 0;
    }
    if (!test_division_and_domains(precision)) {
        return 0;
    }
    if (!test_power_expand_and_abs(precision)) {
        return 0;
    }
    printf("  precision %ld bits: PASS\n", (long)precision);
    return 1;
}

int main(void)
{
    const mpfr_prec_t precisions[] = {53, 500, 5000};
    const size_t precision_count = sizeof(precisions) / sizeof(precisions[0]);
    cert_interval uninitialized_storage;
    size_t i;

    EXPECT(
        cert_interval_init(NULL, 53) == CERT_INTERVAL_INVALID_ARGUMENT,
        "NULL initialization had the wrong status"
    );
    EXPECT(
        cert_interval_init(&uninitialized_storage, 0)
            == CERT_INTERVAL_INVALID_PRECISION,
        "invalid precision had the wrong status"
    );
    EXPECT(
        cert_interval_status_string(CERT_INTERVAL_OK) != NULL,
        "status string was NULL"
    );

    puts("Testing outward-rounded shared interval arithmetic:");
    for (i = 0; i < precision_count; ++i) {
        if (!run_precision(precisions[i])) {
            return EXIT_FAILURE;
        }
    }
    if (!test_unchanged_on_range_failure()) {
        return EXIT_FAILURE;
    }
    puts("certificate_interval tests: PASS");
    return EXIT_SUCCESS;
}
