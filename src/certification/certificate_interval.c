#include "certificate_interval.h"

#include <stddef.h>

/*
 * All arithmetic is performed in temporaries having the destination's
 * precision.  Besides keeping each endpoint rounded in the conservative
 * direction, this makes every operation safe when result aliases an input.
 */

static cert_interval_status destination_precision(
    const cert_interval *result,
    mpfr_prec_t *precision
)
{
    mpfr_prec_t lo_precision;
    mpfr_prec_t hi_precision;

    if (result == NULL || precision == NULL) {
        return CERT_INTERVAL_INVALID_ARGUMENT;
    }
    lo_precision = mpfr_get_prec(result->lo);
    hi_precision = mpfr_get_prec(result->hi);
    if (lo_precision != hi_precision) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    *precision = lo_precision;
    return CERT_INTERVAL_OK;
}

static cert_interval_status commit(
    cert_interval *result,
    mpfr_srcptr lo,
    mpfr_srcptr hi
)
{
    mpfr_prec_t precision;
    mpfr_t rounded_lo;
    mpfr_t rounded_hi;
    cert_interval_status status;

    if (!mpfr_number_p(lo) || !mpfr_number_p(hi) || mpfr_cmp(lo, hi) > 0) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }

    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }

    /*
     * Most callers already supply temporaries at the destination precision.
     * A copy may instead narrow a higher-precision interval.  Round that case
     * into separate temporaries first so an exponent-range failure cannot
     * partially overwrite the destination promised unchanged on failure.
     */
    if (mpfr_get_prec(lo) != precision || mpfr_get_prec(hi) != precision) {
        mpfr_init2(rounded_lo, precision);
        mpfr_init2(rounded_hi, precision);
        mpfr_set(rounded_lo, lo, MPFR_RNDD);
        mpfr_set(rounded_hi, hi, MPFR_RNDU);
        if (!mpfr_number_p(rounded_lo) || !mpfr_number_p(rounded_hi) ||
            mpfr_cmp(rounded_lo, rounded_hi) > 0) {
            mpfr_clear(rounded_lo);
            mpfr_clear(rounded_hi);
            return CERT_INTERVAL_INVALID_INTERVAL;
        }
        mpfr_set(result->lo, rounded_lo, MPFR_RNDD);
        mpfr_set(result->hi, rounded_hi, MPFR_RNDU);
        mpfr_clear(rounded_lo);
        mpfr_clear(rounded_hi);
        return CERT_INTERVAL_OK;
    }

    mpfr_set(result->lo, lo, MPFR_RNDD);
    mpfr_set(result->hi, hi, MPFR_RNDU);
    return CERT_INTERVAL_OK;
}

const char *cert_interval_status_string(cert_interval_status status)
{
    switch (status) {
    case CERT_INTERVAL_OK:
        return "success";
    case CERT_INTERVAL_INVALID_ARGUMENT:
        return "invalid argument";
    case CERT_INTERVAL_INVALID_PRECISION:
        return "invalid MPFR precision";
    case CERT_INTERVAL_INVALID_INTERVAL:
        return "invalid interval";
    case CERT_INTERVAL_PARSE_ERROR:
        return "invalid exact decimal";
    case CERT_INTERVAL_DIVISION_BY_ZERO:
        return "division by an interval containing zero";
    case CERT_INTERVAL_NEGATIVE_RADICAND:
        return "square-root interval is not certified nonnegative";
    case CERT_INTERVAL_NEGATIVE_RADIUS:
        return "interval radius is not certified nonnegative";
    default:
        return "unknown interval status";
    }
}

cert_interval_status cert_interval_init(
    cert_interval *x,
    mpfr_prec_t precision
)
{
    if (x == NULL) {
        return CERT_INTERVAL_INVALID_ARGUMENT;
    }
    if (precision < MPFR_PREC_MIN || precision > MPFR_PREC_MAX) {
        return CERT_INTERVAL_INVALID_PRECISION;
    }
    mpfr_init2(x->lo, precision);
    mpfr_init2(x->hi, precision);
    return CERT_INTERVAL_OK;
}

void cert_interval_clear(cert_interval *x)
{
    if (x != NULL) {
        mpfr_clear(x->lo);
        mpfr_clear(x->hi);
    }
}

int cert_interval_valid(const cert_interval *x)
{
    if (x == NULL) {
        return 0;
    }
    return mpfr_get_prec(x->lo) == mpfr_get_prec(x->hi)
        && mpfr_number_p(x->lo)
        && mpfr_number_p(x->hi)
        && mpfr_cmp(x->lo, x->hi) <= 0;
}

cert_interval_status cert_interval_copy(
    cert_interval *result,
    const cert_interval *x
)
{
    mpfr_prec_t precision;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    if (destination_precision(result, &precision) != CERT_INTERVAL_OK) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    (void)precision;
    return commit(result, x->lo, x->hi);
}

cert_interval_status cert_interval_set_si(cert_interval *x, long value)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    status = destination_precision(x, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_set_si(lo, value, MPFR_RNDD);
    mpfr_set_si(hi, value, MPFR_RNDU);
    status = commit(x, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_set_exact_decimal(
    cert_interval *x,
    const char *decimal
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (decimal == NULL) {
        return CERT_INTERVAL_INVALID_ARGUMENT;
    }
    status = destination_precision(x, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    if (mpfr_set_str(lo, decimal, 10, MPFR_RNDD) != 0
        || mpfr_set_str(hi, decimal, 10, MPFR_RNDU) != 0
        || !mpfr_number_p(lo)
        || !mpfr_number_p(hi)) {
        mpfr_clear(lo);
        mpfr_clear(hi);
        return CERT_INTERVAL_PARSE_ERROR;
    }
    status = commit(x, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_add(
    cert_interval *result,
    const cert_interval *x,
    const cert_interval *y
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x) || !cert_interval_valid(y)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_add(lo, x->lo, y->lo, MPFR_RNDD);
    mpfr_add(hi, x->hi, y->hi, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_sub(
    cert_interval *result,
    const cert_interval *x,
    const cert_interval *y
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x) || !cert_interval_valid(y)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_sub(lo, x->lo, y->hi, MPFR_RNDD);
    mpfr_sub(hi, x->hi, y->lo, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_neg(
    cert_interval *result,
    const cert_interval *x
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_neg(lo, x->hi, MPFR_RNDD);
    mpfr_neg(hi, x->lo, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_mul(
    cert_interval *result,
    const cert_interval *x,
    const cert_interval *y
)
{
    mpfr_prec_t precision;
    mpfr_t lower[4];
    mpfr_t upper[4];
    mpfr_t lo;
    mpfr_t hi;
    mpfr_srcptr x_endpoint[2];
    mpfr_srcptr y_endpoint[2];
    cert_interval_status status;
    int index;
    int i;
    int j;

    if (!cert_interval_valid(x) || !cert_interval_valid(y)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    x_endpoint[0] = x->lo;
    x_endpoint[1] = x->hi;
    y_endpoint[0] = y->lo;
    y_endpoint[1] = y->hi;
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    for (i = 0; i < 4; ++i) {
        mpfr_init2(lower[i], precision);
        mpfr_init2(upper[i], precision);
    }

    index = 0;
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            mpfr_mul(
                lower[index],
                x_endpoint[i],
                y_endpoint[j],
                MPFR_RNDD
            );
            mpfr_mul(
                upper[index],
                x_endpoint[i],
                y_endpoint[j],
                MPFR_RNDU
            );
            ++index;
        }
    }
    mpfr_set(lo, lower[0], MPFR_RNDD);
    mpfr_set(hi, upper[0], MPFR_RNDU);
    for (i = 1; i < 4; ++i) {
        if (mpfr_cmp(lower[i], lo) < 0) {
            mpfr_set(lo, lower[i], MPFR_RNDD);
        }
        if (mpfr_cmp(upper[i], hi) > 0) {
            mpfr_set(hi, upper[i], MPFR_RNDU);
        }
    }
    status = commit(result, lo, hi);

    for (i = 0; i < 4; ++i) {
        mpfr_clear(lower[i]);
        mpfr_clear(upper[i]);
    }
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_mul_si(
    cert_interval *result,
    const cert_interval *x,
    long value
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    if (value >= 0) {
        mpfr_mul_si(lo, x->lo, value, MPFR_RNDD);
        mpfr_mul_si(hi, x->hi, value, MPFR_RNDU);
    } else {
        mpfr_mul_si(lo, x->hi, value, MPFR_RNDD);
        mpfr_mul_si(hi, x->lo, value, MPFR_RNDU);
    }
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_square(
    cert_interval *result,
    const cert_interval *x
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    mpfr_t candidate;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_init2(candidate, precision);
    if (mpfr_sgn(x->lo) >= 0) {
        mpfr_mul(lo, x->lo, x->lo, MPFR_RNDD);
        mpfr_mul(hi, x->hi, x->hi, MPFR_RNDU);
    } else if (mpfr_sgn(x->hi) <= 0) {
        mpfr_mul(lo, x->hi, x->hi, MPFR_RNDD);
        mpfr_mul(hi, x->lo, x->lo, MPFR_RNDU);
    } else {
        mpfr_set_zero(lo, 1);
        mpfr_mul(hi, x->lo, x->lo, MPFR_RNDU);
        mpfr_mul(candidate, x->hi, x->hi, MPFR_RNDU);
        if (mpfr_cmp(candidate, hi) > 0) {
            mpfr_set(hi, candidate, MPFR_RNDU);
        }
    }
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    mpfr_clear(candidate);
    return status;
}

int cert_interval_contains_zero(const cert_interval *x)
{
    if (!cert_interval_valid(x)) {
        return 1;
    }
    return mpfr_sgn(x->lo) <= 0 && mpfr_sgn(x->hi) >= 0;
}

cert_interval_status cert_interval_inv(
    cert_interval *result,
    const cert_interval *x
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    if (cert_interval_contains_zero(x)) {
        return CERT_INTERVAL_DIVISION_BY_ZERO;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_ui_div(lo, 1, x->hi, MPFR_RNDD);
    mpfr_ui_div(hi, 1, x->lo, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_div(
    cert_interval *result,
    const cert_interval *x,
    const cert_interval *y
)
{
    cert_interval inverse;
    cert_interval_status status;
    mpfr_prec_t precision;

    if (!cert_interval_valid(x) || !cert_interval_valid(y)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    status = cert_interval_init(&inverse, precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    status = cert_interval_inv(&inverse, y);
    if (status == CERT_INTERVAL_OK) {
        status = cert_interval_mul(result, x, &inverse);
    }
    cert_interval_clear(&inverse);
    return status;
}

cert_interval_status cert_interval_sqrt(
    cert_interval *result,
    const cert_interval *x
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    if (mpfr_sgn(x->lo) < 0) {
        return CERT_INTERVAL_NEGATIVE_RADICAND;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_sqrt(lo, x->lo, MPFR_RNDD);
    mpfr_sqrt(hi, x->hi, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_add_si(
    cert_interval *result,
    const cert_interval *x,
    long value
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_add_si(lo, x->lo, value, MPFR_RNDD);
    mpfr_add_si(hi, x->hi, value, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_ui_sub(
    cert_interval *result,
    unsigned long value,
    const cert_interval *x
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_ui_sub(lo, value, x->hi, MPFR_RNDD);
    mpfr_ui_sub(hi, value, x->lo, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_div_ui(
    cert_interval *result,
    const cert_interval *x,
    unsigned long value
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    if (value == 0UL) {
        return CERT_INTERVAL_DIVISION_BY_ZERO;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_div_ui(lo, x->lo, value, MPFR_RNDD);
    mpfr_div_ui(hi, x->hi, value, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_pow_ui(
    cert_interval *result,
    const cert_interval *x,
    unsigned long power
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    mpfr_t candidate;
    cert_interval_status status;

    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_init2(candidate, precision);
    if (power == 0UL) {
        mpfr_set_ui(lo, 1, MPFR_RNDD);
        mpfr_set_ui(hi, 1, MPFR_RNDU);
    } else if ((power & 1UL) != 0UL) {
        mpfr_pow_ui(lo, x->lo, power, MPFR_RNDD);
        mpfr_pow_ui(hi, x->hi, power, MPFR_RNDU);
    } else if (mpfr_sgn(x->lo) >= 0) {
        mpfr_pow_ui(lo, x->lo, power, MPFR_RNDD);
        mpfr_pow_ui(hi, x->hi, power, MPFR_RNDU);
    } else if (mpfr_sgn(x->hi) <= 0) {
        mpfr_pow_ui(lo, x->hi, power, MPFR_RNDD);
        mpfr_pow_ui(hi, x->lo, power, MPFR_RNDU);
    } else {
        mpfr_set_zero(lo, 1);
        mpfr_pow_ui(hi, x->lo, power, MPFR_RNDU);
        mpfr_pow_ui(candidate, x->hi, power, MPFR_RNDU);
        if (mpfr_cmp(candidate, hi) > 0) {
            mpfr_set(hi, candidate, MPFR_RNDU);
        }
    }
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    mpfr_clear(candidate);
    return status;
}

cert_interval_status cert_interval_expand(
    cert_interval *result,
    const cert_interval *center,
    const cert_interval *radius
)
{
    mpfr_prec_t precision;
    mpfr_t lo;
    mpfr_t hi;
    cert_interval_status status;

    if (!cert_interval_valid(center) || !cert_interval_valid(radius)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    if (mpfr_sgn(radius->lo) < 0) {
        return CERT_INTERVAL_NEGATIVE_RADIUS;
    }
    status = destination_precision(result, &precision);
    if (status != CERT_INTERVAL_OK) {
        return status;
    }
    mpfr_init2(lo, precision);
    mpfr_init2(hi, precision);
    mpfr_sub(lo, center->lo, radius->hi, MPFR_RNDD);
    mpfr_add(hi, center->hi, radius->hi, MPFR_RNDU);
    status = commit(result, lo, hi);
    mpfr_clear(lo);
    mpfr_clear(hi);
    return status;
}

cert_interval_status cert_interval_abs_upper(
    mpfr_ptr result,
    const cert_interval *x
)
{
    mpfr_prec_t precision;
    mpfr_t lo_abs;
    mpfr_t hi_abs;

    if (result == NULL) {
        return CERT_INTERVAL_INVALID_ARGUMENT;
    }
    if (!cert_interval_valid(x)) {
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    precision = mpfr_get_prec(result);
    mpfr_init2(lo_abs, precision);
    mpfr_init2(hi_abs, precision);
    mpfr_abs(lo_abs, x->lo, MPFR_RNDU);
    mpfr_abs(hi_abs, x->hi, MPFR_RNDU);
    if (!mpfr_number_p(lo_abs) || !mpfr_number_p(hi_abs)) {
        mpfr_clear(lo_abs);
        mpfr_clear(hi_abs);
        return CERT_INTERVAL_INVALID_INTERVAL;
    }
    mpfr_set(
        result,
        mpfr_cmp(lo_abs, hi_abs) >= 0 ? lo_abs : hi_abs,
        MPFR_RNDU
    );
    mpfr_clear(lo_abs);
    mpfr_clear(hi_abs);
    return CERT_INTERVAL_OK;
}
