#ifndef CERTIFICATE_INTERVAL_H
#define CERTIFICATE_INTERVAL_H

/*
 * Small outward-rounded MPFR interval kernel shared by the certificate
 * programs.  Each interval owns its two MPFR endpoints; callers choose the
 * precision when the interval is initialized.
 *
 * None of these functions prints a diagnostic or terminates the process.
 * Operations that can fail return a status, leaving the destination
 * unchanged on failure.  Unless stated otherwise, all interval operands and
 * destinations must first have been initialized with cert_interval_init().
 */

#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    mpfr_t lo;
    mpfr_t hi;
} cert_interval;

typedef enum {
    CERT_INTERVAL_OK = 0,
    CERT_INTERVAL_INVALID_ARGUMENT,
    CERT_INTERVAL_INVALID_PRECISION,
    CERT_INTERVAL_INVALID_INTERVAL,
    CERT_INTERVAL_PARSE_ERROR,
    CERT_INTERVAL_DIVISION_BY_ZERO,
    CERT_INTERVAL_NEGATIVE_RADICAND,
    CERT_INTERVAL_NEGATIVE_RADIUS
} cert_interval_status;

const char *cert_interval_status_string(cert_interval_status status);

cert_interval_status cert_interval_init(
    cert_interval *x,
    mpfr_prec_t precision
);
void cert_interval_clear(cert_interval *x);

/* True exactly when both endpoints are finite, equally precise, and ordered. */
int cert_interval_valid(const cert_interval *x);

cert_interval_status cert_interval_copy(
    cert_interval *result,
    const cert_interval *x
);
cert_interval_status cert_interval_set_si(cert_interval *x, long value);
cert_interval_status cert_interval_set_exact_decimal(
    cert_interval *x,
    const char *decimal
);

cert_interval_status cert_interval_add(
    cert_interval *result,
    const cert_interval *x,
    const cert_interval *y
);
cert_interval_status cert_interval_sub(
    cert_interval *result,
    const cert_interval *x,
    const cert_interval *y
);
cert_interval_status cert_interval_neg(
    cert_interval *result,
    const cert_interval *x
);
cert_interval_status cert_interval_mul(
    cert_interval *result,
    const cert_interval *x,
    const cert_interval *y
);
cert_interval_status cert_interval_mul_si(
    cert_interval *result,
    const cert_interval *x,
    long value
);
cert_interval_status cert_interval_square(
    cert_interval *result,
    const cert_interval *x
);

/* Invalid intervals conservatively return true from this predicate. */
int cert_interval_contains_zero(const cert_interval *x);

cert_interval_status cert_interval_inv(
    cert_interval *result,
    const cert_interval *x
);
cert_interval_status cert_interval_div(
    cert_interval *result,
    const cert_interval *x,
    const cert_interval *y
);
cert_interval_status cert_interval_sqrt(
    cert_interval *result,
    const cert_interval *x
);

cert_interval_status cert_interval_add_si(
    cert_interval *result,
    const cert_interval *x,
    long value
);
cert_interval_status cert_interval_ui_sub(
    cert_interval *result,
    unsigned long value,
    const cert_interval *x
);
cert_interval_status cert_interval_div_ui(
    cert_interval *result,
    const cert_interval *x,
    unsigned long value
);
cert_interval_status cert_interval_pow_ui(
    cert_interval *result,
    const cert_interval *x,
    unsigned long power
);

/* Enclose center + [-radius, radius], using radius->hi as the radius. */
cert_interval_status cert_interval_expand(
    cert_interval *result,
    const cert_interval *center,
    const cert_interval *radius
);

/* Store an upward-rounded upper bound for |x| in an initialized MPFR value. */
cert_interval_status cert_interval_abs_upper(
    mpfr_ptr result,
    const cert_interval *x
);

#ifdef __cplusplus
}
#endif

#endif
