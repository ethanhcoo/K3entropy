/*
 * This is the program behind Appendix C.4. It certifies the nine relative
 * homotopy classes shown in the top panel of Figure 4.
 *
 * First, it uses the outward-rounded interval arithmetic in
 * certificate_interval.c and the orbit data in k3_orbit_data.c to put the ten
 * marked points in disjoint rational boxes and certify their left-to-right
 * order. It also checks that stereographic projection is defined on the
 * required neighborhoods and covers the image of each thickened standard
 * segment by a chain of overlapping rational rectangles.
 *
 * Next, it chooses a rational polygonal path through each rectangle chain.
 * The checks on these polygons are exact: the program verifies containment,
 * embeddedness, and the transversality conditions, computes their arc-data,
 * simplifies it using arc_word.c, and compares the unreduced and reduced
 * words with the rows in arc_certificate_data.c. The reduced rows are used by
 * mclass.c.
 *
 * The orbit decimals in k3_orbit_data.c are treated as exact finite decimals.
 * If an interval or exact-geometric check cannot be certified, the program
 * exits unsuccessfully.
 */

#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include "arc_certificate_data.h"
#include "arc_word.h"
#include "certificate_interval.h"
#include "k3_orbit_data.h"

#ifndef PRECISION
#define PRECISION 500
#endif

#define ORBIT_LENGTH ARC_CERTIFICATE_ORBIT_LENGTH
#define ARC_COUNT ARC_CERTIFICATE_ARC_COUNT
#define INITIAL_PANEL_BITS 8
#define MAX_PANEL_DEPTH 28
#define RADIAL_BISECTIONS 36
#define RADIAL_GROWTH_STEPS 20
#define PUNCTURE_CENTER_BITS 64
#define PUNCTURE_RADIUS_BITS 45
#define INITIAL_IMAGE_WIDTH_BITS 5
#define MAX_POLYGON_ATTEMPTS 32
#define MAX_PANEL_COUNT 200000

typedef cert_interval interval;

typedef struct {
    interval x;
    interval y;
} interval_point2;

typedef struct {
    interval x;
    interval y;
    interval z;
} interval_point3;

/* Value and derivative intervals with respect to the segment parameter. */
typedef struct {
    interval value;
    interval derivative;
} dual;

typedef struct {
    dual x;
    dual y;
    dual z;
} dual_point3;

typedef struct {
    interval radius;
    interval q_radius_derivative;
} radial_certificate;

typedef struct {
    mpq_t x_min;
    mpq_t x_max;
    mpq_t y_min;
    mpq_t y_max;
} rational_box;

typedef struct {
    mpq_t x;
    mpq_t y;
} rational_point;

typedef struct {
    mpq_t t0;
    mpq_t t1;
    rational_box image;
} certified_panel;

typedef struct {
    certified_panel *items;
    size_t count;
    size_t capacity;
} panel_vector;

typedef struct {
    rational_point *vertices;
    size_t count;
} rational_path;

typedef struct {
    unsigned long panel_evaluations;
    unsigned long radial_refinements;
    int maximum_depth;
    mpfr_t minimum_projection_denominator;
    mpfr_t minimum_q_radius_derivative;
} arc_statistics;

static const int expected_tau[ORBIT_LENGTH] = {
    0, 2, 7, 3, 9, 6, 4, 1, 5, 8
};

static int active_arc = -1;
static int active_panel = -1;
static const char *inconclusive_reason = "unspecified interval dependency";
static int dump_certificate = 0;

static _Noreturn void fail(const char *format, ...)
{
    va_list arguments;

    fprintf(stderr, "CERTIFICATION FAILURE");
    if (active_arc >= 0) {
        fprintf(stderr, " (arc %d", active_arc);
        if (active_panel >= 0) {
            fprintf(stderr, ", panel %d", active_panel);
        }
        fputc(')', stderr);
    }
    fputs(": ", stderr);
    va_start(arguments, format);
    vfprintf(stderr, format, arguments);
    va_end(arguments);
    fputc('\n', stderr);
    exit(EXIT_FAILURE);
}

static void require_interval_status(
    cert_interval_status status,
    const char *context
)
{
    if (status != CERT_INTERVAL_OK) {
        fail("%s: %s", context, cert_interval_status_string(status));
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

static void intervals_init(interval *items, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        interval_init(&items[i]);
    }
}

static void intervals_clear(interval *items, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        interval_clear(&items[i]);
    }
}

static void point2_init(interval_point2 *point)
{
    interval_init(&point->x);
    interval_init(&point->y);
}

static void point2_clear(interval_point2 *point)
{
    interval_clear(&point->x);
    interval_clear(&point->y);
}

static void point3_init(interval_point3 *point)
{
    interval_init(&point->x);
    interval_init(&point->y);
    interval_init(&point->z);
}

static void point3_clear(interval_point3 *point)
{
    interval_clear(&point->x);
    interval_clear(&point->y);
    interval_clear(&point->z);
}

static void dual_init(dual *x)
{
    interval_init(&x->value);
    interval_init(&x->derivative);
}

static void dual_clear(dual *x)
{
    interval_clear(&x->value);
    interval_clear(&x->derivative);
}

static void duals_init(dual *items, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        dual_init(&items[i]);
    }
}

static void duals_clear(dual *items, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        dual_clear(&items[i]);
    }
}

static void dual_point3_init(dual_point3 *point)
{
    dual_init(&point->x);
    dual_init(&point->y);
    dual_init(&point->z);
}

static void dual_point3_clear(dual_point3 *point)
{
    dual_clear(&point->x);
    dual_clear(&point->y);
    dual_clear(&point->z);
}

static void radial_certificate_init(radial_certificate *certificate)
{
    interval_init(&certificate->radius);
    interval_init(&certificate->q_radius_derivative);
}

static void radial_certificate_clear(radial_certificate *certificate)
{
    interval_clear(&certificate->radius);
    interval_clear(&certificate->q_radius_derivative);
}

static void interval_require_valid(const interval *x, const char *context)
{
    if (!cert_interval_valid(x)) {
        fail("invalid interval in %s", context);
    }
}

static int interval_contains_zero(const interval *x)
{
    interval_require_valid(x, "zero-containment test");
    return cert_interval_contains_zero(x);
}

static void interval_copy(interval *result, const interval *x)
{
    require_interval_status(
        cert_interval_copy(result, x),
        "interval copy"
    );
}

static void interval_set_si(interval *x, long value)
{
    require_interval_status(
        cert_interval_set_si(x, value),
        "integer interval conversion"
    );
}

static void interval_set_exact_decimal(interval *x, const char *decimal)
{
    cert_interval_status status = cert_interval_set_exact_decimal(x, decimal);

    if (status != CERT_INTERVAL_OK) {
        fail(
            "invalid exact decimal '%s': %s",
            decimal,
            cert_interval_status_string(status)
        );
    }
}

static void interval_set_mpfr(interval *x, const mpfr_t value)
{
    mpfr_set(x->lo, value, MPFR_RNDD);
    mpfr_set(x->hi, value, MPFR_RNDU);
    interval_require_valid(x, "MPFR point-interval conversion");
}

static void interval_add(interval *result, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_add(result, x, y),
        "interval addition"
    );
}

static void interval_sub(interval *result, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_sub(result, x, y),
        "interval subtraction"
    );
}

static void interval_neg(interval *result, const interval *x)
{
    require_interval_status(
        cert_interval_neg(result, x),
        "interval negation"
    );
}

static void interval_mul(interval *result, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_mul(result, x, y),
        "interval multiplication"
    );
}

static void interval_mul_si(interval *result, const interval *x, long value)
{
    require_interval_status(
        cert_interval_mul_si(result, x, value),
        "integer interval multiplication"
    );
}

static void interval_square(interval *result, const interval *x)
{
    require_interval_status(
        cert_interval_square(result, x),
        "interval square"
    );
}

static void interval_inv(interval *result, const interval *x)
{
    require_interval_status(
        cert_interval_inv(result, x),
        "interval reciprocal"
    );
}

static void interval_div(interval *result, const interval *x, const interval *y)
{
    require_interval_status(
        cert_interval_div(result, x, y),
        "interval division"
    );
}

static void interval_sqrt(interval *result, const interval *x)
{
    require_interval_status(
        cert_interval_sqrt(result, x),
        "interval square root"
    );
}

static void interval_one_plus_square(interval *result, const interval *x)
{
    interval square;
    interval one;

    interval_init(&square);
    interval_init(&one);
    interval_square(&square, x);
    interval_set_si(&one, 1);
    interval_add(result, &one, &square);
    interval_clear(&one);
    interval_clear(&square);
}

static void interval_intersect(interval *result, const interval *candidate)
{
    interval_require_valid(result, "interval intersection accumulator");
    interval_require_valid(candidate, "interval intersection candidate");
    if (mpfr_cmp(candidate->lo, result->lo) > 0) {
        mpfr_set(result->lo, candidate->lo, MPFR_RNDD);
    }
    if (mpfr_cmp(candidate->hi, result->hi) < 0) {
        mpfr_set(result->hi, candidate->hi, MPFR_RNDU);
    }
    if (mpfr_cmp(result->lo, result->hi) > 0) {
        fail("two rigorous enclosures had empty intersection");
    }
    interval_require_valid(result, "interval intersection result");
}

static void dual_copy(dual *result, const dual *x)
{
    interval_copy(&result->value, &x->value);
    interval_copy(&result->derivative, &x->derivative);
}

static void dual_set_si(dual *x, long value)
{
    interval_set_si(&x->value, value);
    interval_set_si(&x->derivative, 0);
}

static void dual_add(dual *result, const dual *x, const dual *y)
{
    dual work;

    dual_init(&work);
    interval_add(&work.value, &x->value, &y->value);
    interval_add(&work.derivative, &x->derivative, &y->derivative);
    dual_copy(result, &work);
    dual_clear(&work);
}

static void dual_sub(dual *result, const dual *x, const dual *y)
{
    dual work;

    dual_init(&work);
    interval_sub(&work.value, &x->value, &y->value);
    interval_sub(&work.derivative, &x->derivative, &y->derivative);
    dual_copy(result, &work);
    dual_clear(&work);
}

static void dual_neg(dual *result, const dual *x)
{
    dual work;

    dual_init(&work);
    interval_neg(&work.value, &x->value);
    interval_neg(&work.derivative, &x->derivative);
    dual_copy(result, &work);
    dual_clear(&work);
}

static void dual_mul(dual *result, const dual *x, const dual *y)
{
    dual work;
    interval terms[2];

    dual_init(&work);
    intervals_init(terms, 2);
    interval_mul(&work.value, &x->value, &y->value);
    interval_mul(&terms[0], &x->derivative, &y->value);
    interval_mul(&terms[1], &x->value, &y->derivative);
    interval_add(&work.derivative, &terms[0], &terms[1]);
    dual_copy(result, &work);
    intervals_clear(terms, 2);
    dual_clear(&work);
}

static void dual_mul_si(dual *result, const dual *x, long value)
{
    dual work;

    dual_init(&work);
    interval_mul_si(&work.value, &x->value, value);
    interval_mul_si(&work.derivative, &x->derivative, value);
    dual_copy(result, &work);
    dual_clear(&work);
}

static void dual_square(dual *result, const dual *x)
{
    dual work;
    interval twice_value;

    dual_init(&work);
    interval_init(&twice_value);
    interval_square(&work.value, &x->value);
    interval_mul_si(&twice_value, &x->value, 2);
    interval_mul(&work.derivative, &twice_value, &x->derivative);
    dual_copy(result, &work);
    interval_clear(&twice_value);
    dual_clear(&work);
}

static void dual_inv(dual *result, const dual *x)
{
    dual work;
    interval square;
    interval quotient;

    dual_init(&work);
    interval_init(&square);
    interval_init(&quotient);
    interval_inv(&work.value, &x->value);
    interval_square(&square, &x->value);
    interval_div(&quotient, &x->derivative, &square);
    interval_neg(&work.derivative, &quotient);
    dual_copy(result, &work);
    interval_clear(&quotient);
    interval_clear(&square);
    dual_clear(&work);
}

static void dual_div(dual *result, const dual *x, const dual *y)
{
    dual inverse;
    dual work;

    dual_init(&inverse);
    dual_init(&work);
    dual_inv(&inverse, y);
    dual_mul(&work, x, &inverse);
    dual_copy(result, &work);
    dual_clear(&work);
    dual_clear(&inverse);
}

static void dual_sqrt(dual *result, const dual *x)
{
    dual work;
    interval denominator;

    if (mpfr_sgn(x->value.lo) <= 0) {
        fail("dual square root was not certified positive");
    }
    dual_init(&work);
    interval_init(&denominator);
    interval_sqrt(&work.value, &x->value);
    interval_mul_si(&denominator, &work.value, 2);
    interval_div(&work.derivative, &x->derivative, &denominator);
    dual_copy(result, &work);
    interval_clear(&denominator);
    dual_clear(&work);
}

static void dual_one_plus_square(dual *result, const dual *x)
{
    dual work[2];

    duals_init(work, 2);
    dual_square(&work[0], x);
    dual_set_si(&work[1], 1);
    dual_add(result, &work[1], &work[0]);
    duals_clear(work, 2);
}

static void rational_box_init(rational_box *box)
{
    mpq_inits(box->x_min, box->x_max, box->y_min, box->y_max, NULL);
}

static void rational_box_clear(rational_box *box)
{
    mpq_clears(box->x_min, box->x_max, box->y_min, box->y_max, NULL);
}

static void rational_box_copy(rational_box *result, const rational_box *box)
{
    mpq_set(result->x_min, box->x_min);
    mpq_set(result->x_max, box->x_max);
    mpq_set(result->y_min, box->y_min);
    mpq_set(result->y_max, box->y_max);
}

static void rational_point_init(rational_point *point)
{
    mpq_inits(point->x, point->y, NULL);
}

static void rational_point_clear(rational_point *point)
{
    mpq_clears(point->x, point->y, NULL);
}

static void rational_point_copy(rational_point *result, const rational_point *point)
{
    mpq_set(result->x, point->x);
    mpq_set(result->y, point->y);
}

static void rational_box_to_interval(
    interval_point2 *result,
    const rational_box *box
)
{
    mpfr_set_q(result->x.lo, box->x_min, MPFR_RNDD);
    mpfr_set_q(result->x.hi, box->x_max, MPFR_RNDU);
    mpfr_set_q(result->y.lo, box->y_min, MPFR_RNDD);
    mpfr_set_q(result->y.hi, box->y_max, MPFR_RNDU);
}

static void rational_box_midpoint(
    rational_point *result,
    const rational_box *box
)
{
    mpq_add(result->x, box->x_min, box->x_max);
    mpq_div_2exp(result->x, result->x, 1);
    mpq_add(result->y, box->y_min, box->y_max);
    mpq_div_2exp(result->y, result->y, 1);
}

static int rational_box_contains_point(
    const rational_box *box,
    const rational_point *point
)
{
    return mpq_cmp(box->x_min, point->x) <= 0 &&
           mpq_cmp(point->x, box->x_max) <= 0 &&
           mpq_cmp(box->y_min, point->y) <= 0 &&
           mpq_cmp(point->y, box->y_max) <= 0;
}

static int rational_box_contains_box(
    const rational_box *outer,
    const rational_box *inner
)
{
    return mpq_cmp(outer->x_min, inner->x_min) <= 0 &&
           mpq_cmp(inner->x_max, outer->x_max) <= 0 &&
           mpq_cmp(outer->y_min, inner->y_min) <= 0 &&
           mpq_cmp(inner->y_max, outer->y_max) <= 0;
}

static int rational_boxes_disjoint(
    const rational_box *left,
    const rational_box *right
)
{
    return mpq_cmp(left->x_max, right->x_min) < 0 ||
           mpq_cmp(right->x_max, left->x_min) < 0 ||
           mpq_cmp(left->y_max, right->y_min) < 0 ||
           mpq_cmp(right->y_max, left->y_min) < 0;
}

static void rational_box_hull(
    rational_box *result,
    const rational_box *other
)
{
    if (mpq_cmp(other->x_min, result->x_min) < 0) {
        mpq_set(result->x_min, other->x_min);
    }
    if (mpq_cmp(other->x_max, result->x_max) > 0) {
        mpq_set(result->x_max, other->x_max);
    }
    if (mpq_cmp(other->y_min, result->y_min) < 0) {
        mpq_set(result->y_min, other->y_min);
    }
    if (mpq_cmp(other->y_max, result->y_max) > 0) {
        mpq_set(result->y_max, other->y_max);
    }
}

static int rational_box_intersection(
    rational_box *result,
    const rational_box *a,
    const rational_box *b
)
{
    mpq_set(result->x_min, mpq_cmp(a->x_min, b->x_min) >= 0
                               ? a->x_min : b->x_min);
    mpq_set(result->x_max, mpq_cmp(a->x_max, b->x_max) <= 0
                               ? a->x_max : b->x_max);
    mpq_set(result->y_min, mpq_cmp(a->y_min, b->y_min) >= 0
                               ? a->y_min : b->y_min);
    mpq_set(result->y_max, mpq_cmp(a->y_max, b->y_max) <= 0
                               ? a->y_max : b->y_max);
    return mpq_cmp(result->x_min, result->x_max) <= 0 &&
           mpq_cmp(result->y_min, result->y_max) <= 0;
}

static void rational_box_fprint(FILE *stream, const rational_box *box)
{
    gmp_fprintf(
        stream,
        "[%Qd,%Qd] x [%Qd,%Qd]",
        box->x_min,
        box->x_max,
        box->y_min,
        box->y_max
    );
}

static void rational_point_fprint(FILE *stream, const rational_point *point)
{
    gmp_fprintf(stream, "(%Qd,%Qd)", point->x, point->y);
}

static void certified_panel_init(certified_panel *panel)
{
    mpq_inits(panel->t0, panel->t1, NULL);
    rational_box_init(&panel->image);
}

static void certified_panel_clear(certified_panel *panel)
{
    rational_box_clear(&panel->image);
    mpq_clears(panel->t0, panel->t1, NULL);
}

static void panel_vector_init(panel_vector *vector)
{
    vector->items = NULL;
    vector->count = 0;
    vector->capacity = 0;
}

static void panel_vector_clear(panel_vector *vector)
{
    for (size_t i = 0; i < vector->count; ++i) {
        certified_panel_clear(&vector->items[i]);
    }
    free(vector->items);
    vector->items = NULL;
    vector->count = 0;
    vector->capacity = 0;
}

static certified_panel *panel_vector_append(panel_vector *vector)
{
    certified_panel *new_items;

    if (vector->count >= MAX_PANEL_COUNT) {
        fail("panel count exceeded %d", MAX_PANEL_COUNT);
    }
    if (vector->count == vector->capacity) {
        size_t new_capacity = vector->capacity == 0 ? 512 : 2 * vector->capacity;

        if (new_capacity > MAX_PANEL_COUNT) {
            new_capacity = MAX_PANEL_COUNT;
        }
        new_items = realloc(vector->items, new_capacity * sizeof(*new_items));
        if (new_items == NULL) {
            fail("could not allocate panel vector");
        }
        vector->items = new_items;
        vector->capacity = new_capacity;
    }
    certified_panel_init(&vector->items[vector->count]);
    return &vector->items[vector->count++];
}

static void rational_path_init(rational_path *path, size_t count)
{
    path->vertices = calloc(count, sizeof(*path->vertices));
    if (path->vertices == NULL) {
        fail("could not allocate polygonal path");
    }
    path->count = count;
    for (size_t i = 0; i < count; ++i) {
        rational_point_init(&path->vertices[i]);
    }
}

static void rational_path_clear(rational_path *path)
{
    for (size_t i = 0; i < path->count; ++i) {
        rational_point_clear(&path->vertices[i]);
    }
    free(path->vertices);
    path->vertices = NULL;
    path->count = 0;
}

static void arc_statistics_init(arc_statistics *statistics)
{
    statistics->panel_evaluations = 0;
    statistics->radial_refinements = 0;
    statistics->maximum_depth = 0;
    mpfr_init2(statistics->minimum_projection_denominator, PRECISION);
    mpfr_init2(statistics->minimum_q_radius_derivative, PRECISION);
    mpfr_set_inf(statistics->minimum_projection_denominator, 1);
    mpfr_set_inf(statistics->minimum_q_radius_derivative, 1);
}

static void arc_statistics_clear(arc_statistics *statistics)
{
    mpfr_clear(statistics->minimum_projection_denominator);
    mpfr_clear(statistics->minimum_q_radius_derivative);
}

static void update_minimum(mpfr_t minimum, const mpfr_t value)
{
    if (mpfr_cmp(value, minimum) < 0) {
        mpfr_set(minimum, value, MPFR_RNDD);
    }
}

/* ------------------------------------------------------------------------- */
/* The K3 chart, radial, and automorphism formulas.                           */
/* ------------------------------------------------------------------------- */

static void discriminant(
    interval *result,
    const interval *u,
    const interval *v
)
{
    interval w[9];

    intervals_init(w, 9);
    interval_square(&w[0], u);
    interval_square(&w[1], v);
    interval_one_plus_square(&w[2], u);
    interval_one_plus_square(&w[3], v);
    interval_mul(&w[4], &w[0], &w[1]);
    interval_mul_si(&w[4], &w[4], 100);
    interval_mul(&w[5], &w[2], &w[3]);
    interval_mul_si(&w[5], &w[5], 8);
    interval_square(&w[6], &w[2]);
    interval_square(&w[7], &w[3]);
    interval_mul(&w[8], &w[6], &w[7]);
    interval_mul_si(&w[8], &w[8], 4);
    interval_add(&w[6], &w[4], &w[5]);
    interval_sub(result, &w[6], &w[8]);
    intervals_clear(w, 9);
}

static void chart_dependent_coordinate(
    interval *result,
    const interval *u,
    const interval *v,
    int sign
)
{
    interval w[8];

    intervals_init(w, 8);
    discriminant(&w[0], u, v);
    if (mpfr_sgn(w[0].lo) <= 0) {
        fail("chart discriminant was not certified positive");
    }
    interval_sqrt(&w[1], &w[0]);
    interval_mul(&w[2], u, v);
    interval_mul_si(&w[3], &w[2], -10);
    if (sign == 1) {
        interval_add(&w[4], &w[3], &w[1]);
    } else if (sign == -1) {
        interval_sub(&w[4], &w[3], &w[1]);
    } else {
        fail("invalid chart sign %d", sign);
    }
    interval_one_plus_square(&w[5], u);
    interval_one_plus_square(&w[6], v);
    interval_mul(&w[7], &w[5], &w[6]);
    interval_mul_si(&w[7], &w[7], 2);
    interval_div(result, &w[4], &w[7]);
    intervals_clear(w, 8);
}

static void chart_lift(
    interval_point3 *point,
    const interval *u,
    const interval *v,
    int number,
    int sign
)
{
    interval dependent;

    interval_init(&dependent);
    chart_dependent_coordinate(&dependent, u, v, sign);
    if (number == 1) {
        interval_copy(&point->x, &dependent);
        interval_copy(&point->y, u);
        interval_copy(&point->z, v);
    } else if (number == 2) {
        interval_copy(&point->x, u);
        interval_copy(&point->y, &dependent);
        interval_copy(&point->z, v);
    } else if (number == 3) {
        interval_copy(&point->x, u);
        interval_copy(&point->y, v);
        interval_copy(&point->z, &dependent);
    } else {
        fail("invalid chart number %d", number);
    }
    interval_clear(&dependent);
}

static void q_surface(interval *result, const interval_point3 *point)
{
    interval w[10];

    intervals_init(w, 10);
    interval_one_plus_square(&w[0], &point->x);
    interval_one_plus_square(&w[1], &point->y);
    interval_one_plus_square(&w[2], &point->z);
    interval_mul(&w[3], &w[0], &w[1]);
    interval_mul(&w[4], &w[3], &w[2]);
    interval_mul(&w[5], &point->x, &point->y);
    interval_mul(&w[6], &w[5], &point->z);
    interval_mul_si(&w[7], &w[6], 10);
    interval_add(&w[8], &w[4], &w[7]);
    interval_set_si(&w[9], 2);
    interval_sub(result, &w[8], &w[9]);
    intervals_clear(w, 10);
}

/* Apply R, then stereographic projection through gamma=(0,0,1). */
static int stereographic_project_interval(
    interval_point2 *result,
    interval *denominator_out,
    const interval_point3 *point
)
{
    interval w[10];
    const interval *denominator;

    intervals_init(w, 10);
    interval_square(&w[0], &point->x);
    interval_square(&w[1], &point->y);
    interval_square(&w[2], &point->z);
    interval_add(&w[3], &w[0], &w[1]);
    interval_add(&w[4], &w[3], &w[2]);
    if (mpfr_sgn(w[4].lo) <= 0) {
        inconclusive_reason = "radial norm was not certified positive";
        intervals_clear(w, 10);
        return 0;
    }
    interval_sqrt(&w[5], &w[4]);
    interval_sub(&w[6], &w[5], &point->z);

    /* Avoid cancellation in ||p||-z when z is positive. */
    interval_add(&w[7], &w[5], &point->z);
    denominator = &w[6];
    if (mpfr_sgn(w[3].lo) > 0 && mpfr_sgn(w[7].lo) > 0) {
        interval_div(&w[8], &w[3], &w[7]);
        denominator = &w[8];
    }
    if (mpfr_sgn(denominator->lo) <= 0) {
        inconclusive_reason =
            "stereographic denominator was not certified positive";
        intervals_clear(w, 10);
        return 0;
    }
    interval_div(&result->x, &point->x, denominator);
    interval_div(&result->y, &point->y, denominator);
    if (denominator_out != NULL) {
        interval_copy(denominator_out, denominator);
    }
    intervals_clear(w, 10);
    return 1;
}

static void dual_q_surface(dual *result, const dual_point3 *point)
{
    dual w[10];

    duals_init(w, 10);
    dual_one_plus_square(&w[0], &point->x);
    dual_one_plus_square(&w[1], &point->y);
    dual_one_plus_square(&w[2], &point->z);
    dual_mul(&w[3], &w[0], &w[1]);
    dual_mul(&w[4], &w[3], &w[2]);
    dual_mul(&w[5], &point->x, &point->y);
    dual_mul(&w[6], &w[5], &point->z);
    dual_mul_si(&w[7], &w[6], 10);
    dual_add(&w[8], &w[4], &w[7]);
    dual_set_si(&w[9], 2);
    dual_sub(result, &w[8], &w[9]);
    duals_clear(w, 10);
}

static void inverse_stereographic_direction(
    dual_point3 *direction,
    const dual *u,
    const dual *v
)
{
    dual w[4];

    /*
     * psi^{-1}(u,v) is a positive multiple of
     *     d=(2u,2v,u^2+v^2-1).
     * The omitted denominator is absorbed into the positive radial scale.
     */
    duals_init(w, 4);
    dual_square(&w[0], u);
    dual_square(&w[1], v);
    dual_add(&w[2], &w[0], &w[1]);
    dual_set_si(&w[3], 1);
    dual_mul_si(&direction->x, u, 2);
    dual_mul_si(&direction->y, v, 2);
    dual_sub(&direction->z, &w[2], &w[3]);
    duals_clear(w, 4);
}

/*
 * For d=(dx,dy,dz), evaluate
 * q(r d)=A r^2+10P r^3+B r^4+C r^6-1,
 * where A=dx^2+dy^2+dz^2,
 * B=dx^2 dy^2+dx^2 dz^2+dy^2 dz^2,
 * C=dx^2 dy^2 dz^2, and P=dx dy dz.
 */
static void radial_q_at_scalar(
    interval *result,
    const dual_point3 *direction,
    const mpfr_t radius
)
{
    interval w[25];

    intervals_init(w, 25);
    interval_square(&w[0], &direction->x.value);
    interval_square(&w[1], &direction->y.value);
    interval_square(&w[2], &direction->z.value);
    interval_mul(&w[3], &w[0], &w[1]);
    interval_mul(&w[4], &w[0], &w[2]);
    interval_mul(&w[5], &w[1], &w[2]);
    interval_add(&w[6], &w[3], &w[4]);
    interval_add(&w[7], &w[6], &w[5]);       /* B */
    interval_mul(&w[8], &w[3], &w[2]);       /* C */
    interval_mul(&w[9], &direction->x.value, &direction->y.value);
    interval_mul(&w[10], &w[9], &direction->z.value); /* P */

    /* Here A=dx^2+dy^2+dz^2=(dz+2)^2. */
    interval_set_si(&w[11], 2);
    interval_add(&w[12], &direction->z.value, &w[11]);
    interval_square(&w[13], &w[12]);

    interval_set_mpfr(&w[14], radius);
    interval_square(&w[15], &w[14]);
    interval_mul(&w[16], &w[15], &w[14]);
    interval_square(&w[17], &w[15]);
    interval_mul(&w[18], &w[17], &w[15]);
    interval_mul(&w[19], &w[13], &w[15]);
    interval_mul(&w[20], &w[7], &w[17]);
    interval_mul(&w[21], &w[8], &w[18]);
    interval_mul(&w[22], &w[10], &w[16]);
    interval_mul_si(&w[22], &w[22], 10);
    interval_set_si(&w[23], 1);
    interval_sub(&w[24], &w[19], &w[23]);
    interval_add(&w[24], &w[24], &w[20]);
    interval_add(&w[24], &w[24], &w[21]);
    interval_add(result, &w[24], &w[22]);
    intervals_clear(w, 25);
}

static void radial_q_derivative(
    interval *result,
    const dual_point3 *direction,
    const interval *radius
)
{
    interval w[24];

    intervals_init(w, 24);
    interval_square(&w[0], &direction->x.value);
    interval_square(&w[1], &direction->y.value);
    interval_square(&w[2], &direction->z.value);
    interval_mul(&w[3], &w[0], &w[1]);
    interval_mul(&w[4], &w[0], &w[2]);
    interval_mul(&w[5], &w[1], &w[2]);
    interval_add(&w[6], &w[3], &w[4]);
    interval_add(&w[7], &w[6], &w[5]);
    interval_mul(&w[8], &w[3], &w[2]);
    interval_mul(&w[9], &direction->x.value, &direction->y.value);
    interval_mul(&w[10], &w[9], &direction->z.value);
    interval_set_si(&w[11], 2);
    interval_add(&w[12], &direction->z.value, &w[11]);
    interval_square(&w[13], &w[12]);
    interval_square(&w[14], radius);
    interval_mul(&w[15], &w[14], radius);
    interval_mul(&w[16], &w[15], &w[14]);
    interval_mul(&w[17], &w[13], radius);
    interval_mul_si(&w[17], &w[17], 2);
    interval_mul(&w[18], &w[7], &w[15]);
    interval_mul_si(&w[18], &w[18], 4);
    interval_mul(&w[19], &w[8], &w[16]);
    interval_mul_si(&w[19], &w[19], 6);
    interval_mul(&w[20], &w[10], &w[14]);
    interval_mul_si(&w[20], &w[20], 30);
    interval_add(&w[21], &w[17], &w[18]);
    interval_add(&w[22], &w[21], &w[19]);
    interval_add(&w[23], &w[22], &w[20]);
    interval_copy(result, &w[23]);
    intervals_clear(w, 24);
}

/*
 * Uniformly bracket the positive radial root for every ray in direction.
 * The strict signs at the two endpoints, together with the global uniqueness
 * supplied by the star-shapedness lemma, prove containment. Interval Newton
 * only contracts that proven bracket; it is never used as an unvalidated root
 * finder.
 */
static int certify_radial_root(
    radial_certificate *certificate,
    const dual_point3 *direction,
    arc_statistics *statistics
)
{
    mpfr_t negative_bound;
    mpfr_t positive_bound;
    mpfr_t search_bound;
    mpfr_t midpoint;
    interval q_value;
    interval newton[4];
    int found_positive = 0;

    mpfr_inits2(
        PRECISION,
        negative_bound,
        positive_bound,
        search_bound,
        midpoint,
        (mpfr_ptr) 0
    );
    interval_init(&q_value);
    intervals_init(newton, 4);
    mpfr_set_zero(negative_bound, 1);
    mpfr_set_ui(positive_bound, 1, MPFR_RNDN);

    radial_q_at_scalar(&q_value, direction, negative_bound);
    if (mpfr_sgn(q_value.hi) >= 0) {
        fail("q(0) was not certified negative");
    }
    for (int growth = 0; growth < RADIAL_GROWTH_STEPS; ++growth) {
        radial_q_at_scalar(&q_value, direction, positive_bound);
        if (mpfr_sgn(q_value.lo) > 0) {
            found_positive = 1;
            break;
        }
        mpfr_mul_2ui(positive_bound, positive_bound, 1, MPFR_RNDN);
    }
    if (!found_positive) {
        inconclusive_reason = "could not certify a positive radial sign";
        goto inconclusive;
    }

    /* Refine the lower endpoint while preserving uniform negativity. */
    mpfr_set_zero(negative_bound, 1);
    mpfr_set(search_bound, positive_bound, MPFR_RNDN);
    for (int iteration = 0; iteration < RADIAL_BISECTIONS; ++iteration) {
        mpfr_add(midpoint, negative_bound, search_bound, MPFR_RNDN);
        mpfr_div_2ui(midpoint, midpoint, 1, MPFR_RNDN);
        radial_q_at_scalar(&q_value, direction, midpoint);
        if (mpfr_sgn(q_value.hi) < 0) {
            mpfr_set(negative_bound, midpoint, MPFR_RNDN);
        } else {
            mpfr_set(search_bound, midpoint, MPFR_RNDN);
        }
        ++statistics->radial_refinements;
    }

    /* Refine the upper endpoint while preserving uniform positivity. */
    mpfr_set_zero(search_bound, 1);
    for (int iteration = 0; iteration < RADIAL_BISECTIONS; ++iteration) {
        mpfr_add(midpoint, search_bound, positive_bound, MPFR_RNDN);
        mpfr_div_2ui(midpoint, midpoint, 1, MPFR_RNDN);
        radial_q_at_scalar(&q_value, direction, midpoint);
        if (mpfr_sgn(q_value.lo) > 0) {
            mpfr_set(positive_bound, midpoint, MPFR_RNDN);
        } else {
            mpfr_set(search_bound, midpoint, MPFR_RNDN);
        }
        ++statistics->radial_refinements;
    }

    if (mpfr_cmp(negative_bound, positive_bound) >= 0) {
        inconclusive_reason = "radial sign searches crossed";
        goto inconclusive;
    }
    mpfr_set(certificate->radius.lo, negative_bound, MPFR_RNDD);
    mpfr_set(certificate->radius.hi, positive_bound, MPFR_RNDU);
    radial_q_at_scalar(&q_value, direction, certificate->radius.lo);
    if (mpfr_sgn(q_value.hi) >= 0) {
        fail("retained lower radial sign was not strict");
    }
    radial_q_at_scalar(&q_value, direction, certificate->radius.hi);
    if (mpfr_sgn(q_value.lo) <= 0) {
        fail("retained upper radial sign was not strict");
    }

    radial_q_derivative(
        &certificate->q_radius_derivative,
        direction,
        &certificate->radius
    );
    if (mpfr_sgn(certificate->q_radius_derivative.lo) <= 0) {
        inconclusive_reason = "q_lambda was not positive on radial bracket";
        goto inconclusive;
    }

    for (int iteration = 0; iteration < 4; ++iteration) {
        mpfr_add(
            midpoint,
            certificate->radius.lo,
            certificate->radius.hi,
            MPFR_RNDN
        );
        mpfr_div_2ui(midpoint, midpoint, 1, MPFR_RNDN);
        radial_q_at_scalar(&q_value, direction, midpoint);
        interval_set_mpfr(&newton[0], midpoint);
        radial_q_derivative(&newton[1], direction, &certificate->radius);
        if (interval_contains_zero(&newton[1])) {
            inconclusive_reason = "interval Newton divisor contained zero";
            goto inconclusive;
        }
        interval_div(&newton[2], &q_value, &newton[1]);
        interval_sub(&newton[3], &newton[0], &newton[2]);
        interval_intersect(&certificate->radius, &newton[3]);
    }
    radial_q_derivative(
        &certificate->q_radius_derivative,
        direction,
        &certificate->radius
    );
    if (mpfr_sgn(certificate->q_radius_derivative.lo) <= 0) {
        inconclusive_reason = "q_lambda lost positivity after contraction";
        goto inconclusive;
    }

    update_minimum(
        statistics->minimum_q_radius_derivative,
        certificate->q_radius_derivative.lo
    );
    intervals_clear(newton, 4);
    interval_clear(&q_value);
    mpfr_clears(
        negative_bound,
        positive_bound,
        search_bound,
        midpoint,
        (mpfr_ptr) 0
    );
    return 1;

inconclusive:
    intervals_clear(newton, 4);
    interval_clear(&q_value);
    mpfr_clears(
        negative_bound,
        positive_bound,
        search_bound,
        midpoint,
        (mpfr_ptr) 0
    );
    return 0;
}

static int radial_lift(
    dual_point3 *result,
    radial_certificate *certificate,
    const dual_point3 *direction,
    arc_statistics *statistics
)
{
    dual radius;
    dual_point3 fixed_radius_point;
    dual partial_q;
    interval negative_partial;

    if (!certify_radial_root(certificate, direction, statistics)) {
        return 0;
    }
    dual_init(&radius);
    dual_point3_init(&fixed_radius_point);
    dual_init(&partial_q);
    interval_init(&negative_partial);

    interval_copy(&radius.value, &certificate->radius);
    interval_set_si(&radius.derivative, 0);
    dual_mul(&fixed_radius_point.x, &radius, &direction->x);
    dual_mul(&fixed_radius_point.y, &radius, &direction->y);
    dual_mul(&fixed_radius_point.z, &radius, &direction->z);
    dual_q_surface(&partial_q, &fixed_radius_point);

    /* Differentiate q(r(t)d(t))=0 to obtain r'=-q_t/q_r. */
    interval_neg(&negative_partial, &partial_q.derivative);
    interval_div(
        &radius.derivative,
        &negative_partial,
        &certificate->q_radius_derivative
    );
    dual_mul(&result->x, &radius, &direction->x);
    dual_mul(&result->y, &radius, &direction->y);
    dual_mul(&result->z, &radius, &direction->z);

    interval_clear(&negative_partial);
    dual_clear(&partial_q);
    dual_point3_clear(&fixed_radius_point);
    dual_clear(&radius);
    return 1;
}

static void dual_alpha(dual *result, const dual *u, const dual *v)
{
    dual w[4];

    duals_init(w, 4);
    dual_mul(&w[0], u, v);
    dual_one_plus_square(&w[1], u);
    dual_one_plus_square(&w[2], v);
    dual_mul(&w[3], &w[1], &w[2]);
    dual_div(result, &w[0], &w[3]);
    duals_clear(w, 4);
}

/* Apply f=sigma_3 o sigma_2 o sigma_1 in place. */
static void dual_apply_f(dual_point3 *point)
{
    dual alpha;
    dual correction;
    dual negative;
    dual replacement;

    dual_init(&alpha);
    dual_init(&correction);
    dual_init(&negative);
    dual_init(&replacement);

    dual_alpha(&alpha, &point->y, &point->z);
    dual_mul_si(&correction, &alpha, 10);
    dual_neg(&negative, &point->x);
    dual_sub(&replacement, &negative, &correction);
    dual_copy(&point->x, &replacement);

    dual_alpha(&alpha, &point->x, &point->z);
    dual_mul_si(&correction, &alpha, 10);
    dual_neg(&negative, &point->y);
    dual_sub(&replacement, &negative, &correction);
    dual_copy(&point->y, &replacement);

    dual_alpha(&alpha, &point->x, &point->y);
    dual_mul_si(&correction, &alpha, 10);
    dual_neg(&negative, &point->z);
    dual_sub(&replacement, &negative, &correction);
    dual_copy(&point->z, &replacement);

    dual_clear(&replacement);
    dual_clear(&negative);
    dual_clear(&correction);
    dual_clear(&alpha);
}

static int dual_stereographic_project(
    dual *u,
    dual *v,
    interval *denominator_out,
    const dual_point3 *point
)
{
    dual w[10];
    const dual *denominator;

    duals_init(w, 10);
    dual_square(&w[0], &point->x);
    dual_square(&w[1], &point->y);
    dual_square(&w[2], &point->z);
    dual_add(&w[3], &w[0], &w[1]);
    dual_add(&w[4], &w[3], &w[2]);
    if (mpfr_sgn(w[4].value.lo) <= 0) {
        inconclusive_reason = "image norm was not certified positive";
        duals_clear(w, 10);
        return 0;
    }
    dual_sqrt(&w[5], &w[4]);
    dual_sub(&w[6], &w[5], &point->z);
    dual_add(&w[7], &w[5], &point->z);
    denominator = &w[6];
    if (mpfr_sgn(w[3].value.lo) > 0 && mpfr_sgn(w[7].value.lo) > 0) {
        dual_div(&w[8], &w[3], &w[7]);
        denominator = &w[8];
    }
    if (mpfr_sgn(denominator->value.lo) <= 0) {
        inconclusive_reason =
            "image stereographic denominator was not certified positive";
        duals_clear(w, 10);
        return 0;
    }
    dual_div(u, &point->x, denominator);
    dual_div(v, &point->y, denominator);
    interval_copy(denominator_out, &denominator->value);
    duals_clear(w, 10);
    return 1;
}

static void affine_segment_coordinate(
    dual *result,
    const interval *left,
    const interval *right,
    const interval *parameter
)
{
    interval difference;
    interval product;

    interval_init(&difference);
    interval_init(&product);
    interval_sub(&difference, right, left);
    interval_mul(&product, parameter, &difference);
    interval_add(&result->value, left, &product);
    interval_copy(&result->derivative, &difference);
    interval_clear(&product);
    interval_clear(&difference);
}

/* Enclose Pi((1-t)left+t right) for every t in parameter. */
static int evaluate_pi_panel(
    dual *image_x,
    dual *image_y,
    interval *projection_denominator,
    radial_certificate *radial,
    const interval_point2 *left,
    const interval_point2 *right,
    const interval *parameter,
    arc_statistics *statistics
)
{
    dual u;
    dual v;
    dual_point3 direction;
    dual_point3 surface_point;
    int success = 0;

    ++statistics->panel_evaluations;
    dual_init(&u);
    dual_init(&v);
    dual_point3_init(&direction);
    dual_point3_init(&surface_point);

    affine_segment_coordinate(&u, &left->x, &right->x, parameter);
    affine_segment_coordinate(&v, &left->y, &right->y, parameter);
    inverse_stereographic_direction(&direction, &u, &v);
    if (!radial_lift(&surface_point, radial, &direction, statistics)) {
        goto cleanup;
    }
    dual_apply_f(&surface_point);
    if (!dual_stereographic_project(
            image_x,
            image_y,
            projection_denominator,
            &surface_point
        )) {
        goto cleanup;
    }
    interval_require_valid(&image_x->value, "Pi first coordinate");
    interval_require_valid(&image_y->value, "Pi second coordinate");
    interval_require_valid(&image_x->derivative, "Pi first derivative");
    interval_require_valid(&image_y->derivative, "Pi second derivative");
    update_minimum(
        statistics->minimum_projection_denominator,
        projection_denominator->lo
    );
    success = 1;

cleanup:
    dual_point3_clear(&surface_point);
    dual_point3_clear(&direction);
    dual_clear(&v);
    dual_clear(&u);
    return success;
}

static void sharpen_image_by_mean_value(
    dual *image_x,
    dual *image_y,
    const dual *midpoint_x,
    const dual *midpoint_y,
    const interval *parameter,
    const mpfr_t midpoint
)
{
    interval delta;
    interval variation;
    interval centered;

    interval_init(&delta);
    interval_init(&variation);
    interval_init(&centered);
    mpfr_sub(delta.lo, parameter->lo, midpoint, MPFR_RNDD);
    mpfr_sub(delta.hi, parameter->hi, midpoint, MPFR_RNDU);
    interval_mul(&variation, &image_x->derivative, &delta);
    interval_add(&centered, &midpoint_x->value, &variation);
    interval_intersect(&image_x->value, &centered);
    interval_mul(&variation, &image_y->derivative, &delta);
    interval_add(&centered, &midpoint_y->value, &variation);
    interval_intersect(&image_y->value, &centered);
    interval_clear(&centered);
    interval_clear(&variation);
    interval_clear(&delta);
}

/* ------------------------------------------------------------------------- */
/* Step 1: certified rational puncture boxes.                                 */
/* ------------------------------------------------------------------------- */

static void make_shadowing_coordinate(
    interval *result,
    const char *center_decimal,
    const interval *symmetric_error
)
{
    interval center;

    interval_init(&center);
    interval_set_exact_decimal(&center, center_decimal);
    interval_add(result, &center, symmetric_error);
    interval_clear(&center);
}

static void make_dyadic_axis_box(
    mpq_t lower,
    mpq_t upper,
    const interval *enclosure
)
{
    mpfr_t midpoint;
    mpfr_t scaled;
    mpz_t center_integer;
    mpq_t center;
    mpq_t radius;

    mpfr_init2(midpoint, PRECISION);
    mpfr_init2(scaled, PRECISION);
    mpz_init(center_integer);
    mpq_inits(center, radius, NULL);

    mpfr_add(midpoint, enclosure->lo, enclosure->hi, MPFR_RNDN);
    mpfr_div_2ui(midpoint, midpoint, 1, MPFR_RNDN);
    mpfr_mul_2ui(scaled, midpoint, PUNCTURE_CENTER_BITS, MPFR_RNDN);
    mpfr_get_z(center_integer, scaled, MPFR_RNDN);
    mpq_set_z(center, center_integer);
    mpq_div_2exp(center, center, PUNCTURE_CENTER_BITS);
    mpq_set_ui(radius, 1, 1);
    mpq_div_2exp(radius, radius, PUNCTURE_RADIUS_BITS);
    mpq_sub(lower, center, radius);
    mpq_add(upper, center, radius);

    if (mpfr_cmp_q(enclosure->lo, lower) <= 0 ||
        mpfr_cmp_q(enclosure->hi, upper) >= 0) {
        fail(
            "fixed dyadic puncture radius 2^(-%d) did not contain enclosure",
            PUNCTURE_RADIUS_BITS
        );
    }

    mpq_clears(center, radius, NULL);
    mpz_clear(center_integer);
    mpfr_clear(scaled);
    mpfr_clear(midpoint);
}

static void certify_puncture_boxes(
    rational_box punctures[ORBIT_LENGTH],
    rational_point reference_points[ORBIT_LENGTH],
    int tau[ORBIT_LENGTH],
    int inverse_tau[ORBIT_LENGTH],
    int dynamics[ORBIT_LENGTH]
)
{
    interval_point2 temporal_enclosures[ORBIT_LENGTH];
    rational_box temporal_boxes[ORBIT_LENGTH];
    int left_to_temporal[ORBIT_LENGTH];
    interval eta;
    interval error;
    interval u;
    interval v;
    interval denominator;

    for (int j = 0; j < ORBIT_LENGTH; ++j) {
        point2_init(&temporal_enclosures[j]);
        rational_box_init(&temporal_boxes[j]);
    }
    interval_init(&eta);
    interval_init(&error);
    interval_init(&u);
    interval_init(&v);
    interval_init(&denominator);
    interval_set_exact_decimal(&eta, "5e-27");
    mpfr_neg(error.lo, eta.hi, MPFR_RNDD);
    mpfr_set(error.hi, eta.hi, MPFR_RNDU);

    for (int j = 0; j < ORBIT_LENGTH; ++j) {
        interval_point3 surface;
        interval surface_equation;
        int number;
        int sign;

        point3_init(&surface);
        interval_init(&surface_equation);
        if (!k3_orbit_chart_at((size_t)j, &number, &sign)) {
            fail("missing shared chart data at temporal index %d", j);
        }
        make_shadowing_coordinate(&u, k3_orbit_a_decimal[j], &error);
        make_shadowing_coordinate(&v, k3_orbit_b_decimal[j], &error);
        chart_lift(&surface, &u, &v, number, sign);
        q_surface(&surface_equation, &surface);
        if (!interval_contains_zero(&surface_equation)) {
            fail("chart %d did not enclose q=0", j);
        }
        if (!stereographic_project_interval(
                &temporal_enclosures[j],
                &denominator,
                &surface
            )) {
            fail(
                "chart %d did not certify avoidance of gamma: %s",
                j,
                inconclusive_reason
            );
        }
        make_dyadic_axis_box(
            temporal_boxes[j].x_min,
            temporal_boxes[j].x_max,
            &temporal_enclosures[j].x
        );
        make_dyadic_axis_box(
            temporal_boxes[j].y_min,
            temporal_boxes[j].y_max,
            &temporal_enclosures[j].y
        );
        interval_clear(&surface_equation);
        point3_clear(&surface);
    }

    /* Sort only to propose an order; all strict separations are rechecked. */
    for (int i = 0; i < ORBIT_LENGTH; ++i) {
        left_to_temporal[i] = i;
    }
    for (int i = 1; i < ORBIT_LENGTH; ++i) {
        int candidate = left_to_temporal[i];
        int position = i;

        while (position > 0 &&
               mpq_cmp(
                   temporal_boxes[candidate].x_min,
                   temporal_boxes[left_to_temporal[position - 1]].x_min
               ) < 0) {
            left_to_temporal[position] = left_to_temporal[position - 1];
            --position;
        }
        left_to_temporal[position] = candidate;
    }

    for (int label = 0; label < ORBIT_LENGTH; ++label) {
        int temporal = left_to_temporal[label];

        rational_box_copy(&punctures[label], &temporal_boxes[temporal]);
        rational_box_midpoint(&reference_points[label], &punctures[label]);
        if (mpq_cmp(punctures[label].x_min, reference_points[label].x) >= 0 ||
            mpq_cmp(reference_points[label].x, punctures[label].x_max) >= 0 ||
            mpq_cmp(punctures[label].y_min, reference_points[label].y) >= 0 ||
            mpq_cmp(reference_points[label].y, punctures[label].y_max) >= 0) {
            fail("reference point for B_%d was not strictly interior", label);
        }
        tau[temporal] = label;
        inverse_tau[label] = temporal;
    }
    for (int label = 0; label + 1 < ORBIT_LENGTH; ++label) {
        if (mpq_cmp(
                punctures[label].x_max,
                punctures[label + 1].x_min
            ) >= 0) {
            fail("puncture boxes %d and %d were not strictly ordered",
                 label, label + 1);
        }
    }
    for (int i = 0; i < ORBIT_LENGTH; ++i) {
        for (int j = i + 1; j < ORBIT_LENGTH; ++j) {
            if (!rational_boxes_disjoint(&punctures[i], &punctures[j])) {
                fail("puncture boxes %d and %d were not disjoint", i, j);
            }
        }
    }
    for (int temporal = 0; temporal < ORBIT_LENGTH; ++temporal) {
        if (tau[temporal] != expected_tau[temporal]) {
            fail("certified tau disagrees with expected label at time %d",
                 temporal);
        }
    }
    for (int label = 0; label < ORBIT_LENGTH; ++label) {
        int temporal = inverse_tau[label];
        dynamics[label] = tau[(temporal + 1) % ORBIT_LENGTH];
    }

    printf("Step 1: certified ten rational puncture boxes and gamma avoidance.\n");
    printf("  tau (temporal -> left label):");
    for (int j = 0; j < ORBIT_LENGTH; ++j) {
        printf(" %d", tau[j]);
    }
    fputc('\n', stdout);
    for (int label = 0; label < ORBIT_LENGTH; ++label) {
        printf("  B_%d (time %d) = ", label, inverse_tau[label]);
        rational_box_fprint(stdout, &punctures[label]);
        fputc('\n', stdout);
    }

    interval_clear(&denominator);
    interval_clear(&v);
    interval_clear(&u);
    interval_clear(&error);
    interval_clear(&eta);
    for (int j = 0; j < ORBIT_LENGTH; ++j) {
        rational_box_clear(&temporal_boxes[j]);
        point2_clear(&temporal_enclosures[j]);
    }
}

/* ------------------------------------------------------------------------- */
/* Step 2a: adaptive rational rectangle chains.                               */
/* ------------------------------------------------------------------------- */

static void parameter_interval_from_rationals(
    interval *parameter,
    const mpq_t t0,
    const mpq_t t1
)
{
    if (mpq_cmp(t0, t1) > 0) {
        fail("reversed rational parameter interval");
    }
    mpfr_set_q(parameter->lo, t0, MPFR_RNDD);
    mpfr_set_q(parameter->hi, t1, MPFR_RNDU);
}

static void rational_box_from_image(
    rational_box *box,
    const interval *x,
    const interval *y
)
{
    mpfr_t endpoint;

    interval_require_valid(x, "rational x hull");
    interval_require_valid(y, "rational y hull");
    mpfr_init2(endpoint, PRECISION);

    mpfr_set(endpoint, x->lo, MPFR_RNDD);
    mpfr_nextbelow(endpoint);
    mpfr_get_q(box->x_min, endpoint);
    mpfr_set(endpoint, x->hi, MPFR_RNDU);
    mpfr_nextabove(endpoint);
    mpfr_get_q(box->x_max, endpoint);
    mpfr_set(endpoint, y->lo, MPFR_RNDD);
    mpfr_nextbelow(endpoint);
    mpfr_get_q(box->y_min, endpoint);
    mpfr_set(endpoint, y->hi, MPFR_RNDU);
    mpfr_nextabove(endpoint);
    mpfr_get_q(box->y_max, endpoint);
    mpfr_clear(endpoint);
}

static int interval_width_at_most_power_of_two(
    const interval *x,
    int negative_exponent
)
{
    mpfr_t width;
    int result;

    mpfr_init2(width, PRECISION);
    mpfr_sub(width, x->hi, x->lo, MPFR_RNDU);
    result = mpfr_cmp_ui_2exp(width, 1, -negative_exponent) <= 0;
    mpfr_clear(width);
    return result;
}

static int image_avoids_forbidden_punctures(
    const rational_box *image,
    const rational_box punctures[ORBIT_LENGTH],
    int allowed_label
)
{
    for (int label = 0; label < ORBIT_LENGTH; ++label) {
        if (label != allowed_label &&
            !rational_boxes_disjoint(image, &punctures[label])) {
            return 0;
        }
    }
    return 1;
}

static void evaluate_and_sharpen_panel(
    dual *image_x,
    dual *image_y,
    interval *projection_denominator,
    radial_certificate *radial,
    const interval_point2 *left,
    const interval_point2 *right,
    const interval *parameter,
    arc_statistics *statistics,
    int *success
)
{
    dual midpoint_x;
    dual midpoint_y;
    interval midpoint_parameter;
    interval midpoint_denominator;
    radial_certificate midpoint_radial;
    mpfr_t midpoint;

    dual_init(&midpoint_x);
    dual_init(&midpoint_y);
    interval_init(&midpoint_parameter);
    interval_init(&midpoint_denominator);
    radial_certificate_init(&midpoint_radial);
    mpfr_init2(midpoint, PRECISION);
    *success = 0;

    if (!evaluate_pi_panel(
            image_x,
            image_y,
            projection_denominator,
            radial,
            left,
            right,
            parameter,
            statistics
        )) {
        goto cleanup;
    }
    mpfr_add(midpoint, parameter->lo, parameter->hi, MPFR_RNDN);
    mpfr_div_2ui(midpoint, midpoint, 1, MPFR_RNDN);
    interval_set_mpfr(&midpoint_parameter, midpoint);
    if (!evaluate_pi_panel(
            &midpoint_x,
            &midpoint_y,
            &midpoint_denominator,
            &midpoint_radial,
            left,
            right,
            &midpoint_parameter,
            statistics
        )) {
        goto cleanup;
    }
    sharpen_image_by_mean_value(
        image_x,
        image_y,
        &midpoint_x,
        &midpoint_y,
        parameter,
        midpoint
    );
    *success = 1;

cleanup:
    mpfr_clear(midpoint);
    radial_certificate_clear(&midpoint_radial);
    interval_clear(&midpoint_denominator);
    interval_clear(&midpoint_parameter);
    dual_clear(&midpoint_y);
    dual_clear(&midpoint_x);
}

/*
 * Subdivide until each image enclosure has the required width and avoids
 * every puncture box except the permitted box at an endpoint.
 */
static void certify_panel_recursive(
    panel_vector *panels,
    const interval_point2 *left,
    const interval_point2 *right,
    const rational_box punctures[ORBIT_LENGTH],
    int start_target,
    int end_target,
    const mpq_t t0,
    const mpq_t t1,
    int depth,
    int image_width_bits,
    arc_statistics *statistics
)
{
    interval parameter;
    interval projection_denominator;
    dual image_x;
    dual image_y;
    radial_certificate radial;
    rational_box image_box;
    mpq_t midpoint;
    int success;
    int touches_start = mpq_sgn(t0) == 0;
    int touches_end = mpq_cmp_ui(t1, 1, 1) == 0;
    int allowed_label = touches_start ? start_target
                         : (touches_end ? end_target : -1);
    int acceptable;

    if (depth > statistics->maximum_depth) {
        statistics->maximum_depth = depth;
    }
    interval_init(&parameter);
    interval_init(&projection_denominator);
    dual_init(&image_x);
    dual_init(&image_y);
    radial_certificate_init(&radial);
    rational_box_init(&image_box);
    mpq_init(midpoint);
    parameter_interval_from_rationals(&parameter, t0, t1);
    evaluate_and_sharpen_panel(
        &image_x,
        &image_y,
        &projection_denominator,
        &radial,
        left,
        right,
        &parameter,
        statistics,
        &success
    );

    acceptable = success;
    if (acceptable) {
        rational_box_from_image(
            &image_box,
            &image_x.value,
            &image_y.value
        );
        if (touches_start) {
            rational_box_hull(&image_box, &punctures[start_target]);
        }
        if (touches_end) {
            rational_box_hull(&image_box, &punctures[end_target]);
        }
        acceptable =
            image_avoids_forbidden_punctures(
                &image_box,
                punctures,
                allowed_label
            ) &&
            interval_width_at_most_power_of_two(
                &image_x.value,
                image_width_bits
            ) &&
            interval_width_at_most_power_of_two(
                &image_y.value,
                image_width_bits
            );
    }

    if (acceptable) {
        certified_panel *accepted = panel_vector_append(panels);

        mpq_set(accepted->t0, t0);
        mpq_set(accepted->t1, t1);
        rational_box_copy(&accepted->image, &image_box);
    } else {
        if (depth >= MAX_PANEL_DEPTH) {
            fail(
                "subdivision limit reached (%s)",
                success ? "image rectangle remained too coarse or met a puncture"
                        : inconclusive_reason
            );
        }
        mpq_add(midpoint, t0, t1);
        mpq_div_2exp(midpoint, midpoint, 1);
        if (mpq_cmp(midpoint, t0) <= 0 || mpq_cmp(midpoint, t1) >= 0) {
            fail("could not form a strict rational panel midpoint");
        }
        certify_panel_recursive(
            panels,
            left,
            right,
            punctures,
            start_target,
            end_target,
            t0,
            midpoint,
            depth + 1,
            image_width_bits,
            statistics
        );
        certify_panel_recursive(
            panels,
            left,
            right,
            punctures,
            start_target,
            end_target,
            midpoint,
            t1,
            depth + 1,
            image_width_bits,
            statistics
        );
    }

    mpq_clear(midpoint);
    rational_box_clear(&image_box);
    radial_certificate_clear(&radial);
    dual_clear(&image_y);
    dual_clear(&image_x);
    interval_clear(&projection_denominator);
    interval_clear(&parameter);
}

static void build_rectangle_chain(
    panel_vector *panels,
    const rational_box punctures[ORBIT_LENGTH],
    int input_left,
    int input_right,
    int start_target,
    int end_target,
    int image_width_bits,
    arc_statistics *statistics
)
{
    interval_point2 left;
    interval_point2 right;
    mpq_t t0;
    mpq_t t1;
    const unsigned long initial_count = 1UL << INITIAL_PANEL_BITS;

    point2_init(&left);
    point2_init(&right);
    mpq_inits(t0, t1, NULL);
    rational_box_to_interval(&left, &punctures[input_left]);
    rational_box_to_interval(&right, &punctures[input_right]);

    for (unsigned long index = 0; index < initial_count; ++index) {
        active_panel = (int) index;
        mpq_set_ui(t0, index, initial_count);
        mpq_set_ui(t1, index + 1, initial_count);
        certify_panel_recursive(
            panels,
            &left,
            &right,
            punctures,
            start_target,
            end_target,
            t0,
            t1,
            0,
            image_width_bits,
            statistics
        );
    }
    active_panel = -1;
    if (panels->count < 2) {
        fail("rectangle chain has fewer than two panels");
    }
    if (mpq_sgn(panels->items[0].t0) != 0 ||
        mpq_cmp_ui(panels->items[panels->count - 1].t1, 1, 1) != 0) {
        fail("rectangle chain does not cover [0,1]");
    }
    if (!rational_box_contains_box(
            &panels->items[0].image,
            &punctures[start_target]
        ) ||
        !rational_box_contains_box(
            &panels->items[panels->count - 1].image,
            &punctures[end_target]
        )) {
        fail("endpoint image rectangles do not contain their puncture boxes");
    }
    for (size_t i = 0; i + 1 < panels->count; ++i) {
        rational_box overlap;

        if (mpq_cmp(panels->items[i].t1, panels->items[i + 1].t0) != 0) {
            fail("rectangle-chain parameters are not consecutive");
        }
        rational_box_init(&overlap);
        if (!rational_box_intersection(
                &overlap,
                &panels->items[i].image,
                &panels->items[i + 1].image
            )) {
            rational_box_clear(&overlap);
            fail("consecutive certified image rectangles do not overlap");
        }
        if (mpq_cmp(overlap.x_min, overlap.x_max) >= 0 ||
            mpq_cmp(overlap.y_min, overlap.y_max) >= 0) {
            rational_box_clear(&overlap);
            fail("consecutive rectangles lack an open rational overlap");
        }
        rational_box_clear(&overlap);
    }

    mpq_clears(t0, t1, NULL);
    point2_clear(&right);
    point2_clear(&left);
}

/* ------------------------------------------------------------------------- */
/* Step 2b: exact rational polygonal geometry and arc-data.                   */
/* ------------------------------------------------------------------------- */

/*
 * The reference point bar_z_i is the rational midpoint of B_i. Since the
 * boxes are disjoint and both z_i and bar_z_i lie in int(B_i), the local
 * homeomorphisms taking z_i to bar_z_i have disjoint supports. Their union is
 * the map theta in the manuscript. Only its support is used here.
 */

static void rational_affine_coordinate(
    mpq_t result,
    const mpq_t lower,
    const mpq_t upper,
    unsigned long numerator,
    unsigned long denominator
)
{
    mpq_t width;
    mpq_t fraction;

    if (denominator == 0 || numerator > denominator) {
        fail("invalid rational affine fraction");
    }
    mpq_inits(width, fraction, NULL);
    mpq_sub(width, upper, lower);
    mpq_set_ui(fraction, numerator, denominator);
    mpq_mul(width, width, fraction);
    mpq_add(result, lower, width);
    mpq_clears(width, fraction, NULL);
}

static void choose_overlap_point(
    rational_point *point,
    const rational_box *overlap,
    unsigned int attempt,
    size_t boundary_index
)
{
    unsigned long x_numerator;
    unsigned long y_numerator;
    const unsigned long denominator = 1009;

    if (mpq_cmp(overlap->x_min, overlap->x_max) >= 0 ||
        mpq_cmp(overlap->y_min, overlap->y_max) >= 0) {
        fail("cannot choose a point in a degenerate overlap");
    }
    if (attempt == 0) {
        x_numerator = 1;
        y_numerator = 1;
        rational_affine_coordinate(
            point->x,
            overlap->x_min,
            overlap->x_max,
            x_numerator,
            2
        );
        rational_affine_coordinate(
            point->y,
            overlap->y_min,
            overlap->y_max,
            y_numerator,
            2
        );
        return;
    }

    /* Try deterministic rational points in the interior of the overlap. */
    x_numerator = 1 +
        ((unsigned long) attempt * 421UL +
         (unsigned long) boundary_index * 173UL) % (denominator - 1);
    y_numerator = 1 +
        ((unsigned long) attempt * 769UL +
         (unsigned long) boundary_index * 317UL + 97UL) %
        (denominator - 1);
    rational_affine_coordinate(
        point->x,
        overlap->x_min,
        overlap->x_max,
        x_numerator,
        denominator
    );
    rational_affine_coordinate(
        point->y,
        overlap->y_min,
        overlap->y_max,
        y_numerator,
        denominator
    );
}

static int build_polygon_candidate(
    rational_path *path,
    const panel_vector *panels,
    const rational_point reference_points[ORBIT_LENGTH],
    int start_label,
    int end_label,
    unsigned int attempt
)
{
    rational_box overlap;

    if (path->count != panels->count + 1) {
        fail("polygon and rectangle-chain sizes disagree");
    }
    rational_point_copy(&path->vertices[0], &reference_points[start_label]);
    rational_point_copy(
        &path->vertices[path->count - 1],
        &reference_points[end_label]
    );
    rational_box_init(&overlap);
    for (size_t i = 0; i + 1 < panels->count; ++i) {
        if (!rational_box_intersection(
                &overlap,
                &panels->items[i].image,
                &panels->items[i + 1].image
            )) {
            rational_box_clear(&overlap);
            return 0;
        }
        choose_overlap_point(&path->vertices[i + 1], &overlap, attempt, i);
    }
    rational_box_clear(&overlap);

    /* Each segment then lies in its containing rectangle by convexity. */
    for (size_t i = 0; i < panels->count; ++i) {
        if (!rational_box_contains_point(
                &panels->items[i].image,
                &path->vertices[i]
            ) ||
            !rational_box_contains_point(
                &panels->items[i].image,
                &path->vertices[i + 1]
            )) {
            return 0;
        }
    }
    return 1;
}

static int rational_point_equal(
    const rational_point *a,
    const rational_point *b
)
{
    return mpq_equal(a->x, b->x) && mpq_equal(a->y, b->y);
}

static int orientation_sign(
    const rational_point *a,
    const rational_point *b,
    const rational_point *c
)
{
    mpq_t bax;
    mpq_t bay;
    mpq_t cax;
    mpq_t cay;
    mpq_t first;
    mpq_t second;
    mpq_t determinant;
    int sign;

    mpq_inits(bax, bay, cax, cay, first, second, determinant, NULL);
    mpq_sub(bax, b->x, a->x);
    mpq_sub(bay, b->y, a->y);
    mpq_sub(cax, c->x, a->x);
    mpq_sub(cay, c->y, a->y);
    mpq_mul(first, bax, cay);
    mpq_mul(second, bay, cax);
    mpq_sub(determinant, first, second);
    sign = mpq_sgn(determinant);
    mpq_clears(bax, bay, cax, cay, first, second, determinant, NULL);
    return sign;
}

static int rational_between_closed(
    const mpq_t value,
    const mpq_t endpoint_a,
    const mpq_t endpoint_b
)
{
    if (mpq_cmp(endpoint_a, endpoint_b) <= 0) {
        return mpq_cmp(endpoint_a, value) <= 0 &&
               mpq_cmp(value, endpoint_b) <= 0;
    }
    return mpq_cmp(endpoint_b, value) <= 0 &&
           mpq_cmp(value, endpoint_a) <= 0;
}

static int point_on_segment(
    const rational_point *point,
    const rational_point *a,
    const rational_point *b
)
{
    return orientation_sign(a, b, point) == 0 &&
           rational_between_closed(point->x, a->x, b->x) &&
           rational_between_closed(point->y, a->y, b->y);
}

static int segment_bounding_boxes_overlap(
    const rational_point *a,
    const rational_point *b,
    const rational_point *c,
    const rational_point *d
)
{
    const mpq_t *ab_x_min = mpq_cmp(a->x, b->x) <= 0 ? &a->x : &b->x;
    const mpq_t *ab_x_max = mpq_cmp(a->x, b->x) >= 0 ? &a->x : &b->x;
    const mpq_t *ab_y_min = mpq_cmp(a->y, b->y) <= 0 ? &a->y : &b->y;
    const mpq_t *ab_y_max = mpq_cmp(a->y, b->y) >= 0 ? &a->y : &b->y;
    const mpq_t *cd_x_min = mpq_cmp(c->x, d->x) <= 0 ? &c->x : &d->x;
    const mpq_t *cd_x_max = mpq_cmp(c->x, d->x) >= 0 ? &c->x : &d->x;
    const mpq_t *cd_y_min = mpq_cmp(c->y, d->y) <= 0 ? &c->y : &d->y;
    const mpq_t *cd_y_max = mpq_cmp(c->y, d->y) >= 0 ? &c->y : &d->y;

    return mpq_cmp(*ab_x_max, *cd_x_min) >= 0 &&
           mpq_cmp(*cd_x_max, *ab_x_min) >= 0 &&
           mpq_cmp(*ab_y_max, *cd_y_min) >= 0 &&
           mpq_cmp(*cd_y_max, *ab_y_min) >= 0;
}

static int segments_intersect(
    const rational_point *a,
    const rational_point *b,
    const rational_point *c,
    const rational_point *d
)
{
    int o1;
    int o2;
    int o3;
    int o4;

    if (!segment_bounding_boxes_overlap(a, b, c, d)) {
        return 0;
    }
    o1 = orientation_sign(a, b, c);
    o2 = orientation_sign(a, b, d);
    o3 = orientation_sign(c, d, a);
    o4 = orientation_sign(c, d, b);
    if (o1 == 0 && point_on_segment(c, a, b)) {
        return 1;
    }
    if (o2 == 0 && point_on_segment(d, a, b)) {
        return 1;
    }
    if (o3 == 0 && point_on_segment(a, c, d)) {
        return 1;
    }
    if (o4 == 0 && point_on_segment(b, c, d)) {
        return 1;
    }
    return ((o1 < 0 && o2 > 0) || (o1 > 0 && o2 < 0)) &&
           ((o3 < 0 && o4 > 0) || (o3 > 0 && o4 < 0));
}

static int path_is_embedded(
    const rational_path *path,
    const panel_vector *panels
)
{
    size_t segment_count = path->count - 1;

    for (size_t i = 0; i < segment_count; ++i) {
        if (rational_point_equal(&path->vertices[i], &path->vertices[i + 1]) ||
            mpq_equal(path->vertices[i].x, path->vertices[i + 1].x)) {
            return 0;
        }
        if (i + 1 < segment_count &&
            orientation_sign(
                &path->vertices[i],
                &path->vertices[i + 1],
                &path->vertices[i + 2]
            ) == 0) {
            return 0;
        }
    }

    for (size_t i = 0; i < segment_count; ++i) {
        for (size_t j = i + 2; j < segment_count; ++j) {
            /* Segments in disjoint containing rectangles cannot intersect. */
            if (rational_boxes_disjoint(
                    &panels->items[i].image,
                    &panels->items[j].image
                )) {
                continue;
            }
            if (segments_intersect(
                    &path->vertices[i],
                    &path->vertices[i + 1],
                    &path->vertices[j],
                    &path->vertices[j + 1]
                )) {
                return 0;
            }
        }
    }
    return 1;
}

static int gate_strictly_between(
    const mpq_t gate_x,
    const mpq_t x0,
    const mpq_t x1
)
{
    return (mpq_cmp(x0, gate_x) < 0 && mpq_cmp(gate_x, x1) < 0) ||
           (mpq_cmp(x1, gate_x) < 0 && mpq_cmp(gate_x, x0) < 0);
}

static int append_segment_crossing(
    arc_gate_word *word,
    const rational_point *a,
    const rational_point *b,
    const rational_point reference_points[ORBIT_LENGTH],
    int gate,
    int direction
)
{
    mpq_t dx;
    mpq_t dy;
    mpq_t gate_offset;
    mpq_t crossing_y;
    mpq_t work;
    int comparison;

    if (word->length >= ARC_WORD_MAX_LETTERS) {
        return 0;
    }
    mpq_inits(dx, dy, gate_offset, crossing_y, work, NULL);
    mpq_sub(dx, b->x, a->x);
    mpq_sub(dy, b->y, a->y);
    mpq_sub(gate_offset, reference_points[gate].x, a->x);
    mpq_mul(work, gate_offset, dy);
    mpq_div(work, work, dx);
    mpq_add(crossing_y, a->y, work);
    comparison = mpq_cmp(crossing_y, reference_points[gate].y);
    if (comparison == 0) {
        mpq_clears(dx, dy, gate_offset, crossing_y, work, NULL);
        return 0;
    }
    word->letters[word->length].gate = gate;
    word->letters[word->length].side =
        comparison > 0 ? ARC_SIDE_UPPER : ARC_SIDE_LOWER;
    word->letters[word->length].direction = direction;
    ++word->length;
    mpq_clears(dx, dy, gate_offset, crossing_y, work, NULL);
    return 1;
}

/*
 * Compute the Definition 3.11 arc-data, plus crossing directions, in path
 * order.
 */
static int compute_exact_arc_word(
    arc_gate_word *word,
    const rational_path *path,
    const rational_point reference_points[ORBIT_LENGTH],
    int start_label,
    int end_label
)
{
    word->start = start_label;
    word->end = end_label;
    word->length = 0;

    for (size_t vertex = 1; vertex + 1 < path->count; ++vertex) {
        for (int gate = 0; gate < ORBIT_LENGTH; ++gate) {
            if (mpq_equal(path->vertices[vertex].x,
                          reference_points[gate].x)) {
                return 0;
            }
        }
    }

    for (size_t segment = 0; segment + 1 < path->count; ++segment) {
        const rational_point *a = &path->vertices[segment];
        const rational_point *b = &path->vertices[segment + 1];
        int direction = mpq_cmp(a->x, b->x) < 0 ? 1 : -1;

        if (mpq_equal(a->x, b->x)) {
            return 0;
        }
        if (direction > 0) {
            for (int gate = 0; gate < ORBIT_LENGTH; ++gate) {
                if (gate_strictly_between(
                        reference_points[gate].x,
                        a->x,
                        b->x
                    ) &&
                    !append_segment_crossing(
                        word,
                        a,
                        b,
                        reference_points,
                        gate,
                        direction
                    )) {
                    return 0;
                }
            }
        } else {
            for (int gate = ORBIT_LENGTH - 1; gate >= 0; --gate) {
                if (gate_strictly_between(
                        reference_points[gate].x,
                        a->x,
                        b->x
                    ) &&
                    !append_segment_crossing(
                        word,
                        a,
                        b,
                        reference_points,
                        gate,
                        direction
                    )) {
                    return 0;
                }
            }
        }
    }
    return 1;
}

static void print_arc_word(const char *label, const arc_gate_word *word)
{
    printf("  %s: %d ->", label, word->start);
    for (size_t i = 0; i < word->length; ++i) {
        printf(
            " %c%d%c",
            word->letters[i].side == ARC_SIDE_UPPER ? 'U' : 'L',
            word->letters[i].gate,
            word->letters[i].direction > 0 ? '+' : '-'
        );
    }
    printf(" -> %d\n", word->end);
}

static int polygon_certifies_expected_word(
    const rational_path *path,
    const panel_vector *panels,
    const rational_point reference_points[ORBIT_LENGTH],
    int arc,
    int start_label,
    int end_label,
    arc_gate_word *raw_word,
    arc_gate_word *reduced_word
)
{
    if (!path_is_embedded(path, panels)) {
        return 0;
    }
    if (!compute_exact_arc_word(
            raw_word,
            path,
            reference_points,
            start_label,
            end_label
        )) {
        return 0;
    }
    if (!arc_word_matches_legacy(
            raw_word,
            arc_certificate_unreduced_positions[arc],
            arc_certificate_unreduced_heights[arc],
            (size_t) arc_certificate_unreduced_path_lengths[arc]
        )) {
        return 0;
    }

    /*
     * The comparison above checks every crossing. mclass.c instead uses the
     * reduced word stored in arc_certificate_data.c. Endpoint spurs are
     * removed by a relative homotopy near the endpoint; consecutive opposite
     * crossings of the same half-line bound a bigon in a puncture-free
     * vertical cell.
     */
    *reduced_word = *raw_word;
    if (!arc_word_reduce_relative(reduced_word, ORBIT_LENGTH)) {
        return 0;
    }
    return arc_word_matches_legacy(
        reduced_word,
        arc_certificate_reduced_positions[arc],
        arc_certificate_reduced_heights[arc],
        (size_t) arc_certificate_reduced_path_lengths[arc]
    );
}

static void dump_arc_certificate(
    int arc,
    const panel_vector *panels,
    const rational_path *path
)
{
    if (!dump_certificate) {
        return;
    }
    printf("  exact rectangle chain for arc %d:\n", arc);
    for (size_t i = 0; i < panels->count; ++i) {
        gmp_printf("    I_%zu=[%Qd,%Qd], R_%zu=", i,
                   panels->items[i].t0, panels->items[i].t1, i);
        rational_box_fprint(stdout, &panels->items[i].image);
        fputc('\n', stdout);
    }
    printf("  exact vertices of L_%d:\n", arc);
    for (size_t i = 0; i < path->count; ++i) {
        printf("    v_%td=", (ptrdiff_t) i - 1);
        rational_point_fprint(stdout, &path->vertices[i]);
        fputc('\n', stdout);
    }
}

static void certify_one_arc(
    int arc,
    const rational_box punctures[ORBIT_LENGTH],
    const rational_point reference_points[ORBIT_LENGTH],
    const int dynamics[ORBIT_LENGTH]
)
{
    int start_label = dynamics[arc];
    int end_label = dynamics[arc + 1];
    int certified = 0;

    active_arc = arc;
    for (int refinement_round = 0;
         refinement_round < 4 && !certified;
         ++refinement_round) {
        panel_vector panels;
        rational_path path;
        arc_statistics statistics;
        arc_gate_word raw_word;
        arc_gate_word reduced_word;
        int image_width_bits = INITIAL_IMAGE_WIDTH_BITS + refinement_round;

        panel_vector_init(&panels);
        arc_statistics_init(&statistics);
        build_rectangle_chain(
            &panels,
            punctures,
            arc,
            arc + 1,
            start_label,
            end_label,
            image_width_bits,
            &statistics
        );
        rational_path_init(&path, panels.count + 1);
        for (unsigned int attempt = 0;
             attempt < MAX_POLYGON_ATTEMPTS;
             ++attempt) {
            memset(&raw_word, 0, sizeof(raw_word));
            memset(&reduced_word, 0, sizeof(reduced_word));
            if (build_polygon_candidate(
                    &path,
                    &panels,
                    reference_points,
                    start_label,
                    end_label,
                    attempt
                ) &&
                polygon_certifies_expected_word(
                    &path,
                    &panels,
                    reference_points,
                    arc,
                    start_label,
                    end_label,
                    &raw_word,
                    &reduced_word
                )) {
                printf(
                    "Arc %d: B_%d--B_%d maps from B_%d to B_%d; "
                    "%zu certified panels.\n",
                    arc,
                    arc,
                    arc + 1,
                    start_label,
                    end_label,
                    panels.count
                );
                print_arc_word("exact polygon arc-data", &raw_word);
                print_arc_word("reduced mclass word", &reduced_word);
                printf(
                    "  PASS: rational containment, open overlaps, "
                    "embeddedness, transversality, and unreduced/reduced "
                    "row comparisons.\n"
                );
                printf(
                    "  interval work: %lu Pi evaluations, depth %d, ",
                    statistics.panel_evaluations,
                    statistics.maximum_depth
                );
                mpfr_printf(
                    "min(q_lambda)=%.6RDe, "
                    "min(stereo denominator)=%.6RDe\n",
                    statistics.minimum_q_radius_derivative,
                    statistics.minimum_projection_denominator
                );
                dump_arc_certificate(arc, &panels, &path);
                certified = 1;
                break;
            }
        }
        if (!certified) {
            printf(
                "  Arc %d: refining rectangle diameter from 2^(-%d).\n",
                arc,
                image_width_bits
            );
        }
        rational_path_clear(&path);
        arc_statistics_clear(&statistics);
        panel_vector_clear(&panels);
    }
    if (!certified) {
        fail("could not construct an embedded transverse polygon with the "
             "certified relative word");
    }
    active_arc = -1;
}

static void usage(const char *program)
{
    fprintf(
        stderr,
        "usage: %s [--dump-certificate]\n"
        "  --dump-certificate  print every exact t_r, R_r, and polygon vertex\n",
        program
    );
}

static void validate_expected_row(
    const char *kind,
    int arc,
    int length,
    int capacity,
    const int *positions,
    const int *heights
)
{
    if (length < 2 || length > capacity) {
        fail(
            "stored %s arc row %d has invalid length %d",
            kind,
            arc,
            length
        );
    }
    if (positions[0] < 0 ||
        positions[0] >= ORBIT_LENGTH ||
        positions[length - 1] < 0 ||
        positions[length - 1] >= ORBIT_LENGTH ||
        positions[0] == positions[length - 1]) {
        fail("stored %s arc row %d has invalid endpoints", kind, arc);
    }
    if (heights[0] != ARC_CERTIFICATE_ENDPOINT_HEIGHT ||
        heights[length - 1] != ARC_CERTIFICATE_ENDPOINT_HEIGHT) {
        fail("stored %s arc row %d lacks endpoint sentinels", kind, arc);
    }
    for (int i = 1; i + 1 < length; ++i) {
        if (positions[i] < 0 ||
            positions[i] >= ORBIT_LENGTH ||
            (heights[i] != ARC_SIDE_LOWER &&
             heights[i] != ARC_SIDE_UPPER)) {
            fail(
                "stored %s arc row %d has malformed crossing %d",
                kind,
                arc,
                i
            );
        }
    }
}

static void validate_expected_row_shapes(void)
{
    for (int arc = 0; arc < ARC_COUNT; ++arc) {
        int unreduced_length =
            arc_certificate_unreduced_path_lengths[arc];
        int reduced_length = arc_certificate_reduced_path_lengths[arc];

        validate_expected_row(
            "unreduced",
            arc,
            unreduced_length,
            ARC_CERTIFICATE_UNREDUCED_ROW_CAPACITY,
            arc_certificate_unreduced_positions[arc],
            arc_certificate_unreduced_heights[arc]
        );
        validate_expected_row(
            "reduced",
            arc,
            reduced_length,
            ARC_CERTIFICATE_REDUCED_ROW_CAPACITY,
            arc_certificate_reduced_positions[arc],
            arc_certificate_reduced_heights[arc]
        );
        if (
            arc_certificate_unreduced_positions[arc][0] !=
                arc_certificate_reduced_positions[arc][0] ||
            arc_certificate_unreduced_positions[arc][unreduced_length - 1] !=
                arc_certificate_reduced_positions[arc][reduced_length - 1]
        ) {
            fail("stored arc rows %d have inconsistent endpoints", arc);
        }
    }
}

int main(int argc, char **argv)
{
    rational_box punctures[ORBIT_LENGTH];
    rational_point reference_points[ORBIT_LENGTH];
    int tau[ORBIT_LENGTH];
    int inverse_tau[ORBIT_LENGTH];
    int dynamics[ORBIT_LENGTH];

    setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc == 2 && strcmp(argv[1], "--dump-certificate") == 0) {
        dump_certificate = 1;
    } else if (argc != 1) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (ORBIT_LENGTH != K3_ORBIT_LENGTH || ARC_COUNT != 9 ||
        !k3_orbit_data_are_valid()) {
        fail("shared certificate dimensions changed unexpectedly");
    }
    validate_expected_row_shapes();
    for (int i = 0; i < ORBIT_LENGTH; ++i) {
        rational_box_init(&punctures[i]);
        rational_point_init(&reference_points[i]);
    }

    printf("Appendix C.4 path certification (%d-bit MPFR intervals).\n", PRECISION);
    printf("All analytic intervals are outward-rounded; polygonal geometry "
           "is exact over Q.\n");
    printf("In printed words, U/L is the stored above/below arc-data; +/- is "
           "the auxiliary exact crossing direction used for reduction.\n");
    certify_puncture_boxes(
        punctures,
        reference_points,
        tau,
        inverse_tau,
        dynamics
    );
    printf("Step 2: constructing certified rectangle chains and polygons.\n");
    for (int arc = 0; arc < ARC_COUNT; ++arc) {
        certify_one_arc(
            arc,
            punctures,
            reference_points,
            dynamics
        );
    }

    for (int i = 0; i < ORBIT_LENGTH; ++i) {
        rational_point_clear(&reference_points[i]);
        rational_box_clear(&punctures[i]);
    }
    puts("PATH CERTIFICATION COMPLETE.");
    return EXIT_SUCCESS;
}
