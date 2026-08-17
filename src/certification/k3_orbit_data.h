#ifndef K3_ORBIT_DATA_H
#define K3_ORBIT_DATA_H

#include <stddef.h>

/*
 * The ten entries below are in temporal order.  The decimal strings are
 * intended to be parsed as exact finite decimals by the certification programs.
 */
enum { K3_ORBIT_LENGTH = 10 };

extern const char *const k3_orbit_a_decimal[K3_ORBIT_LENGTH];
extern const char *const k3_orbit_b_decimal[K3_ORBIT_LENGTH];

/*
 * Legacy chart codes used by periodicExact.c and the manuscript tables:
 *
 *     1 = Psi_3^+    2 = Psi_2^+    3 = Psi_1^+
 *     4 = Psi_3^-    5 = Psi_2^-    6 = Psi_1^-
 *
 * These codes are retained as the canonical stored representation because
 * derivativesInterval.c and verifyPeriodicOrbit.c already use them.
 */
extern const int k3_orbit_chart_code[K3_ORBIT_LENGTH];

/* Return nonzero exactly when the supplied value belongs to the domain. */
int k3_orbit_index_is_valid(size_t index);
int k3_chart_code_is_valid(int code);
int k3_chart_number_is_valid(int number);
int k3_chart_sign_is_valid(int sign);

/*
 * Convert without terminating the caller.  On invalid input these functions
 * return zero and leave every non-NULL output unchanged.  A NULL output is
 * permitted when the caller only wants to validate the input.
 */
int k3_chart_decode(int code, int *number, int *sign);
int k3_chart_encode(int number, int sign, int *code);

/* Decode the canonical chart belonging to a temporal orbit index. */
int k3_orbit_chart_at(size_t index, int *number, int *sign);

/* Return a short printable label such as "Psi2-", or NULL for a bad code. */
const char *k3_chart_label(int code);

/* Check the internal table and every stored chart-code round trip. */
int k3_orbit_data_are_valid(void);

#endif
