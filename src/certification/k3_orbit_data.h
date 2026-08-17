#ifndef K3_ORBIT_DATA_H
#define K3_ORBIT_DATA_H

#include <stddef.h>

/*
 * The ten entries below are in temporal order.  The certification programs
 * parse the coordinate strings as exact finite decimals.
 */
enum { K3_ORBIT_LENGTH = 10 };

extern const char *const k3_orbit_a_decimal[K3_ORBIT_LENGTH];
extern const char *const k3_orbit_b_decimal[K3_ORBIT_LENGTH];

/*
 * Chart encoding for the chart choices displayed in Table 1:
 *
 *     1 = Psi_3^+    2 = Psi_2^+    3 = Psi_1^+
 *     4 = Psi_3^-    5 = Psi_2^-    6 = Psi_1^-
 *
 * These codes are the shared representation used by the certification
 * programs.
 */
extern const int k3_orbit_chart_code[K3_ORBIT_LENGTH];

/* Return nonzero for a valid orbit index, chart code, number, or sign. */
int k3_orbit_index_is_valid(size_t index);
int k3_chart_code_is_valid(int code);
int k3_chart_number_is_valid(int number);
int k3_chart_sign_is_valid(int sign);

/*
 * Encode or decode a chart.  Invalid input returns zero and leaves every
 * non-NULL output unchanged.  Outputs may be NULL when only validation is
 * needed.
 */
int k3_chart_decode(int code, int *number, int *sign);
int k3_chart_encode(int number, int sign, int *code);

/* Decode the stored chart belonging to a temporal orbit index. */
int k3_orbit_chart_at(size_t index, int *number, int *sign);

/* Return a short label such as "Psi2-", or NULL for an invalid code. */
const char *k3_chart_label(int code);

/* Check the internal table and every stored chart-code round trip. */
int k3_orbit_data_are_valid(void);

#endif
