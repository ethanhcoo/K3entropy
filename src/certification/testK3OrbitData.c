/* Check the shared pseudo-orbit decimals and chart codes. */

#include "k3_orbit_data.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition)                                                     \
    do {                                                                     \
        if (!(condition)) {                                                  \
            fprintf(                                                         \
                stderr,                                                      \
                "testK3OrbitData: check failed at line %d: %s\n",           \
                __LINE__,                                                    \
                #condition                                                   \
            );                                                               \
            return EXIT_FAILURE;                                             \
        }                                                                    \
    } while (0)

static int check_canonical_orbit(void)
{
    static const char *const expected_a[K3_ORBIT_LENGTH] = {
        "1.041643093944314148360673792017",
        "-0.439586738044637984442175311821",
        "1.111402054756051352317454938205",
        "-0.328869789067645570794396144391",
        "1.488818954806569814700993326668",
        "0.533829900932504729554816817729",
        "-1.004464450964796276608907444033",
        "-0.867950394543647373540816310310",
        "0.926435350008842162121508383319",
        "0.555943953085459715621476373770"
    };
    static const char *const expected_b[K3_ORBIT_LENGTH] = {
        "1.726895448754858426328854724474",
        "0.555943953085459715621476373770",
        "-0.926435350008842162121508383319",
        "0.867950394543647373540816310310",
        "1.004464450964796276608907444033",
        "-0.533829900932504729554816817729",
        "-1.488818954806569814700993326668",
        "-0.328869789067645570794396144391",
        "-1.111402054756051352317454938205",
        "0.439586738044637984442175311821"
    };
    static const int expected_code[K3_ORBIT_LENGTH] = {
        6, 4, 5, 6, 5, 5, 5, 1, 5, 3
    };
    static const int expected_number[K3_ORBIT_LENGTH] = {
        1, 3, 2, 1, 2, 2, 2, 3, 2, 1
    };
    static const int expected_sign[K3_ORBIT_LENGTH] = {
        -1, -1, -1, -1, -1, -1, -1, 1, -1, 1
    };
    size_t index;

    CHECK(K3_ORBIT_LENGTH == 10);
    CHECK(k3_orbit_data_are_valid());

    for (index = 0; index < (size_t)K3_ORBIT_LENGTH; ++index) {
        int number = 0;
        int sign = 0;

        CHECK(k3_orbit_index_is_valid(index));
        CHECK(strcmp(k3_orbit_a_decimal[index], expected_a[index]) == 0);
        CHECK(strcmp(k3_orbit_b_decimal[index], expected_b[index]) == 0);
        CHECK(k3_orbit_chart_code[index] == expected_code[index]);
        CHECK(k3_orbit_chart_at(index, &number, &sign));
        CHECK(number == expected_number[index]);
        CHECK(sign == expected_sign[index]);
    }
    CHECK(!k3_orbit_index_is_valid((size_t)K3_ORBIT_LENGTH));
    return EXIT_SUCCESS;
}

static int check_chart_conversions(void)
{
    static const int expected_number[6] = {3, 2, 1, 3, 2, 1};
    static const int expected_sign[6] = {1, 1, 1, -1, -1, -1};
    static const char *const expected_label[6] = {
        "Psi3+", "Psi2+", "Psi1+", "Psi3-", "Psi2-", "Psi1-"
    };
    int code;

    for (code = 1; code <= 6; ++code) {
        int number = 0;
        int sign = 0;
        int round_trip_code = 0;

        CHECK(k3_chart_code_is_valid(code));
        CHECK(k3_chart_decode(code, &number, &sign));
        CHECK(number == expected_number[code - 1]);
        CHECK(sign == expected_sign[code - 1]);
        CHECK(k3_chart_encode(number, sign, &round_trip_code));
        CHECK(round_trip_code == code);
        CHECK(strcmp(k3_chart_label(code), expected_label[code - 1]) == 0);

        /* NULL outputs are valid for callers performing validation only. */
        CHECK(k3_chart_decode(code, NULL, NULL));
        CHECK(k3_chart_encode(number, sign, NULL));
    }
    return EXIT_SUCCESS;
}

static int check_rejected_inputs(void)
{
    static const int bad_codes[] = {-1, 0, 7, 100};
    static const int bad_numbers[] = {-1, 0, 4, 100};
    static const int bad_signs[] = {-2, 0, 2, 100};
    size_t index;

    for (index = 0; index < sizeof(bad_codes) / sizeof(bad_codes[0]); ++index) {
        int number = 71;
        int sign = 73;

        CHECK(!k3_chart_code_is_valid(bad_codes[index]));
        CHECK(!k3_chart_decode(bad_codes[index], &number, &sign));
        CHECK(number == 71);
        CHECK(sign == 73);
        CHECK(k3_chart_label(bad_codes[index]) == NULL);
    }

    for (index = 0;
         index < sizeof(bad_numbers) / sizeof(bad_numbers[0]);
         ++index) {
        int code = 79;

        CHECK(!k3_chart_number_is_valid(bad_numbers[index]));
        CHECK(!k3_chart_encode(bad_numbers[index], 1, &code));
        CHECK(code == 79);
    }

    for (index = 0; index < sizeof(bad_signs) / sizeof(bad_signs[0]); ++index) {
        int code = 83;

        CHECK(!k3_chart_sign_is_valid(bad_signs[index]));
        CHECK(!k3_chart_encode(1, bad_signs[index], &code));
        CHECK(code == 83);
    }

    {
        int number = 89;
        int sign = 97;

        CHECK(!k3_orbit_chart_at((size_t)K3_ORBIT_LENGTH, &number, &sign));
        CHECK(number == 89);
        CHECK(sign == 97);
    }
    return EXIT_SUCCESS;
}

int main(void)
{
    CHECK(check_canonical_orbit() == EXIT_SUCCESS);
    CHECK(check_chart_conversions() == EXIT_SUCCESS);
    CHECK(check_rejected_inputs() == EXIT_SUCCESS);

    puts("k3 orbit-data tests passed");
    return EXIT_SUCCESS;
}
