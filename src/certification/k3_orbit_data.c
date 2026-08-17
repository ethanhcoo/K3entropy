#include "k3_orbit_data.h"

/*
 * Table 1 pseudo-orbit centers, in temporal order.  Keeping these strings in
 * one translation unit prevents the three rigorous verifiers from silently
 * drifting apart.  They are data, not binary floating-point constants.
 */
const char *const k3_orbit_a_decimal[K3_ORBIT_LENGTH] = {
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

const char *const k3_orbit_b_decimal[K3_ORBIT_LENGTH] = {
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

const int k3_orbit_chart_code[K3_ORBIT_LENGTH] = {
    6, 4, 5, 6, 5, 5, 5, 1, 5, 3
};

int k3_orbit_index_is_valid(size_t index)
{
    return index < (size_t)K3_ORBIT_LENGTH;
}

int k3_chart_code_is_valid(int code)
{
    return code >= 1 && code <= 6;
}

int k3_chart_number_is_valid(int number)
{
    return number >= 1 && number <= 3;
}

int k3_chart_sign_is_valid(int sign)
{
    return sign == -1 || sign == 1;
}

int k3_chart_decode(int code, int *number, int *sign)
{
    int decoded_number;
    int decoded_sign;

    if (!k3_chart_code_is_valid(code)) {
        return 0;
    }

    if (code <= 3) {
        decoded_number = 4 - code;
        decoded_sign = 1;
    } else {
        decoded_number = 7 - code;
        decoded_sign = -1;
    }

    if (number != NULL) {
        *number = decoded_number;
    }
    if (sign != NULL) {
        *sign = decoded_sign;
    }
    return 1;
}

int k3_chart_encode(int number, int sign, int *code)
{
    int encoded_code;

    if (!k3_chart_number_is_valid(number) ||
        !k3_chart_sign_is_valid(sign)) {
        return 0;
    }

    encoded_code = sign > 0 ? 4 - number : 7 - number;
    if (code != NULL) {
        *code = encoded_code;
    }
    return 1;
}

int k3_orbit_chart_at(size_t index, int *number, int *sign)
{
    if (!k3_orbit_index_is_valid(index)) {
        return 0;
    }
    return k3_chart_decode(k3_orbit_chart_code[index], number, sign);
}

const char *k3_chart_label(int code)
{
    static const char *const labels[6] = {
        "Psi3+", "Psi2+", "Psi1+", "Psi3-", "Psi2-", "Psi1-"
    };

    if (!k3_chart_code_is_valid(code)) {
        return NULL;
    }
    return labels[code - 1];
}

int k3_orbit_data_are_valid(void)
{
    size_t index;

    for (index = 0; index < (size_t)K3_ORBIT_LENGTH; ++index) {
        int number;
        int sign;
        int round_trip_code;

        if (k3_orbit_a_decimal[index] == NULL ||
            k3_orbit_a_decimal[index][0] == '\0' ||
            k3_orbit_b_decimal[index] == NULL ||
            k3_orbit_b_decimal[index][0] == '\0' ||
            !k3_orbit_chart_at(index, &number, &sign) ||
            !k3_chart_encode(number, sign, &round_trip_code) ||
            round_trip_code != k3_orbit_chart_code[index]) {
            return 0;
        }
    }
    return 1;
}
