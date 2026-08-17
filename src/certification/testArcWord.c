/*
 * Check the stored arc-data rows and the arc-word simplification.  The tests
 * also reject inconsistent crossing directions.
 */

#include "arc_word.h"
#include "arc_certificate_data.h"

#include <assert.h>
#include <stdio.h>

#define ORBIT_LENGTH 10
#define ARC_COUNT (ORBIT_LENGTH - 1)

static const size_t expected_reduced_length[ARC_COUNT] = {
    13, 3, 8, 12, 11, 11, 2, 7, 7
};

static const int expected_reduced_position[ARC_COUNT][13] = {
    {2, 3, 4, 5, 6, 7, 8, 9, 9, 8, 7, 6, 5},
    {5, 6, 7},
    {7, 6, 5, 5, 6, 7, 8, 9},
    {9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 1},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8},
    {8, 9, 9, 8, 7, 6, 5, 4, 3, 3, 4},
    {4, 3},
    {3, 4, 4, 3, 2, 1, 0},
    {0, 1, 2, 3, 4, 5, 6}
};

static const int expected_reduced_height[ARC_COUNT][13] = {
    {-1000, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, -1000},
    {-1000, 1, -1000},
    {-1000, 1, 1, 0, 1, 0, 0, -1000},
    {-1000, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, -1000},
    {-1000, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1000},
    {-1000, 1, 0, 0, 0, 0, 0, 0, 0, 1, -1000},
    {-1000, -1000},
    {-1000, 0, 1, 1, 0, 1, -1000},
    {-1000, 1, 1, 1, 1, 0, -1000}
};

static void test_expected_rows(void)
{
    for (int i = 0; i < ARC_COUNT; i++) {
        arc_gate_word word;

        assert(arc_word_from_legacy(
            &word,
            expected_reduced_position[i],
            expected_reduced_height[i],
            expected_reduced_length[i],
            ORBIT_LENGTH
        ));
        assert(arc_word_matches_legacy(
            &word,
            expected_reduced_position[i],
            expected_reduced_height[i],
            expected_reduced_length[i]
        ));
    }
}

static void test_shared_certificate_rows(void)
{
    for (int i = 0; i < ARC_COUNT; ++i) {
        assert(
            arc_certificate_reduced_path_lengths[i] ==
            (int) expected_reduced_length[i]
        );
        for (size_t j = 0; j < expected_reduced_length[i]; ++j) {
            assert(
                arc_certificate_reduced_positions[i][j] ==
                expected_reduced_position[i][j]
            );
            assert(
                arc_certificate_reduced_heights[i][j] ==
                expected_reduced_height[i][j]
            );
        }
    }
}

/*
 * Check the stored row format.  The certificate program, not this unit test,
 * compares these rows with exact signed words and performs the certified
 * simplification.
 */
static void test_stored_unreduced_rows(void)
{
    for (int i = 0; i < ARC_COUNT; ++i) {
        int length = arc_certificate_unreduced_path_lengths[i];
        int reduced_length = arc_certificate_reduced_path_lengths[i];

        assert(length >= 2);
        assert(length <= ARC_CERTIFICATE_UNREDUCED_ROW_CAPACITY);
        assert(reduced_length >= 2);
        assert(reduced_length <= length);
        assert(
            arc_certificate_unreduced_heights[i][0] ==
            ARC_CERTIFICATE_ENDPOINT_HEIGHT
        );
        assert(
            arc_certificate_unreduced_heights[i][length - 1] ==
            ARC_CERTIFICATE_ENDPOINT_HEIGHT
        );
        assert(
            arc_certificate_unreduced_positions[i][0] ==
            arc_certificate_reduced_positions[i][0]
        );
        assert(
            arc_certificate_unreduced_positions[i][length - 1] ==
            arc_certificate_reduced_positions[i][reduced_length - 1]
        );
        for (int j = 1; j + 1 < length; ++j) {
            assert(arc_certificate_unreduced_positions[i][j] >= 0);
            assert(
                arc_certificate_unreduced_positions[i][j] < ORBIT_LENGTH
            );
            assert(
                arc_certificate_unreduced_heights[i][j] == ARC_SIDE_LOWER ||
                arc_certificate_unreduced_heights[i][j] == ARC_SIDE_UPPER
            );
        }
    }
}

static void test_endpoint_and_backtracking_reduction(void)
{
    arc_gate_word word = {
        .start = 2,
        .end = 5,
        .length = 6,
        .letters = {
            {2, ARC_SIDE_UPPER, 1},
            {3, ARC_SIDE_UPPER, 1},
            {4, ARC_SIDE_LOWER, 1},
            {4, ARC_SIDE_LOWER, -1},
            {4, ARC_SIDE_LOWER, 1},
            {5, ARC_SIDE_UPPER, -1}
        }
    };

    assert(arc_word_reduce_relative(&word, ORBIT_LENGTH));
    assert(word.length == 2);
    assert(word.letters[0].gate == 3);
    assert(word.letters[0].direction == 1);
    assert(word.letters[1].gate == 4);
    assert(word.letters[1].direction == 1);
}

static void test_same_direction_is_not_backtracking(void)
{
    arc_gate_word word = {
        .start = 2,
        .end = 5,
        .length = 4,
        .letters = {
            {3, ARC_SIDE_UPPER, 1},
            {4, ARC_SIDE_LOWER, 1},
            {4, ARC_SIDE_LOWER, 1},
            {5, ARC_SIDE_UPPER, -1}
        }
    };

    assert(!arc_word_reduce_relative(&word, ORBIT_LENGTH));
}

static void test_certified_direction_must_match_comb_path(void)
{
    arc_gate_word word = {
        .start = 2,
        .end = 5,
        .length = 2,
        .letters = {
            {3, ARC_SIDE_UPPER, -1},
            {4, ARC_SIDE_LOWER, 1}
        }
    };

    assert(!arc_word_reduce_relative(&word, ORBIT_LENGTH));
}

static void test_mixed_direction_metadata_is_rejected(void)
{
    arc_gate_word word = {
        .start = 2,
        .end = 5,
        .length = 2,
        .letters = {
            {3, ARC_SIDE_UPPER, 1},
            {4, ARC_SIDE_LOWER, 0}
        }
    };

    assert(!arc_word_reduce_relative(&word, ORBIT_LENGTH));
}

int main(void)
{
    test_shared_certificate_rows();
    test_expected_rows();
    test_stored_unreduced_rows();
    test_endpoint_and_backtracking_reduction();
    test_same_direction_is_not_backtracking();
    test_certified_direction_must_match_comb_path();
    test_mixed_direction_metadata_is_rejected();
    puts("PASS: stored arc rows, normalization, and legacy conversion");
    return 0;
}
