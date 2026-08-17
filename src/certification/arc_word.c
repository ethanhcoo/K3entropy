/*
 * Simplification and validation of the arc-data used by verifyArcs.c.  This
 * removes endpoint spurs and cancelling crossings, checks crossing directions,
 * and provides comparison with stored position/height rows.
 */

#include "arc_word.h"

#include <string.h>

static int valid_gate_side(
    const arc_gate_letter *letter,
    int puncture_count
)
{
    return
        letter->gate >= 0 &&
        letter->gate < puncture_count &&
        (letter->side == ARC_SIDE_LOWER ||
         letter->side == ARC_SIDE_UPPER) &&
        (letter->direction == -1 ||
         letter->direction == 0 ||
         letter->direction == 1);
}

static void delete_letter(arc_gate_word *word, size_t index)
{
    if (index + 1 < word->length) {
        memmove(
            &word->letters[index],
            &word->letters[index + 1],
            (word->length - index - 1) * sizeof(word->letters[0])
        );
    }
    word->length--;
}

static void delete_pair(arc_gate_word *word, size_t index)
{
    if (index + 2 < word->length) {
        memmove(
            &word->letters[index],
            &word->letters[index + 2],
            (word->length - index - 2) * sizeof(word->letters[0])
        );
    }
    word->length -= 2;
}

int arc_word_reduce_relative(arc_gate_word *word, int puncture_count)
{
    int changed;
    int has_certified_direction = 0;
    int has_legacy_direction = 0;

    if (
        word == NULL ||
        puncture_count < 2 ||
        word->start < 0 ||
        word->start >= puncture_count ||
        word->end < 0 ||
        word->end >= puncture_count ||
        word->start == word->end ||
        word->length > ARC_WORD_MAX_LETTERS
    ) {
        return 0;
    }

    for (size_t i = 0; i < word->length; i++) {
        if (!valid_gate_side(&word->letters[i], puncture_count)) {
            return 0;
        }
        if (word->letters[i].direction == 0) {
            has_legacy_direction = 1;
        } else {
            has_certified_direction = 1;
        }
    }

    /*
     * Either every direction is +/-1, or every direction is zero and will be
     * inferred after simplification.  Mixing the two would make cancellation
     * depend on missing direction data.
     */
    if (has_certified_direction && has_legacy_direction) {
        return 0;
    }

    do {
        changed = 0;

        if (
            word->length > 0 &&
            word->letters[0].gate == word->start
        ) {
            delete_letter(word, 0);
            changed = 1;
            continue;
        }

        if (
            word->length > 0 &&
            word->letters[word->length - 1].gate == word->end
        ) {
            delete_letter(word, word->length - 1);
            changed = 1;
            continue;
        }

        for (size_t i = 0; i + 1 < word->length; i++) {
            if (
                word->letters[i].gate == word->letters[i + 1].gate &&
                word->letters[i].side == word->letters[i + 1].side
            ) {
                /*
                 * Two consecutive crossings of the same ray can cancel only
                 * when their directions are opposite.  A zero direction comes
                 * from a position/height row and is inferred after
                 * simplification.
                 */
                if (
                    word->letters[i].direction != 0 &&
                    word->letters[i + 1].direction != 0 &&
                    word->letters[i].direction ==
                        word->letters[i + 1].direction
                ) {
                    return 0;
                }
                delete_pair(word, i);
                changed = 1;
                break;
            }
        }
    } while (changed);

    return arc_word_infer_directions(word, puncture_count);
}

int arc_word_infer_directions(arc_gate_word *word, int puncture_count)
{
    int strip;

    if (
        word == NULL ||
        puncture_count < 2 ||
        word->start < 0 ||
        word->start >= puncture_count ||
        word->end < 0 ||
        word->end >= puncture_count ||
        word->start == word->end ||
        word->length > ARC_WORD_MAX_LETTERS
    ) {
        return 0;
    }

    if (word->length == 0) {
        return
            word->start + 1 == word->end ||
            word->end + 1 == word->start;
    }

    if (word->letters[0].gate < word->start) {
        strip = word->start;
    } else if (word->letters[0].gate > word->start) {
        strip = word->start + 1;
    } else {
        /* Endpoint-simplified words never start on the start barrier. */
        return 0;
    }

    for (size_t i = 0; i < word->length; i++) {
        arc_gate_letter *letter = &word->letters[i];
        int inferred_direction;

        if (!valid_gate_side(letter, puncture_count)) {
            return 0;
        }

        if (strip == letter->gate) {
            inferred_direction = 1;
            strip = letter->gate + 1;
        } else if (strip == letter->gate + 1) {
            inferred_direction = -1;
            strip = letter->gate;
        } else {
            return 0;
        }
        if (
            letter->direction != 0 &&
            letter->direction != inferred_direction
        ) {
            return 0;
        }
        letter->direction = inferred_direction;
    }

    /*
     * The final strip must be one of the two strips incident to the ending
     * puncture.  A terminal crossing of the end barrier would have been
     * removed by relative simplification.
     */
    return strip == word->end || strip == word->end + 1;
}

int arc_word_from_legacy(
    arc_gate_word *word,
    const int *position,
    const int *height,
    size_t legacy_length,
    int puncture_count
)
{
    if (
        word == NULL ||
        position == NULL ||
        height == NULL ||
        legacy_length < 2 ||
        legacy_length - 2 > ARC_WORD_MAX_LETTERS
    ) {
        return 0;
    }

    word->start = position[0];
    word->end = position[legacy_length - 1];
    word->length = legacy_length - 2;

    for (size_t i = 0; i < word->length; i++) {
        word->letters[i].gate = position[i + 1];
        word->letters[i].side = height[i + 1];
        word->letters[i].direction = 0;
    }

    return arc_word_reduce_relative(word, puncture_count);
}

int arc_word_matches_legacy(
    const arc_gate_word *word,
    const int *position,
    const int *height,
    size_t legacy_length
)
{
    if (
        word == NULL ||
        position == NULL ||
        height == NULL ||
        legacy_length != word->length + 2 ||
        position[0] != word->start ||
        position[legacy_length - 1] != word->end
    ) {
        return 0;
    }

    for (size_t i = 0; i < word->length; i++) {
        if (
            position[i + 1] != word->letters[i].gate ||
            height[i + 1] != word->letters[i].side
        ) {
            return 0;
        }
    }

    return 1;
}
