#ifndef ARC_WORD_H
#define ARC_WORD_H

#include <stddef.h>

#define ARC_WORD_MAX_LETTERS 1024

enum {
    ARC_SIDE_LOWER = 0,
    ARC_SIDE_UPPER = 1
};

typedef struct {
    int gate;
    int side;
    int direction;
} arc_gate_letter;

typedef struct {
    int start;
    int end;
    size_t length;
    arc_gate_letter letters[ARC_WORD_MAX_LETTERS];
} arc_gate_word;

/*
 * Reduce a gate word relative to its endpoint punctures.
 *
 * In addition to ordinary cancellation of adjacent inverse crossings of the
 * same gate on the same side, crossings of the start barrier at the beginning
 * and crossings of the end barrier at the end are endpoint spurs.  They are
 * deleted iteratively.  When directions are supplied, a cancellable pair must
 * have opposite directions.  Direction zero means "not yet inferred", as in a
 * legacy position/height row.  This is the topological normalization
 * implemented by simplify_arrays() in arcs.c and simplify_path() in mclass.c.
 *
 * Returns 1 on success and 0 if the input is malformed.
 */
int arc_word_reduce_relative(arc_gate_word *word, int puncture_count);

/*
 * Infer crossing directions from an endpoint-normalized gate/side sequence.
 * The barriers are numbered from left to right.  direction is set to +1 for
 * left-to-right and -1 for right-to-left.
 *
 * Returns 1 exactly when the sequence is a valid edge path in the comb graph.
 */
int arc_word_infer_directions(arc_gate_word *word, int puncture_count);

/*
 * Convert the legacy position/height encoding to a reduced signed gate word.
 * The first and last positions are endpoint labels; their heights are ignored.
 */
int arc_word_from_legacy(
    arc_gate_word *word,
    const int *position,
    const int *height,
    size_t legacy_length,
    int puncture_count
);

/*
 * Compare a gate word with a row in the legacy position/height encoding.
 * The row may be unreduced or reduced; this function performs no
 * normalization.  Crossing direction and the endpoint height sentinel are
 * deliberately not part of the comparison because neither belongs to the
 * stored position/height encoding.
 */
int arc_word_matches_legacy(
    const arc_gate_word *word,
    const int *position,
    const int *height,
    size_t legacy_length
);

#endif
