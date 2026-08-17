#include "arc_certificate_data.h"

/*
 * Literal arc-data of the rational polygonal paths checked by verifyArcs.c.
 * The first and last entries in each row are endpoint labels; all intervening
 * entries are crossings in path order, before simplification.
 */
const int
arc_certificate_unreduced_path_lengths[ARC_CERTIFICATE_ARC_COUNT] = {
    13, 4, 23, 15, 17, 14, 4, 22, 18
};

const int
arc_certificate_unreduced_positions[ARC_CERTIFICATE_ARC_COUNT]
                                    [ARC_CERTIFICATE_UNREDUCED_ROW_CAPACITY] = {
    { 2, 3, 4, 5, 6, 7, 8, 9, 9, 8, 7, 6, 5 },
    { 5, 5, 6, 7 },
    { 7, 6, 5, 4, 3, 2, 1, 0, 0, 1, 1, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9 },
    { 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 1, 1, 1 },
    { 1, 2, 3, 4, 4, 3, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8 },
    { 8, 9, 9, 8, 7, 6, 5, 4, 3, 2, 2, 3, 4, 4 },
    { 4, 4, 3, 3 },
    { 3, 4, 5, 6, 7, 8, 9, 9, 8, 8, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0 },
    { 0, 0, 1, 2, 3, 4, 4, 3, 2, 1, 0, 0, 1, 2, 3, 4, 5, 6 }
};

const int
arc_certificate_unreduced_heights[ARC_CERTIFICATE_ARC_COUNT]
                                  [ARC_CERTIFICATE_UNREDUCED_ROW_CAPACITY] = {
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      0, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT }
};

/*
 * Reduced rows used by mclass.c.  They are obtained from the
 * unreduced rows above by the relative simplification in arc_word.c.
 */
const int
arc_certificate_reduced_path_lengths[ARC_CERTIFICATE_ARC_COUNT] = {
    13, 3, 8, 12, 11, 11, 2, 7, 7
};

const int
arc_certificate_reduced_positions[ARC_CERTIFICATE_ARC_COUNT]
                                  [ARC_CERTIFICATE_REDUCED_ROW_CAPACITY] = {
    { 2, 3, 4, 5, 6, 7, 8, 9, 9, 8, 7, 6, 5 },
    { 5, 6, 7 },
    { 7, 6, 5, 5, 6, 7, 8, 9 },
    { 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 1 },
    { 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8 },
    { 8, 9, 9, 8, 7, 6, 5, 4, 3, 3, 4 },
    { 4, 3 },
    { 3, 4, 4, 3, 2, 1, 0 },
    { 0, 1, 2, 3, 4, 5, 6 }
};

const int
arc_certificate_reduced_heights[ARC_CERTIFICATE_ARC_COUNT]
                                [ARC_CERTIFICATE_REDUCED_ROW_CAPACITY] = {
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 1, 0, 1, 0, 0,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      0, 0, 1, 0, 1, 1, 1, 1, 1, 0,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 0, 0, 0, 0, 0, 0, 0, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      0, 1, 1, 0, 1,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT },
    { ARC_CERTIFICATE_ENDPOINT_HEIGHT,
      1, 1, 1, 1, 0,
      ARC_CERTIFICATE_ENDPOINT_HEIGHT }
};
