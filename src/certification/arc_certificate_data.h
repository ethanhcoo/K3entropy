#ifndef ARC_CERTIFICATE_DATA_H
#define ARC_CERTIFICATE_DATA_H

enum {
    ARC_CERTIFICATE_ORBIT_LENGTH = 10,
    ARC_CERTIFICATE_ARC_COUNT = 9,
    ARC_CERTIFICATE_UNREDUCED_ROW_CAPACITY = 23,
    ARC_CERTIFICATE_REDUCED_ROW_CAPACITY = 13,
    ARC_CERTIFICATE_ENDPOINT_HEIGHT = -1000
};

/*
 * Both row formats include the starting and ending puncture labels as the
 * first and last positions.  Their height entries use the sentinel
 * ARC_CERTIFICATE_ENDPOINT_HEIGHT; every interior height is 0 (below) or
 * 1 (above).  Row k in each array belongs to arc k.
 *
 * The unreduced rows record the arc-data of the rational polygonal paths
 * constructed and checked by verifyArcs.c.  The reduced rows are obtained by
 * deleting endpoint spurs and adjacent inverse crossings.  Only the reduced
 * rows are input to mclass.c.
 *
 * Crossing direction is not stored because it is not part of the manuscript's
 * definition of arc-data.  verifyArcs.c determines each direction exactly
 * from the rational polygon before simplifying the row.
 */
extern const int
    arc_certificate_unreduced_path_lengths[ARC_CERTIFICATE_ARC_COUNT];
extern const int
    arc_certificate_unreduced_positions[ARC_CERTIFICATE_ARC_COUNT]
                                        [ARC_CERTIFICATE_UNREDUCED_ROW_CAPACITY];
extern const int
    arc_certificate_unreduced_heights[ARC_CERTIFICATE_ARC_COUNT]
                                      [ARC_CERTIFICATE_UNREDUCED_ROW_CAPACITY];

extern const int
    arc_certificate_reduced_path_lengths[ARC_CERTIFICATE_ARC_COUNT];
extern const int
    arc_certificate_reduced_positions[ARC_CERTIFICATE_ARC_COUNT]
                                      [ARC_CERTIFICATE_REDUCED_ROW_CAPACITY];
extern const int
    arc_certificate_reduced_heights[ARC_CERTIFICATE_ARC_COUNT]
                                    [ARC_CERTIFICATE_REDUCED_ROW_CAPACITY];

#endif
