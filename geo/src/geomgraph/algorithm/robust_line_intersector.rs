use super::{Intersection, LineIntersector};
use crate::algorithm::kernels::{Kernel, Orientation, RobustKernel};
use crate::contains::Contains;
use crate::intersects::Intersects;
use crate::num_traits::Zero;
use crate::{Coordinate, GeoFloat, Rect};

// JTS: /**
// JTS:  * A robust version of {@link LineIntersector}.
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class RobustLineIntersector
// JTS:     extends LineIntersector
// JTS: {
/// A robust version of [LineIntersector](traits.LineIntersector).
#[derive(Clone)]
pub(crate) struct RobustLineIntersector<F: GeoFloat> {
    // TODO: JTS captures some state in the LineIntersector. I'm not sure if it's helpful. Roughly it
    // seems to be mid-computation state and result state. Perhaps that could be better modeled and
    // we could leave the LineIntersector less mutable.
    result: Option<Intersection>,
    // TODO: make this [geo_types::Line; 2]? matching JTS for now
    is_proper: bool,
    input_lines: [[Coordinate<F>; 2]; 2],
    int_pt: [Coordinate<F>; 2],
}

impl<F: GeoFloat> RobustLineIntersector<F> {
    // JTS:   public RobustLineIntersector() {
    // JTS:   }
    pub fn new() -> RobustLineIntersector<F> {
        RobustLineIntersector {
            result: None,
            is_proper: false,
            input_lines: [[Coordinate::zero(); 2]; 2],
            int_pt: [Coordinate::zero(); 2],
        }
    }
}

impl<F: GeoFloat> LineIntersector<F> for RobustLineIntersector<F> {
    fn intersection(&self, intersection_index: usize) -> Coordinate<F> {
        self.int_pt[intersection_index]
    }

    fn result(&self) -> Intersection {
        match self.result {
            None => {
                // JTS initializes result to `NoIntersection`, but it seems like in practice it
                // should never be returned before being initialized.
                debug_assert!(false, "result was unexpectedly None");
                Intersection::NoIntersection
            }
            Some(i) => i,
        }
    }

    fn set_result(&mut self, result: Intersection) {
        self.result = Some(result);
    }

    fn input_lines(&self) -> [[Coordinate<F>; 2]; 2] {
        self.input_lines
    }

    fn set_input_lines(&mut self, line: usize, start_or_end: usize, coord: Coordinate<F>) {
        assert!(line <= 2);
        assert!(start_or_end <= 2);
        self.input_lines[line][start_or_end] = coord;
    }

    fn is_proper(&self) -> bool {
        self.has_intersection() && self.is_proper
    }

    // JTS:
    // JTS:   public void computeIntersection(Coordinate p, Coordinate p1, Coordinate p2) {
    // JTS:     isProper = false;
    // JTS:     // do between check first, since it is faster than the orientation test
    // JTS:     if (Envelope.intersects(p1, p2, p)) {
    // JTS:       if ((Orientation.index(p1, p2, p) == 0)
    // JTS:           && (Orientation.index(p2, p1, p) == 0)) {
    // JTS:         isProper = true;
    // JTS:         if (p.equals(p1) || p.equals(p2)) {
    // JTS:           isProper = false;
    // JTS:         }
    // JTS:         result = POINT_INTERSECTION;
    // JTS:         return;
    // JTS:       }
    // JTS:     }
    // JTS:     result = NO_INTERSECTION;
    // JTS:   }
    // JTS:
    // JTS:   protected int computeIntersect(
    // JTS:                 Coordinate p1, Coordinate p2,
    // JTS:                 Coordinate q1, Coordinate q2  ) {
    fn compute_intersect(
        &mut self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Intersection {
        // JTS:     isProper = false;
        self.is_proper = false;
        // JTS:     // first try a fast test to see if the envelopes of the lines intersect
        // JTS:     if (! Envelope.intersects(p1, p2, q1, q2))
        // JTS:       return NO_INTERSECTION;
        // first try a fast test to see if the envelopes of the lines intersect
        if !Rect::new(p1, p2).intersects(&Rect::new(q1, q2)) {
            return Intersection::NoIntersection;
        }

        // JTS:     // for each endpoint, compute which side of the other segment it lies
        // JTS:     // if both endpoints lie on the same side of the other segment,
        // JTS:     // the segments do not intersect
        // JTS:     int Pq1 = Orientation.index(p1, p2, q1);
        // JTS:     int Pq2 = Orientation.index(p1, p2, q2);
        // JTS:
        // JTS:     if ((Pq1>0 && Pq2>0) || (Pq1<0 && Pq2<0)) {
        // JTS:       return NO_INTERSECTION;
        // JTS:     }
        // for each endpoint, compute which side of the other segment it lies
        // if both endpoints lie on the same side of the other segment,
        // the segments do not intersect
        let p_q1 = RobustKernel::orient2d(p1, p2, q1);
        let p_q2 = RobustKernel::orient2d(p1, p2, q2);
        match (p_q1, p_q2) {
            (Orientation::Clockwise, Orientation::Clockwise)
            | (Orientation::CounterClockwise, Orientation::CounterClockwise) => {
                return Intersection::NoIntersection
            }
            _ => (),
        }

        // JTS:     int Qp1 = Orientation.index(q1, q2, p1);
        // JTS:     int Qp2 = Orientation.index(q1, q2, p2);
        // JTS:
        // JTS:     if ((Qp1>0 && Qp2>0) || (Qp1<0 && Qp2<0)) {
        // JTS:         return NO_INTERSECTION;
        // JTS:     }
        let q_p1 = RobustKernel::orient2d(q1, q2, p1);
        let q_p2 = RobustKernel::orient2d(q1, q2, p2);
        match (q_p1, q_p2) {
            (Orientation::Clockwise, Orientation::Clockwise)
            | (Orientation::CounterClockwise, Orientation::CounterClockwise) => {
                return Intersection::NoIntersection
            }
            _ => (),
        }

        // JTS:     boolean collinear = Pq1 == 0
        // JTS:          && Pq2 == 0
        // JTS:          && Qp1 == 0
        // JTS:          && Qp2 == 0;
        // JTS:     if (collinear) {
        // JTS:       return computeCollinearIntersection(p1, p2, q1, q2);
        // JTS:     }
        if let (
            Orientation::Collinear,
            Orientation::Collinear,
            Orientation::Collinear,
            Orientation::Collinear,
        ) = (p_q1, p_q2, q_p1, q_p2)
        {
            return self.compute_collinear_intersection(p1, p2, q1, q2);
        }
        // JTS:     /**
        // JTS:      * At this point we know that there is a single intersection point
        // JTS:      * (since the lines are not collinear).
        // JTS:      */
        // JTS:
        // JTS:     /**
        // JTS:      *  Check if the intersection is an endpoint. If it is, copy the endpoint as
        // JTS:      *  the intersection point. Copying the point rather than computing it
        // JTS:      *  ensures the point has the exact value, which is important for
        // JTS:      *  robustness. It is sufficient to simply check for an endpoint which is on
        // JTS:      *  the other line, since at this point we know that the inputLines must
        // JTS:      *  intersect.
        // JTS:      */
        // At this point we know that there is a single intersection point (since the lines are not
        // collinear).
        //
        // Check if the intersection is an endpoint. If it is, copy the endpoint as the
        // intersection point. Copying the point rather than computing it ensures the point has the
        // exact value, which is important for robustness. It is sufficient to simply check for an
        // endpoint which is on the other line, since at this point we know that the inputLines
        // must intersect.
        // JTS:     if (Pq1 == 0 || Pq2 == 0 || Qp1 == 0 || Qp2 == 0) {
        // JTS:       isProper = false;
        if p_q1 == Orientation::Collinear
            || p_q2 == Orientation::Collinear
            || q_p1 == Orientation::Collinear
            || q_p2 == Orientation::Collinear
        {
            self.is_proper = false;
            // JTS:       /**
            // JTS:        * Check for two equal endpoints.
            // JTS:        * This is done explicitly rather than by the orientation tests
            // JTS:        * below in order to improve robustness.
            // JTS:        *
            // JTS:        * [An example where the orientation tests fail to be consistent is
            // JTS:        * the following (where the true intersection is at the shared endpoint
            // JTS:        * POINT (19.850257749638203 46.29709338043669)
            // JTS:        *
            // JTS:        * LINESTRING ( 19.850257749638203 46.29709338043669, 20.31970698357233 46.76654261437082 )
            // JTS:        * and
            // JTS:        * LINESTRING ( -48.51001596420236 -22.063180333403878, 19.850257749638203 46.29709338043669 )
            // JTS:        *
            // JTS:        * which used to produce the INCORRECT result: (20.31970698357233, 46.76654261437082, NaN)
            // JTS:        *
            // JTS:        */
            // Check for two equal endpoints.
            // This is done explicitly rather than by the orientation tests below in order to improve
            // robustness.
            //
            // [An example where the orientation tests fail to be consistent is the following (where
            // the true intersection is at the shared endpoint
            // POINT (19.850257749638203 46.29709338043669)
            //
            // LINESTRING ( 19.850257749638203 46.29709338043669, 20.31970698357233 46.76654261437082 )
            // and
            // LINESTRING ( -48.51001596420236 -22.063180333403878, 19.850257749638203 46.29709338043669 )
            //
            // which used to produce the INCORRECT result: (20.31970698357233, 46.76654261437082, NaN)

            // JTS:       if (p1.equals2D(q1)
            // JTS:               || p1.equals2D(q2)) {
            // JTS:           intPt[0] = p1;
            // JTS:       }
            // JTS:       else if (p2.equals2D(q1)
            // JTS:               || p2.equals2D(q2)) {
            // JTS:           intPt[0] = p2;
            // JTS:       }
            if p1 == q1 || p1 == q2 {
                self.int_pt[0] = p1;
            } else if p2 == q1 || p2 == q2 {
                self.int_pt[0] = p2;
            // JTS:
            // JTS:       /**
            // JTS:        * Now check to see if any endpoint lies on the interior of the other segment.
            // JTS:        */
            // JTS:       else if (Pq1 == 0) {
            // JTS:         intPt[0] = new Coordinate(q1);
            // JTS:       }
            // JTS:       else if (Pq2 == 0) {
            // JTS:         intPt[0] = new Coordinate(q2);
            // JTS:       }
            // JTS:       else if (Qp1 == 0) {
            // JTS:         intPt[0] = new Coordinate(p1);
            // JTS:       }
            // JTS:       else if (Qp2 == 0) {
            // JTS:         intPt[0] = new Coordinate(p2);
            // JTS:       }
            // JTS:     }
            // Now check to see if any endpoint lies on the interior of the other segment.
            } else if p_q1 == Orientation::Collinear {
                self.int_pt[0] = q1;
            } else if p_q2 == Orientation::Collinear {
                self.int_pt[0] = q2;
            } else if q_p1 == Orientation::Collinear {
                self.int_pt[0] = p1;
            } else if q_p2 == Orientation::Collinear {
                self.int_pt[0] = p2;
            }
        // JTS:     else {
        } else {
            // JTS:       isProper = true;
            // JTS:       intPt[0] = intersection(p1, p2, q1, q2);
            self.is_proper = true;
            self.int_pt[0] = self.intersection(p1, p2, q1, q2);
        }
        // JTS:     }
        // JTS:     return POINT_INTERSECTION;
        // JTS:   }
        Intersection::PointIntersection
    }
}

impl<F: GeoFloat> RobustLineIntersector<F> {
    // JTS:   private int computeCollinearIntersection(Coordinate p1, Coordinate p2,
    // JTS:       Coordinate q1, Coordinate q2) {
    // CLEANUP: take (p: Line<F>, q: Line<F>) instead?
    fn compute_collinear_intersection(
        &mut self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Intersection {
        // JTS:     boolean p1q1p2 = Envelope.intersects(p1, p2, q1);
        // JTS:     boolean p1q2p2 = Envelope.intersects(p1, p2, q2);
        // JTS:     boolean q1p1q2 = Envelope.intersects(q1, q2, p1);
        // JTS:     boolean q1p2q2 = Envelope.intersects(q1, q2, p2);
        let p_rect = Rect::new(p1, p2);
        let p1q1p2 = p_rect.intersects(&q1);
        let p1q2p2 = p_rect.intersects(&q2);

        let q_rect = Rect::new(q1, q2);
        let q1p1q2 = q_rect.intersects(&p1);
        let q1p2q2 = q_rect.intersects(&p2);

        // JTS:     if (p1q1p2 && p1q2p2) {
        // JTS:       intPt[0] = q1;
        // JTS:       intPt[1] = q2;
        // JTS:       return COLLINEAR_INTERSECTION;
        // JTS:     }
        // JTS:     if (q1p1q2 && q1p2q2) {
        // JTS:       intPt[0] = p1;
        // JTS:       intPt[1] = p2;
        // JTS:       return COLLINEAR_INTERSECTION;
        // JTS:     }
        // JTS:     if (p1q1p2 && q1p1q2) {
        // JTS:       intPt[0] = q1;
        // JTS:       intPt[1] = p1;
        // JTS:       return q1.equals(p1) && !p1q2p2 && !q1p2q2 ? POINT_INTERSECTION : COLLINEAR_INTERSECTION;
        // JTS:     }
        // JTS:     if (p1q1p2 && q1p2q2) {
        // JTS:       intPt[0] = q1;
        // JTS:       intPt[1] = p2;
        // JTS:       return q1.equals(p2) && !p1q2p2 && !q1p1q2 ? POINT_INTERSECTION : COLLINEAR_INTERSECTION;
        // JTS:     }
        // JTS:     if (p1q2p2 && q1p1q2) {
        // JTS:       intPt[0] = q2;
        // JTS:       intPt[1] = p1;
        // JTS:       return q2.equals(p1) && !p1q1p2 && !q1p2q2 ? POINT_INTERSECTION : COLLINEAR_INTERSECTION;
        // JTS:     }
        // JTS:     if (p1q2p2 && q1p2q2) {
        // JTS:       intPt[0] = q2;
        // JTS:       intPt[1] = p2;
        // JTS:       return q2.equals(p2) && !p1q1p2 && !q1p1q2 ? POINT_INTERSECTION : COLLINEAR_INTERSECTION;
        // JTS:     }
        // JTS:     return NO_INTERSECTION;
        // JTS:   }
        if p1q1p2 && p1q2p2 {
            self.int_pt[0] = q1;
            self.int_pt[1] = q2;
            return Intersection::CollinearIntersection;
        }
        if q1p1q2 && q1p2q2 {
            self.int_pt[0] = p1;
            self.int_pt[1] = p2;
            return Intersection::CollinearIntersection;
        }
        if p1q1p2 && q1p1q2 {
            self.int_pt[0] = q1;
            self.int_pt[1] = p1;
            if q1 == p1 && !p1q2p2 && !q1p2q2 {
                return Intersection::PointIntersection;
            } else {
                return Intersection::CollinearIntersection;
            }
        }
        if p1q1p2 && q1p2q2 {
            self.int_pt[0] = q1;
            self.int_pt[1] = p2;
            if q1 == p2 && !p1q2p2 && !q1p1q2 {
                return Intersection::PointIntersection;
            } else {
                return Intersection::CollinearIntersection;
            }
        }
        if p1q2p2 && q1p1q2 {
            self.int_pt[0] = q2;
            self.int_pt[1] = p1;
            if q2 == p1 && !p1q1p2 && !q1p2q2 {
                return Intersection::PointIntersection;
            } else {
                return Intersection::CollinearIntersection;
            }
        }
        if p1q2p2 && q1p2q2 {
            self.int_pt[0] = q2;
            self.int_pt[1] = p2;
            if q2 == p2 && !p1q1p2 && !q1p1q2 {
                return Intersection::PointIntersection;
            } else {
                return Intersection::CollinearIntersection;
            }
        }
        return Intersection::NoIntersection;
    }

    // JTS:   /**
    // JTS:    * This method computes the actual value of the intersection point.
    // JTS:    * To obtain the maximum precision from the intersection calculation,
    // JTS:    * the coordinates are normalized by subtracting the minimum
    // JTS:    * ordinate values (in absolute value).  This has the effect of
    // JTS:    * removing common significant digits from the calculation to
    // JTS:    * maintain more bits of precision.
    // JTS:    */
    // JTS:   private Coordinate intersection(
    // JTS:     Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2)
    // JTS:   {
    /// This method computes the actual value of the intersection point.
    /// To obtain the maximum precision from the intersection calculation,
    /// the coordinates are normalized by subtracting the minimum
    /// ordinate values (in absolute value).  This has the effect of
    /// removing common significant digits from the calculation to
    /// maintain more bits of precision.
    fn intersection(
        &self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Coordinate<F> {
        // JTS:     Coordinate intPt = intersectionSafe(p1, p2, q1, q2);
        let mut int_pt = self.intersection_safe(p1, p2, q1, q2);
        // JTS:
        // JTS:     /*
        // JTS:     // TESTING ONLY
        // JTS:     Coordinate intPtDD = CGAlgorithmsDD.intersection(p1, p2, q1, q2);
        // JTS:     double dist = intPt.distance(intPtDD);
        // JTS:     System.out.println(intPt + " - " + intPtDD + " dist = " + dist);
        // JTS:     //intPt = intPtDD;
        // JTS:     */
        // JTS:
        // JTS:     /**
        // JTS:      * Due to rounding it can happen that the computed intersection is
        // JTS:      * outside the envelopes of the input segments.  Clearly this
        // JTS:      * is inconsistent.
        // JTS:      * This code checks this condition and forces a more reasonable answer
        // JTS:      *
        // JTS:      * MD - May 4 2005 - This is still a problem.  Here is a failure case:
        // JTS:      *
        // JTS:      * LINESTRING (2089426.5233462777 1180182.3877339689, 2085646.6891757075 1195618.7333999649)
        // JTS:      * LINESTRING (1889281.8148903656 1997547.0560044837, 2259977.3672235999 483675.17050843034)
        // JTS:      * int point = (2097408.2633752143,1144595.8008114607)
        // JTS:      *
        // JTS:      * MD - Dec 14 2006 - This does not seem to be a failure case any longer
        // JTS:      */
        // JTS:     if (! isInSegmentEnvelopes(intPt)) {
        if !self.is_in_segment_envelopes(&int_pt) {
            // JTS: //      System.out.println("Intersection outside segment envelopes: " + intPt);
            // JTS:
            // JTS:       // compute a safer result
            // JTS:       // copy the coordinate, since it may be rounded later
            // JTS:       intPt = new Coordinate(nearestEndpoint(p1, p2, q1, q2));
            // JTS: //    intPt = CentralEndpointIntersector.getIntersection(p1, p2, q1, q2);
            // compute a safer result
            // copy the coordinate, since it may be rounded later
            int_pt = self.nearest_endpoint(p1, p2, q1, q2);

            // REVIEW: omitting this for now - `checkDD` does not change the output, but does a
            //         high precision accuracy check and logs if a discprency is found.
            // JTS: //      System.out.println("Segments: " + this);
            // JTS: //      System.out.println("Snapped to " + intPt);
            // JTS: //      checkDD(p1, p2, q1, q2, intPt);
            // JTS:     }
        }

        // TODO: do we want to introduce precisionModel?
        // JTS:     if (precisionModel != null) {
        // JTS:       precisionModel.makePrecise(intPt);
        // JTS:     }
        // JTS:     return intPt;
        // JTS:   }
        return int_pt;
    }
    // JTS:
    // JTS:   private void checkDD(Coordinate p1, Coordinate p2, Coordinate q1,
    // JTS:       Coordinate q2, Coordinate intPt)
    // JTS:   {
    // JTS:     Coordinate intPtDD = CGAlgorithmsDD.intersection(p1, p2, q1, q2);
    // JTS:     boolean isIn = isInSegmentEnvelopes(intPtDD);
    // JTS:     System.out.println(   "DD in env = " + isIn + "  --------------------- " + intPtDD);
    // JTS:     if (intPt.distance(intPtDD) > 0.0001) {
    // JTS:       System.out.println("Distance = " + intPt.distance(intPtDD));
    // JTS:     }
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Computes a segment intersection using homogeneous coordinates.
    // JTS:    * Round-off error can cause the raw computation to fail,
    // JTS:    * (usually due to the segments being approximately parallel).
    // JTS:    * If this happens, a reasonable approximation is computed instead.
    // JTS:    *
    // JTS:    * @param p1 a segment endpoint
    // JTS:    * @param p2 a segment endpoint
    // JTS:    * @param q1 a segment endpoint
    // JTS:    * @param q2 a segment endpoint
    // JTS:    * @return the computed intersection point
    // JTS:    */
    // JTS:   private Coordinate intersectionSafe(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2)
    // JTS:   {
    fn intersection_safe(
        &self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Coordinate<F> {
        // JTS:     Coordinate intPt = Intersection.intersection(p1, p2, q1, q2);
        // JTS:     if (intPt == null)
        // JTS:       intPt = nearestEndpoint(p1, p2, q1, q2);
        // JTS:  //     System.out.println("Snapped to " + intPt);
        // JTS:     return intPt;
        // JTS:   }
        match line_intersection(p1, p2, q1, q2) {
            Some(c) => c,
            None => self.nearest_endpoint(p1, p2, q1, q2),
        }
    }

    // JTS:   /**
    // JTS:    * Tests whether a point lies in the envelopes of both input segments.
    // JTS:    * A correctly computed intersection point should return <code>true</code>
    // JTS:    * for this test.
    // JTS:    * Since this test is for debugging purposes only, no attempt is
    // JTS:    * made to optimize the envelope test.
    // JTS:    *
    // JTS:    * @return <code>true</code> if the input point lies within both input segment envelopes
    // JTS:    */
    // JTS:   private boolean isInSegmentEnvelopes(Coordinate intPt)
    // JTS:   {
    ///  Tests whether a point lies in the envelopes of both input segments.
    ///  A correctly computed intersection point should return <code>true</code>
    ///  for this test.
    ///  Since this test is for debugging purposes only, no attempt is
    ///  made to optimize the envelope test.
    ///  
    ///  @return <code>true</code> if the input point lies within both input segment envelopes
    fn is_in_segment_envelopes(&self, int_pt: &Coordinate<F>) -> bool {
        // JTS:     Envelope env0 = new Envelope(inputLines[0][0], inputLines[0][1]);
        // JTS:     Envelope env1 = new Envelope(inputLines[1][0], inputLines[1][1]);
        // JTS:     return env0.contains(intPt) && env1.contains(intPt);
        // JTS:   }
        // CLEANUP: are we using input_lines anywhere else? It'd be nice to
        // get rid of that state and just pass around as params.
        let env0 = Rect::new(self.input_lines[0][0], self.input_lines[0][1]);
        let env1 = Rect::new(self.input_lines[1][0], self.input_lines[1][1]);
        return env0.contains(int_pt) && env1.contains(int_pt);
    }

    // JTS:
    // JTS:   /**
    // JTS:    * Finds the endpoint of the segments P and Q which
    // JTS:    * is closest to the other segment.
    // JTS:    * This is a reasonable surrogate for the true
    // JTS:    * intersection points in ill-conditioned cases
    // JTS:    * (e.g. where two segments are nearly coincident,
    // JTS:    * or where the endpoint of one segment lies almost on the other segment).
    // JTS:    * <p>
    // JTS:    * This replaces the older CentralEndpoint heuristic,
    // JTS:    * which chose the wrong endpoint in some cases
    // JTS:    * where the segments had very distinct slopes
    // JTS:    * and one endpoint lay almost on the other segment.
    // JTS:    *
    // JTS:    * @param p1 an endpoint of segment P
    // JTS:    * @param p2 an endpoint of segment P
    // JTS:    * @param q1 an endpoint of segment Q
    // JTS:    * @param q2 an endpoint of segment Q
    // JTS:    * @return the nearest endpoint to the other segment
    // JTS:    */
    // JTS:   private static Coordinate nearestEndpoint(Coordinate p1, Coordinate p2,
    // JTS:       Coordinate q1, Coordinate q2)
    // JTS:   {
    /// Finds the endpoint of the segments P and Q which is closest to the other segment.  This is
    /// a reasonable surrogate for the true intersection points in ill-conditioned cases (e.g.
    /// where two segments are nearly coincident, or where the endpoint of one segment lies almost
    /// on the other segment).
    ///
    /// This replaces the older CentralEndpoint heuristic, which chose the wrong endpoint in some
    /// cases where the segments had very distinct slopes and one endpoint lay almost on the other
    /// segment.
    ///
    /// @param p1 an endpoint of segment P
    /// @param p2 an endpoint of segment P
    /// @param q1 an endpoint of segment Q
    /// @param q2 an endpoint of segment Q
    /// @return the nearest endpoint to the other segment
    fn nearest_endpoint(
        &self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Coordinate<F> {
        use geo_types::private_utils::line_segment_distance;

        // JTS:     Coordinate nearestPt = p1;
        // JTS:     double minDist = Distance.pointToSegment(p1, q1, q2);
        let mut nearest_pt = p1;
        let mut min_dist = line_segment_distance(p1, q1, q2);

        // JTS:     double dist = Distance.pointToSegment(p2, q1, q2);
        // JTS:     if (dist < minDist) {
        // JTS:       minDist = dist;
        // JTS:       nearestPt = p2;
        // JTS:     }
        let dist = line_segment_distance(p2, q1, q1);
        if dist < min_dist {
            min_dist = dist;
            nearest_pt = p2;
        }
        // JTS:     dist = Distance.pointToSegment(q1, p1, p2);
        // JTS:     if (dist < minDist) {
        // JTS:       minDist = dist;
        // JTS:       nearestPt = q1;
        // JTS:     }
        let dist = line_segment_distance(q1, p1, p2);
        if dist < min_dist {
            min_dist = dist;
            nearest_pt = q1;
        }
        // JTS:     dist = Distance.pointToSegment(q2, p1, p2);
        // JTS:     if (dist < minDist) {
        // JTS:       minDist = dist;
        // JTS:       nearestPt = q2;
        // JTS:     }
        let dist = line_segment_distance(q2, p1, p2);
        if dist < min_dist {
            nearest_pt = q2;
        }
        // JTS:     return nearestPt;
        // JTS:   }
        // JTS:
        // JTS: }
        nearest_pt
    }
}

// From modules/core/src/main/java/org/locationtech/jts/algorithm/Intersection.java
// JTS: /*
// JTS:  * Copyright (c) 2019 martin Davis
// JTS:  *
// JTS:  * All rights reserved. This program and the accompanying materials
// JTS:  * are made available under the terms of the Eclipse Public License 2.0
// JTS:  * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
// JTS:  * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
// JTS:  * and the Eclipse Distribution License is available at
// JTS:  *
// JTS:  * http://www.eclipse.org/org/documents/edl-v10.php.
// JTS:  */
// JTS: package org.locationtech.jts.algorithm;
// JTS:
// JTS: import org.locationtech.jts.geom.Coordinate;
// JTS:
// JTS: /**
// JTS:  * Contains functions to compute intersections between lines.
// JTS:  *
// JTS:  * @author Martin Davis
// JTS:  *
// JTS:  */
// JTS: public class Intersection {
// JTS:
// JTS:   /**
// JTS:    * Computes the intersection point of two lines.
// JTS:    * If the lines are parallel or collinear this case is detected
// JTS:    * and <code>null</code> is returned.
// JTS:    * <p>
// JTS:    * In general it is not possible to accurately compute
// JTS:    * the intersection point of two lines, due to
// JTS:    * numerical roundoff.
// JTS:    * This is particularly true when the input lines are nearly parallel.
// JTS:    * This routine uses numerical conditioning on the input values
// JTS:    * to ensure that the computed value should be very close to the correct value.
// JTS:    *
// JTS:    * @param p1 an endpoint of line 1
// JTS:    * @param p2 an endpoint of line 1
// JTS:    * @param q1 an endpoint of line 2
// JTS:    * @param q2 an endpoint of line 2
// JTS:    * @return the intersection point between the lines, if there is one,
// JTS:    * or null if the lines are parallel or collinear
// JTS:    *
// JTS:    * @see CGAlgorithmsDD#intersection(Coordinate, Coordinate, Coordinate, Coordinate)
// JTS:    */
// JTS:   public static Coordinate intersection(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
fn line_intersection<F: GeoFloat>(
    p1: Coordinate<F>,
    p2: Coordinate<F>,
    q1: Coordinate<F>,
    q2: Coordinate<F>,
) -> Option<Coordinate<F>> {
    // JTS:     // compute midpoint of "kernel envelope"
    // JTS:     double minX0 = p1.x < p2.x ? p1.x : p2.x;
    // JTS:     double minY0 = p1.y < p2.y ? p1.y : p2.y;
    // JTS:     double maxX0 = p1.x > p2.x ? p1.x : p2.x;
    // JTS:     double maxY0 = p1.y > p2.y ? p1.y : p2.y;
    // compute midpoint of "kernel envelope"
    let min_x0 = if p1.x < p2.x { p1.x } else { p2.x };
    let min_y0 = if p1.y < p2.y { p1.y } else { p2.y };
    let max_x0 = if p1.x > p2.x { p1.x } else { p2.x };
    let max_y0 = if p1.y > p2.y { p1.y } else { p2.y };

    // JTS:     double minX1 = q1.x < q2.x ? q1.x : q2.x;
    // JTS:     double minY1 = q1.y < q2.y ? q1.y : q2.y;
    // JTS:     double maxX1 = q1.x > q2.x ? q1.x : q2.x;
    // JTS:     double maxY1 = q1.y > q2.y ? q1.y : q2.y;
    let min_x1 = if q1.x < q2.x { q1.x } else { q2.x };
    let min_y1 = if q1.y < q2.y { q1.y } else { q2.y };
    let max_x1 = if q1.x > q2.x { q1.x } else { q2.x };
    let max_y1 = if q1.y > q2.y { q1.y } else { q2.y };

    // JTS:     double intMinX = minX0 > minX1 ? minX0 : minX1;
    // JTS:     double intMaxX = maxX0 < maxX1 ? maxX0 : maxX1;
    // JTS:     double intMinY = minY0 > minY1 ? minY0 : minY1;
    // JTS:     double intMaxY = maxY0 < maxY1 ? maxY0 : maxY1;
    let int_min_x = if min_x0 > min_x1 { min_x0 } else { min_x1 };
    let int_max_x = if max_x0 < max_x1 { max_x0 } else { max_x1 };
    let int_min_y = if min_y0 > min_y1 { min_y0 } else { min_y1 };
    let int_max_y = if max_y0 < max_y1 { max_y0 } else { max_y1 };

    // JTS:     double midx = (intMinX + intMaxX) / 2.0;
    // JTS:     double midy = (intMinY + intMaxY) / 2.0;
    let mid_x = (int_min_x + int_max_x) / (F::one() + F::one());
    let mid_y = (int_min_y + int_max_y) / (F::one() + F::one());

    // JTS:     // condition ordinate values by subtracting midpoint
    // JTS:     double p1x = p1.x - midx;
    // JTS:     double p1y = p1.y - midy;
    // JTS:     double p2x = p2.x - midx;
    // JTS:     double p2y = p2.y - midy;
    // JTS:     double q1x = q1.x - midx;
    // JTS:     double q1y = q1.y - midy;
    // JTS:     double q2x = q2.x - midx;
    // JTS:     double q2y = q2.y - midy;
    // condition ordinate values by subtracting midpoint
    let p1x = p1.x - mid_x;
    let p1y = p1.y - mid_y;
    let p2x = p2.x - mid_x;
    let p2y = p2.y - mid_y;
    let q1x = q1.x - mid_x;
    let q1y = q1.y - mid_y;
    let q2x = q2.x - mid_x;
    let q2y = q2.y - mid_y;

    // JTS:     // unrolled computation using homogeneous coordinates eqn
    // JTS:     double px = p1y - p2y;
    // JTS:     double py = p2x - p1x;
    // JTS:     double pw = p1x * p2y - p2x * p1y;
    // JTS:
    // JTS:     double qx = q1y - q2y;
    // JTS:     double qy = q2x - q1x;
    // JTS:     double qw = q1x * q2y - q2x * q1y;
    // JTS:
    // JTS:     double x = py * qw - qy * pw;
    // JTS:     double y = qx * pw - px * qw;
    // JTS:     double w = px * qy - qx * py;
    // JTS:
    // JTS:     double xInt = x/w;
    // JTS:     double yInt = y/w;
    // unrolled computation using homogeneous coordinates eqn
    let px = p1y - p2y;
    let py = p2x - p1x;
    let pw = p1x * p2y - p2x * p1y;

    let qx = q1y - q2y;
    let qy = q2x - q1x;
    let qw = q1x * q2y - q2x * q1y;

    let x = py * qw - qy * pw;
    let y = qx * pw - px * qw;
    let w = px * qy - qx * py;

    let x_int = x / w;
    let y_int = y / w;

    // JTS:     // check for parallel lines
    // JTS:     if ((Double.isNaN(xInt)) || (Double.isInfinite(xInt)
    // JTS:         || Double.isNaN(yInt)) || (Double.isInfinite(yInt))) {
    // JTS:       return null;
    // JTS:     }
    // JTS:     // de-condition intersection point
    // JTS:     return new Coordinate(xInt + midx, yInt + midy);
    // JTS:   }
    // check for parallel lines
    if (x_int.is_nan() || x_int.is_infinite()) || (y_int.is_nan() || y_int.is_infinite()) {
        None
    } else {
        // de-condition intersection point
        Some(Coordinate {
            x: x_int + mid_x,
            y: y_int + mid_y,
        })
    }
}
// JTS: }

#[cfg(test)]
mod test {
    use super::*;
    use geo_types::Line;

    // JTS: Tests
    // JTS:  /**
    // JTS:    * Following cases were failures when using the CentralEndpointIntersector heuristic.
    // JTS:    * This is because one segment lies at a significant angle to the other,
    // JTS:    * with only one endpoint is close to the other segment.
    // JTS:    * The CE heuristic chose the wrong endpoint to return.
    // JTS:    * The fix is to use a new heuristic which out of the 4 endpoints
    // JTS:    * chooses the one which is closest to the other segment.
    // JTS:    * This works in all known failure cases.
    // JTS:    *
    // JTS:    * @throws ParseException
    // JTS:    */
    // JTS:     public void testCentralEndpointHeuristicFailure()
    // JTS:     throws ParseException
    // JTS:     {
    // JTS:       checkIntersection(
    // JTS:           "LINESTRING (163.81867067 -211.31840378, 165.9174252 -214.1665075)",
    // JTS:           "LINESTRING (2.84139601 -57.95412726, 469.59990601 -502.63851732)",
    // JTS:           1,
    // JTS:           "POINT (163.81867067 -211.31840378)",
    // JTS:           0);
    // JTS:     }
    // JTS:
    // JTS:   public void testCentralEndpointHeuristicFailure2()
    // JTS:   throws ParseException
    // JTS:   {
    // JTS:     checkIntersection(
    // JTS:         "LINESTRING (-58.00593335955 -1.43739086465, -513.86101637525 -457.29247388035)",
    // JTS:         "LINESTRING (-215.22279674875 -158.65425425385, -218.1208801283 -160.68343590235)",
    // JTS:         1,
    // JTS:         "POINT ( -215.22279674875 -158.65425425385 )",
    // JTS:         0);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Tests a case where intersection point is rounded,
    // JTS:    * and it is computed as a nearest endpoint.
    // JTS:    * Exposed a bug due to aliasing of endpoint.
    // JTS:    *
    // JTS:    * MD 8 Mar 2013
    // JTS:    *
    // JTS:    * @throws ParseException
    // JTS:    */
    // JTS:   public void testRoundedPointsNotAltered()
    // JTS:   throws ParseException
    // JTS:   {
    // JTS:     checkInputNotAltered(
    // JTS:         "LINESTRING (-58.00593335955 -1.43739086465, -513.86101637525 -457.29247388035)",
    // JTS:         "LINESTRING (-215.22279674875 -158.65425425385, -218.1208801283 -160.68343590235)",
    // JTS:         100000 );
    // JTS:   }
    // JTS:
    // JTS:
    // JTS: /**
    // JTS:  * Test from Tomas Fa - JTS list 6/13/2012
    // JTS:  *
    // JTS:  * Fails using original JTS DeVillers determine orientation test.
    // JTS:  * Succeeds using DD and Shewchuk orientation
    // JTS:  *
    // JTS:  * @throws ParseException
    // JTS:  */
    // JTS: public void testTomasFa_1()
    // JTS: throws ParseException
    // JTS: {
    // JTS:   checkIntersectionNone(
    // JTS:       "LINESTRING (-42.0 163.2, 21.2 265.2)",
    // JTS:       "LINESTRING (-26.2 188.7, 37.0 290.7)");
    // JTS: }
    // JTS:
    // JTS: /**
    // JTS:  * Test from Tomas Fa - JTS list 6/13/2012
    // JTS:  *
    // JTS:  * Fails using original JTS DeVillers determine orientation test.
    // JTS:  * Succeeds using DD and Shewchuk orientation
    // JTS:  *
    // JTS:  * @throws ParseException
    // JTS:  */
    // JTS: public void testTomasFa_2()
    // JTS: throws ParseException
    // JTS: {
    // JTS:   checkIntersectionNone(
    // JTS:       "LINESTRING (-5.9 163.1, 76.1 250.7)",
    // JTS:       "LINESTRING (14.6 185.0, 96.6 272.6)");
    // JTS: }
    // JTS:
    // JTS: /**
    // JTS:  * Test involving two non-almost-parallel lines.
    // JTS:  * Does not seem to cause problems with basic line intersection algorithm.
    // JTS:  *
    // JTS:  * @throws ParseException
    // JTS:  */
    // JTS: public void testLeduc_1()
    // JTS: throws ParseException
    // JTS: {
    // JTS:   checkIntersection(
    // JTS:       "LINESTRING (305690.0434123494 254176.46578338774, 305601.9999843455 254243.19999846347)",
    // JTS:       "LINESTRING (305689.6153764265 254177.33102743194, 305692.4999844298 254171.4999983967)",
    // JTS:       1,
    // JTS:       "POINT (305690.0434123494 254176.46578338774)",
    // JTS:       0);
    // JTS: }
    #[test]
    fn test_leduc_1() {
        let line_string_1 = Line::new(
            Coordinate {
                x: 305690.0434123494,
                y: 254176.46578338774,
            },
            Coordinate {
                x: 305601.9999843455,
                y: 254243.19999846347,
            },
        );
        let line_string_2 = Line::new(
            Coordinate {
                x: 305689.6153764265,
                y: 254177.33102743194,
            },
            Coordinate {
                x: 305692.4999844298,
                y: 254171.4999983967,
            },
        );
        let intersection_pt = Coordinate {
            x: 305690.0434123494,
            y: 254176.46578338774,
        };
        check_intersection(
            line_string_1,
            line_string_2,
            Intersection::PointIntersection,
            Some(intersection_pt),
            0f64,
        );
    }

    // JTS:
    // JTS:   /**
    // JTS:    * Test from strk which is bad in GEOS (2009-04-14).
    // JTS:    *
    // JTS:    * @throws ParseException
    // JTS:    */
    // JTS:   public void testGEOS_1()
    // JTS:   throws ParseException
    // JTS:   {
    // JTS:       checkIntersection(
    // JTS:               "LINESTRING (588750.7429703881 4518950.493668233, 588748.2060409798 4518933.9452804085)",
    // JTS:               "LINESTRING (588745.824857241 4518940.742239175, 588748.2060437313 4518933.9452791475)",
    // JTS:               1,
    // JTS:               "POINT (588748.2060416829 4518933.945284994)",
    // JTS:               0);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Test from strk which is bad in GEOS (2009-04-14).
    // JTS:    *
    // JTS:    * @throws ParseException
    // JTS:    */
    // JTS:   public void testGEOS_2()
    // JTS:   throws ParseException
    // JTS:   {
    // JTS:       checkIntersection(
    // JTS:               "LINESTRING (588743.626135934 4518924.610969561, 588732.2822865889 4518925.4314047815)",
    // JTS:               "LINESTRING (588739.1191384895 4518927.235700594, 588731.7854614238 4518924.578370095)",
    // JTS:               1,
    // JTS:               "POINT (588733.8306132929 4518925.319423238)",
    // JTS:               0);
    // JTS:   }
    // JTS:
    // JTS:       /**
    // JTS:        * This used to be a failure case (exception), but apparently works now.
    // JTS:        * Possibly normalization has fixed this?
    // JTS:        *
    // JTS:        * @throws ParseException
    // JTS:        */
    // JTS:   public void testDaveSkeaCase()
    // JTS:       throws ParseException
    // JTS:   {
    // JTS:       checkIntersection(
    // JTS:               "LINESTRING ( 2089426.5233462777 1180182.3877339689, 2085646.6891757075 1195618.7333999649 )",
    // JTS:               "LINESTRING ( 1889281.8148903656 1997547.0560044837, 2259977.3672235999 483675.17050843034 )",
    // JTS:               1,
    // JTS:               new Coordinate[] {
    // JTS:                       new Coordinate(2087536.6062609926, 1187900.560566967),
    // JTS:               }, 0);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Outside envelope using HCoordinate method.
    // JTS:    *
    // JTS:    * @throws ParseException
    // JTS:    */
    // JTS:   public void testCmp5CaseWKT()
    // JTS:   throws ParseException
    // JTS:   {
    // JTS:       checkIntersection(
    // JTS:               "LINESTRING (4348433.262114629 5552595.478385733, 4348440.849387404 5552599.272022122 )",
    // JTS:               "LINESTRING (4348433.26211463  5552595.47838573,  4348440.8493874   5552599.27202212  )",
    // JTS:               1,
    // JTS:               new Coordinate[] {
    // JTS:                       new Coordinate(4348440.8493874, 5552599.27202212),
    // JTS:               },
    // JTS:               0);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Result of this test should be the same as the WKT one!
    // JTS:    * @throws ParseException
    // JTS:    */
    // JTS:   public void testCmp5CaseRaw()
    // JTS:   throws ParseException
    // JTS:   {
    // JTS:       checkIntersection(
    // JTS:               new Coordinate[] {
    // JTS:                       new Coordinate(4348433.262114629, 5552595.478385733),
    // JTS:                       new Coordinate(4348440.849387404, 5552599.272022122),
    // JTS:
    // JTS:                       new Coordinate(4348433.26211463,  5552595.47838573),
    // JTS:                       new Coordinate(4348440.8493874,   5552599.27202212)
    // JTS:               },                1,
    // JTS:               new Coordinate[] {
    // JTS:                       new Coordinate(4348440.8493874, 5552599.27202212),
    // JTS:               },
    // JTS:               0);
    // JTS:   }
    // JTS:
    // JTS: void checkIntersectionNone(String wkt1, String wkt2)
    // JTS:   throws ParseException
    // JTS: {
    // JTS:   LineString l1 = (LineString) reader.read(wkt1);
    // JTS:   LineString l2 = (LineString) reader.read(wkt2);
    // JTS:   Coordinate[] pt = new Coordinate[] {
    // JTS:       l1.getCoordinateN(0), l1.getCoordinateN(1),
    // JTS:       l2.getCoordinateN(0), l2.getCoordinateN(1)
    // JTS:   };
    // JTS:   checkIntersection(pt, 0, null, 0);
    // JTS: }
    // JTS:
    // JTS: void checkIntersection(String wkt1, String wkt2,
    // JTS:     int expectedIntersectionNum,
    // JTS:     Coordinate[] intPt,
    // JTS:     double distanceTolerance)
    // JTS:   throws ParseException
    // JTS: {
    // JTS:   LineString l1 = (LineString) reader.read(wkt1);
    // JTS:   LineString l2 = (LineString) reader.read(wkt2);
    // JTS:   Coordinate[] pt = new Coordinate[] {
    // JTS:       l1.getCoordinateN(0), l1.getCoordinateN(1),
    // JTS:       l2.getCoordinateN(0), l2.getCoordinateN(1)
    // JTS:   };
    // JTS:   checkIntersection(pt, expectedIntersectionNum, intPt, distanceTolerance);
    // JTS: }
    // JTS:
    // JTS:   void checkIntersection(String wkt1, String wkt2,
    // JTS:           int expectedIntersectionNum,
    // JTS:           String expectedWKT,
    // JTS:           double distanceTolerance)
    // JTS:       throws ParseException
    // JTS:   {
    // JTS:       LineString l1 = (LineString) reader.read(wkt1);
    // JTS:       LineString l2 = (LineString) reader.read(wkt2);
    // JTS:       Coordinate[] pt = new Coordinate[] {
    // JTS:               l1.getCoordinateN(0), l1.getCoordinateN(1),
    // JTS:               l2.getCoordinateN(0), l2.getCoordinateN(1)
    // JTS:       };
    // JTS:       Geometry g = reader.read(expectedWKT);
    // JTS:       Coordinate[] intPt = g.getCoordinates();
    // JTS:       checkIntersection(pt, expectedIntersectionNum, intPt, distanceTolerance);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Check that intersection of segment defined by points in pt array
    // JTS:    * is equal to the expectedIntPt value (up to the given distanceTolerance).
    // JTS:    *
    // JTS:    * @param pt
    // JTS:    * @param expectedIntersectionNum
    // JTS:    * @param expectedIntPt the expected intersection points (maybe null if not tested)
    // JTS:    * @param distanceTolerance tolerance to use for equality test
    // JTS:    */
    // JTS:   void checkIntersection(Coordinate[] pt,
    // JTS:           int expectedIntersectionNum,
    // JTS:           Coordinate[] expectedIntPt,
    // JTS:           double distanceTolerance)
    // JTS:   {
    fn check_intersection<F: GeoFloat>(
        line_1: Line<F>,
        line_2: Line<F>,
        intersection: Intersection,
        interesection_pt: Option<Coordinate<F>>,
        distance_tolerance: F,
    ) {
        // JTS:       LineIntersector li = new RobustLineIntersector();
        // JTS:       li.computeIntersection(pt[0], pt[1], pt[2], pt[3]);
        let mut li = RobustLineIntersector::new();
        li.compute_intersection(line_1.start, line_1.end, line_2.start, line_2.end);
        // JTS:
        // JTS:       int intNum = li.getIntersectionNum();
        // JTS:       assertEquals("Number of intersections not as expected", expectedIntersectionNum, intNum);
        assert_eq!(li.result(), intersection);
        // JTS:
        // JTS:       if (expectedIntPt != null) {
        // JTS:           assertEquals("Wrong number of expected int pts provided", intNum, expectedIntPt.length);
        // JTS:           // test that both points are represented here
        // JTS:           boolean isIntPointsCorrect = true;
        // JTS:           if (intNum == 1) {
        // JTS:               checkIntPoints(expectedIntPt[0], li.getIntersection(0), distanceTolerance);
        // JTS:           }
        // JTS:           else if (intNum == 2) {
        // JTS:               checkIntPoints(expectedIntPt[1], li.getIntersection(0), distanceTolerance);
        // JTS:               checkIntPoints(expectedIntPt[1], li.getIntersection(0), distanceTolerance);
        // JTS:
        // JTS:               if (! (equals(expectedIntPt[0],li.getIntersection(0), distanceTolerance)
        // JTS:                       || equals(expectedIntPt[0],li.getIntersection(1), distanceTolerance) )) {
        // JTS:                   checkIntPoints(expectedIntPt[0], li.getIntersection(0), distanceTolerance);
        // JTS:                   checkIntPoints(expectedIntPt[0], li.getIntersection(1), distanceTolerance);
        // JTS:               }
        // JTS:               else if (! (equals(expectedIntPt[1],li.getIntersection(0), distanceTolerance)
        // JTS:                       || equals(expectedIntPt[1],li.getIntersection(1), distanceTolerance) )) {
        // JTS:                   checkIntPoints(expectedIntPt[1], li.getIntersection(0), distanceTolerance);
        // JTS:                   checkIntPoints(expectedIntPt[1], li.getIntersection(1), distanceTolerance);
        // JTS:               }
        // JTS:           }
        // JTS:       }
        // JTS:   }
        if let Some(interesection_pt) = interesection_pt {
            match intersection {
                Intersection::NoIntersection => {
                    panic!("shouldn't specify point if no intersection is expected")
                }
                Intersection::PointIntersection => {
                    check_intersection_points(
                        interesection_pt,
                        LineIntersector::intersection(&li, 0),
                        distance_tolerance,
                    );
                }
                Intersection::CollinearIntersection => todo!(),
            }
        }
    }
    // JTS:
    // JTS:   void checkIntPoints(Coordinate expectedPt, Coordinate actualPt, double distanceTolerance)
    // JTS:   {
    // JTS:       boolean isEqual = equals(expectedPt, actualPt, distanceTolerance);
    // JTS:       assertTrue("Int Pts not equal - "
    // JTS:               + "expected " + WKTWriter.toPoint(expectedPt) + " VS "
    // JTS:               + "actual " + WKTWriter.toPoint(actualPt), isEqual);
    // JTS:   }
    fn check_intersection_points<F: GeoFloat>(
        expected: Coordinate<F>,
        actual: Coordinate<F>,
        distance_tolerance: F,
    ) {
        use crate::algorithm::euclidean_distance::EuclideanDistance;
        let distance = expected.euclidean_distance(&actual);
        assert!(
            distance <= distance_tolerance,
            "expected distance: {:?} < tolerance: {:?}",
            distance,
            distance_tolerance
        );
    }

    // JTS:   public static boolean equals(Coordinate p0, Coordinate p1, double distanceTolerance)
    // JTS:   {
    // JTS:       return p0.distance(p1) <= distanceTolerance;
    // JTS:   }
    // JTS:
    // JTS: void checkInputNotAltered(String wkt1, String wkt2, int scaleFactor) throws ParseException
    // JTS: {
    // JTS:   LineString l1 = (LineString) reader.read(wkt1);
    // JTS:   LineString l2 = (LineString) reader.read(wkt2);
    // JTS:   Coordinate[] pt = new Coordinate[] { l1.getCoordinateN(0),
    // JTS:       l1.getCoordinateN(1), l2.getCoordinateN(0), l2.getCoordinateN(1) };
    // JTS:   checkInputNotAltered(pt, scaleFactor);
    // JTS: }
    // JTS:
    // JTS:   public void checkInputNotAltered(Coordinate[] pt, int scaleFactor)
    // JTS:     {
    // JTS:       // save input points
    // JTS:       Coordinate[] savePt = new Coordinate[4];
    // JTS:       for (int i = 0; i < 4; i++) {
    // JTS:         savePt[i] = new Coordinate(pt[i]);
    // JTS:       }
    // JTS:
    // JTS:       LineIntersector li = new RobustLineIntersector();
    // JTS:       li.setPrecisionModel(new PrecisionModel(scaleFactor));
    // JTS:       li.computeIntersection(pt[0], pt[1], pt[2], pt[3]);
    // JTS:
    // JTS:       // check that input points are unchanged
    // JTS:       for (int i = 0; i < 4; i++) {
    // JTS:         assertEquals("Input point " + i + " was altered - ", savePt[i], pt[i]);
    // JTS:       }
    // JTS:     }
}
