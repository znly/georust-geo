use super::{Intersection, LineIntersector};
use crate::algorithm::kernels::{Kernel, Orientation, RobustKernel};
use crate::contains::Contains;
use crate::intersects::Intersects;
use crate::num_traits::Zero;
use geo_types::{Coordinate, Rect};

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
pub(crate) struct RobustLineIntersector<F: num_traits::Float> {
    // TODO: JTS captures some state in the LineIntersector. I'm not sure if it's helpful. Roughly it
    // seems to be mid-computation state and result state. Perhaps that could be better modeled and
    // we could leave the LineIntersector less mutable.
    result: Option<Intersection>,
    // TODO: make this [geo_types::Line; 2]? matching JTS for now
    is_proper: bool,
    input_lines: [[Coordinate<F>; 2]; 2],
    int_pt: [Coordinate<F>; 2],
}

impl<F: num_traits::Float> RobustLineIntersector<F> {
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

impl<F: num_traits::Float> LineIntersector<F> for RobustLineIntersector<F> {
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
        // first try a fast test to see if the envelopes of the lines intersect
        // TODO: This assumes https://github.com/georust/geo/pull/519 is merged
        // otherwise we need to order points by min/max
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
        // TODO: might be nicer if robust::orient2d could take a Coordinate directly (e.g. similar
        // to https://github.com/georust/proj/pull/41, robust::Coord could be a trait, which
        // geo-types::Coordinate could conform to
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
            // JTS:       		|| p1.equals2D(q2)) {
            // JTS:       	intPt[0] = p1;
            // JTS:       }
            // JTS:       else if (p2.equals2D(q1)
            // JTS:       		|| p2.equals2D(q2)) {
            // JTS:       	intPt[0] = p2;
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
        }
        Intersection::PointIntersection
    }
}

impl<F: num_traits::Float> RobustLineIntersector<F> {
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
fn line_intersection<F: num_traits::Float>(
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
