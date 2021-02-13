use super::{Intersection, LineIntersector};
use crate::algorithm::kernels::{Kernel, Orientation, RobustKernel};
use crate::contains::Contains;
use crate::intersects::Intersects;
use crate::num_traits::Zero;
use crate::{Coordinate, GeoFloat, Rect};

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

    fn compute_intersect(
        &mut self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Intersection {
        self.is_proper = false;
        // first try a fast test to see if the envelopes of the lines intersect
        if !Rect::new(p1, p2).intersects(&Rect::new(q1, q2)) {
            return Intersection::NoIntersection;
        }

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

        let q_p1 = RobustKernel::orient2d(q1, q2, p1);
        let q_p2 = RobustKernel::orient2d(q1, q2, p2);
        match (q_p1, q_p2) {
            (Orientation::Clockwise, Orientation::Clockwise)
            | (Orientation::CounterClockwise, Orientation::CounterClockwise) => {
                return Intersection::NoIntersection
            }
            _ => (),
        }

        if let (
            Orientation::Collinear,
            Orientation::Collinear,
            Orientation::Collinear,
            Orientation::Collinear,
        ) = (p_q1, p_q2, q_p1, q_p2)
        {
            return self.compute_collinear_intersection(p1, p2, q1, q2);
        }
        // At this point we know that there is a single intersection point (since the lines are not
        // collinear).
        //
        // Check if the intersection is an endpoint. If it is, copy the endpoint as the
        // intersection point. Copying the point rather than computing it ensures the point has the
        // exact value, which is important for robustness. It is sufficient to simply check for an
        // endpoint which is on the other line, since at this point we know that the inputLines
        // must intersect.
        if p_q1 == Orientation::Collinear
            || p_q2 == Orientation::Collinear
            || q_p1 == Orientation::Collinear
            || q_p2 == Orientation::Collinear
        {
            self.is_proper = false;
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

            if p1 == q1 || p1 == q2 {
                self.int_pt[0] = p1;
            } else if p2 == q1 || p2 == q2 {
                self.int_pt[0] = p2;
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
        } else {
            self.is_proper = true;
            self.int_pt[0] = self.intersection(p1, p2, q1, q2);
        }
        Intersection::PointIntersection
    }
}

impl<F: GeoFloat> RobustLineIntersector<F> {
    // CLEANUP: take (p: Line<F>, q: Line<F>) instead?
    fn compute_collinear_intersection(
        &mut self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Intersection {
        let p_rect = Rect::new(p1, p2);
        let p1q1p2 = p_rect.intersects(&q1);
        let p1q2p2 = p_rect.intersects(&q2);

        let q_rect = Rect::new(q1, q2);
        let q1p1q2 = q_rect.intersects(&p1);
        let q1p2q2 = q_rect.intersects(&p2);

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
        let mut int_pt = self.intersection_safe(p1, p2, q1, q2);
        if !self.is_in_segment_envelopes(&int_pt) {
            // compute a safer result
            // copy the coordinate, since it may be rounded later
            int_pt = self.nearest_endpoint(p1, p2, q1, q2);

            // REVIEW: omitting this for now - `checkDD` does not change the output, but does a
            //         high precision accuracy check and logs if a discprency is found.
        }

        // TODO: do we want to introduce precisionModel?
        return int_pt;
    }
    fn intersection_safe(
        &self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Coordinate<F> {
        match line_intersection(p1, p2, q1, q2) {
            Some(c) => c,
            None => self.nearest_endpoint(p1, p2, q1, q2),
        }
    }

    ///  Tests whether a point lies in the envelopes of both input segments.
    ///  A correctly computed intersection point should return <code>true</code>
    ///  for this test.
    ///  Since this test is for debugging purposes only, no attempt is
    ///  made to optimize the envelope test.
    ///  
    ///  @return <code>true</code> if the input point lies within both input segment envelopes
    fn is_in_segment_envelopes(&self, int_pt: &Coordinate<F>) -> bool {
        // CLEANUP: are we using input_lines anywhere else? It'd be nice to
        // get rid of that state and just pass around as params.
        let env0 = Rect::new(self.input_lines[0][0], self.input_lines[0][1]);
        let env1 = Rect::new(self.input_lines[1][0], self.input_lines[1][1]);
        return env0.contains(int_pt) && env1.contains(int_pt);
    }

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

        let mut nearest_pt = p1;
        let mut min_dist = line_segment_distance(p1, q1, q2);

        let dist = line_segment_distance(p2, q1, q1);
        if dist < min_dist {
            min_dist = dist;
            nearest_pt = p2;
        }
        let dist = line_segment_distance(q1, p1, p2);
        if dist < min_dist {
            min_dist = dist;
            nearest_pt = q1;
        }
        let dist = line_segment_distance(q2, p1, p2);
        if dist < min_dist {
            nearest_pt = q2;
        }
        nearest_pt
    }
}

// From modules/core/src/main/java/org/locationtech/jts/algorithm/Intersection.java
fn line_intersection<F: GeoFloat>(
    p1: Coordinate<F>,
    p2: Coordinate<F>,
    q1: Coordinate<F>,
    q2: Coordinate<F>,
) -> Option<Coordinate<F>> {
    // compute midpoint of "kernel envelope"
    let min_x0 = if p1.x < p2.x { p1.x } else { p2.x };
    let min_y0 = if p1.y < p2.y { p1.y } else { p2.y };
    let max_x0 = if p1.x > p2.x { p1.x } else { p2.x };
    let max_y0 = if p1.y > p2.y { p1.y } else { p2.y };

    let min_x1 = if q1.x < q2.x { q1.x } else { q2.x };
    let min_y1 = if q1.y < q2.y { q1.y } else { q2.y };
    let max_x1 = if q1.x > q2.x { q1.x } else { q2.x };
    let max_y1 = if q1.y > q2.y { q1.y } else { q2.y };

    let int_min_x = if min_x0 > min_x1 { min_x0 } else { min_x1 };
    let int_max_x = if max_x0 < max_x1 { max_x0 } else { max_x1 };
    let int_min_y = if min_y0 > min_y1 { min_y0 } else { min_y1 };
    let int_max_y = if max_y0 < max_y1 { max_y0 } else { max_y1 };

    let mid_x = (int_min_x + int_max_x) / (F::one() + F::one());
    let mid_y = (int_min_y + int_max_y) / (F::one() + F::one());

    // condition ordinate values by subtracting midpoint
    let p1x = p1.x - mid_x;
    let p1y = p1.y - mid_y;
    let p2x = p2.x - mid_x;
    let p2y = p2.y - mid_y;
    let q1x = q1.x - mid_x;
    let q1y = q1.y - mid_y;
    let q2x = q2.x - mid_x;
    let q2y = q2.y - mid_y;

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

#[cfg(test)]
mod test {
    use super::*;
    use geo_types::Line;

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

    fn check_intersection<F: GeoFloat>(
        line_1: Line<F>,
        line_2: Line<F>,
        intersection: Intersection,
        interesection_pt: Option<Coordinate<F>>,
        distance_tolerance: F,
    ) {
        let mut li = RobustLineIntersector::new();
        li.compute_intersection(line_1.start, line_1.end, line_2.start, line_2.end);
        assert_eq!(li.result(), intersection);
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
}
