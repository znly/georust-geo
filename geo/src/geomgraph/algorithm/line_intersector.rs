use crate::{Coordinate, GeoFloat};

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum Intersection {
    /// Indicates that line segments do not intersect
    NoIntersection,

    /// Indicates that line segments intersect in a single point
    PointIntersection,

    /// Indicates that line segments intersect in a line segment
    CollinearIntersection,
}

pub(crate) trait LineIntersector<F: GeoFloat> {
    fn compute_edge_distance(&self, p: Coordinate<F>, p0: Coordinate<F>, p1: Coordinate<F>) -> F {
        let dx = (p1.x - p0.x).abs();
        let dy = (p1.y - p0.y).abs();

        let mut dist: F;
        if p == p0 {
            dist = F::zero();
        } else if p == p1 {
            if dx > dy {
                dist = dx;
            } else {
                dist = dy;
            }
        } else {
            let pdx = (p.x - p0.x).abs();
            let pdy = (p.y - p0.y).abs();
            if dx > dy {
                dist = pdx;
            } else {
                dist = pdy;
            }
            // hack to ensure that non-endpoints always have a non-zero distance
            if dist == F::zero() && p != p0 {
                dist = pdx.max(pdy);
            }
        }
        assert!(!(dist == F::zero() && p != p0), "Bad distance calculation");
        dist
    }

    fn result(&self) -> Intersection;
    fn set_result(&mut self, intersection: Intersection);

    fn input_lines(&self) -> [[Coordinate<F>; 2]; 2];
    fn set_input_lines(&mut self, line: usize, start_or_end: usize, coord: Coordinate<F>);

    fn compute_intersection(
        &mut self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        p3: Coordinate<F>,
        p4: Coordinate<F>,
    ) {
        // TODO: What do we need input_lines for? This API feels overly stateful, but I'm going to
        //       match JTS for now.
        self.set_input_lines(0, 0, p1);
        self.set_input_lines(0, 1, p2);
        self.set_input_lines(1, 0, p3);
        self.set_input_lines(1, 1, p4);

        // TODO: This API feels overly stateful, but I'm going to match JTS for now.
        let result = self.compute_intersect(p1, p2, p3, p4);
        self.set_result(result);
    }

    fn compute_intersect(
        &mut self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Intersection;

    /// Tests whether the input geometries intersect.
    fn has_intersection(&self) -> bool {
        self.result() != Intersection::NoIntersection
    }

    /// Returns the number of intersection points found.  This will be either 0, 1 or 2.
    ///
    /// @return the number of intersection points found (0, 1, or 2)
    fn intersection_num(&self) -> usize {
        match self.result() {
            Intersection::NoIntersection => 0,
            Intersection::PointIntersection => 1,
            Intersection::CollinearIntersection => 2,
        }
    }

    fn intersection(&self, intersection_index: usize) -> Coordinate<F>;

    /// Test whether a point is a intersection point of two line segments.
    ///
    /// Note that if the intersection is a line segment, this method only tests for equality with
    /// the endpoints of the intersection segment.  
    ///
    /// It does <b>not</b> return true if the input point is internal to the intersection segment.
    ///
    /// @return true if the input point is one of the intersection points.
    fn is_intersection(&self, coord: &Coordinate<F>) -> bool {
        for i in 0..self.intersection_num() {
            if self.intersection(i) == *coord {
                return true;
            }
        }
        return false;
    }

    fn is_proper(&self) -> bool;
    fn edge_distance(&self, segment_index: usize, intersection_index: usize) -> F {
        self.compute_edge_distance(
            self.intersection(intersection_index),
            self.input_lines()[segment_index][0],
            self.input_lines()[segment_index][1],
        )
    }
}
