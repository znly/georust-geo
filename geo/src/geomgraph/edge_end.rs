use super::{Edge, Label, Node, Quadrant};
use crate::{Coordinate, GeoFloat};

use std::cell::RefCell;
use std::fmt;

#[derive(Clone)]
pub(crate) struct EdgeEnd<F>
where
    F: GeoFloat,
{
    // edge: RefCell<Edge<F>>,
    label: Label,
    coord_0: Coordinate<F>,
    coord_1: Coordinate<F>,
    delta: Coordinate<F>,
    quadrant: Quadrant,
}

impl<F: GeoFloat> fmt::Debug for EdgeEnd<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EdgeEnd")
            .field("label", &self.label)
            .field(
                "coords",
                &format!("{:?} -> {:?}", &self.coord_0, &self.coord_1),
            )
            .field("quadrant", &self.quadrant)
            .finish()
    }
}

impl<F> EdgeEnd<F>
where
    F: GeoFloat,
{
    pub fn new(
        //edge: RefCell<Edge<F>>,
        coord_0: Coordinate<F>,
        coord_1: Coordinate<F>,
        label: Label,
    ) -> EdgeEnd<F> {
        let delta = coord_1 - coord_0;
        let quadrant = Quadrant::new(delta.x, delta.y);
        EdgeEnd {
            // edge,
            label,
            coord_0,
            coord_1,
            delta,
            quadrant,
        }
    }

    pub fn label(&self) -> &Label {
        &self.label
    }

    pub fn label_mut(&mut self) -> &mut Label {
        &mut self.label
    }

    pub fn coordinate(&self) -> &Coordinate<F> {
        &self.coord_0
    }

    pub fn directed_coordinate(&self) -> &Coordinate<F> {
        &self.coord_1
    }
}

impl<F> std::cmp::Eq for EdgeEnd<F> where F: GeoFloat {}

impl<F> std::cmp::PartialEq for EdgeEnd<F>
where
    F: GeoFloat,
{
    fn eq(&self, other: &EdgeEnd<F>) -> bool {
        self.delta == other.delta
    }
}

impl<F> std::cmp::PartialOrd for EdgeEnd<F>
where
    F: GeoFloat,
{
    fn partial_cmp(&self, other: &EdgeEnd<F>) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<F> std::cmp::Ord for EdgeEnd<F>
where
    F: GeoFloat,
{
    fn cmp(&self, other: &EdgeEnd<F>) -> std::cmp::Ordering {
        self.compare_direction(other)
    }
}

impl<F> EdgeEnd<F>
where
    F: GeoFloat,
{
    pub(crate) fn compare_direction(&self, other: &EdgeEnd<F>) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        if self.delta == other.delta {
            Ordering::Equal
        } else if (self.quadrant as usize) > (other.quadrant as usize) {
            Ordering::Greater
        } else if (self.quadrant as usize) < (other.quadrant as usize) {
            Ordering::Less
        } else {
            use crate::algorithm::kernels::{Kernel, Orientation, RobustKernel};
            // REVIEW:
            match RobustKernel::orient2d(other.coord_0, other.coord_1, self.coord_1) {
                Orientation::Clockwise => Ordering::Less,
                Orientation::CounterClockwise => Ordering::Greater,
                Orientation::Collinear => Ordering::Equal,
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_ord() {
        let fake_label = Label::empty_line();
        let edge_end_1 = EdgeEnd::new(
            Coordinate::zero(),
            Coordinate { x: 1.0, y: 1.0 },
            fake_label.clone(),
        );
        let edge_end_2 = EdgeEnd::new(
            Coordinate::zero(),
            Coordinate { x: 1.0, y: 1.0 },
            fake_label.clone(),
        );
        assert_eq!(edge_end_1.cmp(&edge_end_2), std::cmp::Ordering::Equal);

        // edge_end_3 is clockwise from edge_end_1
        let edge_end_3 = EdgeEnd::new(
            Coordinate::zero(),
            Coordinate { x: 1.0, y: -1.0 },
            fake_label.clone(),
        );
        assert_eq!(edge_end_1.cmp(&edge_end_3), std::cmp::Ordering::Less);
        assert_eq!(edge_end_3.cmp(&edge_end_1), std::cmp::Ordering::Greater);
    }
}
