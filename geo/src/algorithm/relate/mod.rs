mod edge_end_builder;
mod intersection_matrix;
pub mod relate_computer;
mod relate_node;

pub(crate) use edge_end_builder::EdgeEndBuilder;
pub(crate) use intersection_matrix::IntersectionMatrix;
pub(crate) use relate_node::RelateNodeFactory;

use crate::{GeoFloat, GeometryCow};

pub(crate) trait Relate<F, T> {
    fn relate(&self, other: &T) -> IntersectionMatrix;
}

impl<F: 'static + GeoFloat> Relate<F, GeometryCow<'_, F>> for GeometryCow<'_, F> {
    fn relate(&self, other: &GeometryCow<F>) -> IntersectionMatrix {
        let mut relate_computer = relate_computer::RelateComputer::new(self, other);
        relate_computer.compute_intersection_matrix()
    }
}
