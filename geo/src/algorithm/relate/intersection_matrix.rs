use crate::algorithm::dimensions::Dimensions;
use crate::geomgraph::Location;

// CLEANUP: make internal private?
pub struct IntersectionMatrix(pub [[Dimensions; 3]; 3]);

impl IntersectionMatrix {
    pub fn new() -> IntersectionMatrix {
        IntersectionMatrix([[Dimensions::Empty; 3]; 3])
    }

    pub fn set(&mut self, arg_1: Location, arg_2: Location, dimensionality: Dimensions) {
        self.0[arg_1 as usize][arg_2 as usize] = dimensionality;
    }
}
