use crate::geomgraph::Location;

// CLEANUP: make internal private?
pub struct IntersectionMatrix(pub [[usize; 3]; 3]);

impl IntersectionMatrix {
    pub fn new() -> IntersectionMatrix {
        IntersectionMatrix([[0; 3]; 3])
    }

    pub fn set(&mut self, arg_1: Location, arg_2: Location, dimensionality: usize) {
        self.0[arg_1 as usize][arg_2 as usize] = dimensionality;
    }
}
