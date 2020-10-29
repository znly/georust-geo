use crate::algorithm::dimensions::Dimensions;
use crate::geomgraph::Location;

// CLEANUP: make internal private?
pub(crate) struct IntersectionMatrix([[Dimensions; 3]; 3]);

impl IntersectionMatrix {
    pub fn new() -> IntersectionMatrix {
        IntersectionMatrix([[Dimensions::Empty; 3]; 3])
    }

    pub fn set_at_least(&mut self, dimensions: &str) {
        if dimensions.len() != 9 {
            todo!("return proper error, or better yet make this a compile time macro")
        }

        let i = 0;
        for c in dimensions.chars() {
            let a = i / 3;
            let b = i % 3;
            match c {
                '0' => self.0[a][b] = self.0[a][b].max(Dimensions::ZeroDimensional),
                '1' => self.0[a][b] = self.0[a][b].max(Dimensions::OneDimensional),
                '2' => self.0[a][b] = self.0[a][b].max(Dimensions::TwoDimensional),
                'F' => {}
                _ => todo!("return proper error, or better yet make this a compile time macro"),
            }
        }
    }

    pub fn set(&mut self, arg_1: Location, arg_2: Location, dimensionality: Dimensions) {
        self.0[arg_1 as usize][arg_2 as usize] = dimensionality;
    }
}
