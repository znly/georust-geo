/// Utility functions for working with quadrants, which are labeled as follows:
///          (+)
///        NW ┃ NE
///    (-) ━━━╋━━━━ (+)
///        SW ┃ SE
///          (-)
// CLEANUP: can we remove explicit discriminant? It's used in
#[derive(Debug, Clone, Copy)]
pub enum Quadrant {
    NE = 0,
    NW = 1,
    SW = 2,
    SE = 3,
}

impl Quadrant {
    pub fn new<F>(dx: F, dy: F) -> Quadrant
    where
        F: num_traits::Float,
    {
        if dx.is_zero() && dy.is_zero() {
            todo!("gracefully handle non-quadrant")
        }
        if dx >= F::zero() {
            if dy >= F::zero() {
                Quadrant::NE
            } else {
                Quadrant::SE
            }
        } else {
            if dy >= F::zero() {
                Quadrant::NW
            } else {
                Quadrant::SW
            }
        }
    }
}
