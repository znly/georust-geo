use crate::algorithm::dimensions::Dimensions;
use crate::geomgraph::Location;

/// Models a *Dimensionally Extended Nine-Intersection Model (DE-9IM)* matrix.
///
/// DE-9IM matrix values (such as "212FF1FF2") specify the topological relationship between
/// two [Geometeries](struct.Geometry.html).
///
/// TODO: This class can also represent matrix patterns (such as "T*T******")
/// which are used for matching instances of DE-9IM matrices.
///
/// DE-9IM matrices are 3x3 matrices with integer entries.
/// The matrix indices {0,1,2} represent the topological locations
/// that occur in a geometry (Interior, Boundary, Exterior).
/// These are provided by the enum cases
/// [Location::Interior, Location::Boundary, Location::Exterior](enum.Location.html).
///
/// The matrix entries represent the [Dimensions](enum.Dimension.html) of each intersection.
///
/// For a description of the DE-9IM and the spatial predicates derived from it,
/// see the following references:
/// - [OGC 99-049 OpenGIS Simple Features Specification for SQL](http://portal.opengeospatial.org/files/?artifact_id=829), Section 2.1.13
/// - [OGC 06-103r4 OpenGIS Implementation Standard for Geographic information - Simple feature access - Part 1: Common architecture](http://portal.opengeospatial.org/files/?artifact_id=25355), Section 6.1.15 (which provides some further details on certain predicate specifications).
/// - Wikipedia article on [DE-9IM](https://en.wikipedia.org/wiki/DE-9IM)
///
/// This implementation is heavily based on that from the [JTS project](https://github.com/locationtech/jts/blob/master/modules/core/src/main/java/org/locationtech/jts/geom/IntersectionMatrix.java).
#[derive(PartialEq, Eq)]
pub(crate) struct IntersectionMatrix([[Dimensions; 3]; 3]);

impl std::fmt::Debug for IntersectionMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        fn char_for_dim(dim: &Dimensions) -> &'static str {
            match dim {
                Dimensions::Empty => "F",
                Dimensions::ZeroDimensional => "0",
                Dimensions::OneDimensional => "1",
                Dimensions::TwoDimensional => "2",
            }
        }
        let text = self
            .0
            .iter()
            .flat_map(|r| r.iter().map(char_for_dim))
            .collect::<Vec<&str>>()
            .join("");

        write!(f, "IntersectionMatrix({})", &text)
    }
}

impl IntersectionMatrix {
    pub fn empty() -> Self {
        IntersectionMatrix([[Dimensions::Empty; 3]; 3])
    }

    // CLEANUP: remove?
    #[allow(dead_code)]
    pub fn new(dimensions: [[Dimensions; 3]; 3]) -> Self {
        IntersectionMatrix(dimensions)
    }

    #[cfg(test)]
    pub fn from_str(str: &str) -> Self {
        let mut im = IntersectionMatrix::empty();
        im.set_at_least_from_string(str);
        im
    }

    pub fn set(&mut self, location_a: Location, location_b: Location, dimensionality: Dimensions) {
        self.0[location_a as usize][location_b as usize] = dimensionality;
    }

    pub fn set_at_least(
        &mut self,
        location_a: Location,
        location_b: Location,
        minimum_dimension_value: Dimensions,
    ) {
        if self.0[location_a as usize][location_b as usize] < minimum_dimension_value {
            self.0[location_a as usize][location_b as usize] = minimum_dimension_value;
        }
    }

    pub fn set_at_least_if_valid(
        &mut self,
        location_a: Option<Location>,
        location_b: Option<Location>,
        minimum_dimension_value: Dimensions,
    ) {
        if let Some(location_a) = location_a {
            if let Some(location_b) = location_b {
                self.set_at_least(location_a, location_b, minimum_dimension_value);
            }
        }
    }

    pub fn set_at_least_from_string(&mut self, dimensions: &str) {
        if dimensions.len() != 9 {
            todo!("return proper error, or better yet make this a compile time macro")
        }

        let mut i = 0;
        for c in dimensions.chars() {
            let a = i / 3;
            let b = i % 3;
            i += 1;
            match c {
                '0' => self.0[a][b] = self.0[a][b].max(Dimensions::ZeroDimensional),
                '1' => self.0[a][b] = self.0[a][b].max(Dimensions::OneDimensional),
                '2' => self.0[a][b] = self.0[a][b].max(Dimensions::TwoDimensional),
                'F' => {}
                _ => todo!("return proper error, or better yet make this a compile time macro"),
            }
        }
    }

    pub fn is_contains(&self) -> bool {
        self.0[Location::Interior as usize][Location::Interior as usize] != Dimensions::Empty
            && self.0[Location::Exterior as usize][Location::Interior as usize] == Dimensions::Empty
            && self.0[Location::Exterior as usize][Location::Boundary as usize] == Dimensions::Empty
    }
}
