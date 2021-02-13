use super::{Location, Position, TopologyLocation};

use std::fmt;

/// A `Label` indicates the topological relationship of a component of a topology graph to a given
/// [Geometry].
///
/// This class supports labels for relationships to two `Geometry`s, which is sufficient for
/// algorithms for binary operations.
///
/// Topology graphs support the concept of labeling nodes and edges in the graph.  The label of a
/// node or edge specifies its topological relationship to one or more geometries.  (In fact, since
/// JTS operations have only two arguments labels are required for only two geometries).  A label
/// for a node or edge has one or two elements, depending on whether the node or edge occurs in one
/// or both of the input `Geometry`s.  Elements contain attributes which categorize the topological
/// location of the node or edge relative to the parent `Geometry`; that is, whether the node or
/// edge is in the interior, boundary or exterior of the `Geometry`.  Attributes have a value
/// from the set `{Interior, Boundary, Exterior}`.  In a node each element has  a single attribute
/// `(On)`.  For an edge each element has a triplet of attributes `(Left, On, Right)`.
///
/// It is up to the client code to associate the 0 and 1 [TopologyLocation]s with specific
/// geometries.
#[derive(Clone)]
pub(crate) struct Label {
    // REVIEW: better name? what does this stand for - "element location's topology"?
    elt: [TopologyLocation; 2],
}

impl fmt::Debug for Label {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Label {{ A: {:?}, B: {:?} }}",
            &self.elt[0], &self.elt[1]
        )
    }
}

impl Label {
    pub fn empty_line() -> Label {
        Label {
            elt: [
                TopologyLocation::empty_line(),
                TopologyLocation::empty_line(),
            ],
        }
    }

    pub fn empty_area() -> Self {
        Self {
            elt: [
                TopologyLocation::empty_area(),
                TopologyLocation::empty_area(),
            ],
        }
    }

    pub fn new(geom_index: usize, location: TopologyLocation) -> Self {
        let mut elt = match location {
            TopologyLocation::Line { .. } => [
                TopologyLocation::empty_line(),
                TopologyLocation::empty_line(),
            ],
            TopologyLocation::Area { .. } => [
                TopologyLocation::empty_area(),
                TopologyLocation::empty_area(),
            ],
        };
        elt[geom_index] = location;
        Label { elt }
    }

    pub fn flip(&mut self) {
        self.elt[0].flip();
        self.elt[1].flip();
    }

    pub fn location(&self, geom_index: usize, position: Position) -> Option<Location> {
        self.elt[geom_index].get(position)
    }

    pub fn on_location(&self, geom_index: usize) -> Option<Location> {
        self.elt[geom_index].get(Position::On)
    }

    pub fn set_location(&mut self, geom_index: usize, position: Position, location: Location) {
        self.elt[geom_index].set_location(position, location);
    }

    pub fn set_on_location(&mut self, geom_index: usize, location: Location) {
        self.elt[geom_index].set_location(Position::On, location);
    }

    pub fn set_all_locations(&mut self, geom_index: usize, location: Location) {
        self.elt[geom_index].set_all_locations(location)
    }

    pub fn set_all_locations_if_empty(&mut self, geom_index: usize, location: Location) {
        self.elt[geom_index].set_all_locations_if_empty(location)
    }

    pub fn geometry_count(&self) -> usize {
        self.elt
            .iter()
            .filter(|location| !location.is_empty())
            .count()
    }

    pub fn is_empty(&self, geom_index: usize) -> bool {
        self.elt[geom_index].is_empty()
    }

    pub fn is_any_empty(&self, geom_index: usize) -> bool {
        self.elt[geom_index].is_any_empty()
    }

    pub fn is_area(&self) -> bool {
        self.elt[0].is_area() || self.elt[1].is_area()
    }

    pub fn is_geom_area(&self, geom_index: usize) -> bool {
        self.elt[geom_index].is_area()
    }

    pub fn is_line(&self, geom_index: usize) -> bool {
        self.elt[geom_index].is_line()
    }
}
