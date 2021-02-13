pub trait BoundaryNodeRule {
    fn is_in_boundary(&self, boundary_count: usize) -> bool;
}

/// A [BoundaryNodeRule]() specifies that points are in the boundary of a
/// lineal geometry iff the point lies on the boundary of an odd number of
/// components.
///
/// Under this rule LinearRings and closed LineStrings have an empty boundary.
///
/// This is the rule specified by the _OGC SFS_, and is the default rule
pub(crate) struct Mod2BoundaryNodeRule;
impl BoundaryNodeRule for Mod2BoundaryNodeRule {
    fn is_in_boundary(&self, boundary_count: usize) -> bool {
        // the "Mod-2 Rule"
        return boundary_count % 2 == 1;
    }
}
