// NOTE: JTS implements a EdgeEndStar base class, but the inheritance is hard to model vis a vis
// sharing references in particular, so since the only child class we use for the Relate operation
// is the EdgeEndBundleStar, we shoehorn all functionality into that child.
//
// use super::{EdgeEnd, Float, GeometryGraph};
//
// // TODO: delegate to subclass adapter (or maybe as a trait?)
// #[derive(Clone)]
// pub(crate) struct EdgeEndStar<F>(std::marker::PhantomData<F>)
// where
//     F: Float;
//
// impl<F> EdgeEndStar<F>
// where
//     F: Float,
// {
//     pub(crate) fn insert(&mut self, edge_end: EdgeEnd<F>) {
//         todo!()
//     }
//
//     pub(crate) fn compute_labeling(&self, graph_a: &GeometryGraph<F>, graph_b: &GeometryGraph<F>) {
//         todo!()
//     }
//
// }
