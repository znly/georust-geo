use super::{EdgeEnd, Location, Node, NodeFactory};
use crate::{Coordinate, GeoFloat};

use std::collections::BTreeMap;
use std::fmt;
use std::marker::PhantomData;

/// A map of nodes, indexed by the coordinate of the node
pub(crate) struct NodeMap<F, NF>
where
    F: GeoFloat,
    NF: NodeFactory<F>,
{
    map: BTreeMap<NodeKey<F>, (Node<F>, NF::Edges)>,
    _node_factory: PhantomData<NF>,
}

impl<F, NF> fmt::Debug for NodeMap<F, NF>
where
    F: GeoFloat,
    NF: NodeFactory<F>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("NodeMap")
            .field("map.len()", &self.map.len())
            .finish()
    }
}

#[derive(Clone)]
struct NodeKey<F: GeoFloat>(Coordinate<F>);

impl<F: GeoFloat> std::cmp::Ord for NodeKey<F> {
    fn cmp(&self, other: &NodeKey<F>) -> std::cmp::Ordering {
        // TODO: BTree requires Ord - can we guarantee all coords in the graph are non-NaN?
        debug_assert!(!self.0.x.is_nan());
        debug_assert!(!self.0.y.is_nan());
        debug_assert!(!other.0.x.is_nan());
        debug_assert!(!other.0.y.is_nan());

        crate::utils::lex_cmp(&self.0, &other.0)
    }
}

impl<F: GeoFloat> std::cmp::PartialOrd for NodeKey<F> {
    fn partial_cmp(&self, other: &NodeKey<F>) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<F: GeoFloat> std::cmp::PartialEq for NodeKey<F> {
    fn eq(&self, other: &NodeKey<F>) -> bool {
        // TODO: BTree requires Eq - can we guarantee all coords in the graph are non-NaN?
        debug_assert!(!self.0.x.is_nan());
        debug_assert!(!self.0.y.is_nan());
        debug_assert!(!other.0.x.is_nan());
        debug_assert!(!other.0.y.is_nan());
        return self.0 == other.0;
    }
}

impl<F: GeoFloat> std::cmp::Eq for NodeKey<F> {}

impl<F, NF> NodeMap<F, NF>
where
    F: GeoFloat,
    NF: NodeFactory<F>,
{
    pub fn new() -> Self {
        NodeMap {
            map: BTreeMap::new(),
            _node_factory: PhantomData,
        }
    }
    pub fn add_node_with_coordinate(&mut self, coord: Coordinate<F>) -> &mut (Node<F>, NF::Edges) {
        let key = NodeKey(coord);
        self.map.entry(key).or_insert(NF::create_node(coord))
    }

    /// returns the node if found
    pub fn find(&self, coord: Coordinate<F>) -> Option<&(Node<F>, NF::Edges)> {
        self.map.get(&NodeKey(coord))
    }

    pub fn iter(&self) -> impl Iterator<Item = &(Node<F>, NF::Edges)> {
        self.map.values()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut (Node<F>, NF::Edges)> {
        self.map.values_mut()
    }

    pub fn boundary_nodes(&self, geom_index: usize) -> impl Iterator<Item = &(Node<F>, NF::Edges)> {
        self.map.values().filter(move |(node, _edges)| {
            matches!(
                node.label().on_location(geom_index),
                Some(Location::Boundary)
            )
        })
    }
}
