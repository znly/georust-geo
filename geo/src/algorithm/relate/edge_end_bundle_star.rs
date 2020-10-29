use crate::geomgraph::{Float, GeometryGraph};

// JTS: /**
// JTS:  * An ordered list of {@link EdgeEndBundle}s around a {@link RelateNode}.
// JTS:  * They are maintained in CCW order (starting with the positive x-axis) around the node
// JTS:  * for efficient lookup and topology building.
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class EdgeEndBundleStar
// JTS:   extends EdgeEndStar
// JTS: {
// JTS:   /**
// JTS:    * Creates a new empty EdgeEndBundleStar
// JTS:    */
// JTS:   public EdgeEndBundleStar() {
// JTS:   }
// JTS:
pub(crate) struct EdgeEndBundleStar<F>
where
    F: Float,
{
    _marker: std::marker::PhantomData<F>,
}

impl<F> EdgeEndBundleStar<F>
where
    F: Float,
{
    pub(crate) fn new() -> Self {
        EdgeEndBundleStar {
            _marker: std::marker::PhantomData,
        }
    }

    pub(crate) fn compute_labeling(&self, graph_a: &GeometryGraph<F>, graph_b: &GeometryGraph<F>) {
        todo!()
    }
}
// JTS:   /**
// JTS:    * Insert a EdgeEnd in order in the list.
// JTS:    * If there is an existing EdgeStubBundle which is parallel, the EdgeEnd is
// JTS:    * added to the bundle.  Otherwise, a new EdgeEndBundle is created
// JTS:    * to contain the EdgeEnd.
// JTS:    * <br>
// JTS:    */
// JTS:   public void insert(EdgeEnd e)
// JTS:   {
// JTS:     EdgeEndBundle eb = (EdgeEndBundle) edgeMap.get(e);
// JTS:     if (eb == null) {
// JTS:       eb = new EdgeEndBundle(e);
// JTS:       insertEdgeEnd(e, eb);
// JTS:     }
// JTS:     else {
// JTS:       eb.insert(e);
// JTS:     }
// JTS:   }
// JTS:
// JTS:   /**
// JTS:    * Update the IM with the contribution for the EdgeStubs around the node.
// JTS:    */
// JTS:   void updateIM(IntersectionMatrix im)
// JTS:   {
// JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
// JTS:       EdgeEndBundle esb = (EdgeEndBundle) it.next();
// JTS:       esb.updateIM(im);
// JTS:     }
// JTS:   }
// JTS:
// JTS: }
