//! Internal typestate markers for assembler DAG transitions.

use crate::overlap::{OrientedPairSlices, PairOverlap};

#[doc(hidden)]
pub trait OverlapStateStorage<'pair, 'scratch> {
    type Storage;
}

#[derive(Debug, Clone, Copy)]
pub struct OverlapUnsearched;
#[derive(Debug, Clone, Copy)]
pub struct OverlapFound;
#[derive(Debug, Clone, Copy)]
pub struct NoOverlapFound;

impl<'pair, 'scratch> OverlapStateStorage<'pair, 'scratch> for OverlapUnsearched {
    type Storage = OrientedPairSlices<'pair, 'scratch>;
}

impl<'pair, 'scratch> OverlapStateStorage<'pair, 'scratch> for OverlapFound {
    type Storage = PairOverlap<'pair, 'scratch>;
}

impl OverlapStateStorage<'_, '_> for NoOverlapFound {
    type Storage = ();
}

#[derive(Debug, Clone, Copy)]
pub struct Unvalidated;
#[derive(Debug, Clone, Copy)]
pub struct Validated;

#[derive(Debug, Clone, Copy)]
pub struct Unmerged;
#[derive(Debug, Clone, Copy)]
pub struct Merged;

#[derive(Debug, Clone, Copy)]
pub struct Uncorrected;
#[derive(Debug, Clone, Copy)]
pub struct Corrected;
