//! Internal typestate markers for assembler DAG transitions.

use crate::overlap::PairOverlap;

#[doc(hidden)]
pub trait OverlapStateStorage<'pair> {
    type Storage;
}

#[derive(Debug, Clone, Copy)]
pub struct OverlapUnsearched;
#[derive(Debug, Clone, Copy)]
pub struct OverlapFound;
#[derive(Debug, Clone, Copy)]
pub struct NoOverlapFound;

impl OverlapStateStorage<'_> for OverlapUnsearched {
    type Storage = ();
}

impl<'pair> OverlapStateStorage<'pair> for OverlapFound {
    type Storage = PairOverlap<'pair>;
}

impl OverlapStateStorage<'_> for NoOverlapFound {
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
