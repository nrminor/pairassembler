mod bounds;
mod finder;
mod pair_overlap;
mod params;
mod slices;

pub use pair_overlap::PairOverlap;
pub use params::{OverlapParams, TiePolicy};

pub(crate) use bounds::{OverlapBounds, OverlapSpan};
pub(crate) use finder::OverlapFinder;
pub(crate) use slices::{AssemblyScratch, HasOrientedPairSlices, OrientedPairSlices};

pub(crate) mod private {
    pub(crate) trait Sealed {}
}
