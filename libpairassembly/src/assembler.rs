//! Assembler-centered public API with two usage layers.
//!
//! - Layer A: per-pair fluent DAG transitions (`on_pair(...).find_overlap()...`).
//! - Layer B: collection orchestration (`process_iter`, `process_iter_with`).
//!
//! The default convenience path (`process_pair` / `process_iter`) is the checked
//! found-overlap path (`find_overlap -> validate -> merge -> correct`). Expert
//! paths can reorder the same fluent transitions after the overlap search branch;
//! the receiver typestate records whether each merge consumed validated or
//! unvalidated slices.
//!
//! Transition channels are tracked across four dimensions:
//! - `O`: overlap search state (`OverlapUnsearched`/`OverlapFound`/`NoOverlapFound`)
//! - `V`: overlap validated (`Unvalidated`/`Validated`)
//! - `M`: merged state (`Unmerged`/`Merged`)
//! - `C`: correction state (`Uncorrected`/`Corrected`)
//!

mod capability;
mod config;
mod context;
mod input;
mod ops;
mod process_iter;
mod traits;
mod typestate;

pub use crate::merge::{MergeParams, MergeTiePolicy};
pub use config::{Assembler, AssemblerBuilder, AssemblerConfig};
pub use context::{
    CorrectedContext, CorrectedMergedContext, MergedContext, NoOverlapContext, OverlapContext,
    OverlapOutcome, OverlapSearch, PairReady, ValidatedContext, ValidatedCorrectedContext,
    ValidatedCorrectedMergedContext, ValidatedMergedContext,
};
pub use input::PairInput;
pub use process_iter::ProcessIter;
pub use traits::SeqRecordView;
pub use typestate::{
    Corrected, Merged, NoOverlapFound, OverlapFound, OverlapUnsearched, Uncorrected, Unmerged,
    Unvalidated, Validated,
};

pub(crate) use capability::*;

#[cfg(test)]
pub(crate) use crate::{OverlapParams, OverlapValidator, TiePolicy};

#[cfg(test)]
mod tests;
