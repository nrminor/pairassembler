//! Assembler-centered public API with two usage layers.
//!
//! - Layer A: per-pair fluent DAG transitions (`on_pair(...).overlap()...`).
//! - Layer B: collection orchestration (`process_iter`, `process_iter_with`).
//!
//! The default convenience path (`process_pair` / `process_iter`) is the checked
//! path (`overlap -> validate -> merge -> correct`). Unchecked expert paths are
//! available from [`OverlapContext`] (`merge_unchecked`,
//! `correct_pair_unchecked`) and are intentionally explicit.
//!
//! Transition channels are tracked across four dimensions:
//! - `O`: overlap discovered (`NoOverlap`/`HasOverlap`)
//! - `V`: overlap validated (`Unvalidated`/`Validated`)
//! - `M`: merged state (`Unmerged`/`Merged`)
//! - `C`: correction state (`Uncorrected`/`Corrected`)
//!
//! ```

mod capability;
mod config;
mod context;
mod input;
mod ops;
mod process_iter;
mod traits;
mod typestate;

pub use config::{
    Assembler, AssemblerBuilder, AssemblerConfig, BatchPolicy, ExecutionPolicy, MergeParams,
};
pub use context::{
    MergeContext, MergedContext, OverlapContext, PairContext, PairReady, ValidatedContext,
    ValidatedMergedContext,
};
pub use input::PairInput;
pub use process_iter::ProcessIter;
pub use traits::{
    FromRecordParts, IntoOwnedPairRecordParts, IntoOwnedRecordParts, IntoRecordConversion,
    IntoRecordsConversion, SeqRecordView,
};
pub use typestate::{
    Corrected, HasOverlap, Merged, NoOverlap, Uncorrected, Unmerged, Unvalidated, Validated,
};

pub(crate) use capability::*;

#[cfg(test)]
pub(crate) use crate::{BaseCallValidator, OverlapParams, TiePolicy};

#[cfg(test)]
mod tests;
