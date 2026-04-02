//! Internal capability traits for post-overlap operation contracts.

use crate::{
    PairOverlap, ReadPair, Result,
    correct::{CorrectedMergedRead, CorrectedReadPair, CorrectionWindow},
    errors::OverlapError,
    merge::{MergeView, MergedRead},
    validate::{ValidatedOverlap, ValidationMetrics},
};
use std::borrow::Cow;

use super::{
    PairContext,
    context::{OverlapOutcome, OverlapSnapshot},
};

pub(crate) mod private {
    pub(crate) trait Sealed {}
}

/// Internal marker trait for state/output carriers participating in assembler DAG operations.
pub(crate) trait PairState: private::Sealed {}

/// Capability for materializing overlap evidence.
pub(crate) trait HasPairOverlap: PairState {
    fn materialize_pair_overlap(&self) -> Result<PairOverlap<'_>>;
}

/// Capability for borrowing source read-pair evidence.
pub(crate) trait HasReadPair: PairState {
    fn read_pair(&self) -> &ReadPair<'_>;
}

/// Capability for exposing normalized merge-ready overlap views.
pub(crate) trait HasMergeableOverlap: HasReadPair + HasPairOverlap {
    fn merge_view(&self) -> Result<MergeView<'_>>;
}

impl<R, O, V, M, C> HasMergeableOverlap for PairContext<'_, '_, R, O, V, M, C> {
    fn merge_view(&self) -> Result<MergeView<'_>> {
        let pair = self.read_pair_ref();
        let snapshot = match self.overlap_outcome() {
            OverlapOutcome::Found(snapshot) => snapshot,
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                return Err(OverlapError::NoOverlapFound.into());
            },
        };

        build_merge_view_from_snapshot(pair, snapshot)
    }
}

impl HasMergeableOverlap for ValidatedOverlap<'_> {
    fn merge_view(&self) -> Result<MergeView<'_>> {
        let pair = self.read_pair();
        let overlap = self.overlap();
        MergeView::from_pair_bounds(
            pair,
            overlap.len(),
            overlap.forward_start_offset(),
            overlap.forward_end_offset(),
            overlap.reverse_start_offset(),
            overlap.reverse_end_offset(),
        )
    }
}

fn build_merge_view_from_snapshot<'a>(
    pair: &'a ReadPair<'a>,
    snapshot: OverlapSnapshot,
) -> Result<MergeView<'a>> {
    MergeView::from_pair_bounds(
        pair,
        snapshot.overlap_len(),
        snapshot.fwd_start_offset(),
        snapshot.fwd_end_offset(),
        snapshot.rev_start_offset(),
        snapshot.rev_end_offset(),
    )
}

/// Capability for exposing an aligned overlap-local correction window.
pub(crate) trait HasCorrectionWindow: PairState {
    fn correction_window(&self) -> Result<CorrectionWindow<'_>>;
}

/// Capability for exposing consensus record payload.
pub(crate) trait HasConsensusRecord: PairState {
    fn consensus_id(&self) -> &str;
    fn consensus_seq(&self) -> &[u8];
    fn consensus_qual(&self) -> &[u8];
}

/// Capability for exposing retained validation-stage metrics.
pub(crate) trait HasValidationMetrics: PairState {
    fn validation_metrics(&self) -> &ValidationMetrics;
}

impl<R, O, V, M, C> private::Sealed for PairContext<'_, '_, R, O, V, M, C> {}
impl<R, O, V, M, C> PairState for PairContext<'_, '_, R, O, V, M, C> {}

impl private::Sealed for ValidatedOverlap<'_> {}
impl PairState for ValidatedOverlap<'_> {}

impl private::Sealed for MergedRead {}
impl PairState for MergedRead {}

impl private::Sealed for CorrectedMergedRead {}
impl PairState for CorrectedMergedRead {}

impl private::Sealed for CorrectedReadPair {}
impl PairState for CorrectedReadPair {}

impl<R, O, V, M, C> HasReadPair for PairContext<'_, '_, R, O, V, M, C> {
    fn read_pair(&self) -> &ReadPair<'_> {
        self.read_pair_ref()
    }
}

impl<R, O, V, M, C> HasPairOverlap for PairContext<'_, '_, R, O, V, M, C> {
    fn materialize_pair_overlap(&self) -> Result<PairOverlap<'_>> {
        match self.overlap_outcome() {
            OverlapOutcome::Found(snapshot) => {
                Ok(snapshot.materialize_overlap(self.read_pair_ref()))
            },
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                Err(OverlapError::NoOverlapFound.into())
            },
        }
    }
}

impl HasReadPair for ValidatedOverlap<'_> {
    fn read_pair(&self) -> &ReadPair<'_> {
        ValidatedOverlap::read_pair(self)
    }
}

impl HasPairOverlap for ValidatedOverlap<'_> {
    fn materialize_pair_overlap(&self) -> Result<PairOverlap<'_>> {
        let overlap = self.overlap();
        Ok(PairOverlap::try_new(
            overlap.len(),
            overlap.forward_start_offset(),
            overlap.forward_end_offset(),
            overlap.reverse_start_offset(),
            overlap.reverse_end_offset(),
            overlap.forward_sequence(),
            overlap.forward_qualities(),
            overlap.reverse_sequence().to_vec(),
            overlap.reverse_qualities().to_vec(),
        )
        .expect("validated overlaps should always retain structurally valid overlap windows"))
    }
}

impl HasConsensusRecord for MergedRead {
    fn consensus_id(&self) -> &str {
        self.id()
    }

    fn consensus_seq(&self) -> &[u8] {
        self.sequence()
    }

    fn consensus_qual(&self) -> &[u8] {
        self.qualities()
    }
}

impl HasCorrectionWindow for MergedRead {
    fn correction_window(&self) -> Result<CorrectionWindow<'_>> {
        Ok(CorrectionWindow::new(
            Cow::Borrowed(self.provenance().fwd_overlap_seq()),
            Cow::Borrowed(self.provenance().fwd_overlap_qual()),
            Cow::Borrowed(self.provenance().rev_overlap_seq()),
            Cow::Borrowed(self.provenance().rev_overlap_qual()),
        ))
    }
}

impl<'asm, 'pair, R, V> HasCorrectionWindow
    for PairContext<
        'asm,
        'pair,
        R,
        super::typestate::HasOverlap,
        V,
        super::typestate::Unmerged,
        super::typestate::Uncorrected,
    >
where
    PairContext<
        'asm,
        'pair,
        R,
        super::typestate::HasOverlap,
        V,
        super::typestate::Unmerged,
        super::typestate::Uncorrected,
    >: HasReadPair + HasPairOverlap,
{
    fn correction_window(&self) -> Result<CorrectionWindow<'_>> {
        let overlap = self.materialize_pair_overlap()?;
        Ok(CorrectionWindow::new(
            Cow::Owned(overlap.forward_sequence().to_vec()),
            Cow::Owned(overlap.forward_qualities().to_vec()),
            Cow::Owned(overlap.reverse_sequence().to_vec()),
            Cow::Owned(overlap.reverse_qualities().to_vec()),
        ))
    }
}

impl HasConsensusRecord for CorrectedMergedRead {
    fn consensus_id(&self) -> &str {
        self.id()
    }

    fn consensus_seq(&self) -> &[u8] {
        self.sequence_bytes()
    }

    fn consensus_qual(&self) -> &[u8] {
        self.quality_bytes()
    }
}

impl<R, M, C> HasValidationMetrics
    for PairContext<'_, '_, R, super::typestate::HasOverlap, super::typestate::Validated, M, C>
{
    fn validation_metrics(&self) -> &ValidationMetrics {
        self.validation_metrics_ref()
            .expect("validated contexts must retain validation metrics")
    }
}

impl HasValidationMetrics for ValidatedOverlap<'_> {
    fn validation_metrics(&self) -> &ValidationMetrics {
        ValidatedOverlap::validation_metrics(self)
    }
}
