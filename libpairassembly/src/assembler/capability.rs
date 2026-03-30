//! Internal capability traits for post-overlap operation contracts.

use crate::{
    PairOverlap, ReadPair, Result,
    correct::{CorrectedMergedRead, CorrectedReadPair},
    errors::OverlapError,
    merge::UncorrectedMergedRead,
    validate::ValidatedOverlap,
};

use super::{PairContext, context::OverlapOutcome};

pub(crate) mod private {
    pub(crate) trait Sealed {}
}

/// Internal marker trait for state/output carriers participating in assembler DAG operations.
pub(crate) trait PairState: private::Sealed {}

/// Optional validation diagnostics cache for state carriers.
#[derive(Debug, Clone)]
pub(crate) struct ValidationDiag {
    pub(crate) min_overlap_len: usize,
    pub(crate) observed_error_rate: f32,
    pub(crate) maximum_expected_error_rate: Option<f32>,
}

/// Capability for materializing overlap evidence.
pub(crate) trait HasPairOverlap: PairState {
    fn materialize_pair_overlap(&self) -> Result<PairOverlap<'_>>;
}

/// Capability for borrowing source read-pair evidence.
pub(crate) trait HasReadPair: PairState {
    fn read_pair(&self) -> &ReadPair<'_>;
}

/// Capability for exposing merged-read correction evidence.
pub(crate) trait HasCorrectionEvidence: PairState {
    fn forward_source_seq(&self) -> &[u8];
    fn forward_source_qual(&self) -> &[u8];
    fn reverse_source_seq(&self) -> &[u8];
    fn reverse_source_qual(&self) -> &[u8];
}

/// Capability for exposing consensus record payload.
pub(crate) trait HasConsensusRecord: PairState {
    fn consensus_id(&self) -> &str;
    fn consensus_seq(&self) -> &[u8];
    fn consensus_qual(&self) -> &[u8];
}

/// Capability for optional validation diagnostics.
pub(crate) trait HasValidationDiag: PairState {
    fn validation_diag(&self) -> Option<&ValidationDiag>;
}

impl<'asm, 'pair, R, O, V, M, C> private::Sealed for PairContext<'asm, 'pair, R, O, V, M, C> {}
impl<'asm, 'pair, R, O, V, M, C> PairState for PairContext<'asm, 'pair, R, O, V, M, C> {}

impl<'a> private::Sealed for ValidatedOverlap<'a> {}
impl<'a> PairState for ValidatedOverlap<'a> {}

impl private::Sealed for UncorrectedMergedRead {}
impl PairState for UncorrectedMergedRead {}

impl private::Sealed for CorrectedMergedRead {}
impl PairState for CorrectedMergedRead {}

impl private::Sealed for CorrectedReadPair {}
impl PairState for CorrectedReadPair {}

impl<'asm, 'pair, R, O, V, M, C> HasReadPair for PairContext<'asm, 'pair, R, O, V, M, C> {
    fn read_pair(&self) -> &ReadPair<'_> {
        self.read_pair_ref()
    }
}

impl<'asm, 'pair, R, O, V, M, C> HasPairOverlap for PairContext<'asm, 'pair, R, O, V, M, C> {
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

impl<'pair> HasReadPair for ValidatedOverlap<'pair> {
    fn read_pair(&self) -> &ReadPair<'_> {
        self.read_pair()
    }
}

impl<'pair> HasPairOverlap for ValidatedOverlap<'pair> {
    fn materialize_pair_overlap(&self) -> Result<PairOverlap<'_>> {
        let overlap = self.overlap();
        Ok(PairOverlap::from_components(
            overlap.len(),
            overlap.forward_start_offset(),
            overlap.forward_end_offset(),
            overlap.reverse_start_offset(),
            overlap.reverse_end_offset(),
            overlap.forward_sequence(),
            overlap.forward_qualities(),
            overlap.reverse_sequence().to_vec(),
            overlap.reverse_qualities().to_vec(),
        ))
    }
}

impl HasConsensusRecord for UncorrectedMergedRead {
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

impl HasCorrectionEvidence for UncorrectedMergedRead {
    fn forward_source_seq(&self) -> &[u8] {
        self.forward_source_seq()
    }

    fn forward_source_qual(&self) -> &[u8] {
        self.forward_source_qual()
    }

    fn reverse_source_seq(&self) -> &[u8] {
        self.reverse_source_seq()
    }

    fn reverse_source_qual(&self) -> &[u8] {
        self.reverse_source_qual()
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

impl<'asm, 'pair, R, O, V, M, C> HasValidationDiag for PairContext<'asm, 'pair, R, O, V, M, C> {
    fn validation_diag(&self) -> Option<&ValidationDiag> {
        None
    }
}

impl<'a> HasValidationDiag for ValidatedOverlap<'a> {
    fn validation_diag(&self) -> Option<&ValidationDiag> {
        None
    }
}
