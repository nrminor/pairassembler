use crate::{
    Result,
    errors::{ConversionError, InputOutputError, PairingError},
    prelude::utils::reverse_complement,
};

#[derive(Debug, Clone, Copy)]
pub struct SequenceRead<'read> {
    id: &'read str,
    seq: &'read str,
    qual: &'read str,
}

impl<'read> SequenceRead<'read> {
    #[doc(hidden)]
    #[must_use]
    pub const fn from_literal_parts(id: &'read str, seq: &'read str, qual: &'read str) -> Self {
        Self { id, seq, qual }
    }

    pub(crate) fn from_views(id: &'read str, seq: &'read str, qual: &'read str) -> Self {
        Self::new(id, seq, qual)
    }

    pub(crate) fn new(id: &'read str, seq: &'read str, qual: &'read str) -> Self {
        assert_eq!(seq.len(), qual.len());
        SequenceRead { id, seq, qual }
    }

    /// Construct a read after validating sequence and quality lengths match.
    ///
    /// # Errors
    ///
    /// Returns an error when `seq.len() != qual.len()`.
    pub fn try_new(id: &'read str, seq: &'read str, qual: &'read str) -> Result<Self> {
        if seq.len() != qual.len() {
            return Err(InputOutputError::SequenceQualityLengthMismatch(
                seq.to_string(),
                seq.len(),
                qual.to_string(),
                qual.len(),
            )
            .into());
        }

        Ok(SequenceRead { id, seq, qual })
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn reverse_complement(&self) -> Vec<u8> {
        let rc = reverse_complement(self.seq);
        rc.as_bytes().to_vec()
    }

    #[inline]
    #[must_use]
    pub fn check_for_mate(&self, possible_mate: &SequenceRead) -> bool {
        self.id == possible_mate.id
    }

    #[inline]
    #[must_use]
    pub fn id(&self) -> &'read str {
        self.id
    }

    #[inline]
    #[must_use]
    pub fn quality_scores(&self) -> &'read str {
        self.qual
    }

    #[inline]
    #[must_use]
    pub fn sequence(&self) -> &'read str {
        self.seq
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ReadPair<'mate> {
    fwd_mate: SequenceRead<'mate>,
    rev_mate: SequenceRead<'mate>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedSequenceRead {
    id: String,
    seq: String,
    qual: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedReadPair {
    id: String,
    fwd_seq: String,
    fwd_qual: String,
    rev_seq: String,
    rev_qual: String,
}

#[derive(Debug, Default)]
pub(crate) struct OwnedReadPairBuilder {
    id: Option<String>,
    fwd_seq: Option<Vec<u8>>,
    fwd_qual: Option<Vec<u8>>,
    rev_seq: Option<Vec<u8>>,
    rev_qual: Option<Vec<u8>>,
}

impl<'a> ReadPair<'a> {
    pub(crate) fn from_views(fwd_mate: SequenceRead<'a>, rev_mate: SequenceRead<'a>) -> Self {
        debug_assert_eq!(fwd_mate.id(), rev_mate.id());
        Self { fwd_mate, rev_mate }
    }

    /// Construct a read pair from two reads with matching identifiers.
    ///
    /// # Errors
    ///
    /// Returns an error when `read1.id != read2.id`.
    pub fn from(read1: SequenceRead<'a>, read2: SequenceRead<'a>) -> Result<Self> {
        if read1.id != read2.id {
            return Err(
                PairingError::UnmatchedIds(read1.id.to_string(), read2.id.to_string()).into(),
            );
        }
        let pair = ReadPair {
            fwd_mate: read1,
            rev_mate: read2,
        };
        Ok(pair)
    }

    #[inline]
    #[must_use]
    pub fn forward_read(&self) -> &SequenceRead<'a> {
        &self.fwd_mate
    }

    #[inline]
    #[must_use]
    pub fn reverse_read(&self) -> &SequenceRead<'a> {
        &self.rev_mate
    }

    #[inline]
    #[must_use]
    pub fn fwd_id(&self) -> &'a str {
        self.fwd_mate.id()
    }

    #[inline]
    #[must_use]
    pub fn rev_id(&self) -> &'a str {
        self.rev_mate.id()
    }

    #[inline]
    #[must_use]
    pub fn fwd_id_bytes(&self) -> &'a [u8] {
        self.fwd_mate.id().as_bytes()
    }

    #[inline]
    #[must_use]
    pub fn rev_id_bytes(&self) -> &'a [u8] {
        self.rev_mate.id().as_bytes()
    }

    #[inline]
    #[must_use]
    pub fn fwd_sequence(&self) -> &'a str {
        self.fwd_mate.sequence()
    }

    #[inline]
    #[must_use]
    pub fn rev_sequence(&self) -> &'a str {
        self.rev_mate.sequence()
    }

    #[inline]
    #[must_use]
    pub fn fwd_sequence_bytes(&self) -> &'a [u8] {
        self.fwd_mate.sequence().as_bytes()
    }

    #[inline]
    #[must_use]
    pub fn rev_sequence_bytes(&self) -> &'a [u8] {
        self.rev_mate.sequence().as_bytes()
    }

    #[inline]
    #[must_use]
    pub fn fwd_quality_scores(&self) -> &'a str {
        self.fwd_mate.quality_scores()
    }

    #[inline]
    #[must_use]
    pub fn rev_quality_scores(&self) -> &'a str {
        self.rev_mate.quality_scores()
    }

    #[inline]
    #[must_use]
    pub fn fwd_quality_bytes(&self) -> &'a [u8] {
        self.fwd_mate.quality_scores().as_bytes()
    }

    #[inline]
    #[must_use]
    pub fn rev_quality_bytes(&self) -> &'a [u8] {
        self.rev_mate.quality_scores().as_bytes()
    }
}

impl OwnedSequenceRead {
    pub(crate) fn try_from_ascii_bytes(id: String, seq: Vec<u8>, qual: Vec<u8>) -> Result<Self> {
        let seq = String::from_utf8(seq)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        let qual = String::from_utf8(qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        Ok(Self { id, seq, qual })
    }

    #[must_use]
    pub fn as_read(&self) -> SequenceRead<'_> {
        SequenceRead::from_views(&self.id, &self.seq, &self.qual)
    }

    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn sequence(&self) -> &str {
        &self.seq
    }

    #[must_use]
    pub fn sequence_bytes(&self) -> &[u8] {
        self.seq.as_bytes()
    }

    #[must_use]
    pub fn quality_scores(&self) -> &str {
        &self.qual
    }

    #[must_use]
    pub fn quality_bytes(&self) -> &[u8] {
        self.qual.as_bytes()
    }
}

impl OwnedReadPair {
    pub(crate) fn builder() -> OwnedReadPairBuilder {
        OwnedReadPairBuilder::default()
    }

    #[must_use]
    pub fn as_read_pair(&self) -> ReadPair<'_> {
        let fwd = SequenceRead::from_views(&self.id, &self.fwd_seq, &self.fwd_qual);
        let rev = SequenceRead::from_views(&self.id, &self.rev_seq, &self.rev_qual);
        ReadPair::from_views(fwd, rev)
    }

    #[must_use]
    pub fn into_reads(self) -> (OwnedSequenceRead, OwnedSequenceRead) {
        let fwd = OwnedSequenceRead {
            id: self.id.clone(),
            seq: self.fwd_seq,
            qual: self.fwd_qual,
        };
        let rev = OwnedSequenceRead {
            id: self.id,
            seq: self.rev_seq,
            qual: self.rev_qual,
        };
        (fwd, rev)
    }

    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn forward_sequence(&self) -> &str {
        &self.fwd_seq
    }

    #[must_use]
    pub fn forward_quality_scores(&self) -> &str {
        &self.fwd_qual
    }

    #[must_use]
    pub fn reverse_sequence(&self) -> &str {
        &self.rev_seq
    }

    #[must_use]
    pub fn reverse_quality_scores(&self) -> &str {
        &self.rev_qual
    }
}

impl OwnedReadPairBuilder {
    pub(crate) fn id(mut self, id: String) -> Self {
        self.id = Some(id);
        self
    }

    pub(crate) fn forward(mut self, seq: Vec<u8>, qual: Vec<u8>) -> Self {
        self.fwd_seq = Some(seq);
        self.fwd_qual = Some(qual);
        self
    }

    pub(crate) fn reverse(mut self, seq: Vec<u8>, qual: Vec<u8>) -> Self {
        self.rev_seq = Some(seq);
        self.rev_qual = Some(qual);
        self
    }

    pub(crate) fn build(self) -> Result<OwnedReadPair> {
        let id = Self::required(self.id, "read-pair id")?;
        let fwd_seq = Self::required(self.fwd_seq, "forward sequence")?;
        let fwd_qual = Self::required(self.fwd_qual, "forward quality")?;
        let rev_seq = Self::required(self.rev_seq, "reverse sequence")?;
        let rev_qual = Self::required(self.rev_qual, "reverse quality")?;

        let fwd_seq = String::from_utf8(fwd_seq)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        let fwd_qual = String::from_utf8(fwd_qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        let rev_seq = String::from_utf8(rev_seq)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        let rev_qual = String::from_utf8(rev_qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;

        Ok(OwnedReadPair {
            id,
            fwd_seq,
            fwd_qual,
            rev_seq,
            rev_qual,
        })
    }

    fn required<T>(value: Option<T>, name: &'static str) -> Result<T> {
        value.ok_or_else(|| {
            ConversionError::RecordConstruction(format!("missing {name} for owned read pair"))
                .into()
        })
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

    #[test]
    fn test_try_new_rejects_seq_qual_length_mismatch() {
        let result = SequenceRead::try_new("r1", "ACGT", "III");
        assert!(result.is_err());
    }

    #[test]
    fn test_read_mates_from_rejects_mismatched_ids() {
        let r1 = SequenceRead::new("read1", "ACGT", "IIII");
        let r2 = SequenceRead::new("read2", "ACGT", "IIII");
        let result = ReadPair::from(r1, r2);
        assert!(result.is_err());
    }

    #[test]
    fn test_read_mates_from_accepts_matching_ids() {
        let r1 = SequenceRead::new("read1", "ACGT", "IIII");
        let r2 = SequenceRead::new("read1", "TGCA", "IIII");
        let result = ReadPair::from(r1, r2);
        assert!(result.is_ok());
    }
}
