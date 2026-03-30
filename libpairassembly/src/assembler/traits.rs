//! Boundary traits for assembler input/output interop.

use std::fmt::Display;

use crate::{Result, SequenceRead, errors::ConversionError};

/// Boundary trait for pair records accepted by the assembler API.
///
/// Implement this for external record types to use `Assembler` directly.
pub trait SeqRecordView {
    fn id(&self) -> &str;
    fn seq(&self) -> &str;
    fn qual(&self) -> &str;
}

/// Trait for constructing user-space record types from corrected output parts.
pub trait FromRecordParts: Sized {
    type Error;

    /// Construct a record instance from owned identifier, sequence, and quality parts.
    ///
    /// # Errors
    ///
    /// Returns an error if the target record type cannot be constructed from the
    /// provided parts.
    fn try_from_parts(
        id: String,
        seq: Vec<u8>,
        qual: Vec<u8>,
    ) -> std::result::Result<Self, Self::Error>;
}

/// Trait for extracting owned merged-record parts from terminal pipeline outputs.
pub trait IntoOwnedRecordParts {
    fn into_owned_record_parts(self) -> (String, Vec<u8>, Vec<u8>);
}

/// Trait for extracting owned paired-record parts from terminal pipeline outputs.
pub trait IntoOwnedPairRecordParts {
    fn into_owned_pair_record_parts(self) -> (String, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>);
}

/// Blanket conversion helper for merged outputs.
///
/// By convention, `into_*` methods in this crate are reserved for meaningful
/// conversion across representation boundaries. Identity-shaped `into_*`
/// endpoints are intentionally not provided.
pub trait IntoRecordConversion: IntoOwnedRecordParts {
    /// Convert merged terminal output into a user record type.
    ///
    /// # Errors
    ///
    /// Returns an error if the target record type cannot be constructed from
    /// the owned parts extracted from `self`.
    fn into_record<T>(self) -> Result<T>
    where
        Self: Sized,
        T: FromRecordParts,
        T::Error: Display,
    {
        let (id, seq, qual) = self.into_owned_record_parts();
        T::try_from_parts(id, seq, qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()).into())
    }
}

impl<T> IntoRecordConversion for T where T: IntoOwnedRecordParts {}

/// Blanket conversion helper for paired outputs.
///
/// By convention, `into_*` methods in this crate are reserved for meaningful
/// conversion across representation boundaries. Identity-shaped `into_*`
/// endpoints are intentionally not provided.
pub trait IntoRecordsConversion: IntoOwnedPairRecordParts {
    /// Convert paired terminal output into two user record values.
    ///
    /// # Errors
    ///
    /// Returns an error if either target record value cannot be constructed
    /// from the owned parts extracted from `self`.
    fn into_records<T>(self) -> Result<(T, T)>
    where
        Self: Sized,
        T: FromRecordParts,
        T::Error: Display,
    {
        let (id, fwd_seq, fwd_qual, rev_seq, rev_qual) = self.into_owned_pair_record_parts();
        let left = T::try_from_parts(id.clone(), fwd_seq, fwd_qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        let right = T::try_from_parts(id, rev_seq, rev_qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        Ok((left, right))
    }
}

impl<T> IntoRecordsConversion for T where T: IntoOwnedPairRecordParts {}

impl SeqRecordView for SequenceRead<'_> {
    fn id(&self) -> &str {
        self.id()
    }

    fn seq(&self) -> &str {
        self.sequence()
    }

    fn qual(&self) -> &str {
        self.quality_scores()
    }
}
