use std::result::Result as StdResult;

use crate::{
    Error, Result, assembler::SeqRecordView, errors::InputOutputError, read::OwnedSequenceRead,
};

#[derive(Debug, Clone)]
#[repr(transparent)]
pub(crate) struct TupleRecord((String, String, String));

impl TupleRecord {
    pub(crate) fn try_new(id: String, seq: String, qual: String) -> Result<Self> {
        if seq.len() != qual.len() {
            return Err(InputOutputError::SequenceQualityLengthMismatch(
                seq.clone(),
                seq.len(),
                qual.clone(),
                qual.len(),
            )
            .into());
        }
        Ok(Self((id, seq, qual)))
    }

    pub(crate) fn from_strs(id: &str, seq: &str, qual: &str) -> Self {
        match Self::try_new(id.to_string(), seq.to_string(), qual.to_string()) {
            Ok(record) => record,
            Err(err) => panic!("test fixture tuple record should be valid: {err}"),
        }
    }

    pub(crate) fn id(&self) -> &str {
        &self.0.0
    }

    pub(crate) fn seq(&self) -> &str {
        &self.0.1
    }

    pub(crate) fn qual(&self) -> &str {
        &self.0.2
    }
}

impl SeqRecordView for TupleRecord {
    fn id(&self) -> &str {
        self.id()
    }

    fn seq(&self) -> &str {
        self.seq()
    }

    fn qual(&self) -> &str {
        self.qual()
    }
}

impl TryFrom<OwnedSequenceRead> for TupleRecord {
    type Error = Error;

    fn try_from(read: OwnedSequenceRead) -> StdResult<Self, Self::Error> {
        Self::try_new(
            read.id().to_string(),
            read.sequence().to_string(),
            read.quality_scores().to_string(),
        )
    }
}
