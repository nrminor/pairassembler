use crate::{
    Error, Result,
    assembler::{FromRecordParts, SeqRecordView},
    errors::ConversionError,
};

#[derive(Debug, Clone)]
#[repr(transparent)]
pub(crate) struct TupleRecord((String, String, String));

impl TupleRecord {
    pub(crate) fn try_new(id: String, seq: String, qual: String) -> Result<Self> {
        if seq.len() != qual.len() {
            return Err(crate::errors::SequenceQualityLengthMismatch(
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

impl FromRecordParts for TupleRecord {
    type Error = Error;

    fn try_from_parts(
        id: String,
        seq: Vec<u8>,
        qual: Vec<u8>,
    ) -> std::result::Result<Self, Self::Error> {
        let seq = String::from_utf8(seq)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        let qual = String::from_utf8(qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        Self::try_new(id, seq, qual)
    }
}
