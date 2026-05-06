use std::str;

use color_eyre::eyre::{Result, WrapErr};
use libpairassembly::SeqRecordView;
use noodles::fastq::Record as FastqRecord;

#[derive(Clone, Copy)]
pub(super) struct FastqMateKey<'key>(&'key [u8]);

impl<'key> FastqMateKey<'key> {
    pub(super) fn from_record(record: &'key FastqRecord) -> Self {
        Self::from_header(record.name().as_ref())
    }

    fn from_header(header: &'key [u8]) -> Self {
        let first_token = header
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or(header);

        let key = match first_token {
            [prefix @ .., b'/', b'1' | b'2'] => prefix,
            _ => first_token,
        };

        Self(key)
    }

    pub(super) fn as_bytes(self) -> &'key [u8] {
        self.0
    }

    pub(super) fn as_str(self) -> Result<&'key str> {
        str::from_utf8(self.0).wrap_err("FASTQ read name is not valid UTF-8")
    }
}

pub(super) fn mate_keys_match(r1: &FastqRecord, r2: &FastqRecord) -> bool {
    FastqMateKey::from_record(r1).as_bytes() == FastqMateKey::from_record(r2).as_bytes()
}

pub(super) struct FastqReadView<'read> {
    id: &'read str,
    seq: &'read str,
    qual: &'read str,
}

impl<'read> FastqReadView<'read> {
    pub(super) fn from_record(record: &'read FastqRecord, id: &'read str) -> Result<Self> {
        Ok(Self {
            id,
            seq: str::from_utf8(record.sequence()).wrap_err("FASTQ sequence is not valid UTF-8")?,
            qual: str::from_utf8(record.quality_scores())
                .wrap_err("FASTQ quality string is not valid UTF-8")?,
        })
    }
}

impl SeqRecordView for FastqReadView<'_> {
    fn id(&self) -> &str {
        self.id
    }

    fn seq(&self) -> &str {
        self.seq
    }

    fn qual(&self) -> &str {
        self.qual
    }
}

#[cfg(test)]
mod tests {
    use noodles::fastq::{Record as FastqRecord, record::Definition};

    use super::FastqMateKey;

    #[test]
    fn mate_key_strips_slash_mate_suffix() {
        assert_eq!(
            FastqMateKey::from_header(b"read123/1").as_bytes(),
            b"read123"
        );
        assert_eq!(
            FastqMateKey::from_header(b"read123/2").as_bytes(),
            b"read123"
        );
        assert_eq!(FastqMateKey::from_header(b"read123").as_bytes(), b"read123");
    }

    #[test]
    fn mate_key_uses_first_whitespace_token() {
        assert_eq!(
            FastqMateKey::from_header(b"read123/1 instrument stuff").as_bytes(),
            b"read123"
        );
        assert_eq!(
            FastqMateKey::from_header(b"read123 comment").as_bytes(),
            b"read123"
        );
    }

    #[test]
    fn mate_key_can_be_extracted_from_fastq_record() {
        let record = FastqRecord::new(
            Definition::new("read123/1", "instrument stuff"),
            "AAAA",
            "IIII",
        );

        assert_eq!(FastqMateKey::from_record(&record).as_bytes(), b"read123");
    }
}
