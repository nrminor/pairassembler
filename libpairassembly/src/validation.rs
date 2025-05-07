#![allow(clippy::pedantic, clippy::perf)]
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]

use crate::{Read, ReadMates};

struct EntropyCalc {
    k: usize,
    min_entropy: EntropyStrictness,
}

#[derive(Debug, Clone, Default)]
enum EntropyStrictness {
    Loose,
    #[default]
    Normal,
    Strict,
    Other(u8),
}

impl EntropyStrictness {
    const LOOSE: u8 = 30;
    const NORMAL: u8 = 39;
    const STRICT: u8 = 44;

    fn get(&self) -> u8 {
        match self {
            EntropyStrictness::Loose => Self::LOOSE,
            EntropyStrictness::Normal => Self::NORMAL,
            EntropyStrictness::Strict => Self::STRICT,
            EntropyStrictness::Other(val) => *val,
        }
    }
}

impl Default for EntropyCalc {
    fn default() -> Self {
        let k = 3;
        let min_entropy = EntropyStrictness::default();
        Self { k, min_entropy }
    }
}

impl EntropyCalc {
    fn new() -> Self {
        Self::default()
    }

    fn with_k(self, k: usize) -> Self {
        Self { k, ..self }
    }

    fn with_min_entropy(self, min_entropy: u8) -> Self {
        let min_entropy = match min_entropy {
            30 => EntropyStrictness::Loose,
            39 => EntropyStrictness::Normal,
            44 => EntropyStrictness::Strict,
            _ if min_entropy > 0 => EntropyStrictness::Other(min_entropy),
            _ => EntropyStrictness::Normal,
        };
        Self {
            min_entropy,
            ..self
        }
    }

    fn compute_min_entropy(&self, mates: &ReadMates) -> color_eyre::Result<Option<usize>> {
        todo!()
    }
}

#[derive(Debug)]
struct MateOverlap<'a> {
    r1_offset: usize,
    r1_seq_view: &'a [u8],
    r1_qual_view: &'a [u8],
    r2_seq_view: &'a [u8],
    r2_qual_view: &'a [u8],
}

impl MateOverlap<'_> {
    fn compute_error_rate(&self) {
        todo!()
    }

    fn compute_mismatch_count(&self) {
        todo!()
    }
}

#[derive(Debug)]
pub struct ValidatedOverlap<'read> {
    mates: ReadMates<'read>,
    overlap: MateOverlap<'read>,
}

impl<'read> ValidatedOverlap<'read> {
    pub fn new(mates: ReadMates<'read>, overlap: MateOverlap<'read>) -> Self {
        todo!()
    }

    pub fn merge(self) -> color_eyre::Result<Read<'read>> {
        todo!()
    }
}
