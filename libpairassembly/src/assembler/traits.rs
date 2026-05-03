//! Boundary traits for assembler input interop.

use crate::read::SequenceRead;

/// Boundary trait for pair records accepted by the assembler API.
///
/// Implement this for external record types to use `Assembler` directly.
pub trait SeqRecordView {
    fn id(&self) -> &str;
    fn seq(&self) -> &str;
    fn qual(&self) -> &str;
}

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
