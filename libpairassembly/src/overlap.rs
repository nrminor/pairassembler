use crate::{ReadMates, ValidatedOverlap};

mod methods {
    mod bbmerge;
    mod fastp;
    mod vsearch;
}

impl ReadMates<'_> {
    pub async fn find_overlap(&self) -> color_eyre::Result<Option<ValidatedOverlap>> {
        todo!()
    }
}
