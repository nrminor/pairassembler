//! Internal typestate markers for assembler DAG transitions.

#[derive(Debug, Clone, Copy)]
pub struct NoOverlap;
#[derive(Debug, Clone, Copy)]
pub struct HasOverlap;

#[derive(Debug, Clone, Copy)]
pub struct Unvalidated;
#[derive(Debug, Clone, Copy)]
pub struct Validated;

#[derive(Debug, Clone, Copy)]
pub struct Unmerged;
#[derive(Debug, Clone, Copy)]
pub struct Merged;

#[derive(Debug, Clone, Copy)]
pub struct Uncorrected;
#[derive(Debug, Clone, Copy)]
pub struct Corrected;
