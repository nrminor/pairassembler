use crate::assembler::{Assembler, ExecutionPolicy};

#[test]
fn test_builder_with_defaults() {
    let asm = Assembler::builder()
        .build()
        .expect("default assembler builder should produce a valid configuration");
    assert!(matches!(asm.config().execution, ExecutionPolicy::Record));
}
