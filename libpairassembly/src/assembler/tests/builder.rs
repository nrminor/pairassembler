use crate::assembler::Assembler;

#[test]
fn test_builder_with_defaults() {
    let asm = Assembler::builder()
        .build()
        .expect("default assembler builder should produce a valid configuration");
    let _config = asm.config();
}
