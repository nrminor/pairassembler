// ROADMAP:
// For I/O to be useful in this library, I'll need to support two or three FASTQ I/O crates, each
// of which must implement some means of extracting the read ID. To make this work, I'll need to
// do the following:
//
// 1. Use noodles by default, making that the first feature, and making noodles installed with the crate by default.
// 2. Make a trait, e.g., `PairAssemble` or `MateOverlap` or something like that, that exposes an `.id()` method (and maybe some other stuff? Maybe an extension trait?) for custom `Record` types or `NewType` pattern wrappers of pre-existing libraries' record types.
// 3. Write a derive macro for the above trait, which can only be derived on structs that have an `id()` method, or something like that.
// 4. Write a bunch of `From<>` and `AsRef<>` impls to make `libpairassembly` useable with noodles and whatever else I decide to support.
// 5. Decide whether I should leave the actual reader and writer stuff to the importing crate.
// 6. Pairing up read mates could be as simple as a blanket implementation of `PartialEq` on all types that implement `MateOverlap`. This could also be a derive macro that is usable when a type first implements `MateOverlap`.
//
