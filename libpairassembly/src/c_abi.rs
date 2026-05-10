//! Minimal C ABI for pair assembly.
//!
//! This module intentionally exposes a smaller surface than the Rust API. C callers get an opaque
//! assembler handle, borrowed input read views, an owned merged-read output, and matching free
//! functions for Rust-allocated memory.

use std::{ffi::CString, ffi::c_char, panic::AssertUnwindSafe, ptr, slice, str};

use crate::{Assembler, PairInput, SequenceRead};

/// Opaque assembler handle owned by Rust and passed through C as a pointer.
pub enum PairasmAssembler {}

struct AssemblerHandle {
    assembler: Assembler,
}

/// Borrowed byte slice passed across the C ABI.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct PairasmBytes {
    pub ptr: *const u8,
    pub len: usize,
}

/// Borrowed read record view passed across the C ABI.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct PairasmReadView {
    pub id: PairasmBytes,
    pub sequence: PairasmBytes,
    pub quality: PairasmBytes,
}

/// Owned merged read returned across the C ABI.
#[repr(C)]
#[derive(Debug)]
pub struct PairasmMergedRead {
    pub id: *mut u8,
    pub id_len: usize,
    pub sequence: *mut u8,
    pub sequence_len: usize,
    pub quality: *mut u8,
    pub quality_len: usize,
}

impl PairasmMergedRead {
    fn empty() -> Self {
        Self {
            id: ptr::null_mut(),
            id_len: 0,
            sequence: ptr::null_mut(),
            sequence_len: 0,
            quality: ptr::null_mut(),
            quality_len: 0,
        }
    }
}

/// Status code returned by C ABI functions.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairasmStatus {
    Ok = 0,
    NoOverlap = 1,
    InvalidInput = 2,
    AssemblyError = 3,
    Panic = 255,
}

/// Allocate a default assembler handle.
///
/// The caller owns the returned handle and must release it with [`pairasm_assembler_free`].
///
/// # Safety
///
/// `out` must be either null or a valid writable pointer to storage for one assembler pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn pairasm_assembler_new(out: *mut *mut PairasmAssembler) -> PairasmStatus {
    let result = std::panic::catch_unwind(AssertUnwindSafe(|| {
        if out.is_null() {
            return PairasmStatus::InvalidInput;
        }

        let Ok(assembler) = Assembler::builder().build() else {
            return PairasmStatus::AssemblyError;
        };
        let handle = Box::new(AssemblerHandle { assembler });
        let raw = Box::into_raw(handle).cast::<PairasmAssembler>();

        unsafe {
            *out = raw;
        }

        PairasmStatus::Ok
    }));

    result.unwrap_or(PairasmStatus::Panic)
}

/// Free an assembler handle allocated by [`pairasm_assembler_new`].
///
/// # Safety
///
/// `assembler` must be null or a pointer previously returned by [`pairasm_assembler_new`] that has
/// not already been freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn pairasm_assembler_free(assembler: *mut PairasmAssembler) {
    if assembler.is_null() {
        return;
    }

    unsafe {
        drop(Box::from_raw(assembler.cast::<AssemblerHandle>()));
    }
}

/// Process one read pair with a Rust-owned assembler.
///
/// On [`PairasmStatus::Ok`], `out` owns a merged read and must be freed with
/// [`pairasm_merged_read_free`]. On [`PairasmStatus::NoOverlap`], `out` is empty. On error,
/// `error_out` receives a Rust-allocated NUL-terminated message when `error_out` is non-null.
///
/// # Safety
///
/// `assembler` must be a valid pointer returned by [`pairasm_assembler_new`]. Each non-empty input
/// byte view must point to readable memory for its declared length for the duration of this call.
/// `out` must be a valid writable pointer. `error_out` may be null; when non-null, it must be a
/// valid writable pointer to storage for one error-string pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn pairasm_process_pair(
    assembler: *mut PairasmAssembler,
    forward: PairasmReadView,
    reverse: PairasmReadView,
    out: *mut PairasmMergedRead,
    error_out: *mut *mut c_char,
) -> PairasmStatus {
    let result = std::panic::catch_unwind(AssertUnwindSafe(|| unsafe {
        if !error_out.is_null() {
            *error_out = ptr::null_mut();
        }
        if !out.is_null() {
            *out = PairasmMergedRead::empty();
        }

        if assembler.is_null() || out.is_null() {
            set_error(error_out, "assembler and output pointers must not be null");
            return PairasmStatus::InvalidInput;
        }

        let handle = &mut *assembler.cast::<AssemblerHandle>();
        let forward = match read_from_view(forward) {
            Ok(read) => read,
            Err(message) => {
                set_error(error_out, &message);
                return PairasmStatus::InvalidInput;
            },
        };
        let reverse = match read_from_view(reverse) {
            Ok(read) => read,
            Err(message) => {
                set_error(error_out, &message);
                return PairasmStatus::InvalidInput;
            },
        };

        let pair = PairInput::new(forward, reverse);
        match handle.assembler.process_pair(&pair) {
            Ok(Some(merged)) => {
                *out = PairasmMergedRead {
                    id: into_owned_bytes(merged.id().as_bytes()),
                    id_len: merged.id().len(),
                    sequence: into_owned_bytes(merged.sequence_bytes()),
                    sequence_len: merged.sequence_bytes().len(),
                    quality: into_owned_bytes(merged.quality_bytes()),
                    quality_len: merged.quality_bytes().len(),
                };
                PairasmStatus::Ok
            },
            Ok(None) => PairasmStatus::NoOverlap,
            Err(error) => {
                set_error(error_out, &error.to_string());
                PairasmStatus::AssemblyError
            },
        }
    }));

    result.unwrap_or_else(|_| unsafe {
        set_error(
            error_out,
            "libpairassembly panicked while processing read pair",
        );
        PairasmStatus::Panic
    })
}

/// Free a merged read allocated by [`pairasm_process_pair`].
///
/// # Safety
///
/// `read` must be null or a valid writable pointer to a [`PairasmMergedRead`] initialized by this
/// library. Each non-null field pointer in the struct must still be owned by the struct.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn pairasm_merged_read_free(read: *mut PairasmMergedRead) {
    if read.is_null() {
        return;
    }

    unsafe {
        let read = &mut *read;
        free_owned_bytes(read.id, read.id_len);
        free_owned_bytes(read.sequence, read.sequence_len);
        free_owned_bytes(read.quality, read.quality_len);
        *read = PairasmMergedRead::empty();
    }
}

/// Free an error string allocated by this library.
///
/// # Safety
///
/// `error` must be null or a pointer previously returned through an ABI `error_out` parameter that
/// has not already been freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn pairasm_error_free(error: *mut c_char) {
    if error.is_null() {
        return;
    }

    unsafe {
        drop(CString::from_raw(error));
    }
}

unsafe fn read_from_view<'a>(view: PairasmReadView) -> Result<SequenceRead<'a>, String> {
    let id = unsafe { bytes_to_str(view.id, "read id")? };
    let sequence = unsafe { bytes_to_str(view.sequence, "read sequence")? };
    let quality = unsafe { bytes_to_str(view.quality, "read quality")? };

    SequenceRead::try_new(id, sequence, quality).map_err(|error| error.to_string())
}

unsafe fn bytes_to_str<'a>(bytes: PairasmBytes, name: &'static str) -> Result<&'a str, String> {
    if bytes.ptr.is_null() && bytes.len != 0 {
        return Err(format!(
            "{name} pointer must not be null when length is nonzero"
        ));
    }

    let bytes = if bytes.len == 0 {
        &[]
    } else {
        unsafe { slice::from_raw_parts(bytes.ptr, bytes.len) }
    };

    str::from_utf8(bytes).map_err(|error| format!("{name} must be valid UTF-8: {error}"))
}

fn into_owned_bytes(bytes: &[u8]) -> *mut u8 {
    let mut bytes = bytes.to_vec().into_boxed_slice();
    let ptr = bytes.as_mut_ptr();
    std::mem::forget(bytes);
    ptr
}

unsafe fn free_owned_bytes(ptr: *mut u8, len: usize) {
    if ptr.is_null() {
        return;
    }

    unsafe {
        drop(Box::from_raw(ptr::slice_from_raw_parts_mut(ptr, len)));
    }
}

unsafe fn set_error(error_out: *mut *mut c_char, message: &str) {
    if error_out.is_null() {
        return;
    }

    let message = message.replace('\0', "\\0");
    let Ok(message) = CString::new(message) else {
        return;
    };

    unsafe {
        *error_out = message.into_raw();
    }
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use super::*;

    fn bytes(value: &str) -> PairasmBytes {
        PairasmBytes {
            ptr: value.as_ptr(),
            len: value.len(),
        }
    }

    fn read(id: &str, sequence: &str, quality: &str) -> PairasmReadView {
        PairasmReadView {
            id: bytes(id),
            sequence: bytes(sequence),
            quality: bytes(quality),
        }
    }

    #[test]
    fn c_abi_process_pair_returns_merged_read() {
        let mut assembler = ptr::null_mut();
        let status = unsafe { pairasm_assembler_new(&raw mut assembler) };
        assert_eq!(status, PairasmStatus::Ok);

        let forward = read(
            "read-1",
            "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        );
        let reverse = read(
            "read-1",
            "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        );
        let mut merged = PairasmMergedRead::empty();
        let mut error = ptr::null_mut();

        let status = unsafe {
            pairasm_process_pair(assembler, forward, reverse, &raw mut merged, &raw mut error)
        };

        assert_eq!(status, PairasmStatus::Ok);
        assert!(error.is_null());
        assert!(!merged.sequence.is_null());
        assert_eq!(merged.sequence_len, merged.quality_len);

        unsafe {
            pairasm_merged_read_free(&raw mut merged);
            pairasm_assembler_free(assembler);
        }
    }

    #[test]
    fn c_abi_process_pair_reports_no_overlap() {
        let mut assembler = ptr::null_mut();
        let status = unsafe { pairasm_assembler_new(&raw mut assembler) };
        assert_eq!(status, PairasmStatus::Ok);

        let forward = read(
            "read-1",
            "AAAAAAAAAAAAAAAAAAAAAAAA",
            "IIIIIIIIIIIIIIIIIIIIIIII",
        );
        let reverse = read(
            "read-1",
            "CCCCCCCCCCCCCCCCCCCCCCCC",
            "IIIIIIIIIIIIIIIIIIIIIIII",
        );
        let mut merged = PairasmMergedRead::empty();

        let status = unsafe {
            pairasm_process_pair(
                assembler,
                forward,
                reverse,
                &raw mut merged,
                ptr::null_mut(),
            )
        };

        assert_eq!(status, PairasmStatus::NoOverlap);
        assert!(merged.sequence.is_null());

        unsafe {
            pairasm_assembler_free(assembler);
        }
    }

    #[test]
    fn c_abi_process_pair_rejects_invalid_input() {
        let mut assembler = ptr::null_mut();
        let status = unsafe { pairasm_assembler_new(&raw mut assembler) };
        assert_eq!(status, PairasmStatus::Ok);

        let forward = read("read-1", "ACGT", "III");
        let reverse = read("read-1", "ACGT", "IIII");
        let mut merged = PairasmMergedRead::empty();
        let mut error = ptr::null_mut();

        let status = unsafe {
            pairasm_process_pair(assembler, forward, reverse, &raw mut merged, &raw mut error)
        };

        assert_eq!(status, PairasmStatus::InvalidInput);
        assert!(!error.is_null());

        unsafe {
            pairasm_error_free(error);
            pairasm_assembler_free(assembler);
        }
    }
}
