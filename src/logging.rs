//! Logging infrastructure for ripopt.
//!
//! All solver output goes through `log_to_current`. By default this writes to
//! stderr via `eprintln!`. Frontends (GAMS link, custom C code) can redirect
//! output to their own logging system by installing a callback with
//! `ripopt_set_log_callback` (see `c_api.rs`).
//!
//! The callback is stored in a thread-local so that concurrent solves on
//! different threads each get their own channel.

use std::cell::Cell;
use std::ffi::CString;
use std::os::raw::{c_char, c_void};

/// C-callable log callback: receives a NUL-terminated message and a user_data pointer.
pub type LogCb = unsafe extern "C" fn(msg: *const c_char, user_data: *mut c_void);

thread_local! {
    /// Active log callback for the current thread, if any.
    static LOG_CALLBACK: Cell<Option<(LogCb, *mut c_void)>> = Cell::new(None);
}

/// Install a log callback for the current thread.
/// Pass `None` to revert to `eprintln!`.
pub(crate) fn set_log_callback(cb: Option<(LogCb, *mut c_void)>) {
    LOG_CALLBACK.with(|cell| cell.set(cb));
}

/// Emit a log message, routing it through the active callback or `eprintln!`.
pub(crate) fn log_to_current(msg: &str) {
    LOG_CALLBACK.with(|cell| {
        if let Some((cb, user_data)) = cell.get() {
            if let Ok(cstr) = CString::new(msg) {
                unsafe { cb(cstr.as_ptr(), user_data) };
            }
        } else {
            eprintln!("{}", msg);
        }
    });
}

/// Convenience macro — same syntax as `eprintln!` but routes through the
/// active log callback (or stderr when no callback is installed).
macro_rules! rip_log {
    ($($arg:tt)*) => {
        $crate::logging::log_to_current(&format!($($arg)*))
    };
}
pub(crate) use rip_log;
