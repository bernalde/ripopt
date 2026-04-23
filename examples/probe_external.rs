//! Diagnostic: dlopen an AMPL external-function library and list its
//! registered functions.
//!
//! Usage: cargo run --example probe_external -- <path/to/library.dylib>

use ripopt::nl::external::ExternalLibrary;
use std::path::PathBuf;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = if args.len() >= 2 {
        PathBuf::from(&args[1])
    } else {
        PathBuf::from(std::env::var_os("HOME").expect("HOME not set"))
            .join(".idaes/bin/general_helmholtz_external.dylib")
    };
    let lib = match ExternalLibrary::load(&path) {
        Ok(l) => l,
        Err(e) => {
            eprintln!("load failed: {e}");
            std::process::exit(1);
        }
    };
    let mut names: Vec<String> = lib.function_names().map(|s| s.to_string()).collect();
    names.sort();
    println!("{} functions registered in {}:", names.len(), path.display());
    for n in &names {
        let rf = lib.get(n).unwrap();
        println!("  {:<24} type={} nargs={}", n, rf.ty, rf.nargs);
    }
}
