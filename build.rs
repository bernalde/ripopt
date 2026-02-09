fn main() {
    // Only link ipopt when the ipopt-native feature is enabled
    if std::env::var("CARGO_FEATURE_IPOPT_NATIVE").is_ok() {
        let output = std::process::Command::new("pkg-config")
            .args(["--libs-only-L", "ipopt"])
            .output()
            .expect("pkg-config not found; install ipopt via homebrew");
        let lib_path = String::from_utf8(output.stdout).unwrap();
        let lib_path = lib_path.trim().trim_start_matches("-L");
        println!("cargo:rustc-link-search=native={}", lib_path);
        println!("cargo:rustc-link-lib=dylib=ipopt");

        let output = std::process::Command::new("pkg-config")
            .args(["--cflags-only-I", "ipopt"])
            .output()
            .unwrap();
        let inc_path = String::from_utf8(output.stdout).unwrap();
        let inc_path = inc_path.trim().trim_start_matches("-I");
        println!("cargo:rustc-env=IPOPT_INCLUDE={}", inc_path);
    }
}
