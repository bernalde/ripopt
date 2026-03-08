fn main() {
    // Link the MUMPS C wrapper and MUMPS libraries for benchmarks
    println!("cargo:rustc-link-search=native=/tmp");
    println!("cargo:rustc-link-lib=static=mumps_wrapper");

    println!("cargo:rustc-link-search=native=/opt/homebrew/Cellar/ipopt/3.14.19/lib");
    println!("cargo:rustc-link-lib=dylib=dmumps");
    println!("cargo:rustc-link-lib=dylib=mumps_common");
    println!("cargo:rustc-link-lib=dylib=mpiseq");
    println!("cargo:rustc-link-lib=dylib=pord");

    println!("cargo:rustc-link-search=native=/opt/homebrew/opt/gcc/lib/gcc/current");
    println!("cargo:rustc-link-lib=dylib=gfortran");
}
