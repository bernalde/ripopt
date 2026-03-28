using Test

# Point RIPOPT_LIBRARY_PATH to the built library
const RIPOPT_LIB_DIR = get(
    ENV, "RIPOPT_LIBRARY_PATH",
    joinpath(@__DIR__, "..", "..", "target", "release"),
)
ENV["RIPOPT_LIBRARY_PATH"] = RIPOPT_LIB_DIR

using Ripopt

@testset "Ripopt.jl" begin
    include("test_C_wrapper.jl")
    include("test_MOI_wrapper.jl")
end
