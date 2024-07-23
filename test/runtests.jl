using TestItemRunner

include("general.jl")
include("filtering.jl")

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(GeoGrids; ambiguities=false)
    Aqua.test_ambiguities(GeoGrids)
end

@testitem "JET" begin
    using JET
    report_package("GeoGrids")
end

@run_package_tests verbose = true