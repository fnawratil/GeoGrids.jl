@testitem "Plots Plotly Base" tags = [:general] begin
    using PlotlyBase
    @test plot_unitarysphere(GeoGrids._icogrid(100; coord=:cart)) isa Plot
    @test plot_geo_points(icogrid(sepAng=5°)) isa Plot
    @test plot_geo_points(icogrid(sepAng=deg2rad(4) * rad); camera=:threedim) isa Plot
    @test plot_geo_points(rectgrid(5)) isa Plot
    @test plot_geo_points(rectgrid(5°); camera=:threedim) isa Plot
end