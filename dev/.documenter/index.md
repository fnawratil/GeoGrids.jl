


# GeoGrids

Documentation for [GeoGrids](https://gitlab.esa.int/tec-esc-tools/GeoGrids.jl).
- [`GeoGrids.EO`](#GeoGrids.EO)
- [`GeoGrids.GeoRegion`](#GeoGrids.GeoRegion)
- [`GeoGrids.GeoRegionOffset`](#GeoGrids.GeoRegionOffset)
- [`GeoGrids.GlobalRegion`](#GeoGrids.GlobalRegion)
- [`GeoGrids.H3`](#GeoGrids.H3)
- [`GeoGrids.HEX`](#GeoGrids.HEX)
- [`GeoGrids.ICO`](#GeoGrids.ICO)
- [`GeoGrids.LatBeltRegion`](#GeoGrids.LatBeltRegion)
- [`GeoGrids.MultiBorder`](#GeoGrids.MultiBorder)
- [`GeoGrids.PolyBorder`](#GeoGrids.PolyBorder)
- [`GeoGrids.PolyRegion`](#GeoGrids.PolyRegion)
- [`GeoGrids.PolyRegionOffset`](#GeoGrids.PolyRegionOffset)
- [`GeoGrids._adapted_icogrid`](#GeoGrids._adapted_icogrid-Tuple{Number})
- [`GeoGrids._add_angular_offset`](#GeoGrids._add_angular_offset-Tuple{Any,%20Any})
- [`GeoGrids._fibonaccisphere_classic_partial`](#GeoGrids._fibonaccisphere_classic_partial-Tuple{Any})
- [`GeoGrids._find_min_separation_angle`](#GeoGrids._find_min_separation_angle-Tuple{Any})
- [`GeoGrids._gen_regular_lattice`](#GeoGrids._gen_regular_lattice-Union{Tuple{T},%20Tuple{T,%20Any,%20Any}}%20where%20T)
- [`GeoGrids._get_theta_phi`](#GeoGrids._get_theta_phi-Tuple{Number,%20Number})
- [`GeoGrids._hex_tesselation_centroids`](#GeoGrids._hex_tesselation_centroids-Tuple{Point{Meshes.🌐,%20<:CoordRefSystems.GeodeticLatLon{WGS84Latest}},%20Number})
- [`GeoGrids._icogrid`](#GeoGrids._icogrid-Tuple{Int64})
- [`GeoGrids._offset_ring`](#GeoGrids._offset_ring-Union{Tuple{T},%20Tuple{Meshes.Ring{Meshes.🌐,%20CoordRefSystems.GeodeticLatLon{WGS84Latest,%20Unitful.Quantity{T,%20NoDims,%20Unitful.FreeUnits{(°,),%20NoDims,%20nothing}}},%20CircularArrays.CircularArray{Point{Meshes.🌐,%20CoordRefSystems.GeodeticLatLon{WGS84Latest,%20Unitful.Quantity{T,%20NoDims,%20Unitful.FreeUnits{(°,),%20NoDims,%20nothing}}}},%201,%20Array{Point{Meshes.🌐,%20CoordRefSystems.GeodeticLatLon{WGS84Latest,%20Unitful.Quantity{T,%20NoDims,%20Unitful.FreeUnits{(°,),%20NoDims,%20nothing}}}},%201}}},%20Any}}%20where%20T)
- [`GeoGrids._points_required_for_separation_angle`](#GeoGrids._points_required_for_separation_angle-Tuple{Union{Real,%20Unitful.Quantity{<:Real,%20<:Any,%20<:Union{Unitful.FreeUnits{(°,),%20NoDims,%20nothing},%20Unitful.FreeUnits{(rad,),%20NoDims,%20nothing}}}}})
- [`GeoGrids._tesselate`](#GeoGrids._tesselate)
- [`GeoGrids._wrap_latlon`](#GeoGrids._wrap_latlon-Tuple{Number,%20Number})
- [`GeoGrids.cartesian_geometry`](#GeoGrids.cartesian_geometry-Tuple{PolyArea{Meshes.🌐,%20var"#s12",%20R,%20V}%20where%20{var"#s12"<:(CoordRefSystems.GeodeticLatLon{WGS84Latest,%20Unitful.Quantity{T,%20NoDims,%20Unitful.FreeUnits{(°,),%20NoDims,%20nothing}}}%20where%20T),%20R<:(Meshes.Ring{Meshes.🌐,%20var"#s12",%20V}%20where%20V<:(CircularArrays.CircularArray{Point{Meshes.🌐,%20var"#s12"},%201,%20A}%20where%20A<:AbstractArray{Point{Meshes.🌐,%20var"#s12"},%201})),%20V<:AbstractVector{R}}})
- [`GeoGrids.fibonaccisphere_alternative1`](#GeoGrids.fibonaccisphere_alternative1-Tuple{Int64})
- [`GeoGrids.fibonaccisphere_optimization1`](#GeoGrids.fibonaccisphere_optimization1-Tuple{Int64})
- [`GeoGrids.filter_points`](#GeoGrids.filter_points-Tuple{AbstractVector{<:Union{CoordRefSystems.GeodeticLatLon,%20Point{Meshes.🌐,%20<:CoordRefSystems.GeodeticLatLon{WGS84Latest}}}},%20Union{LatBeltRegion,%20PolyRegion,%20PolyRegionOffset}})
- [`GeoGrids.gen_circle_pattern`](#GeoGrids.gen_circle_pattern-Tuple{AbstractVector{<:Point{Meshes.🌐,%20<:CoordRefSystems.GeodeticLatLon{WGS84Latest}}},%20Number})
- [`GeoGrids.gen_hex_lattice`](#GeoGrids.gen_hex_lattice)
- [`GeoGrids.gen_hex_pattern`](#GeoGrids.gen_hex_pattern-Tuple{AbstractVector{<:Point{Meshes.🌐,%20<:CoordRefSystems.GeodeticLatLon{WGS84Latest}}},%20AbstractVector{<:Number},%20Meshes.SimpleMesh})
- [`GeoGrids.generate_tesselation`](#GeoGrids.generate_tesselation-Tuple{GlobalRegion,%20Number,%20ICO})
- [`GeoGrids.generate_tesselation`](#GeoGrids.generate_tesselation-Tuple{Union{GeoRegion,%20GeoRegionOffset,%20PolyRegion,%20PolyRegionOffset},%20Number,%20HEX})
- [`GeoGrids.icogrid`](#GeoGrids.icogrid-Tuple{})
- [`GeoGrids.latlon_geometry`](#GeoGrids.latlon_geometry-Tuple{PolyArea{Meshes.𝔼{2},%20var"#s12",%20R,%20V}%20where%20{var"#s12"<:(Cartesian{WGS84Latest,%202,%20Unitful.Quantity{T,%20𝐋,%20Unitful.FreeUnits{(m,),%20𝐋,%20nothing}}}%20where%20T),%20R<:(Meshes.Ring{Meshes.𝔼{2},%20var"#s12",%20V}%20where%20V<:(CircularArrays.CircularArray{Point{Meshes.𝔼{2},%20var"#s12"},%201,%20A}%20where%20A<:AbstractArray{Point{Meshes.𝔼{2},%20var"#s12"},%201})),%20V<:AbstractVector{R}}})
- [`GeoGrids.offset_region`](#GeoGrids.offset_region-Union{Tuple{P},%20Tuple{D},%20Tuple{GeoRegion{D,%20P},%20Any}}%20where%20{D,%20P})
- [`GeoGrids.rectgrid`](#GeoGrids.rectgrid-Tuple{Union{Real,%20Unitful.Quantity{<:Real,%20<:Any,%20<:Union{Unitful.FreeUnits{(°,),%20NoDims,%20nothing},%20Unitful.FreeUnits{(rad,),%20NoDims,%20nothing}}}}})
- [`GeoGrids.vecgrid`](#GeoGrids.vecgrid-Tuple{Union{Real,%20Unitful.Quantity{<:Real,%20<:Any,%20<:Union{Unitful.FreeUnits{(°,),%20NoDims,%20nothing},%20Unitful.FreeUnits{(rad,),%20NoDims,%20nothing}}}}})


## Public
<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.EO' href='#GeoGrids.EO'><span class="jlbinding">GeoGrids.EO</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
EO
```


Struct used to create function methods that return more than one output. Used within multiple methods of the GeoGrids API, usually given as last optional argument.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L208-L213)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.GeoRegion' href='#GeoGrids.GeoRegion'><span class="jlbinding">GeoGrids.GeoRegion</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GeoRegion{D,P} <: AbstractRegion
```


Type representing a geographical region based on CountriesBorders.

Fields:
- `name::String`: Name of the region
- `continent::String`: Continent of the region
- `subregion::String`: Subregion within the continent
- `admin::String`: Administrative area
- `domain::D`: Domain of the region
- `convexhull::PolyBorder{P}`: Convex hull of the region

Where `D` is the domain type and `P` is the precision type for coordinates.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L64-L78)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.GeoRegionOffset' href='#GeoGrids.GeoRegionOffset'><span class="jlbinding">GeoGrids.GeoRegionOffset</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GeoRegionOffset{D,P} <: AbstractRegion
```


Type representing an enlarged geographical region based on a GeoRegion.

Fields:
- `original::GeoRegion{D,P}`: The original GeoRegion
- `name::String`: Name of the enlarged region
- `domain::MultiBorder{P}`: Domain of the enlarged region
- `convexhull::PolyBorder{P}`: Convex hull of the enlarged region

**Constructors**

```
GeoRegionOffset(; name="offset_georegion", continent="", subregion="", admin="", delta::Number, resolution=110, refRadius=constants.Re_mean, magnitude=3, precision=7)
GeoRegionOffset(gr::GeoRegion, delta::Number; name="offset_georegion", refRadius=constants.Re_mean, magnitude=3, precision=7)
```


Create an enlarged GeoRegion either from scratch or from an existing GeoRegion.

**Arguments**
- `delta`: Distance to enlarge the region by, in meters
- `gr::GeoRegion`: The original GeoRegion to enlarge (for the second constructor)

**Keyword Arguments**
- `name::String="enlarged_region"`: Name of the enlarged region
- `continent::String=""`: Continent of the region (only for the first constructor)
- `subregion::String=""`: Subregion within the continent (only for the first constructor)
- `admin::String=""`: Administrative area (only for the first constructor)
- `resolution::Int=110`: Resolution of the geographical data (only for the first constructor)
- `refRadius::Float64=constants.Re_mean`: Reference radius of the Earth
- `magnitude::Int=3`: Magnitude for polygon offsetting
- `precision::Int=7`: Precision for polygon offsetting

**Returns**
- `GeoRegionOffset`: The enlarged geographical region


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/offset_types.jl#L1-L35)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.GlobalRegion' href='#GeoGrids.GlobalRegion'><span class="jlbinding">GeoGrids.GlobalRegion</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GlobalRegion <: AbstractRegion
```


Type representing a global region.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L57-L61)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.H3' href='#GeoGrids.H3'><span class="jlbinding">GeoGrids.H3</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
H3 <: AbstractTiling
```


Struct representing an H3 tiling.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L201-L205)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.HEX' href='#GeoGrids.HEX'><span class="jlbinding">GeoGrids.HEX</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
HEX <: AbstractTiling
```


Struct representing a hexagonal tiling.

Fields:
- `direction::Symbol`: Default direction of hexagons in the tiling (:pointy or :flat)
- `pattern::Symbol`: Default pattern shape to be used with this type of tiling (:circ or :hex)


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L179-L187)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.ICO' href='#GeoGrids.ICO'><span class="jlbinding">GeoGrids.ICO</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ICO <: AbstractTiling
```


Struct representing an icosahedral tiling.

Fields:
- `correction::Number`: Default correction factor for the icosahedral cell grid partial overlap
- `pattern::Symbol`: Default pattern shape to be used with this type of tiling (:circ or :hex)


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L158-L166)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.LatBeltRegion' href='#GeoGrids.LatBeltRegion'><span class="jlbinding">GeoGrids.LatBeltRegion</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
LatBeltRegion <: AbstractRegion
```


Type representing a latitude belt region.

Fields:
- `name::String`: Name of the region
- `lim::Tuple{ValidAngle,ValidAngle}`: Latitude limits of the belt in degrees


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L119-L127)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.MultiBorder' href='#GeoGrids.MultiBorder'><span class="jlbinding">GeoGrids.MultiBorder</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MultiBorder{P} <: Geometry{🌐,LATLON}
```


Struct representing a Multi in both LatLon and Cartesian coordinates.

Fields:
- `latlon::MULTI_LATLON{P}`: The borders in LatLon CRS
- `cart::MULTI_CART{P}`: The borders in Cartesian2D CRS

Where `P` is the precision type (e.g., Float32, Float64) for the coordinates.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L34-L44)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.PolyBorder' href='#GeoGrids.PolyBorder'><span class="jlbinding">GeoGrids.PolyBorder</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PolyBorder{P} <: Geometry{🌐,LATLON}
```


Struct representing a PolyArea in both LatLon and Cartesian coordinates.

Fields:
- `latlon::POLY_LATLON{P}`: The borders in LatLon CRS
- `cart::POLY_CART{P}`: The borders in Cartesian2D CRS

Where `P` is the precision type (e.g., Float32, Float64) for the coordinates.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L13-L23)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.PolyRegion' href='#GeoGrids.PolyRegion'><span class="jlbinding">GeoGrids.PolyRegion</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PolyRegion{P} <: AbstractRegion
```


Type representing a polygonal region based on PolyArea.

Fields:
- `name::String`: Name of the region
- `domain::PolyBorder{P}`: Domain of the region as a PolyBorder

Where `P` is the precision type for coordinates.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/basic_types.jl#L100-L110)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.PolyRegionOffset' href='#GeoGrids.PolyRegionOffset'><span class="jlbinding">GeoGrids.PolyRegionOffset</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PolyRegionOffset{P} <: AbstractRegion
```


Struct representing an enlarged polygonal region.

Fields:
- `original::PolyRegion{P}`: The original PolyRegion
- `name::String`: Name of the enlarged region
- `domain::MultiBorder{P}`: Domain of the enlarged region as a MultiBorder

**Constructors**

```
PolyRegionOffset(delta::Number; kwargs...)
PolyRegionOffset(pr::PolyRegion, delta::Number; kwargs...)
```


Create an enlarged PolyRegion either from scratch or from an existing PolyRegion.

**Arguments**
- `deltaDist`: Distance to enlarge the region by, in meters
- `pr::PolyRegion`: The original PolyRegion to enlarge (for the second constructor)

**Keyword Arguments**
- `name::String="offset_polyregion"`: Name of the enlarged region
- `domain`: Domain of the region (only for the first constructor)
- `refRadius::Float64=constants.Re_mean`: Reference radius of the Earth
- `magnitude::Int=3`: Magnitude for polygon offsetting
- `precision::Int=7`: Precision for polygon offsetting

**Returns**
- `PolyRegionOffset`: The enlarged polygonal region


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/offset_types.jl#L62-L92)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._tesselate' href='#GeoGrids._tesselate'><span class="jlbinding">GeoGrids._tesselate</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
_tesselate(points::AbstractVector{<:Point{🌐,<:LatLon{WGS84Latest}}}, method::TesselationMethod=VoronoiTesselation()) -> TesselationResult
_tesselate(point::Point{🌐,<:LatLon{WGS84Latest}}; kwargs...) -> TesselationResult
```


The function `_tesselate` uses tesselate from Meshes.jl is used to create a tasselation starting from a vector of geographical points (latitude and longitude) and tesselates them according to the specified tesselation method (default is `VoronoiTesselation()`). In this function, latitude is treated as the y-coordinate and longitude as the x-coordinate. The single point version of the function (`_tesselate(point::Point{🌐,<:LatLon{WGS84Latest}}; kwargs...)`) is a convenience method that allows you to tesselate a single point by internally converting it to a vector containing just that point.

**Arguments**
- `points::AbstractVector{<:Point{🌐,<:LatLon{WGS84Latest}}}`: A vector of points defined in the WGS84 coordinate system. Each point represents a geographical location with latitude and longitude.
- `method::TesselationMethod=VoronoiTesselation()`: The method used for tesselation. The default is `VoronoiTesselation()`, but other methods can be specified.

**Keyword Arguments (for the second method signature)**
- `kwargs...`: Additional keyword arguments that will be passed to the first method.

**Returns**
- `TesselationResult`: The result of the tesselation, which could be a set of polygons or other geometrical structures, depending on the tesselation method used.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L329-L351)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.cartesian_geometry-Tuple{PolyArea{Meshes.🌐, var"#s12", R, V} where {var"#s12"<:(CoordRefSystems.GeodeticLatLon{WGS84Latest, Unitful.Quantity{T, NoDims, Unitful.FreeUnits{(°,), NoDims, nothing}}} where T), R<:(Meshes.Ring{Meshes.🌐, var"#s12", V} where V<:(CircularArrays.CircularArray{Point{Meshes.🌐, var"#s12"}, 1, A} where A<:AbstractArray{Point{Meshes.🌐, var"#s12"}, 1})), V<:AbstractVector{R}}}' href='#GeoGrids.cartesian_geometry-Tuple{PolyArea{Meshes.🌐, var"#s12", R, V} where {var"#s12"<:(CoordRefSystems.GeodeticLatLon{WGS84Latest, Unitful.Quantity{T, NoDims, Unitful.FreeUnits{(°,), NoDims, nothing}}} where T), R<:(Meshes.Ring{Meshes.🌐, var"#s12", V} where V<:(CircularArrays.CircularArray{Point{Meshes.🌐, var"#s12"}, 1, A} where A<:AbstractArray{Point{Meshes.🌐, var"#s12"}, 1})), V<:AbstractVector{R}}}'><span class="jlbinding">GeoGrids.cartesian_geometry</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cartesian_geometry(poly::PolyArea{🌐,<:LATLON})
cartesian_geometry(multi::Multi{🌐,<:LATLON})
```


Convert geometries from LatLon to Cartesian coordinate systems.

**Arguments**
- `poly::PolyArea{🌐,<:LATLON}`: A polygon in LatLon coordinates.
- `multi::Multi{🌐,<:LATLON}`: A multi-geometry in LatLon coordinates.

**Returns**
- `PolyArea` or `Multi`: The converted geometry in Cartesian coordinate system.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/helper_func.jl#L95-L107)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.filter_points-Tuple{AbstractVector{<:Union{CoordRefSystems.GeodeticLatLon, Point{Meshes.🌐, <:CoordRefSystems.GeodeticLatLon{WGS84Latest}}}}, Union{LatBeltRegion, PolyRegion, PolyRegionOffset}}' href='#GeoGrids.filter_points-Tuple{AbstractVector{<:Union{CoordRefSystems.GeodeticLatLon, Point{Meshes.🌐, <:CoordRefSystems.GeodeticLatLon{WGS84Latest}}}}, Union{LatBeltRegion, PolyRegion, PolyRegionOffset}}'><span class="jlbinding">GeoGrids.filter_points</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
filter_points(points::AbstractVector{<:LatLon}, domain::Union{GeoRegion, PolyRegion, LatBeltRegion, GeoRegionOffset, PolyRegionOffset}) -> Vector{Input Type}
```


Filters a list of points based on whether they fall within a specified geographical domain.

**Arguments**
- `points`: An array of points. The points are `LatLon`.
- `domain`: A geographical domain which can be of type `GeoRegion` or `PolyRegion`, `LatBeltRegion`, `GeoRegionOffset`, or `PolyRegionOffset` in alternative a `Meshes.Domain` of type `GeometrySet` or `PolyArea`.
- `::EO`: An `EO` object for additional output containing the indices of the filtered points (wrt the input).

**Returns**
- A vector of points that fall within the specified domain, subsection of the input vector. The output is of the same type as the input.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/filtering_func.jl#L1-L14)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.gen_circle_pattern-Tuple{AbstractVector{<:Point{Meshes.🌐, <:CoordRefSystems.GeodeticLatLon{WGS84Latest}}}, Number}' href='#GeoGrids.gen_circle_pattern-Tuple{AbstractVector{<:Point{Meshes.🌐, <:CoordRefSystems.GeodeticLatLon{WGS84Latest}}}, Number}'><span class="jlbinding">GeoGrids.gen_circle_pattern</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
gen_circle_pattern(centers::AbstractVector{Point{🌐,<:LatLon{WGS84Latest}}}, radius::Number; refRadius::Number=constants.Re_mean, n::Int=20) -> Vector{Vector{Point{🌐,<:LatLon{WGS84Latest}}}
gen_circle_pattern(c::Point{🌐,<:LatLon{WGS84Latest}}, radius::Number; kwargs...) -> Vector{Point{🌐,<:LatLon{WGS84Latest}}
gen_circle_pattern(centers::AbstractVector{<:LatLon}, radius::Number; kwargs...) -> Vector{Vector{Point{🌐,<:LatLon{WGS84Latest}}}
gen_circle_pattern(c::LatLon, radius::Number; kwargs...) -> Vector{Point{🌐,<:LatLon{WGS84Latest}}
```


The `gen_circle_pattern` function generates circles of geographical points centered at each point in the `centers` vector. The points are generated on the Earth&#39;s surface using a spherical approximation, where latitude and longitude are converted to spherical coordinates (theta-phi), and then an angular offset is applied to generate the circle. The single point versions of the function (`gen_circle_pattern(c::Point{🌐,<:LatLon{WGS84Latest}}, radius::Number; kwargs...)` and `gen_circle_pattern(c::LatLon, radius::Number; kwargs...)`) are convenience methods that allow you to generate a circle pattern around a single center point. Tis function is used to create a plottable patter od the circles around a center point.

**Arguments**
- `centers::AbstractVector{Point{🌐,<:LatLon{WGS84Latest}}}`: A vector of geographical points in the WGS84 coordinate system. Each point represents the center of a circle.
- `radius::Number`: The radius of the circles to generate, in the same units as the reference radius (`refRadius`).
- `refRadius::Number=constants.Re_mean`: The reference radius for the spherical approximation, defaulting to the mean Earth radius (`Re_mean`).
- `n::Int=20`: The number of points to generate along each circle&#39;s circumference.

**Keyword Arguments**
- `kwargs...`: Additional keyword arguments passed to other variations of the function.

**Returns**
- `Vector{Vector{Point{🌐,<:LatLon{WGS84Latest}}}}`: A vector where each element is a vector of `LatLon` points representing a circle.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L368-L389)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.gen_hex_lattice' href='#GeoGrids.gen_hex_lattice'><span class="jlbinding">GeoGrids.gen_hex_lattice</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
gen_hex_lattice(spacing, direction = :pointy; kwargs...)
```


Generate a hexagonal lattice of points with equal distance `spacing` between neighboring points.

The generated hexagonal lattice will have distance between points on the same row/column that will depend on the second argument `direction`:

If `direction = :pointy`, neighboring points on the same row (points which have the same `y` coordinate) will be at a distance `spacing` from one another, while points on the same column (sharing the `x` coordinate) will have a distance equivalent to `√3 * spacing`.

If `direction = :flat`, the distance will be reversed, so points on the same column will have a distance equivalent to `spacing` while points on the same row will have a distance equivalent to `√3 * spacing`.

**Arguments**
- `spacing`: spacing between neighboring points
- `direction`: specifies the direction of minimum distance between neighboringpoints. Defaults to `:pointy`.

See also: [`_gen_regular_lattice`](/index#GeoGrids._gen_regular_lattice-Union{Tuple{T},%20Tuple{T,%20Any,%20Any}}%20where%20T)


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L189-L211)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.gen_hex_pattern-Tuple{AbstractVector{<:Point{Meshes.🌐, <:CoordRefSystems.GeodeticLatLon{WGS84Latest}}}, AbstractVector{<:Number}, Meshes.SimpleMesh}' href='#GeoGrids.gen_hex_pattern-Tuple{AbstractVector{<:Point{Meshes.🌐, <:CoordRefSystems.GeodeticLatLon{WGS84Latest}}}, AbstractVector{<:Number}, Meshes.SimpleMesh}'><span class="jlbinding">GeoGrids.gen_hex_pattern</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
gen_hex_pattern(filtered::AbstractVector{Point{🌐,<:LatLon{WGS84Latest}}}, idxs::AbstractVector{<:Number}, mesh::SimpleMesh) -> Vector{Vector{LatLon}}
gen_hex_pattern(p::Point{🌐,<:LatLon{WGS84Latest}}, idx::Number, mesh::SimpleMesh) -> Vector{LatLon}
gen_hex_pattern(filtered::AbstractVector{<:LatLon}, idxs::AbstractVector{<:Number}, mesh::SimpleMesh) -> Vector{Vector{LatLon}}
gen_hex_pattern(p::LatLon, idx::Number, mesh::SimpleMesh) -> Vector{LatLon}
```


The `gen_hex_pattern` function generates patterns of hexagons (or other polygons) around a set of geographical points using a provided mesh. The mesh is expected to contain polygons that are indexed by the `idxs` argument, and each polygon is converted into a set of geographical points in latitude and longitude. The function iterates over each polygon in the mesh corresponding to the indices in `idxs`, converting the vertices of the polygon into `LatLon` points that represent the corners of the hexagon (or other polygon) around the corresponding center point. The single point versions of the function (`gen_hex_pattern(p::Point{🌐,<:LatLon{WGS84Latest}}, idx::Number, mesh::SimpleMesh)` and `gen_hex_pattern(p::LatLon, idx::Number, mesh::SimpleMesh)`) are convenience methods that allow you to generate a pattern around a single center point.

**Arguments**
- `filtered::AbstractVector{Point{🌐,<:LatLon{WGS84Latest}}}`: A vector of geographical points in the WGS84 coordinate system that represent the centers of the hexagons or polygons.
- `idxs::AbstractVector{<:Number}`: A vector of indices corresponding to the polygons within the mesh.
- `mesh::SimpleMesh`: A mesh object containing the polygons (typically hexagons) used to generate the patterns.

**Keyword Arguments**
- `kwargs...`: Additional keyword arguments passed to other variations of the function.

**Returns**
- `Vector{Vector{LatLon}}`: A vector where each element is a vector of `LatLon` points representing the vertices of the polygons (typically hexagons) for each center point.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L414-L443)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.generate_tesselation-Tuple{GlobalRegion, Number, ICO}' href='#GeoGrids.generate_tesselation-Tuple{GlobalRegion, Number, ICO}'><span class="jlbinding">GeoGrids.generate_tesselation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
generate_tesselation(region::GlobalRegion, radius::Number, type::ICO; refRadius::Number=constants.Re_mean) -> AbstractVector{<:LatLon}
generate_tesselation(region::GlobalRegion, radius::Number, type::ICO, ::EO; refRadius::Number=constants.Re_mean) -> AbstractVector{<:LatLon}, AbstractVector{<:AbstractVector{<:LatLon}}
generate_tesselation(region::Union{LatBeltRegion, GeoRegion, PolyRegion, GeoRegionOffset, PolyRegionOffset}, radius::Number, type::ICO; refRadius::Number=constants.Re_mean) -> AbstractVector{<:LatLon}
generate_tesselation(region::Union{LatBeltRegion, GeoRegion, PolyRegion, GeoRegionOffset, PolyRegionOffset}, radius::Number, type::ICO, ::EO; refRadius::Number=constants.Re_mean) -> AbstractVector{<:LatLon}, AbstractVector{<:AbstractVector{<:LatLon}}
```


The `generate_tesselation` function generates a cell layout using an icosahedral grid for a given geographical region. The function adapts the grid based on the specified radius and applies a correction factor (see `_adapted_icogrid`). The radius as to be intended as the semparation angle of the point on the icosahedral grid. It then filters the grid points to include only those within the specified region.

**Arguments**
- `region::Union{LatBeltRegion, GeoRegion, PolyRegion}`: The geographical region for which the cell layout is generated. This can be a `LatBeltRegion`, `GeoRegion`, or `PolyRegion`.
- `radius::Number`: The radius used to adapt the icosahedral grid.
- `type::ICO`: An object specifying the type of icosahedral grid and its correction factor.
- `::EO`: an extra parameter enabling a `Vector{Ngon}` the contours of each cell. The mesh originating these contours is obtained using `VoronoiTesselation`.

**Returns**
- `Vector{Point{🌐,<:LatLon{WGS84Latest}}}`: A Vecotr of points (`LatLon`) representing the cell centers within the specified region.

See also: [`_adapted_icogrid()`](/index#GeoGrids._adapted_icogrid-Tuple{Number}), [`icogrid()`](/index#GeoGrids.icogrid-Tuple{}), [`filter_points()`](/index#GeoGrids.filter_points-Tuple{AbstractVector{<:Union{CoordRefSystems.GeodeticLatLon,%20Point{Meshes.🌐,%20<:CoordRefSystems.GeodeticLatLon{WGS84Latest}}}},%20Union{LatBeltRegion,%20PolyRegion,%20PolyRegionOffset}}), [`GeoRegion`](/index#GeoGrids.GeoRegion), [`LatBeltRegion`](/index#GeoGrids.LatBeltRegion), [`PolyRegion`](/index#GeoGrids.PolyRegion), [`GlobalRegion`](/index#GeoGrids.GlobalRegion), [`ICO`](/index#GeoGrids.ICO)


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L68-L93)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.generate_tesselation-Tuple{Union{GeoRegion, GeoRegionOffset, PolyRegion, PolyRegionOffset}, Number, HEX}' href='#GeoGrids.generate_tesselation-Tuple{Union{GeoRegion, GeoRegionOffset, PolyRegion, PolyRegionOffset}, Number, HEX}'><span class="jlbinding">GeoGrids.generate_tesselation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
generate_tesselation(region::Union{GeoRegion, PolyRegion, GeoRegionOffset, PolyRegionOffset}, radius::Number, type::HEX; refRadius::Number=constants.Re_mean, kwargs_lattice...) -> AbstractVector{<:Point{🌐,<:LatLon{WGS84Latest}}}
generate_tesselation(region::Union{GeoRegion, PolyRegion, GeoRegionOffset, PolyRegionOffset}, radius::Number, type::HEX, ::EO; refRadius::Number=constants.Re_mean, kwargs_lattice...) -> AbstractVector{<:Point{🌐,<:LatLon{WGS84Latest}}}, AbstractVector{<:AbstractVector{<:Point{🌐,<:LatLon{WGS84Latest}}}}
```


The `generate_tesselation` function generates a hexagonal cell layout for a given geographical region. It calculates the cell grid layout centered around the centroid of the main area of the region and returns the points within the specified region. A method for using hexagonal tesselation for GlobalRegion and LatBeltRegion is not provided because of the large size of surface to be covered, the local hexagonal tasellation would be inaccurate.

**Arguments**
- `region::Union{GeoRegion, PolyRegion}`: The geographical region for which the

cell layout is generated. Larger regions like global and LatBeltRegions are not supported because of the problem of regular tassellation of the sphere.
- `radius::Number`: The radius of each hexagonal cell. Has to be intended as the circumscribed circumference.
- `type::HEX`: A parameter indicating the type of lattice (only HEX is supported).
- `refRadius::Number`: The radius of the Earth in meters (default is `constants.Re_mean`).
- `::EO`: an extra parameter enabling a `Vector{Ngon}` the contours of each cell. The mesh originating these contours is obtained using `VoronoiTesselation`.
- `kwargs_lattice...`: Additional keyword arguments passed to the `gen_hex_lattice` function.

**Returns**
- `Vector{Point{🌐,<:LatLon{WGS84Latest}}}`: A Vecotr of points (`LatLon`) representing the cell centers within the specified region.

See also: [`gen_hex_lattice`](/index#GeoGrids.gen_hex_lattice), [`_generate_tesselation`](@ref), [`_hex_tesselation_centroids`](/index#GeoGrids._hex_tesselation_centroids-Tuple{Point{Meshes.🌐,%20<:CoordRefSystems.GeodeticLatLon{WGS84Latest}},%20Number}), [`_tesselate`](/index#GeoGrids._tesselate), [`HEX`](/index#GeoGrids.HEX), [`GeoRegion`](/index#GeoGrids.GeoRegion), [`PolyRegion`](/index#GeoGrids.PolyRegion)


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L1-L29)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.icogrid-Tuple{}' href='#GeoGrids.icogrid-Tuple{}'><span class="jlbinding">GeoGrids.icogrid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
icogrid(;N=nothing, sepAng=nothing, unit=:rad, type=:lla) -> Vector{Point{🌐,<:LatLon{WGS84Latest}}}
```


This function returns a vector `Nx2` of LAT, LON values for a `N` points grid built with the Fibonacci Spiral method.

**Arguments:**
- `N`: The number of points to generate.
- `sepAng`: The separation angle for the grid of points to be generated [rad].
- `unit`: `:rad` or `:deg`
- `type`: `:lla` or `:point`. Output type either `LLA` or `Point2`

**Output:**
- `out`: an Vector{Point{🌐,&lt;:LatLon{WGS84Latest}}} of points in the icosahedral grid.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/ico_func.jl#L1-L15)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.latlon_geometry-Tuple{PolyArea{Meshes.𝔼{2}, var"#s12", R, V} where {var"#s12"<:(Cartesian{WGS84Latest, 2, Unitful.Quantity{T, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}} where T), R<:(Meshes.Ring{Meshes.𝔼{2}, var"#s12", V} where V<:(CircularArrays.CircularArray{Point{Meshes.𝔼{2}, var"#s12"}, 1, A} where A<:AbstractArray{Point{Meshes.𝔼{2}, var"#s12"}, 1})), V<:AbstractVector{R}}}' href='#GeoGrids.latlon_geometry-Tuple{PolyArea{Meshes.𝔼{2}, var"#s12", R, V} where {var"#s12"<:(Cartesian{WGS84Latest, 2, Unitful.Quantity{T, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}} where T), R<:(Meshes.Ring{Meshes.𝔼{2}, var"#s12", V} where V<:(CircularArrays.CircularArray{Point{Meshes.𝔼{2}, var"#s12"}, 1, A} where A<:AbstractArray{Point{Meshes.𝔼{2}, var"#s12"}, 1})), V<:AbstractVector{R}}}'><span class="jlbinding">GeoGrids.latlon_geometry</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
latlon_geometry(poly::PolyArea{𝔼{2},<:CART})
latlon_geometry(multi::Multi{𝔼{2},<:CART})
```


Convert geometries from Cartesian to LatLon coordinate systems.

**Arguments**
- `poly::PolyArea{𝔼{2},<:CART}`: A polygon in Cartesian coordinates.
- `multi::Multi{𝔼{2},<:CART}`: A multi-geometry in Cartesian coordinates.

**Returns**
- `PolyArea` or `Multi`: The converted geometry in LatLon coordinate system.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/helper_func.jl#L115-L127)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.offset_region-Union{Tuple{P}, Tuple{D}, Tuple{GeoRegion{D, P}, Any}} where {D, P}' href='#GeoGrids.offset_region-Union{Tuple{P}, Tuple{D}, Tuple{GeoRegion{D, P}, Any}} where {D, P}'><span class="jlbinding">GeoGrids.offset_region</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
offset_region(originalRegion::GeoRegion, deltaDist; refRadius=constants.Re_mean, magnitude=3, precision=7)
offset_region(originalRegion::PolyRegion, deltaDist; refRadius=constants.Re_mean, magnitude=3, precision=7)
```


Offset a GeoRegion or PolyRegion by a given distance. This function offsets each polygon in the region separately and combines the results into a Multi geometry.

**Arguments**
- `originalRegion::Union{GeoRegion,PolyRegion}`: The original region to be offset.
- `deltaDist`: The distance to offset the region by, in meters. Positive for enlargement, negative for shrinking.
- `refRadius::Float64=constants.Re_mean`: The reference radius to use for the Earth.
- `magnitude::Int=3`: The number of integer digits for IntPoint conversion.
- `precision::Int=7`: The total number of digits to be considered for each coordinate in IntPoint conversion.

**Returns**
- `Multi`: A Multi geometry containing the offset polygons.

**Notes**
- For GeoRegion, only the outer ring of each polygon is considered for offsetting.
- For PolyRegion, if multiple outer rings are produced, inner rings are ignored and separate PolyAreas are created for each outer ring.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/offsetting_func.jl#L1-L21)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.rectgrid-Tuple{Union{Real, Unitful.Quantity{<:Real, <:Any, <:Union{Unitful.FreeUnits{(°,), NoDims, nothing}, Unitful.FreeUnits{(rad,), NoDims, nothing}}}}}' href='#GeoGrids.rectgrid-Tuple{Union{Real, Unitful.Quantity{<:Real, <:Any, <:Union{Unitful.FreeUnits{(°,), NoDims, nothing}, Unitful.FreeUnits{(rad,), NoDims, nothing}}}}}'><span class="jlbinding">GeoGrids.rectgrid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rectgrid(xRes::ValidAngle; yRes::ValidAngle=xRes) -> Array{Point{🌐,<:LatLon{WGS84Latest}}, 2}
```


Create a rectangular grid of latitude and longitude points with a specified grid resolutions. The function validates the input resolutions to ensure they are within the range of `-180°` to `180°`. If a negative resolution is provided, it is converted to a positive value with a warning. The resolution values are used to create a rectangular grid covering the entire range of latitudes `[-90°, 90°]` and longitudes `[-180°, 180°)`. 

**Arguments**
- `xRes::ValidAngle`: The resolution for the latitude grid spacing. This can be a real number (interpreted as degrees) or a `ValidAngle`.
- `yRes::ValidAngle`: The resolution for the longitude grid spacing. This is optional and defaults to `xRes` if not provided. This can be a real number (interpreted as degrees) or a `ValidAngle`.

**Returns**
- A 2D array of `Point{🌐,<:LatLon{WGS84Latest}}` objects representing the grid of latitude and longitude points.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/rect_func.jl#L1-L17)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.vecgrid-Tuple{Union{Real, Unitful.Quantity{<:Real, <:Any, <:Union{Unitful.FreeUnits{(°,), NoDims, nothing}, Unitful.FreeUnits{(rad,), NoDims, nothing}}}}}' href='#GeoGrids.vecgrid-Tuple{Union{Real, Unitful.Quantity{<:Real, <:Any, <:Union{Unitful.FreeUnits{(°,), NoDims, nothing}, Unitful.FreeUnits{(rad,), NoDims, nothing}}}}}'><span class="jlbinding">GeoGrids.vecgrid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
vecgrid(gridRes::ValidAngle) -> Vector{Point{🌐,<:LatLon{WGS84Latest}}}
```


Generate a vector of latitude points from the equator to the North Pole with a specified resolution. The function validates the input resolution to ensure it is within the range of `-90°` to `90°`. If a negative resolution is provided, it is converted to a positive value with a warning. The resolution value is then used to create a vector of latitude points ranging from `0°` to `90°` (the North Pole). Each latitude point is represented as a `LatLon` object with a fixed longitude of `0°`.

**Arguments**
- `gridRes::ValidAngle`: The resolution for the latitude grid spacing. This can be a real number (interpreted as degrees) or a `ValidAngle`.

**Returns**
- A vector of `Point{🌐,<:LatLon{WGS84Latest}}` objects representing latitude points from the equator (0°) to the North Pole (90°) with the specified resolution.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/rect_func.jl#L68-L84)

</details>


## Private
<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._adapted_icogrid-Tuple{Number}' href='#GeoGrids._adapted_icogrid-Tuple{Number}'><span class="jlbinding">GeoGrids._adapted_icogrid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_adapted_icogrid(radius::Number; correctionFactor=3/2)
```


The `_adapted_icogrid` function generates an icosahedral grid with a specified radius. It defines the separation angle for the icosahedral grid using a correction factor to adapt the cell centers&#39; distances, ensuring the grid is appropriate for the desired scale.

**Arguments**
- `radius::Number`: The radius used to define the separation angle for the icosahedral grid. This radius helps determine the distance between the grid points.

**Keyword Arguments**
- `correctionFactor=1.2`: The correction factor used to adapt the cell centers&#39; distances to ensure the grid is appropriate for the desired scale.

**Returns**
- `grid`: The generated icosahedral grid based on the calculated separation angle. The specific structure and format of the returned grid depend on the `icogrid` function being used.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L229-L245)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._add_angular_offset-Tuple{Any, Any}' href='#GeoGrids._add_angular_offset-Tuple{Any, Any}'><span class="jlbinding">GeoGrids._add_angular_offset</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_add_angular_offset(inputθϕ, offsetθϕ) -> NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}
```


Add an angular offset to given spherical coordinates.

**Arguments**
- `inputθϕ::NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}`: The input spherical coordinates with components `θ` (polar angle) and `ϕ` (azimuthal angle) in radians.
- `offsetθϕ::NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}`: The offset spherical coordinates with components `θ` (polar angle) and `ϕ` (azimuthal angle) in radians.

**Returns**
- `NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}`: The new spherical coordinates after applying the angular offset, with components `θ` and `ϕ` in radians.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/helper_func.jl#L39-L50)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._fibonaccisphere_classic_partial-Tuple{Any}' href='#GeoGrids._fibonaccisphere_classic_partial-Tuple{Any}'><span class="jlbinding">GeoGrids._fibonaccisphere_classic_partial</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_fibonaccisphere_classic_partial(N; spheRadius=1.0, pointsToCheck::Int=50)
```


**Arguments:**
- `N`: an integer representing the number of points to generate on the surface of the sphere.
- `spheRadius`: (optional) a float representing the radius of the sphere.
- `pointsToCheck`: (optional) an integer representing the number of points to return starting from the first generated.

**Output:**
- `points`: an array of 3D points on the surface of the sphere represented as SVector{3}.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/ico_func.jl#L194-L204)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._find_min_separation_angle-Tuple{Any}' href='#GeoGrids._find_min_separation_angle-Tuple{Any}'><span class="jlbinding">GeoGrids._find_min_separation_angle</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_find_min_separation_angle(points)
```


This function takes an array of 3D Cartesian points as input and computes the smallest angle between any two points in the array. It does this by iterating over all unique pairs of points in the array and computing the angle between them using the `angle``function. The smallest angle encountered during the iteration is stored in the variable`sep` and returned as the output of the function.

**Arguments:**
- `points`: an array of 3D points in the Cartesian plane represented as Tuples, Arrays, SVectors.

**Output:**
- `sep`: the smallest angle between any two points in the input array, returned as Uniful.Quantity in degrees.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/ico_func.jl#L157-L172)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._gen_regular_lattice-Union{Tuple{T}, Tuple{T, Any, Any}} where T' href='#GeoGrids._gen_regular_lattice-Union{Tuple{T}, Tuple{T, Any, Any}} where T'><span class="jlbinding">GeoGrids._gen_regular_lattice</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_gen_regular_lattice(dx::T, dy, ds; x0=zero(T), y0=zero(T), M::Int=100, N::Int=M) where {T}
```


The `_gen_regular_lattice` function generates a regular lattice of points in a two-dimensional space. The lattice is defined by its spacing in the x and y directions (`dx` and `dy`), an additional offset (`ds`), and optional starting positions (`x0` and `y0`). The lattice spans `2M + 1` rows and `2N + 1` columns centered around the origin.

**Arguments**
- `dx::T`: The spacing between points in the x direction.
- `dy::T`: The spacing between points in the y direction.
- `ds::T`: The additional offset in the x direction per row.
- `x0::T`: The x-coordinate of the starting position. Default is `zero(T)`.
- `y0::T`: The y-coordinate of the starting position. Default is `zero(T)`.
- `M::Int`: The number of points in the x direction from the center. Default is 100.
- `N::Int`: The number of points in the y direction from the center. Default is equal to `M`.

**Returns**
- `Array{SVector{2,T},2}`: A 2D array of points represented as static vectors (`SVector{2,T}`) from the `StaticArrays` package. Each point is in the form `(x, y)`.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L152-L172)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._get_theta_phi-Tuple{Number, Number}' href='#GeoGrids._get_theta_phi-Tuple{Number, Number}'><span class="jlbinding">GeoGrids._get_theta_phi</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_get_theta_phi(k::Number, N::Number) -> Tuple{Number, Number}
```


Calculate the spherical coordinates θ (theta) and ϕ (phi) for a given index `k` and total number of points `N` using the Golden Ratio method.

**Arguments**
- `k::Number`: The index of the point for which the spherical coordinates are to be calculated. This should typically be an integer between `0` and `N-1`.
- `N::Number`: The total number of points for which the spherical coordinates are to be calculated.

**Returns**
- A tuple `(θ, ϕ)` where:
    - `θ::Number`: The longitude angle in rad, ranging from `[0, 360]`.
    - `ϕ::Number`: The latitude angle in rad, ranging from `[0, 180]` from the North Pole.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/ico_func.jl#L220-L234)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._hex_tesselation_centroids-Tuple{Point{Meshes.🌐, <:CoordRefSystems.GeodeticLatLon{WGS84Latest}}, Number}' href='#GeoGrids._hex_tesselation_centroids-Tuple{Point{Meshes.🌐, <:CoordRefSystems.GeodeticLatLon{WGS84Latest}}, Number}'><span class="jlbinding">GeoGrids._hex_tesselation_centroids</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_hex_tesselation_centroids(origin::Point{🌐,<:LatLon{WGS84Latest}}, radius::Number; direction::Symbol=:pointy, Re::Number=constants.Re_mean, kwargs_lattice...)
```


This function generates the centroids of a hexagonal tessellation on the Earth&#39;s surface, centered at the origin. The tessellation is created based on a given radius and direction. The function converts the offsets of the hexagonal grid to latitude and longitude coordinates.

**Arguments**
- `origin::Point{🌐,<:LatLon{WGS84Latest}}`: The lat-lon coordinates of the center of the tessellation. 
- `radius::Number`: The radius of the hexagons in the tessellation in meters.
- `direction::Symbol`: The direction of the hexagons, either `:pointy` (default) or `:flat`.
- `Re::Number`: The mean radius of the Earth in meters (default is `constants.Re_mean`).
- `kwargs_lattice...`: Additional keyword arguments for the hexagonal lattice generation.

**Returns**
- `Vector{Point{🌐,<:LatLon{WGS84Latest}}}`: A vector of `Point{🌐,<:LatLon{WGS84Latest}}` objects representing the centroids of the hexagonal tessellation in latitude and longitude.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/tessellation_func.jl#L261-L278)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._icogrid-Tuple{Int64}' href='#GeoGrids._icogrid-Tuple{Int64}'><span class="jlbinding">GeoGrids._icogrid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_icogrid(N::Int)
```


This function generates `N` uniformly distributed points on the surface of a unitary sphere using the classic Fibonacci Spiral method as described in [1]. Contrary to the Ichosahedral grid generation process, with the Fibonacci Spiral method it is possible to generate a grid of points uniformly distributed in the area for a generic `N` value. As a drawback, the structure of the points do not follow a &quot;perfect simmetry&quot; however, the density of points in the area is preserved quasi-constant.

**Arguments:**
- `N::Int`: The number of points to generate.
- `coord::Symbol`: The type of coordinates of generated points (`:sphe` | `:cart`).
- `radius`: the sphere radius in meters (unitary as default)

**Output:**
- `pointsVec`: a `Vector{SVector}` of the generated points. Each element corresponds to a point on the surface of the sphere, the SVector contains either the x, y, and z (:cart) or lat, lon (:sphe) (LAT=x, LON=y) in rad coordinates of the point.

**References**
1. http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/ico_func.jl#L48-L69)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._offset_ring-Union{Tuple{T}, Tuple{Meshes.Ring{Meshes.🌐, CoordRefSystems.GeodeticLatLon{WGS84Latest, Unitful.Quantity{T, NoDims, Unitful.FreeUnits{(°,), NoDims, nothing}}}, CircularArrays.CircularArray{Point{Meshes.🌐, CoordRefSystems.GeodeticLatLon{WGS84Latest, Unitful.Quantity{T, NoDims, Unitful.FreeUnits{(°,), NoDims, nothing}}}}, 1, Array{Point{Meshes.🌐, CoordRefSystems.GeodeticLatLon{WGS84Latest, Unitful.Quantity{T, NoDims, Unitful.FreeUnits{(°,), NoDims, nothing}}}}, 1}}}, Any}} where T' href='#GeoGrids._offset_ring-Union{Tuple{T}, Tuple{Meshes.Ring{Meshes.🌐, CoordRefSystems.GeodeticLatLon{WGS84Latest, Unitful.Quantity{T, NoDims, Unitful.FreeUnits{(°,), NoDims, nothing}}}, CircularArrays.CircularArray{Point{Meshes.🌐, CoordRefSystems.GeodeticLatLon{WGS84Latest, Unitful.Quantity{T, NoDims, Unitful.FreeUnits{(°,), NoDims, nothing}}}}, 1, Array{Point{Meshes.🌐, CoordRefSystems.GeodeticLatLon{WGS84Latest, Unitful.Quantity{T, NoDims, Unitful.FreeUnits{(°,), NoDims, nothing}}}}, 1}}}, Any}} where T'><span class="jlbinding">GeoGrids._offset_ring</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_offset_ring(ring::Ring{🌐,<:LatLon{WGS84Latest}}, delta; magnitude=3, precision=7)
```


Offset a ring by a given delta value. This function uses the Clipper library for polygon offsetting. It may return multiple rings even when starting from a single Ring.

**Arguments**
- `ring::Ring{🌐,<:LatLon{WGS84Latest}}`: The ring to be offset.
- `delta`: The distance to offset the ring by. Positive for enlargement, negative for shrinking.
- `magnitude::Int=3`: The number of integer digits for IntPoint conversion.
- `precision::Int=7`: The total number of digits to be considered for each coordinate in IntPoint conversion.

**Returns**
- Vector of `Ring{🌐,<:LatLon{WGS84Latest}}`: The resulting offset rings.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/offsetting_func.jl#L93-L108)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._points_required_for_separation_angle-Tuple{Union{Real, Unitful.Quantity{<:Real, <:Any, <:Union{Unitful.FreeUnits{(°,), NoDims, nothing}, Unitful.FreeUnits{(rad,), NoDims, nothing}}}}}' href='#GeoGrids._points_required_for_separation_angle-Tuple{Union{Real, Unitful.Quantity{<:Real, <:Any, <:Union{Unitful.FreeUnits{(°,), NoDims, nothing}, Unitful.FreeUnits{(rad,), NoDims, nothing}}}}}'><span class="jlbinding">GeoGrids._points_required_for_separation_angle</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_points_required_for_separation_angle(angle)
```


This function computes the minimum number of points required on the surface of a sphere to achieve a desired separation angle between any two adjacent points. The function uses the bisection method to find the minimum number of points needed and returns the higher end of the precision. Instead of checking all the possible pairs of the N points generated with the Fibonacci Spiral, only the first `pointsToCheck` are checked in order to evaluate the minimum separation angle.

**Arguments:**
- `sepAng`: a float representing the desired separation angle between two adjacent points on the surface of the sphere.
- `spheRadius`: an optional float representing the radius of the sphere. If not provided, it defaults to 1.0.
- `pointsToCheck`: an optional integer representing the number of points to generate on the surface of the sphere. If not provided, it defaults to 50.
- `maxPrec`: an optional integer representing the maximum precision for the number of points generated on the surface of the sphere. If not provided, it defaults to 10^7.
- `tol`: an optional integer representing the tolerance for the bisection method used to find the minimum number of points needed to achieve the desired separation angle. If not provided, it defaults to 10.

**Output:**
- `Ns[2]`: an integer representing the minimum number of points required on the surface of the sphere to achieve the desired separation angle.
- `thisSep`: an Uniful.Quantity in degrees representing the separation angle between two adjacent points on the surface of the sphere.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/ico_func.jl#L91-L112)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids._wrap_latlon-Tuple{Number, Number}' href='#GeoGrids._wrap_latlon-Tuple{Number, Number}'><span class="jlbinding">GeoGrids._wrap_latlon</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_wrap_latlon(lat::Number, lon::Number)
```


The `_wrap_latlon` function normalizes and wraps geographic coordinates, latitude (`lat`) and longitude (`lon`). It ensures that the latitude is within the range [-90, 90] degrees and the longitude is within the range [-180, 180) degrees. This function is useful for handling geographic data where coordinates might exceed their typical bounds.

**Arguments**
- `lat::Number`: The latitude value to be normalized and wrapped, expressed in degrees.
- `lon::Number`: The longitude value to be normalized and wrapped, expressed in degrees.

**Returns**
- `Tuple{Number, Number}`: A tuple `(lat, lon)` in degrees where `lat` is in the range [-90, 90] and `lon` is in the range [-180, 180).


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/helper_func.jl#L1-L16)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.fibonaccisphere_alternative1-Tuple{Int64}' href='#GeoGrids.fibonaccisphere_alternative1-Tuple{Int64}'><span class="jlbinding">GeoGrids.fibonaccisphere_alternative1</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fibonaccisphere_alternative1(N::Int)
```


This function generates points on the surface of a unit sphere using the Fibonacci spiral method. The function takes an integer `N` as an input, which specifies the number of points to be generated.

**Arguments**
- `N::Int`: The number of points to generate. This is an integer value.

**Output**
- `N x 3` matrix containing the generated points. Each row of the matrix corresponds to a point on the surface of the sphere, and the columns correspond to the x, y, and z coordinates of the point.


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/ico_func.jl#L307-L318)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoGrids.fibonaccisphere_optimization1-Tuple{Int64}' href='#GeoGrids.fibonaccisphere_optimization1-Tuple{Int64}'><span class="jlbinding">GeoGrids.fibonaccisphere_optimization1</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fibonaccisphere_optimization1(N)
```


Create a set of `N` uniformly distributed points on the surface of a sphere using the Fibonacci spiral method optimized as described in [1].

This method is called Offset Fibonacci Lattice which is one method to optimize the minimum nearest-neighbor distance. We need to move (offset) all the points slightly farther away from the poles. This of course means, that almost all of them become slightly closer together. Offsetting the points of the Fibonacci lattice slightly away from the poles produces a packing that is up to 8.3% tighter than the canonical Fibonacci lattice.

For `n>100`, an improvement can be made beyond this, by initially placing a point at each pole, and then placing the remaining `n-2` points. This not only (very sightly) improves minimal nearest packing, but it also prevents a large gap at each pole.

**Arguments**
- `N::Int`: The number of points to generate.

**Output**
- `points::Matrix{Float64}`: A `N`x`3` matrix where each row corresponds to a point `(x,y,z)` on the surface of the unitary sphere.

**References**
1. http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/


[source](https://github.com/https://blob/f2f9584d448cd6444b5ca41ab237fa732af12ede/src/ico_func.jl#L248-L274)

</details>

