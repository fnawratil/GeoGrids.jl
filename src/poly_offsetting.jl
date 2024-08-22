mutable struct GeoRegionEnlarged <: AbstractRegion
    standardRegion::GeoRegion
    enlargedRegion::GeoRegion

    function GeoRegionEnlarged(standardRegion::GeoRegion, delta::Real)
        # Add code for enlargement of the region

        # new(standardRegion, enlargedRegion)
    end
end

function _offset_polygon(polygon, delta)
    # co = ClipperOffset()
    # add_path!(co, polygon, JoinTypeMiter, EndTypeClosedPolygon)
    # offset_polygons = execute(co, delta)
    # return offset_polygons

    # get_lat(vertices(r.domain[1].latlon.geoms[2])[1])


    # intPolygon = map(vertices(originalPolygon)) do vertex
    #     IntPoint(vertex.x, vertex.y, 3, 6)
    # end
    # polygon = IntPoint[]
	# push!(polygon, IntPoint(0.0,0.0,3,6)) # 3 integer, 6 total values preserved
	# push!(polygon, IntPoint(0.0,1.0,3,6))
	# push!(polygon, IntPoint(1.0,1.0,3,6))
	# push!(polygon, IntPoint(1.0,0.0,3,6))
	# push!(polygon, IntPoint(0.0,0.0,3,6))
	
	# co = ClipperOffset()
	# add_path!(co, polygon, JoinTypeMiter, EndTypeClosedPolygon)
	# # Need to use quite high number of digit to preserve in order to fine tune the sizing.
	# # Use JoinTypeMiter for more accurate tipe of representation (less distortion of polygon shape, visible especially for simple polygons)
	# # Convert back values using tofloat()
	# offset_polygons = execute(co, 500.0)
end

function offset_region(originalRegion::GeoRegion, delta_km; refRadius=constants.Re_mean, magnitude=3, precision=7)
    # `magnitude` represents the number of integer digits while `precision` the
    # total number of digits that will be considered for each of the coordinates
    # for the `IntPoint` conversion. Look at Clipper documentation for more
    # details. `delta` is the distance to offset the polygon by, it is a
    # positive value for enlargement or negative number for shrinking. The value
    # should be expressed in km.

    # Compute the delta value to be used in the offsetting process
    delta = delta_km / refRadius
    intDelta = IntPoint(delta, delta, magnitude, precision)

    numCountries = length(originalRegion.domain) # Number of Countries in GeoRegion
    for idxCountry in 1:numCountries
        country = originalRegion.domain[idxCountry]
        (; admin, latlon, resolution, table_idx, valid_polyareas) = country

        offsetPolyAreas = map(latlon.geoms) do geom
            _offset_polygon(geom, intDelta)
        end
    end

 

end


# //NOTE: NExt Step: 
# [x] Use Clipping.jl for offsetting the polygon
# [x] Starting from the example notebook `polygon_offset.jl`, understand the relation between the value δ and the increase/decrease in polygon area
# [] Test the package with CountriesBorders to see if it works with those domains
# [] Map the value of δ to a Lat-Lon quantity for enlargement/decrease in area
# [] Write an interface function between CountriesBorders/Meshes and Clipping.jl

# [] Investigate alternative packages for polygon offsetting:
#    - GeometryBasics.jl: Provides basic geometric types and operations
#    - Meshes.jl: Offers advanced mesh processing capabilities
#    - LazySets.jl: Includes polygon operations and may support offsetting
#    - Clipper.jl: A Julia wrapper for the Clipper library, which supports polygon offsetting
# [] Compare the features and performance of these packages with Clipping.jl
# [] Choose the most suitable package for our polygon offsetting needs




