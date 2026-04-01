struct Clim_obj{T <: Real}
    air_temp :: Array{T,3}               # air temperature (:,:,12)
    precip :: Array{T,3}               # precipitation (:,:,12)
    γₚ :: T
    h_ref :: Array{T,2}
end

function Clim_obj(;
                  air_temp=nothing,
                  precip=nothing,
                  γₚ=0.75,
                  h_ref=nothing,
                  )
    
    ~(air_temp === nothing) || throw(ArgumentError("You must input a initial air_temp such as zeros(nx,ny,12)"))
    ~(precip === nothing) || throw(ArgumentError("You must input a initial precip such as zeros(nx,ny,12)"))
    ~(h_ref === nothing) || throw(ArgumentError("You must input a initial h_ref such as your 'bed_elevation' or 'bed_elevation + ice_thickness' "))

    return Clim_obj(air_temp, precip, γₚ, h_ref)
end

# function update_Clim_obj!(clim::Clim_obj,fields)
#     @unpack gh = fields
#     println("precip_max_before: $(maximum(clim.precip))")
#     @. clim.precip *= exp(-1 * clim.γₚ * max(0, gh.s - clim.h_ref))  # adjust the precipitation according to current elevation
#     println("precip_max: $(maximum(clim.precip))")
#     println("Precipitation is corrected")
#     return nothing
# end