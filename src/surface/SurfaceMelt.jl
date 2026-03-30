include("./Clim.jl")
include("./PDD.jl")


function update_surface_melt!(surface_melt::PDD_obj, clim::Clim_obj, params, fields)
    @unpack gh = fields
    if !surface_melt.PDD_scheme
        gh.accumulation .= params.accumulation_rate
        println("PDD_scheme is off")
    else
        # update_Clim_obj!(clim, fields) # correct the precipitation
        update_PDD_obj!(surface_melt, clim, params,fields) # change the SMB
        gh.accumulation .= surface_melt.SMB
    end
end