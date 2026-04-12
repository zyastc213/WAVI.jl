using SpecialFunctions
# Distributions

struct PDD_obj{T <: Real, X} <: AbstractSurfaceProcess
    T_min :: T                    # air temperature all precipitation as snowfall
    T_max :: T                    # air temperature all precipitation as rainfall
    σ_T :: X                      # standard deviation of daily temperature variation
    σ_lapse_lat_base :: T         # the latitude where = σ_T 
    σ_lapse_lat_rate :: T         # the lapse rate of σ_T
    θr :: T                       # refreeze fraction
    Fi :: T                       # factor ice 
    Fs :: T                       # factor snow
    SMB_max :: T
    PDD_scheme :: Bool            # choose if you want to use PDD scheme
    SMB :: Array{T,2}             # surface mass balance
    IceMelt :: Bool
end 

function PDD_obj(;
            T_min=273.15,
            T_max=275.15,
            σ_T=5.0,
            σ_lapse_lat_base=72.0,
            σ_lapse_lat_rate=0.0,
            θr=0.6,
            Fi=0.00879120879120879,
            Fs=0.0032967032967033,
            SMB_max=999.0,
            PDD_scheme = true,
            SMB = nothing,
            IceMelt = false
            )
 
    @assert T_max > T_min
    @assert any(σ_T .> 0)
    @assert 0 <= θr <= 1 
    @assert Fi > 0 && Fs >0
    ~(SMB === nothing) || throw(ArgumentError("You must input a initial SMB such as zeros(nx,ny)"))

    return PDD_obj(T_min, T_max, σ_T, σ_lapse_lat_base, σ_lapse_lat_rate, θr, Fi, Fs, SMB_max, PDD_scheme, SMB, IceMelt)
end

# calculate the number of PDDs (in:T out:PDD) 
function compute_PDD_CalovGreve05(σ_T,
                                  T)
    N = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    PDD = zeros(size(T[:,:,1]))
    for i=1:12
        Z = T[:,:,i] ./ (sqrt(2) .* σ_T[:,:,i])
        PDD .+= (σ_T[:,:,i] ./ sqrt(2*π) .* exp.(-Z.^2) + T[:,:,i] ./ 2 .* erfc.(-Z)) .* N[i]
    end
    println("PDD_max: $(maximum(PDD))")
    return PDD
end

# update snow depth (in:T,P out:h_snow)
function compute_snow_depth(T_min,
                            T_max,
                            density_ice,
                            T,
                            P,
                            P_lapse_rate)
    h_snow = zeros(size(P[:,:,1]))
    P_tmp = zeros(size(P[:,:,1]))
    N = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for i in 1:12
        # @views P[:,:,i] ./= density_ice/(N[i]*24*60*60)                    # convert "kg m-2 second-1" to "m month-1" (snow depth rate)
        # # use linear calculation to determine how much precipitation interpreted as snow fall
        # @views @. P[:,:,i] = ifelse(T[:,:,i] <= T_min, P[:,:,i],
        #                      ifelse(T[:,:,i] >= T_max, 0.0,
        #                           P[:,:,i] * (T_max - T[:,:,i]) / (T_max - T_min)))
        # h_snow .+= P[:,:,i]
        P_tmp .= P[:,:,i] ./ density_ice .* (N[i]*24*60*60) .* P_lapse_rate       # convert "kg m-2 second-1" to "m month-1" (snow depth rate)
        # use linear calculation to determine how much precipitation interpreted as snow fall
        @. P_tmp = ifelse(T[:,:,i] <= T_min, P_tmp,
                             ifelse(T[:,:,i] >= T_max, 0.0,
                                  P_tmp * (T_max - T[:,:,i]) / (T_max - T_min)))
        h_snow .+= P_tmp
    end
    println("hsnow_max: $(maximum(h_snow))")
    return h_snow
end

function compute_runoff(PDD,
                        Fs,
                        Fi,
                        θr,
                        h_snow,
                        IceMelt)
    snow_melt = min.(h_snow, PDD .* Fs)
    ice_melt = (PDD .- snow_melt / Fs) .* Fi
    if !IceMelt
        refreeze = snow_melt .* θr
    else
        refreeze = (snow_melt .+ ice_melt) .* θr
    end
    runoff = snow_melt .+ ice_melt .- refreeze
    return runoff
end

function compute_surface_mass_balance(T_min, T_max, σ_T, θr, Fi, Fs, SMB_max, IceMelt, density_ice, air_temp, precip, P_lapse_rate)
    PDD = compute_PDD_CalovGreve05(σ_T, air_temp.-273.15) # convert Kelvin to °C
    h_snow = compute_snow_depth(T_min, T_max, density_ice, air_temp, precip, P_lapse_rate)
    runoff = compute_runoff(PDD, Fs, Fi, θr, h_snow, IceMelt)
    SMB = h_snow - runoff # check if the formula is correct???? 
    if SMB_max != 999.0
        SMB .= ifelse.(SMB.>0, SMB_max .* tanh.(SMB./SMB_max), SMB)
    end
    println("SMBmax: $(maximum(SMB)), SMBmin: $(minimum(SMB))")
    return SMB
end

function update_PDD_obj!(surface_melt::PDD_obj, clim::Clim_obj, params, fields)
    @unpack air_temp, precip, h_ref, γₚ = clim
    @unpack T_min, T_max, σ_T, θr, Fi, Fs, SMB_max, IceMelt=surface_melt
    P_lapse_rate = zeros(size(h_ref))
    P_lapse_rate .= exp.(-1 * γₚ .* max.(0, (fields.gh.s .- h_ref)))
    println("P_lapse_rate_max: $(minimum(P_lapse_rate))")
    surface_melt.SMB .= compute_surface_mass_balance(T_min, T_max, σ_T, θr, Fi, Fs, SMB_max, IceMelt, params.density_ice, air_temp, precip, P_lapse_rate)
    println("Surface mass balance is updated")
    return nothing
end

"""
refreeze_fraction = 0.6
melt = snow_melt + firn_melt + ice_melt
PDD_excess = PDD * snow_melt / Fₛ
ice_melt = PDD_excess * Fᵢ
refreeze = (snow_melt + ice_melt) * refreeze_fraction
runoff = snow_melt + ice_melt - refreeze
smb = accumulation - runoff
H_snow .-= snow_melt
"""
