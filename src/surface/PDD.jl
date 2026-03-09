using SpecialFunctions,Distributions

struct PDD_obj{T <: Real} <: AbstractSurfaceProcess
    T_min :: T                    # air temperature all precipitation as snowfall
    T_max :: T                    # air temperature all precipitation as rainfall
    σ_T :: T                      # standard deviation of daily temperature variation
    σ_lapse_lat_base :: T         # the latitude where = σ_T 
    σ_lapse_lat_rate :: T         # the lapse rate of σ_T
    θr :: T                       # refreeze fraction
    Fi :: T                       # factor ice 
    Fs :: T                       # factor snow
    T :: Array{T,3}               # air temperature (:,:,12)
    P :: Array{T,3}               # precipitation (:,:,12)
    SMB :: Array{T,2}             # surface mass balance
end 

function PDD_obj(;
            T_min=275.15,
            T_max=273.15,
            σ_T=5.0,
            σ_lapse_lat_base=72.0,
            σ_lapse_lat_rate=0.0,
            θr=0.6,
            Fi=0.00879120879120879,
            Fs=0.0032967032967033,
            T=air_temp,
            P=precip)
    SMB = compute_surface_mass_balance()
    return PDD_obj(T_min, T_max, σ_T, σ_lapse_lat_base, σ_lapse_lat_rate, θr, Fi, Fs, T, P, SMB)
end

# calculate the number of PDDs (in:T out:PDD) 
function compute_PDD_CalovGreve05(sigma_T,
                                  T_ac)
    N = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for i=1:12
        Z = T_ac[:,:,i] ./ (sqrt(2) * sigma_T)
        PDD .+= (sigma_T / sqrt(2*π) .* exp.(-Z.^2) + T_ac[:,:,i] ./ 2 .* erfc.(-Z)) .* N[i]
    end
    return PDD
end

# update snow depth (in:T,P out:h_snow)
function compute_snow_depth(T_min,
                            T_max,
                            T,
                            P)
    @unpack density_ice=params.density_ice
    for i in 1:12
        @views P[:,:,i] ./= density_ice                    # convert "kg m-2 second-1" to "m second-1" (snow depth rate)
        # use linear calculation to determine how much precipitation interpreted as snow fall
        @views @. P[:,:,i] = ifelse(T[:,:,i] <= T_min, P[:,:,i],
                             ifelse(T[:,:,i] >= T_max, 0.0,
                                  P[:,:,i] * (T_max - T[:,:,i]) / (T_max - T_min)))
        h_snow .+= P[:,:,i]
    end
    return h_snow
end

function compute_runoff(PDD,
                        Fs,
                        Fi,
                        θr,
                        h_snow)
    snow_melt = max.(h_snow, PDD .* Fs)
    ice_melt = (PDD .- snow_melt / Fs) .* Fi
    refreeze = (snow_melt .+ ice_melt) .* θr
    runoff = snow_melt .+ ice_melt .- refreeze
    return runoff
end

function compute_surface_mass_balance(PDD_obj,T,P)
    @unpack σ_T = PDD_obj
    @unpack h = model.fields.gh

    PDD = compute_PDD_CalovGreve05(σ_T, T)
    h_snow = compute_snow_depth(T_min, T_max, T, P)
    runoff = compute_runoff(PDD, Fs, Fi, θr, h_snow)
end

refreeze_fraction = 0.6
melt = snow_melt + firn_melt + ice_melt
PDD_excess = PDD * snow_melt / Fₛ
ice_melt = PDD_excess * Fᵢ
refreeze = (snow_melt + ice_melt) * refreeze_fraction
runoff = snow_melt + ice_melt - refreeze
smb = accumulation - runoff
H_snow .-= snow_melt

"""
OUTLINE

SMB(accumulation(P(:,:,(12)),T(:,:,(12)),T_max,T_min), runoff(PDD(T(:,:,12),sigma_T(:,:,12)), Fₛ, Fᵢ,θᵣ))

accumulation(P(:,:,12),T(:,:,12),T_max,T_min) => update_snow_depth!
PDD(T(:,:,12),sigma_T(:,:,12)) => update_PDD_CalovGreve05!


"""

function update_accumulation_rate!()
    @unpack accumulation_rate = params.accumulation_rate
    accumulation_rate = A
end

"""
function calculate_PDD(PDD::Float64,T_d::Float64)
    N = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for i=1:12
        T_d = T_mm[i] + rand(Normal(0,sigma_T))
        PDD[i] = N[i] * (T_d)


end

function PDD_scheme()
    potential_melt_snow = Fₛ * PDD_days
end
"""


"""
function pdd_model()


end

function pdd_model(T_daily::Vector{Float64}, P_daily::Vector{Float64};
                   T_snow_threshold::Float64 = 0.0,
                   PDD_factor_snow::Float64 = 0.003,
                   PDD_factor_ice::Float64 = 0.008,
                   refreeze_fraction::Float64 = 0.6,
                   sigma_T::Float64 = 5.0)
    
    n_days = length(T_daily)
    @assert length(P_daily) == n_days  # check the length of T and P
    
    # 初始化结果数组
    snowfall = zeros(n_days)      # 降雪量 [m w.e.]
    rainfall = zeros(n_days)      # 降雨量 [m w.e.]
    melt_snow = zeros(n_days)     # 雪的融化量 [m w.e.]
    melt_ice = zeros(n_days)      # 冰的融化量 [m w.e.]
    accumulation = zeros(n_days)  # 积累量 [m w.e.]
    ablation = zeros(n_days)      # 消融量 [m w.e.]
    snow_depth = zeros(n_days+1)  # 雪深 [m w.e.]
    
    # 初始雪深
    snow_depth[1] = 0.0
    
    for day in 1:n_days
        T = T_daily[day]
        P = P_daily[day]
        
        # 1. 降水相态分离
        if T <= T_snow_threshold
            snowfall[day] = P
            rainfall[day] = 0.0
        else
            snowfall[day] = 0.0
            rainfall[day] = P
        end
        
        # 2. 计算正度日（考虑温度波动）
        PDD = calculate_pdd(T, sigma_T)
        
        # 3. 融化计算（先融化雪，再融化冰）
        if PDD > 0.0
            # 先计算雪的融化潜力
            potential_melt_snow = PDD_factor_snow * PDD
            
            if snow_depth[day] > 0.0
                # 有积雪：先融化雪
                melt_snow[day] = min(potential_melt_snow, snow_depth[day])
                remaining_melt_potential = potential_melt_snow - melt_snow[day]
                
                # 剩余融化潜力用于融化冰
                melt_ice[day] = PDD_factor_ice * remaining_melt_potential
            else
                # 无积雪：直接融化冰
                melt_snow[day] = 0.0
                melt_ice[day] = PDD_factor_ice * PDD
            end
        else
            melt_snow[day] = 0.0
            melt_ice[day] = 0.0
        end
        
        # 4. 积累和消融
        accumulation[day] = snowfall[day]
        ablation[day] = melt_snow[day] + melt_ice[day] - refreeze_fraction * rainfall[day]
        
        # 5. 更新雪深
        snow_depth[day+1] = max(0.0, snow_depth[day] + snowfall[day] - melt_snow[day])
    end

"""