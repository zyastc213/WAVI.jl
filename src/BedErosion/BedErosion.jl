struct BedErosion_obj{T <: Real} <: AbstractBedErosion
    Kg :: T 
    l :: T 
    erosion_rate :: Array{T,2}
    erosion_accumulation :: Array{T,2}
    change_bedrock :: Bool       
end

function BedErosion_obj(;
                        Kg = 5e-5,
                        l = 1.0,
                        erosion_rate = nothing,
                        erosion_accumulation = nothing,
                        change_bedrock = true
                        )
    ~(erosion_rate === nothing) || throw(ArgumentError("You must input a initial erosion_rate such as zeros(nx,ny)"))
    ~(erosion_accumulation === nothing) || throw(ArgumentError("You must input a initial erosion_accumulation such as zeros(nx,ny)"))
    return BedErosion_obj(Kg, l, erosion_rate, erosion_accumulation, change_bedrock)
end

function update_BedErosion_obj!(bed_erosion::BedErosion_obj, fields)
    @unpack bed_speed = fields.gh
    @unpack Kg, l, erosion_rate, erosion_accumulation = bed_erosion
    erosion_rate .=  Kg * bed_speed .^ l
    erosion_accumulation .+= erosion_rate
    println("bed erosion is calculated")
end

function update_bed_elevation!(bed_erosion::BedErosion_obj, fields)
    @unpack gh=fields
    gh.b .-= bed_erosion.erosion_rate
end
