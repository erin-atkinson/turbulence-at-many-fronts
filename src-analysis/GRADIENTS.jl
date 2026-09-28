include("terms/terms.jl")

fields = (:u, :v, :w, :b)

for ξ in fields
    @eval $ξ = input_fields.$ξ
end

vorticity_x = Field(VorticityX(u, v, w))
vorticity_y = Field(VorticityY(u, v, w))
vorticity_z = Field(VorticityZ(u, v, w))
vorticity = (; vorticity_x, vorticity_y, vorticity_z)

M² = Field(∂x(b))
N² = Field(∂z(b))
buoyancy = (; M², N²)

S² = Field(∂z(u)^2 + ∂z(v)^2)
Ri = Field(Richardson(u, v, b))
Rib = Field(BalancedRichardson(b, sp))
shear = (; S², Ri, Rib)

q = Field(PotentialVorticity(u, v, w, b, sp))
potential_vorticity = (; q)

skip_update = filter(a->a ∉ fields, keys(input_fields))
dependency_fields = merge(vorticity, buoyancy, shear, potential_vorticity)
output_fields = dependency_fields
