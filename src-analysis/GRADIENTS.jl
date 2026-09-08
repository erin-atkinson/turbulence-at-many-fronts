include("terms/terms.jl")

fields = (:u, :v, :w, :b)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

vorticity_x = Field(VorticityX(u_bar, v_bar, w_bar))
vorticity_y = Field(VorticityY(u_bar, v_bar, w_bar))
vorticity_z = Field(VorticityZ(u_bar, v_bar, w_bar))
vorticity = (; vorticity_x, vorticity_y, vorticity_z)

M² = Field(∂x(b_bar))
N² = Field(∂z(b_bar))
buoyancy = (; M², N²)

S² = Field(∂z(u_bar)^2 + ∂z(v_bar)^2)
Ri = Field(Richardson(u_bar, v_bar, b_bar))
Rib = Field(BalancedRichardson(b_bar, sp))
shear = (; S², Ri, Rib)

q = Field(PotentialVorticity(u_bar, v_bar, w_bar, b_bar, sp))
potential_vorticity = (; q)

skip_update = filter(a->a ∉ fields, keys(input_fields))
dependency_fields = merge(vorticity, buoyancy, shear, potential_vorticity)
output_fields = dependency_fields
