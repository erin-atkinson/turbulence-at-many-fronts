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


vorticity_x = Field(-∂z(v_bar))
vorticity_y = Field(∂z(u_bar) - ∂x(w_bar))
vorticity_z = Field(∂x(v_bar))
vorticity = (; vorticity_x, vorticity_y, vorticity_z)


M² = Field(∂x(b_bar))
N² = Field(∂z(b_bar))
buoyancy = (; M², N²)

S² = Field(∂z(u_bar)^2 + ∂z(v_bar)^2)
Ri = Field(N² / S²)
shear = (; S², Ri)

skip_update = filter(a->a ∉ fields, keys(input_fields))
dependency_fields = merge(vorticity, buoyancy, shear)
output_fields = dependency_fields