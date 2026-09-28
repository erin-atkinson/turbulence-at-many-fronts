include("terms/terms.jl")

fields = (
    :u, :v, :w, :b,
    :uu, :uv, :uw, 
    :vu, :vv, :vw, 
    :wu, :wv, :ww,
    :ub, :vb, :wb,
)

for ξ in fields
    @eval $ξ = input_fields.$ξ
end

println("Turbulent fluxes")
u′u′ = Field(uu - UuFlux(centered, u, u))
u′v′ = Field(uv - UvFlux(centered, u, v))
u′w′ = Field(uw - UwFlux(centered, u, w))
u′b′ = Field(ub - UcFlux(centered, u, b))

v′u′ = Field(vu - VuFlux(centered, v, u))
v′v′ = Field(vv - VvFlux(centered, v, v))
v′w′ = Field(vw - VwFlux(centered, v, w))
v′b′ = Field(vb - VcFlux(centered, v, b))

w′u′ = Field(wu - WuFlux(centered, w, u))
w′v′ = Field(wv - WvFlux(centered, w, v))
w′w′ = Field(ww - WwFlux(centered, w, w))
w′b′ = Field(wb - WcFlux(centered, w, b))

turbulent_fluxes = (;
    uu = u′u′, uv = u′v′, uw = u′w′, ub = u′b′,
    vu = v′u′, vv = v′v′, vw = v′w′, vb = v′b′,
    wu = w′u′, wv = w′v′, ww = w′w′, wb = w′b′,
)
dependency_fields = merge(dependency_fields, turbulent_fluxes)

println("Turbulent forcing")
𝔉u = Field(-(∂x(u′u′) + ∂z(w′u′))) # Field(-(∂x(u′u′) + ∂y(v′u′) + ∂z(w′u′)))
𝔉v = Field(-(∂x(u′v′) + ∂z(w′v′))) # Field(-(∂x(u′v′) + ∂y(v′v′) + ∂z(w′v′)))
𝔉w = Field(-(∂x(u′w′) + ∂z(w′w′))) # Field(-(∂x(u′w′) + ∂y(v′w′) + ∂z(w′w′)))
𝔇 = Field(-(∂x(u′b′) + ∂z(w′b′))) # Field(-(∂x(u′b′) + ∂y(v′b′) + ∂z(w′b′)))

turbulent_forcing = (; 𝔉u, 𝔉v, 𝔉w, 𝔇)
dependency_fields = merge(dependency_fields, turbulent_forcing)

println("Potential vorticity")
q_mean_fields = (; u=u, v=v, w=w, b=b)
q_turbulent_forcing = (; u=𝔉u, v=𝔉v, w=𝔉w, b=𝔇)

q = Field(PV(u, v, w, b, sp))

Jx = Field(PVFluxX(q_mean_fields, sp))
Jy = Field(PVFluxY(q_mean_fields, sp))
Jz = Field(PVFluxZ(q_mean_fields, sp))

# Since we want a gauge where total Jy = 0
Jb_inf = Field(PVFluxBackground(q_mean_fields, input_fields.U, sp))
Jb = Field(RefToZero(Jb_inf))

𝔍x = Field(TurbulentPVFluxX(q_mean_fields, q_turbulent_forcing, sp))
𝔍y = Field(TurbulentPVFluxY(q_mean_fields, q_turbulent_forcing, sp))
𝔍z = Field(TurbulentPVFluxZ(q_mean_fields, q_turbulent_forcing, sp))

potential_vorticity = (; q, Jx, Jy, Jz, Jb_inf, Jb, 𝔍x, 𝔍y, 𝔍z)
dependency_fields = merge(dependency_fields, potential_vorticity)

skip_update = filter(a->a ∉ fields, keys(input_fields))
output_fields = (; q, Jx, Jy, Jz, Jb, 𝔍x, 𝔍y, 𝔍z)
