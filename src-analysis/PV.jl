include("terms/terms.jl")
include("terms/vorticity/vorticity.jl")

include("terms/pv/q.jl")
include("terms/pv/flux.jl")
include("terms/pv/turbulence.jl")
include("terms/pv/background.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

fields = (
    :u, :v, :w, :b,
    :uu, :uv, :uw, 
    :vu, :vv, :vw, 
    :wu, :wv, :ww,
    :ub, :vb, :wb,
)

mean_fields = NamedTuple()
model_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end
dependency_fields = mean_fields

println("Turbulent fluxes")
u′u′ = Field(uu_bar - UuFlux(centered, u_bar, u_bar))
u′v′ = Field(uv_bar - UvFlux(centered, u_bar, v_bar))
u′w′ = Field(uw_bar - UwFlux(centered, u_bar, w_bar))
u′b′ = Field(ub_bar - UcFlux(centered, u_bar, b_bar))

v′u′ = Field(vu_bar - VuFlux(centered, v_bar, u_bar))
v′v′ = Field(vv_bar - VvFlux(centered, v_bar, v_bar))
v′w′ = Field(vw_bar - VwFlux(centered, v_bar, w_bar))
v′b′ = Field(vb_bar - VcFlux(centered, v_bar, b_bar))

w′u′ = Field(wu_bar - WuFlux(centered, w_bar, u_bar))
w′v′ = Field(wv_bar - WvFlux(centered, w_bar, v_bar))
w′w′ = Field(ww_bar - WwFlux(centered, w_bar, w_bar))
w′b′ = Field(wb_bar - WcFlux(centered, w_bar, b_bar))

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
q_mean_fields = (; u=u_bar, v=v_bar, w=w_bar, b=b_bar)
q_turbulent_forcing = (; u=𝔉u, v=𝔉v, w=𝔉w, b=𝔇)

q = Field(PV(u_bar, v_bar, w_bar, b_bar, sp))

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
