include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

fields = (
    :u, :v, :w, :b,
    :uu, :uv, :uw, :vu, :vv, :vw, :wu, :wv, :ww,
    :ub, :vb, :wb,
    :u_prev, :v_prev, :w_prev, :b_prev,
    :ke
)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

velocities = (; u=u_bar, v=v_bar, w=w_bar)
velocities_prev = (; u=u_prev_bar, v=v_prev_bar, w=w_prev_bar)

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

println("Kinetic energy terms...")
include("terms/energy/dsp_density.jl")
include("terms/energy/lsp_density.jl")
include("terms/energy/vsp_density.jl")
include("terms/energy/buoyancy_density.jl")
include("terms/energy/mke_density.jl")
include("terms/energy/sponge_mke_density.jl")
include("terms/energy/strain_mke_density.jl")
include("terms/energy/stress.jl")

dsp_density = Field(DSPDensity(clock, velocities, velocities_prev, sp))
lsp_density = Field(LSPDensity(velocities, velocities_prev, turbulent_fluxes))
vsp_density = Field(VSPDensity(velocities, velocities_prev, turbulent_fluxes))
buoyancy_density = Field(BUOYANCYDensity(w_bar, w_prev_bar, b_bar))
mke_density = Field(MKEDensity(velocities))
sponge_mke_density = Field(SPONGEMKEDensity(velocities, velocities_prev, sp))
strain_mke_density = Field(STRAINMKEDensity(clock, velocities, velocities_prev, sp))
stress = Field(STRESS(clock, velocities, velocities_prev, sp))
tke_density = Field(ke_bar - mke_density)

mke_production_density = (;
    dsp_density,
    lsp_density,
    vsp_density,
    buoyancy_density,
    mke_density,
    sponge_mke_density,
    strain_mke_density,
    stress,
    tke_density
)
dependency_fields = merge(dependency_fields, mke_production_density)
@info mke_production_density

dsp = Field(Integral(dsp_density))
lsp = Field(Integral(lsp_density))
vsp = Field(Integral(vsp_density))
buoyancy = Field(Integral(buoyancy_density))
mke = Field(Integral(mke_density))
sponge_mke = Field(Integral(sponge_mke_density))
strain_mke = Field(Integral(strain_mke_density))
wind = Field(Integral(stress))
tke = Field(Integral(tke_density))

mke_production = (;
    dsp,
    lsp,
    vsp,
    buoyancy,
    mke,
    sponge_mke,
    strain_mke,
    wind,
    tke
)
dependency_fields = merge(dependency_fields, mke_production)
@info mke_production

println("Potential energy terms...")
include("terms/energy/mld.jl")
include("terms/energy/mpe_density.jl")
include("terms/energy/mixed_density.jl")
include("terms/energy/strain_mpe_density.jl")
include("terms/energy/cooling.jl")
include("terms/energy/bflux_density.jl")

b_profile = Field(Average(b_bar; dims=1))
b_prev_profile = Field(Average(b_bar; dims=1))

h_ml = Field(MLD(b_profile, 2sp.Δb / sp.H))
h_ml_prev = Field(MLD(b_prev_profile, 2sp.Δb / sp.H))

mixed_layer = (; b_profile, b_prev_profile, h_ml, h_ml_prev)
dependency_fields = merge(dependency_fields, mixed_layer)

mpe_density = Field(MPEDensity(b_bar, h_ml))
mixed_density = Field(MIXEDDensity(clock, b_bar, b_prev_bar, h_ml, h_ml_prev))
strain_mpe_density = Field(STRAINMPEDensity(clock, b_bar, h_ml, h_ml_prev, sp))
bflux_density = Field(BFLUXDensity(w′b′))

mpe_production_density = (;
    mpe_density,
    mixed_density,
    strain_mpe_density,
    bflux_density
)
dependency_fields = merge(dependency_fields, mpe_production_density)
@info mpe_production_density

mpe = Field(Integral(mpe_density))
mixed = Field(Integral(mixed_density))
strain_mpe = Field(Integral(strain_mpe_density))
bflux = Field(Integral(bflux_density))
cooling = Field(COOLING(clock, h_ml, h_ml_prev, sp))

mpe_production = (;
    h_ml,
    mpe,
    mixed,
    strain_mpe,
    bflux,
    cooling,
)
dependency_fields = merge(dependency_fields, mpe_production)
@info mpe_production

mke_total = Field(dsp + wind + buoyancy + sponge_mke - lsp - vsp + strain_mke)
mpe_total = Field(-buoyancy - bflux + mixed + cooling + strain_mpe)

totals = (; mke_total, mpe_total)
dependency_fields = merge(dependency_fields, totals)

output_fields = merge(
    mke_production_density, 
    mke_production, 
    mpe_production_density, 
    mpe_production,
    totals
)
skip_update = (:p, )
