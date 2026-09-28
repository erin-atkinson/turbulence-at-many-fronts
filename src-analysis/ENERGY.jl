include("terms/terms.jl")

fields = (
    :u, :v, :w, :b,
    :uu, :uv, :uw, :vu, :vv, :vw, :wu, :wv, :ww,
    :ub, :vb, :wb,
    :u_prev, :v_prev, :w_prev, :b_prev,
    :ke
)

for ξ in fields
    @eval $ξ = input_fields.$ξ
end

velocities = (; u=u, v=v, w=w)
velocities_prev = (; u=u_prev, v=v_prev, w=w_prev)

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
dependency_fields = turbulent_fluxes

println("Kinetic energy terms...")

dsp_density = Field(DSPDensity(clock, velocities, velocities_prev, sp))
lsp_density = Field(LSPDensity(velocities, velocities_prev, turbulent_fluxes))
vsp_density = Field(VSPDensity(velocities, velocities_prev, turbulent_fluxes))
buoyancy_density = Field(BUOYANCYDensity(w, w_prev, b))
mke_density = Field(MKEDensity(velocities))
sponge_mke_density = Field(SPONGEMKEDensity(velocities, velocities_prev, sp))
strain_mke_density = Field(STRAINMKEDensity(clock, velocities, velocities_prev, sp))
stress = Field(STRESS(clock, velocities, velocities_prev, sp))
tke_density = Field(ke - mke_density)

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

b_profile = Field(Average(b; dims=1))
b_prev_profile = Field(Average(b_prev; dims=1))

h_ml = Field(BuoyancyBasedMixedLayer(b_profile, sp.Δb))
h_ml_prev = Field(BuoyancyBasedMixedLayer(b_prev_profile, sp.Δb))

mixed_layer = (; b_profile, b_prev_profile, h_ml, h_ml_prev)
dependency_fields = merge(dependency_fields, mixed_layer)

mpe_density = Field(MPEDensity(b, h_ml))
mixed_density = Field(MIXEDDensity(clock, b, b_prev, h_ml, h_ml_prev))
strain_mpe_density = Field(STRAINMPEDensity(clock, b, h_ml, h_ml_prev, sp))
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
