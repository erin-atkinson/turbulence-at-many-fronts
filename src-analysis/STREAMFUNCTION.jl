include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/streamfunction/streamfunction.jl")
include("terms/streamfunction/sce.jl")
include("terms/streamfunction/ageostrophic.jl")
include("terms/streamfunction/strain.jl")
include("terms/streamfunction/turbulence.jl")
include("terms/streamfunction/sponge.jl")

fields = (:u, :v, :w, :b, :uu, :wu, :uw, :ww)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

# Streamfunction
ψ = Field(Streamfunction(u_bar))
streamfunction = (; ψ)

println("Turbulent fluxes")
u′u′ = Field(uu_bar - UuFlux(centered, u_bar, u_bar))
u′w′ = Field(uw_bar - UwFlux(centered, u_bar, w_bar))

w′u′ = Field(wu_bar - WuFlux(centered, w_bar, u_bar))
w′w′ = Field(ww_bar - WwFlux(centered, w_bar, w_bar))

turbulent_fluxes = (;
    uu = u′u′, uw = u′w′,
    wu = w′u′, ww = w′w′,
)

sce_density = Field(SCEDensity(ψ))
turbulence_density = Field(TURBULENCESCEDensity(ψ, turbulent_fluxes))
strain_density = Field(STRAINSCEDensity(clock, ψ, sp))
ageostrophic_density = Field(AGEOSTROPHICDensity(ψ, v_bar, b_bar, sp))
sponge_density = Field(SPONGESCEDensity(ψ, sp))

sce_density_terms = (; 
    sce_density,
    turbulence_density, 
    strain_density, 
    ageostrophic_density, 
    sponge_density
)

sce = Field(Integral(sce_density)) 
turbulence = Field(Integral(turbulence_density)) 
strain = Field(Integral(strain_density)) 
ageostrophic = Field(Integral(ageostrophic_density)) 
sponge = Field(Integral(sponge_density)) 

sce_terms = (;
    sce,
    turbulence,
    strain,
    ageostrophic,
    sponge
)

total = Field(turbulence + strain + ageostrophic + sponge)

skip_update = filter(a->a ∉ fields, keys(input_fields))
dependency_fields = merge(mean_fields, streamfunction, turbulent_fluxes, sce_density_terms, sce_terms)
output_fields = merge(streamfunction, sce_density_terms, sce_terms)
