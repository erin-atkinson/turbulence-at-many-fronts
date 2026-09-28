include("terms/terms.jl")

fields = (:u, :v, :w, :b, :uu, :wu, :uw, :ww)

for ξ in fields
    @eval $ξ = input_fields.$ξ
end

# Streamfunction
ψ = Field(Streamfunction(u))
streamfunction = (; ψ)

println("Turbulent fluxes")
u′u′ = Field(uu - UuFlux(centered, u, u))
u′w′ = Field(uw - UwFlux(centered, u, w))

w′u′ = Field(wu - WuFlux(centered, w, u))
w′w′ = Field(ww - WwFlux(centered, w, w))

turbulent_fluxes = (;
    uu = u′u′, uw = u′w′,
    wu = w′u′, ww = w′w′,
)

sce_density = Field(SCEDensity(ψ))
turbulence_density = Field(TURBULENCESCEDensity(ψ, turbulent_fluxes))
strain_density = Field(STRAINSCEDensity(clock, ψ, sp))
ageostrophic_density = Field(AGEOSTROPHICDensity(ψ, v, b, sp))
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
