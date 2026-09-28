include("terms/terms.jl")

fields = (
    :u, :v, :w, :b, :p, 
    :uu, :uv, :uw, 
    :vu, :vv, :vw, 
    :wu, :wv, :ww, 
    :ub, :vb, :wb, 
    :ke
)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ=$ξ_bar)
    end
end

ψ = Field(Streamfunction(u_bar))
H = Field(BuoyancyBasedMixedLayer(b_bar, sp.Δb))

mean_fields = (; mean_fields..., ψ, H)

skip_update = filter(a->a ∉ fields, keys(input_fields))

dependency_fields = mean_fields
output_fields = dependency_fields
