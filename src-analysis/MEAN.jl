# Get mean state using an along-front average
include("terms/terms.jl")

fields = (:u, :v, :w, :b, :p, :uu, :uv, :uw, :vu, :vv, :vw, :wu, :wv, :ww)

for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (mean_fields..., $ξ_bar)
    end
end

ψ = Field(CumulativeIntegral(-u_bar; dims=3))

skip_update = filter(a->a ∉ fields, keys(input_fields))

dependency_fields = mean_fields
output_fields = dependency_fields
