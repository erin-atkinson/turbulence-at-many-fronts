# balances.jl
# The balance equations have quite a lot of repetition so it makes sense to combine them

import Base.getindex

u_balance_terms = (
    :advection_x, :advection_background, :advection_z,
    :mixing_x, :mixing_z,
    :coriolis, :pressure, :strain, :sponge, :surface
)

u_balance_term_labels = (;
    advection_x = L"-\overline{u}\frac{\partial \overline{u}}{\partial x}",
    advection_background = L"-U\frac{\partial \overline{u}}{\partial x}",
    advection_z = L"-w\frac{\partial \overline{u}}{\partial z}",
    mixing_x = L"-\frac{\partial }{\partial x}\overline{u'u'}",
    mixing_z = L"-\frac{\partial }{\partial z}\overline{w'u'}",
    coriolis = L"f\overline{v}",
    pressure = L"\frac{\partial \overline{p}}{\partial x}",
    strain = L"-\frac{\partial U}{\partial x}\overline{u}",
    sponge = L"-S_u",
    surface = L"\tau_x\delta (z)"
)

v_balance_terms = (
    :advection_x, :advection_background, :advection_z,
    :mixing_x, :mixing_z,
    :coriolis, :strain, :sponge, :surface
)

v_balance_term_labels = (;
    advection_x = L"-\overline{u}\frac{\partial \overline{v}}{\partial x}",
    advection_background = L"-U\frac{\partial \overline{v}}{\partial x}",
    advection_z = L"-w\frac{\partial \overline{v}}{\partial z}",
    mixing_x = L"-\frac{\partial }{\partial x}\overline{u'v'}",
    mixing_z = L"-\frac{\partial }{\partial z}\overline{w'v'}",
    coriolis = L"-f\overline{u}",
    strain = L"\frac{\partial U}{\partial x}\overline{v}",
    sponge = L"-S_v",
    surface = L"\tau_y\delta (z)"
)

b_balance_terms = (
    :advection_x, :advection_background, :advection_z,
    :mixing_x, :mixing_z,
    :surface
)

b_balance_term_labels = (;
    advection_x = L"-\overline{u}\frac{\partial \overline{v}}{\partial x}",
    advection_background = L"-U\frac{\partial \overline{v}}{\partial x}",
    advection_z = L"-w\frac{\partial \overline{v}}{\partial z}",
    mixing_x = L"-\frac{\partial }{\partial x}\overline{u'v'}",
    mixing_z = L"-\frac{\partial }{\partial z}\overline{w'v'}",
    surface = L"-B\delta (z)"
)

MKE_terms = (
    :dsp,
    :lsp,
    :vsp,
    :buoyancy,
    :sponge_mke,
    :strain_mke,
    :wind
)

MKE_term_labels = (;
    dsp = L"\text{DSP}",
    lsp = L"-\text{LSP}",
    vsp = L"-\text{VSP}",
    buoyancy = L"\text{BUOYANCY}",
    sponge_mke = L"\text{SPONGE}_\text{MKE}",
    strain_mke = L"\text{STRAIN}_\text{MKE}",
    wind = L"\text{WIND}"
)

MKE_term_signs = (;
    dsp = 1,
    lsp = -1,
    vsp = -1,
    buoyancy = 1,
    sponge_mke = 1,
    strain_mke = 1,
    wind = 1,
)

MPE_terms = (;
    :buoyancy,
    :bflux,
    :cooling,
    :strain_mpe,
    :mixed
)

MPE_term_labels = (;
    buoyancy = L"-\text{BUOYANCY}",
    bflux = L"-\text{BFLUX}",
    cooling = L"\text{COOLING}",
    strain_mpe = L"\text{STRAIN}_\text{MPE}",
    mixed = L"\text{MIXED}",
)

MPE_term_signs = (;
    buoyancy = -1,
    bflux = -1,
    cooling = 1,
    strain_mpe = 1,
    mixed = 1
)

abstract type AbstractBalance end

function scriptname(::AbstractBalance) end
function terms(::AbstractBalance) end
function termlabels(::AbstractBalance) end
function targetterm(::AbstractBalance) end
function totalterm(::AbstractBalance) end
permittedterms(balance::AbstractBalance) = Tuple(terms..., targetterm(balance), totalterm(balance))

windowlength(balance::AbstractBalance) = balance.N
filepath(balance::AbstractBalance, altname=nothing) = joinpath(scratchpath, balance.run_id, filename(balance, altname))

function filename(balance::AbstractBalance, altname=nothing)
    name = isnothing(altname) ? scriptname(balance) : altname
    N = windowlength(balance)
    N == 1 && return name * ".jld2"
    return name * "-$N.jld2"
end

iterations_times(balance::AbstractBalance) = iterations_times(filepath(balance))

function Base.getindex(balance::AbstractBalance, term::Symbol)
    term ∉ permittedterms(balance) && throw(FieldError(balance, term))
    termsign = term ∈ terms(balance) ? termsigns(balance)[term] : 1

    iterations, _ = iterations_times(balance)
    return timeseries_of(identity, filepath(balance), string(term), iterations) .* termsign
end

termsigns(balance::AbstractBalance) = NamedTuple(k=>1 for k in terms(balance))

@doc raw"""
    target(balance::AbstractBalance)
Return a timeseries of the target tendency of a balance equation
"""
function target(balance::AbstractBalance)
    _, times = iterations_times(balance)
    timeseries = balance[targetterm(balance)]
    result = similar(timeseries)

    for i in axes(timeseries, 1)
        result[i] = (timeseries[i] - timeseries[max(i, 1)]) / (times[i] - times[max(i, 1)])
    end

    return result
end

total(balance::AbstractBalance) = balance[totalterm(balance)]

@doc raw"""
    alltimeseries(balance::AbstractBalance)
Construct a named tuple of timeseries that sum to target(balance)
"""
function alltimeseries(balance::AbstractBalance)
    NamedTuple(term => balance.term for term in terms(balance))
end

struct MKEBalance <: AbstractBalance
    run_id
    N
end

scriptname(balance::MKEBalance) = "ENERGY"
terms(::MKEBalance) = MKE_terms
termlabels(::MKEBalance) = MKE_term_labels
termsigns(::MKEBalance) = MKE_term_signs
targetterm(::MKEBalance) = :mke
totalterm(::MKEBalance) = :mke_total

struct MPEBalance <: AbstractBalance
    run_id
    N
end

scriptname(balance::MPEBalance) = "ENERGY"
terms(::MKEBalance) = MPE_terms
termlabels(::MPEBalance) = MPE_term_labels
termsigns(::MPEBalance) = MPE_term_signs
targetterm(::MPEBalance) = :mpe
totalterm(::MPEBalance) = :mpe_total

abstract type QuadraticBalance <: AbstractBalance end

function Base.getindex(balance::QuadraticBalance, term::Symbol)
    term ∉ permittedterms(balance) && throw(FieldError(balance, term))
    termstr = term ∈ terms(balance) ? "quadratic_" * string(term) : string(term)

    iterations, _ = iterations_times(balance)
    return timeseries_of(identity, filepath(balance), termstr, iterations)
end

totalterm(::QuadraticBalance) = :quadratic_total

struct UBalance <: QuadraticBalance
    run_id
    N
end

scriptname(balance::UBalance) = "UBALANCE"
terms(::UBalance) = u_balance_terms
termlabels(::UBalance) = u_balance_term_labels
targetterm(::UBalance) = :u²

struct VBalance <: QuadraticBalance
    run_id
    N
end

scriptname(balance::VBalance) = "VBALANCE"
terms(::VBalance) = v_balance_terms
termlabels(::VBalance) = v_balance_term_labels
targetterm(::VBalance) = :v²

struct BBalance <: QuadraticBalance
    run_id
    N
end

scriptname(balance::BBalance) = "BBALANCE"
terms(::BBalance) = b_balance_terms
termlabels(::BBalance) = b_balance_term_labels
targetterm(::BBalance) = :b²

@doc raw"""
    plot_balance!(ax, times, balance::AbstractBalance; kwargs...)
Return a named tuple of lines produced by plotting lines! for each term in balance

See also [`plot_target!`](@ref), [`plot_total!`](@ref)
"""
function plot_balance!(ax, times, balance::AbstractBalance, unit=1; kwargs...)
    lns = NamedTuple(term => lines!(ax, times, balance[term] ./ unit; kwargs...) for term in terms(balance))
    return lns
end
plot_target!(ax, times, balance::AbstractBalance, unit=1.0; color=:black, linestyle=:dash, kwargs...) = lines!(ax, times, target(balance) ./ unit; color, linestyle, kwargs...)
plot_total!(ax, times, balance::AbstractBalance, unit=1.0; color=:black, kwargs...) = lines!(ax, times, total(balance) ./ unit; color, linestyle, kwargs...)

make_legend!(gl, lns, balance; kwargs...) = Legend(gl, [lns...], termlabels(balance); kwargs...)
