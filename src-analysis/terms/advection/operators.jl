using Oceananigans: location
using Oceananigans.Fields: ZeroField

@inline swap_loc(loc::Face) = Center()
@inline swap_loc(loc::Center) = Face()
@inline swap_loc(loc::Nothing) = nothing

@inline swap_loc(loc::Type{Face}) = Center
@inline swap_loc(loc::Type{Center}) = Face
@inline swap_loc(loc::Type{Nothing}) = Nothing

for ξ in (:u, :v, :w), χ in (:U, :V, :W)
    advective_momentum_flux_density_χξ = Symbol(:advective_momentum_flux_density_, χ, ξ)
    χξ_func = Symbol(χ, ξ, :_func)
    @eval begin
        @inline function $χξ_func(i, j, k, grid, advection, U, u, Ub=ZeroField(eltype(grid)))
            U_tot = SumOfArrays{2}(U, Ub)
            $advective_momentum_flux_density_χξ(i, j, k, grid, advection, U_tot, u)
        end
    end
end

for (χ, l) in zip((:U, :V, :W), (:x, :y, :z))
    advective_tracer_flux_density_l = Symbol(:advective_tracer_flux_density_, l)
    χc_func = Symbol(χ, :c_func)
    @eval begin
        @inline function $χc_func(i, j, k, grid, advection, U, c, Ub=ZeroField(eltype(grid)))
            U_tot = SumOfArrays{2}(U, Ub)
            $advective_tracer_flux_density_l(i, j, k, grid, advection, U_tot, c)
        end
    end
end

for ξ in (:u, :v, :w, :c), χ in (:U, :V, :W)
    χξ_func = Symbol(χ, ξ, :_func)
    χξFlux = Symbol(χ, ξ, :Flux)
    @eval begin
        function $χξFlux(advection, U, u; background=ZeroField(eltype(grid)))
            grid = U.grid
            U_loc = location(U)
            u_loc = location(u)
            Uu_loc = [a isa Type{Face} ? swap_loc(b) : b for (a, b) in zip(U_loc, u_loc)]
            
            KernelFunctionOperation{Uu_loc...}($χξ_func, grid, advection, U, u, background)
        end
    end
end
