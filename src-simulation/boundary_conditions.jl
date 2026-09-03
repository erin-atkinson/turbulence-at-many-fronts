# boundary_conditions.jl

# ---------------------------------------
# Cooling turns on slowly
@inline function b_flux_func(x, y, t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return sp.B * turnon
end

# θ: angle relative to a down-front wind
# We only include wind in the central region
@inline function u_flux_func(x, y, t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return -sp.τ * turnon * sin(sp.θτ) * exp(-x^2 / 4sp.L^2)
end

@inline function v_flux_func(x, y, t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return -sp.τ * turnon * cos(sp.θτ) * exp(-x^2 / 4sp.L^2)
end
# ---------------------------------------

# ---------------------------------------
# Boundary conditions have no flow outside of the domain
b_bcs = FieldBoundaryConditions(;
    bottom=GradientBoundaryCondition(sp.N₀²),
    top=FluxBoundaryCondition(b_flux_func; parameters=(; sp.B, sp.f, sp.start_time))
)
u_bcs=FieldBoundaryConditions(;
    top=FluxBoundaryCondition(u_flux_func; parameters=(; sp.θτ, sp.f, sp.τ, sp.start_time, sp.L))
)
v_bcs=FieldBoundaryConditions(;
    east=ValueBoundaryCondition(0),
    west=ValueBoundaryCondition(0),
    top=FluxBoundaryCondition(v_flux_func; parameters=(; sp.θτ, sp.f, sp.τ, sp.start_time, sp.L))
)
w_bcs=FieldBoundaryConditions(;
    east=ValueBoundaryCondition(0),
    west=ValueBoundaryCondition(0)
)
# ---------------------------------------


boundary_conditions = (; u=u_bcs, v=v_bcs, b=b_bcs, w=w_bcs)
@info "Created boundary conditions"