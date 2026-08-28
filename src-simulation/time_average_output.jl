using Oceananigans.Operators
using Oceananigans.Advection: advective_momentum_flux_Uu,
                              advective_momentum_flux_Vu,
                              advective_momentum_flux_Wu,
                              advective_momentum_flux_Uv,
                              advective_momentum_flux_Vv,
                              advective_momentum_flux_Wv,
                              advective_momentum_flux_Uw,
                              advective_momentum_flux_Vw,
                              advective_momentum_flux_Ww,
                              advective_tracer_flux_x,
                              advective_tracer_flux_y,
                              advective_tracer_flux_z

using Oceananigans.Utils: SumOfArrays
using Oceananigans.OutputWriters: AveragedSpecifiedTimes

@info "Constructing time averaged outputs..."

u, v, w = model.velocities
p = PressureField(model)
b = model.tracers.b

advection = model.advection
centered = Centered(; order=2)

# Across-front
function uu_func(i, j, k, grid, advection, centered, u, U)
    u_tot = SumOfArrays{2}(u, U)
    adv = advective_momentum_flux_Uu(i, j, k, grid, advection, u_tot, u)
    back = advective_momentum_flux_Uu(i, j, k, grid, centered, U, u)
    
    return (adv - back) / Axᶜᶜᶜ(i, j, k, grid)
end

function uv_func(i, j, k, grid, advection, centered, u, U, v)
    u_tot = SumOfArrays{2}(u, U)
    adv = advective_momentum_flux_Uv(i, j, k, grid, advection, u_tot, v)
    back = advective_momentum_flux_Uv(i, j, k, grid, centered, U, v)
    
    return  (adv - back) / Axᶠᶠᶜ(i, j, k, grid)
end

function uw_func(i, j, k, grid, advection, centered, u, U, w)
    u_tot = SumOfArrays{2}(u, U)
    adv = advective_momentum_flux_Uw(i, j, k, grid, advection, u_tot, w)
    back = advective_momentum_flux_Uw(i, j, k, grid, centered, U, w)
    
    return (adv - back) / Axᶠᶜᶠ(i, j, k, grid)
end

function ub_func(i, j, k, grid, advection, centered, u, U, b)
    u_tot = SumOfArrays{2}(u, U)
    adv = advective_tracer_flux_x(i, j, k, grid, advection, u_tot, b)
    back = advective_tracer_flux_x(i, j, k, grid, centered, U, b)
    
    return (adv - back) / Axᶠᶜᶜ(i, j, k, grid)
end

# Along-front
function vu_func(i, j, k, grid, advection, u, v)
    return advective_momentum_flux_Vu(i, j, k, grid, advection, v, u) / Ayᶠᶠᶜ(i, j, k, grid)
end

function vv_func(i, j, k, grid, advection, v)
    return advective_momentum_flux_Vv(i, j, k, grid, advection, v, v) / Ayᶜᶜᶜ(i, j, k, grid)
end

function vw_func(i, j, k, grid, advection, v, w)
    return advective_momentum_flux_Vw(i, j, k, grid, advection, v, w) / Ayᶜᶠᶠ(i, j, k, grid)
end

function vb_func(i, j, k, grid, advection, v, b)
    return advective_tracer_flux_y(i, j, k, grid, advection, v, b) / Ayᶜᶠᶜ(i, j, k, grid)
end

# Vertical
function wu_func(i, j, k, grid, advection, u, w)
    return advective_momentum_flux_Wu(i, j, k, grid, advection, w, u) / Azᶠᶜᶠ(i, j, k, grid)
end

function wv_func(i, j, k, grid, advection, v, w)
    return advective_momentum_flux_Wv(i, j, k, grid, advection, w, v) / Azᶜᶠᶠ(i, j, k, grid)
end

function ww_func(i, j, k, grid, advection, w)
    return advective_momentum_flux_Ww(i, j, k, grid, advection, w, w) / Azᶜᶜᶜ(i, j, k, grid)
end

function wb_func(i, j, k, grid, advection, w, b)
    return advective_tracer_flux_z(i, j, k, grid, advection, w, b) / Azᶜᶜᶠ(i, j, k, grid)
end

@inline u_sq(i, j, k, grid, u) = @inbounds u[i, j, k]^2
function ke_func(i, j, k, grid, u, v, w)
    u² = ℑxᶜᵃᵃ(i, j, k, grid, u_sq, u)
    v² = ℑyᵃᶜᵃ(i, j, k, grid, u_sq, v)
    w² = ℑzᵃᵃᶜ(i, j, k, grid, u_sq, w)

    return (u² + v² + w²) / 2
end

uu = KernelFunctionOperation{Center, Center, Center}(uu_func, grid, advection, centered, u, U)
uv = KernelFunctionOperation{Face, Face, Center}(uv_func, grid, advection, centered, u, U, v)
uw = KernelFunctionOperation{Face, Center, Face}(uw_func, grid, advection, centered, u, U, w)
ub = KernelFunctionOperation{Face, Center, Center}(ub_func, grid, advection, centered, u, U, b)

vu = KernelFunctionOperation{Face, Face, Center}(vu_func, grid, advection, u, v)
vv = KernelFunctionOperation{Center, Center, Center}(vv_func, grid, advection, v)
vw = KernelFunctionOperation{Center, Face, Face}(vw_func, grid, advection, v, w)
vb = KernelFunctionOperation{Center, Face, Center}(vb_func, grid, advection, v, b)

wu = KernelFunctionOperation{Face, Center, Face}(wu_func, grid, advection, u, w)
wv = KernelFunctionOperation{Center, Face, Face}(wv_func, grid, advection, v, w)
ww = KernelFunctionOperation{Center, Center, Center}(ww_func, grid, advection, w)
wb = KernelFunctionOperation{Center, Center, Face}(wb_func, grid, advection, w, b)

ke = KernelFunctionOperation{Center, Center, Center}(ke_func, grid, u, v, w)

output_fields = (; u, v, w, b, uu, uv, uw, ub, vu, vv, vw, vb, wu, wv, ww, wb, p, ke)

writing_times_pos = filter(t-> t > prev_time, 0:sp.save_time:sp.stop_time)
writing_times_neg = filter(t-> t > prev_time, (sp.start_time:sp.save_time:0)[1:end-1])
writing_times = [writing_times_neg; writing_times_pos]

output_symbol = Symbol(:fields, prev_iteration)
snapshot_symbol = Symbol(:snapshot, prev_iteration)
checkpointer_symbol = Symbol(:checkpointer, prev_iteration)

simulation.output_writers[output_symbol] = JLD2Writer(model, output_fields; 
    filename = "$output_folder/AVERAGE.jld2", 
    schedule = AveragedSpecifiedTimes(writing_times, sp.save_time),
    overwrite_existing = false,
    with_halos = true,
    init = init_jld2!
)

simulation.output_writers[snapshot_symbol] = JLD2Writer(model, (; u, v, w, b); 
    filename = "$output_folder/OUTPUT.jld2", 
    schedule = SpecifiedTimes(writing_times),
    overwrite_existing = false,
    with_halos = true,
    init = init_jld2!
)

simulation.output_writers[checkpointer_symbol] = Checkpointer(model;
    schedule = SpecifiedTimes(writing_times),
    dir = output_folder, 
    overwrite_existing = true,
    verbose = true,
    cleanup = true
)
