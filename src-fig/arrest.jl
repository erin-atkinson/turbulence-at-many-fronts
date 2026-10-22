@inline function dSdx_func(i, j, k, grid, mld)
    h = @inbounds mld[i, j, k]
    
    return sqrt(1 / h) / 4
end

function SCoordinate(mld)
    grid = mld.grid
    loc = location(mld)
    
    return CumulativeIntegral(KernelFunctionOperation{loc...}(dSdx_func, grid, mld); dims=1)
end

@inline function s_func(i, j, k, grid, S, x₁)
    S₁ = Oceananigans.interpolate((x₁, ), S, (Face(), nothing, nothing), grid)
    
    Sx = ℑxᶜᵃᵃ(i, j, k, grid, S)
    x = xnode(i, grid, Center())
    
    return (Sx - S₁)
end

function sCoordinate(S, x₁)
    grid = S.grid
    loc = location(S)
    
    return KernelFunctionOperation{loc...}(s_func, grid, S, x₁)
end

function isvalidbounds(data_points, bounds)
    x = data_points.x
    i_linear = bounds.i_linear
    i_ttw = bounds.i_ttw

    return i_linear >= i_ttw
end

function updatebounds(data_points, bounds, i_linear=5, i_ttw=length(data_points.x)-4)
    x = data_points.x
    b = data_points.b
    s = data_points.s
    h = data_points.h
    
    x_linear = x[i_linear]
    b_linear = b[i_linear]
    s_linear = s[i_linear]
    h_linear = h[i_linear]
    x_ttw = x[i_ttw]
    b_ttw = b[i_ttw]
    s_ttw = s[i_ttw]
    h_ttw = h[i_ttw]
    
    new_bounds = (; 
        i_linear, i_ttw, 
        x_linear, x_ttw,
        b_linear, b_ttw,
        s_linear, s_ttw,
        h_linear, h_ttw,
    )
    
    return merge(bounds, new_bounds)
end

"""
    front_fit_coordinates(b, h, sp)
Prepare data for fitting 
"""
function front_fit_coordinates(b, h, sp; front_region=(-sp.Lh/2, sp.Lh/2), cutoff=0.1)
    grid = h.grid

    # Edges of the frontal region
    b_left = Oceananigans.interpolate(front_region[1], b)
    b_right = Oceananigans.interpolate(front_region[2], b)
    Δb = b_right - b_left

    # Arrays
    x_data = xnodes(b)
    b_data = interior(b, :, 1, 1)[:]
    h_data = interior(h, :, 1, 1)[:]
    
    i_left = findfirst(x -> x >= front_region[1], x_data)
    i_right = findlast(x -> x <= front_region[2], x_data)
    
    # Position of minimum buoyancy
    i_min = argmin(b_data)
    x_min = x_data[i_min]
    b_min = b_data[i_min]
    h_min = h_data[i_min]
    
    # We cut out the bottom, top 10% for the fit
    b₀ = b_left + cutoff * (b_right - b_left)
    b₁ = b_left + (1 - cutoff) * (b_right - b_left)
    i₀ = findfirst(ib -> ib[2] >= b₀ && ib[1] >= i_left, collect(enumerate(b_data)))
    i₁ = findlast(ib -> ib[2] <= b₁ && ib[1] <= i_right, collect(enumerate(b_data)))

    x_data = x_data[i₀:i₁]
    b_data = b_data[i₀:i₁]
    h_data = h_data[i₀:i₁]
    
    x₀ = x_data[1]
    x₁ = x_data[end]
    b₀ = b_data[1]
    b₁ = b_data[end]
    h₀ = h_data[1]
    h₁ = h_data[end]
    
    S = Field(SCoordinate(h))
    fill_halo_regions!(S)
    
    s = Field(sCoordinate(S, x₁))
    fill_halo_regions!(s)

    s_data = interior(s, :, 1, 1)[:]
    s_min = s_data[i_min]
    
    s_data = s_data[i₀:i₁]
    s₀ = s_data[1]
    s₁ = s_data[end]
    
    data_points = (; x=x_data, s=s_data, b=b_data, h=h_data)
    bounds = (; x₀, x₁, s₀, s₁, b₀, b₁, h₀, h₁, x_min, s_min, b_min, h_min, Δb, b_left, b_right)
    bounds = updatebounds(data_points, bounds)
    
    return data_points, bounds
end

# -----------------------------------------------------------------------------------
function linear_weight_func(x, bounds) 
    x₀ = bounds.x₀
    x₁ = bounds.x₁
    x_linear = bounds.x_linear
    x_ttw = bounds.x_ttw
    
    x < x_linear && return one(x)
    x > x_ttw && return zero(x)
    return (x_ttw - x) / abs(x_ttw - x_linear)
end

function ttw_weight_func(x, bounds) 
    x₀ = bounds.x₀
    x₁ = bounds.x₁
    x_linear = bounds.x_linear
    x_ttw = bounds.x_ttw
    
    x > x_ttw && return one(x)
    x < x_linear && return zero(x)
    return (x - x_linear) / abs(x_ttw - x_linear)
end

x_fit_func(x, bounds, coeffs) = coeffs.Γa * (x - bounds.x_min) + coeffs.b_min
χ_fit_func(χ, bounds, coeffs) = coeffs.b_ref - coeffs.A²f² * χ^2
s_fit_func(s, bounds, coeffs) = χ_fit_func(s + coeffs.χ₁, bounds, coeffs)

function front_fit_func(x, s, bounds, coeffs)
    linear_weight = linear_weight_func(x, bounds)
    ttw_weight = ttw_weight_func(x, bounds)
    total_weight = linear_weight + ttw_weight
    
    linear_part = linear_weight * x_fit_func(x, bounds, coeffs)
    ttw_part = ttw_weight * s_fit_func(s, bounds, coeffs)
    return (linear_part + ttw_part) / total_weight
end

function front_fit_data(data_points, bounds, coeffs)
    x = data_points.x
    s = data_points.s
    b = data_points.b
    
    b_linear = similar(x)
    b_ttw = similar(s)
    b_front = similar(x)

    for i in eachindex(x)
        b_linear[i] = x_fit_func(x[i], bounds, coeffs)
        b_ttw[i] = s_fit_func(s[i], bounds, coeffs)
        b_front[i] = front_fit_func(x[i], s[i], bounds, coeffs)
    end
    
    return (; b_linear, b_ttw, b_front)
end

function front_fit_error(data_points, bounds, coeffs)
    x_data = data_points.x
    s_data = data_points.s
    b_data = data_points.b
    
    square_error = mapreduce(+, x_data[1:bounds.i_linear], b_data[1:bounds.i_linear]) do x, b
        (b - x_fit_func(x, bounds, coeffs))^2
    end
    
    square_error += mapreduce(+, s_data[bounds.i_ttw:end], b_data[bounds.i_ttw:end]) do s, b
        (b - s_fit_func(s, bounds, coeffs))^2
    end
    
    return square_error
end
# -----------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Functions to fit each part
"""
    linear_solution(data_points, bounds)
Try to fit a linear front up to `bounds.i_linear`
"""
function linear_solution(data_points, bounds)
    i = bounds.i_linear
    x = data_points.x[1:i]
    b = data_points.b[1:i]
    
    prob = CurveFitProblem(x, b)
    sol = solve(prob, LinearCurveFitAlgorithm())
    
    b_min = sol(bounds.x_min)
    Γa = sol.u[1]
    
    return (; b_min, Γa)
end

"""
    ttw_solution(data_points, bounds)
Try to fit a locally-ttw front down to `bounds.i_ttw`
"""
function ttw_solution(data_points, bounds)
    i = bounds.i_ttw
    s = data_points.s[i:end]
    b = data_points.b[i:end]
    
    prob = CurveFitProblem(s, b)
    sol = solve(prob, PolynomialFitAlgorithm(2))

    a₀ = sol.u[1]
    a₁ = sol.u[2]
    a₂ = sol.u[3]
    
    χ₁ = a₁ / 2a₂
    A²f² = -a₂
    b_ref = a₀ + A²f² * χ₁^2
    
    return (; χ₁, A²f², b_ref)
end

"""
    get_front_fit(data_points, bounds)
Try to fit the whole front with specific wall and linear region locations
Note that the data should be only in the frontal region.
"""
function get_front_fit(data_points, bounds)
    # Solve linear and TTW problem
    linear_coeffs = linear_solution(data_points, bounds)
    ttw_coeffs = ttw_solution(data_points, bounds)
    coeffs = merge(linear_coeffs, ttw_coeffs)
    
    return coeffs
end

function objective_func(data_points, bounds, coeffs, λ_linear, λ_transition, λ_ttw)
    b₀ = bounds.b₀
    b₁ = bounds.b₁
    x₀ = bounds.x₀
    x₁ = bounds.x₁
    x_linear = bounds.x_linear
    x_ttw = bounds.x_ttw
    
    mse_front = front_fit_error(data_points, bounds, coeffs) / (b₁ - b₀)^2
    
    linear_reg = λ_linear * (x_linear - x₀)^2 / (x₁ - x₀)^2
    transition_reg = λ_transition * (x_ttw - x_linear)^2 / (x₁ - x₀)^2
    ttw_reg = λ_ttw * (x₁ - x_ttw)^2 / (x₁ - x₀)^2
    
    return mse_front + transition_reg + linear_reg + ttw_reg
end

function maybegetobjective!(objective_cache, data_points, bounds, λ1, λ2)
    i_linear = bounds.i_linear
    i_ttw = bounds.i_ttw
    L = objective_cache[i_linear, i_ttw]
    !iszero(L) && return L
    
    coeffs = get_front_fit(data_points, bounds)
    L = objective_func(data_points, bounds, coeffs, λ1, λ2)
    objective_cache[i_linear, i_ttw] = L
    return L
end

@doc raw"""
    get_best_front_fit(data_points, bounds)
Try to fit a piecewise front. 

The objective function is
```math
\mathcal{L} = \sum \left (\frac{b_i - b(x_i)}{b_1 - b_0}\right )^2 + \lambda_1\left(\frac{x_\text{ttw} - x_\text{linear}}{x_1 - x_0}\right )^2 + \lambda_2\left(\frac{x_\text{linear} - x_0}{x_1 - x_0}\right )^2
```
i.e. we optimise the fit, and regularise via looking for the smallest linear region, and smallest difference between TTW and linear profile. 

Looking globally isn't that big of a deal
"""
function get_best_front_fit(data_points, bounds, λ_linear, λ_transition, λ_ttw)
    N = length(data_points.x)
    objective_cache = zeros(eltype(data_points.x), N, N)
    
    obj_min = Inf
    minimiser = [1, N]
    
    for i_linear in 5:(N-4), i_ttw in (i_linear+1):(N-4)
        #i_ttw = i_linear + 2
        bounds = updatebounds(data_points, bounds, i_linear, i_ttw)
        coeffs = get_front_fit(data_points, bounds)
        obj = objective_func(data_points, bounds, coeffs, λ_linear, λ_transition, λ_ttw)
    
        if obj < obj_min
            obj_min = obj
            minimiser[1] = i_linear
            minimiser[2] = i_ttw
        end
    end

    bounds = updatebounds(data_points, bounds, minimiser[1], minimiser[2])
    coeffs = get_front_fit(data_points, bounds)
    
    return bounds, coeffs
end
# -----------------------------------------------------------------------------------------------------------------------
