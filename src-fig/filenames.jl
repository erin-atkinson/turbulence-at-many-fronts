function make_filename(sp, ext=nothing, pre=""; βH=sp.βH, βα=sp.βα, βB=sp.βB, βτ=sp.βτ, θτ=sp.θτ)
    strs = map([βH, βα, βB, βτ, θτ]) do β
        replace(string(β), "."=>"_")
    end
    ext = ext == nothing ? "" : ".$ext"
    return joinpath(pre, join(strs, "-") * ext)
end

function θτ_from_str(ip)
    ip.βτ == 0 && return "C"
    ip.θτ == 0 && return "N"
    ip.θτ == π/2 && return "E"
end

default_ip = (;
    run_time = 8e5, start_time = -2e5, save_time = 1e4, max_time = -1.0,
    f = 1e-4, L = 1e3,
    βx = 20, βh = 3,
    Nx = 1400, Nh = 1024, Ny = 128, Nz = 64,
    βℓ = 6, βH = 0.1,
    βα = 0.1, βB = 0.03, βτ = 0.00, θτ = 0.0, β₀ = 6,
    comment = "Cooling, depth initialisation"
)

ips = [
    (; default_ip..., Nz=39, βH=0.06, βB=0.01, run_time=4.3e5, start_time=-4.3e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.03, run_time=3e5, start_time=-3e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.05, run_time=2.5e5, start_time=-2.5e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.1, run_time=2e5, start_time=-2e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.2, run_time=1.6e5, start_time=-1.6e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.01, run_time=6.0e5, start_time=-6.0e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.03, run_time=4.2e5, start_time=-4.2e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.05, run_time=3.5e5, start_time=-3.5e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.1, run_time=2.8e5, start_time=-2.8e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.2, run_time=2.2e5, start_time=-2.2e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.01, run_time=7.5e5, start_time=-7.5e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.03, run_time=5.2e5, start_time=-5.2e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.05, run_time=4.4e5, start_time=-4.4e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.1, run_time=3.5e5, start_time=-3.5e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.2, run_time=2.8e5, start_time=-2.8e5),
]

run_ids = [map(ip -> make_filename(ip; θτ=θτ_from_str(ip), βα), ips) for βα in (0.05, 0.1, 0.2)][2]
n_runs = length(run_ids)