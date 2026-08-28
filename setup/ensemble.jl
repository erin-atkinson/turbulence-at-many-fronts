# Default
default_ip = (;
    stop_time = 10e5, start_time = -2e5, save_time = 1e4, max_time = -1.0,
    f = 1e-4, L = 1e3,
    βx = 20, βh = 3,
    Nx = 1024, Nh = 768, Ny = 128, Nz = 64,
    βℓ = 6, βH = 0.1,
    βα = 0.1, βB = 0.03, βτ = 0.00, θτ = 0.0, β₀ = 6,
    comment = ""
)

# Some tests for central region size
size_test = (;
    ips = [
        (; default_ip..., βh=2, comment="Central region size test 2"),
        (; default_ip..., βh=3, comment="Central region size test 3"),
        (; default_ip..., βh=4, comment="Central region size test 4"),
        (; default_ip..., βh=5, comment="Central region size test 5"),
    ],
    filenames = [
        "size-test-h2",
        "size-test-h3",
        "size-test-h4",
        "size-test-h5",
    ]
)

# Some tests for domain resolution
resolution_test = (;
    ips = [
        (; default_ip..., βh = 4, Nx=1400, Nh=1024, Nz=96, comment="Resolution test x1400 z96"),
        (; default_ip..., βh = 4, Nx=1400, Nh=1024, Nz=64, comment="Resolution test x1400 z64"),
        (; default_ip..., βh = 4, Nx=1024, Nh=768, Nz=96, comment="Resolution test x1024 z96"),
        (; default_ip..., βh = 4, Nx=1024, Nh=768, Nz=32, comment="Resolution test x1024 z32"),
        (; default_ip..., βh = 4, Nx=800, Nh=600, Nz=64, comment="Resolution test x800 z64"),
        (; default_ip..., βh = 4, Nx=800, Nh=600, Nz=32, comment="Resolution test x800 z32"),
    ],
    filenames = [
        "resolution-test-x1400-z96",
        "resolution-test-x1400-z64",
        "resolution-test-x1024-z96",
        "resolution-test-x1024-z32",
        "resolution-test-x800-z64",
        "resolution-test-x800-z32",
    ]
)

# Some tests for outer region size
outer_test = (; 
    ips = [
        (; default_ip..., βh = 3, Nx=800, Nh=700, Nz=64, comment="Outer test x800"),
        (; default_ip..., βh = 3, Nx=900, Nh=700, Nz=64, comment="Outer test x900"),
        (; default_ip..., βh = 3, Nx=1000, Nh=700, Nz=64, comment="Outer test x1000"),
        (; default_ip..., βh = 3, Nx=1100, Nh=700, Nz=64, comment="Outer test x1100"),
    ],
    filenames = [
        "outer-test-x800",
        "outer-test-x900",
        "outer-test-x1000",
        "outer-test-x1100",
    ]
)

default_ip = (;
    stop_time = 10e5, start_time = -2e5, save_time = 1e4, max_time = -1.0,
    f = 1e-4, L = 1e3,
    βx = 20, βh = 3,
    Nx = 1400, Nh = 1024, Ny = 128, Nz = 64,
    βℓ = 6, βH = 0.1,
    βα = 0.1, βB = 0.03, βτ = 0.00, θτ = 0.0, β₀ = 60,
    comment = "Cooling, depth initialisation"
)

ips = [
    (; default_ip..., Nz=39, βH=0.06, βB=0.01, stop_time=0.0, start_time=-4.3e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.03, stop_time=0.0, start_time=-3e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.05, stop_time=0.0, start_time=-2.5e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.1, stop_time=0.0, start_time=-2e5),
    (; default_ip..., Nz=39, βH=0.06, βB=0.2, stop_time=0.0, start_time=-1.6e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.01, stop_time=0.0, start_time=-6.0e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.03, stop_time=0.0, start_time=-4.2e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.05, stop_time=0.0, start_time=-3.5e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.1, stop_time=0.0, start_time=-2.8e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.2, stop_time=0.0, start_time=-2.2e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.01, stop_time=0.0, start_time=-7.5e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.03, stop_time=0.0, start_time=-5.2e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.05, stop_time=0.0, start_time=-4.4e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.1, stop_time=0.0, start_time=-3.5e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.2, stop_time=0.0, start_time=-2.8e5),
]

filenames = map(ips) do ip
    make_filename(ip; θτ=θτ_from_str(ip), βα="init")
end

destfilenamess = map(ips) do ip
    [
        make_filename(ip; θτ=θτ_from_str(ip), βα=0.1), 
        make_filename(ip; θτ=θτ_from_str(ip), βα=0.2), 
        make_filename(ip; θτ=θτ_from_str(ip), βα=0.05)
    ]
end

cooling_depth_init = (;
    ips,
    filenames,
    destfilenamess
)

ips_01 = map(ips) do ip
    (; ip..., stop_time=10e5, βα=0.1)
end

ips_02 = map(ips[7:9]) do ip
    (; ip..., stop_time=5e5, βα=0.2)
end

ips_005 = map(ips[7:9]) do ip
    (; ip..., stop_time=20e5, βα=0.05)
end

cooling_depth_01 = (;
    ips = ips_01,
    filenames = map(ip -> make_filename(ip; θτ=θτ_from_str(ip)), ips_01)
)

cooling_depth_02 = (;
    ips = ips_02,
    filenames = map(ip -> make_filename(ip; θτ=θτ_from_str(ip)), ips_02)
)

cooling_depth_005 = (;
    ips = ips_005,
    filenames = map(ip -> make_filename(ip; θτ=θτ_from_str(ip)), ips_005)
)

test_set = (;
    ips = cooling_depth_01.ips[7:9],
    filenames = cooling_depth_01.filenames[7:9],
)

test_set_init = (;
    ips = cooling_depth_init.ips[7:9],
    filenames = cooling_depth_init.filenames[7:9],
)

#= I probably won't use these
# Cooling only
ips = [
    (; default_ip..., Nz=64, βH=0.1, βB=0.01, stop_time=0.0, start_time=-6.0e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.02, stop_time=0.0, start_time=-4.8e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.03, stop_time=0.0, start_time=-4.2e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.04, stop_time=0.0, start_time=-3.8e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.05, stop_time=0.0, start_time=-3.5e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.06, stop_time=0.0, start_time=-3.3e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.07, stop_time=0.0, start_time=-3.1e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.08, stop_time=0.0, start_time=-3.0e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.09, stop_time=0.0, start_time=-2.9e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.1, stop_time=0.0, start_time=-2.8e5),
    (; default_ip..., Nz=64, βH=0.1, βB=0.11, stop_time=0.0, start_time=-2.7e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.12, stop_time=0.0, start_time=-2.6e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.13, stop_time=0.0, start_time=-2.6e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.14, stop_time=0.0, start_time=-2.5e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.15, stop_time=0.0, start_time=-2.4e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.2, stop_time=0.0, start_time=-2.2e5),
]
ips = map(ips) do ip
    (; ip..., stop_time=10e5, βα=0.1)
end

filenames = map(ips) do ip
    make_filename(ip; θτ=θτ_from_str(ip))
end

cooling_only = (;
    ips,
    filenames,
)

# Depth only
ips = [
    (; default_ip..., Nz=39, βH=0.06, βB=0.03, stop_time=0.0, start_time=-3e5),
    (; default_ip..., Nz=51, βH=0.08, βB=0.03, stop_time=0.0, start_time=-3.6e5), 
    (; default_ip..., Nz=58, βH=0.09, βB=0.03, stop_time=0.0, start_time=-3.9e5), 
    (; default_ip..., Nz=64, βH=0.1, βB=0.03, stop_time=0.0, start_time=-4.2e5),
    (; default_ip..., Nz=77, βH=0.12, βB=0.03, stop_time=0.0, start_time=-4.7e5), 
    (; default_ip..., Nz=83, βH=0.13, βB=0.03, stop_time=0.0, start_time=-5.0e5),
    (; default_ip..., Nz=90, βH=0.14, βB=0.03, stop_time=0.0, start_time=-5.2e5)
]

ips = map(ips) do ip
    (; ip..., stop_time=10e5, βα=0.1)
end

filenames = map(ips) do ip
    make_filename(ip; θτ=θτ_from_str(ip))
end

depth_only = (;
    ips,
    filenames,
)
=#
