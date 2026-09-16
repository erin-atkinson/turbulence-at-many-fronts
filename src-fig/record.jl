function prettyrecord(observable, fig, filename, frames::AbstractVector; record_kw...)
    N = length(frames)
    N == 1 && return prettyrecord(observable, fig, filename, frames[1]; record_kw...)
    filename = joinpath(videopath, filename * ".mp4")

    mkpath(splitdir(filename)[1])
    t0 = time()
    record(fig, filename, 1:N; record_kw...) do i
        observable[] = frames[i]
        
        t = time() - t0
        total_time = N * t / i
        
        str = @sprintf "%ds / %ds (%d / %d)" t total_time i N
        nl = i == N ? "\n" : "\r"
        
        print(str * nl)
    end
    return nothing
end

function prettyrecord(observable, fig, filename, frame::Number; record_kw...)
    isnothing(filename) && return nothing
    observable[] = frame
    filename = joinpath(videopath, filename * ".png")

    mkpath(splitdir(filename)[1])
    save(filename, fig; record_kw...)
    return nothing
end
