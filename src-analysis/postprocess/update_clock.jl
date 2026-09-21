function update_clock!(clock, iterations, times, frame)
    clock.time = times[frame]
    clock.iteration = iterations[frame]
    clock.last_Δt = frame > 1 ? times[frame] - times[frame-1] : Inf
    return nothing
end
