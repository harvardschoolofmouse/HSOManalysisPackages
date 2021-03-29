function render_trace(trace; show_data=true, limit_y=true, Alpha=1, ax=nothing, suppressAlpha=false)
    xs = get_args(trace)[1]
    xmin = minimum(xs)
    xmax = maximum(xs)
    if show_data
        ys = Vector{Float64}(undef, length(xs));
        for i=1:length(xs)
            yi = trace[:data => i => :y]
            ys[i] = yi;
        end
    end
    ymin = minimum(ys)
    ymax = maximum(ys)
    if isnothing(ax)
        ax = gca()
    end
    try
        slope = trace[:slope]
        intercept = trace[:intercept]
        noiseband = trace[:noise]
        noiseband_min = slope *  [xmin, xmax] .+ intercept .- noiseband
        noiseband_max = slope *  [xmin, xmax] .+ intercept .+ noiseband
        if suppressAlpha
            ax.fill_between([xmin, xmax],noiseband_min, noiseband_max, color="grey")
            plot([xmin, xmax], slope *  [xmin, xmax] .+ intercept, color="black")
        else
            ax.fill_between([xmin, xmax],noiseband_min, noiseband_max, color="grey",alpha=0.1)
            plot([xmin, xmax], slope *  [xmin, xmax] .+ intercept, color="black", alpha=0.5)
        end
    catch
        try
            intercept_flat = trace[:intercept_flat]
            noiseband = trace[:noise_flat]
            noiseband_min = intercept_flat .- noiseband
            noiseband_max = intercept_flat .+ noiseband
            if suppressAlpha
                ax.fill_between([xmin, xmax],noiseband_min, noiseband_max, color="grey")
                plot([xmin, xmax], [intercept_flat, intercept_flat], color="black")
            else
                ax.fill_between([xmin, xmax],noiseband_min, noiseband_max, color="grey",alpha=0.1)
                plot([xmin, xmax], [intercept_flat, intercept_flat], color="black", alpha=0.5)
            end
        catch
            left_segment_amp = trace[:left_segment_amp]
            right_segment_amp = trace[:right_segment_amp]
            step_position = trace[:step_position]
            noiseband = trace[:noise_step]
            ls = left_segment_amp.*ones(1,length(xs[1:step_position]))
            rs = right_segment_amp.*ones(1, length(xs[step_position+1:end]))
            mean_segment = vec(hcat(ls,rs))
            noiseband_min = mean_segment .- noiseband
            noiseband_max = mean_segment .+ noiseband
            if suppressAlpha
                ax.fill_between(xs, noiseband_min, noiseband_max, color="grey")
                plot(xs, mean_segment, color="black")
                plot(xs, noiseband_min, color="grey")
                plot(xs, noiseband_max, color="grey")
            else
                ax.fill_between(xs, noiseband_min, noiseband_max, color="grey",alpha=0.1)
                plot(xs, mean_segment, color="black", alpha=0.5)
                plot(xs, noiseband_min, color="grey", alpha=0.1)
                plot(xs, noiseband_max, color="grey", alpha=0.1)
            end
        end
    end
    if suppressAlpha
        scatter(xs, ys, c="blue")
    else
        scatter(xs, ys, c="blue", alpha=Alpha)
    end
    ax.set_xlim((xmin, xmax))
    if limit_y
        ax.set_ylim((ymin, ymax))
    end
end;

function grid(renderer::Function, traces; ncols=6, nrows=3)
    figure(figsize=(16, 8))
    for (i, trace) in enumerate(traces)
        subplot(nrows, ncols, i)
        renderer(trace, Alpha=1)
    end
end;
function overlay(renderer, traces; ax=nothing, suppressAlpha=false)
    for i=1:length(traces)
        renderer(traces[i], Alpha=1/length(traces), ax=ax,suppressAlpha=suppressAlpha)
    end
end;
function render_distribution(trace, data; yl="p", xl="data", ax=nothing, xlims=nothing)
    (xs,_)=Gen.get_args(trace)
    if isnothing(ax)
        figure(figsize=(2,2))
        ax = gca()
    end
    if isnothing(xlims)
        xlims = [0.7,maximum(xs)]
    end
        
    weights = ones(size(data))./length(data)
    ax.hist(data, weights=weights)#density=true)
    ax.set_ylabel(yl)
    ax.set_xlabel(xl)
    ax.set_xlim(xlims)
end;
function plot_summary(traces; param="slope", traceNames=[])
    n = length(traces)
    figure(figsize=(2,n))
    ax = []
    for i=1:n
        push!(ax, subplot(n,1,i))
        tr = traces[i]
        output_params = extract_mean_parameters(tr)
        if param == "slope"
            data = filter(!isnan,output_params.slopes)
        elseif param == "intercept"
            data = filter(!isnan,output_params.intercepts)
        elseif param == "step_position"
            data = output_params.step_positions
            (xs,_) = get_args(tr[1])
            for i = 1:length(tr)
                if !isnan(data[i])
                    data[i] = get_step_time_s(nothing; xs=xs,step_position=data[i])
                end
            end
            data = filter(!isnan,data)
        end
        if !isempty(data)
            xlims=[minimum([0.,minimum(data)]),maximum(data)]
            render_distribution(tr[1], data; yl="p", xl=param, ax=ax[i], xlims=xlims)
        else
            ax[i].set_xlim([0.,0.00001])
        end
        if !isempty(traceNames)
            ax[i].set_ylabel(traceNames[i])
        end
    end
    set_xaxes_same_scale(ax)
end
function render_proposal(xs, ys, proposed_choices; axs=nothing)
    slope = proposed_choices[:slope]
    intercept = proposed_choices[:intercept]
    noise = proposed_choices[:noise]
    noise_step = proposed_choices[:noise_step]
    left_segment_amp = proposed_choices[:left_segment_amp]
    right_segment_amp = proposed_choices[:right_segment_amp]
    step_position = round(Int,proposed_choices[:step_position])
    xmin = minimum(xs)
    xmax = maximum(xs)
    
    if isnothing(axs)
        axs=[]
        push!(axs,subplot(1,2,1))
    end
    noiseband = noise
    noiseband_min = slope *  [xmin, xmax] .+ intercept .- noiseband
    noiseband_max = slope *  [xmin, xmax] .+ intercept .+ noiseband
    axs[1].fill_between([xmin, xmax],noiseband_min, noiseband_max, color="grey",alpha=0.1)
    axs[1].plot([xmin, xmax], slope *  [xmin, xmax] .+ intercept, color="black", alpha=0.5)
    axs[1].plot(xs, ys)
    axs[1].set_title("ramp model")
    
    if length(axs)==1
        push!(axs,subplot(1,2,2))
    end
    noiseband = noise_step
    # Draw the line
    ls = left_segment_amp.*ones(1,length(xs[1:step_position]))
    rs = right_segment_amp.*ones(1, length(xs[step_position+1:end]))
    mean_segment = vec(hcat(ls,rs))
    noiseband_min = mean_segment .- noiseband
    noiseband_max = mean_segment .+ noiseband
    axs[2].fill_between(xs, noiseband_min, noiseband_max, color="grey",alpha=0.1)
    axs[2].plot(xs, mean_segment, color="black", alpha=0.5)
    axs[2].plot(xs, noiseband_min, color="grey", alpha=0.1)
    axs[2].plot(xs, noiseband_max, color="grey", alpha=0.1)
    axs[2].plot(xs, ys)
    axs[2].set_title("step model")
    return axs
end;