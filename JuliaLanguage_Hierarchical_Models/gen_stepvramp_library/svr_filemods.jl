#
# Includes tools for autoloading dependencies
#
struct Result
   traces
   flags
end
struct MR3
    sessionCode::String
    modelID::String
    amount_of_computation::Int
    ntraces_per_trial::Int
    p::Vector{Float64}
    traces::Array{Vector,1}
    slopes::Array{Vector{Float64},1}
    intercepts::Array{Vector{Float64},1}
    l::Array{Vector{Float64},1}
    r::Array{Vector{Float64},1}
    st::Array{Vector{Float64},1}
    flags::Array{Vector{Bool},1}
end
struct MR4
    sessionCode::String
    modelID::String
    amount_of_computation::Int
    ntraces_per_trial::Int
    p::Vector{Float64}
    traces::Array{Vector,1}
    slopes::Array{Vector{Float64},1}
    intercepts::Array{Vector{Float64},1}
    l::Array{Vector{Float64},1}
    r::Array{Vector{Float64},1}
    st::Array{Vector{Float64},1}
    flags::Array{Vector{Bool},1}
    tNo::Vector{Int}
    lick_time::Vector{Float64}
end

function saveMR4asMAT(data::MR4, path, trim_cue_s, trim_lick_s)
	headsup("refreshed")
	retdir = pwd()
	cd(path)
	filename = join([data.sessionCode[1], "_amtcomp", data.amount_of_computation[1], "_ntr", data.ntraces_per_trial[1], timestamp_now(), ".mat"])
	file = matopen(filename, "w")
	write(file, "sessionCode", data.sessionCode)
	write(file, "modelID", data.modelID)
	write(file, "amount_of_computation", data.amount_of_computation)
	write(file, "ntraces_per_trial", data.ntraces_per_trial)
	write(file, "p", data.p)
	write(file, "slopes", data.slopes)
	write(file, "intercepts", data.intercepts)
	write(file, "left_segment", data.l)
	write(file, "right_segment", data.r)
	write(file, "step_time", data.st)
	write(file, "flags", data.flags)
	write(file, "tNo", data.tNo)
	write(file, "lick_time", data.lick_time)
	write(file, "trim_cue_s", trim_cue_s)
	write(file, "trim_lick_s", trim_lick_s)
	close(file)
	cd(retdir)
end


function logmeanexp(scores::Vector{Float64})
    logsumexp(scores) - log(length(scores))
end;
function make_constraints(ys::Vector{Float64})
    constraints = Gen.choicemap()
    for i=1:length(ys)
        constraints[:data => i => :y] = ys[i]
    end
    constraints
end;
function run_generic_inference(inference_fxn, generative_fxn, genfxn_args, rendring_fxn, xs::Vector{Float64}, ys::Vector{Float64}; n=12, Iters=500, long=false)
    #genfxn_args = (xs, noiseCenter)
    observations = make_constraints(ys);
    scores = Vector{Float64}(undef, n)
    traces = [];
    for i=1:n
        if ~long
            @time tr = inference_fxn(generative_fxn, genfxn_args, xs, ys, iters=Iters)
        else
            progressbar(i,n)
            tr = inference_fxn(generative_fxn, genfxn_args, xs, ys, iters=Iters)
        end
        scores[i] = get_score(tr)
        push!(traces, tr)

    end
    println("Mean model score: ", logmeanexp(scores))
    get_mean_slope(traces, verbose=true)
    get_mean_intercept(traces, verbose=true)
    get_break_points(traces, verbose=true)
    traces
end;
function get_rand_trial_data(data; idx=nothing)
    if isnothing(idx)
        idx = rand(1:length(data.xdata))
    end
    xs = data.xdata[idx]
    ys = data.ydata[idx]
    tNo = data.trialNo[idx]
    lt = data.lickTime_s[idx]
    path = data.path
    seshcode = data.sessionCode;
    trial_data = TrialData([xs],[ys],[tNo],[lt],path,seshcode)
    return trial_data
end;