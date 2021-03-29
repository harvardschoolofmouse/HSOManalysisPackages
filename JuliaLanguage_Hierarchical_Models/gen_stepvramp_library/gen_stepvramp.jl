function get_mean_slope(traces; verbose=false)
    try
        ii = [traces[i][:slope] for i = 1:length(traces)]
        if verbose
            println("Mean slope: ", round(sum(ii)/length(ii),2), "(",minimum(ii),",", maximum(ii), ")")
        end
        (sum(ii)/length(ii), minimum(ii), maximum(ii))
    catch
    end
end;
function get_mean_intercept(traces; verbose=false)
    try
        ii = [traces[i][:intercept] for i = 1:length(traces)]
        if verbose
            println("Mean intercept: ", round(sum(ii)/length(ii),2), "(",minimum(ii),",", maximum(ii), ")")
            render_distribution(traces[1],ii, xl="slope")
        end
        (sum(ii)/length(ii), minimum(ii), maximum(ii))
    catch
        try
            ii = [traces[i][:intercept_flat] for i = 1:length(traces)]
            if verbose
                println("Mean intercept: ", round(sum(ii)/length(ii),2), "(",minimum(ii),",", maximum(ii), ")")
                render_distribution(traces[1],ii, xl="intercept") 
            end
            (sum(ii)/length(ii), minimum(ii), maximum(ii))
        catch
        end
    end
end;
function get_break_points(traces; verbose=false)
   try
        ii = [get_args(traces[i])[1][traces[i][:step_position]] for i = 1:length(traces)]
        if verbose
            println("Mean break point: ", round(sum(ii)/length(ii),digits=2), "s (",minimum(ii),",", maximum(ii), ")")
            render_distribution(traces[1], ii, xl="break point (s)") 
        end
        (sum(ii)/length(ii), minimum(ii), maximum(ii))
    catch
    end 
end;


function get_step_time_s(trace; xs,step_position)
    if !isnothing(trace)
        (xs,_) = Gen.get_args(trace)
        step_position = trace[:step_position]
    end
    step_position = round(Int,step_position)
    step_time_s = (xs[step_position]+xs[step_position])/2
    step_time_s
end;


struct RANSACParams
    """the number of random subsets to try"""
    iters::Int

    """the number of points to use to construct a hypothesis"""
    subset_size::Int

    """the error threshold below which a datum is considered an inlier"""
    eps::Float64
    
    function RANSACParams(iters, subset_size, eps)
        if iters < 1
            error("iters < 1")
        end
        new(iters, subset_size, eps)
    end
end


function ransac(xs::Vector{Float64}, ys::Vector{Float64}, params::RANSACParams)
    best_num_inliers::Int = -1
    best_slope::Float64 = NaN
    best_intercept::Float64 = NaN
    for i=1:params.iters
        # select a random subset of points
        rand_ind = StatsBase.sample(1:length(xs), params.subset_size, replace=false)
        subset_xs = xs[rand_ind]
        subset_ys = ys[rand_ind]
        
        # estimate slope and intercept using least squares
        A = hcat(subset_xs, ones(length(subset_xs)))
        slope, intercept = A\subset_ys
        
        ypred = intercept .+ slope * xs

        # count the number of inliers for this (slope, intercept) hypothesis
        inliers = abs.(ys - ypred) .< params.eps
        num_inliers = sum(inliers)

        if num_inliers > best_num_inliers
            best_slope, best_intercept = slope, intercept
            best_num_inliers = num_inliers
        end
    end

    # return the hypothesis that resulted in the most inliers
    (best_slope, best_intercept)
end;

function ransac_dynamic_noise(xs::Vector{Float64}, ys::Vector{Float64}, params::RANSACParams, noise)
    best_num_inliers::Int = -1
    best_slope::Float64 = NaN
    best_intercept::Float64 = NaN
    for i=1:params.iters
        # select a random subset of points
        rand_ind = StatsBase.sample(1:length(xs), params.subset_size, replace=false)
        subset_xs = xs[rand_ind]
        subset_ys = ys[rand_ind]
        
        # estimate slope and intercept using least squares
        A = hcat(subset_xs, ones(length(subset_xs)))
        slope, intercept = A\subset_ys
        
        ypred = intercept .+ slope * xs

        # count the number of inliers for this (slope, intercept) hypothesis
        inliers = abs.(ys - ypred) .< noise #params.eps
        num_inliers = sum(inliers)

        if num_inliers > best_num_inliers
            best_slope, best_intercept = slope, intercept
            best_num_inliers = num_inliers
        end
    end

    # return the hypothesis that resulted in the most inliers
    (best_slope, best_intercept)
end;

@gen function ransac_proposal_dynamic_noise(prev_trace, xs, ys)
    noise = prev_trace[:noise]
    base_slope = 0.08
    base_intercept = 0.05
    subset_size = ceil(length(xs)/10)
    (slope_guess, intercept_guess) = ransac_dynamic_noise(xs, ys, RANSACParams(10, subset_size, 1.),noise)
    slope ~ normal(slope_guess, base_slope/10)
    intercept ~ normal(intercept_guess, base_intercept/10)
end;
@gen function ransac_dn_update(tr, ys, u, iters)    
    (xs,) = get_args(tr)
    n = length(xs)
    # Use RANSAC to (potentially) jump to a better line
    # from wherever we are
    (tr, _) = mh(tr, ransac_proposal_dynamic_noise, (xs, ys))
    # Spend a while refining the parameters, using Gaussian drift
    # to tune the slope and intercept, and resimulation for the noise
    # and outliers.
    for j=1:iters
        (tr, _) = mh(tr, select(:noise))
        (tr, _) = mh(tr, line_proposal, ())
    end
    tr
end;    
function ransac_inference_dn(generative_fxn, genfxn_args, xs, ys; iters=20)
    observations = make_constraints(ys)
    #     
    #     Get initial ransac proposal
    #     
    base_slope = 0.08
    base_intercept = 0.05
    subset_size = ceil(length(xs)/10)
    (slope, intercept) = ransac(xs, ys, RANSACParams(10, subset_size, 1.))
    slope_intercept_init = choicemap()
    slope_intercept_init[:slope] = slope
    slope_intercept_init[:intercept] = intercept
    (tr_start, _) = generate(generative_fxn, genfxn_args, merge(observations, slope_intercept_init))
    tr = tr_start;
    for iter=1:5
        tr = ransac_dn_update(tr, ys, iter, iters)
    end
    tr
end;

function flat_ransac(xs::Vector{Float64}, ys::Vector{Float64}, params::RANSACParams)
    best_num_inliers::Int = -1
    best_intercept::Float64 = NaN
    for i=1:params.iters
        # select a random subset of points
        rand_ind = StatsBase.sample(1:length(xs), params.subset_size, replace=false)
        subset_xs = xs[rand_ind]
        subset_ys = ys[rand_ind]
        
        # estimate intercept
        intercept = sum(subset_ys)/length(subset_ys)
        ypred = intercept

        # count the number of inliers for this (slope, intercept) hypothesis
        inliers = abs.(ys .- ypred) .< params.eps
        num_inliers = sum(inliers)

        if num_inliers > best_num_inliers
            best_intercept = intercept
            best_num_inliers = num_inliers
        end
    end

    # return the hypothesis that resulted in the most inliers
    best_intercept
end;
function flat_ransac_dynamic_noise(xs::Vector{Float64}, ys::Vector{Float64}, params::RANSACParams, noise)
    best_num_inliers::Int = -1
    best_intercept::Float64 = NaN
    for i=1:params.iters
        # select a random subset of points
        rand_ind = StatsBase.sample(1:length(xs), params.subset_size, replace=false)
        subset_xs = xs[rand_ind]
        subset_ys = ys[rand_ind]
        
        # estimate intercept
        intercept = sum(subset_ys)/length(subset_ys)
        ypred = intercept

        # count the number of inliers for this (slope, intercept) hypothesis
        inliers = abs.(ys .- ypred) .< noise
        num_inliers = sum(inliers)

        if num_inliers > best_num_inliers
            best_intercept = intercept
            best_num_inliers = num_inliers
        end
    end

    # return the hypothesis that resulted in the most inliers
    best_intercept
end;
function flat_ransac_inference_dn(generative_fxn, genfxn_args, xs, ys; iters=20, verbose=false)
    observations = make_constraints(ys)
    #     
    #     Get initial ransac proposal
    #     
    base_intercept = 0.05
    subset_size = ceil(length(xs)/10)
    intercept_flat = flat_ransac(xs, ys, RANSACParams(10, subset_size, 1.))
    intercept_init = choicemap()
    intercept_init[:intercept_flat] = intercept_flat
    (tr_start, _) = generate(generative_fxn, genfxn_args, merge(observations, intercept_init))
    tr = tr_start;
    for iter=1:5
        tr = flat_ransac_dn_update(tr, ys, iter, iters)
    end
    tr
end;
function step_ransac_dynamic_noise(xs::Vector{Float64}, ys::Vector{Float64}, params::RANSACParams, noise, step_position)
    best_num_inliers_LHS::Int = -1
    best_LHS::Float64 = NaN
    best_num_inliers_RHS::Int = -1
    best_RHS::Float64 = NaN
    for i=1:params.iters
        #
        # Get best setting for the LHS
        #
        # select a random subset of points
        if params.subset_size > length(xs[1:step_position])
            subset_size = length(xs[1:step_position])
        else
            subset_size = params.subset_size
        end
        rand_ind_LHS = StatsBase.sample(1:length(xs[1:step_position]), subset_size, replace=false)
        subset_xs_LHS = xs[rand_ind_LHS]
        subset_ys_LHS = ys[rand_ind_LHS]
        
        # estimate amplitude of LHS
        amp_LHS = sum(subset_ys_LHS)/length(subset_ys_LHS)
        ypred = amp_LHS

        # count the number of inliers for this hypothesis
        inliers_LHS = abs.(ys .- ypred) .< noise
        num_inliers_LHS = sum(inliers_LHS)

        if num_inliers_LHS > best_num_inliers_LHS
            best_LHS = amp_LHS
            best_num_inliers_LHS = num_inliers_LHS
        end        
        #
        # Get best setting for the RHS
        #
        # select a random subset of points
        if params.subset_size > length(xs[step_position+1:end])
            subset_size = length(xs[step_position+1:end])
        else
            subset_size = params.subset_size
        end
        rand_ind_RHS = StatsBase.sample(1:length(xs[step_position+1:end]), subset_size, replace=false)
        subset_xs_RHS = xs[rand_ind_RHS]
        subset_ys_RHS = ys[rand_ind_RHS]
        
        # estimate intercept
        amp_RHS = sum(subset_ys_RHS)/length(subset_ys_RHS)
        ypred = amp_RHS

        # count the number of inliers for this hypothesis
        inliers_RHS = abs.(ys .- ypred) .< noise
        num_inliers_RHS = sum(inliers_RHS)

        if num_inliers_RHS > best_num_inliers_RHS
            best_RHS = amp_RHS
            best_num_inliers_RHS = num_inliers_RHS
        end
    end
    # return the hypothesis that resulted in the most inliers
    (best_LHS,best_RHS)
end;
function step_model_inference(generative_fxn, genfxn_args, xs, ys; iters=20)
    observations = make_constraints(ys)
    #     
    #     Get initial ransac proposal
    #     
    (proposed_choices, _, _) = propose(step_proposal_init, (xs, ys))
    (tr_start, _) = generate(generative_fxn, genfxn_args, merge(proposed_choices, observations))
    tr = tr_start;
    for iter=1:5
        tr = step_update_course(tr, ys, iter, iters)
    end
    tr
end;

@gen function regression_variable_noise(xs::Vector{Float64}, noiseCenter)
    # default noiseCenter = 1
    #
    # First, generate some parameters of the model. We make these
    # random choices, because later, we will want to infer them
    # from data. The distributions we use here express our assumptions
    # about the parameters: we think the slope and intercept won't be
    # too far from 0; that the noise is relatively small; and that
    # the proportion of the dataset that don't fit a linear relationship
    # (outliers) could be anything between 0 and 1.
    slope ~ normal(0, 0.1)
    intercept ~ normal(0, 0.05)
    noise ~ gamma(noiseCenter,noiseCenter)#gamma(0.1,0.1) #gamma(0.025,0.05) #mean 0.00125 gamma(0.1,0.1) #mean 0.01, gamma(1, 1)-mean1
    
    # Next, we generate the actual y coordinates.
    n = length(xs)
    ys = Vector{Float64}(undef, n)
    
    for i = 1:n
        (mu, std) = (xs[i] * slope + intercept, noise)
        # Sample a y value for this point
        ys[i] = ({:data => i => :y} ~ normal(mu, std))
    end
    ys
end;
@gen function line_proposal(current_trace)
    baseslope = 0.02;
    slope ~ normal(current_trace[:slope], baseslope/2)
    intercept ~ normal(current_trace[:intercept], baseslope/2)
end;


@gen function flat_regression_variable_noise(xs::Vector{Float64}, noiseCenter)
    # default noiseCenter = 1
    #
    # Just as before, only we are only allowing an offset term
    #     
    intercept_flat ~ normal(0, 0.05)
    noise_flat ~ gamma(noiseCenter,noiseCenter)
    
    # Next, we generate the actual y coordinates.
    n = length(xs)
    ys = Vector{Float64}(undef, n)
    
    for i = 1:n
        (mu, std) = (intercept_flat, noise_flat)
        # Sample a y value for this point
        ys[i] = ({:data => i => :y} ~ normal(mu, std))
    end
    ys
end;

@gen function flat_ransac_proposal_dynamic_noise(prev_trace, xs, ys)
    noise = prev_trace[:noise_flat]
    base_intercept = 0.02
    subset_size = ceil(length(xs)/10)
    (intercept_guess) = flat_ransac_dynamic_noise(xs, ys, RANSACParams(10, subset_size, 1.),noise)
    intercept_flat ~ normal(intercept_guess, base_intercept/2)
end;
@gen function flat_line_proposal(current_trace)
    base_intercept = 0.02;
    intercept_flat ~ normal(current_trace[:intercept_flat], base_intercept/2)
end;

@gen function flat_ransac_dn_update(tr, ys, u, iters)    
    (xs,) = get_args(tr)
    n = length(xs)
    # Use RANSAC to (potentially) jump to a better line
    # from wherever we are
    (tr, _) = mh(tr, flat_ransac_proposal_dynamic_noise, (xs, ys))
    # Spend a while refining the parameters, using Gaussian drift
    # to tune the slope and intercept, and resimulation for the noise
    # and outliers.
    for j=1:iters
        (tr, _) = mh(tr, select(:noise_flat))
        (tr, _) = mh(tr, flat_line_proposal, ())
    end
    tr
end;  


@gen function generate_step_model(xs::Vector{Float64}, noiseCenter)
    step_position ~ uniform_discrete(1, length(xs)-1)
    left_segment_amp ~ normal(0, 0.25)
    right_segment_amp ~ normal(0, 0.25)
    noise_step ~ gamma(noiseCenter,noiseCenter)
    
    n = length(xs)
    ys = Vector{Float64}(undef, n)
    
    for i = 1:n
        if i <= step_position
            mu = left_segment_amp
        else
            mu = right_segment_amp
        end
        std = noise_step
        ys[i] = ({:data => i => :y} ~ normal(mu, std))
    end
    (ys, get_step_time_s(nothing, xs=xs, step_position=step_position))
end;

@gen function step_ransac_proposal_dynamic_noise(prev_trace, xs, ys)
    noise = prev_trace[:noise_step]
    step_position = prev_trace[:step_position]
    base_amp = 0.02
    subset_size = ceil(length(xs)/20) # div by twenty because we are doing for 2 subsets...
    (guess_LHS,guess_RHS) = step_ransac_dynamic_noise(xs, ys, RANSACParams(10, subset_size, 1.),noise,step_position)
    left_segment_amp ~ normal(guess_LHS, base_amp)
    right_segment_amp ~ normal(guess_RHS, base_amp)    
end;
@gen function gd_step_proposal(current_trace)
    base_amp = 0.02
    step_position ~ normal(prev_trace[:step_position], prev_trace[:step_position]/2)
    left_segment_amp ~ normal(prev_trace[:left_segment_amp], base_amp)
    right_segment_amp ~ normal(prev_trace[:right_segment_amp], base_amp)
    noise_step ~ normal(prev_trace[:noise_step], prev_trace[:noise_step]/2)
end;
@gen function step_update_course(tr, ys, u, iters)   
    #
    #  Test out the hypothesis space coursely, updating each param randomly, but using RANSAC to guide
    #
    (xs,) = get_args(tr)
    n = length(xs)
    # Use RANSAC to (potentially) jump to a better parameterization of the split-point
    (tr, _) = mh(tr, step_ransac_proposal_dynamic_noise, (xs, ys))
    # Spend a while refining the parameters, using Gaussian drift
    # to tune the slope and intercept, and resimulation for the noise
    # and outliers.
    for j=1:iters
        (tr, _) = mh(tr, select(:step_position))
        (tr, _) = mh(tr, select(:left_segment_amp))
        (tr, _) = mh(tr, select(:right_segment_amp))
        (tr, _) = mh(tr, select(:noise_step))
    end
    #
    tr
end;  
@gen function step_update_fine(tr, ys, u, iters) 
    # Use RANSAC to (potentially) jump to a better parameterization of the split-point
    (tr, _) = mh(tr, step_ransac_proposal_dynamic_noise, (xs, ys))
    # Spend a while refining the parameters, using Gaussian drift
    for j=1:iters
        (tr, _) = mh(tr, gd_step_proposal, ())
    end
    #
    tr
end;  
@gen function step_proposal_init(xs, ys, step_position=nothing)
    #     
    #     we will use what we know about the dataset--it seems to ramp on average.
    #     as such, we will propose a LHS value that's near the min of the LHS,
    #     and a RHS value near the max of the RHS
    #     
    #     Start by randomly sampling the step_position
    if isnothing(step_position)
        step_position = round(Int,0.5*length(xs))
    end
    lhs_min = minimum(ys[1:step_position])
    rhs_max = maximum(ys[step_position+1:end])
    left_segment_amp ~ normal(lhs_min, 0.01)
    right_segment_amp ~ normal(rhs_max, 0.01)   
end;


@gen function generate_hierarchical_model(xs, noiseCenter)
    if (model_choice_is_slope = {:model_choice_is_slope} ~ bernoulli(0.5))
        {*} ~ regression_variable_noise(xs, noiseCenter)
    else
        {*} ~ generate_step_model(xs, noiseCenter)
    end
end;

struct OutputParameters
    slopes::Vector{Float64}
    intercepts::Vector{Float64}
    step_positions::Vector{Float64}
end

function extract_mean_parameters(traces; verbose=true)
    output_params = OutputParameters(Vector{Float64}(undef,length(traces)),
        Vector{Float64}(undef,length(traces)),
        Vector{Float64}(undef,length(traces)))
    
    for i = 1:length(traces)
        try
            output_params.slopes[i] = traces[i][:slope]
            output_params.intercepts[i] = traces[i][:intercept]
        catch
            output_params.slopes[i] = NaN
            output_params.intercepts[i] = NaN
        end
        try 
            output_params.step_positions[i] = traces[i][:step_position]
        catch
            output_params.step_positions[i] = NaN
        end
    end
    return output_params
end;


struct OutputParameters
    slopes::Vector{Float64}
    intercepts::Vector{Float64}
    step_positions::Vector{Float64}
end

function extract_mean_parameters(traces; verbose=true)
    output_params = OutputParameters(Vector{Float64}(undef,length(traces)),
        Vector{Float64}(undef,length(traces)),
        Vector{Float64}(undef,length(traces)))
    
    for i = 1:length(traces)
        try
            output_params.slopes[i] = traces[i][:slope]
            output_params.intercepts[i] = traces[i][:intercept]
        catch
            output_params.slopes[i] = NaN
            output_params.intercepts[i] = NaN
        end
        try 
            output_params.step_positions[i] = traces[i][:step_position]
        catch
            output_params.step_positions[i] = NaN
        end
    end
    return output_params
end;

function derivative_constraint(ys::Vector{Float64})
    dy = ys[2:end] - ys[1:end-1]
    negs = sum(y -> y < 0, dy)
    pos = sum(y -> y > 0, dy)
    if negs > pos
        return (dy, -1)
    else
        return (dy, +1)
    end
end;

function test_proposal(proposal, xs::Vector{Float64}, ys::Vector{Float64}; n::Int=10)
    figure(figsize=(8,3))
    axs=[]
    push!(axs, subplot(1,2,1))
    push!(axs, subplot(1,2,2))
    for i = 1:n
        (pc, _, _) = Gen.propose(proposal, (xs, ys))
        render_proposal(xs, ys, pc, axs=axs)
    end
    return nothing
end;



@gen function ransac_assisted_model_selection_proposal(xs::Vector{Float64}, ys::Vector{Float64})
    # characterize the data
    (global_deriv, global_deriv_sign) = derivative_constraint(ys)   
    ys_std = StatsBase.std(ys)
    globalmin = minimum(ys)
    globalmax = maximum(ys)
    xn = minimum(xs)
    xx = maximum(xs)
    maxslope = (globalmax-globalmin)/(xx-xn)
    # randomly sample x and y values based on y amplitude, then we can trace around this with some jitter
    (y_idx_Q1, y_Q1_std) = inverse_sample(ys, consideration_interval=[0.,0.25])
    (y_idx_Q2, y_Q2_std) = inverse_sample(ys, consideration_interval=[0.25,0.5])
    (y_idx_Q3, y_Q3_std) = inverse_sample(ys, consideration_interval=[0.5,0.75])
    (y_idx_Q4, y_Q4_std) = inverse_sample(ys, consideration_interval=[0.75,1.0])
    y_Q1 = ys[y_idx_Q1]
    y_Q2 = ys[y_idx_Q2]
    y_Q3 = ys[y_idx_Q3]
    y_Q4 = ys[y_idx_Q4]
    # randomly sample x positions based on derivative dy
    (x_idx_Q4, x_Q4_std) = inverse_sample(global_deriv, consideration_interval=[0.95,1.])
    if isnan(x_idx_Q4) || isnan(x_Q4_std)
        # warning("encoountered unexpected nan. If in std, is because we had only one timepoint to choose from, so setting to zero")
        # warning(join(["x_idx_Q4=", x_idx_Q4]))
        # warning(join(["x_Q4_std=", x_Q4_std]))
        x_Q4_std = 0;
    end
    # sample a step_position
    step_position ~ uniform_discrete(round(Int, x_idx_Q4-x_Q4_std), round(Int, x_idx_Q4+x_Q4_std))
    # sample noise, noise_step  
    (a,b) = get_beta_params(std(ys))#StatsBase.std(global_deriv)) # want to high-pass filter the slope model to get noise
    noise ~ beta(a,b)
    if step_position<=1
        (a,b) = get_beta_params(mean([StatsBase.std(ys[1:step_position+1]),StatsBase.std(ys[step_position+2:end])])) # want to use a more local estimate of noise
    elseif step_position>=length(ys)-1
        (a,b) = get_beta_params(mean([StatsBase.std(ys[1:step_position-1]),StatsBase.std(ys[step_position:end])])) # want to use a more local estimate of noise
    else #standard case
        (a,b) = get_beta_params(mean([StatsBase.std(ys[1:step_position]),StatsBase.std(ys[step_position+1:end])])) # want to use a more local estimate of noise
    end
    if isnan(a) || isnan(b)
        # warning(join(["a=", a]))
        # warning(join(["b=", b]))
        # warning(join(["setting noise-step to noise."]))
        noise_step = noise
    else
        noise_step ~ beta(a,b)#sum([y_Q1_std,y_Q4_std])/2)
    end
    
    # sample the step qualities with similar to generative model with ransac
    subset_size = ceil(length(xs)/20) # div by twenty because we are doing for 2 subsets...
    left_segment_amp ~ normal(mean(ys[1:step_position]), std(ys[1:step_position]))
    right_segment_amp ~ normal(mean(ys[step_position+1:end]), std(ys[step_position+1:end]))
    
    # sample slope with ransac
    base_slope = maxslope/2 # intermediate data-driven slope
    subset_size = ceil(length(xs)/10)
    (slope_guess, intercept_guess) = ransac_dynamic_noise(xs, ys, RANSACParams(10, subset_size, 1.),noise)
    slope ~ normal(slope_guess, base_slope/10)
    
    base_intercept = (globalmin-slope*xn + mean(ys))/2 # intermediate data-driven intercept
    intercept ~ normal(intercept_guess, base_intercept/10)
    
    return nothing
end;

@gen function data_driven_inference(model::GenerativeFunction, model_args::Tuple, 
        ys::Vector{Float64}, proposal::GenerativeFunction, proposal_args::Tuple, 
        amount_of_computation::Int)
    observations = make_constraints(ys)
    # Call importance_resampling to obtain a likely trace consistent
    # with our observations and data-driven proposal
    (traces, WT) = importance_resampling(model,
        model_args, observations,
        proposal, proposal_args,
        amount_of_computation, verbose=false)
    
#     println("WT:", WT)
    iters = 0
    flag = false
    while isinf(WT) || WT < -100 || isnan(WT)
        if isnan(WT) || iters >=3
#             println("FLAG: Did not converge")
            flag = true
            break
        end
        amount_of_computation = amount_of_computation+10
        println("WARNING: did not converge. Trying again with amt_of_computation=", amount_of_computation)
        (traces, WT) = importance_resampling(model,
            model_args, observations,
            proposal, proposal_args,
            amount_of_computation, verbose=false)
        iters += 1
    end
#     println("Flag: ", flag)
    result = Result(traces,flag)
    return result
end;
function data_driven_probability_of_model_class(model::GenerativeFunction, 
        model_args::Tuple, ys::Vector{Float64}, proposal::GenerativeFunction, 
        proposal_args::Tuple; amount_of_computation::Int=10000, ntraces::Int=100)
    result = [data_driven_inference(model, model_args, ys, proposal, proposal_args, 
        amount_of_computation) for _=1:ntraces];
    traces = [result[i].traces for i=1:length(result)]
    flags = [result[i].flags for i=1:length(result)]
    function percent_slope(traces=traces, flags=flags)
        c = 0
        ntrue = length(findall(x->x==false, flags))
        if ntrue == 0
            return NaN
        else
            for ii = 1:length(traces)
                if traces[ii][:model_choice_is_slope] && flags[ii] == false
                    c=c+1
                end
            end
            return c/ntrue
        end
    end
    return (percent_slope(traces), traces, flags)
end;
function test_out_inf(xx,yy)
    figure(figsize=(1,1))
    plot(xx, yy)
    @time (percent_slope, traces) =data_driven_probability_of_model_class(generate_hierarchical_model, 
            (xx,1), yy, ransac_assisted_model_selection_proposal, 
            (xx,yy); amount_of_computation=10, ntraces=12)

    println(string("With 10 iters, p(slope NOT step): ", percent_slope))
    @time (percent_slope, traces2) =data_driven_probability_of_model_class(generate_hierarchical_model, 
            (xx,1), yy, ransac_assisted_model_selection_proposal, 
            (xx,yy); amount_of_computation=100, ntraces=12)
    println(string("With 100 iters, p(slope NOT step): ", percent_slope))
    
    grid(render_trace, traces[1:12])
    grid(render_trace, traces2[1:12])
    percent_slope_test = percent_slope;
    traces_test = traces2;
    println(" ")
    return (percent_slope_test, traces_test)
end;