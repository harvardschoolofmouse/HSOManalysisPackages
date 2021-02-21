#
#  The following functions to do hazard fitting processing
#
#. 500 ms hazard normalized comparison:
# (hazfit250_2_50, p_by_trialfit250_2_50, ltfit250_2_50) = get_fit_hazard(results_200, seshNo=1, p_prior=0.00386*2, ndp_per_sample=50, pre_normalize_p=false, normalize=true);
#


function average_df(predictorIDs, idxs, df, yID)
    # extract the part of the dataframe
    df2 = df[idxs,predictorIDs]
    df2[yID]=df[idxs,yID]
    meanpred = []
    # take the mean
    newdf = DataFrame()
    for i = 1:length(predictorIDs)
        mm = mean(df2[predictorIDs[i]])
        push!(meanpred, mm)
        newdf[predictorIDs[i]] = [mm]
    end
    newdf[yID] = [maximum(df2[yID]) == 1.]
   return newdf 
end
function get_prior(collated_results, seshNo, yID)
    df = collated_results.results[seshNo].modelData[1].df[1]
    t = countmap(df[yID])[true]
    f = countmap(df[yID])[false]
    prior = t / (t+f)
end
function extract_behavior_distribution(og_df::DataFrame; verbose=true)
    #
    # If the licktime was normalized, need to correct this back to time in sec...
    #
    
    trialIDs = unique(og_df.TrialNo)
    
    lick_times = []
    dps = []
    for i in trialIDs
        # get the datapoints corresponding to this trial
        dp = findfirst(x->x==i, og_df.TrialNo)
        push!(dps, dp)
        push!(lick_times, og_df.LickTime[dp])
    end
    min_time = minimum(og_df[dps[2:end].-1, :X])
    max_time = maximum(og_df[dps[2:end].-1, :X])
    lick_times = (lick_times .* max_time) .+ min_time
    if verbose
        println("min licktime=", min_time, " max licktime=", max_time)
        figure(figsize=(3,3))
        render_distribution(lick_times, "lick time (s)", t="true normalized distribution", bins=50, ax=gca())
        xticks(0:17)
        xlim([0,17])
    end
    return lick_times
end
function IRT_byOpportunity(d; edges=0:0.25:17, verbose=true)
    #
    # Given a distribution, calculate a hazard fxn
    #
    result = fit(Histogram, Float64.(d), edges)
    N = result.weights
    IRTbyOP = nanmat(length(N),1);
    for nn = 1:length(IRTbyOP)
        IRTbyOP[nn] = N[nn]/(sum(N[nn:end]));
    end
    if verbose
        figure(figsize=(8,3))
        suptitle("IRT by Opportunity")
        subplot(1,2,1)
        plot(edges[1:end-1], IRTbyOP, "k-")
        ylabel("Hazard Rate")
        xlabel("Time relative to cue (s)")
        
        subplot(1,2,2)
        plot(edges[1:end-1-round(Int,length(edges)/2)], IRTbyOP[1:end-round(Int,length(edges)/2)], "k-")
        xticks(0:7)
        xlim([0,7])

        xlabel("Time relative to cue (s)")
    end
    return IRTbyOP
end
function predict_logit(yID, predictors, th, df; plotOn=true)
    th = hcat(Float64.(th),)
    yfit = []
    for i =1:nrow(df)
        xx = df[i, predictors[1]]
        v = vec(convert(Array, xx))
        xx = hcat(vcat([1],v),)
        x = th'*xx
        yfiti = (exp.(x))./(1 .+ exp.(x)) 
        push!(yfit, yfiti[1])
    end
    
    sortedidx = sortperm(df[:Y])
    
    yactual = df[yID]
    if plotOn
        figure(figsize=(3,3))
        plot(vec(df[sortedidx,:Y]), vec(yactual[sortedidx]), "k.")
        plot(vec(df[sortedidx,:Y]), vec(yfit[sortedidx]), "r.")
        xlabel("DA level")
        ylabel("p(move|predictors)")
        
#         figure(figsize=(3,3))
#         plot_mean_fit_vs_predictor2(df, yfit, yactual, :Y, :Hx10, :LickState)
    end
    return (yfit, yactual)
end

function get_fit_hazard(collated_results; seshNo=1, p_prior=1.0, ndp_per_sample=1, pre_normalize_p=false, normalize=false)
    #
    pprior_star = get_prior(collated_results, seshNo, :LickState)
    println("pprior_star=", pprior_star)
    yID = collated_results.results[seshNo].modelData[1].yID[1]
    th = collated_results.results[seshNo].modelData[1].th_means[1]
    df = collated_results.results[seshNo].modelData[1].df[1]
    X = df.X
    predictors = collated_results.results[seshNo].modelData[1].predictors # must be a vector
    println(predictors)
    # Let's get the range of probabilities the model can do on the whole dataset...
    all_p,_ = predict_logit(yID, predictors, th, df,plotOn=false)
    min_p = minimum(all_p)
    max_p = maximum(all_p)
    println("min p_fit=", min_p)
    println("max p_fit=", max_p)
    

    trialIDs = unique(df.TrialNo)
    p_by_trial = nanmat(length(trialIDs), 1701)
    for ii =1:length(trialIDs)
        i = trialIDs[ii]
        # get the datapoints corresponding to this trial
        dps = findall(x->x==i, df.TrialNo)
        nj = ceil(Int,length(dps)/ndp_per_sample) 

        foundlick=false;
        for j = 1:nj
            mn = ndp_per_sample*(j-1)+1
            mx = round(Int,minimum([ndp_per_sample*(j), length(dps)]))
            # Average the signals within this time window...
            dpsrange = dps[mn:mx]

            # make an average dataframe of this set
            avedf = average_df(predictors[1], dpsrange, df, yID)

            # look at the DA signal and predictors at that time
            # calculate its p(move | dopamine, factors) from the model 
            prediction_p,_ = predict_logit(yID, predictors, th, avedf,plotOn=false)
            p = prediction_p[1]

            # Let's correct this p for our uncertainty. 
            #. If it's close to min prediction on our dataset, let's set to zero. 
            #. If it's close to max set close to 1
            if pre_normalize_p
                p_corrected = (p-min_p)/max_p
            else 
                p_corrected = p
            end
            #
            # Scale everything by the prior... (p_prior)
            #
            p_corrected = p_corrected*p_prior

            p_by_trial[ii,j] = p_corrected
        end
    end
    haz = nanmean_mat(p_by_trial)
    lt = extract_behavior_distribution(df)
    haz_results(haz, lt, ndp_per_sample, normalize=normalize)
    return (haz, p_by_trial, lt, ndp_per_sample)
end
function haz_results(haz, lt, ndp_per_sample; normalize=false)
    # get the nanmean of the vectors
    IRTbyOP = IRT_byOpportunity(lt, edges=0:0.01*ndp_per_sample:17)
    figure(figsize=(8,3))
    xs = range(0.01, step=0.01*ndp_per_sample, stop=7)
    
    if normalize
        IRTbyOP = (IRTbyOP .- minimum(IRTbyOP[1:length(xs)])) ./ maximum((IRTbyOP .- minimum(IRTbyOP[1:length(xs)]))[1:length(xs)])
        haz = (haz .- minimum(haz[1:length(xs)])) ./ maximum((haz .- minimum(haz[1:length(xs)]))[1:length(xs)])
    end
    suptitle(join(["Rsq=", Rsq(IRTbyOP[1:length(xs)], haz[1:length(xs)])]))
    xs = range(0.01, step=0.01*ndp_per_sample, stop=17)
    subplot(1,2,1)
    plot(xs, haz[1:length(xs)], "k-")
    xlim([0,maximum(xs)])
    xticks(0:17)
    title("logit fit hazard")
    xlabel("time (s)") 
    
    subplot(1,2,2)
    xs = range(0.01, step=0.01*ndp_per_sample, stop=7)
    plot(xs, haz[1:length(xs)], "k-")
    xlim([0,maximum(xs)])
    xticks(0:7)
    title("logit fit hazard")
    xlabel("time (s)") 
    
    figure(figsize=(3,3))
    suptitle(join(["Rsq=", round(Rsq(IRTbyOP[1:length(xs)], haz[1:length(xs)]), digits=3)]))
    xs = range(0.01, step=0.01*ndp_per_sample, stop=7)
    
    edges=0:0.01*ndp_per_sample:17
    plot(edges[1:length(xs)], IRTbyOP[1:length(xs)], "k-", label="true hazard")
    plot(xs, haz[1:length(xs)], "r-", label="fit hazard")
    xlabel("time (s)") 
    ylabel("normalized hazard") 
    legend(loc="right", bbox_to_anchor=(1.7, 0.5))
    xticks(0:1:7)
    xlim([0,7])
    
    
    return haz
end