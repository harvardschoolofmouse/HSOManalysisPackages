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
    
    
    
    yactual = df[yID]
    if plotOn
        sortedidx = sortperm(df[:Y])
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

function haz_results_composite(hazs, lts, seshCodes; ndp_per_sample=50, normalize=true, figname="", figurePath="")
    # fig 1 is te normalized hazard true, fit and overlay
    fig=figure(figsize=(12,3))
    axIRT = subplot(1,3,1)
    axHaz = subplot(1,3,2)
    axOverlay = subplot(1,3,3)
    # fig 2 is a misguided attempt at CI. This will never work
    fig2=figure(figsize=(12,3))
    ax2IRT = subplot(1,3,1)
    ax2Haz = subplot(1,3,2)
    ax2Overlay = subplot(1,3,3)
    # fig 3 is going to be where we plot the by-animal averages normalized and the grandave.
    fig3=figure(figsize=(12,3))
    ax3IRT = subplot(1,3,1)
    ax3Haz = subplot(1,3,2)
    ax3Overlay = subplot(1,3,3)


    meanIRT = []
    meanHaz = []
    goodidx = findall(x->!isnan(x[1]), hazs)
    failidx = findall(x->isnan(x[1]), hazs)

    # we will store these and save properly to have a matlab-compatible variable
    xs = range(0.01, step=0.01*ndp_per_sample, stop=7)
    normIRTs = nanmat(length(hazs), length(1:length(xs)))
    normhazs = nanmat(length(hazs), length(hazs[1]))
    println("normIRTs=", size(normIRTs))
    println("normhazs=", size(normhazs))
    for ii=1:length(hazs)#length(goodidx)
        if ii in goodidx
            # ii = goodidx[i]
            i = findall(x->x==ii, goodidx)
            sesh=seshCodes[ii]
            # get the nanmean of the vectors
            IRTbyOP = IRT_byOpportunity(lts[ii], edges=0:0.01*ndp_per_sample:17, verbose=false)   
            xs = range(0.01, step=0.01*ndp_per_sample, stop=7)
            haz=hazs[ii]
            if normalize
                IRTbyOP = (IRTbyOP .- nanmin(IRTbyOP[1:length(xs)])) ./ nanmax((IRTbyOP .- nanmin(IRTbyOP[1:length(xs)]))[1:length(xs)])
                haz = (haz .- nanmin(haz[1:length(xs)])) ./ nanmax((haz .- nanmin(haz[1:length(xs)]))[1:length(xs)])
            end
            println("IRTbyOP=", size(IRTbyOP))
            println("haz=", size(haz))
            # need to keep track to plot mean by mouse...
            normIRTs[ii, :] = IRTbyOP
            normhazs[ii, :] = haz

            println(join([seshCodes[ii], " Rsq=", round(Rsq(IRTbyOP[1:length(xs)], haz[1:length(xs)]), digits=3)]))
            xs = range(0.01, step=0.01*ndp_per_sample, stop=7)

            edges=0:0.01*ndp_per_sample:17
            axIRT.plot(edges[1:length(xs)], IRTbyOP[1:length(xs)], linewidth=0.5, "k-", label=sesh)
            axHaz.plot(xs, haz[1:length(xs)], "r-", linewidth=0.5,label=sesh)
            axOverlay.plot(edges[1:length(xs)], IRTbyOP[1:length(xs)], linewidth=0.5, "k-", label=sesh)
            axOverlay.plot(xs, haz[1:length(xs)], "r-",linewidth=0.5, label=sesh)
            
            if i==1
                meanIRT = IRTbyOP
                meanHaz = haz
            else
                meanIRT = hcat(meanIRT, IRTbyOP)
                meanHaz = hcat(meanHaz, haz)
        
            end
        else

        end
    end





    allIRT = meanIRT
    allHaz = meanHaz
    # sort these out to get 95% CI
    for i = 1:size(meanIRT)[1] # for each timepoint
        allIRT[i,:] = sort(allIRT[i,:])
        allHaz[i,:] = sort(allHaz[i,:])
    end
    # get means:
    xs = range(0.01, step=0.01*ndp_per_sample, stop=7)
    edges=0:0.01*ndp_per_sample:17
    meanIRT = nanmean_mat(meanIRT, 2)
    meanHaz = nanmean_mat(meanHaz, 2)



    # now we need to get the by animal situation...
    names = ["B1", "B2", "B3", "b5", "B6", "H3", "H4", "H5", "h6", "h7", "H14", "H15"]
    IRTbymouse = nanmat(length(names), size(normIRTs)[2])
    hazbymouse = nanmat(length(names), size(normhazs)[2])
    for i = 1:length(names)
        name = names[i]
        ix = findall( x -> occursin(name, x), seshcodes)

        IRT_thismouse = normIRTs[ix, :]
        IRTbymouse[i,:] = IRT_thismouse
        ax3IRT.plot(edges[1:length(xs)],nanmean_mat(IRT_thismouse, 1)[1:length(xs)], label=name)
        ax3Overlay.plot(edges[1:length(xs)], nanmean_mat(IRT_thismouse, 1)[1:length(xs)], label=name)


        haz_thismouse = normhazs[ix, :]
        hazbymouse[i,:] = haz_thismouse
        ax3Haz.plot(edges[1:length(xs)], nanmean_mat(haz_thismouse, 1)[1:length(xs)], label=name)
        ax3Overlay.plot(edges[1:length(xs)], nanmean_mat(haz_thismouse, 1)[1:length(xs)], label=name)

    end
    ax3IRT.plot(edges[1:length(xs)], meanIRT[1:length(xs)], "k-", linewidth=3, label="MEAN")
    ax3Haz.plot(xs, meanHaz[1:length(xs)], "r-", linewidth=3, label="MEAN")
    ax3Overlay.plot(edges[1:length(xs)], meanIRT[1:length(xs)], "k-", linewidth=3, label="MEAN")
    ax3Overlay.plot(xs, meanHaz[1:length(xs)], "r-", linewidth=3, label="MEAN")
    ax3Haz.legend()








    idxx_min = round(Int,0.025*size(allIRT)[2])
    idxx_max = round(Int,0.975*size(allIRT)[2])
    if idxx_min < 1
        idxx_min = 1
    end
    println("CIlower index: ",idxx_min, " | CIupper index: ", idxx_max)
    println(size(allIRT))
    ax2IRT.plot(edges[1:length(xs)], allIRT[1:length(xs), idxx_min], linewidth=0.5, "k-", label="CI95-lower")
    ax2Haz.plot(xs, allHaz[1:length(xs), idxx_min], "r-", linewidth=0.5,label="CI95-lower")
    ax2Overlay.plot(edges[1:length(xs)], allIRT[1:length(xs), idxx_min], linewidth=0.5, "k-", label="CI95-lower")
    ax2Overlay.plot(xs, allHaz[1:length(xs), idxx_min], "r-",linewidth=0.5, label="CI95-lower")
    ax2IRT.plot(edges[1:length(xs)], allIRT[1:length(xs), idxx_max], linewidth=0.5, "k-", label="CI95-upper")
    ax2Haz.plot(xs, allHaz[1:length(xs), idxx_max], "r-", linewidth=0.5,label="CI95-upper")
    ax2Overlay.plot(edges[1:length(xs)], allIRT[1:length(xs), idxx_max], linewidth=0.5, "k-", label="CI95-upper")
    ax2Overlay.plot(xs, allHaz[1:length(xs), idxx_max], "r-",linewidth=0.5, label="CI95-upper")
    ax2IRT.plot(edges[1:length(xs)], meanIRT[1:length(xs)], "k-", linewidth=3, label="MEAN")
    ax2Haz.plot(xs, meanHaz[1:length(xs)], "r-", linewidth=3, label="MEAN")
    ax2Overlay.plot(edges[1:length(xs)], meanIRT[1:length(xs)], "k-", linewidth=3, label="MEAN")
    ax2Overlay.plot(xs, meanHaz[1:length(xs)], "r-", linewidth=3, label="MEAN")


    
    
    h = fig.suptitle(join(["Mean Rsq=", round(Rsq(meanIRT[1:length(xs)], meanHaz[1:length(xs)]), digits=3)]), y=1.15)
    h2 = fig2.suptitle(join(["Mean Rsq=", round(Rsq(meanIRT[1:length(xs)], meanHaz[1:length(xs)]), digits=3)]), y=1.15)
    axIRT.plot(edges[1:length(xs)], meanIRT[1:length(xs)], "k-", linewidth=3, label="MEAN")
    axHaz.plot(xs, meanHaz[1:length(xs)], "r-", linewidth=3, label="MEAN")
    axOverlay.plot(edges[1:length(xs)], meanIRT[1:length(xs)], "k-", linewidth=3, label="MEAN")
    axOverlay.plot(xs, meanHaz[1:length(xs)], "r-", linewidth=3, label="MEAN")
    
    
    axIRT.set_title("True Hazard")
    axHaz.set_title("Fit Hazard")
    axOverlay.set_title("Overlay")

    ax2IRT.set_title("True Hazard")
    ax2Haz.set_title("Fit Hazard")
    ax2Overlay.set_title("Overlay")
    
    
#     axIRT.legend(loc="bottom", bbox_to_anchor=(0, 1.3))
#     axHaz.legend(loc="bottom", bbox_to_anchor=(0, 1.3))
#     axOverlay.legend(loc="bottom", bbox_to_anchor=(0, 1.3))
    axIRT.set_xticks(0:1:7)
    axHaz.set_xticks(0:1:7)
    axOverlay.set_xticks(0:1:7)
    axIRT.set_xlim([0,7])
    axHaz.set_xlim([0,7])
    axOverlay.set_xlim([0,7])
    axIRT.set_xlabel("time (s)") 
    axHaz.set_xlabel("time (s)") 
    axOverlay.set_xlabel("time (s)") 
    axIRT.set_ylabel("normalized hazard") 

    ax2IRT.set_xticks(0:1:7)
    ax2Haz.set_xticks(0:1:7)
    ax2Overlay.set_xticks(0:1:7)
    ax2IRT.set_xlim([0,7])
    ax2Haz.set_xlim([0,7])
    ax2Overlay.set_xlim([0,7])
    ax2IRT.set_xlabel("time (s)") 
    ax2Haz.set_xlabel("time (s)") 
    ax2Overlay.set_xlabel("time (s)") 
    ax2IRT.set_ylabel("normalized hazard") 
    
    # diagnostics:
    println("Good fits included: ")
    print("     ")
    pretty_print_list(seshCodes[goodidx], orient="horizontal")
    println(" ")
    println("Bad fits excluded: ")
    print("     ")
    pretty_print_list(seshCodes[failidx], orient="horizontal")
    println(" ")
     
    printFigure(join(["haz_results_composite_", figname]); fig=fig,figurePath=figurePath, suptitle=true, h_suptitle=h)
    printFigure(join(["haz_results_CI_composite_", figname]); fig=fig2,figurePath=figurePath, suptitle=true, h_suptitle=h)

    ret_dir = pwd()

    cd(figurePath)
    # save the new variables...
    writeMATLAB(allIRT, allIRT)
    writeMATLAB(allHaz, allHaz)
    writeMATLAB(meanIRT, meanIRT)
    writeMATLAB(meanHaz, meanHaz)
    writeMATLAB(xs, 1:length(xs))
    writeMATLAB(IRTbymouse, IRTbymouse)
    writeMATLAB(hazbymouse, hazbymouse)
    writeMATLAB(names, names)
    


    




    cd(ret_dir)
end
function composite_fit_hazard(collated_results; seshNos=[], p_prior=0., ndp_per_sample=50, pre_normalize_p=false, normalize=true, filename="", savepath="")
    
    
    hazs = []
    lts = []
    seshCodes = []
    for i = 1:length(seshNos)
        seshNo = seshNos[i]
        progressbar(i,length(seshNos))
        pprior_star = get_prior(collated_results, seshNo, :LickState)
#         println("pprior_star=", pprior_star)
        yID = collated_results.results[seshNo].modelData[1].yID[1]
        th = collated_results.results[seshNo].modelData[1].th_means[1]
        df = collated_results.results[seshNo].modelData[1].df[1]
        push!(seshCodes, df.SessionCode[1])
        X = df.X
        predictors = collated_results.results[seshNo].modelData[1].predictors # must be a vector
#         println(predictors)
        # Let's get the range of probabilities the model can do on the whole dataset...
        # println(yID)
        # println(predictors)
        all_p,_ = predict_logit(yID, predictors, th, df,plotOn=false)
        min_p = minimum(all_p)
        max_p = maximum(all_p)
#         println("min p_fit=", min_p)
#         println("max p_fit=", max_p)


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
                p_corrected = p_corrected*pprior_star

                p_by_trial[ii,j] = p_corrected
            end
        end
        haz = nanmean_mat(p_by_trial)
        lt = extract_behavior_distribution(df, verbose=false)
        push!(hazs, haz)
        push!(lts, lt)
    end
    print("Finished:   ")
    pretty_print_list(seshCodes, orient="horizontal")
    cc=pwd()
    cd(savepath)
    println("saving... \n")

    try
        writeMATLAB(hazs, hazs)
        writeMATLAB(lts, lts)
        writeMATLAB(seshCodes, seshCodes)
    catch
        warning("Was not able to write matlab file...check this. Writing CSV instead")
        CSV.write(join([filename, "_hazs.csv"]), DataFrame(hazs=hazs))
        CSV.write(join([filename, "_lts.csv"]), DataFrame(lts=lts))
        CSV.write(join([filename, "_seshCodes.csv"]), DataFrame(seshCodes=seshCodes))
    end
    
    cd(cc)
    return (hazs,lts,seshCodes)
end
function bootCI_excludenan(Arr, nboot; a=0.05)
    # replicates are in DIM 1
    #
    #  sample x datapoints, e.g., 98 animals x 5000 ms hazard
    #
    n = size(Arr)[2]
    reps = size(Arr)[1]
    bootData = nanmat(nboot, n)

    # for iboot = 1:nboot
    idxs = sample(1:reps,n)
    println(idxs)
    bootData = Arr[idxs,:]
    # end
    ix_min = round(Int, nboot * a/2)
    ix_max = round(Int, nboot - nboot * a/2)
    sortedBoot = sort(bootData, dims=1)
    CImin = sortedBoot[ix_min, :]
    CImax = sortedBoot[ix_max, :]

    figure(figsize=(5,3))
    plot(nanmean_mat(Arr, 1),"k-", linewidth=3)
    plot(nanmean_mat(bootData, 1),"r-", linewidth=1)
    plot(CImin,"g-", linewidth=1)
    plot(CImax,"g-", linewidth=1)


end