function run_hierarchy_on_session(inference_function, trial_data::TrialData; ntraces_per_trial::Int=20, amount_of_computation::Int=50, figpath=pwd(),generate_hierarchical_model_fxn=generate_hierarchical_model)
    sessionCode = trial_data.sessionCode
    modelID = "(percent_slope, traces2) =data_driven_probability_of_model_class(generate_hierarchical_model_fxn, 
            (xx,1), yy, ransac_assisted_model_selection_proposal, 
            (xx,yy); amount_of_computation=amount_of_computation, ntraces=ntraces)"
    p = []
    traces = []
    slopes = []
    intercepts = []
    ls = []
    rs = []
    sps = []
    flags = []

    function get_t(t)
        println("starting trial ", t)
    end
    
    for t = 1:length(trial_data.xdata)
#         println(t)
        get_t(t)

        progressbar(t,length(trial_data.xdata))
        tNo = trial_data.trialNo[t]
        lt = trial_data.lickTime_s[t]
        xx = trial_data.xdata[t]
        yy = trial_data.ydata[t]
#         println(size(xx))

        percent_slope=NaN
        tr=[];
        flgs=[];
        try
            (percent_slope, tr, flgs) =inference_function(generate_hierarchical_model_fxn, 
                    (xx,1), yy, ransac_assisted_model_selection_proposal, 
                    (xx,yy); amount_of_computation=amount_of_computation, ntraces=ntraces_per_trial)
            println("finished inference_function")
        catch
            warning(join(["failed to complete inference function on trial ", t]))
            badnews(join([e.msg, "!"]))
            continue
        end
        # # Plotting takes a long time to render. So let's implement tools to save the figures but not leave them open.

        f = figure(figsize=(3,3))
        overlay(render_trace,tr[1:end], suppressAlpha=true)
        ax = gca()
        ax.set_title(join(["t=",tNo, " lt=", lt, " p(slope)=", percent_slope]))
        name = join([sessionCode,"_t",tNo])

        println("finished figure")
        
        printFigure(name; fig=f, figurePath=figpath, verbose=false, suptitle=false, h_suptitle=[])
        close(f)

        println("printed figure")


        push!(p,percent_slope)
        push!(traces, tr)
        """Here, we will identify the fit slopes and other params from each trace"""
        slope = []
        intercept = []
        l = []
        r = []
        sp = []
        for itrace = 1:length(tr)
            if tr[itrace][:model_choice_is_slope]
                push!(slope,tr[itrace][:slope])
                push!(intercept,tr[itrace][:intercept])
            else
                push!(l,tr[itrace][:left_segment_amp])
                push!(r,tr[itrace][:right_segment_amp])
                push!(sp,tr[itrace][:step_position])
            end
        end
        push!(slopes, slope)
        push!(intercepts, intercept)
        push!(ls, l)
        push!(rs, r)
        push!(sps, sp)
        push!(flags, flgs)
        println("reached end of loop")
    end

    results = MR4(sessionCode, modelID, amount_of_computation, ntraces_per_trial, 
        p,traces,slopes,intercepts,ls,rs,sps, flags,trial_data.trialNo,trial_data.lickTime_s)
    return results
end;



function hierarchy_v1(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
    #
    # Use this to build a new analysis
    #
    # name the package and runID
    packagename = join(["hierarchy_v1_",runID])
    if getpackagename
        return packagename
    end

    # do the business of the package on this session
    data_single_trial = extract_data(joinpath(path, "singletrial"), blmode=false, LOImode=false)
    lb = -0.15
    println("truncating data back: ",lb, "s")
    data_baseline_trial = extract_data(joinpath(path, "baseline"), blmode=true, LOImode=false);
    data_LOI_trial = extract_data(joinpath(path, "LOI"), blmode=false, LOImode=true);
    (data_single_trial,data_baseline_trial,data_LOI_trial) = truncate_at_cue_and_lickp250(data_single_trial; cbuffer_s=0.5, lbuffer_s=lb, bl_data=data_baseline_trial, loi_data=data_LOI_trial);
    data = data_single_trial;


    # Try to enter the results folder
    # savepath = joinpath(path, join(["results_", packagename]))
    savepath = joinpath(path, packagename)
    try 
        cd(savepath)
    catch
        mkdir(savepath)
        cd(savepath)
    end
    data_fullmodel = run_hierarchy_on_session(data_driven_probability_of_model_class, 
        data; 
        ntraces_per_trial=20, 
        amount_of_computation=50, 
        figpath=savepath, 
        generate_hierarchical_model_fxn=generate_hierarchical_model);
    saveMR4asMAT(data_fullmodel, savepath)


    # Save each variable to our results folder
    

    # make a working result df with all the results to keep in workspace
    result = DataFrame(savepath=savepath)
    return result
end
