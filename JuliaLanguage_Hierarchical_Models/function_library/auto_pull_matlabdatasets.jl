#  auto_pull_matlabdatasets.jl
# 
#   Puts together the tools needed to parse a collation folder for its datasets
# 
#. collatedPath = the path to the collated CSV folder for all the sessions to run
#.     This folder has a folder for each session, within
#			within this folder is singletrial, baseline, LOI so that we can import data for each file
#
using GLM

function run_collated_model(collatedPath::String, modelpackagefunction::Function; 
	runID=0, pathIDx = [], runFails=false, failDirs=[], 
	postprocessingfunction::Function=template_postprocessingfunction, 
	compositesavepath::String="", suppressFigures=false)
	start_dir = pwd()
	if isempty(compositesavepath)
	#
	# Make a results folder ABOVE our collated path
	#
		cd(collatedPath)
		cd("..")
		try
			mkdir(join(["Collated_results_", modelpackagefunction(""; runID=runID, getpackagename=true)]))
		catch
		end
		cd(join(["Collated_results_", modelpackagefunction(""; runID=runID, getpackagename=true)]))
		compositesavepath = pwd()
		cd(start_dir)
	end
	# Can rerun an ID with runID !=0
	if runID==0
		runID = rand(1:10000)
	end
	# Find all the sessions in the folder
	# pathIDx allows you to select only some of the sessions in the path. only use this for forward mode (NOT FAIL MODE)
	println("-----------------------------------------------")
    println("	")
    println("Initializing run_collated_model for ", modelpackagefunction(""; runID=runID, getpackagename=true))
    println("	")

    
    try
    	root = collatedPath
    	dirs = String[]
    	files = String[]
    	a = readdir(collatedPath)
    	for i in a
		    if i == ".DS_Store"
		        #skip
		    elseif isdir(joinpath(root, i))
		        push!(dirs, i)
		    elseif isfile(joinpath(root, i))
		        push!(files, i)
		    end
		end
    	# for (rooti, dirsi, filesi) in walkdir(collatedPath)
    	# 	root = rooti
    	# 	push!(dirs, dirsi)
    	# 	push!(files, filesi)
    	# end
        if !runFails
	        println("Found ", length(dirs), " sessions in $root:")
	        pretty_print_list(dirs, orient="vertical", enum=true)
	        if isempty(pathIDx)
	        	println("	Using all these sessions.")
	        	dirs = dirs
        	else
				println("	Using only a subset of these sessions:")
				dirs = dirs[pathIDx]
				pretty_print_list(dirs, orient="vertical", enum=false)
        	end
        else
        	println("Re-running failures in $root:")
        	pretty_print_list(failDirs, orient="vertical", enum=true)
        	dirs = failDirs	
        end
        sessionIDs = String[]
        sessionPaths = String[]
        sessionIdx = Vector{Int}()
        results = []#Vector{DataFrame}()
        failDirs = String[]
        result_df = DataFrame([Int,String, DataFrame, String], [:sessionIdx, :sessionIDs, :results, :sessionPaths])
        println("	")

	    for i = 1:length(dirs)
	    	dir = dirs[i]
	        # println(joinpath(root, dir)) # path to directories
	        try
	        	println("-----------------------------------------------")
	        	println("Processing ", dir, "...(", timestamp_now(), ")")
	        	print("   ")
	        	progressbar(i, length(dirs))
	        	println("	")
	        	push!(sessionIdx, i)
		        push!(sessionIDs, dir)
		        push!(sessionPaths, joinpath(root, dir))

		        result = modelpackagefunction(joinpath(root, dir), sessionID = dir, runID = runID, suppressFigures=suppressFigures)
		        push!(results, result)
		        
	        	result_df = DataFrame(sessionIdx=sessionIdx, sessionIDs=sessionIDs, results = results, sessionPaths=sessionPaths)
	        	# println(typeof(result_df.results))
	        	# println(typeof(result_df.results[1]))
	        	if suppressFigures
	        		close() # closes all the open figures
        		end
	        catch
	        	println("	!********* Encountered error! Skipping this directory")
	        	println("	")
	        	push!(failDirs, dir)
	        	rethrow()
	        end
	    end
	    println("-----------------------------------------------")
	    println("	")
	    println("Completed procesing of ", length(dirs) - length(failDirs), " sessions. (", timestamp_now(), ") ~")
	    println("	")
	    println("	Find results in each sessions' folder in: ", modelpackagefunction(""; runID=runID, getpackagename=true))
	    println(collatedPath)
	    println("	")
	    println("Initiating post-modeling collation of results...")
	    post_processing_results = postprocessingfunction(result_df, compositesavepath, modelpackagefunction; runID=runID)
	    println("Post-modeling collation of results complete and variables saved to:")
	    println(compositesavepath)
	    println("	")
	    println("-----------------------------------------------")
	    println("	")
	    cd(start_dir)
	    return failDirs, result_df, post_processing_results
    catch
        cd(start_dir)
        try
        	typeof(result_df)
    	catch
    		result_df = []
		end
        # return failDirs, result_df
        rethrow()
    end
end

function template_modelpackage(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	#
	# Use this to build a new analysis
	#
	# name the package and runID
	packagename = join(["templatepkg_",runID])
	if getpackagename
		return packagename
	end

	# do the business of the package on this session
	println("	good!")
	
	# Try to enter the results folder
	savepath = joinpath(path, join(["results_", packagename]))
	try 
		cd(savepath)
	catch
		mkdir(savepath)
		cd(savepath)
	end

	# Save each variable to our results folder
	CSV.write("result_good.csv", DataFrame(good = ["good", "good"]))
	CSV.write("result_othervar.csv", DataFrame(othervar = ["var1", "var2"]))

	# make a working result df with all the results to keep in workspace
	result = DataFrame(
		good = ["good", "good"], 
		othervar = ["var1", "var2"]
		)
	return result
end
function template_postprocessingfunction(results::DataFrame, compositesavepath, modelpackagefunction; runID=0)
	#
	# Use this to compile the analysis
	#
	# name the package and runID
	packagename = modelpackagefunction(""; sessionID ="", getpackagename=true, runID=runID)
	

	# do the business of the package on this session
	println("	This is a template postprocessingfunction. To do work, you need to 
		implement what you want to collate from the results here!") #"
	
	# Try to enter the composite results folder
	cd(compositesavepath)

	# Save each variable to our results folder
	CSV.write(join([packagename, "_result_good_",".csv"]), DataFrame(good = ["good", "good"]))
	CSV.write(join([packagename, "_result_othervar_",".csv"]), DataFrame(othervar = ["var1", "var2"]))
	return Nothing
end


function bootlogit_modelpackage1(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["bootlogit_modelpackage1_",runID])
	if getpackagename
		return packagename
	end
# Try to enter the results folder
	savepath = joinpath(path, join(["results_", packagename]))
	figurePath = joinpath(path, join(["figures_", packagename]))
	try 
		cd(savepath)
	catch
		mkdir(savepath)
		mkdir(figurePath)
		cd(savepath)
	end
# do the business of the package on this session
	#
	# first, we extract the relevant data: the singletrial, baseline and LOI sets and make a df
	# (expecting singletrial, baseline, and LOI folders for each dataset with CSV files from matlab)
	#
	ndf = extract_data_with_baselineandLOI(path; normalize=true)
	#
	# next, we want to specify our model formulae, including the nested models
	#
	formulas = [
	    @formula(LickState ~ Y),
	    @formula(LickState ~ Mean_Baseline),
	    @formula(LickState ~ Median_Baseline),
	    @formula(LickState ~ Mean_LOI),
	    @formula(LickState ~ Median_LOI),
	    @formula(LickState ~ Mean_Baseline + Mean_LOI + Y),
	    @formula(LickState ~ Median_Baseline + Median_LOI + Y),
	    @formula(LickState ~ LickTime_1back),
	    @formula(LickState ~ LickTime_2back),
	    @formula(LickState ~ LickTime_2back + LickTime_1back),
	    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back),
	    @formula(LickState ~ Rxn_1back + Early_1back + Reward_1back + ITI_1back),
	    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + 
	    	Rxn_1back + Early_1back + Reward_1back + ITI_1back),

	    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + 
	    	Rxn_1back + Early_1back + Reward_1back + ITI_1back + 
	    	Y),
	    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + 
	    	Rxn_1back + Early_1back + Reward_1back + ITI_1back +
	        Mean_Baseline + Mean_LOI + Y),
	    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + 
	    	Rxn_1back + Early_1back + Reward_1back + ITI_1back + 
	        Median_Baseline + Median_LOI + Y),

	    @formula(LickState ~ LickTime_2back + LickTime_1back + 
	    	Rxn_2back + Early_2back + Reward_2back + ITI_2back +
	    	Rxn_1back + Early_1back + Reward_1back + ITI_1back + 
	        Y),
	    @formula(LickState ~ LickTime_2back + LickTime_1back + 
	    	Rxn_2back + Early_2back + Reward_2back + ITI_2back +
	    	Rxn_1back + Early_1back + Reward_1back + ITI_1back + 
	        Mean_Baseline + Mean_LOI + Y),
	    @formula(LickState ~ LickTime_2back + LickTime_1back + 
	    	Rxn_2back + Early_2back + Reward_2back + ITI_2back + 
	    	Rxn_1back + Early_1back + Reward_1back + ITI_1back + 
	        Median_Baseline + Median_LOI + Y),
		]

	modelNames = [
	    "DA-only",
	    "μBl-only",
	    "medBl-only",
	    "μLOI-only",
	    "medLOI-only",
	    "DA-μBl-μLOI",
	    "DA-medBl-medLOI",
	    "Lt1b-only",
	    "Lt2b-only",
	    "Lt1b-Lt2b",
	    "oc1b-only",
	    "oc1b-oc2b",
	    "DA-oc1b-oc2b",
	    "DA-μBl-μLOI-oc1b-oc2b",
	    "DA-medBl-medLOI-oc1b-oc2b",
	    "DA-Lt1b-Lt2b-oc1b-oc2b",
	    "DA-μBl-μLOI-Lt1b-Lt2b-oc1b-oc2b",
	    "DA-medBl-medLOI-Lt1b-Lt2b-oc1b-oc2b",
	]

	
	results = modelSelectionByAICBICxval(ndf, :LickState, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# this is already handled by the modelSelectionByAICBICxval function

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end


function bootlogit_modelpackage2(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
#
# This one only runs time in session as predictor
#
# name the package and runID
	packagename = join(["bootlogit_modelpackage2_",runID])
	if getpackagename
		return packagename
	end
# Try to enter the results folder
	savepath = joinpath(path, join(["results_", packagename]))
	figurePath = joinpath(path, join(["figures_", packagename]))
	try 
		cd(savepath)
	catch
		mkdir(savepath)
		mkdir(figurePath)
		cd(savepath)
	end
# do the business of the package on this session
	#
	# first, we extract the relevant data: the singletrial, baseline and LOI sets and make a df
	# (expecting singletrial, baseline, and LOI folders for each dataset with CSV files from matlab)
	#
	ndf = extract_data_with_baselineandLOI(path; normalize=true)
	#
	# next, we want to specify our model formulae, including the nested models
	#
	formulas = [
	    @formula(LickState ~ X),
	    @formula(LickState ~ X + Y),
		]

	modelNames = [
	    "elapsedtime",
	    "elapsedtime+DA"
	]

	
	results = modelSelectionByAICBICxval(ndf, :LickState, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# this is already handled by the modelSelectionByAICBICxval function

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end

function extract_data_with_baselineandLOI(path; normalize=true)
	ret_dir = pwd()
	cd(path)
	cd("..")
	data_single_trial = extract_data(joinpath(path, "singletrial"), blmode=false, LOImode=false)
	lb = -0.15
	println("truncating data back: ",lb, "s")
	
	# data_single_trial.sessionCode = [sessionID for _=1:length(data_single_trial.sessionCode)]
	data_baseline_trial = extract_data(joinpath(path, "baseline"), blmode=true, LOImode=false)
	data_LOI_trial = extract_data(joinpath(path, "LOI"), blmode=false, LOImode=true)
	
	(data_single_trial,data_baseline_trial,data_LOI_trial) = truncate_at_cue_and_lickp250(data_single_trial; cbuffer_s=0., lbuffer_s=lb, bl_data=data_baseline_trial, loi_data=data_LOI_trial)
	# check_imported_data(data_single_trial; idx=11)
	# printFigure(join(["datacheck_",data_single_trial.sessionCode,"_"]); fig=gcf(), figurePath=pwd())
	# check_imported_data(data_baseline_trial; idx=11)
	# printFigure(join(["datacheck_bl_",data_baseline_trial.sessionCode,"_"]); fig=gcf(), figurePath=pwd())
	# check_imported_data(data_LOI_trial; idx=11)
	# printFigure(join(["datacheck_loi_",data_LOI_trial.sessionCode,"_"]); fig=gcf(), figurePath=pwd())
	# println("saved figs to ", pwd())
	df = makeSessionDataFrame(data_single_trial; normalize=normalize, includeBL_LOI=true, baseline_data=data_baseline_trial, LOI_data=data_LOI_trial)
	cd(ret_dir)
	return df
end
# function extract_data_with_baselineandLOI(path; normalize=true)
# 	data_single_trial = extract_data(joinpath(path, "singletrial"), blmode=false, LOImode=false)
# 	# data_single_trial.sessionCode = [sessionID for _=1:length(data_single_trial.sessionCode)]
# 	data_baseline_trial = extract_data(joinpath(path, "baseline"), blmode=true, LOImode=false)
# 	data_LOI_trial = extract_data(joinpath(path, "LOI"), blmode=false, LOImode=true)
# 	df = makeSessionDataFrame(data_single_trial; normalize=normalize, includeBL_LOI=true, baseline_data=data_baseline_trial, LOI_data=data_LOI_trial)
# 	return df
# end

function bootlogit_postprocessingfunction1(results::DataFrame, compositesavepath, modelpackagefunction; runID=0)
	#
	# Use this to compile the analysis
	#
	# name the package and runID
	packagename = modelpackagefunction(""; sessionID ="", getpackagename=true, runID=runID)
	

	# do the business of the package on this session
	combine_th_across_sessions(results, compositesavepath, runID, packagename)
	combine_AICBIC_across_sessions(results, compositesavepath, runID, packagename)
	
	# Try to enter the composite results folder
	# cd(compositesavepath)

	# Save each variable to our results folder
	# 		This is handled by the combine functions ABOVE
	# CSV.write(join([packagename, "_result_good_",".csv"]), DataFrame(good = ["good", "good"]))
	# CSV.write(join([packagename, "_result_othervar_",".csv"]), DataFrame(othervar = ["var1", "var2"]))
	return Nothing
end


#
# 	Collation
#
function combine_th_across_sessions(results, compositesavepath, runID,packagename)
    ret_dir = pwd()
    cd(compositesavepath)
    modelNames = unique(results.results[1].th_summary[1].modelName)
    n_sesh = length(results.results)
    i_sesh=1
    i_model=1
    i_iter=1
    n_models = length(results.results[i_sesh].AICs)
    n_iters = length(results.results[i_sesh].AICs[i_model])
    
    axs = []
    fs = []
    th_summary = DataFrame(modelName=[], composite_th=[], composite_se=[], composite_CImin=[], composite_CImax=[], composite_mdof = [])
    for model = 1:n_models

        result_df = DataFrame(train_dof = [results.results[1].dofs[model][1]], train_ths = [[results.results[1].ths[model][1]]], train_se_ths=[[results.results[1].se_ths[model][1]]])
        for i_sesh = 1:n_sesh
        	n_iters = length(results.results[i_sesh].AICs[model])
      #   	println("	model:", model, " session: ", i_sesh)
	    	# println("n_iters=", n_iters)
		    # println("size BIC=", length(results.results[i_sesh].BICs[model]))
		    # println("size dofs=", length(results.results[i_sesh].dofs[model]))
		    # println("size ths=", length(results.results[i_sesh].ths[model]))
		    # println("size se_ths=", length(results.results[i_sesh].se_ths[model]))
            if i_sesh==1
                iterstart=2
            else
                iterstart=1
            end

            for i=iterstart:n_iters
                append!(result_df, 
                DataFrame(train_dof = [results.results[i_sesh].dofs[model][i]], 
                        train_ths = [[results.results[i_sesh].ths[model][i]]], 
                        train_se_ths=[[results.results[i_sesh].se_ths[model][i]]],
                    ))   
            end
        end
        
        (composite_th, composite_se, composite_CImin, composite_CImax, ax, f, composite_mdof) = theta_summary(result_df; Mode = "sparseFit", result_df=result_df)
        title(modelNames[model])
        push!(axs, ax)
        push!(fs, f)

        append!(th_summary, 
            DataFrame(
                modelName = modelNames[model],
                composite_th=composite_th, 
                composite_se=composite_se, 
                composite_CImin=composite_CImin, 
                composite_CImax=composite_CImax,
                composite_mdof = composite_mdof,
                ))
        
        #
        # Save the composite variables
        #
        # CSV.write(join(["MODELno",model, "_composite_ths_nboot", n_iters, "_nsesh", n_sesh, 
        #             "_", packagename, ".csv"]),DataFrame(train_ths = result_df.train_ths))
        # CSV.write(join(["MODELno",model, "_composite_se_ths_nboot", n_iters, "_nsesh", n_sesh, 
        #             "_", packagename, ".csv"]),DataFrame(train_se_ths = result_df.train_se_ths))
        # CSV.write(join(["MODELno",model, "_composite_dofs_nboot", n_iters, "_nsesh", n_sesh, 
        #             "_", packagename, ".csv"]),DataFrame(train_dof = result_df.train_dof))
        # CSV.write(join(["MODELno",model, "_composite_th_summary_nboot", n_iters, "_nsesh", n_sesh, 
        #             "_", packagename, ".csv"]),DataFrame(th_summary = th_summary))
        # println(pwd())
        # CSV.write("test.csv",DataFrame(train_ths = result_df.train_ths))

        CSV.write(join([modelNames[model], "_composite_ths_nboot", n_iters, "_nsesh", n_sesh, 
                    ".csv"]),DataFrame(train_ths = result_df.train_ths))
        CSV.write(join([modelNames[model], "_composite_se_ths_nboot", n_iters, "_nsesh", n_sesh, 
                    ".csv"]),DataFrame(train_se_ths = result_df.train_se_ths))
        CSV.write(join([modelNames[model], "_composite_dofs_nboot", n_iters, "_nsesh", n_sesh, 
                    ".csv"]),DataFrame(train_dof = result_df.train_dof))
        CSV.write(join([modelNames[model], "_composite_th_summary_nboot", n_iters, "_nsesh", n_sesh, 
                    ".csv"]),DataFrame(th_summary = th_summary))
    end

    set_xaxes_same_scale(axs)
    set_yaxes_same_scale(axs)
    println("Figs saved to:", compositesavepath)
    cd(compositesavepath)
    for i=1:length(fs)
        # println(i)
        # printFigure(join(["composite_", modelNames[i], "_theta_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=fs[i],figurePath=compositesavepath)

        printFigure(join(["composite_", modelNames[i], "_theta_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=fs[i],figurePath=compositesavepath)
    end
    cd(ret_dir)  
    return th_summary  
end


function combine_AICBIC_across_sessions(results, compositesavepath, runID, packagename)
    n_sesh = length(results.results)
    i_sesh=1
    i_model=1
    i_iter=1
    n_models = length(results.results[i_sesh].AICs)
    n_iters = length(results.results[i_sesh].AICs[i_model])
    allAICs = [[] for _=1:n_models]
    allAICcs = [[] for _=1:n_models]
    allBICs = [[] for _=1:n_models]
    all_Sn_accuracy = [[] for _=1:n_models]
    all_test_accuracy = [[] for _=1:n_models]
    
    for i_sesh = 1:n_sesh
        for i_model = 1:n_models
        	n_iters = length(results.results[i_sesh].AICs[i_model])
      #   	println("	model:", i_model, " session: ", i_sesh)
	    	# println("n_iters=", n_iters)
		    # println("size BIC=", length(results.results[i_sesh].BICs[i_model]))
		    # println("size dofs=", length(results.results[i_sesh].dofs[i_model]))
		    # println("size ths=", length(results.results[i_sesh].ths[i_model]))
		    # println("size se_ths=", length(results.results[i_sesh].se_ths[i_model]))
            AICs = [results.results[i_sesh].AICs[i_model][i_iter] for i_iter = 1:n_iters]
            append!(allAICs[i_model], AICs)
            AICcs = [results.results[i_sesh].AICcs[i_model][i_iter] for i_iter = 1:n_iters]
            append!(allAICcs[i_model], AICcs)
            BICs = [results.results[i_sesh].BICs[i_model][i_iter] for i_iter = 1:n_iters]
            append!(allBICs[i_model], BICs)
            Sn_accuracies = [results.results[i_sesh].Sn_accuracy[i_model][i_iter] for i_iter = 1:n_iters]
            append!(all_Sn_accuracy[i_model], Sn_accuracies)
            test_accuracies = [results.results[i_sesh].test_accuracy[i_model][i_iter] for i_iter = 1:n_iters]
            append!(all_test_accuracy[i_model], test_accuracies)
        end
    end
    mean_all_AICs = [mean(allAICs[i]) for i=1:n_models]
    mean_all_AICcs = [mean(allAICcs[i]) for i=1:n_models]
    mean_all_BICs = [mean(allBICs[i]) for i=1:n_models]
    mean_all_Sn_accuracy = [mean(all_Sn_accuracy[i]) for i=1:n_models]
    mean_all_test_accuracy = [mean(all_test_accuracy[i]) for i=1:n_models]
    f = figure(figsize=(20,3))
    ax1=subplot(1,3,1)
    compare_AICBIC(mean_all_AICs, allAICs; yl=join(["AIC nsesh=",n_sesh]), iters=n_iters, ax=ax1, minmax="min")
    ax2=subplot(1,3,2)
    compare_AICBIC(mean_all_AICcs, allAICcs; yl=join(["AICc nsesh=",n_sesh]), iters=n_iters, ax=ax2, minmax="min")
    ax3=subplot(1,3,3)
    compare_AICBIC(mean_all_BICs, allBICs; yl=join(["BIC nsesh=",n_sesh]), iters=n_iters, ax=ax3, minmax="min")
    # printFigure(join(["compositeAICBIC_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=f, figurePath=compositesavepath)
    printFigure(join(["compositeAICBIC_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=f, figurePath=compositesavepath)
    
    f = figure(figsize=(20,3))
    ax1=subplot(1,3,1)
    compare_AICBIC(mean_all_Sn_accuracy, all_Sn_accuracy; yl=join(["Train Accuracy nsesh=",n_sesh]), iters=n_iters, ax=ax1, minmax="max")
    ax1.set_ylim([0., 1.])
    ax2=subplot(1,3,2)
    compare_AICBIC(mean_all_test_accuracy, all_test_accuracy; yl=join(["Test Accuracy nsesh=",n_sesh]), iters=n_iters, ax=ax2, minmax="max")
    ax2.set_ylim([0., 1.])
    ax3=subplot(1,3,3)
    # printFigure(join(["composite_Accuracy_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=f, figurePath=compositesavepath)
    printFigure(join(["composite_Accuracy_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=f, figurePath=compositesavepath)
    #
    # Save the variables to the composite folder
    #
    ret_dir = pwd()
    cd(compositesavepath)

    # CSV.write(join(["composite_AICs_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(allAICs = allAICs))
    # CSV.write(join(["composite_meanAIC_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(mean_all_AICs = mean_all_AICs))
    # CSV.write(join(["composite_AICcs_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(allAICcs = allAICcs))
    # CSV.write(join(["composite_meanAICc_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(mean_all_AICcs = mean_all_AICcs))
    # CSV.write(join(["composite_BICs_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(allBICs = allBICs))
    # CSV.write(join(["composite_meanBIC_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(mean_all_BICs = mean_all_BICs))
    
    # CSV.write(join(["composite_Sn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(all_Sn_accuracy = all_Sn_accuracy))
    # CSV.write(join(["composite_meanSn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(mean_all_Sn_accuracy = mean_all_Sn_accuracy))
    # CSV.write(join(["composite_test_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(all_test_accuracy = all_test_accuracy))
    # CSV.write(join(["composite_meantest_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
    #             "_", packagename, ".csv"]),DataFrame(mean_all_test_accuracy = mean_all_test_accuracy))
    CSV.write(join(["composite_AICs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(allAICs = allAICs))
    CSV.write(join(["composite_meanAIC_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(mean_all_AICs = mean_all_AICs))
    CSV.write(join(["composite_AICcs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(allAICcs = allAICcs))
    CSV.write(join(["composite_meanAICc_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(mean_all_AICcs = mean_all_AICcs))
    CSV.write(join(["composite_BICs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(allBICs = allBICs))
    CSV.write(join(["composite_meanBIC_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(mean_all_BICs = mean_all_BICs))
    
    CSV.write(join(["composite_Sn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(all_Sn_accuracy = all_Sn_accuracy))
    CSV.write(join(["composite_meanSn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(mean_all_Sn_accuracy = mean_all_Sn_accuracy))
    CSV.write(join(["composite_test_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(all_test_accuracy = all_test_accuracy))
    CSV.write(join(["composite_meantest_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(mean_all_test_accuracy = mean_all_test_accuracy))
    
    cd(ret_dir)
end


function bootlogit_timeslice_modelpackage1(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["bootlogit_timeslice_modelpackage1_",runID])
	if getpackagename
		return packagename
	end
# Try to enter the results folder
	savepath = joinpath(path, join(["results_", packagename]))
	figurePath = joinpath(path, join(["figures_", packagename]))
	try 
		cd(savepath)
	catch
		mkdir(savepath)
		mkdir(figurePath)
		cd(savepath)
	end
# do the business of the package on this session
	#
	# first, we extract the relevant data: the singletrial, baseline and LOI sets and make a df
	# (expecting singletrial, baseline, and LOI folders for each dataset with CSV files from matlab)
	#
	ndf = extract_data_with_baselineandLOI(path; normalize=true)
	#
	# next, we need to specify binning for our timeslices
	#. as a start to match the hazard analysis, I'll use 250 ms bins
	#
	slice_width_ms = 250.#3000.#250.
	println("slice_width_ms: ", slice_width_ms)
	(binned_ndfs, binEdges) = slice_dataframe_into_timebins(ndf, slice_width_ms)
	# name the timeslices
	timeslice_names = []
	for i = 1:length(binned_ndfs)
		push!(timeslice_names, join(["_", round(Int,1000*binEdges[i]), "-", round(Int,1000*binEdges[i+1])]))
	end
	#
	# next, we want to specify our model formulae, including the nested models
	#
	formulas = [
	    @formula(LickState ~ Y),
	    @formula(LickState ~ Mean_Baseline),
	    @formula(LickState ~ Median_Baseline),
	    @formula(LickState ~ Mean_LOI),
	    @formula(LickState ~ Median_LOI),
	    @formula(LickState ~ Mean_Baseline + Mean_LOI + Y),
	    @formula(LickState ~ Median_Baseline + Median_LOI + Y),
		]

	modelNames = [
	    "DA-only",
	    "μBl-only",
	    "medBl-only",
	    "μLOI-only",
	    "medLOI-only",
	    "DA-μBl-μLOI",
	    "DA-medBl-medLOI",
	]
	nbins = length(binEdges) - 1
	results = Vector{DataFrame}()
	for slice = 1:nbins-1
		progressbar(slice,nbins)
		println("----------slice: ", timeslice_names[slice])
		#
		# Specify the timeslice
		#
		modelNames_slice = [join(hcat(modelNames[x], timeslice_names[slice])) for x=1:length(modelNames)]
		result = modelSelectionByAICBICxval(binned_ndfs[slice], :LickState, formulas, modelNames_slice, "logit"; 
	    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures, slice=timeslice_names[slice])
		result[:TimeSlice] = [timeslice_names[slice] for _ = 1:nrow(result)]
		push!(results, result)
	end

	

# make a working result list of dfs with all the results to keep in workspace
	return results#, ndf
end
function bootlogit_timeslice_postprocessingfunction1(results::DataFrame, compositesavepath, modelpackagefunction; runID=0)
	#
	# Use this to compile the analysis when we have a LIST of result dfs
	#
	# name the package and runID
	packagename = modelpackagefunction(""; sessionID ="", getpackagename=true, runID=runID)
	

	# do the business of the package on this session
	# 
	# We need to create a working list for each model from each session
	ret_dir = pwd()
	cd(compositesavepath)
	sdir = join(["Composite_", packagename])
	try 
		mkdir(sdir)
	catch
	end
	cd(sdir)
	sdd = pwd()

	#
	# results_df[session (df)].results[timeslice# (array)].TimeSlice
	#
	# results.results[session][timeslice][:COLUMN(e.g., TimeSlice)][modlNo]

	# we want to gather the results across all sessions at same timeslice...
	
	# println("ntimeslices=", nrow(results.results[1][1]))
	# println("n-sesh=", nrow(results))
	by_slice_dfs = []
	by_slice_composite_ths = []
	for tslice = 1:length(results.results[1])#.results) # for each timeslice
		slicepath = results.results[1][tslice][:TimeSlice][1] # the timeslice name from session 1 model 1
		try 
			mkdir(join([runID,slicepath]))
		catch
		end
		cd(join([runID,slicepath]))
		sav_dir = pwd()

		# println("saving model results for this timeslice across sessions to : ",pwd())
		# normally, we have a results_df = DataFrame(sessionIdx=sessionIdx, sessionIDs=sessionIDs, results = results, sessionPaths=sessionPaths)
		results_df_this_slice = DataFrame(sessionIdx = results[1,:sessionIdx], sessionIDs = results[1,:sessionIdx], results = results.results[1][tslice], sessionPaths = results[1,:sessionPaths])
		if nrow(results) > 1
			for sesh = 2:nrow(results)
				rr = DataFrame(sessionIdx = results[sesh,:sessionIdx], sessionIDs = results[sesh,:sessionIdx], results = results.results[sesh][tslice], sessionPaths = results[sesh,:sessionPaths])
				append!(results_df_this_slice,rr)
			end
		end
		push!(by_slice_dfs,results_df_this_slice)
		# println("		")
		# println("	tslice=", tslice)
		th_summary = combine_th_across_sessions(results_df_this_slice, sav_dir, join([runID,slicepath]), packagename)
		combine_AICBIC_across_sessions(results_df_this_slice, sav_dir, join([runID,slicepath]), packagename)
		cd("..")
		th_summary[:sliceID] = [slicepath for _=1:nrow(th_summary)]
		push!(by_slice_composite_ths,th_summary)
	end
	plot_th_vs_timeslice(by_slice_composite_ths,savedir=sdd)
	cd(ret_dir)
	return by_slice_dfs,by_slice_composite_ths
end

function slice_dataframe_into_timebins(df::DataFrame, slice_width_ms::Float64=250.)
	#
	# We have a dataframe whose :X column is the time in the interval.
	#   We want to divide this into timebins
	#
	# println("Unique of df LickState:", unique(df[:LickState]))
	nbins = ceil(Int, 7000. /slice_width_ms)+1
	#
	# Now we want to gather the bin edges
	#
	binEdges = collect(range(0., length=nbins, step=slice_width_ms/1000))
	#
	# For each left bin edge, find all the indicies from df that go in that bin and make a new df
	#
	binned_dfs = []
	for i = 1:nbins-1
	    # get indicies
	    idx = findall(x-> x >= binEdges[i] && x < binEdges[i+1], df.X)
	    # get all the unique trials from these indicies and get an average of the fields of interest
	    trials = unique(df[idx,:TrialNo])
	    # for each trial, get the mean of the fields...
	    newData = similar(df,0)
	    for t = 1:length(trials)
	        # find the idxs that match the trial
	        tidx = idx[findall(x->x==trials[t], df[idx, :TrialNo])]
	        global __TrialData = DataFrame()
	        coldata = []#[[] for _= 1:ncol(df)]
	        for col = 1:ncol(df)
	            if names(df)[col] == "LickState"
	            	# what does this actually do? We want to preseve a bool! What the heck?
	                # cc = maximum(df[tidx, col])
	                # let's instead return a bool
	                # println("Number of lick states this trial: ", length(findall(x->x, df[tidx, col])))
	                cc = length(findall(x->x, df[tidx, col]))
	            elseif typeof(df[!,col][1])<:Number
	                cc = mean(df[tidx, col])
	            else 
	                cc = df[tidx[1], col]
	            end
	            ss = names(df)[col]
	            if typeof(cc) == String
	            	if Sys.iswindows()
	            		cc = join(split(cc, "\\"), "\\\\")
	            		eval(Meta.parse(join(["__TrialData.",ss,"=[\"", cc ,"\"]"])))
	            	else
	                	eval(Meta.parse(join(["__TrialData.",ss,"=[\"", cc ,"\"]"])))
                	end
	            else
	                eval(Meta.parse(join(["__TrialData.",ss,"=[",cc,"]"])))
	            end
	            # println(__TrialData)
	            
	        end
	        newData = vcat(newData, __TrialData)
	        # println("Bin No:", i, "Unique of df LickState:", unique(newData[:LickState]))
	    end
	    push!(binned_dfs, newData)
	end
	# for i = 1:nbins-1
	# 	# get indicies
	# 	idx = findall(x-> x >= binEdges[i] && x < binEdges[i+1], df.X)
	# 	# get all the unique trials from these indicies and get an average of the fields of interest
	# 	trials = unique(df[idx,:TrialNo])
	# 	# for each trial, get the mean of the fields...
	# 	newData = similar(df,0)
	# 	for t = 1:length(trials)
	# 		# find the idxs that match the trial
	# 		tidx = idx[findall(df[idx, :TrialNo])]
	# 		trialData = similar(df,0)
	# 		coldata = []#[[] for _= 1:ncol(df)]
	# 		for col = 1:ncol(df)
	# 			if names(df)[col] == :LickState
	# 				cc = maximum(df[tidx, col])
	# 			elseif typeof(ndf[!,col][1])<:Number
	# 				cc = mean(df[tidx, col])
	# 			else 
	# 				cc = df[tidx[1], col]
	# 			end
	# 			push(coldata, cc)
	# 		end
	# 		push!(newData, coldata)
	# 	end

	#     # add the new data row by row
	#     # for j = 1:length(idx)
	#     #     push!(newData, df[idx[j],:])
	#     # end
	#     push!(binned_dfs, newData)
	# end


	return (binned_dfs, binEdges)
end


function plot_th_vs_timeslice(by_slice_composite_ths; savedir=pwd())
    ret_dir=pwd()
    cd(savedir)
    try
        mkdir(join(["th_vs_timeslice", timestamp_now()]))
    catch
    end
    cd(join(["th_vs_timeslice", timestamp_now()]))
   # by_slice_composite_ths is an array of th_summary's for each timeslice
    #
    #  extract the thetas for each model
    #             by_slice_dfs2[iTimeSlice][:results][sesh][:dofs][mdl]
    #   Get model idxs
    nTimeSlices = length(by_slice_composite_ths)
    nmodels = length(unique(by_slice_composite_ths[1][:,:modelName]))
    modelIdxs = [findall(x->x==ii, by_slice_composite_ths[1][:,:modelName]) for ii in unique(by_slice_composite_ths[1][:,:modelName])]
    ths = [[] for _=1:nmodels]
    se_ths = [[] for _=1:nmodels]
    mdofs = [[] for _=1:nmodels]
    CImin = [[] for _=1:nmodels]
    CImax = [[] for _=1:nmodels]
    sliceID = [[] for _=1:nmodels]
    composite_th_summary = [[] for _=1:nmodels]
    ax = []
    axs2 = []
    for imodel = 1:nmodels
        for iTimeSlice = 1:nTimeSlices
            push!(ths[imodel], by_slice_composite_ths[iTimeSlice][:composite_th][modelIdxs[imodel]])
            push!(se_ths[imodel], by_slice_composite_ths[iTimeSlice][:composite_se][modelIdxs[imodel]])
            push!(CImin[imodel], by_slice_composite_ths[iTimeSlice][:composite_CImin][modelIdxs[imodel]])
            push!(CImax[imodel], by_slice_composite_ths[iTimeSlice][:composite_CImax][modelIdxs[imodel]])
            push!(mdofs[imodel], by_slice_composite_ths[iTimeSlice][:composite_mdof][modelIdxs[imodel]])
            push!(sliceID[imodel], by_slice_composite_ths[iTimeSlice][:sliceID][modelIdxs[imodel][1]])
        end
        
        d_model = length(ths[imodel][1])
        f = figure(figsize=(d_model*5,2))
        my_suptitle=suptitle(by_slice_composite_ths[1][:modelName][modelIdxs[imodel][1]], y=1.2)
        composite_th_summaries = DataFrame(thID=[], th = [], se_th=[], CImin=[], CImax=[], mdof=[], nslices=[])
        for iTh = 1:d_model

            subplot(1,d_model,iTh)
            thi = [ths[imodel][i_timeslice][iTh] for i_timeslice=1:nTimeSlices]
            sliceIDi = [sliceID[imodel][i_timeslice] for i_timeslice=1:nTimeSlices]
            se_th_i = [se_ths[imodel][i_timeslice][iTh] for i_timeslice=1:nTimeSlices]
            CImin_i = [CImin[imodel][i_timeslice][iTh] for i_timeslice=1:nTimeSlices]
            CImax_i = [CImax[imodel][i_timeslice][iTh] for i_timeslice=1:nTimeSlices]
            dof_i = [mdofs[imodel][i_timeslice][iTh] for i_timeslice=1:nTimeSlices]
            plot(0:nTimeSlices+1, zeros(length(0:nTimeSlices+1)), "k-")    
            plot(1:nTimeSlices, thi, "r.-")
            plot_with_CI(thi, CImin_i, CImax_i; dataLabels="", ylab="wt", ax=gca())
            gca().set_xticks(1:nTimeSlices)
            xlabel("Timeslice")
            if !isempty(findall(x->!isnan(x), thi))
	            mn = maximum([minimum(thi[findall(x->!isnan(x), thi)]), -10]) - 1
	            mx = minimum([maximum(thi[findall(x->!isnan(x), thi)]), 10]) + 1
            else
            	mn = -10
            	mx = 10
            end
            ylim([mn,mx])
            if iTh==1
                ylabel("weight")
            end
            xticks(rotation = 45)
            gca().set_xticklabels(sliceIDi)
            title(join(["model: ", imodel, " th", iTh-1]))
            push!(ax, gca())
            this_model_this_th = DataFrame(thi=thi, sliceIDi=sliceIDi, CImin_i=CImin_i, CImax_i=CImax_i)
            this_model_this_th[:iTh] = [iTh for _=1:nrow(this_model_this_th)]
            this_model_this_th[:model] = [imodel for _=1:nrow(this_model_this_th)]
            CSV.write(join(["thbytimeslice_df_th", iTh-1, "_model_", imodel, ".csv"]),this_model_this_th)
            #             
            #  Get a composite th across slices
            #             
            #    Don't include bad models!!
            #
            goodidx = findall(x->abs(x)<10 && !isnan(x), CImin_i)
            (meanTh_c, propagated_se_th_c, CImin_c, CImax_c, mdf_c) = getCompositeTheta(thi[goodidx], se_th_i[goodidx], dof_i[goodidx])
            composite_th_summary = DataFrame(thID=join(["th", iTh]), th = meanTh_c, se_th=propagated_se_th_c, CImin=CImin_c, CImax=CImax_c, mdof=mdf_c, nslices=length(goodidx))
            append!(composite_th_summaries,composite_th_summary)
            
        end
        figname = join(["thbytimeslice_Model_", imodel, "_", timestamp_now(), ".eps"])
        f.savefig(figname, transparent=true, format="eps", bbox_inches="tight",bbox_extra_artists=[my_suptitle])
        
        
        #
        # Don't use any bad models, ie where SE is > 10
        #
        thiii = composite_th_summaries.th
        CIminiii = composite_th_summaries.CImin
        CImaxiii = composite_th_summaries.CImax
        CSV.write(join(["th_across_timeslices_df_th_summary__model_", imodel, ".csv"]),composite_th_summaries)
        f2 = figure(figsize=(d_model,2))
        plot(1:d_model, zeros(length(1:d_model)), "k-")    
        plot_with_CI(thiii, CIminiii, CImaxiii; dataLabels="", ylab="wt", ax=gca())
        gca().set_xticks(1:d_model)
        xlabel("th")
        if !isempty(findall(x->!isnan(x), thiii))
	        mn = maximum([minimum(thiii[findall(x->!isnan(x), thiii)]), -10]) - 1
    	    mx = minimum([maximum(thiii[findall(x->!isnan(x), thiii)]), 10]) + 1
	    else
	    	mn = -10
	    	mx = 10
	    end
        ylim([mn,mx])
        ylabel("weight")
        title(join(["model: ", imodel, " nslices=", composite_th_summaries.nslices]))
        push!(axs2, gca())
        figname = join(["th_across_timeslices_Model_", imodel, "_", timestamp_now(), ".eps"])
        my_suptitle=suptitle(by_slice_composite_ths[1][:modelName][modelIdxs[imodel][1]], y=1.2)
        f2.savefig(figname, transparent=true, format="eps", bbox_inches="tight",bbox_extra_artists=[my_suptitle])
    end
    set_yaxes_same_scale(ax)
    set_yaxes_same_scale(axs2)
    cd(ret_dir)
end
