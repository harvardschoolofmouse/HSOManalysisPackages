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
        results = Vector{DataFrame}()
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
	        	if suppressFigures
	        		close() # closes all the open figures
        		end
	        catch
	        	println("	!********* Encountered error! Skipping this directory")
	        	println("	")
	        	push!(failDirs, dir)
	        	# rethrow()
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
	    postprocessingfunction(result_df, compositesavepath, modelpackagefunction; runID=runID)
	    println("Post-modeling collation of results complete and variables saved to:")
	    println(compositesavepath)
	    println("	")
	    println("-----------------------------------------------")
	    println("	")
	    cd(start_dir)
	    return failDirs, result_df
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

function extract_data_with_baselineandLOI(path; normalize=true)
	data_single_trial = extract_data(joinpath(path, "singletrial"), blmode=false, LOImode=false)
	# data_single_trial.sessionCode = [sessionID for _=1:length(data_single_trial.sessionCode)]
	data_baseline_trial = extract_data(joinpath(path, "baseline"), blmode=true, LOImode=false)
	data_LOI_trial = extract_data(joinpath(path, "LOI"), blmode=false, LOImode=true)
	df = makeSessionDataFrame(data_single_trial; normalize=normalize, includeBL_LOI=true, baseline_data=data_baseline_trial, LOI_data=data_LOI_trial)
	return df
end

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
    th_summary = DataFrame(modelName=[], composite_th=[], composite_se=[], composite_CImin=[], composite_CImax=[])
    for model = 1:n_models
        result_df = DataFrame(train_dof = [results.results[1].dofs[model][1]], train_ths = [[results.results[1].ths[model][1]]], train_se_ths=[[results.results[1].se_ths[model][1]]])
        for i_sesh = 1:n_sesh
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
        
        (composite_th, composite_se, composite_CImin, composite_CImax, ax, f) = theta_summary(result_df; Mode = "sparseFit", result_df=result_df)
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
                ))
        
        #
        # Save the composite variables
        #
        CSV.write(join(["MODELno",model, "_composite_ths_nboot", n_iters, "_nsesh", n_sesh, 
                    "_", packagename, ".csv"]),DataFrame(train_ths = result_df.train_ths))
        CSV.write(join(["MODELno",model, "_composite_se_ths_nboot", n_iters, "_nsesh", n_sesh, 
                    "_", packagename, ".csv"]),DataFrame(train_se_ths = result_df.train_se_ths))
        CSV.write(join(["MODELno",model, "_composite_dofs_nboot", n_iters, "_nsesh", n_sesh, 
                    "_", packagename, ".csv"]),DataFrame(train_dof = result_df.train_dof))
        CSV.write(join(["MODELno",model, "_composite_th_summary_nboot", n_iters, "_nsesh", n_sesh, 
                    "_", packagename, ".csv"]),DataFrame(th_summary = th_summary))
    end

    set_xaxes_same_scale(axs)
    set_yaxes_same_scale(axs)
    for i=1:length(fs)
        # println(i)
        printFigure(join(["composite_", modelNames[i], "_theta_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=fs[i],figurePath=compositesavepath)
    end
    cd(ret_dir)    
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
    printFigure(join(["compositeAICBIC_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=f, figurePath=compositesavepath)
    
    f = figure(figsize=(20,3))
    ax1=subplot(1,3,1)
    compare_AICBIC(mean_all_Sn_accuracy, all_Sn_accuracy; yl=join(["Train Accuracy nsesh=",n_sesh]), iters=n_iters, ax=ax1, minmax="max")
    ax1.set_ylim([0., 1.])
    ax2=subplot(1,3,2)
    compare_AICBIC(mean_all_test_accuracy, all_test_accuracy; yl=join(["Test Accuracy nsesh=",n_sesh]), iters=n_iters, ax=ax2, minmax="max")
    ax2.set_ylim([0., 1.])
    ax3=subplot(1,3,3)
    printFigure(join(["composite_Accuracy_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=f, figurePath=compositesavepath)
    #
    # Save the variables to the composite folder
    #
    ret_dir = pwd()
    cd(compositesavepath)

    CSV.write(join(["composite_AICs_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(allAICs = allAICs))
    CSV.write(join(["composite_meanAIC_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(mean_all_AICs = mean_all_AICs))
    CSV.write(join(["composite_AICcs_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(allAICcs = allAICcs))
    CSV.write(join(["composite_meanAICc_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(mean_all_AICcs = mean_all_AICcs))
    CSV.write(join(["composite_BICs_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(allBICs = allBICs))
    CSV.write(join(["composite_meanBIC_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(mean_all_BICs = mean_all_BICs))
    
    CSV.write(join(["composite_Sn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(all_Sn_accuracy = all_Sn_accuracy))
    CSV.write(join(["composite_meanSn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(mean_all_Sn_accuracy = mean_all_Sn_accuracy))
    CSV.write(join(["composite_test_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(all_test_accuracy = all_test_accuracy))
    CSV.write(join(["composite_meantest_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
                "_", packagename, ".csv"]),DataFrame(mean_all_test_accuracy = mean_all_test_accuracy))
    
    cd(ret_dir)
end

