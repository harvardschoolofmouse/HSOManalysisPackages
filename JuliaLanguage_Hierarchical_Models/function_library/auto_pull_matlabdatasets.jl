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

function extract_data_with_baselineandLOI(path; normalize=true, useHx=false, history_spacing_s=0.1, omit_cue=false)
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
	if useHx
		df = makeSessionDataFrame(data_single_trial; normalize=normalize, includeBL_LOI=true, baseline_data=data_baseline_trial, 
    		LOI_data=data_LOI_trial, include_history=true, history_spacing_s=history_spacing_s, n_hx_terms = 10, cut_out_cue=omit_cue)
	else
		df = makeSessionDataFrame(data_single_trial; normalize=normalize, includeBL_LOI=true, baseline_data=data_baseline_trial, 
			LOI_data=data_LOI_trial, cut_out_cue=omit_cue)
	end
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
# 	pwd()
	warning("Using new combine_th_across_sessions 3/2/2021")
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
        model_indicies = findall(x->x==modelNames[model], results.results[1].th_summary[1].modelName)
        d = length(model_indicies)
        ths_this_model = [Vector{Float64}(undef, 0) for i=1:length(model_indicies)]
        se_ths_this_model = [Vector{Float64}(undef, 0) for i=1:length(model_indicies)]
        mdofs_this_model = [Vector{Float64}(undef, 0) for i=1:length(model_indicies)]
        

        # only propagate error ACROSS sessions, not within a session... 3/2/2021
        # start by pulling out the composite ths for the model for sesh 1
        for i_sesh = 1:n_sesh
            th_summary_this_sesh = results.results[i_sesh].th_summary[1][model_indicies, :]
            for i_th = 1:d
                push!(ths_this_model[i_th], th_summary_this_sesh.composite_th[i_th])
                push!(se_ths_this_model[i_th], th_summary_this_sesh.composite_se[i_th])
                push!(mdofs_this_model[i_th], th_summary_this_sesh.composite_mdof[i_th])
            end
        end
        # now, get the composite ths for each theta
        composite_th = Vector{Float64}(undef, 0)
        composite_se = Vector{Float64}(undef, 0)
        composite_CImin = Vector{Float64}(undef, 0)
        composite_CImax = Vector{Float64}(undef, 0)
        composite_mdof = Vector{Float64}(undef, 0)
        for i_th = 1:length(model_indicies)
            result_df_this_theta = DataFrame(train_dof = [mdofs_this_model[i_th]], 
                        train_ths = [[ths_this_model[i_th]]], 
                        train_se_ths = [[se_ths_this_model[i_th]]])
            
            
            (meanTh, propagated_se_th, CImin, CImax,mdf) = getCompositeTheta(ths_this_model[i_th], se_ths_this_model[i_th], mdofs_this_model[i_th], propagate=true)
        
            push!(composite_th, meanTh)
            push!(composite_se, propagated_se_th)
            push!(composite_CImin, CImin)
            push!(composite_CImax, CImax)
            push!(composite_mdof, mdf)
        end
        
        f = figure(figsize=(5,3))
	    ax = subplot(1,1,1)
	    plot_with_CI(composite_th, composite_CImin, composite_CImax, ax=ax)
	    ax.set_title(modelNames[model])
	    ax.set_xticks(collect(1:d))
        
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
        
        CSV.write(join([modelNames[model], "_composite_ths_nboot", n_iters, "_nsesh", n_sesh, 
                    ".csv"]),DataFrame(composite_th = composite_th))
        CSV.write(join([modelNames[model], "_composite_se_ths_nboot", n_iters, "_nsesh", n_sesh, 
                    ".csv"]),DataFrame(composite_se = composite_se))
        CSV.write(join([modelNames[model], "_composite_dofs_nboot", n_iters, "_nsesh", n_sesh, 
                    ".csv"]),DataFrame(composite_mdof = composite_mdof))
        CSV.write(join([modelNames[model], "_composite_th_summary_nboot", n_iters, "_nsesh", n_sesh, 
                    ".csv"]),DataFrame(th_summary = th_summary))
    end
    set_xaxes_same_scale(axs)
    set_yaxes_same_scale(axs)

    cd(compositesavepath)
    for i=1:length(fs)
        println(modelNames[i])
        printFigure(join(["composite_", modelNames[i], "_theta_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=fs[i],figurePath=compositesavepath)
    end
    println("Figs saved to:", compositesavepath)
    cd(ret_dir)  
    return th_summary
end;

# function combine_th_across_sessions(results, compositesavepath, runID,packagename)
# 	pwd()
#     ret_dir = pwd()
#     cd(compositesavepath)
#     modelNames = unique(results.results[1].th_summary[1].modelName)
#     n_sesh = length(results.results)
#     i_sesh=1
#     i_model=1
#     i_iter=1
#     n_models = length(results.results[i_sesh].AICs)
#     n_iters = length(results.results[i_sesh].AICs[i_model])
    
#     axs = []
#     fs = []
#     th_summary = DataFrame(modelName=[], composite_th=[], composite_se=[], composite_CImin=[], composite_CImax=[], composite_mdof = [])
#     for model = 1:n_models

#         # result_df = DataFrame(train_dof = [results.results[1].dofs[model][1]], train_ths = [[results.results[1].ths[model][1]]], train_se_ths=[[results.results[1].se_ths[model][1]]])
# 		# only propagate error ACROSS sessions, not within a session... 3/2/2021
# 		result_df_thissesh = 
# 				DataFrame(train_dof = [results.results[1].dofs[model][1]], 
# 					train_ths = [[results.results[1].ths[model][1]]], 
# 					train_se_ths=[[results.results[1].se_ths[model][1]]])
# 		(composite_th, composite_se, composite_CImin, composite_CImax, ax, f, composite_mdof) = 
# 					theta_summary(result_df_thissesh; Mode = "sparseFit", result_df=result_df_thissesh, propagate=false) 
# 		result_df_acrosssesh = 
# 				DataFrame(train_dof = [composite_mdof], 
# 					train_ths = [[composite_th]], 
# 					train_se_ths=[[composite_se]])


#         for i_sesh = 1:n_sesh
#         	n_iters = length(results.results[i_sesh].AICs[model])
#       #   	println("	model:", model, " session: ", i_sesh)
# 	    	# println("n_iters=", n_iters)
# 		    # println("size BIC=", length(results.results[i_sesh].BICs[model]))
# 		    # println("size dofs=", length(results.results[i_sesh].dofs[model]))
# 		    # println("size ths=", length(results.results[i_sesh].ths[model]))
# 		    # println("size se_ths=", length(results.results[i_sesh].se_ths[model]))
#             if i_sesh==1
#                 iterstart=2
#             else
#                 iterstart=1
#             end

#             for i=iterstart:n_iters
#             	# append!(result_df, 
#              #    	DataFrame(train_dof = [results.results[i_sesh].dofs[model][i]], 
#              #            train_ths = [[results.results[i_sesh].ths[model][i]]], 
#              #            train_se_ths=[[results.results[i_sesh].se_ths[model][i]]],
#              #        ))  

#             	# only propagate error ACROSS sessions, not within a session... 3/2/2021
#             	result_df_thissesh = 
# 					DataFrame(train_dof = [results.results[i_sesh].dofs[model][i]], 
#                         train_ths = [[results.results[i_sesh].ths[model][i]]], 
#                         train_se_ths=[[results.results[i_sesh].se_ths[model][i]]],
#                     )
# 				(composite_th, composite_se, composite_CImin, composite_CImax, ax, f, composite_mdof) = 
# 						theta_summary(result_df_thissesh; Mode = "sparseFit", result_df=result_df_thissesh, propagate=false) 
#             	append!(result_df_acrosssesh, 
#                 	DataFrame(train_dof = [composite_mdof], 
#                         train_ths = [[composite_th]], 
#                         train_se_ths=[[composite_se]],
#                     ))  
                 
#             end
#         end
#         result_df = result_df_acrosssesh
        
#         (composite_th, composite_se, composite_CImin, composite_CImax, ax, f, composite_mdof) = theta_summary(result_df; Mode = "sparseFit", result_df=result_df, propagate=true)
#         title(modelNames[model])
#         push!(axs, ax)
#         push!(fs, f)

#         append!(th_summary, 
#             DataFrame(
#                 modelName = modelNames[model],
#                 composite_th=composite_th, 
#                 composite_se=composite_se, 
#                 composite_CImin=composite_CImin, 
#                 composite_CImax=composite_CImax,
#                 composite_mdof = composite_mdof,
#                 ))
        
#         #
#         # Save the composite variables
#         #
#         # CSV.write(join(["MODELno",model, "_composite_ths_nboot", n_iters, "_nsesh", n_sesh, 
#         #             "_", packagename, ".csv"]),DataFrame(train_ths = result_df.train_ths))
#         # CSV.write(join(["MODELno",model, "_composite_se_ths_nboot", n_iters, "_nsesh", n_sesh, 
#         #             "_", packagename, ".csv"]),DataFrame(train_se_ths = result_df.train_se_ths))
#         # CSV.write(join(["MODELno",model, "_composite_dofs_nboot", n_iters, "_nsesh", n_sesh, 
#         #             "_", packagename, ".csv"]),DataFrame(train_dof = result_df.train_dof))
#         # CSV.write(join(["MODELno",model, "_composite_th_summary_nboot", n_iters, "_nsesh", n_sesh, 
#         #             "_", packagename, ".csv"]),DataFrame(th_summary = th_summary))
#         # println(pwd())
#         # CSV.write("test.csv",DataFrame(train_ths = result_df.train_ths))

#         CSV.write(join([modelNames[model], "_composite_ths_nboot", n_iters, "_nsesh", n_sesh, 
#                     ".csv"]),DataFrame(train_ths = result_df.train_ths))
#         CSV.write(join([modelNames[model], "_composite_se_ths_nboot", n_iters, "_nsesh", n_sesh, 
#                     ".csv"]),DataFrame(train_se_ths = result_df.train_se_ths))
#         CSV.write(join([modelNames[model], "_composite_dofs_nboot", n_iters, "_nsesh", n_sesh, 
#                     ".csv"]),DataFrame(train_dof = result_df.train_dof))
#         CSV.write(join([modelNames[model], "_composite_th_summary_nboot", n_iters, "_nsesh", n_sesh, 
#                     ".csv"]),DataFrame(th_summary = th_summary))
#     end

#     set_xaxes_same_scale(axs)
#     set_yaxes_same_scale(axs)
#     println("Figs saved to:", compositesavepath)
#     cd(compositesavepath)
#     for i=1:length(fs)
#         # println(i)
#         # printFigure(join(["composite_", modelNames[i], "_theta_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=fs[i],figurePath=compositesavepath)

#         printFigure(join(["composite_", modelNames[i], "_theta_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=fs[i],figurePath=compositesavepath)
#     end
#     cd(ret_dir)  
#     return th_summary  
# end


# function combine_AICBIC_across_sessions(results, compositesavepath, runID, packagename)
#     n_sesh = length(results.results)
#     i_sesh=1
#     i_model=1
#     i_iter=1
#     n_models = length(results.results[i_sesh].AICs)
#     n_iters = length(results.results[i_sesh].AICs[i_model])
#     allAICs = [[] for _=1:n_models]
#     allAICcs = [[] for _=1:n_models]
#     allBICs = [[] for _=1:n_models]
#     all_Sn_accuracy = [[] for _=1:n_models]
#     all_test_accuracy = [[] for _=1:n_models]
    
#     for i_sesh = 1:n_sesh
#         for i_model = 1:n_models
#         	n_iters = length(results.results[i_sesh].AICs[i_model])
#       #   	println("	model:", i_model, " session: ", i_sesh)
# 	    	# println("n_iters=", n_iters)
# 		    # println("size BIC=", length(results.results[i_sesh].BICs[i_model]))
# 		    # println("size dofs=", length(results.results[i_sesh].dofs[i_model]))
# 		    # println("size ths=", length(results.results[i_sesh].ths[i_model]))
# 		    # println("size se_ths=", length(results.results[i_sesh].se_ths[i_model]))
#             AICs = [results.results[i_sesh].AICs[i_model][i_iter] for i_iter = 1:n_iters]
#             append!(allAICs[i_model], AICs)
#             AICcs = [results.results[i_sesh].AICcs[i_model][i_iter] for i_iter = 1:n_iters]
#             append!(allAICcs[i_model], AICcs)
#             BICs = [results.results[i_sesh].BICs[i_model][i_iter] for i_iter = 1:n_iters]
#             append!(allBICs[i_model], BICs)
#             Sn_accuracies = [results.results[i_sesh].Sn_accuracy[i_model][i_iter] for i_iter = 1:n_iters]
#             append!(all_Sn_accuracy[i_model], Sn_accuracies)
#             test_accuracies = [results.results[i_sesh].test_accuracy[i_model][i_iter] for i_iter = 1:n_iters]
#             append!(all_test_accuracy[i_model], test_accuracies)
#         end
#     end
#     mean_all_AICs = [mean(allAICs[i]) for i=1:n_models]
#     mean_all_AICcs = [mean(allAICcs[i]) for i=1:n_models]
#     mean_all_BICs = [mean(allBICs[i]) for i=1:n_models]
#     mean_all_Sn_accuracy = [mean(all_Sn_accuracy[i]) for i=1:n_models]
#     mean_all_test_accuracy = [mean(all_test_accuracy[i]) for i=1:n_models]
#     f = figure(figsize=(20,3))
#     ax1=subplot(1,3,1)
#     compare_AICBIC(mean_all_AICs, allAICs; yl=join(["AIC nsesh=",n_sesh]), iters=n_iters, ax=ax1, minmax="min")
#     ax2=subplot(1,3,2)
#     compare_AICBIC(mean_all_AICcs, allAICcs; yl=join(["AICc nsesh=",n_sesh]), iters=n_iters, ax=ax2, minmax="min")
#     ax3=subplot(1,3,3)
#     compare_AICBIC(mean_all_BICs, allBICs; yl=join(["BIC nsesh=",n_sesh]), iters=n_iters, ax=ax3, minmax="min")
#     # printFigure(join(["compositeAICBIC_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=f, figurePath=compositesavepath)
#     printFigure(join(["compositeAICBIC_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=f, figurePath=compositesavepath)
    
#     f = figure(figsize=(20,3))
#     ax1=subplot(1,3,1)
#     compare_AICBIC(mean_all_Sn_accuracy, all_Sn_accuracy; yl=join(["Train Accuracy nsesh=",n_sesh]), iters=n_iters, ax=ax1, minmax="max")
#     ax1.set_ylim([0., 1.])
#     ax2=subplot(1,3,2)
#     compare_AICBIC(mean_all_test_accuracy, all_test_accuracy; yl=join(["Test Accuracy nsesh=",n_sesh]), iters=n_iters, ax=ax2, minmax="max")
#     ax2.set_ylim([0., 1.])
#     ax3=subplot(1,3,3)
#     # printFigure(join(["composite_Accuracy_summary_nboot", n_iters, "_nsesh", n_sesh, "_", packagename]); fig=f, figurePath=compositesavepath)
#     printFigure(join(["composite_Accuracy_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=f, figurePath=compositesavepath)
#     #
#     # Save the variables to the composite folder
#     #
#     ret_dir = pwd()
#     cd(compositesavepath)

#     # CSV.write(join(["composite_AICs_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(allAICs = allAICs))
#     # CSV.write(join(["composite_meanAIC_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(mean_all_AICs = mean_all_AICs))
#     # CSV.write(join(["composite_AICcs_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(allAICcs = allAICcs))
#     # CSV.write(join(["composite_meanAICc_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(mean_all_AICcs = mean_all_AICcs))
#     # CSV.write(join(["composite_BICs_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(allBICs = allBICs))
#     # CSV.write(join(["composite_meanBIC_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(mean_all_BICs = mean_all_BICs))
    
#     # CSV.write(join(["composite_Sn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(all_Sn_accuracy = all_Sn_accuracy))
#     # CSV.write(join(["composite_meanSn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(mean_all_Sn_accuracy = mean_all_Sn_accuracy))
#     # CSV.write(join(["composite_test_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(all_test_accuracy = all_test_accuracy))
#     # CSV.write(join(["composite_meantest_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
#     #             "_", packagename, ".csv"]),DataFrame(mean_all_test_accuracy = mean_all_test_accuracy))
#     CSV.write(join(["composite_AICs_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(allAICs = allAICs))
#     CSV.write(join(["composite_meanAIC_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(mean_all_AICs = mean_all_AICs))
#     CSV.write(join(["composite_AICcs_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(allAICcs = allAICcs))
#     CSV.write(join(["composite_meanAICc_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(mean_all_AICcs = mean_all_AICcs))
#     CSV.write(join(["composite_BICs_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(allBICs = allBICs))
#     CSV.write(join(["composite_meanBIC_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(mean_all_BICs = mean_all_BICs))
    
#     CSV.write(join(["composite_Sn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(all_Sn_accuracy = all_Sn_accuracy))
#     CSV.write(join(["composite_meanSn_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(mean_all_Sn_accuracy = mean_all_Sn_accuracy))
#     CSV.write(join(["composite_test_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(all_test_accuracy = all_test_accuracy))
#     CSV.write(join(["composite_meantest_accuracy_nboot", n_iters, "_nsesh", n_sesh, 
#                 ".csv"]),DataFrame(mean_all_test_accuracy = mean_all_test_accuracy))
    
#     cd(ret_dir)
# end
function combine_AICBIC_across_sessions(results, compositesavepath, runID, packagename)
    # this new version will propagate error across sessions rather than using the range on all sesh...
    modelNames = unique(results.results[1].th_summary[1].modelName)
    warning("using new combine_AICBIC_across_sessions2 that propagates error across sessions")
    n_sesh = length(results.results)
    i_sesh=1
    i_model=1
    i_iter=1
    n_models = length(results.results[i_sesh].AICs)
    n_iters = length(results.results[i_sesh].AICs[i_model])
    
    compositeAICs = nanmat(n_models, 1)
    composite_stdAICs = nanmat(n_models, 1)
    compositeAICcs = nanmat(n_models, 1)
    composite_stdAICcs = nanmat(n_models, 1)
    compositeBICs = nanmat(n_models, 1)
    composite_stdBICs = nanmat(n_models, 1)
    composite_Sn_accuracys = nanmat(n_models, 1)
    composite_std_Sn_accuracys = nanmat(n_models, 1)
    composite_accuracys_test = nanmat(n_models, 1)
    composite_std_accuracys_test = nanmat(n_models, 1)
    composite_Sn_dev_explained = nanmat(n_models, 1)
    composite_std_Sn_dev_explained = nanmat(n_models, 1)
    composite_dev_explained = nanmat(n_models, 1)
    composite_std_dev_explained = nanmat(n_models, 1)
    
    CImin_stdAICs = nanmat(n_models, 1)
    CImin_stdAICcs = nanmat(n_models, 1)
    CImin_stdBICs = nanmat(n_models, 1)
    CImin_std_Sn_accuracys = nanmat(n_models, 1)
    CImin_std_accuracys_test = nanmat(n_models, 1)
    CImin_std_Sn_dev_explained = nanmat(n_models, 1)
    CImin_std_dev_explained = nanmat(n_models, 1)
    
    CImax_stdAICs = nanmat(n_models, 1)
    CImax_stdAICcs = nanmat(n_models, 1)
    CImax_stdBICs = nanmat(n_models, 1)
    CImax_std_Sn_accuracys = nanmat(n_models, 1)
    CImax_std_accuracys_test = nanmat(n_models, 1)
    CImax_std_Sn_dev_explained = nanmat(n_models, 1)
    CImax_std_dev_explained = nanmat(n_models, 1)

    
    for i_model = 1:n_models
        model_indicies = findall(x->x==modelNames[i_model], results.results[1].th_summary[1].modelName)
        meanAICs = Vector{Float64}(undef, 0)
        stdAICs = Vector{Float64}(undef, 0)
        meanAICcs = Vector{Float64}(undef, 0)
        stdAICcs = Vector{Float64}(undef, 0)
        meanBICs = Vector{Float64}(undef, 0)
        stdBICs = Vector{Float64}(undef, 0)
        mean_Sn_accuracys = Vector{Float64}(undef, 0)
        std_Sn_accuracys = Vector{Float64}(undef, 0)
        mean_accuracys_test = Vector{Float64}(undef, 0)
        std_accuracys_test = Vector{Float64}(undef, 0)
        mean_Sn_dev_explained = Vector{Float64}(undef, 0)
        std_Sn_dev_explained = Vector{Float64}(undef, 0)
        mean_dev_explained = Vector{Float64}(undef, 0)
        std_dev_explained = Vector{Float64}(undef, 0)

        dof = []
        for i_sesh = 1:n_sesh
            th_summary_this_sesh = results.results[i_sesh].th_summary[1][model_indicies, :]
            push!(dof, th_summary_this_sesh.composite_mdof[1])
            mdofs_this_model = [Vector{Float64}(undef, 0) for i=1:length(model_indicies)]
            
            push!(meanAICs, results.results[i_sesh].meanAIC[i_model])
            push!(meanAICcs, results.results[i_sesh].meanAICc[i_model])
            push!(meanBICs, results.results[i_sesh].meanBIC[i_model])
            push!(mean_Sn_accuracys, results.results[i_sesh].meanAccuracy_Sn[i_model])
            push!(mean_accuracys_test, results.results[i_sesh].meanAccuracy_test[i_model])
            push!(mean_Sn_dev_explained, results.results[i_sesh].mean_deviance_explained_Sn[i_model])
            push!(mean_dev_explained, results.results[i_sesh].mean_deviance_explained[i_model])
            
            push!(stdAICs, results.results[i_sesh].stdAIC[i_model])
            push!(stdAICcs, results.results[i_sesh].stdAICc[i_model])
            push!(stdBICs, results.results[i_sesh].stdBIC[i_model])
            push!(std_Sn_accuracys, results.results[i_sesh].std_Sn_accuracy[i_model])
            push!(std_accuracys_test, results.results[i_sesh].std_accuracy_test[i_model])
            push!(std_Sn_dev_explained, results.results[i_sesh].std_Sn_dev_explained[i_model])
            push!(std_dev_explained, results.results[i_sesh].std_dev_explained[i_model])            

        end
        # We can now propagate the std across sessions using the error propagation method, same as th...
        (compositeAICs[i_model], composite_stdAICs[i_model], CImin_stdAICs[i_model], CImax_stdAICs[i_model],mdf) = getCompositeTheta(meanAICs, stdAICs, dof, propagate=true)
        (compositeAICcs[i_model], composite_stdAICcs[i_model], CImin_stdAICcs[i_model], CImax_stdAICcs[i_model],mdf) = getCompositeTheta(meanAICcs, stdAICcs, dof, propagate=true)
        (compositeBICs[i_model], composite_stdBICs[i_model], CImin_stdBICs[i_model], CImax_stdBICs[i_model],mdf) = getCompositeTheta(meanBICs, stdBICs, dof, propagate=true)
        (composite_Sn_accuracys[i_model], composite_std_Sn_accuracys[i_model], CImin_std_Sn_accuracys[i_model], CImax_std_Sn_accuracys[i_model],mdf) = getCompositeTheta(mean_Sn_accuracys, std_Sn_accuracys, dof, propagate=true)
        (composite_accuracys_test[i_model], composite_std_accuracys_test[i_model], CImin_std_accuracys_test[i_model], CImax_std_accuracys_test[i_model],mdf) = getCompositeTheta(mean_accuracys_test, std_accuracys_test, dof, propagate=true)
        (composite_Sn_dev_explained[i_model], composite_std_Sn_dev_explained[i_model], CImin_std_Sn_dev_explained[i_model], CImax_std_Sn_dev_explained[i_model],mdf) = getCompositeTheta(mean_Sn_dev_explained, std_Sn_dev_explained, dof, propagate=true)
        (composite_dev_explained[i_model], composite_std_dev_explained[i_model], CImin_std_dev_explained[i_model], CImax_std_dev_explained[i_model],mdf) = getCompositeTheta(mean_dev_explained, std_dev_explained, dof, propagate=true)
        
    end
    
    
    compositeAICs = vec(compositeAICs)
    compositeAICcs = vec(compositeAICcs)
    compositeBICs = vec(compositeBICs)
    composite_Sn_accuracys = vec(composite_Sn_accuracys)
    composite_accuracys_test = vec(composite_accuracys_test)
    composite_Sn_dev_explained = vec(composite_Sn_dev_explained)
    composite_dev_explained = vec(composite_dev_explained)
    
    CImin_stdAICs = vec(CImin_stdAICs)
    CImin_stdAICcs = vec(CImin_stdAICcs)
    CImin_stdBICs = vec(CImin_stdBICs)
    CImin_std_Sn_accuracys = vec(CImin_std_Sn_accuracys)
    CImin_std_accuracys_test = vec(CImin_std_accuracys_test)
    CImin_std_Sn_dev_explained = vec(CImin_std_Sn_dev_explained)
    CImin_std_dev_explained = vec(CImin_std_dev_explained)
    
    CImax_stdAICs =vec(CImax_stdAICs)
    CImax_stdAICcs =vec(CImax_stdAICcs)
    CImax_stdBICs =vec(CImax_stdBICs)
    CImax_std_Sn_accuracys =vec(CImax_std_Sn_accuracys)
    CImax_std_accuracys_test =vec(CImax_std_accuracys_test)
    CImax_std_Sn_dev_explained =vec(CImax_std_Sn_dev_explained)
    CImax_std_dev_explained =vec(CImax_std_dev_explained)
    
    

    f = figure(figsize=(20,3))
    ax1=subplot(1,3,1)
    plot_with_CI_min(compositeAICs, CImin_stdAICs, CImax_stdAICs; 
        ax=ax1, ylab=join(["AIC"]), tit=join(["meanAIC nsesh=",n_sesh, " iters=", n_iters]), xl="Model #")
    ax2=subplot(1,3,2)
    plot_with_CI_min(compositeAICcs, CImin_stdAICcs, CImax_stdAICcs; 
        ax=ax2, ylab=join(["AICc"]), tit=join(["meanAICc nsesh=",n_sesh, " iters=", n_iters]), xl="Model #")
    ax3=subplot(1,3,3)
    plot_with_CI_min(compositeBICs, CImin_stdBICs, CImax_stdBICs; 
        ax=ax3, ylab=join(["BIC"]), tit=join(["meanBIC nsesh=",n_sesh, " iters=", n_iters]), xl="Model #")
    printFigure(join(["compositeAICBIC_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=f, figurePath=compositesavepath)
    
    f = figure(figsize=(20,3))
    ax1=subplot(1,3,1)
    plot_with_CI_max(composite_Sn_accuracys, CImin_std_Sn_accuracys, CImax_std_Sn_accuracys; 
        ax=ax1, ylab=join(["% correct"]), tit=join(["Train Accuracy nsesh=",n_sesh, " iters=", n_iters]), xl="Model #")
    ax1.set_ylim([0., 1.])
    ax2=subplot(1,3,2)
    plot_with_CI_max(composite_accuracys_test, CImin_std_accuracys_test, CImax_std_accuracys_test; 
        ax=ax2, ylab=join(["% correct"]), tit=join(["Test Accuracy nsesh=",n_sesh, " iters=", n_iters]), xl="Model #")
    ax2.set_ylim([0., 1.])
    ax3=subplot(1,3,3)
    printFigure(join(["composite_Accuracy_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=f, figurePath=compositesavepath)
    
    f = figure(figsize=(20,3))
    ax1=subplot(1,3,1)
    plot_with_CI_max(composite_Sn_dev_explained, CImin_std_Sn_dev_explained, CImax_std_Sn_dev_explained; 
        ax=ax1, ylab=join(["R^2"]), tit=join(["Train Dev Explained nsesh=",n_sesh, " iters=", n_iters]), xl="Model #")

    warning("is new")
    
    if sum(isnan.(composite_Sn_dev_explained)) == length(composite_Sn_dev_explained) 
    	ax1.set_ylim([-0.1, 0.2])
    else
	    ax1.set_ylim([-0.1, nanmax(composite_Sn_dev_explained) + 0.1])
    end
    ax2=subplot(1,3,2)
    plot_with_CI_max(composite_dev_explained, CImin_std_dev_explained, CImax_std_dev_explained; 
        ax=ax2, ylab=join(["R^2"]), tit=join(["Test Dev Explained nsesh=",n_sesh, " iters=", n_iters]), xl="Model #")
    # ax2.set_ylim([-0.1, 0.2])
    if sum(isnan.(composite_dev_explained)) == length(composite_dev_explained) 
    	ax2.set_ylim([-0.1, 0.2])
    else
	    ax2.set_ylim([-0.1, nanmax(composite_dev_explained) + 0.1])
    end
    ax3=subplot(1,3,3)
    printFigure(join(["composite_DevExpl_summary_nboot", n_iters, "_nsesh", n_sesh]); fig=f, figurePath=compositesavepath)
    #
    # Save the variables to the composite folder
    #
    ret_dir = pwd()
    cd(compositesavepath)
    CSV.write(join(["compositeAICs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(compositeAICs = compositeAICs))
    CSV.write(join(["compositeAICcs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(compositeAICcs = compositeAICcs))
    CSV.write(join(["compositeBICs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(compositeBICs = compositeBICs))
    CSV.write(join(["composite_Sn_accuracys_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(composite_Sn_accuracys = composite_Sn_accuracys))
    CSV.write(join(["composite_accuracys_test_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(composite_accuracys_test = composite_accuracys_test))
    CSV.write(join(["composite_Sn_dev_explained_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(composite_Sn_dev_explained = composite_Sn_dev_explained))
    CSV.write(join(["composite_dev_explained_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(composite_dev_explained = composite_dev_explained))
    
    
    CSV.write(join(["CI_AICs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(CImin_stdAICs = CImin_stdAICs, CImax_stdAICs = CImax_stdAICs))
    CSV.write(join(["CI_AICcs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(CImin_stdAICcs = CImin_stdAICcs, CImax_stdAICcs = CImax_stdAICcs))
    CSV.write(join(["CI_BICs_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(CImin_stdBICs = CImin_stdBICs, CImax_stdBICs = CImax_stdBICs))
    CSV.write(join(["CI_Sn_accuracys_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(CImin_std_Sn_accuracys = CImin_std_Sn_accuracys, CImax_std_Sn_accuracys = CImax_std_Sn_accuracys))
    CSV.write(join(["CI_accuracys_test_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(CImin_std_accuracys_test = CImin_std_accuracys_test, CImax_std_accuracys_test = CImax_std_accuracys_test))
    CSV.write(join(["CI_Sn_dev_explained_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(CImin_std_Sn_dev_explained = CImin_std_Sn_dev_explained, CImax_std_Sn_dev_explained = CImax_std_Sn_dev_explained))
    CSV.write(join(["CI_dev_explained_nboot", n_iters, "_nsesh", n_sesh, 
                ".csv"]),DataFrame(CImin_std_dev_explained = CImin_std_dev_explained, CImax_std_dev_explained = CImax_std_dev_explained))
    

    cd(ret_dir)
    # let's return the best model ID for each assessment:
    return(compositeAICs,
	compositeAICcs,
	compositeBICs,
	composite_Sn_accuracys,
	composite_accuracys_test,
	composite_Sn_dev_explained,
	composite_dev_explained)
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
	# formulas = [
	#     @formula(LickState ~ Y),
	#     @formula(LickState ~ Mean_Baseline),
	#     @formula(LickState ~ Median_Baseline),
	#     @formula(LickState ~ Mean_LOI),
	#     @formula(LickState ~ Median_LOI),
	#     @formula(LickState ~ Mean_Baseline + Mean_LOI + Y),
	#     @formula(LickState ~ Median_Baseline + Median_LOI + Y),
	# 	]

	# modelNames = [
	#     "DA-only",
	#     "μBl-only",
	#     "medBl-only",
	#     "μLOI-only",
	#     "medLOI-only",
	#     "DA-μBl-μLOI",
	#     "DA-medBl-medLOI",
	# ]

	formulas = [
		    @formula(LickState ~ Y),
		    @formula(LickState ~ Mean_Baseline),
		    @formula(LickState ~ Median_Baseline),
		    @formula(LickState ~ Mean_LOI),
		    @formula(LickState ~ Median_LOI),
		    @formula(LickState ~ Mean_Baseline + Mean_LOI + Y),
		    @formula(LickState ~ Median_Baseline + Median_LOI + Y),
		    @formula(LickState ~ LickTime_2back),
		    @formula(LickState ~ LickTime_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back),
		    @formula(LickState ~ Rxn_1back + Early_1back + Reward_1back + ITI_1back),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back + Y),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back +Mean_Baseline + Mean_LOI + Y),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back + Median_Baseline + Median_LOI + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Rxn_2back + Early_2back + Reward_2back + ITI_2back +Rxn_1back + Early_1back + Reward_1back + ITI_1back + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Rxn_2back + Early_2back + Reward_2back + ITI_2back +Rxn_1back + Early_1back + Reward_1back + ITI_1back + Mean_Baseline + Mean_LOI + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back + Median_Baseline + Median_LOI + Y),
			]

		modelNames = [
		    "DA",
		    "μBl",
		    "medBl",
		    "μLOI",
		    "medLOI",
		    "DA-μBl-μLOI",
		    "DA-medBl-medLOI",
		    "Lt2b",
		    "Lt1b",
		    "Lt1b-Lt2b",
		    "oc2b",
		    "oc1b",
		    "oc1b-oc2b",
		    "DA-oc1b-oc2b",
		    "DA-μBlLOI-oc12b",
		    "DA-mdBlLOI-oc12b",
		    "DA-Ltoc12b",
		    "DA-μBlLO-Ltoc12b",
		    "DA-mdBlLO-Ltoc12b",
		]

	if length(modelNames) != length(formulas)
		error("Model names and formulas not matched...")
	end

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


function bootlogit_timeslice_modelpackage2(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	warning("timeslice model updated 3-3-2021 for new Hx predictors and updates to error prop...")
# name the package and runID
	packagename = join(["bl_ts_pkg2_",runID])
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s = 0.2)
	ndf[:NoReward_1back] = (ndf[:Rxn_1back] .+ ndf[:Early_1back] .+ ndf[:Reward_1back] .+ ndf[:ITI_1back]) .== 1 
	ndf[:NoReward_2back] = (ndf[:Rxn_2back] .+ ndf[:Early_2back] .+ ndf[:Reward_2back] .+ ndf[:ITI_2back]) .== 1 
	#
	# next, we need to specify binning for our timeslices
	#. as a start to match the hazard analysis, I'll use 250 ms bins
	#
	slice_width_ms = 250.#3000.#250.
	println("slice_width_ms: ", slice_width_ms)
	(binned_ndfs, binEdges) = slice_dataframe_into_timebins(ndf, slice_width_ms)
	# name the timeslices
	timeslice_names = []
	for i = 1:round(Int, 7000 ./slice_width_ms)#length(binned_ndfs)
		push!(timeslice_names, join(["_", round(Int,1000*binEdges[i]), "-", round(Int,1000*binEdges[i+1])]))
	end
	


	formulas = [
			@formula(LickState ~ LickTime_2back),
		    @formula(LickState ~ LickTime_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Y),
		    # @formula(LickState ~ LickTime_2back + LickTime_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ NoReward_2back),
		    @formula(LickState ~ NoReward_1back),
		    @formula(LickState ~ NoReward_2back + NoReward_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + NoReward_2back + NoReward_1back),
		    @formula(LickState ~ NoReward_2back + NoReward_1back + Y),
		    # @formula(LickState ~ NoReward_2back + NoReward_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + NoReward_2back + NoReward_1back + Y),
		    # @formula(LickState ~ LickTime_2back + LickTime_1back + NoReward_2back + NoReward_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ Y),
		   #  @formula(LickState ~ Hx10),
		   #  @formula(LickState ~ Hx10 + Hx9),
		   #  @formula(LickState ~ Hx10 + Hx9 + Hx8),
		   #  @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7),
		   #  @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6),
		   #  @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5),
		  	# @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4),
		  	# @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3),
		  	# @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2),
		  	# @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
			]

		modelNames = [
			"Lt2b",
		    "Lt1b",
		    "Lt1b-Lt2b",
		    "Lt1b-Lt2b-DA",
		    # "Lt12b-Hx0-2-DA",
		    "oc2b",
		    "oc1b",
		    "oc1b-oc2b",
		    "Lt12b-oc12b",
		    "oc1b-oc2b-DA",
		    # "oc12b-Hx0-2-DA",
		    "Lt12b-oc12b-DA",
		    # "Lt12b-oc12b-Hx0-2-DA",
		    "DA",
		    # "Hx2s",
		    # "Hx2s_1-6s",
		    # "Hx2s_1-4s",
		    # "Hx2s_1-2s",
		    # "Hx2s_1s",
		    # "Hx2s_-8s",
		    # "Hx2s_-6s",
		    # "Hx2s_-4s",
		    # "Hx2s_-2s",
		    # "Hx2s_-2s_DA",
		]

		predictors = [
			[:LickTime_2back],
		    [:LickTime_1back],
		    [:LickTime_2back, :LickTime_1back],
		    [:LickTime_2back, :LickTime_1back, :Y],
		    # [:LickTime_2back, :LickTime_1back, :Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2, :Y],
		    [:NoReward_2back],
		    [:NoReward_1back],
		    [:NoReward_2back, :NoReward_1back],
		    [:LickTime_2back, :LickTime_1back, :NoReward_2back, :NoReward_1back],
		    [:NoReward_2back, :NoReward_1back, :Y],
		    # [:NoReward_2back, :NoReward_1back, :Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2, :Y],
		    [:LickTime_2back, :LickTime_1back, :NoReward_2back, :NoReward_1back, :Y],
		    # [:LickTime_2back, :LickTime_1back, :NoReward_2back, :NoReward_1back, :Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2, :Y],
		    [:Y],
		   #  [:Hx10],
		   #  [:Hx10, :Hx9],
		   #  [:Hx10, :Hx9, :Hx8],
		   #  [:Hx10, :Hx9, :Hx8, :Hx7],
		   #  [:Hx10, :Hx9, :Hx8, :Hx7, :Hx6],
		   #  [:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5],
		  	# [:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4],
		  	# [:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3],
		  	# [:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2],
		  	# [:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2, :Y],
			]

	if length(modelNames) != length(formulas)
		error("Model names and formulas not matched...")
	end

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
		result[:predictors] = predictors
		push!(results, result)
	end

	

# make a working result list of dfs with all the results to keep in workspace
	return results#, ndf
end








function extract_testables_bootlogit_timeslice_modelpackage2(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	# THIS JUST LETS US PULL VARIABLES FOR TESTING!!!!!!!!!!!!!



	warning("timeslice model updated 3-3-2021 for new Hx predictors and updates to error prop...")
# name the package and runID
	packagename = join(["bl_ts_pkg2_",runID])
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s = 0.2)
	ndf[:NoReward_1back] = (ndf[:Rxn_1back] .+ ndf[:Early_1back] .+ ndf[:Reward_1back] .+ ndf[:ITI_1back]) .== 1 
	ndf[:NoReward_2back] = (ndf[:Rxn_2back] .+ ndf[:Early_2back] .+ ndf[:Reward_2back] .+ ndf[:ITI_2back]) .== 1 
	#
	# next, we need to specify binning for our timeslices
	#. as a start to match the hazard analysis, I'll use 250 ms bins
	#
	slice_width_ms = 250.#3000.#250.
	println("slice_width_ms: ", slice_width_ms)
	(binned_ndfs, binEdges) = slice_dataframe_into_timebins(ndf, slice_width_ms)
	# name the timeslices
	timeslice_names = []
	for i = 1:round(Int, 7000 ./slice_width_ms)#length(binned_ndfs)
		push!(timeslice_names, join(["_", round(Int,1000*binEdges[i]), "-", round(Int,1000*binEdges[i+1])]))
	end
	


	formulas = [
			@formula(LickState ~ LickTime_2back),
		    @formula(LickState ~ LickTime_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ NoReward_2back),
		    @formula(LickState ~ NoReward_1back),
		    @formula(LickState ~ NoReward_2back + NoReward_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + NoReward_2back + NoReward_1back),
		    @formula(LickState ~ NoReward_2back + NoReward_1back + Y),
		    @formula(LickState ~ NoReward_2back + NoReward_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + NoReward_2back + NoReward_1back + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + NoReward_2back + NoReward_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ Y),
		    @formula(LickState ~ Hx10),
		    @formula(LickState ~ Hx10 + Hx9),
		    @formula(LickState ~ Hx10 + Hx9 + Hx8),
		    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7),
		    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6),
		    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5),
		  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4),
		  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3),
		  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2),
		  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
			]

		modelNames = [
			"Lt2b",
		    "Lt1b",
		    "Lt1b-Lt2b",
		    "Lt1b-Lt2b-DA",
		    "Lt12b-Hx0-2-DA",
		    "oc2b",
		    "oc1b",
		    "oc1b-oc2b",
		    "Lt12b-oc12b",
		    "oc1b-oc2b-DA",
		    "oc12b-Hx0-2-DA",
		    "Lt12b-oc12b-DA",
		    "Lt12b-oc12b-Hx0-2-DA",
		    "DA",
		    "Hx2s",
		    "Hx2s_1-6s",
		    "Hx2s_1-4s",
		    "Hx2s_1-2s",
		    "Hx2s_1s",
		    "Hx2s_-8s",
		    "Hx2s_-6s",
		    "Hx2s_-4s",
		    "Hx2s_-2s",
		    "Hx2s_-2s_DA",
		]

		predictors = [
			[:LickTime_2back],
		    [:LickTime_1back],
		    [:LickTime_2back, :LickTime_1back],
		    [:LickTime_2back, :LickTime_1back, :Y],
		    [:LickTime_2back, :LickTime_1back, :Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2, :Y],
		    [:NoReward_2back],
		    [:NoReward_1back],
		    [:NoReward_2back, :NoReward_1back],
		    [:LickTime_2back, :LickTime_1back, :NoReward_2back, :NoReward_1back],
		    [:NoReward_2back, :NoReward_1back, :Y],
		    [:NoReward_2back, :NoReward_1back, :Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2, :Y],
		    [:LickTime_2back, :LickTime_1back, :NoReward_2back, :NoReward_1back, :Y],
		    [:LickTime_2back, :LickTime_1back, :NoReward_2back, :NoReward_1back, :Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2, :Y],
		    [:Y],
		    [:Hx10],
		    [:Hx10, :Hx9],
		    [:Hx10, :Hx9, :Hx8],
		    [:Hx10, :Hx9, :Hx8, :Hx7],
		    [:Hx10, :Hx9, :Hx8, :Hx7, :Hx6],
		    [:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5],
		  	[:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4],
		  	[:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3],
		  	[:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2],
		  	[:Hx10, :Hx9, :Hx8, :Hx7, :Hx6, :Hx5, :Hx4, :Hx3, :Hx2, :Y],
			]

	if length(modelNames) != length(formulas)
		error("Model names and formulas not matched...")
	end

	nbins = length(binEdges) - 1
	results = DataFrame(ndf = [ndf],
		slice_width_ms = [slice_width_ms],
		binned_ndfs = [binned_ndfs],
		binEdges = [binEdges],
		timeslice_names = [timeslice_names],
		formulas = [formulas],
		modelNames = [modelNames],
		predictors = [predictors],
		nbins = [nbins],
		figurePath = [figurePath],
		savepath = [savepath])

	
	# extract from with:
	# ndf = variables_r[:ndf][1]
	# slice_width_ms = variables_r[:slice_width_ms][1]
	# binned_ndfs = variables_r[:binned_ndfs][1]
	# binEdges = variables_r[:binEdges][1]
	# timeslice_names = variables_r[:timeslice_names][1]
	# formulas = variables_r[:formulas][1]
	# modelNames = variables_r[:modelNames][1]
	# predictors = variables_r[:predictors][1]
	# nbins = variables_r[:nbins][1]
	# figurePath = variables_r[:figurePath][1]
	# savepath = variables_r[:savepath][1]


# make a working result list of dfs with all the results to keep in workspace
	return results#, ndf
end









function bootlogit_timeslice_postprocessingfunction1(results::DataFrame, compositesavepath, modelpackagefunction; runID=0, inclusionthresh=10.)
	#
	# Use this to compile the analysis when we have a LIST of result dfs
	#
	# name the package and runID
	packagename = modelpackagefunction(""; sessionID ="", getpackagename=true, runID=runID)
	println(packagename)

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

	compositeAICs = []
	compositeAICcs = []
	compositeBICs = []
	composite_Sn_accuracys = []
	composite_accuracy_tests = []
	composite_Sn_dev_explaineds = []
	composite_dev_explaineds = []
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
		(compositeAIC,
	compositeAICc,
	compositeBIC,
	composite_Sn_accuracy,
	composite_accuracy_test,
	composite_Sn_dev_explained,
	composite_dev_explained) = combine_AICBIC_across_sessions(results_df_this_slice, sav_dir, join([runID,slicepath]), packagename)
		cd("..")
		println(names(th_summary))
		th_summary[:sliceID] = [slicepath for _=1:nrow(th_summary)]
		push!(by_slice_composite_ths,th_summary)

		push!(compositeAICs, compositeAIC)
		push!(compositeAICcs, compositeAICc)
		push!(compositeBICs, compositeBIC)
		push!(composite_Sn_accuracys, composite_Sn_accuracy)
		push!(composite_accuracy_tests, composite_accuracy_test)
		push!(composite_Sn_dev_explaineds, composite_Sn_dev_explained)
		push!(composite_dev_explaineds, composite_dev_explained)
	end
	plot_th_vs_timeslice(by_slice_composite_ths,savedir=sdd, inclusionthresh=inclusionthresh)

	plot_composite_AIC_slice(compositeAICs, "AIC", sdd)
	plot_composite_AIC_slice(compositeAICcs, "AICc", sdd)
	plot_composite_AIC_slice(compositeBICs, "BIC", sdd)
	plot_composite_AIC_slice(composite_Sn_accuracys, "Train Accuracy", sdd)
	plot_composite_AIC_slice(composite_accuracy_tests, "Test Accuracy", sdd)
	plot_composite_AIC_slice(composite_Sn_dev_explaineds, "Train Dev Explained", sdd)
	plot_composite_AIC_slice(composite_dev_explaineds, "Test Dev Explained", sdd)

	cd(ret_dir)
	return by_slice_dfs,by_slice_composite_ths, compositeAICs,
			compositeAICcs,
			compositeBICs,
			composite_Sn_accuracys,
			composite_accuracy_tests,
			composite_Sn_dev_explaineds,
			composite_dev_explaineds
end

function plot_composite_AIC_slice(slices, metricname, compositesavepath)
	ret_dir = pwd()
	# pull out all the data for each slice and then get an average
	modeldata = [[] for _=1:length(slices[1])]
	for slice in 1:length(slices)
		for model in 1:length(slices[slice])
			push!(modeldata[model], slices[slice][model])
		end
	end
	f = figure(figsize=(5,3))
	for model = 1:length(slices[1])
		plot(model.*vec(ones(size(modeldata[model]))), modeldata[model], "k.", markersize=10)
		plot(model, nanmean(modeldata[model]), "r.", markersize=20)
	end
	title(join(["composite ", metricname, " across slices"]))
	xlabel("model #")
	ylabel(metricname)
	xticks(1:length(modeldata))
	xlim([0, length(modeldata)+1])
	cd(compositesavepath)
    printFigure(join(["composite_", metricname, "_across_slices_and_sesssions"]); fig=f,figurePath=compositesavepath)
    cd(ret_dir)  
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
	            	# println(join(["SHOULD ALWAYS GET THIS for LickState. names(df)[col] =", names(df)[col], " unique(df[tidx, col])=", unique(df[tidx, col])]))
	            	# what does this actually do? We want to preseve a bool! What the heck?
	                # cc = maximum(df[tidx, col])
	                # let's instead return a bool
	                # println("Number of lick states this trial: ", length(findall(x->x, df[tidx, col])))
	                cc = length(findall(x->x==true, df[tidx, col])) == 1
	                # println("returning: ", cc)
                elseif typeof(df[!,col][1])<:Bool
                	# need to try to convert to float.......
                	if length(unique(df[tidx, col])) !=1 # if there is a 0 and a 1 we better be looking at LickState. I can't get this to handle properly...
	                	# warning(join([" This should be LickState, returning true. Name: ", names(df)[col], " found unique(df[tidx, col])=", unique(df[tidx, col])]))
	                	# println("length:", length(findall(x->x==true, df[tidx, col])), "(should be 1)")
	                	if names(df)[col] != :LickState
	                		warning(join([" This should be LickState, returning true. Name: ", names(df)[col], " found unique(df[tidx, col])=", unique(df[tidx, col])]))
	                	end
	                	cc = true
                	elseif names(df)[col] == :LickState
                		cc = length(findall(x->x==true, df[tidx, col])) == 1
                	else
                		# warning(join([" Bool is unique, handling: ", names(df)[col], " as entry #1 - unique(df[tidx, col])=", unique(df[tidx, col])]))
                		cc = Float64(df[tidx[1], col])
                	end
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
	warning(join(["unique of LickState:", unique(binned_dfs[1][:LickState])]))

	return (binned_dfs, binEdges)
end


function plot_th_vs_timeslice(by_slice_composite_ths; savedir=pwd(), inclusionthresh=10.)
	include_poor_models = false
	propagate_error = false
	if include_poor_models
		warning("including poor model fits in ths")
	else
		warning(join(["excluding ths with more than ", inclusionthresh, " CI"]))
	end
	if propagate_error
		warning("propagating error across timeslices. This may not be appropriate (3/6/2021)")
	else
		warning("averaging standard error across timeslices to get the overall ths - this seems approp (3/6/2021)")
	end

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
    f2s = []
    my_suptitles = []
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

        if !include_poor_models
	    	# first check for all the poorly fitting models. don't include thetas from any where one of the std coeffs is > 10
	    	goodidx = vec(ones(nTimeSlices,1))
	    	for iTh = 1:d_model

	    		CImin_i = [CImin[imodel][i_timeslice][iTh] for i_timeslice=1:nTimeSlices]
	    		# println("model: ", imodel, " iTh", iTh," goodidx: ", goodidx)
	    		# println("	 CImin_i: ", CImin_i)

	    		for i_timeslice = 1:nTimeSlices
		    		if abs(CImin_i[i_timeslice]) > inclusionthresh || isnan(CImin_i[i_timeslice])
		    			goodidx[i_timeslice] = 0.
	    			end
    			end
    			# println(" goodidx now: ", goodidx)
	    	end
	    	goodidx = goodidx .== 1.
	    	# println("model: ", imodel, "goodidx final: ", goodidx)
        	goodidx = findall(x->x==true, goodidx)#findall(x->abs(x)<10 && !isnan(x), goodidx)\
        	# println("model: ", imodel, "goodidx final: ", goodidx)
        	# println(" ")
	    	# println(". ")
    	else
        	# warning("here")
            goodidx = 1:length(nTimeSlices)#findall(x->abs(x)<10 && !isnan(x), CImin_i)
            
    	end

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
            #. Let's try with bad models -- a little funny not to. We are doing ave error so propagation shhould correct this
            #
            

            # println("model=", imodel, " iTh=", iTh, " goodidx=", goodidx, " length(goodidx)=", length(goodidx))

            # trying without propagation of error
            if propagate_error
            	(meanTh_c, propagated_se_th_c, CImin_c, CImax_c, mdf_c) = getCompositeTheta(thi[goodidx], se_th_i[goodidx], dof_i[goodidx],; propagate=true)
            else
	            (meanTh_c, propagated_se_th_c, CImin_c, CImax_c, mdf_c) = getCompositeTheta(thi[goodidx], se_th_i[goodidx], dof_i[goodidx],; propagate=false)
            end
            composite_th_summary = DataFrame(thID=join(["th", iTh]), th = meanTh_c, se_th=propagated_se_th_c, CImin=CImin_c, CImax=CImax_c, mdof=mdf_c, nslices=length(goodidx))
            append!(composite_th_summaries,composite_th_summary)
            
        end
        figname = join(["thbytimeslice_Model_", imodel, "_", timestamp_now(), ".eps"])
        f.savefig(figname, transparent=true, format="eps", bbox_inches="tight",bbox_extra_artists=[my_suptitle])
        
        
        #
        # Don't use any bad models, ie where SE is > 10
        #	Actually, trying now with those bad models
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
        push!(f2s, f2)
        push!(my_suptitles, my_suptitle)
    end
    set_yaxes_same_scale(ax)
    set_yaxes_same_scale(axs2)
    
    for i=1:length(f2s)
    	figname = join(["th_across_timeslices_Model_", i, "_", timestamp_now(), ".eps"])
    	f2s[i].savefig(figname, transparent=true, format="eps", bbox_inches="tight",bbox_extra_artists=my_suptitles[i])
	end
	cd(ret_dir)
end












# models that use HISTORY
function btlogit_100ms_hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["btlogit_100ms_hx_pkg",runID])
	history_spacing_s = 0.1
	if getpackagename
		return packagename
	end
	modelNames = [
	    "DA",
	    "DA_hx0-100ms",
	    "DA_hx0-200ms",
	    "DA_hx0-300ms",
	    "DA_hx0-400ms",
	    "DA_hx0-500ms",
	    "DA_hx0-600ms",
	    "DA_hx0-700ms",
	    "DA_hx0-800ms",
	    "DA_hx0-900ms",
	    "DA_hx0-1000ms",
	]
# Call the runner
	result = btlogit_hxrunner(path; sessionID =sessionID, getpackagename=getpackagename, 
		runID=runID, suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename)
	return result
end
function btlogit_200ms_hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["btlogit_200ms_hx_pkg",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "DA",
	    "DA_hx0-200ms",
	    "DA_hx0-400ms",
	    "DA_hx0-600ms",
	    "DA_hx0-800ms",
	    "DA_hx0-1000ms",
	    "DA_hx0-1200ms",
	    "DA_hx0-1400ms",
	    "DA_hx0-1600ms",
	    "DA_hx0-1800ms",
	    "DA_hx0-2000ms",
	]
# Call the runner
	result = btlogit_hxrunner(path; sessionID =sessionID, getpackagename=getpackagename, runID=runID, 
		suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename)
	return result
end
function btlogit_250ms_hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["btlogit_250ms_hx_pkg",runID])
	history_spacing_s = 0.25
	if getpackagename
		return packagename
	end
	modelNames = [
	    "DA",
	    "DA_hx0-250ms",
	    "DA_hx0-500ms",
	    "DA_hx0-750ms",
	    "DA_hx0-1000ms",
	    "DA_hx0-1250ms",
	    "DA_hx0-1500ms",
	    "DA_hx0-1750ms",
	    "DA_hx0-2000ms",
	    "DA_hx0-2250ms",
	    "DA_hx0-2500ms",
	]
# Call the runner
	result = btlogit_hxrunner(path; sessionID =sessionID, getpackagename=getpackagename, 
		runID=runID, suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename)
	return result
end

function btlogit_hxrunner(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false, history_spacing_s=0.0, modelNames=[], packagename="")
	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s)
	#
	# next, we want to specify our model formulae, including the nested models
	#
	formulas = [
	    @formula(LickState ~ Y),
	    @formula(LickState ~ Y + Hx1),
	    @formula(LickState ~ Y + Hx1 + Hx2),
	    @formula(LickState ~ Y + Hx1 + Hx2 + Hx3),
	    @formula(LickState ~ Y + Hx1 + Hx2 + Hx3 + Hx4),
	    @formula(LickState ~ Y + Hx1 + Hx2 + Hx3 + Hx4 + Hx5),
	    @formula(LickState ~ Y + Hx1 + Hx2 + Hx3 + Hx4 + Hx5 + Hx6),
	    @formula(LickState ~ Y + Hx1 + Hx2 + Hx3 + Hx4 + Hx5 + Hx6 + Hx7),
	  	@formula(LickState ~ Y + Hx1 + Hx2 + Hx3 + Hx4 + Hx5 + Hx6 + Hx7 + Hx8),
	  	@formula(LickState ~ Y + Hx1 + Hx2 + Hx3 + Hx4 + Hx5 + Hx6 + Hx7 + Hx8 + Hx9),
	  	@formula(LickState ~ Y + Hx1 + Hx2 + Hx3 + Hx4 + Hx5 + Hx6 + Hx7 + Hx8 + Hx9 + Hx10),

		]

	

	
	results = modelSelectionByAICBICxval(ndf, :LickState, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# this is already handled by the modelSelectionByAICBICxval function

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end



function btlogit_noHx1_hxrunner(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false, history_spacing_s=0.0, modelNames=[], packagename="")
	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s)
	#
	# next, we want to specify our model formulae, including the nested models
	#
	formulas = [
	    @formula(LickState ~ Y),
	    @formula(LickState ~ Y + Hx2),
	    @formula(LickState ~ Y + Hx2 + Hx3),
	    @formula(LickState ~ Y + Hx2 + Hx3 + Hx4),
	    @formula(LickState ~ Y + Hx2 + Hx3 + Hx4 + Hx5),
	    @formula(LickState ~ Y + Hx2 + Hx3 + Hx4 + Hx5 + Hx6),
	    @formula(LickState ~ Y + Hx2 + Hx3 + Hx4 + Hx5 + Hx6 + Hx7),
	  	@formula(LickState ~ Y + Hx2 + Hx3 + Hx4 + Hx5 + Hx6 + Hx7 + Hx8),
	  	@formula(LickState ~ Y + Hx2 + Hx3 + Hx4 + Hx5 + Hx6 + Hx7 + Hx8 + Hx9),
	  	@formula(LickState ~ Y + Hx2 + Hx3 + Hx4 + Hx5 + Hx6 + Hx7 + Hx8 + Hx9 + Hx10),

		]

	

	
	results = modelSelectionByAICBICxval(ndf, :LickState, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# this is already handled by the modelSelectionByAICBICxval function

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end


function btlogit_50ms_hx_noHx1_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["btlogit_50ms_hx_noHx1",runID])
	history_spacing_s = 0.05
	if getpackagename
		return packagename
	end
	modelNames = [
	    "DA",
	    "DA_hx50-100ms",
	    "DA_hx50-150ms",
	    "DA_hx50-200ms",
	    "DA_hx50-250ms",
	    "DA_hx50-300ms",
	    "DA_hx50-350ms",
	    "DA_hx50-400ms",
	    "DA_hx50-450ms",
	    "DA_hx50-500ms",
	]
# Call the runner
	result = btlogit_noHx1_hxrunner(path; sessionID =sessionID, getpackagename=getpackagename, 
		runID=runID, suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename)
	return result
end
function btlogit_100ms_hx_noHx1_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["btlogit_100ms_hx_noHx1",runID])
	history_spacing_s = 0.1
	if getpackagename
		return packagename
	end
	modelNames = [
	    "DA",
	    "DA_hx100-200ms",
	    "DA_hx100-300ms",
	    "DA_hx100-400ms",
	    "DA_hx100-5000ms",
	    "DA_hx100-600ms",
	    "DA_hx100-700ms",
	    "DA_hx100-800ms",
	    "DA_hx100-900ms",
	    "DA_hx100-1000ms",
	]
# Call the runner
	result = btlogit_noHx1_hxrunner(path; sessionID =sessionID, getpackagename=getpackagename, runID=runID, 
		suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename)
	return result
end
function btlogit_200ms_hx_noHx1_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["btlogit_200ms_hx_noHx1",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "DA",
	    "DA_hx200-400ms",
	    "DA_hx200-600ms",
	    "DA_hx200-800ms",
	    "DA_hx200-1000ms",
	    "DA_hx200-1200ms",
	    "DA_hx200-1400ms",
	    "DA_hx200-1600ms",
	    "DA_hx200-1800ms",
	    "DA_hx200-2000ms",
	]
# Call the runner
	result = btlogit_noHx1_hxrunner(path; sessionID =sessionID, getpackagename=getpackagename, 
		runID=runID, suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename)
	return result
end









# nested history models....
function nestlogit_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["nestlogit_200hx",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "DA",
	    "Hx2s",
	    "Hx2s_1-6s",
	    "Hx2s_1-4s",
	    "Hx2s_1-2s",
	    "Hx2s_1s",
	    "Hx2s_-8s",
	    "Hx2s_-6s",
	    "Hx2s_-4s",
	    "Hx2s_-2s",
	    "Hx2s_-2s_DA",
	]
	formulas = [
		@formula(LickState ~ Y),
	    @formula(LickState ~ Hx10),
	    @formula(LickState ~ Hx10 + Hx9),
	    @formula(LickState ~ Hx10 + Hx9 + Hx8),
	    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7),
	    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6),
	    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5),
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4),
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3),
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2),
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		]

		println("HERE")
# Call the runner
	result = nestlogitrunner(path; sessionID =sessionID, getpackagename=getpackagename, 
		runID=runID, suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename, formulas=formulas)

	return result
end
function nestlogit_allpred_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	warning("package version 3-3-2021 for use with updated error prop and hazard type analyses")
# name the package and runID
	packagename = join(["nestlogit_allpred_200hx",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	

	formulas = [
			@formula(LickState ~ LickTime_2back),
		    @formula(LickState ~ LickTime_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back),
		    @formula(LickState ~ Rxn_1back + Early_1back + Reward_1back + ITI_1back),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back + Y),
		    @formula(LickState ~ Rxn_2back + Early_2back + Reward_2back + ITI_2back + Rxn_1back + Early_1back + Reward_1back + ITI_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Rxn_2back + Early_2back + Reward_2back + ITI_2back +Rxn_1back + Early_1back + Reward_1back + ITI_1back + Y),
		    @formula(LickState ~ LickTime_2back + LickTime_1back + Rxn_2back + Early_2back + Reward_2back + ITI_2back +Rxn_1back + Early_1back + Reward_1back + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		    @formula(LickState ~ Y),
		    @formula(LickState ~ Hx10),
		    @formula(LickState ~ Hx10 + Hx9),
		    @formula(LickState ~ Hx10 + Hx9 + Hx8),
		    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7),
		    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6),
		    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5),
		  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4),
		  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3),
		  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2),
		  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
			]

		modelNames = [
			"Lt2b",
		    "Lt1b",
		    "Lt1b-Lt2b",
		    "Lt1b-Lt2b-DA",
		    "Lt12b-Hx0-2-DA",
		    "oc2b",
		    "oc1b",
		    "oc1b-oc2b",
		    "Lt12b-oc12b",
		    "oc1b-oc2b-DA",
		    "oc12b-Hx0-2-DA",
		    "Lt12b-oc12b-DA",
		    "Lt12b-oc12b-Hx0-2-DA",
		    "DA",
		    "Hx2s",
		    "Hx2s_1-6s",
		    "Hx2s_1-4s",
		    "Hx2s_1-2s",
		    "Hx2s_1s",
		    "Hx2s_-8s",
		    "Hx2s_-6s",
		    "Hx2s_-4s",
		    "Hx2s_-2s",
		    "Hx2s_-2s_DA",
		]


# Call the runner
	result = nestlogitrunner(path; sessionID =sessionID, getpackagename=getpackagename, 
		runID=runID, suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename, formulas=formulas)
	return result
end









function nestlogitrunner(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false, 
	history_spacing_s=0.0, modelNames=[], packagename="", formulas=[], omit_cue=false)	
	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s,omit_cue=omit_cue)	
	
	results = modelSelectionByAICBICxval(ndf, :LickState, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)

# Save each variable to our results folder
	# Now let's gather up the coefficients of our model and the formulas
	# println(typeof(results))
	# WE DONT DO THIS FOR NESTED -- TOO MANY MODELS TO CARRY AROUND
	# yID = :LickState
	# modelData = DataFrame(yID = [yID], predictors=[predictors], th_means=[results.th_summary[1].composite_th], df=[ndf])	
	# modelData_result = [DataFrame() for _=1:nrow(results)]
	# modelData_result[1] = modelData
	# results[:modelData] = modelData_result
# Save each variable to our results folder
	# this is already handled by the modelSelectionByAICBICxval function

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end



function fullmodel_logit_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	#
	#. This model runs just one model type and returns the df and model data so we can sim data from it.
	#
# name the package and runID
	packagename = join(["fullmodel_logit_200hx",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "Hx2s_-2s_DA",
	]
	yID = :LickState
	predictors = [:Hx10,
					:Hx9,
					:Hx8, 
					:Hx7, 
					:Hx6,
					:Hx5,
					:Hx4, 
					:Hx3,
					:Hx2,
					:Y,
					]
	formulas = [
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		]

	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s)	
	
	results = modelSelectionByAICBICxval(ndf, yID, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# Now let's gather up the coefficients of our model and the formulas
	println(typeof(results))
	modelData = DataFrame(yID = [yID], predictors=[predictors], th_means=[results.th_summary[1].composite_th], df=[ndf])	
	modelData_result = [DataFrame() for _=1:nrow(results)]
	modelData_result[1] = modelData
	results[:modelData] = modelData_result

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end


function fullmodel_cue700_logit_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	#
	#. This model runs just one model type and returns the df and model data so we can sim data from it.
	#
# name the package and runID
	packagename = join(["fullmodel_c0-7_logit_200hx",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "Hx2s_-2s_DA",
	]
	yID = :LickState
	predictors = [:Hx10,
					:Hx9,
					:Hx8, 
					:Hx7, 
					:Hx6,
					:Hx5,
					:Hx4, 
					:Hx3,
					:Hx2,
					:Y,
					]
	formulas = [
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		]

	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s, omit_cue=true)	
	
	results = modelSelectionByAICBICxval(ndf, yID, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# Now let's gather up the coefficients of our model and the formulas
	println(typeof(results))
	modelData = DataFrame(yID = [yID], predictors=[predictors], th_means=[results.th_summary[1].composite_th], df=[ndf])	
	modelData_result = [DataFrame() for _=1:nrow(results)]
	modelData_result[1] = modelData
	results[:modelData] = modelData_result

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end



function nestlogit_cue700_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
# name the package and runID
	packagename = join(["nestlogit_c0-7_200hx",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "Hx2s",
	    "Hx2s_1-6s",
	    "Hx2s_1-4s",
	    "Hx2s_1-2s",
	    "Hx2s_1s",
	    "Hx2s_-8s",
	    "Hx2s_-6s",
	    "Hx2s_-4s",
	    "Hx2s_-2s",
	    "Hx2s_-2s_DA",
	]
	formulas = [
	    @formula(LickState ~ Hx10),
	    @formula(LickState ~ Hx10 + Hx9),
	    @formula(LickState ~ Hx10 + Hx9 + Hx8),
	    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7),
	    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6),
	    @formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5),
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4),
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3),
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2),
	  	@formula(LickState ~ Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		]

# Call the runner
	result = nestlogitrunner(path; sessionID =sessionID, getpackagename=getpackagename, 
		runID=runID, suppressFigures=suppressFigures, history_spacing_s = history_spacing_s, 
			modelNames=modelNames, packagename=packagename, formulas=formulas, omit_cue=true)

	return result
end


function fullmodel_cue700_timeinsesh_logit_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	#
	#. This model runs just one model type and returns the df and model data so we can sim data from it.
	#
# name the package and runID
	packagename = join(["fullmodel_cue700_timeinsesh_logit_200hx",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "TinSesh_Hx2s_-2s_DA",
	]
	yID = :LickState
	predictors = [:DataID,
					:Hx10,
					:Hx9,
					:Hx8, 
					:Hx7, 
					:Hx6,
					:Hx5,
					:Hx4, 
					:Hx3,
					:Hx2,
					:Y,
					]
	formulas = [
	  	@formula(LickState ~ DataID + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		]

	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s, omit_cue=true)	
	
	results = modelSelectionByAICBICxval(ndf, yID, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# Now let's gather up the coefficients of our model and the formulas
	println(typeof(results))
	modelData = DataFrame(yID = [yID], predictors=[predictors], th_means=[results.th_summary[1].composite_th], df=[ndf])	
	modelData_result = [DataFrame() for _=1:nrow(results)]
	modelData_result[1] = modelData
	results[:modelData] = modelData_result

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end

function fullmodel_timeinsesh_logit_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	#
	#. This model runs just one model type and returns the df and model data so we can sim data from it.
	#
# name the package and runID
	packagename = join(["fullmodel_timeinsesh_logit_200hx",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "TinSesh_Hx2s_-2s_DA",
	]
	yID = :LickState
	predictors = [:DataID,
					:Hx10,
					:Hx9,
					:Hx8, 
					:Hx7, 
					:Hx6,
					:Hx5,
					:Hx4, 
					:Hx3,
					:Hx2,
					:Y,
					]
	formulas = [
	  	@formula(LickState ~ DataID + Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		]

	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s, omit_cue=false)	
	
	results = modelSelectionByAICBICxval(ndf, yID, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# Now let's gather up the coefficients of our model and the formulas
	println(typeof(results))
	modelData = DataFrame(yID = [yID], predictors=[predictors], th_means=[results.th_summary[1].composite_th], df=[ndf])	
	modelData_result = [DataFrame() for _=1:nrow(results)]
	modelData_result[1] = modelData
	results[:modelData] = modelData_result

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end


function fullmodel_ocs_logit_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	#
	#. This model runs just one model type and returns the df and model data so we can sim data from it.
	#
# name the package and runID
	packagename = join(["fm_ocs_logit_200hx",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "OCS_Hx2s_-2s_DA",
	]
	yID = :LickState
	predictors = [:LickTime_2back,
					:Rxn_2back,
					:Early_2back,
					:Reward_2back,
					:ITI_2back,
					:LickTime_1back,
					:Rxn_1back,
					:Early_1back,
					:Reward_1back,
					:ITI_1back,
					:Hx10,
					:Hx9,
					:Hx8, 
					:Hx7, 
					:Hx6,
					:Hx5,
					:Hx4, 
					:Hx3,
					:Hx2,
					:Y,
					]
	formulas = [
	  	@formula(LickState ~ LickTime_2back + 
					Rxn_2back + 
					Early_2back +
					Reward_2back +
					ITI_2back +
					LickTime_1back +
					Rxn_1back +
					Early_1back +
					Reward_1back +
					ITI_1back + 
					Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		]

	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s, omit_cue=false)	
	
	results = modelSelectionByAICBICxval(ndf, yID, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# Now let's gather up the coefficients of our model and the formulas
	println(typeof(results))
	modelData = DataFrame(yID = [yID], predictors=[predictors], th_means=[results.th_summary[1].composite_th], df=[ndf])	
	modelData_result = [DataFrame() for _=1:nrow(results)]
	modelData_result[1] = modelData
	results[:modelData] = modelData_result

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end


function DAHxmodel_logit_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	#
	#. This model runs just one model type and returns the df and model data so we can sim data from it.
	#
# name the package and runID
	packagename = join(["DAHxmodel_logit_200hx_pkg",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "DAHx_-2s_DA",
	]
	yID = :LickState
	predictors = [:Hx10,
					:Hx9,
					:Hx8, 
					:Hx7, 
					:Hx6,
					:Hx5,
					:Hx4, 
					:Hx3,
					:Hx2,
					:Y,
					]
	formulas = [
	  	@formula(LickState ~ 
					Hx10 + Hx9 + Hx8 + Hx7 + Hx6 + Hx5 + Hx4 + Hx3 + Hx2 + Y),
		]

	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s, omit_cue=false)	
	
	results = modelSelectionByAICBICxval(ndf, yID, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# Now let's gather up the coefficients of our model and the formulas
	println(typeof(results))
	modelData = DataFrame(yID = [yID], predictors=[predictors], th_means=[results.th_summary[1].composite_th], df=[ndf])	
	modelData_result = [DataFrame() for _=1:nrow(results)]
	modelData_result[1] = modelData
	results[:modelData] = modelData_result

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end
function SHUFFLE_DAHxmodel_logit_200hx_pkg(path; sessionID ="", getpackagename=false, runID=0, suppressFigures=false)
	#
	#. This model runs just one model type and returns the df and model data so we can sim data from it.
	#
# name the package and runID
	packagename = join(["SHUFFLE_DAHxmodel_logit_200hx_pkg",runID])
	history_spacing_s = 0.2
	if getpackagename
		return packagename
	end
	modelNames = [
	    "SHUFFLE_DAHx_-2s_DA",
	]
	yID = :LickState
	predictors = [:Hx10Shuffle,
					:Hx9Shuffle,
					:Hx8Shuffle, 
					:Hx7Shuffle, 
					:Hx6Shuffle,
					:Hx5Shuffle,
					:Hx4Shuffle, 
					:Hx3Shuffle,
					:Hx2Shuffle,
					:Yshuffle,
					]
	formulas = [
	  	@formula(LickState ~ 
					Hx10Shuffle + Hx9Shuffle + Hx8Shuffle + Hx7Shuffle + Hx6Shuffle + Hx5Shuffle + Hx4Shuffle + Hx3Shuffle + Hx2Shuffle + Yshuffle),
		]

	if history_spacing_s == 0.0
		error("Improper history_spacing_s! Check the calling fxn, it should specify nonzero")
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
	ndf = extract_data_with_baselineandLOI(path; normalize=true, useHx=true, history_spacing_s=history_spacing_s, omit_cue=false)
	#
	# Add shuffle rows to the dataframe
	#
	ndf[:Hx10Shuffle] = Random.shuffle(ndf[:Hx10])
	ndf[:Hx9Shuffle] = Random.shuffle(ndf[:Hx9])
	ndf[:Hx8Shuffle] = Random.shuffle(ndf[:Hx8])
	ndf[:Hx7Shuffle] = Random.shuffle(ndf[:Hx7])
	ndf[:Hx6Shuffle] = Random.shuffle(ndf[:Hx6])
	ndf[:Hx5Shuffle] = Random.shuffle(ndf[:Hx5])
	ndf[:Hx4Shuffle] = Random.shuffle(ndf[:Hx4])
	ndf[:Hx3Shuffle] = Random.shuffle(ndf[:Hx3])
	ndf[:Hx2Shuffle] = Random.shuffle(ndf[:Hx2])

	
	results = modelSelectionByAICBICxval(ndf, yID, formulas, modelNames, "logit"; 
    		n_iters=100,updownsampleYID=true, figurePath=figurePath, savePath = savepath, suppressFigures=suppressFigures)
# Save each variable to our results folder
	# Now let's gather up the coefficients of our model and the formulas
	println(typeof(results))
	modelData = DataFrame(yID = [yID], predictors=[predictors], th_means=[results.th_summary[1].composite_th], df=[ndf])	
	modelData_result = [DataFrame() for _=1:nrow(results)]
	modelData_result[1] = modelData
	results[:modelData] = modelData_result

# make a working result df with all the results to keep in workspace
	result = results
	return result#, ndf
end

