#
#  The following functions allow us to use GLMs and report them
#
#	sigmoid(x)
#	normalize_vector_0_1(v; mx=nothing, mn=nothing) -- vector of vectors
# 	zero_to_one(signal) -- single vector
# 	normalize_to_ref_signal(signal, ref_signal)
#	build_and_report_logit_model(Formula, train, test; modelName="model", modelClass="logit")
# 	build_and_report_linear_model(Formula, train, test; modelName="model")
#.  plot_fit_results(X::Vector, Y::Vector, Yfit, correct; modelClass="logit")
# 	plot_residuals_by_trial
#	select_single_timepoint_per_trial_dataset
#	get_model_stats -- now works for both Julia GLM and for my own ridge regression 
#
# 	nx_kcomparisons_GLM
# 	plot_k_xval_stat
# 	wrap
# 	collect_theta
# 	unwrap_theta
# 	theta_summary
#
#	kfold_xval
#	plot_xval_loss
# 	modelSelectionByAICBICxval
# 	compare_AICBIC
# 	getBootCI
#
using ClassImbalance
import Random
using MLBase
using Lathe
using GLM
using ROCAnalysis
using Lathe.preprocess: TrainTestSplit
using LinearAlgebra
using Distributions
using HypothesisTests

# function sigmoid(x)
#     return (exp.(x))./(1 .+ exp.(x))
# end

function normalize_vector_0_1(v; mx=nothing, mn=nothing)
    # get the vector into one long vector
    #
    if isnothing(mx)
        v = reduce(vcat, v)
        mx = maximum(v)
        mn = minimum(v)
    end
    return (v.-mn)./(mx-mn)
end
function normalize_to_ref_signal(signal, ref_signal)
    println(size(ref_signal))
    l = minimum(ref_signal);
    u = maximum(ref_signal);
    # normalize signal
    zsig = zero_to_one(signal);
    nsig = zsig.*u.+l
end;
function zero_to_one(signal)
	#
	#	for a single vector
	#
    l = minimum(signal);
    u = maximum(signal);
    normd = (signal.-l)./u
end;

function MSELoss(y,yFit)
	return 1/length(y)*sum((y .- yFit).^2)
end

# function build_and_report_logit_model(Formula, train, test; modelName="model", modelClass="logit")
# !!!!!!!!!! SEEMS TO HAVE AN ERROR??

#     # use @formula(Y_name ~ X_1, X_2...) to make the formula
#     # train the model
#     if modelClass == "logit"
#         logit_model = glm(Formula, train, Binomial(), ProbitLink())
#     elseif modelClass == "normal"
#         logit_model = glm(Formula, train, Normal())
#     end
    
   
#     # Predict the target variable on training data 
#     prediction_Sn = predict(logit_model,train)
#     # Convert probability score to class
#     prediction_class_Sn = [if x < 0.5 0 else 1 end for x in prediction_Sn];
#     prediction_df_Sn = DataFrame(y_actual = train.LickState, y_predicted = prediction_class_Sn, prob_predicted = prediction_Sn);
#     prediction_df_Sn.correctly_classified = prediction_df_Sn.y_actual .== prediction_df_Sn.y_predicted
#     # Accuracy Score
#     accuracy_Sn = mean(prediction_df_Sn.correctly_classified)
#     println("Accuracy of the training model is : ",accuracy_Sn)
#     # Confusion Matrix
#     confusion_matrix_Sn = MLBase.roc(prediction_df_Sn.y_actual, prediction_df_Sn.y_predicted)
#     println(confusion_matrix_Sn)
#     # Predict the target variable on test data 
#     prediction = predict(logit_model,test)
#     # Convert probability score to class
#     prediction_class = [if x < 0.5 0 else 1 end for x in prediction];
#     prediction_df = DataFrame(y_actual = test.LickState, y_predicted = prediction_class, prob_predicted = prediction);
#     prediction_df.correctly_classified = prediction_df.y_actual .== prediction_df.y_predicted
#     # Accuracy Score
#     accuracy = mean(prediction_df.correctly_classified)
#     println("Accuracy of the test model is : ",accuracy)

#     # Confusion Matrix
#     confusion_matrix = MLBase.roc(prediction_df.y_actual, prediction_df.y_predicted)
#     println(confusion_matrix)


# 	if modelClass == "logit"
#         logit_model_test = glm(Formula, test, Binomial(), ProbitLink()) 
#     elseif modelClass == "normal"
#         logit_model_test = glm(Formula, test, Normal()) 
#     end

#     model_results = DataFrame(fit = ["train", "test"], 
#         model = [logit_model, logit_model], 
#         formula = [formula, formula],
#         X = [logit_model.mf.data[2:end], logit_model_test.mf.data[2:end]],
#         y = [logit_model.mf.data[1], logit_model_test.mf.data[1]],
#         prediction = [prediction_Sn, prediction], 
#         prediction_df = [prediction_df_Sn, prediction_df],
#         accuracy = [accuracy_Sn, accuracy],
#         confusion_matrix = [confusion_matrix_Sn, confusion_matrix])
    

    
#     for i = 2:length(logit_model.mf.data)
#         outputname = Formula.lhs
#         if length(logit_model.mf.data) > 2
#             predictorname = Formula.rhs[i-1]
#         else
#             predictorname = Formula.rhs
#         end
#         figure(figsize=(8,3))
#         X = logit_model.mf.data[i]
#         Y = logit_model.mf.data[1]
#         Yfit = prediction_Sn
#         correct = prediction_df_Sn.correctly_classified
#         Accuracy = accuracy_Sn
#         subplot(1,2,1)
#         plot_fit_results(X,Y,Yfit,correct)
#         title(join(["training fit, accuracy=", round(Accuracy, digits=2)]))
#         ax=gca()
#         ax.set_xlabel(predictorname)
#         ax.set_ylabel(outputname)


#         X = logit_model_test.mf.data[i]
#         Y = logit_model_test.mf.data[1]
#         Yfit = prediction
#         correct = prediction_df.correctly_classified
#         Accuracy = accuracy
#         subplot(1,2,2)
#         plot_fit_results(X,Y,Yfit,correct)
#         title(join(["test fit, accuracy=", round(Accuracy, digits=2)]))
#         ax=gca()
#         ax.set_xlabel(predictorname)
#         printFigure(join([modelName, "Results_", i-1]); fig=gcf())
#     end
#     println("TRAINING")
#     println(logit_model)
#     println("TEST")
#     println(logit_model_test)
#     return (logit_model, logit_model_test, model_results)
# end
function build_and_report_logit_model(Formula, train, test; modelName="model", modelClass="logit", verbose=true, figurePath=".")
    # use @formula(Y_name ~ X_1, X_2...) to make the formula
    # train the model
    if modelClass == "logit"
        logit_model = glm(Formula, train, Binomial(), ProbitLink())
    elseif modelClass == "normal"
        logit_model = glm(Formula, train, Normal())
    end
    
   
    # Predict the target variable on training data 
    prediction_Sn = StatsBase.predict(logit_model,train)
    # Convert probability score to class
    prediction_class_Sn = [if x < 0.5 0 else 1 end for x in prediction_Sn];
    prediction_df_Sn = DataFrame(y_actual = train.LickState, y_predicted = prediction_class_Sn, prob_predicted = prediction_Sn);
    prediction_df_Sn.correctly_classified = prediction_df_Sn.y_actual .== prediction_df_Sn.y_predicted
    # Accuracy Score
    accuracy_Sn = mean(prediction_df_Sn.correctly_classified)
    if verbose
	    println("Accuracy of the training model is : ",accuracy_Sn)
    end
    # Confusion Matrix
    confusion_matrix_Sn = MLBase.roc(prediction_df_Sn.y_actual, prediction_df_Sn.y_predicted)
    if verbose
    	println(confusion_matrix_Sn)
	end
    # Predict the target variable on test data 
    prediction = StatsBase.predict(logit_model,test)
    # Convert probability score to class
    prediction_class = [if x < 0.5 0 else 1 end for x in prediction];
    prediction_df = DataFrame(y_actual = test.LickState, y_predicted = prediction_class, prob_predicted = prediction);
    prediction_df.correctly_classified = prediction_df.y_actual .== prediction_df.y_predicted
    # Accuracy Score
    accuracy = mean(prediction_df.correctly_classified)
    if verbose
    	println("Accuracy of the test model is : ",accuracy)
	end

    # Confusion Matrix
    confusion_matrix = MLBase.roc(prediction_df.y_actual, prediction_df.y_predicted)
    if verbose
    	println(confusion_matrix)
	end

    model_results = DataFrame(fit = ["train", "test"], 
        model = [logit_model, logit_model], 
        formula = [formula, formula],
        prediction = [prediction_Sn, prediction], 
        prediction_df = [prediction_df_Sn, prediction_df],
        accuracy = [accuracy_Sn, accuracy],
        confusion_matrix = [confusion_matrix_Sn, confusion_matrix])
    

    if modelClass == "logit"
        logit_model_test = glm(Formula, test, Binomial(), ProbitLink()) 
    elseif modelClass == "normal"
        logit_model_test = glm(Formula, test, Normal()) 
    end
    if verbose
	    for i = 2:length(logit_model.mf.data)
	        outputname = Formula.lhs
	        if length(logit_model.mf.data) > 2
	            predictorname = Formula.rhs[i-1]
	        else
	            predictorname = Formula.rhs
	        end
	        figure(figsize=(8,3))
	        X = logit_model.mf.data[i]
	        Y = logit_model.mf.data[1]
	        Yfit = prediction_Sn
	        correct = prediction_df_Sn.correctly_classified
	        Accuracy = accuracy_Sn
	        subplot(1,2,1)
	        plot_fit_results(X,Y,Yfit,correct=correct)
	        title(join(["training fit, accuracy=", round(Accuracy, digits=2)]))
	        ax=gca()
	        ax.set_xlabel(predictorname)
	        ax.set_ylabel(outputname)


	        X = logit_model_test.mf.data[i]
	        Y = logit_model_test.mf.data[1]
	        Yfit = prediction
	        correct = prediction_df.correctly_classified
	        Accuracy = accuracy
	        subplot(1,2,2)
	        plot_fit_results(X,Y,Yfit,correct=correct)
	        title(join(["test fit, accuracy=", round(Accuracy, digits=2)]))
	        ax=gca()
	        ax.set_xlabel(predictorname)
	        printFigure(join([modelName, "Results_", i-1]); fig=gcf(), figurePath=figurePath)
	    end
    end
    if verbose
	    println("TRAINING")
	    println(logit_model)
	    println("TEST")
	    println(logit_model_test)
    end
    return (logit_model, logit_model_test, accuracy_Sn, accuracy)
end


function build_and_report_linear_model(Formula, train, test; modelName="model", verbose=true)
    # use @formula(Y_name ~ X_1, X_2...) to make the formula
    # train the model
    model_train = glm(Formula, train, Normal()) 
    model_test = glm(Formula, test, Normal())    
   
    # Predict the target variable on training data 
    prediction_Sn = StatsBase.predict(model_train,train)
    y_actual_Sn = model_train.mf.data[1]
    prediction_df_Sn = DataFrame(y_actual = y_actual_Sn, y_predicted = prediction_Sn);
    # Rsq
    Rsq_Sn = Rsq(prediction_df_Sn.y_actual, prediction_df_Sn.y_predicted)
    # println("Rsq of the training model is : ",Rsq_Sn)


    # Predict the target variable on test data 
    prediction_test = StatsBase.predict(model_train,test)
	y_actual_test = model_test.mf.data[1]
    prediction_df_test = DataFrame(y_actual = y_actual_test, y_predicted = prediction_test);
    # Rsq
    Rsq_test = Rsq(prediction_df_test.y_actual, prediction_df_test.y_predicted)
    # println("Rsq of the test set is : ",Rsq_test)

    outputname = Formula.lhs
    predictornames = Formula.rhs
    # if length(model_train.mf.data) > 2
    #     predictorname = Formula.rhs[i-1]
    # else
    #     predictorname = Formula.rhs
    # end

    model_results = DataFrame(fit = ["train", "test"],
    	modelName = modelName,
        model = [model_train, model_test], 
        formula = [formula, formula],
        outputname = outputname,
        predictornames = predictornames,
        # X = [model.mf.data[2:end], model_test.mf.data[2:end]],
        y = [y_actual_Sn, y_actual_test],
        prediction = [prediction_Sn, prediction_test], 
        prediction_df = [prediction_df_Sn, prediction_df_test],
        Rsq = [Rsq_Sn, Rsq_test],
        )
    
    if verbose
	    report_linear_model(model_results)
    end
    
    
    return (model_train, model_test, model_results)
end

function report_linear_model(model_results)
	model_train = model_results.model[1]
	model_test = model_results.model[2]
    Formula = model_results.formula[1]
    y_actual_Sn = model_results.y[1]
    y_actual_test = model_results.y[2]
    prediction_Sn = model_results.prediction[1]
    prediction_test = model_results.prediction[2] 
    Rsq_Sn = model_results.Rsq[1]
    Rsq_test = model_results.Rsq[2]
    outputname = model_results.outputname[1]
    predictornames = model_results.predictornames
    modelName = model_results.modelName
	for i = 2:length(model_train.mf.data)
        if length(model_train.mf.data) > 2
	        predictorname = predictornames[1][i-1]
	    else
	        predictorname = predictornames[1]
	    end
        figure(figsize=(16,6))
        X = model_train.mf.data[i]
        Y = y_actual_Sn
        Yfit = prediction_Sn
        ax1 = subplot(2,2,1)
        render_xy(X, Y; t=join(["training fit, Rsq=", round(Rsq_Sn, digits=4)]), xl=predictorname, yl=outputname, ax=ax1)
        ax1.plot(X, Yfit, "r-")
        ax3 = subplot(2,2,3)
		ax3.plot(X, vec(zeros(size(X))), "k-")
		render_xy(X, Y-Yfit; t="", xl=predictorname, yl="Y-Yfit", ax=ax3)
		
        X = model_test.mf.data[i]
        Y = y_actual_test
        Yfit = prediction_test
        ax2 = subplot(2,2,2)
		render_xy(X, Y; t=join(["test fit, Rsq=", round(Rsq_test, digits=4)]), xl=predictorname, yl=outputname, ax=ax2)
        plot(X, Yfit, "r-")        
        title(join(["test fit, Rsq=", round(Rsq_test, digits=2)]))
        set_xaxes_same_scale([ax1, ax2])
        set_yaxes_same_scale([ax1, ax2])
        ax4 = subplot(2,2,4)
		ax4.plot(X, vec(zeros(size(X))), "k-")
		render_xy(X, Y-Yfit; t="", xl=predictorname, yl="Y-Yfit", ax=ax4)

		printFigure(join([modelName, "Results_", i-1]); fig=gcf())
    end
    println("TRAINING")
    println(model_train)
    println("TEST")
    println(model_test)
end


function plot_xval_fit_results(result_df, model_summary; modelClass="logit")
	# extract the relevant data from the result struct
	X = result_df.X[1]
	Y = vec(vcat(result_df.y[1],))
	Yfit = vec(vcat(result_df.yFit[1],))
	lam = result_df.lam[1]
	predictors = result_df.predictors[1]
	predicted = result_df.predicted[1]
	k = result_df.n_kfold_xval[1]
	for i = 2:length(predictors)
		figure(figsize=(2,2))
		ax=gca()
		if modelClass == "logit"
			plot_fit_results(vec(X[i, :]), Y, Yfit; correct=correct, modelClass="logit", ax=ax)
		else
			plot_fit_results(vec(X[i, :]), Y, Yfit; modelClass="linear", ax=ax)
		end
		ax.set_xlabel(predictors[i])
		ax.set_ylabel(predicted)
		ax.set_title(join(["lam=", lam, " k=",k]))
	end
	model_summary
end

function plot_fit_results(X::Vector, Y::Vector, Yfit; correct=nothing, modelClass="logit", ax=gca(), alph=nothing, figurePath=".")
	# for logit model currently
	if isnothing(alph)
		if length(Yfit) > 10000
			alph = 1/(20*log(length(Yfit)))
		else
			alph = 0.2
		end
		# alph = 1/length(x)
	end
    idxs = sortperm(X)

    if modelClass == "logit"
	    correct_idxs = findall(x->x != 0,correct)
	    jitter_correct = (rand(length(correct_idxs)).-0.5)./2
	    incorrect_idxs = findall(x->x == 0, correct)
	    jitter_incorrect = (rand(length(incorrect_idxs)).-0.5)./2

	    ax.plot(X[idxs], Yfit[idxs])
	    ax.plot(X[correct_idxs],Y[correct_idxs].+jitter_correct, "g.", alpha=alph)
	    ax.plot(X[incorrect_idxs],Y[incorrect_idxs].+jitter_incorrect, "r.", alpha=alph)

    elseif modelClass == "linear"
    	if length(Y) > 10000
			alph = 1/(20*log(length(Y)))
		else
			alph = 0.2
		end
	    ax.plot(X[idxs],Y[idxs], "k.", alpha=alph)
	    ax.plot(X[idxs], Yfit[idxs], "r-")

    end
end;


function Rsq(y, yfit)
	ESS = sum((yfit .- mean(y)).^2);
	RSS = sum((yfit .- y).^2);
	return Rsq = ESS/(RSS+ESS)
end


function updownsample(df::DataFrame, refvect::Vector{Bool}, n::Int, downOnly=false)
    # refvect is the bool vector of states so we can up/downsample them
    # n is the number of desired resamples with replacement
    
    idxs_true = findall(x->x, refvect)
    idxs_false = findall(x->!x, refvect)
    
    # resample with replacement the true and false idxs
    if downOnly
    	resampled_idx_true = idxs_true
    else
	    resampled_idx_true = [sample(idxs_true) for x = 1:n]
    end
    resampled_idx_false = [sample(idxs_false) for x = 1:n]
    all_idx = vcat(resampled_idx_true,resampled_idx_false)
    # warning("RBF: check that trial idxs not shuffled...")
    # all_idx = shuffle!(all_idx)
#     Names = names(df)
    newData = similar(df,0)
    # add the new data row by row
    for i = 1:length(all_idx)
        push!(newData, df[all_idx[i],:])
    end
#     for i = 1:length(names(df))
#         og_data = df[!,names(norm_df)[i]]
#         newData[:,i] = vcat(og_data[idxs_true],og_data[idxs_false])
#     end
    return newData
end



# % PURPOSE: computes Cochrane-Orcutt ols Regression for AR1 errors 
	# %--------------------------------------------------- 
	# % USAGE: results = olsc(y,x) 
	# % where: y = dependent variable vector (nobs x 1) 
	# %        x = independent variables matrix (nobs x nvar) 
	# %--------------------------------------------------- 
	# % RETURNS: a structure 
	# %        results.meth  = 'olsc' 
	# %        results.beta  = bhat estimates 
	# %        results.rho   = rho estimate 
	# %        results.tstat = t-stats 
	# %        results.trho  = t-statistic for rho estimate 
	# %        results.yhat  = yhat 
	# %        results.resid = residuals 
	# %        results.sige  = e'*e/(n-k) 
	# %        results.rsqr  = rsquared 
	# %        results.rbar  = rbar-squared 
	# %        results.iter  = niter x 3 matrix of [rho converg iteration#] 
	# %        results.nobs  = nobs 
	# %        results.nvar  = nvars 
	# %        results.y     = y data vector 
	# % -------------------------------------------------- 
	# % SEE ALSO: prt_reg(results), plt_reg(results) 
	# %--------------------------------------------------- 
	 
	# % written by: 
	# % James P. LeSage, Dept of Economics 
	# % University of Toledo 
	# % 2801 W. Bancroft St, 
	# % Toledo, OH 43606 
	# % jpl@jpl.econ.utoledo.edu 
	# Modified for Julia Language by ahamilos, 1/16/2021
function olsc(y,x; lam::Float64=0.0, predictors::Vector{String}=["b1"], predicted="y") 
	# convert to arrays
	x = hcat(x)
	y=hcat(y)
 
	# % do error checking on inputs 
	# if (nargin ~= 2); error('Wrong # of arguments to olsc'); end; 
	warning("check that NO ones as a predictor in x")
	x = hcat(x, ones(size(x)[1]))
	(nobs, nvar) = size(x)[1], size(x)[2]; 
	if size(y)[1] == 1 
		y = y'
	end
	nobs2 = length(y); 
	 
	if nobs != nobs2
		error("x and y must have same # obs in Cochrane-Orcutt ols Regression for AR1 errors (olsc)")
	end 
	
	# check our predictor names match nvar
	if length(predictors)+1 != nvar
		predictors = Vector{String}(undef, 0)
		for i = 1:nvar
			push!(predictors, [join(["th", i])])
		end
		predictors = push!(predictors,"th0")
	else
		predictors = Vector{String}(undef, 0)
		predictors = push!(predictors,"th0")
	end

	 
	# % ----- setup parameters 
	ITERMAX = 100; 
	converg = 1.0; 
	rho = 0.0; 
	iter = 1; 
	# xtmp = x(1:end-1,:)#lag(x,1); 
	# ytmp = y(1:end-1)#lag(y,1); 
	 
	# # % truncate 1st observation to feed the lag 
	# xlag = xtmp(2:nobs,:); 
	# ylag = ytmp(2:nobs,1); 
	# yt = y(2:nobs,1); 
	# xt = x(2:nobs,:); 
	yt = y[2:nobs,1]; 
	ylag = y[1:nobs-1,1];
	xt = x[2:nobs,:]; 
	xlag = x[1:nobs-1,:]; 
	         
	# % setup storage for iteration results 
	iterout = zeros(ITERMAX,3); 
	 
	while (converg > 0.0001) & (iter < ITERMAX) 
		println(" ")
		progressbar(iter,ITERMAX)
		println("Try #", iter)
		println("	OLS of data")
		# % step 1, using intial rho = 0, do OLS to get bhat 
		ystar = yt .- rho.*ylag; 
		xstar = xt .- rho.*xlag; 
		
		(res, _) = ridge_ols(ystar,xstar,lam=lam, predictors=predictors, predicted=predicted); 
		println(join(["		OLS Coeffs this round: ", res.th]))
		println(join(["		se this round: ", res.se_th]))


		e = y .- x*res.th[1]; #'
		# elag = lag(e); 
		elag = e[1:nobs-1,1];
		 
		# % truncate 1st observation to account for the lag 
		et = e[2:nobs,1]; 
		et = hcat(et)
		# elagt = elag(2:nobs,1); 
		elagt = elag;
		elagt = hcat(elag, ones(length(elag)[1]))
		 
		# % step 2, update estimate of rho using residuals 
		# %         from step 1 
		println("	AR[1] OLS of errors")
		(res_rho, _) = ridge_ols(et,elagt); 
		rho_last = rho; 
		println("		Error AR[1] Coeffs this round: ",res_rho.th)
		rho = res_rho.th[1][1]; 
		println("		Error AR[1] slope this round: ",rho)
		println(join(["		se AR[1] this round: ", res_rho.se_th]))
		converg = broadcast(abs, rho .- rho_last);
		 
		iterout[iter,1] = rho; 
		iterout[iter,2] = converg; 
		iterout[iter,3] = iter; 
		 
		iter = iter + 1; 
	 
	end; #% end of while loop 
	 
	if iter == ITERMAX 
		error("ols_corc did not converge in 100 iterations"); 
	end
	 
	# result.iter = iterout[1:iter-1,:]; 
	 
	# % after convergence produce a final set of estimates using rho-value 
	if length(rho) != 1
		warning("I think size of rho should be 1,1: ", size(rho))
	end
	ystar = yt .- rho.*ylag; 
	xstar = xt .- rho.*xlag; 
	 
	(result, model_summary) = ridge_ols(ystar,xstar, lam=lam, predictors=predictors, predicted=predicted); 

	meth = "ridge_olsc"; 
	rho = rho; 
	iter = iterout[1:iter-1,:]; 
	# % compute t-statistic for rho 
	varrho = (1. - rho*rho)/(nobs-2.); 
	trho = rho/sqrt(varrho); 

	df = DataFrame(predicted=[predicted], predictors=[predictors],th=result.th, X=result.X, y=result.y, 
		lam=result.lam, yFit=result.yFit, CVmat=result.CVmat, 
		se_model=result.se_model, se_th=result.se_th, 
		signifCoeff=result.signifCoeff,
		meth=meth,
		rho=rho,
		iter=[iter],
		varrho=varrho,
		trho=trho,
		ESS=result.ESS,RSS=result.RSS,Rsq=result.Rsq, model_summary=[model_summary])

	model_summary[:,:]
	
	return (df, model_summary)
end

function ridge_ols(y,X;lam=0.0, predictors=["b1", "b2"], predicted="y")
	# ensure that X = (dxn), y=(1xn) to start
	if size(X)[2] != length(y)
		X = X' #'
	end
	if size(y)[1] != 1
		y = y' #'
	end
	
	d = size(X)[1]
    th = (X*X'+lam.*I)\X*y'
    yFit = th'*X; #'
    XtX = X*X' #'
	se_model = sqrt(sum((y .- yFit).^2. ./size(yFit)[2]));

	if d > 1
		CVmat = (XtX + (lam * (size(XtX)[1] * I)))^-1*XtX*(XtX + (lam* (size(XtX)[1] * I)))^-1;
		if size(XtX)[1] != length(th)
			error("CVmat unexpected dimensions, RBF")
		end
		se_th = se_model*diag(CVmat).^0.5; 
	else
	    th = [th]
	    XtX = hcat([XtX],)
		CVmat = (XtX + hcat([lam],))^-1*XtX*(XtX + hcat([lam],))^-1;
		se_th = se_model*diag(CVmat).^0.5; 
	end

	
	distFromZero = broadcast(abs, th) .- 2. .*broadcast(abs, se_th);
	signifCoeff = distFromZero .> 0.;

	Resid = y .- yFit;
	std_Resid = sqrt(sum(Resid.^2.)./(length(y) - length(th)));
	std_yActual = std(y);
	# explainedVarianceR2 = 1. - std_Resid^2. /std_yActual^2.; # THIS IS WRONG!!!!!!!!!
	# % 
	# %  Check consistent
	# % 	
	ESS = sum((yFit .- mean(y)).^2.);
	RSS = sum((yFit .- y).^2.);
	Rsq = ESS/(RSS+ESS);
	# if round(explainedVarianceR2, digits=3) != round(Rsq, digits=3)
	# 	warning(join(["expected Rsq: ", Rsq, "\n    Got: ",explainedVarianceR2]))
 #    end

    nSim = 10000
    n = length(y)
    # (CI_th, p_th, sims) = simulate_CI_th(th, se_model, CVmat, n; nSim=nSim)
    stats_df = get_model_stats([]; n=n, th=th, se_th=se_th, dof=n-length(th), th_names=predictors)
    # DataFrame(th_names=th_names, n=n.*ones(length(th)), th=th, se_th=se_th, std_th=std_th, dof=dof*ones(length(th)), 
    #    p_t=p_t, CImin_t=CImin_t, CImax_t=CImax_t, p_z=p_z, CImin_z=CImin_z, CImax_z=CImax_z)

    result_df = DataFrame(predictors=[predictors], predicted=[predicted], th=[th], X=[X], y=[y], lam=lam, yFit=[yFit], 
    	CVmat=[CVmat], se_model=se_model, se_th=[se_th], 
    	signifCoeff=[signifCoeff], 
    	ESS=ESS,RSS=RSS,Rsq=Rsq, n_kfold_xval=NaN)

    if length(predictors) != d
		predictors = []
		for i = 1:d-1
			push!(predictors, [join(["th", i])])
		end
		predictors = push!(predictors,"th0")
	end
	

    # model_summary = DataFrame(predictors=predictors, Coeff=vec(th), StdError=se_th, nBoot=[sims for i=1:d], p=p_th, CI95_lower=CI_th[1:end,1], CI95_upper=CI_th[1:end,2])
    model_summary = DataFrame(predictors=predictors, Coeff=vec(th), StdError=se_th, p=stats_df.p_t, CI95_lower=stats_df.CImin_t, CI95_upper=stats_df.CImax_t)

	return (result_df, model_summary, stats_df)
end

function run_regression(data::DataFrame, predictorIDs::Array, yID; lam=0.0, add_th0=true)
    #
    # predictorIDs and yID should be SYMBOLS matching the column name.
    #     run_regression(ksets[1], [:Y], :Time2Lick)
    #
    # predictorIDs = [x->Symbol(x) for x in predictorIDs]
    yID = Symbol(yID)
    y = data[yID]' #'
    if !isempty(predictorIDs)
	    X = data[predictorIDs[1]]' #'
	    if length(predictorIDs)>1
	        for x = 2:length(predictorIDs)
	            X = vcat(X, data[predictorIDs[x]]') #'
	        end
	    end
    else
    	X = ones(length(y))'
    	add_th0 = false
    	predictorIDs = [:th0]
    end
	predictors = Vector{String}(undef, 0)
	if add_th0
    	push!(predictors,"th0")
        X = vcat(ones(length(X[1,:]))',X)#'
    else
        X = hcat(X,)
	end
    for i=1:length(predictorIDs)
    	push!(predictors,string(predictorIDs[i]))
	end
	
    (result_df, model_summary, stats_df) = run_regression_run(data, X, y; lam=lam, predictors=predictors,add_th0=add_th0, predicted=string(yID))
    return (result_df, model_summary, stats_df)
end
function run_regression(data::DataFrame, predictorID::Symbol, yID; lam=0.0, add_th0=true)
    #
    # predictorIDs and yID should be SYMBOLS matching the column name.
    #     run_regression(ksets[1], [:Y], :Time2Lick)
    #
    predictorIDs = Symbol(predictorID)
    yID = Symbol(yID)
    y = data[yID]' #'
    X = data[predictorID]' #'
    predictors = string(predictorID)
    (result_df, model_summary, stats_df) = run_regression_run(data, X, y; lam=lam, predictors=predictors,add_th0=add_th0, predicted=string(yID))
    return (result_df, model_summary, stats_df)
end
function run_regression_run(data::DataFrame, X, y; lam=0.0, predictors=["b1", "b2"], add_th0=true, predicted="y")
    #
    #  Here we run the actual regression!
    #
    X = hcat(X,)
    y = hcat(y,)

    (result_df, model_summary, stats_df) = ridge_ols(y,X;lam=lam, predictors=predictors, predicted=predicted)
    #result_df = DataFrame(th=[th], X=[X], y=[y], lam=lam, yFit=[yFit], 
#     	CVmat=[CVmat], se_model=se_model, se_th=[se_th], 
#     	signifCoeff=[signifCoeff], 
#     	ESS=ESS,RSS=RSS,Rsq=Rsq)
    #model_summary = DataFrame(predictors=predictors, Coeff=vec(th), 
#         StdError=se_th, nBoot=[sims for i=1:d], p=p_th, CI95_lower=CI_th[1:end,1], 
#         CI95_upper=CI_th[1:end,2])
	# println("Loss: ", MSELoss(y,result_df.yFit[1]))
    return (result_df, model_summary, stats_df)
end



function simulate_CI_th(th, se_model, CVmat, n; nSim=10000)
	# % 
 #    % 	We will redraw coeffs and calculate model error/get distribution of consistent coeffs
 #    % 
 #    % 	sim_th = (d x nSim) matrix
 #    % 
 #    % ----------------------------------------------------
 #    %  	th is a d x 1 vector of coefficients to the model, where model predicts n datapoints (length of y_actual)
 #    % 	sim_th is a dxnSim matrix of simulated coefficients
 	# NB! When the model has very little error, we get CV's very close to zero...
 	# we need to multiply everything out to avoid hitting zero Non-Hermitian errors
 	fact = 1#10^6#(round(log(abs(minimum(CVmat)))))
 	println("multiplier: ", fact)
 	th = fact*th
 	se_model = se_model*fact
 	CVmat = CVmat*fact

   d = length(th);
   sim_th = nanmat(d, nSim);
   sim_se_model = nanmat(1, nSim);
   # for iSim = 1:nSim
   hits = 0
   tries = 0
   while hits < nSim && tries < 2.5*nSim
   		tries = tries + 1
		X = rand(Chisq(n-d))
		sim_se_model[hits+1] = se_model .* sqrt((n-d)/X);
		try
			mvnorm = MvNormal(th, sim_se_model[hits+1]*CVmat)
			sim_th[1:d, hits+1] = rand(mvnorm)/fact
			hits = hits + 1
		catch e
			tries = tries + 1
		end
	end
	if hits < nSim
		println(join(["		Expected ", nSim, " sim models but only found ", hits]))
	end
	ix5 = round(Int,0.025 * hits)
   	ix95 = round(Int,0.975 * hits)
   	CI = nanmat(d, 2);
	p = nanmat(d, 1);
   	if ix5 != 0 && ix95 !=0
		for ith = 1:d
			ths = sim_th[ith,:]
	   		nonnanth = ths[findall(x->!isnan(x), ths)]
			ths = sort(nonnanth)
			lb = ths[ix5]
			ub = ths[ix95]
			CI[ith,1] = lb
			CI[ith,2] = ub
			p[ith,1] = sum(x->x<=0, ths)/length(nonnanth)
	   	end
   	else
   		println(join(["		simulated fit did not converge! Can't simulate coeff p values."]))
	end
	return (CI, vec(p), hits)
end

function plot_OLSC_model(X,y,yFit;xIdx=1)
	# 
	# 	We have lots of predictors, so we want to pick an informative one to plot.
	# 
	plot(X[xIdx,:],y, "k.", label="data", alpha=0.2)
	plot(X[xIdx,:],yFit, "r.", label="fit")
	gca().legend()
end
function plot_OLSC_model(df::DataFrame;xIdx=1)
	# 
	# 	We have lots of predictors, so we want to pick an informative one to plot.
	# 
	X = df.X[1]
	y = df.y[1]
	yFit = df.yFit[1]
	plot_OLSC_model(X,y,yFit,xIdx=xIdx)
	gca().set_title(join(["Rsq=", round(df.Rsq[1], digits=3)]))
	gca().set_xlabel(df.predictors[1][xIdx])
	gca().set_ylabel(df.predicted[1])
end


function getCompositeTheta(ths, se_ths, dofs)
# 	% 
# 	% 	Translated from STAT_Collate photometry Matlab
# 			Difference: only one theta (e.g., b1) calculated at a time rather than all at once... 
#				will use a wrapper to do for each other theta (e.g., b0, b2...)
# 	% 	
	# Only use non-nan indicies
	idxs_nonan = findall(x->!isnan(x), dofs)

	N = length(ths[idxs_nonan]); # the number of datasets, N degrees of freedom
	# NN = N.*ones(1, size(ths, 2)); # I think the number of ths. We only do one th at a time, so NN is ignored
	# NN = [number of b0, number of b1, number of b2...] this was needed because not all sets had tdt. We just need N
	
	
	meanTh = 1/N .* nansum(ths[idxs_nonan]);
	propagated_se_th = 1/N .* sqrt(nansum(se_ths[idxs_nonan].^2))
	mdf = sum(dofs[idxs_nonan])#.*ones(1, size(meanTh,2)); # m is the number of degrees of freedom across all the models
	# % 
	# % 	Now, calculate the CI = b +/- t(0.025, n(m-1))*se
	# % 
	# for nn = 1:length(meanTh, 2)
		# CImin(nn) = meanTh(nn) - abs(tinv(.025,numel(NN(nn))*(mdf(nn) - 1))).*propagated_se_th(nn);
		# CImax(nn) = meanTh(nn) + abs(tinv(.025,numel(NN(nn))*(mdf(nn) - 1))).*propagated_se_th(nn);
	if isnan(mdf)
		CImin = NaN
		CImax = NaN
	else
		try
			CImin = meanTh - abs( quantile(TDist(N*(mdf - 1)),0.025) ).*propagated_se_th;
			CImax = meanTh + abs( quantile(TDist(N*(mdf - 1)),0.025) ).*propagated_se_th;
		catch
			warning(join(["Got error at TDist... N*(mdf - 1))=", N*(mdf - 1), "\n Returning NaN for CI..."]))
			CImin = NaN
			CImax = NaN
		end
	end
	# CImin = meanTh - abs( tinv(.025,N*(mdf - 1)) ).*propagated_se_th(nn);
	# CImax = meanTh + abs( tinv(.025,N*(mdf - 1)) ).*propagated_se_th(nn);
	
# %                 Tried below, too, but yields same result. Not different
# %                 and I think above is correct version
# %                 CIminA(nn) = meanTh(nn) - abs(tinv(.025,numel(mdf(nn))*(NN(nn) - 1))).*propagated_se_th(nn);
# % 				CImaxA(nn) = meanTh(nn) + abs(tinv(.025,numel(mdf(nn))*(NN(nn) - 1))).*propagated_se_th(nn);
	# end


	return (meanTh, propagated_se_th, CImin, CImax, mdf)
end



function plot_residuals_by_trial(model_df::DataFrame, train_df::DataFrame; trialIDs=[], predictor="predictor", predicted="predicted")
#     We will extract the residuals of the TRAINING model and plot them by trial...
    if isempty(trialIDs)
        trialIDs = unique(train_df.TrialNo)
        println(trialIDs)
    end
    # take only the training data to look at...
    y = vec(model_df.y[1])
    yfit = model_df.prediction[1]
    TrialNo = train_df.TrialNo
    X=train_df.X
    plot_residuals_by_trial(y, yfit, data_df=train_df, TrialNo=TrialNo, trialIDs=trialIDs, predictor="predictor", predicted="predicted")
end    

function plot_residuals_by_trial(oscf_df::DataFrame;trialNo=[], data_df::DataFrame=DataFrame(), trialIDs=[], predictor="predictor", predicted="predicted")
#     We will extract the residuals of the TRAINING model and plot them by trial...

    if isempty(trialIDs)
        trialIDs = unique(data_df.TrialNo)
    end
    TrialNo = data_df.TrialNo
    # take only the training data to look at...
    y = oscf_df.y[1]
    yfit = oscf_df.yFit[1]
    plot_residuals_by_trial(vec(y), yfit, data_df=data_df, TrialNo=TrialNo, trialIDs=trialIDs, predictor=predictor, predicted=predicted)
end    
    
function plot_residuals_by_trial(y::AbstractArray{Float64,1}, yfit; data_df::DataFrame=DataFrame(), TrialNo=[], trialIDs=[], predictor="predictor", predicted="predicted")
#     We will extract the residuals of the TRAINING model and plot them by trial...
#     if isempty(trialIDs)
#         trialIDs = unique(train_df.TrialNo)
#         println(trialIDs)
#     end

    
#     # take only the training data to look at...
#     y = model_df.y[1]
#     yfit = model_df.prediction[1]
    for i = 1:length(trialIDs)
        idxs = findall(x->x==i, TrialNo)
        resid = yfit[idxs].-y[idxs]
        figure(figsize=(8,2))
        ax1 = subplot(1,3,1)
        ax1.plot(data_df.X[idxs], resid)
        ax1.set_title(join(["t#", i, " ", predicted]))
        ax1.set_xlabel("time")
        ax1.set_ylabel(join(["residuals (error) \n actual_",predicted," - fit"]))
        if !isempty(resid) && length(resid)>1
            ax1.set_ylim([minimum(resid), maximum(resid[1:end-1])])
        end
        
        ax2 = subplot(1,3,2)
        ax2.plot(data_df.X[idxs], data_df.Y[idxs])
        ax2.set_xlabel("time")
        ax2.set_ylabel(predictor)
        ax2.set_title(join([predictor, " vs time"]))
        ax3=subplot(1,3,3)
        ax3.plot(data_df.X[idxs], y[idxs], label="actual")
        ax3.plot(data_df.X[idxs], yfit[idxs], label="fit")
        ax3.set_ylabel(join([predicted]))
        ax3.set_xlabel("time")
        ax3.legend()
        if !isempty(resid) && length(resid)>1
            ax3.set_ylim([minimum(yfit[idxs]), maximum(yfit[idxs[1:end-1]])])
        end
        
#         y1 = vcat(minimum(yfit[1:end-2]), minimum(y[1:end-2]))
#         y2 = vcat(maximum(yfit[1:end-2]), maximum(y[1:end-2]))
#         ax.set_ylim([minimum(y1), maximum(y2)])
    end
end
 

function select_single_timepoint_per_trial_dataset(df::DataFrame)
    trialIDs = unique(df.TrialNo) 
    one_idx_per_trial_idxs = []
    for i = trialIDs
        ixs = findall(x->x==i, df.TrialNo)
        ix = rand(ixs)
        push!(one_idx_per_trial_idxs, ix)
    end
    newData = similar(df,0)
    # add the new data row by row
    for i = 1:length(one_idx_per_trial_idxs)
        push!(newData, df[one_idx_per_trial_idxs[i],:])
    end
    reducedPool = df[setdiff(1:end, one_idx_per_trial_idxs), :]
    return (newData, reducedPool)
end



function get_model_stats(model; n::Int=0, th::Vector{Float64}=Vector{Float64}(undef, 0), se_th::Vector{Float64}=Vector{Float64}(undef, 0), dof::Int=0, th_names::Vector{String}=[""])
	# if using Julia GLM.jl, input model. Otherwise, input the params from my ridge model
    if isempty(th)
        th = StatsBase.coef(model)
        se_th = stderror(model)
        dof = dof_residual(model)
        n = floor(Int,dof_residual(model)+length(th))
    end
    if isempty(th_names[1])
        th_names = Vector{String}(undef, 0)
        for i=0:length(th)-1
            push!(th_names, join(["th", i]))
        end
    end
    std_th = se_th .* sqrt(n)
    p_z = Vector{Float64}(undef, 0)
    CImin_z = Vector{Float64}(undef, 0)
    CImax_z = Vector{Float64}(undef, 0)
    p_t = Vector{Float64}(undef, 0)
    CImin_t = Vector{Float64}(undef, 0)
    CImax_t = Vector{Float64}(undef, 0)
    for i = 1:length(th)
        z = OneSampleZTest(th[i], std_th[i], n, 0)
        push!(CImin_z, confint(z, level=0.95, tail=:both)[1])
        push!(CImax_z, confint(z, level=0.95, tail=:both)[2])
        push!(p_z, pvalue(z; tail = :both))
        t = OneSampleTTest(th[i], std_th[i], n, 0)
        # if round(confint(t; level = 0.95, tail = :both)[1],digits=2) != round(th[i] - abs(quantile(TDist(dof-1),0.025)).*se_th[i], digits=2)
        #     warning(join(["Expected CImin=", confint(t; level = 0.95, tail = :both)[1], " but got ", th[i] - abs(quantile(TDist(dof-1),0.025)).*se_th[i]]))
        # end
        push!(CImin_t, th[i] - abs(quantile(TDist(dof-1),0.025)).*se_th[i])
        push!(CImax_t, th[i] + abs(quantile(TDist(dof-1),0.025)).*se_th[i])
        push!(p_t, pvalue(t; tail = :both))
    end
    return DataFrame(th_names=th_names, n=n.*ones(length(th)), th=th, se_th=se_th, std_th=std_th, dof=dof*ones(length(th)), 
        p_t=p_t, CImin_t=CImin_t, CImax_t=CImax_t, p_z=p_z, CImin_z=CImin_z, CImax_z=CImax_z)
end


function nx_kcomparisons_GLM(model_fxn; custom_GLM=false, loss_fxn=nothing, Formula, k_sets, modelName::String)
	if custom_GLM == false
		println("using GLM.jl, not custom methods")
	else
		println("using custom methods, no formula needed")
	end
    if isnothing(loss_fxn)
        println("Loss function not specified, using MSE Loss (linear model)")
        loss_fxn = MSELoss
    end
    # 
    #  Data collected as vectors for each kth training set. Pushed within are the results of xval
    #     ths = [ [[kths1-1], [kths1-2]...], [[kths2-1], [kths2-2]...], ...]
    #
    #
    train_ths = [[] for _=1:length(k_sets)]
    train_se_ths = [[] for _=1:length(k_sets)]
    train_dof =  [[] for _=1:length(k_sets)]
    train_CImin =  [[] for _=1:length(k_sets)]
    train_CImax =  [[] for _=1:length(k_sets)]
    train_rsqs = [[] for _=1:length(k_sets)]
    train_loss = [[] for _=1:length(k_sets)]
    train_aic = []
    train_aicc = []
    train_bic = []
    
    test_rsqs = [[] for _=1:length(k_sets)]
    test_loss = [[] for _=1:length(k_sets)]
    
    train_models_df = []
    

    kidxs = 1:length(k_sets)
    for i = kidxs
        println("----- Working on k-set no. ", i) 
        idxs = kidxs[1:end .!= i]
        for j = 1:length(idxs)
        	if !custom_GLM
	            trainmodel, testmodel, df = model_fxn(Formula, k_sets[i], k_sets[idxs[j]], modelName=modelName, verbose=false);
            else
            	trainmodel, testmodel, df = model_fxn(k_sets[i], k_sets[idxs[j]], modelName=modelName, verbose=false);
        	end
            # if is loop 1, get the training model stats:
            if j==1
                y_actual_Sn = df.y[1]
                prediction_Sn = df.prediction[1]
                (th_Sn, se_th_Sn, dof_Sn, CImin_Sn, CImax_Sn) = get_model_stats(trainmodel)
                loss_Sn = loss_fxn(y_actual_Sn,prediction_Sn)
                Rsq_Sn = df.Rsq[1]
                push!(train_ths[i], th_Sn)
                push!(train_se_ths[i], se_th_Sn)
                push!(train_dof[i], dof_Sn)
                push!(train_CImin[i], CImin_Sn)
                push!(train_CImax[i], CImax_Sn)
                push!(train_rsqs[i], Rsq_Sn)
                push!(train_loss[i], loss_Sn)
                push!(train_models_df, df)
                push!(train_aic, aic(trainmodel))
                push!(train_aicc, aicc(trainmodel))
                push!(train_bic, bic(trainmodel))
                # print summary of the training model
                println(trainmodel)
            end
            y_actual_test = df.y[2]
            prediction_test = df.prediction[2] 
            loss_test = loss_fxn(y_actual_test,prediction_test)
            Rsq_test = df.Rsq[2]
            push!(test_rsqs[i], loss_test)
            push!(test_loss[i], Rsq_test)
        end
    end
    # collect summary stats:
    mean_train_rsq_k = [mean(x) for x in train_rsqs]
    mean_test_rsq_k = [mean(x) for x in test_rsqs]
    mean_train_loss_k = [mean(x) for x in train_loss]
    mean_test_loss_k = [mean(x) for x in test_loss]
    results = DataFrame(train_ths = train_ths,
            train_se_ths = train_se_ths,
            train_dof = train_dof,
            train_CImin = train_CImin,
            train_CImax = train_CImax,
            train_rsqs = train_rsqs,
            train_loss = train_loss,
            test_rsqs = test_rsqs,
            test_loss = test_loss,
            mean_train_rsq_k = mean_train_rsq_k,
            mean_test_rsq_k = mean_test_rsq_k,
            mean_train_loss_k = mean_train_loss_k,
            mean_test_loss_k = mean_test_loss_k,
            train_models_df = train_models_df,
            train_aic = train_aic,
            train_aicc = train_aicc,
            train_bic = train_bic,
            )
    
    figure(figsize=(14,3))
    ax1 = subplot(1,4,1)
    plot_k_xval_stat(train_rsqs; statname="train Rsq")
    ax2 = subplot(1,4,2)
    plot_k_xval_stat(test_rsqs; statname="test Rsq")
    ax3 = subplot(1,4,3)
    plot_k_xval_stat(train_loss; statname="train Loss")
    ax4 = subplot(1,4,4)
    plot_k_xval_stat(test_loss; statname="test Loss")
    set_yaxes_same_scale([ax1, ax2])
    set_yaxes_same_scale([ax3, ax4])
    
    figure(figsize=(14,3))
    ax1 = subplot(1,3,1)
    plot_k_xval_stat(wrap(train_aic); statname="train AIC")
    ax2 = subplot(1,3,2)
    plot_k_xval_stat(wrap(train_aicc); statname="test AICc")
    ax3 = subplot(1,3,3)
    plot_k_xval_stat(wrap(train_bic); statname="train BIC")

    
    return results
end
function plot_k_xval_stat(statvector; statname="stat")
    mean_stat_k = [mean(x) for x in statvector]
    overallmean = mean(mean_stat_k)
    k = 1:length(mean_stat_k)
    jittervec = (rand(statvector[1]) .- 0.5) ./ 3
    plot([0,length(k)], [0, 0], "k-")
    for i in k
        plot(i, mean_stat_k[i], "r.", markersize=20)
        ee = i.*ones(length(statvector[i])) .+ jittervec
        plot(ee, statvector[i], "k.")
    end
    title(join([statname, " mean=", round(overallmean, digits=2)]))
    gca().set_xlabel("k-set")
    gca().set_ylabel(statname)
end
function wrap(vect)
    [[x] for x in vect]
end
function collect_theta(thetas, th_idx)
    # collects and wraps for plotting
   [[x[1][th_idx]] for x in thetas] 
   # results_all.train_ths[1][1][2]
end


function unwrap_theta(thetas)
   # unwraps the nested []
    [x[1] for x in thetas];
end

function theta_summary(stats_df::DataFrame; Mode = "sparseFit", result_df::DataFrame)
	# 
	# Modes = "spaseFit", "oneFit"
	# 	Use the stats_df for my own ridge_regression
	# 	use kfold_result for the sparseFit case
	#
	if Mode=="sparseFit"
		# for propapagted error across many fits
	    # get the dof for each model fit
	    dofs = [x[1] for x in result_df.train_dof];
	    # get all the thetas organized
	    # println(result_df.train_ths[1][1])
	    d = length(result_df.train_ths[1][1])
	    th_by_k_dataset = result_df.train_ths
	    se_by_k_dataset = result_df.train_se_ths
	    composite_th = []
	    composite_se = []
	    composite_CImin = []
	    composite_CImax = []
	    composite_mdof = []
	    for i = 1:d
	        # extract the proper theta and se_th
	        th_plot = collect_theta(th_by_k_dataset, i)
	        th = unwrap_theta(th_plot)
	        se_th_plot = collect_theta(se_by_k_dataset, i);
	        se_th = unwrap_theta(se_th_plot)
	        # get composite theta
	        if isnan(dofs[1])
	        	meanTh = NaN#vec(nanmat(d,1))
	        	propagated_se_th = NaN #vec(nanmat(d,1))
	        	CImin = NaN#vec(nanmat(d,1))
	        	CImax = NaN#vec(nanmat(d,1))
	        	mdf = NaN
	        else
		        (meanTh, propagated_se_th, CImin, CImax,mdf) = getCompositeTheta(th, se_th, dofs)
	        end
	        push!(composite_th, meanTh)
	        push!(composite_se, propagated_se_th)
	        push!(composite_CImin, CImin)
	        push!(composite_CImax, CImax)
	        push!(composite_mdof, mdf)
	    end
	    k = length(dofs)
    elseif Mode == "oneFit"
	    dofs = [x[1] for x in stats_df.dof];
	    # get all the thetas organized
	    d = length(stats_df.th)
	    composite_th = stats_df.th
	    composite_se = stats_df.se_th
	    composite_CImin = stats_df.CImin_t
		composite_CImax = stats_df.CImax_t
		k = result_df.n_kfold_xval[1]
    end
    f = figure(figsize=(5,3))
    ax = subplot(1,1,1)
    plot_with_CI(composite_th, composite_CImin, composite_CImax, ax=ax)
    ax.set_title(join(["coefficients, k=", k]))
    ax.set_xticks(collect(1:d))
    # xticks(collect(1:d), labels=stats_df.)
    return (composite_th, composite_se, composite_CImin, composite_CImax, ax, f, composite_mdof)
end



function kfold_xval(model_fxn=run_regression; yID::Symbol, predictornames=[:A], loss_fxn=nothing, k_sets, modelName::String, verbose=false)
    println("========== k-folds xvalidation ", timestamp_now(), "==========")
    if isnothing(loss_fxn)
        println("	Loss function not specified, using MSE Loss (linear model)")
        loss_fxn = MSELoss
    end
    # 
    #  Data collected as vectors for each kth training set. Pushed within are the results of xval
    #     ths = [ [[kths1-1], [kths1-2]...], [[kths2-1], [kths2-2]...], ...]
    #
    #
    # train_ths = [[] for _=1:length(k_sets)]
    # train_se_ths = [[] for _=1:length(k_sets)]
    # train_dof =  [[] for _=1:length(k_sets)]
    # train_CImin =  [[] for _=1:length(k_sets)]
    # train_CImax =  [[] for _=1:length(k_sets)]
    # train_rsqs = [[] for _=1:length(k_sets)]
    # train_loss = [[] for _=1:length(k_sets)]
    # train_aic = []
    # train_aicc = []
    # train_bic = []
    
    # test_rsqs = [[] for _=1:length(k_sets)]
    # test_loss = [[] for _=1:length(k_sets)]
    
    # train_models_df = []

    allset = k_sets[1]
    idxset = [collect(1:nrow(k_sets[1]))]
    s = nrow(k_sets[1])
    for i = 2:length(k_sets)
    	allset = vcat(allset, k_sets[i])
    	xx=collect(s+1:s+nrow(k_sets[i]))
    	push!(idxset,xx)
    	s = s+nrow(k_sets[i])
	end
    # println(idxset[2])
    #
    #	We will test a range of lambdas
    #
    lambdas = [0.0, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4] # this one will be our running pool
    all_lambdas = lambdas # this one will be where we keep track of everything we have tried
    # println(" lambdas=", lambdas)
    # println(" ")

    train_rsqs = []
    test_rsqs = []
    train_loss = []
    test_loss = []
    mean_test_loss_ea_lam = Vector{Float64}(undef, 0)
    mean_train_loss_ea_lam = Vector{Float64}(undef, 0)
    mean_train_Rsq_ea_lam = Vector{Float64}(undef, 0)
    mean_test_Rsq_ea_lam = Vector{Float64}(undef, 0)

    kidxs = 1:length(k_sets)
    tries = 1;
    repeat = true
    lam_star = 0.0
    while repeat == true
    	print("Testing lam = [")
    	pretty_print_list(lambdas, orient="horizontal")
    	println("] ", timestamp_now())
	    for lam in lambdas
	    	if verbose
				println("-- lambda=", lam, " (", timestamp_now(), "):")
			end
			push!(train_rsqs, Vector{Float64}(undef, 0))
			push!(test_rsqs, Vector{Float64}(undef, 0))
			push!(train_loss, Vector{Float64}(undef, 0))
			push!(test_loss, Vector{Float64}(undef, 0))
		    for i = kidxs
		        # println("	Working on k-set no. ", i) 
		        # We will use ALL the data except the k-set to fit

		        idxs = kidxs[1:end .!= i]
		        fitset = copy(allset)
		        fitset = deleterows!(fitset, idxset[i])
		        testset = k_sets[i]
		        if verbose
		        	println("the allset has size ", size(allset))
		        	println("the fitset has size ", size(fitset))
					println("the testset has size ", size(testset))
				end

		        
		        (result_df, model_summary, stats_df) = model_fxn(fitset, predictornames, yID; lam=lam, add_th0=true)
		        y_actual_Sn = result_df.y[1]
		        y_fit_Sn = result_df.yFit[1]
		        th_Sn = stats_df.th[1]
		        # println(size(y_actual_Sn))
		        # println(size(y_fit_Sn))
		        loss_Sn = loss_fxn(y_actual_Sn,y_fit_Sn)
                Rsq_Sn = result_df.Rsq[1]
                push!(train_rsqs[end], Rsq_Sn)
                push!(train_loss[end], loss_Sn)
                if verbose
	                println(model_summary)
                end
                
                # Now get the loss and Rsq on the test set
                # it's convenient to grab our X and y from the fit function...
                (result_df_test, _, _) = model_fxn(testset, predictornames, yID; lam=lam, add_th0=true)
                y_actual_test = result_df_test.y[1]
                X_test = result_df_test.X[1]
                y_predicted_test = th_Sn'*X_test #'
                loss_test = loss_fxn(y_actual_test,y_predicted_test)
                Rsq_test = Rsq(y_actual_test,y_predicted_test)
	            push!(test_rsqs[end], loss_test)
	            push!(test_loss[end], Rsq_test)
	            # println(" 		Loss_Sn=", loss_Sn, " | Rsq_Sn=", Rsq_Sn)
	            # println(" 		:::Loss_test=", loss_test, " | Rsq_test=", Rsq_test)		        
		    end
		    
		    mean_train_Rsq = mean(train_rsqs[end])
		    push!(mean_train_Rsq_ea_lam, mean_train_Rsq)
		    mean_test_Rsq = mean(test_rsqs[end])
		    push!(mean_test_Rsq_ea_lam, mean_test_Rsq)
		    mean_train_loss = mean(train_loss[end])
			push!(mean_train_loss_ea_lam, mean_train_loss)
		    mean_test_loss = mean(test_loss[end])
            push!(mean_test_loss_ea_lam, mean_test_loss)
            if verbose
	            print("	Loss_Sn: [Mean=", round(mean_train_loss, digits=3),"] - ")
			    pretty_print_list(train_loss[end], orient="horizontal", digits=3)
			    println(" ")
			    print(" 	Loss_test: [Mean=", round(mean_test_loss, digits=3),"] - ")
			    pretty_print_list(test_loss[end], orient="horizontal", digits=3)
			    println(" ")
			    print("	Rsq_Sn: [Mean=", round(mean_train_Rsq, digits=3),"] - ")
			    pretty_print_list(train_rsqs[end], orient="horizontal", digits=3)
			    println(" ")
			    print("	Rsq_test: [Mean=", round(mean_test_Rsq, digits=3),"] - ")
			    pretty_print_list(test_rsqs[end], orient="horizontal", digits=3)
			    println(" ")
		    end
		    
	    end
	    # find the best lambda and expand in that region if the difference in loss between lambdas is > 1%
	    # println("mean_test_loss_ea_lam=", mean_test_loss_ea_lam)
	    lambest = findall(x-> x==minimum(mean_test_loss_ea_lam), mean_test_loss_ea_lam)[1]
	    # println("lambest=", lambest)
	    if tries == 1
	    	if all_lambdas[lambest] == 0.0
		    	# we should explore smaller lams
		    	if (mean_test_loss_ea_lam[2] - mean_test_loss_ea_lam[1])/(mean_test_loss_ea_lam[2] + mean_test_loss_ea_lam[1]) > 0.01
		    		repeat = true
		    		lambdas = [0,10^-10,10^-9,10^-8,10^-7,10^-6,10^-5]
		    		all_lambdas = vcat(all_lambdas, lambdas)
	    		else
	    			println(" ")
					println("Best lam was 0, and difference in loss with next lam is < 1%. Breaking")
	    			repeat = false
	    		end

			elseif all_lambdas[lambest]==lambdas[end]
				# we should explore bigger lams
				if (mean_test_loss_ea_lam[end] - mean_test_loss_ea_lam[end-1])/(mean_test_loss_ea_lam[end] + mean_test_loss_ea_lam[end-1]) > 0.01
		    		repeat = true
		    		lambdas = [0,10^5,10^6,10^7,10^8,10^9,10^10]
		    		all_lambdas = vcat(all_lambdas, lambdas)
	    		else
	    			println(" ")
					println("Best lam was ", all_lambdas[lambest], " and difference in loss with next lower lam is < 1%. Breaking")
	    			repeat = false
	    		end
    		else
    			nextbest = findall(x->x==minimum([mean_test_loss_ea_lam[lambest-1], mean_test_loss_ea_lam[lambest+1]]), [mean_test_loss_ea_lam[lambest-1],mean_test_loss_ea_lam[lambest+1]])[1]
    			# println("nextbest", nextbest)
    			if nextbest == 1
	    			lambdas = collect(range(all_lambdas[lambest-1], stop=all_lambdas[lambest], length=10))
    			elseif nextbest == 2
    				lambdas = collect(range(all_lambdas[lambest], stop=all_lambdas[lambest+1], length=10))
				end
				all_lambdas = vcat(all_lambdas, lambdas)
				if verbose
					println("~~Testing fine-grain lambdas: ", lambdas)
				end
			end

		elseif tries == 2
			if all_lambdas[lambest] == 0.0 || all_lambdas[lambest]==lambdas[end]
				println(" ")
				println("Best lam was ", lambdas[lambest], " after 2 tries. Breaking")
    			repeat = false
			else
				nextbest = findall(x->x==minimum([mean_test_loss_ea_lam[lambest-1], mean_test_loss_ea_lam[lambest+1]]), [mean_test_loss_ea_lam[lambest-1],mean_test_loss_ea_lam[lambest+1]])[1]
				# println("nextbest", nextbest)
				if nextbest == 1
	    			lambdas = collect(range(all_lambdas[lambest-1], stop=all_lambdas[lambest], length=10))
				elseif nextbest == 2
					lambdas = collect(range(all_lambdas[lambest], stop=all_lambdas[lambest+1], length=10))
				end
				all_lambdas = vcat(all_lambdas, lambdas)
				if verbose
					println("~~Testing fine-grain lambdas: ", lambdas)
				end
			end
		else
			println(" ")
			println("Best lam was ", all_lambdas[lambest], " after 3 tries. Breaking")
			repeat = false
		end
		tries = tries + 1
		lam_star = all_lambdas[lambest]
	end
	
	
	println("Fitting with best lambda and getting final model...")
	println(" 	Lam=", lam_star)
    (result_df, model_summary, stats_df) = model_fxn(allset, predictornames, yID; lam=lam_star, add_th0=true)
    y_actual_final = result_df.y[1]
    y_fit_final = result_df.yFit[1]
	loss = loss_fxn(y_actual_final,y_fit_final)
	Rsq_final = result_df.Rsq[1]
    println("	Loss=", round(loss, digits=3))
    println("	Rsq=", round(Rsq_final, digits=3))
    result_df[:n_kfold_xval][1] = length(ksets)

    xval_data = DataFrame(k=length(ksets), all_lambdas=[all_lambdas], lam_star=lam_star, 
    	mean_train_loss_ea_lam=[mean_train_loss_ea_lam], mean_test_loss_ea_lam=[mean_test_loss_ea_lam], 
    	mean_train_Rsq_ea_lam=[mean_train_Rsq_ea_lam], mean_test_Rsq_ea_lam=[mean_test_Rsq_ea_lam])
    # plot_xval_loss(xval_data)
    # plot_xval_fit_results(result_df,model_summary; modelClass="linear")
    return (result_df, model_summary, stats_df, xval_data)
end

function no_xval_control(model_fxn; 
	predictornames, yID, k_sets, modelName="", lam=0.0, add_th0=true)
	
	allset = k_sets[1]
    idxset = [collect(1:nrow(k_sets[1]))]
    s = nrow(k_sets[1])
    for i = 2:length(k_sets)
    	allset = vcat(allset, k_sets[i])
    	xx=collect(s+1:s+nrow(k_sets[i]))
    	push!(idxset,xx)
    	s = s+nrow(k_sets[i])
	end

	loss_fxn=MSELoss
	
	(result_df, model_summary, stats_df) = model_fxn(allset, predictornames, yID; lam=lam, add_th0=true)
    y_actual_final = result_df.y[1]
    y_fit_final = result_df.yFit[1]
	loss = loss_fxn(y_actual_final,y_fit_final)
	Rsq_final = result_df.Rsq[1]
    println("	Loss=", round(loss, digits=3))
    println("	Rsq=", round(Rsq_final, digits=3))
    result_df[:n_kfold_xval][1] = 0

    return (result_df, model_summary, stats_df, allset)
end

function plot_xval_loss(xval_data::DataFrame)
	lam_star = xval_data.lam_star[1]
	all_lambdas = xval_data.all_lambdas[1]
	mean_train_loss_ea_lam = xval_data.mean_train_loss_ea_lam[1]
	mean_test_loss_ea_lam = xval_data.mean_test_loss_ea_lam[1]
	mean_train_Rsq_ea_lam = xval_data.mean_train_Rsq_ea_lam[1]
	mean_test_Rsq_ea_lam = xval_data.mean_test_Rsq_ea_lam[1]
	# Plot the results
	#
	lamix = findall(x->x==lam_star, all_lambdas)[1]
	ii = sortperm(all_lambdas)
	# Loss:
	figure(figsize=(8,3))
	ax1=subplot(1,2,1)
	ax1.plot(all_lambdas[ii], mean_train_loss_ea_lam[ii], "k-", label="train")
	# ax1.set_title("training loss")
	ax1.set_xlabel("loss")
	ax1.set_ylabel("lambda")

	# ax2=subplot(1,2,1)
	ax1.plot(all_lambdas[ii], mean_test_loss_ea_lam[ii], "r-", label="test")
	ax1.set_title(join(["test loss*=",round(mean_test_loss_ea_lam[lamix], digits=3)]))
	ax1.set_ylabel("lambda")
	plotsignif(lam_star,maximum([mean_train_loss_ea_lam[lamix], mean_test_loss_ea_lam[lamix]])-0.025; ax=ax1)

	ax3=subplot(1,2,2)
	ax3.plot(all_lambdas[ii], mean_train_Rsq_ea_lam[ii], "k-", label="train")
	ax3.set_title(join(["test Rsq*=",round(mean_test_Rsq_ea_lam[lamix], digits=3)]))
	# ax3.set_title("train Rsq")
	# ax3.set_xlabel("lambda")
	ax3.set_ylabel("Rsq")
	# ax3=subplot(1,4,3)
	ax3.plot(all_lambdas[ii], mean_test_Rsq_ea_lam[ii], "r-", label="test")
	# ax3.set_title("test Rsq")
	ax3.set_xlabel("lambda")
	plotsignif(lam_star,maximum([mean_train_Rsq_ea_lam[lamix], mean_test_Rsq_ea_lam[lamix]])-0.025; ax=ax3)
	ax3.legend()
end


function get_k_1pt_datasets(full_df::DataFrame; threshold_retention=0.8)
    #
    # note that threshold retention can be a fract or an int for the 
    # min number of trials to include. We actually probably do better to 
    # some extent after eliminating the really short ones...
    #
    reducedPool = full_df
    one_per_trial_dfs = []
    trialsAvail_init = unique(reducedPool.TrialNo)
    trialsAvail = trialsAvail_init
    println(length(trialsAvail), " trials avail initially.")
    println("     ")
    if threshold_retention < 1
        threshold_retention = threshold_retention*length(trialsAvail_init)
    end
    while length(trialsAvail) > threshold_retention
        print(length(trialsAvail), " ")
        (one_per_trial_df, reducedPool) = select_single_timepoint_per_trial_dataset(reducedPool)
        trialsAvail = unique(reducedPool.TrialNo)
        push!(one_per_trial_dfs, one_per_trial_df)
    end
    println(" ")
    println(" ")
    println("Found ", length(one_per_trial_dfs), " one-point-per-trial datasets with at least ", round(threshold_retention, digits=2), " trials")
    return one_per_trial_dfs
end


function modelSelectionByAICBICxval(all_df::DataFrame, yID::Symbol, formulas, modelNames, modelClass="logit"; n_iters=10,updownsampleYID=false, figurePath=".", savePath=".", suppressFigures=suppressFigures, slice="")
	maxattempts = 10
    #
    # Here, we do model selection by AIC, BIC criteria and also find fit coeffs for bootstrapped models.
    #  this seems useful for logit models where we only fit on a subset of the data at a given time...
    #
    AICs = [[] for _=1:length(modelNames)]
    AICcs = [[] for _=1:length(modelNames)]
    BICs = [[] for _=1:length(modelNames)]
    testloss = [[] for _=1:length(modelNames)]
    Sn_accuracy = [[] for _=1:length(modelNames)]
    test_accuracy = [[] for _=1:length(modelNames)]
    th_names = [[] for _=1:length(modelNames)]
    ths = [[] for _=1:length(modelNames)]
    se_ths = [[] for _=1:length(modelNames)]
    dofs = [[] for _=1:length(modelNames)]


    #
    #	Check the df isn't empty
    #
    if !isempty(all_df) && sum(all_df[yID])>15
    	try 
	    	npercat = sum(all_df[yID])#length(findall(x->x, all_df[yID]))
    	catch
    		println("\x1b[31m\"!!!!!!!!!!!!!!!! ERROR!\"\x1b[0m")
    		println("\x1b[31m\"      sum(all_df[yID])=",sum(all_df[yID]), "\"\x1b[0m")
    		println("\x1b[31m\"      yID=",yID, "\"\x1b[0m")
    		println("\x1b[31m\"      all_df[yID]=",all_df[yID], "\"\x1b[0m")
		end
    	println(join(["Downsampling the no-lick state => using ", npercat, "=n, the number of first-licks"]))
	 #    #
	 #    # Get the number of up/down samples by querying smote -- Jan thinks is no good
	 #    #
	 #    X2, y2 =smote(all_df[!,[:Y, :X]], all_df[yID], k = 5, pct_under = 150, pct_over = 200)
		# df_balanced = X2
		# df_balanced[yID] = y2;
		# a = countmap(df_balanced[yID])
		# if a[true] != a[false]
		# 	println(a)
		# 	error("we didn't get an even number.")
		# end
	 #    npercat = a[true]

	    # println(join(["Used smote to estimate the number of up/down sampling needed => using ", npercat, "=n"]))
	    retry = true # we should try to sample enough times to get cholesky factorizable matrix. I'll give it 10 shots
	    attempts = 1
	    nomodelfits = false
	    for i = 1:n_iters
	    	if !suppressFigures
		        progressbar(i,n_iters)
	        end
	        retry=true
	        while retry
		        try # we might have a bad subsample and not be able to fit the model
			        if updownsampleYID

			            working_df = updownsample(all_df, all_df[yID], npercat, downOnly=true); # I chose this based on the 
			            # sample n that was picked by smote. In the next version for across models, be sure to adjust for the dataset...
			        else
			            working_df = all_df
			        end
			        train, test = TrainTestSplit(working_df, 0.75)
			        for model = 1:length(modelNames)
			            if modelClass == "logit"
			                # get AIC, BIC with the original model and an average coefficient with propagated se
			                if length(AICs[model]) != n_iters
				                (logit_model, _, accuracy_Sn, accuracy_test) = build_and_report_logit_model(formulas[model], train, test; modelName=modelNames[model], modelClass="logit", verbose=false)
				                push!(AICs[model], aic(logit_model))
				                push!(AICcs[model], aicc(logit_model))
				                push!(BICs[model], bic(logit_model))
				                push!(Sn_accuracy[model], accuracy_Sn)
				                push!(test_accuracy[model], accuracy_test)
				                stats_df = get_model_stats(logit_model)
				                push!(th_names[model], stats_df.th_names)
				                push!(ths[model], stats_df.th)
				                push!(se_ths[model], stats_df.se_th)
				                push!(dofs[model], stats_df.dof)
			                else
			                	# println("*************went past number of iters desired, I think this happens when one model fits ok for a sampled dataset, when a later model does not.")# break
			                end
			            else
			                error("not implemented for non-logit yet")
			            end
			        end
			        retry = false
		        catch e 
		        	if isa(e, PosDefException)
			        	if attempts < maxattempts
			        		if !suppressFigures
			        			println("***PosDefException...retrying. (", attempts, "/", maxattempts, ")")
		        			end
			        		retry = true
			        		attempts = attempts + 1
			    		else
			    			println(join(["***PosDefException -- Matrix could not fit after ", attempts, " attempts at up/down sampling. Ignoring this set!"]))
			    			retry = false
			    			# rethrow()
			    			nomodelfits = true
			        	end
			        else
			        	println(typeof(e))
			        	println("	")
			        	warning("YIKES! WEIRD ERROR - could be fit never converged!")
			        	# rethrow()
		        	end
		        end
	        end
	    end
    else
    	nomodelfits = true
    	npercat = 0
    end
    th_summary = DataFrame(modelName=[], composite_th=[], composite_se=[], composite_CImin=[], composite_CImax=[],composite_mdof = [])
    meanAIC = []
    meanAICc = []
    meanBIC = []
    meanAccuracy_Sn = []
    meanAccuracy_test = []
    axs = []
    fs = []
    if nomodelfits
    	for model = 1:length(modelNames)
    		#
    		#  get d for each model
    		#
    		d = NaN
    		try
	    		d = length(formulas[model].rhs) + 1
    		catch
    			d = 2
			end
			NAMES = String[]
			for i=0:d-1
				push!(NAMES, join(["th",i]))
			end
			push!(th_names[model], NAMES)
            push!(ths[model], vec(nanmat(d,1)))
            push!(se_ths[model], vec(nanmat(d,1)))
            push!(dofs[model], vec(nanmat(d,1))) 
	        result_df = DataFrame(train_dof = [vec(nanmat(d,1))], train_ths = [[vec(nanmat(d,1))]], train_se_ths=[[vec(nanmat(d,1))]])
	        
	        for i=2:n_iters
	        	push!(th_names[model], NAMES)
	            push!(ths[model], vec(nanmat(d,1)))
	            push!(se_ths[model], vec(nanmat(d,1)))
	            push!(dofs[model], vec(nanmat(d,1)))  
	            
	            xxx = DataFrame(train_dof = [dofs[model][i]], 
	                    train_ths = [[ths[model][i]]], 
	                    train_se_ths=[[se_ths[model][i]]],
	                )
	            
	            append!(result_df, xxx)   
                   
	        end

	        (composite_th, composite_se, composite_CImin, composite_CImax, ax, f, composite_mdof) = theta_summary(result_df; Mode = "sparseFit", result_df=result_df)
	        title(modelNames[model])
	        push!(axs, ax)
	        push!(fs, f)
	#         printFigure(join([modelNames[model], "_theta_summary_nboot", n_iters, "_npercat", npercat]); fig=fs[model])
	        
	        append!(th_summary, 
	            DataFrame(
	                modelName = modelNames[model],
	                composite_th=composite_th, 
	                composite_se=composite_se, 
	                composite_CImin=composite_CImin, 
	                composite_CImax=composite_CImax,
	                composite_mdof = composite_mdof,
	                ))
	        push!(meanAIC,NaN)
	        push!(meanAICc,NaN)
	        push!(meanBIC,NaN)
	        push!(meanAccuracy_Sn,NaN)
	        push!(meanAccuracy_test,NaN)
	        AICs[model] = vec(nanmat(n_iters,1))
			AICcs[model] = vec(nanmat(n_iters,1))
			BICs[model] = vec(nanmat(n_iters,1))
	        Sn_accuracy[model] = vec(nanmat(n_iters,1))
	        test_accuracy[model] = vec(nanmat(n_iters,1))

	    end
    else
	    for model = 1:length(modelNames)
	        result_df = DataFrame(train_dof = [dofs[model][1]], train_ths = [[ths[model][1]]], train_se_ths=[[se_ths[model][1]]])
	        for i=2:n_iters
	            append!(result_df, 
	            DataFrame(train_dof = [dofs[model][i]], 
	                    train_ths = [[ths[model][i]]], 
	                    train_se_ths=[[se_ths[model][i]]],
	                ))      
	        end
	        

	        (composite_th, composite_se, composite_CImin, composite_CImax, ax, f, composite_mdof) = theta_summary(result_df; Mode = "sparseFit", result_df=result_df)
	        title(modelNames[model])
	        push!(axs, ax)
	        push!(fs, f)
	#         printFigure(join([modelNames[model], "_theta_summary_nboot", n_iters, "_npercat", npercat]); fig=fs[model])
	        
	        append!(th_summary, 
	            DataFrame(
	                modelName = modelNames[model],
	                composite_th=composite_th, 
	                composite_se=composite_se, 
	                composite_CImin=composite_CImin, 
	                composite_CImax=composite_CImax,
	                composite_mdof = composite_mdof,
	                ))
	        push!(meanAIC,mean(AICs[model]))
	        push!(meanAICc,mean(AICcs[model]))
	        push!(meanBIC,mean(BICs[model]))
	        push!(meanAccuracy_Sn,mean(Sn_accuracy[model]))
	        push!(meanAccuracy_test,mean(test_accuracy[model]))
	        
	    end
    end
    set_xaxes_same_scale(axs)
    set_yaxes_same_scale(axs)
    for i=1:length(fs)
        # println(i)
        printFigure(join([modelNames[i], "_theta_summary_nboot", n_iters, "_npercat", npercat, slice]); fig=fs[i],figurePath=figurePath)
    end
    if suppressFigures
    	for i=1:length(fs)
    		close()
		end
	end
    
    f = figure(figsize=(20,3))
    ax1=subplot(1,3,1)
    compare_AICBIC(meanAIC, AICs; yl="AIC", iters=n_iters, ax=ax1, minmax="min")
    ax2=subplot(1,3,2)
    compare_AICBIC(meanAICc, AICcs; yl="AICc", iters=n_iters, ax=ax2, minmax="min")
    ax3=subplot(1,3,3)
    compare_AICBIC(meanBIC, BICs; yl="BIC", iters=n_iters, ax=ax3, minmax="min")
    printFigure(join(["AICBIC_summary_nboot", n_iters, "_npercat", npercat, slice]); fig=f, figurePath=figurePath)
    if suppressFigures
    	close()
	end
    
    f = figure(figsize=(20,3))
    ax1=subplot(1,3,1)
    compare_AICBIC(meanAccuracy_Sn, Sn_accuracy; yl="Train Accuracy", iters=n_iters, ax=ax1, minmax="max")
    ax1.set_ylim([0., 1.])
    ax2=subplot(1,3,2)
    compare_AICBIC(meanAccuracy_test, test_accuracy; yl="Test Accuracy", iters=n_iters, ax=ax2, minmax="max")
    ax2.set_ylim([0., 1.])
    printFigure(join(["Accuracy_summary_nboot", n_iters, "_npercat", npercat, slice]); fig=f, figurePath=figurePath)
    if suppressFigures
    	close()
	end
        
    results = DataFrame(
        AICs = AICs,
        meanAIC=meanAIC,
        AICcs = AICcs,
        meanAICc=meanAICc,
        BICs = BICs,
        meanBIC=meanBIC,
        testloss = testloss,
        Sn_accuracy = Sn_accuracy,
        meanAccuracy_Sn = meanAccuracy_Sn,
        test_accuracy = test_accuracy,
        meanAccuracy_test = meanAccuracy_test,
        th_names = th_names,
        ths = ths,
        se_ths = se_ths,
        dofs = dofs,
        th_summary = th_summary,
        )
    wd = pwd()
    try
        cd(savePath)
        CSV.write(join(["AICs_", slice, ".csv"]),DataFrame(AICs = AICs))
        CSV.write(join(["meanAIC_", slice, ".csv"]),DataFrame(meanAIC = meanAIC))
        CSV.write(join(["AICcs_", slice, ".csv"]),DataFrame(AICcs = AICcs))
        CSV.write(join(["meanAICc_", slice, ".csv"]),DataFrame(meanAICc = meanAICc))
        CSV.write(join(["BICs_", slice, ".csv"]),DataFrame(BICs = BICs))
        CSV.write(join(["meanBIC_", slice, ".csv"]),DataFrame(meanBIC = meanBIC))
        CSV.write(join(["testloss_", slice, ".csv"]),DataFrame(testloss = testloss))
        CSV.write(join(["Sn_accuracy_", slice, ".csv"]),DataFrame(Sn_accuracy = Sn_accuracy))
        CSV.write(join(["meanAccuracy_Sn_", slice, ".csv"]),DataFrame(meanAccuracy_Sn = meanAccuracy_Sn))
        CSV.write(join(["test_accuracy_", slice, ".csv"]),DataFrame(test_accuracy = test_accuracy))
        CSV.write(join(["meanAccuracy_test_", slice, ".csv"]),DataFrame(meanAccuracy_test = meanAccuracy_test))
        CSV.write(join(["th_names_", slice, ".csv"]),DataFrame(th_names = th_names))
        CSV.write(join(["ths_", slice, ".csv"]),DataFrame(ths = ths))
        CSV.write(join(["se_ths_", slice, ".csv"]),DataFrame(se_ths = se_ths))
        CSV.write(join(["dofs_", slice, ".csv"]),DataFrame(dofs = dofs))
        saveDataFrame(th_summary, join(["th_summary", slice,"_"]); path=savePath)
    catch
        warning("didn't work")
        cd(wd)
        rethrow()
    end
    cd(wd)
    saveDataFrame(results, join(["results_df",slice,"_"]); path=savePath)
    return results
end
function compare_AICBIC(meanAIC, AICs; yl="AIC", iters=0, ax=gca(), minmax="min")
#     ax.plot(1:length(meanAIC), meanAIC, "k.")
    CIls = []
    CIus = []
    for i=1:length(meanAIC)
        (CIl, CIu)=getBootCI(AICs[i]; alph=0.05)
        push!(CIls, CIl)
        push!(CIus, CIu)
    end
    plot_with_CI(meanAIC, CIls, CIus; ax=ax)
    
    # the best is the min
    if minmax=="min"
        bestix = findall(x->minimum(meanAIC)==x, meanAIC)
    else
        bestix = findall(x->maximum(meanAIC)==x, meanAIC)
    end
    plot(bestix, meanAIC[bestix], "g*", markersize=20)

    ax.set_xlabel("Model #")
    ax.set_xticks(collect(1:length(meanAIC)))
    ax.set_ylabel(yl)
    ax.set_title(join(["mean ", yl, " iters=", iters]))
end
function getBootCI(Vec::Vector; alph=0.05)
	Vec = sort!(Vec)
    n = length(Vec)
    minix = round(Int,(alph/2)*n)
    maxix = round(Int,(1-(alph/2))*n)
    if minix == 0 #|| isnan(minix)
    	minix = 1
    end
  #   if isnan(maxix)
		# maxix = 1
  #   end
    CImin = Vec[minix]
    CImax = Vec[maxix]
    # println("mean=",mean(Vec), " (", CImin, ", ", CImax, ")")
    return (CImin, CImax)
end