#
#  The following functions use GLMNet to do xval...
#	I wrote to be compatible and comparable to my own xvalidation
#
#
#
using GLMNet
using LaTeXStrings

function package_df_for_GLMNet(data, predictorIDs::Vector{Symbol}, yID::Symbol; modelClass="logit", add_th0=true)
	#
	#  modelClass=linear, logit
	#
    y = data[yID]
    if !isempty(predictorIDs)
	    X = data[predictorIDs[1]] #'
	    if length(predictorIDs)>1
	        for x = 2:length(predictorIDs)
	            X = hcat(X, data[predictorIDs[x]]) #'
	        end
	    end
    else
    	X = ones(length(y))
    	add_th0 = false
    	predictorIDs = [:th0]
    end
	predictors = Vector{String}(undef, 0)
	if add_th0
    	push!(predictors,"th0")
    	X = hcat(X,)
        # X = hcat(ones(length(X[1,:])),X)
    else
        X = hcat(X,)
	end
    for i=1:length(predictorIDs)
    	push!(predictors,string(predictorIDs[i]))
	end
	if modelClass == "logit"
		y = [string(x) for x in y]
	end
	predicted = string(yID)
	return (X, y, predictors, predicted)
end
function GLMNet_xval(all_set::DataFrame,train_set::DataFrame, validation_set::DataFrame, predictorIDs::Vector{Symbol}, yID::Symbol; modelClass="logit", add_th0=true, k=10)
	#
	#  modelClass=linear, logit
	#
	(X_all, y_all, predictors, predicted) = package_df_for_GLMNet(all_set, predictorIDs, yID; modelClass=modelClass, add_th0=add_th0)
	(X_Sn, y_Sn, predictors, predicted) = package_df_for_GLMNet(train_set, predictorIDs, yID; modelClass=modelClass, add_th0=add_th0)
	(X_val, y_val, _, _) = package_df_for_GLMNet(validation_set, predictorIDs, yID; modelClass=modelClass, add_th0=add_th0)
	#
	#	Do the xval 
	#
	cv = glmnetcv(X_Sn, y_Sn, standardize=false, intercept=add_th0, nfolds=k)
	lam_star = lambdamin(cv)
	if lam_star == 1
		cv = glmnetcv(X_Sn, y_Sn, lambda=[0, 10^-10, 10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0], standardize=false, intercept=add_th0, nfolds=k)
		lam_star = lambdamin(cv)
	end
	lambdas = cv.lambda
	meanloss = cv.meanloss
	stdloss = cv.stdloss
	betas = cv.path.betas
	
	if add_th0
		betas = vcat(cv.path.a0'.*ones(1,size(betas)[2]), betas)
	end
	k = cv.nfolds
	Rsqs = cv.path.dev_ratio


	lamix = findall(x->x==lam_star, lambdas)[1]
	beta_star = betas[:, lamix]
	val_loss_star = meanloss[lamix]
	val_std_loss_star = stdloss[lamix]
	val_Rsq = cv.path.dev_ratio[lamix]


	if modelClass == "logit"
		yFit_Sn = GLMNet.predict(cv, X_Sn, outtype = :prob);
		yFit_xval = GLMNet.predict(cv, X_val, outtype = :prob);
		yFit_all = GLMNet.predict(cv, X_all, outtype = :prob);
		# get prediction and correct/incorrect
		correct_Sn= train_set[yID].==[round(x) for x in yFit_Sn]
		accuracy_Sn = sum(correct_Sn)/length(y_Sn)
		println("Training accuracy: ", accuracy_Sn)

		correct_validation_set= validation_set[yID] .==[round(x) for x in yFit_xval]
		accuracy_validation_set = sum(correct_validation_set)/length(y_val)
		# println(y_val[1:10])
		# println(yFit_xval[1:10])
		# println(correct_validation_set[1:10])
		# println(length(y_val))
		println("Validation accuracy: ", accuracy_validation_set)
		correct_all = all_set[yID] .==[round(x) for x in yFit_all]
		accuracy_all = sum(correct_all)/length(y_all)
		println("Accuracy with regularized model on all data: ", accuracy_all)
		# We need to get the covar matrix to get coeff stats:
		#
		# 	However, it should be noted that there is no closed-form solution to the CVmat for 
		#    regularized logistic regression. As such, this will almost certainly be wrong...
		# 	 Not sure how we should deal with this. We could, I suppose, give an estimate of the
		# 	 parameter across many fits without giving the coefficient error.
		#.  Alternatively, I could use something like AIC/BIC and just report the fit model with
		#.   its own standard errors.    
		#.  For sufficiently small regularization, it seems you could use regularization to get 
		#.   the model selection, and then calculate the coefficient error on the unregularized model...
		#.  Could also bootstrap a coefficient I suppose by fitting a lot of models

		#. https://stats.stackexchange.com/questions/224796/why-are-confidence-intervals-and-p-values-not-reported-as-default-for-penalized
		# 		"Best answer: have a look at section 6 of the vignette for the penalized R package ("L1 and L2 Penalized Regression Models" Jelle Goeman, Rosa Meijer, Nimisha Chaturvedi, Package version 0.9-47), https://cran.r-project.org/web/packages/penalized/vignettes/penalized.pdf.
		# We don't get CIs or standard errors on the coefficients when we use penalized regression because they aren't meaningful. Ordinary linear regression, or logistic regression, or whatever, provides unbiased estimates of the coefficients. A CI around that point estimate, then, can give some indication of how point estimates will be distributed around the true value of the coefficient. Penalized regression, though, uses the bias-variance tradeoff to give us coefficient estimates with lower variance, but with bias. Reporting a CI around a biased estimate will give an unrealistically optimistic indication of how close the true value of the coefficient may be to the point estimate.
		# ("Penalized Regression, Standard Errors, and Bayesian Lassos" Minjung Kyung, Jeff Gilly, Malay Ghosh, and George Casella, Bayesian Analysis (2010) pages 369 - 412, discusses non-parametric (bootstrapped) estimates of p values for penalized regression and, if I understand correctly, they are not impressed. http://www.stat.ufl.edu/archived/casella/Papers/BL-Final.pdf)"
		#
		# Xcov = X
		# if add_th0
		# 	Xcov = hcat(ones(length(Xcov[1,:])),Xcov)
		# end
		# V = Diagonal(yFit_xval.*(1 .- yFit_xval))
		# CVmat = X'*V*X #'
	else
		yFit_train = GLMNet.predict(cv, X_Sn);
		yFit_xval = GLMNet.predict(cv, X_val);
		yFit_all = GLMNet.predict(cv, X_all);
		correct_validation_set = NaN
		accuracy_validation_set = NaN
		correct_Sn = NaN
		accuracy_Sn = NaN
		correct_all = NaN
		accuracy_all = NaN
		# CVmat = (XtX + (lam * (size(XtX)[1] * I)))^-1*XtX*(XtX + (lam* (size(XtX)[1] * I)))^-1;
	end

	d = length(beta_star)
	
	stats_df = DataFrame(
		th_names=predictors, 
		n=vec(length(y_val).*ones(d)), 
		th=vec(beta_star), 
		se_th=vec(nanmat(d,1)), 
		std_th=vec(nanmat(d,1)), 
		dof=vec(nanmat(d,1)), 
	    p_t=vec(nanmat(d,1)), 
	    CImin_t=vec(nanmat(d,1)),
	    CImax_t=vec(nanmat(d,1)), 
	    p_z=vec(nanmat(d,1)), 
	    CImin_z=vec(nanmat(d,1)), 
	    CImax_z=vec(nanmat(d,1)),
	    note=["xval model has no CIs. n=length of validation set" for _=1:d])

	result_df = DataFrame(
		predicted=[predicted], 
		predictors=[predictors],
		th=[beta_star], 
		X=[X_all], 
		y=[float([y=="true" for y in y_all])], 
		lam=[lam_star], 
		yFit=[yFit_all], 
	    CVmat=[nanmat(d,d)],
	    se_model=[NaN],
	    se_th=[nanmat(d,1)],
	    signifCoeff=[nanmat(d,1)],
	    ESS=[NaN],RSS=[NaN],
	    Rsq=[val_Rsq],
	    n_kfold_xval=[k],
	    correct_validation_set=[correct_validation_set],
		accuracy_validation_set=[accuracy_validation_set],
		correct_Sn = [correct_Sn],
		accuracy_Sn = [accuracy_Sn],
		correct_all = [correct_all],
		accuracy_all = [accuracy_all])
	
	model_summary = DataFrame(
		predictors=predictors, 
		Coeff=vec(beta_star), 
		StdError=vec(nanmat(d,1)),
		p=vec(nanmat(d,1)),
		CI95_lower=vec(nanmat(d,1)),
		CI95_upper=vec(nanmat(d,1)),
		xval = [true for _=1:d])

	xval_data = DataFrame(
		k=k, 
		all_lambdas=[lambdas], 
		lam_star=lam_star, 
    	mean_train_loss_ea_lam=[NaN], 
    	mean_test_loss_ea_lam=[meanloss], 
    	std_test_loss_ea_lam = [stdloss],
    	mean_train_Rsq_ea_lam=[NaN], 
    	mean_test_Rsq_ea_lam=[Rsqs],
    	lamix = lamix,
		beta_star = [beta_star],
		val_loss_star = [val_loss_star],
		val_std_loss_star = [val_std_loss_star],
		val_Rsq = [val_Rsq],
    	)

	# println(size(float([y=="true" for y in y_all])))
	# println(vec(X_all))
	figure(figsize=(10,3))
	ax1 = subplot(1,3,1)
	plot_fit_results(vec(X_all), float([y=="true" for y in y_all]), yFit_all; correct=correct_all, modelClass=modelClass, ax=ax1)
	title("All Data")
	ax2 = subplot(1,3,2)
	plot_fit_results(vec(X_Sn), float([y=="true" for y in y_Sn]), yFit_Sn; correct=correct_Sn, modelClass=modelClass, ax=ax2)
	title("Training/K-Fold Pool")
	ax3 = subplot(1,3,3)
	plot_fit_results(vec(X_val), float([y=="true" for y in y_val]), yFit_xval; correct=correct_validation_set, modelClass=modelClass, ax=ax3)
	title("Validation Pool")
	return (result_df, model_summary, stats_df, xval_data, X_all, float([y=="true" for y in y_all]),yFit_all,correct_all)
end


	