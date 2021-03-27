function zero_to_one(signal)
    l = minimum(signal);
    u = maximum(signal);
    normd = (signal.-l)./u
end;
function normalize_to_ref_signal(signal, ref_signal)
    println(size(ref_signal))
    l = minimum(ref_signal);
    u = maximum(ref_signal);
    # normalize signal
    zsig = zero_to_one(signal);
    nsig = zsig.*u.+l
end;
function inverse_sample(ys::Vector{Float64}; consideration_interval=[0.,1.])
	# need to handle nans...
	idxs = findall(x-> !isnan(x), ys)
	if length(idxs) != length(ys)
		warning(join(["inverse_sample of a vector with nans. Found ", length(idxs), " usable idxs out of ", length(ys), "total in vector, ys."]))
	end
	# ys = ys[idxs]
    ys_idxs_sort = sortperm(ys)
    ys_cdf_fx = StatsBase.ecdf(ys)
    ys_cdf = ys_cdf_fx(ys[ys_idxs_sort])
    xi = uniform(consideration_interval[1],consideration_interval[2])
    yi = findfirst(y -> y >= xi, ys_cdf)
    y0 = findfirst(y -> y >= consideration_interval[1], ys_cdf)
    y1 = findfirst(y -> y >= consideration_interval[2], ys_cdf)
    consideration_interval_std = StatsBase.std(ys[ys_idxs_sort[y0:y1]])
    selected_ys_index = ys_idxs_sort[yi]
    if isnan(consideration_interval_std)
    	# warning("ooops, our consideration_interval_std was nan. I suppose this happens if trial too short...")
    	# warning(join(["length(ys)=", length(ys)]))
    	# warning(join(["ys[ys_idxs_sort[y0:y1]]=", ys[ys_idxs_sort[y0:y1]]]))
	end
    return (selected_ys_index, consideration_interval_std)
end
function get_beta_params(desired_mode)#mode, spread
    m = desired_mode
    a = 10
    b = (a*(1-m) + 2*m - 1)/m
    return (a, b)
end
@dist rand_element(myList) = myList[uniform_discrete(1, length(myList))];