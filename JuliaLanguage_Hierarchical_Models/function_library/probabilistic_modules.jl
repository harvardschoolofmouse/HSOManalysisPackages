#
#  The following functions allow us to more slickly do probmods (not including any Gen functions!)
#
#	flip(p=0.5)
# 	flip(output1::Any, output2::Any, p::Float64)
# 	inverse_sample(ys::Vector{Float64}; consideration_interval=[0.,1.])
# 	get_beta_params(desired_mode)#mode, spread
# 
#  NOT INCLUDING ANY GEN FILES
#
using Distributions

function logmeanexp(scores::Vector{Float64})
    logsumexp(scores) - log(length(scores))
end;

function flip(p::Float64=0.5)
    rand() < p;
end;

# @gen function flip(output1::Any, output2::Any, p::Float64)
#     if ({:flip} ~ bernoulli(p))
#         return output1
#     else
#         return output2
#     end
# end;

function inverse_sample(ys::Vector{Float64}; consideration_interval=[0.,1.])
    ys_idxs_sort = sortperm(ys)
    ys_cdf_fx = StatsBase.ecdf(ys)
    ys_cdf = ys_cdf_fx(ys[ys_idxs_sort])
    xi = uniform(consideration_interval[1],consideration_interval[2])
    yi = findfirst(y -> y >= xi, ys_cdf)
    y0 = findfirst(y -> y >= consideration_interval[1], ys_cdf)
    y1 = findfirst(y -> y >= consideration_interval[2], ys_cdf)
    consideration_interval_std = StatsBase.std(ys[ys_idxs_sort[y0:y1]])
    selected_ys_index = ys_idxs_sort[yi]
    return (selected_ys_index, consideration_interval_std)
end

function get_beta_params(desired_mode)#mode, spread
    m = desired_mode
    a = 10
    b = (a*(1-m) + 2*m - 1)/m
    return (a, b)
end