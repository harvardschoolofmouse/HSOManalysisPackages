#
#  The following functions allow us to do math on vectors and NaNs
# 
#   dependencies: StatsBase
#
#   nanmean(x::Vector{Float64})
#   nansum(x::Vector{Float64})
#   nanmat(r, c)
#   nan_elementwise_multiply(x::Vector{Float64}, y::Vector{Float64})
#   nan_elementwise_sum(x::Vector{Float64}, y::Vector{Float64})
#   mean(y::Vector)
#   std(x::Vector; force_nan_zero=false)
#   checksign(p1, p2)

using StatsBase

function nanmean_mat(x::Array{Float64}, dim=1)
    if dim == 1
        iter = 2
    else
        iter = 1
    end
    nanmm = nanmat(size(x)[iter],1)
    for i = 1:size(x)[iter]
        if dim == 1
            v = nanmean(x[:,i])
        else
            v = nanmean(x[i,:])
        end    
        nanmm[i] = v
    end
    return nanmm
end
function nanmean(x::Vector{Float64})
    xx = filter(!isnan, x)
    sum(xx)/length(xx)
end
function nansum(x::Vector{Float64})
    xx = filter(!isnan, x)
    sum(xx)
end;
function nan_elementwise_multiply(x::Vector{Float64}, y::Vector{Float64})
    xx = findall(x->isnan(x), x)
    x[xx] .= 1;
    yy = findall(y->isnan(y), y)
    y[yy] .= 1;
    x.*y
end
function nan_elementwise_sum(x::Vector{Float64}, y::Vector{Float64})
    xx = findall(x->isnan(x), x)
    x[xx] .= 0;
    yy = findall(y->isnan(y), y)
    y[yy] .= 0;
    x.+y
end
function nanmat(r,c)
    NaN.*ones(r,c)
end
function nan2zero(vec)
    if !(typeof(vec) <: Array)
        if isnan(vec)
            return 0.
        else
            return vec
        end
    else
        ix = findall(x->isnan(x), vec)
        vec[ix] = zeros(length(ix))
        return vec
    end
end
# function mean(y)
#     sum(y)/length(y)
# end
function nanstd(x; force_nan_zero=false)
    y = StatsBase.std(x)
    if isnan(y) && force_nan_zero
        y=0
    else
        return y
    end
end

function checksign(d1,d2)
   sign(d1)==sign(d2)
end
