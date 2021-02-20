#
#  The following functions allow us to plot data in useful ways
#
#	render_distribution(data, xl; bins=nothing, t="data", condition="", ax=nothing, normalization="prob", histtype="bar")
# 	extract_behavior_distribution(og_df::DataFrame)
# 	set_xaxes_same_scale(axs)
#   set_yaxes_same_scale(axs)
# 	plot_with_CI
#	drawLine
#	plotsignif
# 
using PyPlot

function render_xy(x, y; t="title", xl="xlabel", yl="ylabel", alph=nothing, ax=nothing)
	if isnothing(alph)
		if length(x) > 10000
			alph = 1/(20*log(length(x)))
		else
			alph = 0.2
		end
		# alph = 1/length(x)
	end
	if isnothing(ax)
		try
			ax = gca()
		catch
			figure
			ax = gca()
		end
	end
	ax.plot(x, y, "k.", alpha=alph)
	ax.set_title(t)
	ax.set_xlabel(xl)
	ax.set_ylabel(yl)
	return ax
end

function render_distribution(data, xl; bins=nothing, t="data", condition="", ax=nothing, normalization="prob", histtype="bar")
    if isnothing(ax)
        figure(figsize=(2,2))
        ax = gca()
    end
    if normalization=="prob"
        weights = ones(size(data))./length(data)
    else
        weights = ones(size(data))
    end
    if isnothing(bins)
        (n, binedges, _) = ax.hist(data, weights=weights,histtype=histtype)
    else
        (n, binedges, _) = ax.hist(data, bins=bins, weights=weights,histtype=histtype)
    end
    ax.set_xlabel(xl)
    ax.set_ylabel(join(["p(", t, "|", condition, ")"]))
    ax.set_xticks(minimum(data):10*10^(floor(log(abs(maximum(data))))):maximum(data))
    title(t)
    return n, binedges
end;


# function extract_behavior_distribution(og_df::DataFrame)
#     trialIDs = unique(og_df.TrialNo)
#     lick_times = []
#     for i in trialIDs
#         # get the datapoints corresponding to this trial
#         dp = findfirst(x->x==i, og_df.TrialNo)
#         push!(lick_times, og_df.LickTime[dp])
#     end
#     render_distribution(lick_times, [0,1], t="true normalized distribution", bins=50)
#     return lick_times
# end
# function extract_behavior_distribution(og_df::DataFrame)
# 	# NOW IN HAZARD MODULES
#     #
#     # If the licktime was normalized, need to correct this back to time in sec...
#     #
    
#     trialIDs = unique(og_df.TrialNo)
    
#     lick_times = []
#     dps = []
#     for i in trialIDs
#         # get the datapoints corresponding to this trial
#         dp = findfirst(x->x==i, og_df.TrialNo)
#         push!(dps, dp)
#         push!(lick_times, og_df.LickTime[dp])
#     end
#     min_time = minimum(og_df[dps[2:end].-1, :X])
#     max_time = maximum(og_df[dps[2:end].-1, :X])
#     lick_times = (lick_times .* max_time) .+ min_time
#     println("min licktime=", min_time, " max licktime=", max_time)
#     figure(figsize=(3,3))
#     render_distribution(lick_times, "lick time (s)", t="true normalized distribution", bins=50, ax=gca())
#     xticks(0:17)
#     xlim([0,17])
#     return lick_times
# end

function set_xaxes_same_scale(axs)
    M=[axs[1].get_xlim()[1] , axs[1].get_xlim()[2]]
    for i = 2:length(axs)
        m=axs[i].get_xlim()
        M[1] = minimum([M[1],m[1]])
        M[2] = maximum([M[2],m[2]])
    end
    for i=1:length(axs)
        axs[i].set_xlim(M)
    end
end;

function set_yaxes_same_scale(axs)
    M=[axs[1].get_ylim()[1] , axs[1].get_ylim()[2]]
    for i = 2:length(axs)
        m=axs[i].get_ylim()
        M[1] = minimum([M[1],m[1]])
        M[2] = maximum([M[2],m[2]])
    end
    for i=1:length(axs)
        axs[i].set_ylim(M)
    end
end;

function plot_with_CI(data, CImin, CImax; dataLabels="", ylab="wt", ax=gca())
    # data, CImin and CImax are vectors
    # println(data)
    if isempty(dataLabels)
        dataLabels = [join(["d", i]) for i=0:length(data)-1]
    end
    ax.plot([0,length(data)+1], [0,0], "k-")
    for i = 1:length(data)
    	# pretty_print_list(data[i], orient="horizontal")
    	# println("")
    	dd = nan2zero(data[i])
    	# println(".	")
    	# dd=data[i]
    	# pretty_print_list(dd, orient="horizontal")
        ax.plot(i, dd, "r.", markersize=30)
        drawLine([i,i], [CImin[i],CImax[i]])
        drawLine([i-0.2,i+0.2], [CImin[i],CImin[i]])
        drawLine([i-0.2,i+0.2], [CImax[i],CImax[i]])
        if checksign(CImin[i],CImax[i])
            plotsignif(i,CImax[i], ax=ax)
        end
    end
    # println(dataLabels)
#     xticks(ticks=[1:length(data)])
#     ax.set_xticklabels(dataLabels)
#     ax.set_ylabel(ylab)
end
function drawLine(xx,yy;style="k-",ax=gca())
    # p1 = [x1, y1], p2=[x2,y2]
#     xx = [p1[1], p2[1]]
#     yy = [p1[2], p2[2]]
    ax.plot(xx, yy, style)
end
function plotsignif(x,CImax; ax=gca())
    plot(x, CImax+0.025, "*", color="gold")
end