#
#  The following functions allow us to parse files, write files, etc
# 
#	using DelimitedFiles, DataFrames
#
# 	Classes:
# 		TrialData
# 
# 
#.  extract_data(path)
#	timestamp_now()
# 	printFigure(name; fig=gcf())
#.  check_imported_data(data::TrialData; idx=nothing)
# 	truncate_at_cue_and_lickp250(td::TrialData; cbuffer_s=0., lbuffer_s=0)

using DelimitedFiles
using DataFrames
using Dates
using CSV

function refresh_tools(path; exact=false)
    # Extracts dependencies from folder
    start_dir = pwd()
    try
        cd(path)
        a = readdir()
        print("Refreshed (", timestamp_now(), ") ")
        for i = 1:length(a)
            if occursin(".jl", a[i])
            	if exact
	                include(joinpath(path, a[i]))
                else
                	include(a[i])
                end
                print(a[i], " ")
            end
        end
        cd(start_dir)
    catch
        cd(start_dir)
        rethrow()
    end
end;


struct TrialData
    xdata::Array{Vector{Float64},1}
    ydata::Array{Vector{Float64},1}
    trialNo::Vector{Int}
    lickTime_s::Vector{Float64}
    path::String
    sessionCode::String
end


function extract_data(path; blmode=false, LOImode=false)
	#
	# Extracts single trial CSV data from folder
	#
    start_dir = pwd()
    cd(path)
    a = readdir()
    nx2 = length(a)
    xstart::Int=0
    ystart::Int=floor(Int, nx2/2)
    xs = []
    ys = []
    trialNo = Vector{Int}(undef, ystart)
    lickTime_s = Vector{Float64}(undef, ystart)
    for i = 1:ystart
    	xx=[]
    	yy=[]
    	tt=[]
    	try
	        xx = vec(readdlm(a[xstart+i], ',', Float64, '\n'));
	        yy = vec(readdlm(a[ystart+i], ',', Float64, '\n'));
	        tt = split(a[xstart+i], "_")
        catch
        	# explain the error
        	println("	")
        	println("\x1b[31m\"!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERROR\"\x1b[0m")
        	println("\x1b[31m\"	On path ", path, "...\"\x1b[0m")
        	println("\x1b[31m\"blmode=", blmode, "  LOImode=", LOImode, "\"\x1b[0m")
        	println("\x1b[31m\"i=", i, "  ystart=", ystart, "  xstart=", xstart, "\"\x1b[0m")
        	println("\x1b[31m\"Could not vec(readdlm(a[xstart+i], ',', Float64, '\n'));\"\x1b[0m")
        	println("\x1b[31m\"a[xstart+i-1]=", a[xstart+i-1], "\"\x1b[0m")
        	println("\x1b[31m\"a[ystart+i-1]=", a[ystart+i-1], "\"\x1b[0m")
        	println("\x1b[31m\"a[xstart+i]=", a[xstart+i], "\"\x1b[0m")
        	println("\x1b[31m\"a[ystart+i]=", a[ystart+i], "\"\x1b[0m")
        	# println(a)
        	rethrow()
        end
        if blmode # if baseline, we have extra baseline text on the filename
        	tNo = parse(Int, tt[3][2:end])
	        lt = parse(Float64, tt[4][3:end-1])
        elseif LOImode
        	tNo = parse(Int, tt[2][2:end])
	        lt = parse(Float64, tt[3][3:end-1])
        else 
	        tNo = parse(Int, tt[2][2:end])
	        lt = parse(Float64, tt[3][3:end-1])
        end
        push!(xs,xx)
        push!(ys,yy)
        trialNo[i] = tNo
        lickTime_s[i] = lt
    end
    if blmode
    	tt = split(a[xstart+1], "_")
    	sessionCode = join(tt[5:end])
	    sessionCode = sessionCode[1:end-4]

    else
    	tt = split(a[xstart+1], "_")
	    sessionCode = join(tt[4:end])
	    sessionCode = sessionCode[1:end-4]    	
    end
    data = TrialData(xs,ys,trialNo,lickTime_s,path, sessionCode)
    cd(start_dir)
    return data
end;




function timestamp_now()
    dform = Dates.DateFormat("yyyy-mm-dd_HHMM")
    return Dates.format(now(), dform)
end

function printFigure(name; fig=gcf(), figurePath=figurePath, verbose=false)
	timestamp = timestamp_now()
	ret_dir = pwd()
	if Sys.iswindows()
		# dlim="\\" #"
		# println("here")
		cd(figurePath)
		# println("here2")
		figname = join([name, "_", timestamp, ".eps"])
		# println("figname=", figname)
	else
		dlim="/"
    	figname = join([figurePath,dlim, name, "_", timestamp, ".eps"])
	end

    # gcf() = fig  
    # fig
    
    if verbose
	    println(figname)
    end
    # println("here3")
    fig.savefig(figname, transparent=true, format="eps")
    # println("here4")
    cd(ret_dir)
end

function check_imported_data(data::TrialData; idx=nothing)
    if length(data.xdata[end]) != length(data.ydata[end])
        error("data fields not matching")
    end
    if isnothing(idx)
        idx = rand(1:length(data.xdata))
    else
        idx = findall(x->x==idx, data.trialNo)
        if !isempty(idx)
	        idx = idx[1]
        else
        	idx = rand(1:length(data.xdata))
    	end
        println(idx)
    end
    xs = data.xdata[idx]
    ys = data.ydata[idx]
    tNo = data.trialNo[idx]
    lt = data.lickTime_s[idx]
    path = data.path
    seshcode = data.sessionCode
    figure(figsize=(5,3))
    plot(xs, ys, color="red", linewidth=2)
    xlabel("time (s)")
    xticks(collect(range(0, stop=maximum(xs), length=5)))
    ylabel("dF/F")
    title(join([seshcode," t=", tNo, " lt=", lt, "s"]));
    println("mean: ", mean(ys))
    println("median: ", median(ys))
    println(path)
    return idx
end;

function makeSessionDataFrame(data::TrialData; normalize=false, includeBL_LOI=false, baseline_data=[], LOI_data=[], include_history=false, history_spacing_s=0.25, n_hx_terms = 10)
	if n_hx_terms != 10
		error("not yet implemented for other than 10 blocks of hx... We have to define variables for each col in the table, so to make easier I am hard coding 10 for now")
	end
	# We will always keep 10 hx terms to make things easier...

	# useful data from the TrialData struct to add to dataframe
    xs = data.xdata
    ys = data.ydata

    bl_xs = baseline_data.xdata
    bl_ys = baseline_data.ydata

    loi_xs = LOI_data.xdata
    loi_ys = LOI_data.ydata

    tNo = data.trialNo
    lt = data.lickTime_s
    path = data.path
    seshcode = data.sessionCode
    if includeBL_LOI
    	meanbl = [mean(x) for x in bl_ys]
    	medianbl = [median(x) for x in bl_ys]
    	meanloi = [mean(x) for x in loi_ys]
    	medianloi = [median(x) for x in loi_ys]
	end

	# We would also like to include history
	#   The tricky thing will be having history before the cue...will need to draw on the LOI and baseline data for this...
	#.  NB we are assuming no trimming of the cue in the dataset!
	if include_history && !includeBL_LOI
		# we must have LOI and baseline to do full history... so throw error
		error("Must include baseline and LOI to do history analysis")
	elseif include_history
		# we need to concatenate and keep track of our time vector...
		xs_hx = [Vector{Float64}(undef,0) for _=1:length(xs)]
		ys_hx = [Vector{Float64}(undef,0) for _=1:length(xs)]

		for i = 1:length(ys)
			if xs[i][1] != 0.00
				warning(join(["Expected trial start at 0.00, but got at ", xs[i][1], " on iter ", i]))
			end
			# concatenate the baseline, LOI and ys...
			# check to see if the edges exactly equal, if so warn bc we have redundant timepoint
			if i==1
				if bl_ys[i][end]==loi_ys[i][1]#baseline_data.ydata[i][end] == LOI_data.ydata[i][1]
					warning(join(["The end of baseline is same as begin of LOI on tNo ", tNo[i], " iter=", i]))
					println(bl_ys[i][end])#baseline_data.ydata[i][end])
					println(loi_ys[i][1])#LOI_data.ydata[i][1])
				elseif loi_ys[i][end] == ys[i][1]#LOI_data.ydata[i][end] == ys[i][1]
					warning(join(["The end of LOI is same as begin of trial on tNo ", tNo[i], " iter=", i]))
					println(loi_ys[i][end])#LOI_data.ydata[i][end])
					println(ys[i][1])
				end
			end

			
			# println("	")
			# println("size of xbaseline=", length(baseline_data.xdata[i]))
			# println("size of ybaseline=", length(baseline_data.ydata[i]))

			yyy = bl_ys[i]#baseline_data.ydata[i]
			append!(ys_hx[i], yyy)
			# println("with baseline, ylen=", length(ys_hx[i]))
			yyy = loi_ys[i]#LOI_data.ydata[i]
			append!(ys_hx[i], yyy)
			# println("with LOI, ylen=", length(ys_hx[i]))
			yyy = ys[i]
			append!(ys_hx[i], yyy)
			# println("with trial, ylen=", length(ys_hx[i]))


			# handle the LOI etc. timestamps
			xxxx=loi_xs[i]#LOI_data.xdata[i]
			loistamps = -1 .*reverse(xxxx) .- 0.01
			
			xxxx=bl_xs[i]#baseline_data.xdata[i]
			blstamps = minimum(loistamps) .+ xxxx .- 0.01
			# println("bl_xs: [", blstamps[1], ",", blstamps[end], "]")
			# println("loi_xs: [", loistamps[1], ",", loistamps[end], "]")
			# println("cue_xs: [", xs[i][1], ",", xs[i][end], "]")
			append!(xs_hx[i], blstamps)
			# println("with baseline, xlen=", length(xs_hx[i]))
			xxxx=loistamps
			append!(xs_hx[i], xxxx)
			# println("with loi, xlen=", length(xs_hx[i]))
			xxxx=xs[i]
			append!(xs_hx[i], xxxx)
			# println("with trial, xlen=", length(xs_hx[i]))
			# println("nunique diffs=", length(unique(round.(xs_hx[i][2:end], digits=3) .- round.(xs_hx[i][1:end-1], digits=3) )), ", ntimes=", length(xs_hx[i]))
			# println(unique(xs_hx[i][2:end] .- xs_hx[i][1:end-1]))

			# Now, for purposes of normalization, we only want to include the times that could be in the history...
			# Thus, trim back from the left side by the maximum time that could be considered...
			
			# println(length(xs_hx[i]))
			# println(length(ys_hx[i]))
			mintime = -1. * history_spacing_s*n_hx_terms -0.01 # a little edge buffer...
			minidx = findall(x->x==mintime, xs_hx[i])[1]
			xs_hx[i] = xs_hx[i][minidx:end]#filter(e->e>mintime,xs_hx[i])
			ys_hx[i] = ys_hx[i][minidx:end]
			# println(length(xs_hx[i]))
			# println(length(ys_hx[i]))
			# println(xs_hx[i][1])
		end
	end

    
    
    # start by sorting the trial numbers so we can get prev trial info...
    sorted_trial_idxs = sortperm(tNo)
    tNo_sorted = tNo[sorted_trial_idxs]
    lt_sorted_nonorm = lt[sorted_trial_idxs]
    xs_sorted = xs[sorted_trial_idxs]
    
    if normalize
    	if include_history
    		warning("NB that we are normalizing history separate from in-trial dopamine...")
			v = reduce(vcat, ys_hx)
        	ys_hx = [normalize_vector_0_1(x, mx=maximum(v), mn=minimum(v)) for x in ys_hx]
		end
        v = reduce(vcat, ys)
        ys = [normalize_vector_0_1(x, mx=maximum(v), mn=minimum(v)) for x in ys]
        v = reduce(vcat, lt)
        lt = [normalize_vector_0_1(x, mx=maximum(v), mn=minimum(v)) for x in lt]
        if includeBL_LOI
        	v = reduce(vcat, meanbl)
        	meanbl = [normalize_vector_0_1(x, mx=maximum(v), mn=minimum(v)) for x in meanbl]
        	v = reduce(vcat, medianbl)
        	medianbl = [normalize_vector_0_1(x, mx=maximum(v), mn=minimum(v)) for x in medianbl]
        	v = reduce(vcat, meanloi)
        	meanloi = [normalize_vector_0_1(x, mx=maximum(v), mn=minimum(v)) for x in meanloi]
        	v = reduce(vcat, medianloi)
        	medianloi = [normalize_vector_0_1(x, mx=maximum(v), mn=minimum(v)) for x in medianloi]
		end
    end

    lt_sorted = lt[sorted_trial_idxs]
    ys_sorted = ys[sorted_trial_idxs]
    if include_history
    	xs_hx_sorted = xs_hx[sorted_trial_idxs]
    	ys_hx_sorted = ys_hx[sorted_trial_idxs]
    end
    if includeBL_LOI
		meanbl_sorted = meanbl[sorted_trial_idxs]
		medianbl_sorted = medianbl[sorted_trial_idxs]
		meanloi_sorted = meanloi[sorted_trial_idxs]
		medianloi_sorted = medianloi[sorted_trial_idxs]
    end
   

    yss = Array{Float64}(undef, 0)
    xss = Array{Float64}(undef, 0)
    time2lick = Array{Float64}(undef, 0)
    pcInterval = Array{Float64}(undef, 0)
    tNos = Array{Int64}(undef, 0)
    lts = Array{Float64}(undef, 0)
    # baseline info
    # if includeBL_LOI
    	meanbls = Array{Float64}(undef, 0)
    	medianbls = Array{Float64}(undef, 0)
    	meanLOIs = Array{Float64}(undef, 0)
    	medianLOIs = Array{Float64}(undef, 0)

    if include_history
    	hx1 = Array{Float64}(undef, 0)
    	hx2 = Array{Float64}(undef, 0)
    	hx3 = Array{Float64}(undef, 0)
    	hx4 = Array{Float64}(undef, 0)
    	hx5 = Array{Float64}(undef, 0)
    	hx6 = Array{Float64}(undef, 0)
    	hx7 = Array{Float64}(undef, 0)
    	hx8 = Array{Float64}(undef, 0)
    	hx9 = Array{Float64}(undef, 0)
    	hx10 = Array{Float64}(undef, 0)
    end
    # end
    #previous trial info
    ltm1s = Array{Float64}(undef, 0)
    rxnm1s = Array{Bool}(undef, 0)
    earlym1s = Array{Bool}(undef, 0)
    rewm1s = Array{Bool}(undef, 0)
    itim1s = Array{Bool}(undef, 0)
    # 2 trials back...
    ltm2s = Array{Float64}(undef, 0)
    rxnm2s = Array{Bool}(undef, 0)
    earlym2s = Array{Bool}(undef, 0)
    rewm2s = Array{Bool}(undef, 0)
    itim2s = Array{Bool}(undef, 0)
    #
    lickstate=Array{Bool}(undef, 0)
    #
    paths = Array{String}(undef, 0)
    sessioncodes = Array{String}(undef, 0)

    for i = 1:length(tNo_sorted)
        # get previous trial info. We will have some issues if there's an exclusion near begin of session, e.g., exc trial 1 or 2... OR if exclude multiple in a row...
        if i-1 <= 0#tNo_sorted[i] == 1 
         	# no previous trial! we have no outcome info...
            if normalize
                ltm1 = 0.0
            else
                ltm1 = 0.0
            end
            outcomes = [false, false, false, false]
            if normalize
                ltm2 = 0.0
            else
                ltm2 = 0.0
            end
            outcomes2 = [false, false, false, false]
        elseif i-2 <=0#tNo_sorted[i] == 2 
        	# no 2 trials back!
            if normalize
                ltm2 = 0.0
            else
                ltm2 = 0.0
            end
            outcomes2 = [false, false, false, false]
        	if i-1 >0 && tNo_sorted[i] - tNo_sorted[i-1] == 1#tNo_sorted[1] == 1 
        		# But we do have 1 trial back
                ltm1 = lt_sorted[i-1]
                outcomes = [false, false, false, false]
                if lt_sorted_nonorm[i-1] < 0.7
                    outcomes[1] = true
                elseif lt_sorted_nonorm[i-1] >= 0.7 && lt_sorted_nonorm[i-1] < 3.333
                    outcomes[2] = true
                elseif lt_sorted_nonorm[i-1] >= 3.333 && lt_sorted_nonorm[i-1] < 7
                    outcomes[3] = true
                else
                    outcomes[4] = true
                end
            else # no prev trials!
				if normalize
                    ltm1 = 0.0
                else
                    ltm1 = 0.0
                end
                outcomes = [false, false, false, false]
            end
        
		else # all other cases
            if tNo_sorted[i] - tNo_sorted[i-2] == 2
                # we have 2 trials back!
                ltm2 = lt_sorted[i-2]
                outcomes2 = [false, false, false, false]
                if lt_sorted_nonorm[i-2] < 0.7
                    outcomes2[1] = true
                elseif lt_sorted_nonorm[i-2] >= 0.7 && lt_sorted_nonorm[i-2] < 3.333
                    outcomes2[2] = true
                elseif lt_sorted_nonorm[i-2] >= 3.333 && lt_sorted_nonorm[i-2] < 7
                    outcomes2[3] = true
                else
                    outcomes2[4] = true
                end
            else
                # no 2 trials back! we have no outcome info...
                if normalize
                    ltm2 = 0.0
                else
                    ltm2 = 0.0
                end
                outcomes2 = [false, false, false, false]
            end
            
            # previous trial data...
            if tNo_sorted[i] - tNo_sorted[i-1] == 1
                # we have a prev trial!
                ltm1 = lt_sorted[i-1]
                outcomes = [false, false, false, false]
                if lt_sorted_nonorm[i-1] < 0.7
                    outcomes[1] = true
                elseif lt_sorted_nonorm[i-1] >= 0.7 && lt_sorted_nonorm[i-1] < 3.333
                    outcomes[2] = true
                elseif lt_sorted_nonorm[i-1] >= 3.333 && lt_sorted_nonorm[i-1] < 7
                    outcomes[3] = true
                else
                    outcomes[4] = true
                end
            elseif tNo_sorted[i] - tNo_sorted[i-1] == 2
                # no previous trial! But we do have 2 trials back...
                if normalize
                    ltm1 = 0.0
                else
                    ltm1 = 0.0
                end
                outcomes = [false, false, false, false]
                # we have 2 trials back!
                ltm2 = lt_sorted[i-1]
                outcomes2 = [false, false, false, false]
                if lt_sorted_nonorm[i-1] < 0.7
                    outcomes2[1] = true
                elseif lt_sorted_nonorm[i-1] >= 0.7 && lt_sorted_nonorm[i-1] < 3.333
                    outcomes2[2] = true
                elseif lt_sorted_nonorm[i-1] >= 3.333 && lt_sorted_nonorm[i-1] < 7
                    outcomes2[3] = true
                else
                    outcomes2[4] = true
                end
            else
                # no previous trial! we have no outcome info...
                if normalize
                    ltm1 = 0.0
                else
                    ltm1 = 0.0
                end
                outcomes = [false, false, false, false]
            end         
        end

        for j = 1:length(ys_sorted[i])
            # extract each timepoint and signal value and relevant data
            push!(yss, ys_sorted[i][j])
            push!(xss, xs_sorted[i][j])
            if !normalize
                push!(time2lick, lt_sorted[i]-xs_sorted[i][j])
                if !isnan(xs_sorted[i][j]/lt_sorted[i])
                    push!(pcInterval, xs_sorted[i][j]/lt_sorted[i])
                else
                    push!(pcInterval, 1.)
                end
            else
                push!(time2lick, lt_sorted_nonorm[i]-xs_sorted[i][j])
                if !isnan(xs_sorted[i][j]/lt_sorted_nonorm[i])
                    push!(pcInterval, xs_sorted[i][j]/lt_sorted_nonorm[i])
                else
                    push!(pcInterval, 1.)
                end
                
            end
            push!(tNos, tNo_sorted[i])
            push!(lts, lt_sorted[i])
            push!(paths, path)
            push!(sessioncodes, seshcode)

            if includeBL_LOI
		    	push!(meanbls, meanbl_sorted[i])
		    	push!(medianbls, medianbl_sorted[i])
		    	push!(meanLOIs, meanloi_sorted[i])
		    	push!(medianLOIs, medianloi_sorted[i])
	    	else
	    		push!(meanbls, NaN)
		    	push!(medianbls, NaN)
		    	push!(meanLOIs, NaN)
		    	push!(medianLOIs, NaN)
	    	end

	    	if include_history
	    		debugg = false
	    		# here we need to take a mean over the history vector at a timepoint hx_n*history_spacing_s back in time...
	    		#
	    		# 	Gather the mean hx terms...
	    		#
	    		current_x = xs_sorted[i][j]
	    		if debugg && i==1 && j==1
	    			println("trialNo: ", tNo_sorted[i])
	    			println("current_x=",current_x)
	    			println("current_idx=",findfirst(x->x>=current_x, xs_hx_sorted[i]))
	    			println("current_y=",ys_sorted[i][j])
	    			println("current_y_hxnormd=",ys_hx_sorted[i][findfirst(x->x>=current_x, xs_hx_sorted[i])])
    			end
	    		hx_vect = []
	    		for k=1:n_hx_terms
	    			left_x = current_x - k*history_spacing_s

	    			# find center idx
	    			left_idx = findfirst(x->x>=left_x, xs_hx_sorted[i])
	    			nsampsdiff = floor(Int,history_spacing_s/(xs_sorted[1][2]-xs_sorted[1][1]) - (xs_sorted[1][2]-xs_sorted[1][1]))
	    			range_idx = collect(left_idx:left_idx+nsampsdiff)
					if debugg && i==1 && j==1
	    				println("size of ys_hx[i]", size(ys_hx_sorted[i]))
    				end
    				try 
		    			hxterm = mean(ys_hx_sorted[i][range_idx])
	    			catch
	    				warning("got error calculating hxterm!!!!!!")
	    				println("trialNo: ", tNo_sorted[i])
		    			println("current_x=",current_x)
		    			println("current_idx=",findfirst(x->x>=current_x, xs_hx_sorted[i]))
		    			println("current_y=",ys_sorted[i][j])
		    			println("current_y_hxnormd=",ys_hx_sorted[i][findfirst(x->x>=current_x, xs_hx_sorted[i])])
		    			println("left_x", k,"=",left_x)
		    			println("	left_idx=", left_idx)
		    			println("	range_idxs=[", range_idx[1], ":", range_idx[end],"]")
		    			# println("	hxterm_min=", minimum(ys_hx[i][range_idx]))
		    			# println("	hxterm_max=", maximum(ys_hx[i][range_idx]))
		    			rethrow()
		    			# println("	hxterm=", mean(ys_hx[i][range_idx]))
	    			end
	    			hxterm = mean(ys_hx_sorted[i][range_idx])
	    			push!(hx_vect, hxterm)
					if debugg && i==1 && j==1
		    			println("left_x", k,"=",left_x)
		    			println("	left_idx=", left_idx)
		    			println("	range_idxs=[", range_idx[1], ":", range_idx[end],"]")
		    			println("	hxterm_min=", minimum(ys_hx_sorted[i][range_idx]))
		    			println("	hxterm_max=", maximum(ys_hx_sorted[i][range_idx]))
		    			println("	hxterm=", mean(ys_hx_sorted[i][range_idx]))
	    			end
	    		end
	    		if debugg && i==1 && j==1
	    			println("hx_vect=", round.(hx_vect, digits=2))
	    		end
	    		push!(hx1, hx_vect[1])
	    		push!(hx2, hx_vect[2])
	    		push!(hx3, hx_vect[3])
	    		push!(hx4, hx_vect[4])
	    		push!(hx5, hx_vect[5])
	    		push!(hx6, hx_vect[6])
	    		push!(hx7, hx_vect[7])
	    		push!(hx8, hx_vect[8])
	    		push!(hx9, hx_vect[9])
	    		push!(hx10, hx_vect[10])
	    	end
            
            push!(ltm1s, ltm1)
            push!(rxnm1s, outcomes[1])
            push!(earlym1s, outcomes[2])
            push!(rewm1s, outcomes[3])
            push!(itim1s, outcomes[4])
            
            push!(ltm2s, ltm2)
            push!(rxnm2s, outcomes2[1])
            push!(earlym2s, outcomes2[2])
            push!(rewm2s, outcomes2[3])
            push!(itim2s, outcomes2[4])
            if j==length(ys_sorted[i])
                push!(lickstate, true)
            else
                push!(lickstate, false)
            end
        end     
       
    end
    
    # normalize the time2lick variable
    if normalize
	    time2lick = zero_to_one(time2lick)
    end

    y_shuffle = Random.shuffle(yss)
    
    if include_history
    	df = DataFrame(DataID=1:length(xss), 
	        TrialNo = tNos,
	        X=xss, 
	        Y=yss, 
	        LickTime=lts,
	        Time2Lick=time2lick,
	        PcInterval=pcInterval,
	        LickTime_1back=ltm1s, 
	        LickTime_2back=ltm2s, 
	        Mean_Baseline = meanbls,
			Median_Baseline = medianbls,
			Mean_LOI = meanLOIs,
			Median_LOI = medianLOIs,
	        Rxn_1back=rxnm1s,
	        Early_1back=earlym1s,
	        Reward_1back=rewm1s,
	        ITI_1back=itim1s,
	        Rxn_2back=rxnm2s,
	        Early_2back=earlym2s,
	        Reward_2back=rewm2s,
	        ITI_2back=itim2s,
	        LickState = lickstate,
	        Yshuffle = y_shuffle,
	        Hx1 = hx1,
	        Hx2 = hx2,
	        Hx3 = hx3,
	        Hx4 = hx4,
	        Hx5 = hx5,
	        Hx6 = hx6,
	        Hx7 = hx7,
	        Hx8 = hx8,
	        Hx9 = hx9,
	        Hx10 = hx10,
	        history_spacing_s = [history_spacing_s for _=1:length(hx1)],
	        SessionCode=sessioncodes, 
	        Path=paths)
    else
	    df = DataFrame(DataID=1:length(xss), 
	        TrialNo = tNos,
	        X=xss, 
	        Y=yss, 
	        LickTime=lts,
	        Time2Lick=time2lick,
	        PcInterval=pcInterval,
	        LickTime_1back=ltm1s, 
	        LickTime_2back=ltm2s, 
	        Mean_Baseline = meanbls,
			Median_Baseline = medianbls,
			Mean_LOI = meanLOIs,
			Median_LOI = medianLOIs,
	        Rxn_1back=rxnm1s,
	        Early_1back=earlym1s,
	        Reward_1back=rewm1s,
	        ITI_1back=itim1s,
	        Rxn_2back=rxnm2s,
	        Early_2back=earlym2s,
	        Reward_2back=rewm2s,
	        ITI_2back=itim2s,
	        LickState = lickstate,
	        Yshuffle = y_shuffle,
	        SessionCode=sessioncodes, 
	        Path=paths)
    end
    return df
end;

function truncate_at_cue_and_lickp250(td::TrialData; cbuffer_s=0., lbuffer_s=0,bl_data=[], loi_data=[])
    newxs = []
    newys = []
    trialNo = []
    lickTime_s = []

    if typeof(loi_data) == TrialData
    	newxs_bl = []
	    newys_bl = []
	    trialNo_bl = []
	    lickTime_s_bl = []
	    newxs_loi = []
	    newys_loi = []
	    trialNo_loi = []
	    lickTime_s_loi = []
    end
    
    for i = 1:length(td.xdata)
        xx = td.xdata[i]
        cidx = findfirst(x->x>=0. +cbuffer_s, xx)
        yy = td.ydata[i]
        lidx = findlast(x->x<=(td.lickTime_s[i]+lbuffer_s), xx)

        if typeof(loi_data) == TrialData
	    	xx_bl = bl_data.xdata[i]
	    	yy_bl = bl_data.ydata[i]
	    	xx_loi = loi_data.xdata[i]
	    	yy_loi = loi_data.ydata[i]
	    end
        if !isnothing(lidx)
	        xx = xx[1:lidx]
	        xx = xx[cidx:end]
	        yy = yy[1:lidx]
	        yy = yy[cidx:end]
	        push!(newxs, xx)
	        push!(newys, yy)
	        push!(trialNo, td.trialNo[i])
	        push!(lickTime_s, td.lickTime_s[i])
	        if typeof(loi_data) == TrialData
	        	push!(newxs_bl, xx_bl)
	        	push!(newxs_loi, xx_loi)
		        push!(newys_bl, yy_bl)
		        push!(newys_loi, yy_loi)
		        push!(trialNo_bl, bl_data.trialNo[i])
		        push!(trialNo_loi, loi_data.trialNo[i])
		        push!(lickTime_s_bl, bl_data.lickTime_s[i])
		        push!(lickTime_s_loi, loi_data.lickTime_s[i])
	        end
        end
    end
    trial_data = TrialData(newxs, newys, trialNo, lickTime_s, td.path, td.sessionCode)
    if typeof(loi_data) == TrialData
    	new_bl_data = TrialData(newxs_bl, newys_bl, trialNo_bl, lickTime_s_bl, bl_data.path, bl_data.sessionCode)
    	new_loi_data = TrialData(newxs_loi, newys_loi, trialNo_loi, lickTime_s_loi, loi_data.path, loi_data.sessionCode)
    	return (trial_data, new_bl_data, new_loi_data)
	else
	    return trial_data 
    end
end;

function saveDataFrame(df, filename; path=nothing)
	if isnothing(path)
		path = pwd()
	end
	if Sys.iswindows()
		dlim="\\" #"
	else
		dlim="/"
	end
	file = join([path, dlim, filename, timestamp_now(), ".csv"])
	println("saving... \n", file)
	CSV.write(file, df)
	return file
end

function importDataFrame(filename; suppresswarn=false)
	# NB! If you have nested DFs or arrays in the df, you need to instead use parse_imported_modelresult(filename)
	if !suppresswarn
		warning("NB! If you have nested DFs or arrays in the df, you need to instead use parse_imported_modelresult(filename)")
		println("Loading... \n", filename)
	end
	df = CSV.read(filename)
	return df
end

function unwrap_CSV_modelresult(filename)
	println("Loading... \n", filename)
	# empty_similar_df is a function to construct an empty dataframe

	# [eval(Meta.parse((f.train_ths[i][4:end]))) for i = 1:length(f.train_ths)]
	f = CSV.File(filename);
	colnames = propertynames(f[1,:][1])
	df_new = DataFrame()
	for i = 1:length(colnames)
	    # get the column data from CSV
	    coldata = [f[ii][i] for ii=1:length(f)]
	    if typeof(coldata[i]) == String
	        try
	            df_new[colnames[i]] = [eval(Meta.parse((coldata[j][4:end]))) for j = 1:length(coldata)] #[1,2]
	        catch
	#             println([eval(Meta.parse((coldata[j][4:end]))) for j = 1:length(coldata)])
	            df_new[colnames[i]] = coldata
	        end
	    else
       	 	df_new[colnames[i]] = coldata
	    end
    end
    return df_new
end



function fixStringForOS(s::String)
	
end
