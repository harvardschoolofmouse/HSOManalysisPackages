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
    	try
	        xx = vec(readdlm(a[xstart+i], ',', Float64, '\n'));
	        yy = vec(readdlm(a[ystart+i], ',', Float64, '\n'));
	        tt = split(a[xstart+i], "_")
        catch
        	# explain the error
        	println("	")
        	println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERROR")
        	println("	On path ", path, "...")
        	println("blmode=", blmode, "  LOImode=", LOImode)
        	println("i=", i, "  ystart=", ystart, "  xstart=", xstart)
        	println("Could not vec(readdlm(a[xstart+i], ',', Float64, '\n'));")
        	println("a:")
        	println(a)
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
	if Sys.iswindows()
		dlim="\\" #"
	else
		dlim="/"
	end

    # gcf() = fig  
    # fig
    timestamp = timestamp_now()
    figname = join([figurePath,dlim, name, "_", timestamp, ".eps"])
    if verbose
	    println(figname)
    end
    fig.savefig(figname, transparent=true, format="eps")
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

function makeSessionDataFrame(data::TrialData; normalize=false, includeBL_LOI=false, baseline_data=[], LOI_data=[])
    # useful data from the TrialData struct to add to dataframe
    xs = data.xdata
    ys = data.ydata
    tNo = data.trialNo
    lt = data.lickTime_s
    path = data.path
    seshcode = data.sessionCode
    if includeBL_LOI
    	meanbl = [mean(x) for x in baseline_data.ydata]
    	medianbl = [median(x) for x in baseline_data.ydata]
    	meanloi = [mean(x) for x in LOI_data.ydata]
    	medianloi = [median(x) for x in LOI_data.ydata]
	end
    
    # start by sorting the trial numbers so we can get prev trial info...
    sorted_trial_idxs = sortperm(tNo)
    tNo_sorted = tNo[sorted_trial_idxs]
    lt_sorted_nonorm = lt[sorted_trial_idxs]
    xs_sorted = xs[sorted_trial_idxs]
    
    if normalize
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
        # get previous trial info
        if tNo_sorted[i] != 1 && tNo_sorted[1]==1
            if tNo_sorted[i] != 2
                # 2 trials back data...
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
                        ltm2 = 1.0
                    else
                        ltm2 = 17.1
                    end
                    outcomes2 = [false, false, false, false]
                end
            else
                # no 2 trials back! we have no outcome info...
                if normalize
                    ltm2 = 1.0
                else
                    ltm2 = 17.1
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
                    ltm1 = 1.0
                else
                    ltm1 = 17.1
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
                    ltm1 = 1.0
                else
                    ltm1 = 17.1
                end
                outcomes = [false, false, false, false]
            end
        else
            # no previous trial! we have no outcome info...
            if normalize
                ltm1 = 1.0
            else
                ltm1 = 17.1
            end
            outcomes = [false, false, false, false]
            if normalize
                ltm2 = 1.0
            else
                ltm2 = 17.1
            end
            outcomes2 = [false, false, false, false]
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
    return df
end;

function truncate_at_cue_and_lickp250(td::TrialData; cbuffer_s=0., lbuffer_s=0)
    newxs = []
    newys = []
    trialNo = []
    lickTime_s = []
    
    for i = 1:length(td.xdata)
        xx = td.xdata[i]
        cidx = findfirst(x->x>=0. +cbuffer_s, xx)
        yy = td.ydata[i]
        lidx = findlast(x->x<=(td.lickTime_s[i]+lbuffer_s), xx)
        if !isnothing(lidx)
	        xx = xx[1:lidx]
	        xx = xx[cidx:end]
	        yy = yy[1:lidx]
	        yy = yy[cidx:end]
	        push!(newxs, xx)
	        push!(newys, yy)
	        push!(trialNo, td.trialNo[i])
	        push!(lickTime_s, td.lickTime_s[i])
        end
    end
    trial_data = TrialData(newxs, newys, trialNo, lickTime_s, td.path, td.sessionCode)
    return trial_data 
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
