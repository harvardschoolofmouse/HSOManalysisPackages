#
#  The following functions allow us to give printed feedback to the user
#
#	pretty_print_list(myList)
# 	warning(msg::String)
# 

function pretty_print_list(myList; orient="vertical", digits=3, enum=false)
	# 
	if orient=="horizontal"
		for ii = 1:length(myList)
			i = myList[ii]
			if typeof(i)<:Number
				i = round(i, digits=digits)
			end
			if enum
				print(ii, ". ")
			end
	    	print(i)
	    	print("  ") 
	    end
	else
	    for ii = 1:length(myList)
	    	i = myList[ii]
	    	if typeof(i)<:Number
				i = round(i, digits=digits)
			end
			if enum
				print(ii, ". ")
			end
	    	println(i) 
	    end
    end
end;
function warning(msg::String)
    red = "\033[1m\033[31m"
    # println("\x1b[31m\"********************************************************************************************\"\x1b[0m")
    println(join(["\x1b[31m\"     WARNING: ", msg, "\"\x1b[0m"]))
    # println("\x1b[31m\"********************************************************************************************\"\x1b[0m")
end
function progressbar(iter,total)
    done = ["=", "=", "=", "=", "=", "=", "=", "=", "=", "="]
    incomplete = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-"]
    if mod(iter,total*0.1) == 0    
        ndone = round(Int,iter/total * 10);
        nincomp = round(Int, (1 - iter/total) * 10);
        println("   *", join(done[1:ndone]), join(incomplete[1:nincomp]), " (", iter, "/", total, ") ",Dates.format(now(), "mm/dd/yy HH:MM:ss") )
    elseif iter==1
        ndone = 0;
        nincomp = 10;
        println("   *", join(done[1:ndone]), join(incomplete[1:nincomp]), " (", iter, "/", total, ") ",Dates.format(now(), "mm/dd/yy HH:MM:ss") )
    end
end;

