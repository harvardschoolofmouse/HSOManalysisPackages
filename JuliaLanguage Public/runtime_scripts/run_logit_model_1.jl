# run_logit_model_1.jl

# USER MUST SET THESE:
collatedPath = "/Volumes/Neurobio/Assad Lab/allison/Collated Analyses/JULIA_CSVs/All_SNc_nosplitsesh"
compositesavepath="" # the save directory will be formed automagically!
runID = 1 # set runID nonZero to specify which run it is.

# Use the latest versions of our files and dependencies!
cd("/Users/lilis/Documents/GitHub/HSOManalysisPackages/JuliaLanguage_Hierarchical_Models")
include("/Users/lilis/Documents/GitHub/HSOManalysisPackages/JuliaLanguage_Hierarchical_Models/function_library/file_modules.jl")
refresh_tools("/Users/lilis/Documents/GitHub/HSOManalysisPackages/JuliaLanguage_Hierarchical_Models/function_library", exact=true)
println("	")
# CANNED FOR LOGIT:
modelpackagefunction = bootlogit_modelpackage1
postprocessingfunction = bootlogit_postprocessingfunction1

warning("suppressing figures from individual sessions, but they will still be saved to host folders for that session.")
println("	")
(fails, results) = run_collated_model(collatedPath, modelpackagefunction; pathIDx = [],runFails=false, failDirs=[], postprocessingfunction=postprocessingfunction,
    compositesavepath=compositesavepath, runID=runID, suppressFigures=true); 