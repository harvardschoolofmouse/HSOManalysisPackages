cd("/Users/ahamilos/Documents/GitHub/HSOManalysisPackages/JuliaLanguage_Hierarchical_Models")
include(join([pwd(), "/function_library/file_modules.jl"]))
refresh_tools("function_library")

runID = "dwnsample_trim-150_250smsbins_allSNc"
collatedPath = "/Volumes/Neurobio/Assad Lab/allison/Collated Analyses/JULIA_CSVs/All_SNc_nosplitsesh"
modelpackagefunction = bootlogit_timeslice_modelpackage1
postprocessingfunction = template_postprocessingfunction#bootlogit_timeslice_postprocessingfunction1

(fails150_b5_2, results150_b5_2,postprocessing150_b5_2) = run_collated_model(collatedPath, modelpackagefunction; pathIDx = [75,77],runFails=false, failDirs=[], postprocessingfunction=postprocessingfunction,
    compositesavepath="", runID=runID, suppressFigures=true); 