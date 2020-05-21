# HSOManalysisPackages
Contains MATLAB 2018B analysis files used in Hamilos et al., 2020 with a sample dataset, instructions, and all necessary custom dependencies.

**Please note!! This supporting instructional documentation is under construction.**

Last Update on May 21, 2020. Please contact ahamilos{at}g.harvard.edu if you have questions or something isn't working. 

  > For example datasets to use with the sample code, please see:
  > https://www.dropbox.com/sh/wdotym743hmo4jc/AABKfTVxkH2JVkZXJ-7UpLyCa?dl=0

-------------------------------------------------
INSTRUCTIONS: Contents:
  
  0. Importing repo code
  1. Loading sample photometry object
  2. Code to reproduce figures from Hamilos et al., 2020
  3. Instructions for generating new objects from raw CED datasets

-------------------------------------------------
0. Importing repo code (designed for Windows 7 and MacOS 10.14)

    1. Clone or download the HSOManalysisPackages repository and add subfolders files to path in MATLAB 2018B or higher
    2. Make sure the following toolboxes are enabled in MATLAB: 
    
          - Curve Fitting Toolbox v3.5.8+
          - Statistics and Machine Learning Toolbox v11.4+
          - Signal Processing Toolbox v8.1+
          - Optimization Toolbox v8.2+
          - Image Processing Toolbox v10.3+
          - DSP System Toolbox v9.7+
          - Control System Toolbox v10.5+
    
-------------------------------------------------
1. Loading single-session sample photometry object:

  - Using the sample datasets at https://www.dropbox.com/sh/wdotym743hmo4jc/AABKfTVxkH2JVkZXJ-7UpLyCa?dl=0
  
    - Navigate to the sample datafolder for animal B5, SNc day 13 (B5_SNc_13)
  
    - Open the file. It will initialize to the workspace as sObj or obj, depending on the dataset
    
Loading composte session object:

  - Navigate to the sample Composite Datasets folder and open the file with the analysis of choice

-------------------------------------------------
2. Code to reproduce figures from Hamilos et al., 2020

**Warning: under construction: 5/21/2020 11:21 EST -- Please check back tomorrow**

  - FIGURE 1: Use the Figure 1 dataset filefolder from dropbox
  
    - 1b: Histograms
            
            load('CollatedStatAnalysisObj_cdf_20191107_18_10_runIDno4486.mat')
            obj.plot('hxg-xsesh')
    - 1c-f: CLTA (cue/lick triggered averages)
    
         NB: Example dataset has 2 mice recorded at SNc and shows basic pattern found more cleanly in the full average across all animals (12-mouse full SNc/tdt dataset and 4-mouse DA2m datasets are too large to put on Dropbox, available on reasonable request to Allison Hamilos, ahamilos{at}g.harvard.edu). The all-mice figures are saved to the 1C-F folders
          
          To plot CLTAs of composite datasets with 500 ms timebins for averaging:
          
          obj = CLASS_HSOM_photometry('v3x', 'times', 68, {'box', 200000}, 30000, [], [], 'off')
          %   Select either SNc or tdt to make a composite SNc or tdt object. 
          %     Then select the corresponding HOST folder from the 
          %     Dropbox sample data.
          
          %     Plotting: NB: you can use the plot function on single-session 
          %       datasets (e.g., B5_SNc_13) and composite datasets
          obj.plot('CLTA', [6, 9, 13, 15], false, 10, 'last-to-first', 1) 
          %
          %     NB: [6,9,13,15] are the timebins -- you can plot any between 1-68
          %
          xlim([-1.5,5])
          %
          % Example figure files included in Dropbox folder
          %
          % See CLASS_HSOM_photometry > plot documentation for instructions 
          %   to parameterize the plot function, there are several 
          %   additional options and functionalities you can try!
          
          
   - FIGURE 2:
     

-------------------------------------------------
3. To generate new objects from CED datasets:

    i. Save CED/Spike2 dataset as .mat file. Ensure that the CED filename has the same name as the .mat file to avoid errors. Ensure waveforms and times are saved and use channel names.
  
    ii. Put raw .mat files into directories in the following way:

    HOSTdirectory >
    
                MOUSENAME_SIGNAL_DAY#
                e.g.:
                - B5_SNc_13
                - B5_SNc_15
                
    Each session folder must contain a CED .mat file and an exclusions text file (.txt) with a filename containing the word "exclusions." Exclusions (for trials) may be written as numbers delimited by any character EXCEPT a dash (-). A dash character denotes a range of trials (e.g., 4-6 = 4,5,6). Any non-existant trials will not affect the file, e.g., excluding trials 400-1000 on a session with 420 trials will be interpreted as 400-420.

    For a standard photometry object, use:
    
        obj = CLASS_HSOM_photometry('v3x', MODE, NBINS, {'GFITMODE', GFITPARAM}, NSAMPLEBUFFER, [], [], 'off')
    MODE: either 'times' or 'trials' -- indicates whether signals should be binned by time in trial or for bins with even numbers of trials in each bin
  
    NBINS: a number. 
    
          For 'times' MODE: Using 17 for the standard task with 'times' results in 17 1-second timebins.   
                  
          For 'trials' MODE: enter the number of bins you wish. If number exceeds number of trials, object will not intialize.
                  
  GFITMODE, GFITPARAM: selects the signal dF/F method (or filtering for movement channels). In general:
  
        {'box', 200000} -- uses the moving average dF/F method with 200000 sample window (200 s at 1kHz sampling)
        {'EMG', []} -- applies default rectification of EMG signal
        {'X', []} -- applies default bandpass filtering followed by rectification for accelerometer signals
        {'CamO', 30} -- processes camera data sampled at 30 Hz
        
  NSAMPLEBUFFER: the number of samples buffered into the ITI. Standard is 30000 for photometry data, 100000 for movement control data


  Standard init for photometry:
  
    obj = CLASS_HSOM_photometry('v3x', 'times', 17, {'box', 200000}, 30000, [], [], 'off')
  Standard init for EMG:
  
    obj = CLASS_HSOM_photometry('v3x', 'times', 17, {'EMG', []}, 100000, [], [], 'off')
  
   NB: data can be rebinned anytime for single session objects. Composite objects (averaging across sessions) cannot be changed after the initialization step.


--------------------------------------

License Notes:

   1. Note that linspecer v1.4.0.0+ by Jonathan C. Lansey is available from MATLAB FileExchange and included in dependencies with its original license: 
    
        https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap?focused=5372538&tab=function

