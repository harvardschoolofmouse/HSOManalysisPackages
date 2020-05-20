classdef CLASS_HSOM_optogenetics_legacy_obj < handle
	% 	
	% 	GITHUB HSOM VERSION
	% 
	% 	Created 	10/25/19	Allison Hamilos (ahamilos{at}g.harvard.edu)  
	% 	Modified 	10/29/19	Allison Hamilos (ahamilos{at}g.harvard.edu)  
	% 
	% 	VERSION CODE: ['CLASS_HSOM_optogenetics_legacy_obj v2.0 Modified 10-29-19 18:36 | obj created: ' datestr(now)];
	% 
	% 	VERSION 1.x
	% 	1	- 5/19/20:	Version used for figures in Hamilos et al., 2020
	% 
	% ----------------------------------------------------------------
	% 
	% 	Compiles and analyzes pre-processed optogenetics datasets: both single and compiled across sessions
	%
	% 
	properties
		iv
		time
		lick
		h
		m		
		lick_data_struct
		valid_blocks
		bootStat
	end

	
	methods
		%-------------------------------------------------------
		%		Methods: Initialization
		%-------------------------------------------------------
		function obj = CLASS_HSOM_optogenetics_legacy_obj(Mode, stim_protocol, stimtype, ChR2data, Xdata, collateKey)
			% 
			% 	Mode:				default: for a single session object
			% 						batch: generate a bunch of single session objects at once
			% 						collate: using a bunch of existing stat objs, opens them, computes an analysis, and then saves them and the analysis results
			% 
			% 	collateKey:			specifies the kind of analysis to capture in each object
			% 						bootEmery -- computes the Emery Bootstrap results for all objects in the HOSTobject path
			% 						ntRaster -- combines raster data across sessions
			% 						ecdf -- saves ecdf figures and does ks test on each one
			% 							obj = CLASS_HSOM_optogenetics_legacy_obj('collate', 'A', 'A', [], [], 'ecdf')
			% 
			% 	stim_protocols:
			% 						ntj = naive test with juice
			% 						default
			% 	
			% 	stim_type:
			% 						'NaiveTest-No Juice'
			% 						'NaiveTest+Juice'
			% 
			% 	ChR2data: 			for all except ntj, this is the post-processed legacy ChR2 structure OR its path
			% 						But for 'NaiveTest+Juice', this is the folder with the CED file and exclusions
            %   MELIE USE:
			%       obj = CLASS_HSOM_optogenetics_legacy_obj('batch', 'A', 'A', [], [], 'emery')
			%	obj = CLASS_HSOM_optogenetics_legacy_obj('collate', 'A', 'A', [], [], 'bootAUC')
			% 
			% 	obj = CLASS_HSOM_optogenetics_legacy_obj('collate', 'nt', 'NaiveTest-No Juice', [], [], 'ntRaster')
			% 	obj = CLASS_HSOM_optogenetics_legacy_obj('collate', 'ntj', 'NaiveTest+Juice', [], [], 'ntjRaster')
			% 	obj = CLASS_HSOM_optogenetics_legacy_obj('batch', 'ntj', 'NaiveTest+Juice', [], [], 'ntjRaster')
			% 
			obj.iv.runID = randi(10000);
			obj.iv.versionCode = ['CLASS_HSOM_optogenetics_legacy_obj v2.0 Modified 10-29-19 18:36 | obj created: ' datestr(now)];
			

			if nargin < 1 || isempty(Mode)
				obj.iv.Mode = 'default';
			else
				obj.iv.Mode = Mode;
			end
			if nargin < 2  || isempty(stim_protocol)
				obj.iv.stim_protocol = 'default';
			else
				obj.iv.stim_protocol = stim_protocol;
			end
			if nargin < 3 || isempty(stimtype)
				% obj.iv.stimtype = 'NaiveTest-No Juice';
				obj.iv.stimtype = 'A';
			else
				obj.iv.stimtype = stimtype;
			end
			if nargin < 4
				ChR2data = [];
			end
			if nargin < 5
				Xdata = [];
			end
			if nargin < 6
				obj.iv.collateKey = '';
			else
				obj.iv.collateKey = collateKey;
			end

			obj.getDataset(ChR2data, Xdata);
			
			if strcmp(obj.iv.Mode, 'batch') 
				disp('Batch processing is complete! The shell object is not saved.')
				alert = ['Optogenetics Legacy Batch Obj Complete: ' num2str(obj.iv.runID)]; 
			    mailAlert(alert);
		    elseif strcmp(obj.iv.Mode, 'collate')
		    	% 
		    	% 	Save the collated results
		    	% 
		    	obj.save;
		    	disp('Batch analysis is complete! The shell object is saved to Host Folder.')
				alert = ['Optogenetics Legacy Analysis Obj Complete: ' num2str(obj.iv.runID)]; 
			    mailAlert(alert);
		    else
				obj.getParsibleLicks;
				obj.resetBootResults;
				obj.save;
			end


		end
		function getDataset(obj, ChR2data, Xdata)
			% 
			% 	INIT >> getDataset
			% 
			if strcmpi(obj.iv.Mode, 'batch')
				if strcmpi(obj.iv.stim_protocol, 'ntj')
		    		% 
		    		% 	NAIVE TEST WITH JUICE MUST BE HANDLED DIFFERENTLY! We need to initialize in a whole new way that gathers variables from CED
		    		% 
		    		instructions = sprintf('Batch Initialization of Optogentics NAIVE TEST PLUS JUICE Objects \n\n **ONLY USE THIS FOR STUFF WITH CED-FILES FOR NTJ FILES -- Won''t work with timing bc CED fields are different \n\n 1. Set Up HOSTObj folder > NTJ > then folders within this for each day. Each folder is named MouseID_ntj_Day#. Put CED file and exclusions.txt in here. \n 2. Select the NTJ Host Folder \n 3. Each folder will be processed and OptoObj saved to the day''s folder.')
		    		hhh = msgbox(instructions);
			    	
	                hostFolder = uigetdir('', 'Select the HOSTObj > NTJ Host Folder for this batch');
	                cd(hostFolder)
	                if exist('hhh', 'var')
		                close(hhh);
	                end
                	disp('====================================================')
				    disp('			Batch Optogentics NTJ Processing 	 	  ')
				    disp('====================================================')
				    disp(' ')
				    disp(['Started: ' datestr(now)])
				    disp(' ')
	                hostFiles = dir(hostFolder);
					dirFlags = [hostFiles.isdir];
					subFolders = hostFiles(dirFlags);
					folderNames = {subFolders(3:end).name};
					folderPaths = {subFolders(3:end).folder};
					obj.iv.files = folderNames;
					disp(char(['Loading the following files...' folderNames]))
					disp(' ')
			    	disp('-----------------------------------------------------------------')
			

					for iset = 1:numel(folderNames)
						fprintf(['Working on file #' num2str(iset) ': ' folderNames{iset} '(' num2str(iset) '/' num2str(numel(folderNames)) ' ' datestr(now,'HH:MM AM') ') \n'])
                		cd(folderNames{iset})
                		% 
						% 	Check what info is available to us in the subfolder. If we want a box200 gfit, we need to load the gfit. If exclusions are present we will add them
						% 
						dirFiles = dir;
						% 
						% 	First, check if a statObj is already present:
						% 
						sObjpos = find(contains({dirFiles.name},'OptoLegacyObj_'));
						
						if ~isempty(sObjpos) && numel(sObjpos) < 2
							disp('		Detected OptoLegacyObj in folder--skip!')
							% sObj = load(obj.correctPathOS([dirFiles(sObjpos).folder, '\' dirFiles(sObjpos).name]));
		     %                sObjfield = fieldnames(sObj);
		     %                eval(['sObj = sObj.' sObjfield{1} ';']);
						else
							% 
							% 	NEED TO LOAD CED DATA AND PROCESS AHEAD OF TIME...do this in default option ABOVE
							% 
							CLASS_HSOM_optogenetics_legacy_obj('default', obj.iv.stim_protocol, obj.iv.stimtype);
						end	
						disp(' ')
						disp('	File saved and complete.')
			    		disp('-----------------------------------------------------------------')
                        cd('..')
					end
		    	else
		    		instructions = sprintf('Batch Initialization of Optogentics Legacy Objects \n\n **ONLY USE THIS FOR STUFF PREPROCESSED WITH roadmap_v1_4_init.m -- Won''t work with naive test + juice because CED fields are different \n\n 1. Select the ChR2 and X legacy files \n 2. Press cancel when finished selecting. You can do multiple files at once. \n 3. The OptoObjs will be saved in the folder where the legacy .mat files are saved.')
		    		hhh = msgbox(instructions);
			    	go = true;
			    	ifile = 1;
			    	while go
	                    [file, path] = uigetfile('*.mat',['Select ChR2 File #' num2str(ifile)], 'MultiSelect','on');
	                    if length(file) > 1
	                        file_paths(ifile).dir = path;
	                        file_paths(ifile).ChR2 = [path, file{1}];
	                        file_paths(ifile).X = [path, file{2}];
	                        cd(path);
	                    else
	                        if file == 0
	                            fprintf(['Cancelled selection of file # ' num2str(ifile) '\n\n']);
	                            break
	                        end
	                        file_paths(ifile).ChR2 = [path, file];
	                        cd(path);
	                        [file, path] = uigetfile('*.mat',['Select X File #' num2str(ifile)]);
	                        if file == 0
	                            fprintf(['Cancelled selection of file #' num2str(ifile) ' \n Deleting ChR2 data... \n\n']);
	                            file_paths(ifile).dir = path;
	                            file_paths(ifile).ChR2 = [];
	                            file_paths(ifile).X = [];
	                        end
	                        file_paths(ifile).X = [path, file];
	                        fprintf([log_data, 'X #', num2str(ifile), '= ' file ' \n']);
	                    end
	                    ifile = ifile + 1;
	                end
	                if exist('hhh', 'var')
		                close(hhh);
	                end

				    returndir = pwd;

				    disp('====================================================')
				    disp('			Batch Optogentics Legacy Processing 	  ')
				    disp('====================================================')
				    disp(' ')
				    disp(['Started: ' datestr(now)])
				    disp(' ')
				    disp('Loading the following files...')
				    for ufile = 1:ifile-1
				    	disp(['	#' num2str(ufile)	': ' file_paths(ufile).ChR2])
			    	end
			    	disp(' ')
			    	disp('-----------------------------------------------------------------')
				    for ufile = 1:ifile-1
				    	cd(file_paths(ufile).dir)
				    	disp(['Working on file #' num2str(ufile) ': ' file_paths(ufile).ChR2])
				    	disp(['	' '(' num2str(ufile) '/' num2str(numel(file_paths)) ' ' datestr(now) ')...'])
				    	CLASS_HSOM_optogenetics_legacy_obj('default', obj.iv.stim_protocol, obj.iv.stimtype, file_paths(ufile).ChR2, file_paths(ufile).X);
				    	disp('	File saved and complete.')
				    	disp('-----------------------------------------------------------------')
			    	end
			    	cd(returndir);
		    	
	    		end
    		elseif strcmpi(obj.iv.Mode, 'collate')
    			obj.runCollate;
			else
				% 
				% 	SINGLE NTJ
				% 
				if strcmpi(obj.iv.stim_protocol, 'ntj')
					dirFiles = dir;
					excPos = find(contains({dirFiles.name},'exclusions'));
					if ~isempty(excPos)
						obj.iv.excludedtrials_  = fileread(obj.correctPathOS([dirFiles(excPos).folder, '\', dirFiles(excPos).name]));
					else
						error('Need to put in exclusion file for autoloader!')
					end
					s7spos = find(~contains({dirFiles(3:end).name},'excl') & ~contains({dirFiles(3:end).name},'Obj') & ~contains({dirFiles(3:end).name},'obj')) + 2;
					s7s = load(obj.correctPathOS([dirFiles(s7spos).folder, '\', dirFiles(s7spos).name]));
					signals = fieldnames(s7s);
					% 
					% 	Once in the spike 2 file, extract anything relevant for actual signals...
					% 
					fieldIdx = contains(signals, 'StimON');
					eval(cell2mat(['obj.time.s.stimOn_s = s7s.' signals(fieldIdx) '.times;']));
					fieldIdx = contains(signals, 'Lick');
					eval(cell2mat(['obj.time.s.lick_s = s7s.' signals(fieldIdx) '.times;']));
					fieldIdx = contains(signals, 'StimOFF');
					eval(cell2mat(['obj.time.s.stimOff_s = s7s.' signals(fieldIdx) '.times;']));
					fieldIdx = contains(signals, 'RefrxOFF');
					eval(cell2mat(['obj.time.s.refrxOff_s = s7s.' signals(fieldIdx) '.times;']));
					fieldIdx = contains(signals, 'Juice');
					eval(cell2mat(['obj.time.s.juice_s = s7s.' signals(fieldIdx) '.times;']));
					fieldIdx = contains(signals, 'ChR2');
					eval(cell2mat(['obj.h.ChR2values = s7s.' signals(fieldIdx) '.values;']));
					eval(cell2mat(['obj.h.ChR2times = s7s.' signals(fieldIdx) '.times;']));
					fieldIdx = contains(signals, 'X');
					eval(cell2mat(['obj.m.X = s7s.' signals(fieldIdx) '.values;']));
					eval(cell2mat(['obj.m.Xtimes = s7s.' signals(fieldIdx) '.times;']));
					obj.m.X = [0;abs(obj.m.X(2:end)-obj.m.X(1:end-1))];
					obj.m.XgfitMode = 'Abs-Xderivative';
					fieldIdx = contains(signals, 'EMG');
					eval(cell2mat(['obj.m.EMG = s7s.' signals(fieldIdx) '.values;']));
					eval(cell2mat(['obj.m.EMGtimes = s7s.' signals(fieldIdx) '.times;']));
					obj.m.EMG = abs(obj.m.EMG);
					obj.m.EMGgfitMode = 'Abs-EMG';
					fieldIdx = contains(signals, 'IRtrig');
					eval(cell2mat(['obj.m.IRtrig = s7s.' signals(fieldIdx) '.times;']));
					fieldIdx = contains(signals, 'CamO');
                    if numel(fieldIdx) > 1
                        fieldIdx = endsWith(signals, 'CamO');
                        eval(cell2mat(['obj.m.CamOtimes = s7s.' signals(fieldIdx) '.times;']));
                        fieldIdx = endsWith(signals, 'CamO2');
                        if sum(fieldIdx)>0
                            eval(cell2mat(['obj.m.CamO2times = s7s.' signals(fieldIdx) '.times;']));
                        end
                    else
                        eval(cell2mat(['obj.m.CamOtimes = s7s.' signals(fieldIdx) '.times;']));
                    end
					
					obj.iv.init_variables_ChR2.num_trials = numel(obj.time.s.stimOn_s);
                    obj.iv.num_trials = obj.iv.init_variables_ChR2.num_trials;
		            obj.parceExclusions;
					
					

				    obj.iv.samples_per_ms_ChR2 = mode(obj.h.ChR2times(2:end)-obj.h.ChR2times(1:end-1))*1000;
				    obj.iv.samples_per_ms_X = mode(obj.m.Xtimes(2:end)-obj.m.Xtimes(1:end-1))*1000;
				    obj.iv.samples_per_ms_EMG = mode(obj.m.EMGtimes(2:end)-obj.m.EMGtimes(1:end-1))*1000;
				    obj.iv.samples_per_ms_CamO = mode(obj.m.CamOtimes(2:end)-obj.m.CamOtimes(1:end-1))*1000;
				    % 
				    % 	Set up time field
				    % 
					% 	- note that zeros indicate a timepoint excluded for a rxn (within 700ms of cue)
					% 	- note that NaNs indicate this was not a stim or unstim trial, depending on the case
					% 
				    obj.time.Plot.stim_on_time_s 	= 0;
				    obj.time.Plot.stim_off_time_s 	= obj.time.s.stimOff_s(1) - obj.time.s.stimOn_s(1);
				    obj.time.Plot.refrx_time_s 		= obj.time.s.refrxOff_s(1) - obj.time.s.stimOn_s(1);
				    obj.time.Plot.max_time_s		= max(obj.time.s.stimOn_s(2:end)-obj.time.s.stimOn_s(1:end-1));

					obj.getBinnedLicks('stimOn', 30000, 30000)
					obj.lick.swrtStimOn.lick_s = obj.time.s.binnedLicks.lick_s;
					obj.lick.swrtStimOn.lick_ex_s = obj.time.s.binnedLicks.lick_ex_s;
					obj.lick.swrtStimOn.f_lick_ex_s = obj.time.s.binnedLicks.f_lick_ex_s_wrtref;
					

					obj.getBinnedLicks('juice', 30000, 30000)
					obj.lick.swrtJuice.lick_s = obj.time.s.binnedLicks.lick_s;
					obj.lick.swrtJuice.lick_ex_s = obj.time.s.binnedLicks.lick_ex_s;
					obj.lick.swrtJuice.f_lick_ex_s = obj.time.s.binnedLicks.f_lick_ex_s_wrtref;	

					% 
					% 	From now on, consideration bound is a non-static property of obj.boot
					% 
					obj.bootStat.notes{1} = 'Consideration boundary is non-static property changed in obj.boot fxn';
					obj.setConsiderationBound(0);		
				% 
				% 	SINGLE Standard OptoObj
				% 
		    	else
		    		if isempty(ChR2data)
	                    [file, path] = uigetfile('*.mat','Select ChR2 File', 'MultiSelect','on');
	                    if length(file) > 1
	                        file_paths(1).ChR2 = [path, file{1}];
	                        file_paths(1).X = [path, file{2}];
	                        cd(path);
	                    else
	                        if file == 0
	                            fprintf('Cancelled selection of file \n\n');
	                            error('No ChR2 file')
	                        end
	                        file_paths(1).ChR2 = [path, file];
	                        cd(path);
	                        [file, path] = uigetfile('*.mat','Select X File');
	                        if file == 0
	                            fprintf('Cancelled selection of file \n Deleting ChR2 data... \n\n');
	                            file_paths(1).ChR2 = [];
	                            file_paths(1).X = [];
	                        end
	                        file_paths(1).X = [path, file];
	                        fprintf(['X #', num2str(1), '= ' file ' \n']);
	                    end
	                end

				    if isempty(ChR2data)
					    optogenetics_data_struct = load(obj.correctPathOS(file_paths(1).ChR2));
		                F = fieldnames(optogenetics_data_struct);
		                optogenetics_data_struct = optogenetics_data_struct.(F{1});
				    elseif isnumeric(ChR2data)
				    	optogenetics_data_struct = ChR2data;
	                else
			    		optogenetics_data_struct = load(obj.correctPathOS(ChR2data));
		                F = fieldnames(optogenetics_data_struct);
		                optogenetics_data_struct = optogenetics_data_struct.(F{1});
			    	end
				    

				    init_variables_ChR2 = optogenetics_data_struct.init_variables;
				    h = optogenetics_data_struct.stim_struct;
				    lick_data_struct = optogenetics_data_struct.lick_data_struct;
				    samples_per_ms_ChR2 = init_variables_ChR2.time_parameters.samples_per_ms;
				    obj.iv.ExcludedTrials = optogenetics_data_struct.exclusions_struct.Excluded_Trials;
				    % 
				    % 	Set up time field
				    % 
					% 	- note that zeros indicate a timepoint excluded for a rxn (within 700ms of cue)
					% 	- note that NaNs indicate this was not a stim or unstim trial, depending on the case
					% 
				    obj.time.ms.cue_on_time_ms 	= init_variables_ChR2.time_parameters.ms.cue_on_time_ms;
				    obj.time.ms.rxn_time_ms 	= init_variables_ChR2.time_parameters.ms.rxn_time_ms;
				    obj.time.ms.rxn_ok_ms 		= init_variables_ChR2.time_parameters.ms.rxn_ok_ms;
				    obj.time.ms.buffer_ms 		= init_variables_ChR2.time_parameters.ms.buffer_ms;
				    obj.time.ms.target_time_ms	= init_variables_ChR2.time_parameters.ms.target_time_ms;
				    obj.time.ms.op_rew_open_ms	= init_variables_ChR2.time_parameters.ms.op_rew_open_ms;
				    obj.time.ms.ITI_time_ms 	= init_variables_ChR2.time_parameters.ms.ITI_time_ms;
				    obj.time.ms.total_time_ms 	= init_variables_ChR2.time_parameters.ms.total_time_ms;

				    obj.time.s.cue_on_time_s 	= init_variables_ChR2.time_parameters.ms.cue_on_time_ms/1000;
				    obj.time.s.rxn_time_s 		= init_variables_ChR2.time_parameters.ms.rxn_time_ms/1000;
				    obj.time.s.rxn_ok_s 		= init_variables_ChR2.time_parameters.ms.rxn_ok_ms/1000;
				    obj.time.s.buffer_s 		= init_variables_ChR2.time_parameters.ms.buffer_ms/1000;
				    obj.time.s.target_time_s	= init_variables_ChR2.time_parameters.ms.target_time_ms/1000;
				    obj.time.s.op_rew_open_s	= init_variables_ChR2.time_parameters.ms.op_rew_open_ms/1000;
				    obj.time.s.ITI_time_s 		= init_variables_ChR2.time_parameters.ms.ITI_time_ms/1000;
				    obj.time.s.total_time_s 	= init_variables_ChR2.time_parameters.ms.total_time_ms/1000;
				    
				    obj.time.pos_ChR2.cue_on_time_pos_ChR2 	= init_variables_ChR2.time_parameters.ms.cue_on_time_ms*samples_per_ms_ChR2;
				    obj.time.pos_ChR2.rxn_time_pos_ChR2		= init_variables_ChR2.time_parameters.ms.rxn_time_ms*samples_per_ms_ChR2;
				    obj.time.pos_ChR2.rxn_ok_pos_ChR2		= init_variables_ChR2.time_parameters.ms.rxn_ok_ms*samples_per_ms_ChR2;
				    obj.time.pos_ChR2.buffer_pos_ChR2		= init_variables_ChR2.time_parameters.ms.buffer_ms*samples_per_ms_ChR2;
				    obj.time.pos_ChR2.target_time_pos_ChR2	= init_variables_ChR2.time_parameters.ms.target_time_ms*samples_per_ms_ChR2;
				    obj.time.pos_ChR2.op_rew_open_pos_ChR2	= init_variables_ChR2.time_parameters.ms.op_rew_open_ms*samples_per_ms_ChR2;
				    obj.time.pos_ChR2.ITI_time_pos_ChR2 	= init_variables_ChR2.time_parameters.ms.ITI_time_ms*samples_per_ms_ChR2;
				    obj.time.pos_ChR2.total_time_pos_ChR2 	= init_variables_ChR2.time_parameters.ms.total_time_ms*samples_per_ms_ChR2;
				    % 
				    % 	Set up the legacy lick field -- we should never use this in analysis, only for debugging
				    %
				    obj.lick.legacy.s.note = {'Only use this field for debugging, not analysis'}; 
				    obj.lick.legacy.s.lick_times_by_trial_s 		= lick_data_struct.lick_times_by_trial;
				    obj.lick.legacy.s.f_lick_rxn_s 					= lick_data_struct.f_lick_rxn;
				    obj.lick.legacy.s.f_lick_operant_no_rew_s 		= lick_data_struct.f_lick_operant_no_rew;
				    obj.lick.legacy.s.f_lick_operant_rew_s 			= lick_data_struct.f_lick_operant_rew;
				    obj.lick.legacy.s.f_lick_ITI_s 					= lick_data_struct.f_lick_ITI;
				    obj.lick.legacy.s.all_first_licks_s 			= lick_data_struct.all_first_licks;
				    obj.lick.legacy.s.lick_ex_times_by_trial_s 		= lick_data_struct.lick_ex_times_by_trial;
				    obj.lick.legacy.s.f_ex_lick_rxn_s 				= lick_data_struct.f_ex_lick_rxn;
				    obj.lick.legacy.s.f_ex_lick_operant_no_rew_s 	= lick_data_struct.f_ex_lick_operant_no_rew;
				    obj.lick.legacy.s.f_ex_lick_operant_rew_s 		= lick_data_struct.f_ex_lick_operant_rew;
				    obj.lick.legacy.s.f_ex_lick_ITI_s 				= lick_data_struct.f_ex_lick_ITI;
				    obj.lick.legacy.s.all_ex_first_licks_s 			= lick_data_struct.all_ex_first_licks;
					% 
					% 	Now handle the non-legacy lick variables
					% 
					obj.iv.init_variables_ChR2 = init_variables_ChR2;
					obj.getParsibleLicks
					obj.getBinnedLicks
					ll = obj.time.s.binnedLicks.f_lick_s_wrtref;
					ll(isnan(ll)) = 0;
					ll(obj.iv.ExcludedTrials) = 0;
					obj.lick.swrtc.stim_first_licks_swrtc = nan(size(ll));
					obj.lick.swrtc.stim_first_licks_swrtc(h.stimTrials) = ll(h.stimTrials);
					obj.lick.swrtc.unstim_first_licks_swrtc = nan(size(ll));
					obj.lick.swrtc.unstim_first_licks_swrtc(h.nostim_trials) = ll(h.nostim_trials);
					% 
					obj.lick.swrtc.lick_ex_times_by_trial_swrtc = obj.time.s.binnedLicks.lick_s;
	                for iExcl = obj.iv.ExcludedTrials
	    				obj.lick.swrtc.lick_ex_times_by_trial_swrtc{iExcl} = [];
	                end


				    if isempty(Xdata)
					    movement_data_struct = load(obj.correctPathOS(file_paths(1).X));
					    F = fieldnames(movement_data_struct);
						movement_data_struct = movement_data_struct.(F{1});
				    elseif isnumeric(Xdata)
				    	movement_data_struct = Xdata;
			    	else
			    		movement_data_struct = load(obj.correctPathOS(Xdata));
					    F = fieldnames(movement_data_struct);
						movement_data_struct = movement_data_struct.(F{1});
			    	end
				    init_variables_X = movement_data_struct.init_variables;	
				    samples_per_ms_X = init_variables_X.time_parameters.samples_per_ms;
				    obj.time.pos_X.cue_on_time_posX 	= init_variables_ChR2.time_parameters.ms.cue_on_time_ms*samples_per_ms_ChR2;
				    obj.time.pos_X.rxn_time_posX 		= init_variables_ChR2.time_parameters.ms.rxn_time_ms*samples_per_ms_ChR2;
				    obj.time.pos_X.rxn_ok_posX 			= init_variables_ChR2.time_parameters.ms.rxn_ok_ms*samples_per_ms_ChR2;
				    obj.time.pos_X.buffer_posX 			= init_variables_ChR2.time_parameters.ms.buffer_ms*samples_per_ms_ChR2;
				    obj.time.pos_X.target_time_posX		= init_variables_ChR2.time_parameters.ms.target_time_ms*samples_per_ms_ChR2;
				    obj.time.pos_X.op_rew_open_posX		= init_variables_ChR2.time_parameters.ms.op_rew_open_ms*samples_per_ms_ChR2;
				    obj.time.pos_X.ITI_time_posX 		= init_variables_ChR2.time_parameters.ms.ITI_time_ms*samples_per_ms_ChR2;
				    obj.time.pos_X.total_time_posX 		= init_variables_ChR2.time_parameters.ms.total_time_ms*samples_per_ms_ChR2;
				    %
				    %	Initialize the movement structure
				    %
				    m = {};
				    %
				    % 
				    for itrial = 1:init_variables_X.num_trials
				        %
				        %	Get smoothed movement data
				        %
				        m.smooth_vbt(itrial, :) = gausssmooth(movement_data_struct.signal_ex_values_by_trial(itrial,:), 50, 'gauss');
				        %
				        %	Get low pass movement data - 1 sec
				        %
				        m.lowPass_vbt(itrial, :) = gausssmooth(movement_data_struct.signal_ex_values_by_trial(itrial,:), 500, 'gauss');
				        % 
				        % 	Get high (band) pass movement data - 1 sec
				        % 
				        m.bandPass_vbt(itrial, :) = gausssmooth(movement_data_struct.signal_ex_values_by_trial(itrial,:) - m.lowPass_vbt(itrial, :), 50, 'gauss');
				    end		
					% 
					% 	Stim position and in seconds wrt cue
					% 
					h.stim_on_position = h.stim_on_positions{1};
					h.stim_on_s_wrt_cue = h.stim_on_position/samples_per_ms_ChR2/1000;
					% 
					% 	From now on, consideration bound is a non-static property of obj.boot
					% 
					obj.bootStat.notes{1} = 'Consideration boundary is non-static property changed in obj.boot fxn';
					obj.setConsiderationBound(0);		
				    % 
				    % 	Gather variables
				    % 
		            if exist('file_paths', 'var')
		    		    obj.iv.file_paths = file_paths;
		            end
				    obj.h = h;
				    obj.m = m;
				    obj.iv.init_variables_ChR2 = init_variables_ChR2;
				    obj.iv.samples_per_ms_ChR2 = samples_per_ms_ChR2;
				    obj.iv.init_variables_X = init_variables_X;
				    obj.iv.samples_per_ms_X = samples_per_ms_X;
				    obj.lick_data_struct = lick_data_struct;
			    end		    
			end
		end
		function runCollate(obj)
			% 
			% 	INIT >> runCollate
			% 
    		% 	Collating results for objects in HOSTObj folder
    		% 
    		instructions = sprintf('Collating Results: \n\n 1. Set Up HOSTObj folder > stimType > then folders within this for each day. Each folder is named MouseID_stimType_Day#. \n 2. Select the Host Folder \n 3. Each folder will be processed -- existing OptoObjs will be opened and then resaved to the day''s folder. \n 4. Results of the collation will be saved in the Obj')
    		hhh = msgbox(instructions);
	    	
            hostFolder = uigetdir('', 'Select the HOSTObj > Stim Host Folder for this collation batch');
            cd(hostFolder)
            if exist('hhh', 'var')
                close(hhh);
            end
        	disp('====================================================')
		    disp(['			Collating Results - ' obj.iv.collateKey ' 	 	  '])
		    disp('====================================================')
		    disp(' ')
		    disp(['Started: ' datestr(now)])
		    disp(' ')
            hostFiles = dir(hostFolder);
			dirFlags = [hostFiles.isdir];
			subFolders = hostFiles(dirFlags);
			folderNames = {subFolders(3:end).name};
			folderPaths = {subFolders(3:end).folder};
			obj.iv.files = folderNames;
			disp(char(['Loading the following files...' folderNames]))
			disp(' ')
	    	disp('-----------------------------------------------------------------')

	    	if strcmpi(obj.iv.collateKey, 'ecdf') || strcmpi(obj.iv.collateKey, 'bootAUC') || strcmpi(obj.iv.collateKey, 'emery')
	    		disp('Select fodler to save figures')
	    		obj.iv.figuresDir = uigetdir;
    		end
	

			for iset = 1:numel(folderNames)
				fprintf(['Working on file #' num2str(iset) ': ' folderNames{iset} '(' num2str(iset) '/' num2str(numel(folderNames)) ' ' datestr(now,'HH:MM AM') ') \n'])
        		cd(folderNames{iset})
        		% 
				% 	Check what info is available to us in the subfolder. If we want a box200 gfit, we need to load the gfit. If exclusions are present we will add them
				% 
				dirFiles = dir;
				% 
				% 	First, check if a statObj is already present:
				% 
				sObjpos = find(contains({dirFiles.name},'OptoLegacyObj_'));
				try
					if ~isempty(sObjpos)
						% 
						% 	Find the newest version of the obj
						% 
						idxdates = [dirFiles(sObjpos).datenum];
						newestObj = find(idxdates == max([dirFiles(sObjpos).datenum]));
						sObjpos = sObjpos(newestObj);

						sObj = load(obj.correctPathOS([dirFiles(sObjpos).folder, '\' dirFiles(sObjpos).name]));
	                    sObjfield = fieldnames(sObj);
	                    eval(['sObj = sObj.' sObjfield{1} ';']);
	                    sObj.resetBootResults;
	                    sObj.analyze(obj.iv.collateKey, obj);
	                    obj.bootStat.collateResult(iset).analysisType = obj.iv.collateKey;
	                    obj.bootStat.collateResult(iset).sessionID = folderNames{iset};
	                    if strcmpi(obj.iv.collateKey, 'bootAUC')
	                    	obj.bootStat.collateResult(iset).p = sObj.bootStat.resultLick.p;
	                    	obj.bootStat.collateResult(iset).delAUC_EOT = sObj.bootStat.resultLick.ref;
                    	elseif strcmpi(obj.iv.collateKey, 'Emery')
                    		obj.bootStat.collateResult(iset).median = sObj.bootStat.resultLick(1).median;
                    		obj.bootStat.collateResult(iset).mean = sObj.bootStat.resultLick(1).mean;
                    		obj.bootStat.collateResult(iset).CI95 = sObj.bootStat.resultLick(1).CI95;
	                	end
	                    obj.bootStat.collateResult(iset).bootStat = sObj.bootStat;
	                    sObj.save;
					else
						error('No OptoLegacyObj in the folder!')
					end	
					disp(' ')
					disp('	File analyzed, resaved and complete.')
				catch ex
					EE = getReport(ex)
					disp(EE)
					warning(['Error while processing this file. The message was:' EE])
					disp('	Skipping this file. Add it to the analysis obj later')

                    obj.bootStat.collateResult(iset).analysisType = obj.iv.collateKey;
                    obj.bootStat.collateResult(iset).sessionID = folderNames{iset};
                    obj.bootStat.collateResult(iset).bootStat = ['Error Encountered:' EE];

				end
	    		disp('-----------------------------------------------------------------')
                cd(hostFolder)
			end
		end

		

		%------------------------------------------------------------------------------------------------------------------- 
		% 	SINGLE SESSION VISUALIZATION METHODS
		%-------------------------------------------------------------------------------------------------------------------
		function ax = plot(obj, Mode, Blocking, division, ref, rxnwin_s)
			% 
			% 	Mode options:
			% 		raster4 -- plots a raster of stimulated and unstimulated licks. top row is flick only, bottom is all licks
			% 		ecdf 	-- plots the full ecdf
			% 		hxg 	-- plots histogram and cdf of flicks and all licks (normalized by total # of licks in each category)
			% 		hxg-counts plots histogram without normalizing for # of trials
			% 		X 		-- plots concatenated X for cue-to-lick for all trials where lick occurs before EOT
			% 
			% 	Blocking: select either valid or even
			% 			NB: for NT, Blocking == ax, the axis for overlaying plots!
			% 
			% 	division: the block # to plot for.
			% 			NB: for NT, division == the collated set idxs in bootStat to use in plot. a vector of numbers
			% 
			% 	IF plotting collated NT --> breaks and sends to the proper obj.plotCollatedNT regardless of input
			% 
			if strcmpi(obj.iv.Mode, 'collate') && strcmpi(obj.iv.stim_protocol, 'nt')
				if nargin < 2
					ax = obj.plotCollatedNT;
				elseif nargin < 3
					ax = obj.plotCollatedNT(Mode);
				elseif nargin < 4
					ax = obj.plotCollatedNT(Mode, Blocking);
				else
					ax = obj.plotCollatedNT(Mode, Blocking, division);
				end
				return
			elseif strcmpi(obj.iv.Mode, 'collate') && strcmpi(obj.iv.stim_protocol, 'ntj')
				if nargin < 2
					ax = obj.plotCollatedNTJ;
				elseif nargin < 3
					ax = obj.plotCollatedNTJ(Mode);
				elseif nargin < 4
					ax = obj.plotCollatedNTJ(Mode, Blocking);
				else
					ax = obj.plotCollatedNTJ(Mode, Blocking, division);
				end
				return
			end
			if isempty(obj.valid_blocks) && ~strcmpi(obj.iv.stim_protocol, 'ntj')
				obj.getvalidblocks;
			end
			if nargin < 2
				Mode = 'raster4';
			end
			if nargin < 3
				Blocking = 'valid';
			end
			if nargin < 4
				division = 1;
			end
			if nargin < 5 && ~strcmpi(obj.iv.stim_protocol, 'ntj')
				ref = 'cue';
			end
			if nargin < 6
				rxnwin_s = 0;
			end
			if strcmpi(obj.iv.stim_protocol, 'ntj')
				if nargin < 5
					ref = 'stimOn';
				end
				ax = obj.plotNTJ(Mode, ref, rxnwin_s);
			else				
				obj.getBinnedLicks(ref, 30000, 30000, rxnwin_s);

				i_div = division;
				if i_div ~=1 || ~(strcmpi(obj.valid_blocks.Mode,'fullday') || strcmpi(obj.valid_blocks.Mode,'full-day'))
					error('Not implemented at obj.valid_blocks.all_trials_this_block @getvalidblocks')
				end


				% 
				% 	Collect var pntrs
				% 
				fs_ms = obj.iv.init_variables_ChR2.time_parameters.samples_per_ms;
				rb_ms = obj.time.ms.op_rew_open_ms;
				eot = obj.time.s.ITI_time_s;
				cue = obj.time.ms.cue_on_time_ms;
				considbound = obj.bootStat.considerationBoundary_ms/1000;
				
				f = figure;
				if strcmpi(Mode, 'raster4')
					ax11nostim = subplot(2, 2,  1, 'Parent', f);
					hold(ax11nostim, 'on')
					ax11stim = subplot(2, 2,  2, 'Parent', f);
					hold(ax11stim, 'on')
					ax11ONLYnostim = subplot(2, 2,  3, 'Parent', f);
					hold(ax11ONLYnostim, 'on')
					ax11ONLYstim = subplot(2, 2,  4, 'Parent', f);
					hold(ax11ONLYstim, 'on')
					% legend(ax11nostim, 'show');
					% 
					% 	Plot raster of all licks with first licks overlaid
					% 
					plot(ax11nostim, zeros(size(obj.valid_blocks.all_trials_this_block)), obj.valid_blocks.all_trials_this_block,'r-', 'DisplayName', 'Cue')
				    plot(ax11nostim, [considbound, considbound], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'b-', 'DisplayName', 'Consideration Limit')
					plot(ax11nostim, [rb_ms/1000, rb_ms/1000], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'k-', 'DisplayName', 'Reward Boundary')
					plot(ax11nostim, [eot, eot], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'k-', 'DisplayName', 'ITI Start')
					plot(ax11nostim, [rxnwin_s, rxnwin_s], [1, numel(obj.time.s.binnedLicks.refevents)],'g-', 'DisplayName', 'Permitted Reaction Window')
					set(ax11nostim,  'YDir','reverse')
					title(ax11nostim, ['No Stim Lick Raster Aligned to ' obj.time.s.binnedLicks.ref])
					xlim(ax11nostim, [0, obj.time.ms.total_time_ms/1000])
			        if numel(obj.valid_blocks.all_trials_this_block(end)) > 1
			    		ylim(ax11nostim, [1, obj.valid_blocks.all_trials_this_block(end)])
			        end
					ax11nostim.XLabel.String = ['Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax11nostim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];

					plot(ax11stim, zeros(size(obj.valid_blocks.all_trials_this_block)), obj.valid_blocks.all_trials_this_block,'r-', 'DisplayName', 'Cue')
				    plot(ax11stim, [considbound, considbound], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'b-', 'DisplayName', 'Consideration Limit')
					plot(ax11stim, [rb_ms/1000, rb_ms/1000], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'k-', 'DisplayName', 'Reward Boundary')
					plot(ax11stim, [eot, eot], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'k-', 'DisplayName', 'ITI Start')
					
					set(ax11stim,  'YDir','reverse')
					title(ax11stim, '+ stimulation')
					xlim(ax11stim, [0, obj.time.ms.total_time_ms/1000])
			        if numel(obj.valid_blocks.all_trials_this_block(end)) > 1
			    		ylim(ax11stim, [1, obj.valid_blocks.all_trials_this_block(end)])
			        end
					ax11stim.XLabel.String = ['Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax11stim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];


					plot(ax11ONLYnostim, zeros(size(obj.valid_blocks.all_trials_this_block)), obj.valid_blocks.all_trials_this_block,'r-', 'DisplayName', 'Cue')
				    plot(ax11ONLYnostim, [considbound, considbound], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'b-', 'DisplayName', 'Consideration Limit')
					plot(ax11ONLYnostim, [rb_ms/1000, rb_ms/1000], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'k-', 'DisplayName', 'Reward Boundary')
					plot(ax11ONLYnostim, [eot, eot], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'k-', 'DisplayName', 'ITI Start')
					
					set(ax11ONLYnostim,  'YDir','reverse')
					title(ax11ONLYnostim, 'ONLY Consec Trials, No stim')
					xlim(ax11ONLYnostim, [0, obj.time.ms.total_time_ms/1000])
			        if numel(obj.h.nostim_trials) > 1
			    		ylim(ax11ONLYnostim, [1, numel(obj.h.nostim_trials)])
			        end
					ax11ONLYnostim.XLabel.String = ['Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax11ONLYnostim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];


					plot(ax11ONLYstim, zeros(size(obj.valid_blocks.all_trials_this_block)), obj.valid_blocks.all_trials_this_block,'r-', 'DisplayName', 'Cue')
				    plot(ax11ONLYstim, [considbound, considbound], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'b-', 'DisplayName', 'Consideration Limit')
					plot(ax11ONLYstim, [rb_ms/1000, rb_ms/1000], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'k-', 'DisplayName', 'Reward Boundary')
					plot(ax11ONLYstim, [eot, eot], [obj.valid_blocks.all_trials_this_block(1), obj.valid_blocks.all_trials_this_block(end)],'k-', 'DisplayName', 'ITI Start')
					
					set(ax11ONLYstim,  'YDir','reverse')
					title(ax11ONLYstim, 'ONLY Consec Trials, + stim')
					xlim(ax11ONLYstim, [0, obj.time.ms.total_time_ms/1000])
			        if numel(obj.h.stimTrials) > 1
			    		ylim(ax11ONLYstim, [1, numel(obj.h.stimTrials)])
			        end
					ax11ONLYstim.XLabel.String = ['Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax11ONLYstim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];

				    unstimTrialTracker = 1;
				    stimTrialTracker = 1;
				    all_lick_times_ex_swrtc = obj.lick.swrtc.lick_ex_times_by_trial_swrtc;
				    for iexc = obj.iv.ExcludedTrials
				    	all_lick_times_ex_swrtc{iexc} = [];
			    	end
					plot(ax11stim, obj.valid_blocks.f_lick_s_wrtcue(i_div).stim, obj.valid_blocks.all_trials_this_block, 'mo', 'DisplayName', 'First Lick', 'markerfacecolor', 'm')
					plot(ax11nostim, obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim, obj.valid_blocks.all_trials_this_block, 'mo', 'DisplayName', 'First Lick', 'markerfacecolor', 'm')
					for itrial = obj.valid_blocks.all_trials_this_block(1):obj.valid_blocks.all_trials_this_block(end)
						obj.valid_blocks.plotpnts = all_lick_times_ex_swrtc{itrial};
						if ismember(itrial, obj.h.stimTrials)
							if ~isempty(obj.valid_blocks.plotpnts) && numel(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim)>=itrial
								plot(ax11stim, obj.valid_blocks.plotpnts, itrial.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
								plot(ax11ONLYstim, obj.valid_blocks.plotpnts, stimTrialTracker.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
								plot(ax11ONLYstim, obj.valid_blocks.f_lick_s_wrtcue(i_div).stim(itrial), stimTrialTracker, 'mo', 'DisplayName', 'First Lick')
							end
							stimTrialTracker = stimTrialTracker +1;
						elseif ismember(itrial, obj.h.nostim_trials) && numel(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim)>=itrial
							if ~isempty(obj.valid_blocks.plotpnts)
								plot(ax11nostim, obj.valid_blocks.plotpnts, itrial.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
								plot(ax11ONLYnostim, obj.valid_blocks.plotpnts, unstimTrialTracker.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
								plot(ax11ONLYnostim, obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim(itrial), unstimTrialTracker, 'mo', 'DisplayName', 'First Lick')
							end
							unstimTrialTracker = unstimTrialTracker +1;
						else
							error('shouldn''t get here');
						end					
					end
					yy = get(ax11ONLYnostim, 'ylim');
					ylim(ax11ONLYstim, yy)
                    xx = get(ax11nostim, 'xlim');
                    plot(ax11nostim, xx, [numel(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim)+0.5, numel(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim)+0.5], 'r-', 'displayName', 'Consecutive Trials limit for AUC')
                    plot(ax11stim, xx, [numel(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim)+0.5, numel(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim)+0.5], 'r-', 'displayName', 'Consecutive Trials limit for AUC')
                    
				elseif strcmpi(Mode, 'ecdf')
					ax2 = subplot(1, 2,  1, 'Parent', f);
					hold(ax2, 'on')
					ax3 = subplot(1, 2,  2, 'Parent', f);
					hold(ax3, 'on')
					legend(ax3, 'show')

					plot(ax2, obj.valid_blocks.eCDF(i_div).unstim_x, obj.valid_blocks.eCDF(i_div).unstim_f, 'LineWidth', 3, 'DisplayName', 'no stim')
					hold(ax2, 'on')
					plot(ax2, obj.valid_blocks.eCDF(i_div).stim_x, obj.valid_blocks.eCDF(i_div).stim_f, 'LineWidth', 3, 'DisplayName', '+ stim')

					plot(ax2, [obj.valid_blocks.Median(i_div).unstim, obj.valid_blocks.Median(i_div).unstim], [0,1],  'b-', 'LineWidth', 2, 'DisplayName', 'median no stim')
					plot(ax2, [obj.valid_blocks.Median(i_div).stim, obj.valid_blocks.Median(i_div).stim], [0,1], 'r-', 'LineWidth', 2, 'DisplayName', 'median + stim')

					title(ax2, ['eCDF of First Licks wrt ' obj.time.s.binnedLicks.ref])
					xlim(ax2, [0, obj.time.ms.total_time_ms/1000])
					ylim(ax2, [0,1])
					ax2.XLabel.String = ['Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax2.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];
					% 
					% 	EOT AUC plot (eg up to 7s)
					% 
					plot(ax3, obj.valid_blocks.eCDF(i_div).unstim_x_EOT, obj.valid_blocks.eCDF(i_div).unstim_f_EOT, 'LineWidth', 3, 'DisplayName', 'no stim')
					hold(ax3, 'on')
					plot(ax3, obj.valid_blocks.eCDF(i_div).stim_x_EOT, obj.valid_blocks.eCDF(i_div).stim_f_EOT, 'LineWidth', 3, 'DisplayName', '+ stim')

					plot(ax3, [obj.valid_blocks.Median(i_div).unstim, obj.valid_blocks.Median(i_div).unstim], [0,1],  'b-', 'LineWidth', 2, 'DisplayName', 'median no stim')
					plot(ax3, [obj.valid_blocks.Median(i_div).stim, obj.valid_blocks.Median(i_div).stim], [0,1], 'r-', 'LineWidth', 2, 'DisplayName', 'median + stim')

					title(ax3, ['Trials #', num2str(obj.valid_blocks.all_trials_this_block(1)), '-', num2str(obj.valid_blocks.all_trials_this_block(end)) ' - EOT'])
					xlim(ax3, [0, obj.time.ms.total_time_ms/1000])
					ylim(ax3, [0,1])
					ax3.XLabel.String = ['Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax3.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];


				elseif strcmpi(Mode, 'hxg') || strcmpi(Mode, 'hxg-counts')
					ax10 = subplot(1, 3,  2, 'Parent', f);
					hold(ax10, 'on')
					ax103 = subplot(1, 3,  1, 'Parent', f);
					hold(ax103, 'on')
					ax102 = subplot(1, 3,  3, 'Parent', f);
					hold(ax102, 'on')

					if strcmpi(Mode, 'hxg')
						histogram(ax10, obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim(~isnan(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim)), 30, 'DisplayName', 'No Stimulation', 'Normalization', 'probability')
						histogram(ax10, obj.valid_blocks.f_lick_s_wrtcue(i_div).stim(~isnan(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim)), 30,'DisplayName', '+ Stimulation', 'Normalization', 'probability')
						ax10.YLabel.String = 'Percentage of First Licks Across Category';
					else
						histogram(ax10, obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim(~isnan(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim)), 30, 'DisplayName', 'No Stimulation')
						histogram(ax10, obj.valid_blocks.f_lick_s_wrtcue(i_div).stim(~isnan(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim)), 30,'DisplayName', '+ Stimulation')
						ax10.YLabel.String = '# of First Licks In Category';
					end
					yy = get(ax10, 'ylim');
					plot(ax10, [0, 0], [0,obj.valid_blocks.num_trials_per_block],'r-', 'DisplayName', 'Cue')
				    plot(ax10, [considbound, considbound], [0,obj.valid_blocks.num_trials_per_block],'b-', 'DisplayName', 'Consideration Limit')
					plot(ax10, [rb_ms/1000, rb_ms/1000], [0,obj.valid_blocks.num_trials_per_block],'k-', 'DisplayName', 'Reward Boundary')
					plot(ax10, [obj.time.s.ITI_time_s, obj.time.s.ITI_time_s], [0,obj.valid_blocks.num_trials_per_block],'k-', 'DisplayName', 'ITI Start')
					plot(ax10, [rxnwin_s, rxnwin_s], [0, numel(obj.time.s.binnedLicks.refevents)],'g-', 'DisplayName', 'Permitted Reaction Window')
					ylim(ax10, yy)
					title(ax10, ['Trials #', num2str(obj.valid_blocks.all_trials_this_block(1)), '-', num2str(obj.valid_blocks.all_trials_this_block(end))])
					legend(ax10, 'show')
					xlim(ax10, [-.4, obj.time.ms.total_time_ms/1000])
					ax10.XLabel.String = 'First Lick Time (s wrt cue)';
					


					histogram(ax103, obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim(~isnan(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim)), 30, 'DisplayName', 'No Stimulation', 'Normalization', 'cdf')
					histogram(ax103, obj.valid_blocks.f_lick_s_wrtcue(i_div).stim(~isnan(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim)), 30,'DisplayName', '+ Stimulation', 'Normalization', 'cdf')
					yy = get(ax103, 'ylim');
					plot(ax103, [0, 0], [0,obj.valid_blocks.num_trials_per_block],'r-', 'DisplayName', 'Cue')
				    plot(ax103, [considbound, considbound], [0,obj.valid_blocks.num_trials_per_block],'b-', 'DisplayName', 'Consideration Limit')
					plot(ax103, [rb_ms/1000, rb_ms/1000], [0,obj.valid_blocks.num_trials_per_block],'k-', 'DisplayName', 'Reward Boundary')
					plot(ax103, [obj.time.s.ITI_time_s, obj.time.s.ITI_time_s], [0,obj.valid_blocks.num_trials_per_block],'k-', 'DisplayName', 'ITI Start')
					plot(ax103, [rxnwin_s, rxnwin_s], [0, numel(obj.time.s.binnedLicks.refevents)],'g-', 'DisplayName', 'Permitted Reaction Window')
					ylim(ax103, yy)
					title(ax103, ['Histogram of first licks Aligned to ' obj.time.s.binnedLicks.ref])
					xlim(ax103, [-.4, obj.time.ms.total_time_ms/1000])
					ax103.XLabel.String = ['First Lick Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax103.YLabel.String = 'CDF of First Licks Across Category';


					warning('Hxg of all licks does not incorporate exclusions. Also haven''t divvied up by blocks yet')
					alllicktimes_stimwrtc = cell2mat(obj.time.s.binnedLicks.lick_s(obj.h.stimTrials)');
					alllicktimes_unstimwrtc = cell2mat(obj.time.s.binnedLicks.lick_s(obj.h.nostim_trials)');
					histogram(ax102, alllicktimes_unstimwrtc, 30000, 'DisplayName', 'No Stimulation', 'Normalization', 'probability')
				    hold(ax102, 'on');
					histogram(ax102, alllicktimes_stimwrtc, 30000, 'DisplayName', '+ Stimulation', 'Normalization', 'probability')
					yy = get(ax102, 'ylim');
					plot(ax102, [0, 0], [0,obj.valid_blocks.num_trials_per_block],'r-', 'DisplayName', 'Cue')
				    plot(ax102, [considbound, considbound], [0,obj.valid_blocks.num_trials_per_block],'b-', 'DisplayName', 'Consideration Limit')
					plot(ax102, [rb_ms/1000, rb_ms/1000], [0,obj.valid_blocks.num_trials_per_block],'k-', 'DisplayName', 'Reward Boundary')
					plot(ax102, [obj.time.s.ITI_time_s, obj.time.s.ITI_time_s], [0,obj.valid_blocks.num_trials_per_block],'k-', 'DisplayName', 'ITI Start')
					plot(ax102, [rxnwin_s, rxnwin_s], [0, numel(obj.time.s.binnedLicks.refevents)],'g-', 'DisplayName', 'Permitted Reaction Window')
					ylim(ax102, yy)
					title(ax102, ['Trials #', num2str(obj.valid_blocks.all_trials_this_block(1)), '-', num2str(obj.valid_blocks.all_trials_this_block(end))])
					xlim(ax102, [-5,18])
					ax102.XLabel.String = ['All Lick Times (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax102.YLabel.String = 'Percentage of Licks Across Session';



				elseif strcmpi(Mode, 'X')
					warning('Possible Problem? -- X data for trials with no-lick is entered as nans to the concat, and the whole cue-to-lick interval is included')
					ax9 = subplot(1, 2,  1, 'Parent', f);
					hold(ax9, 'on')
					ax9_2 = subplot(1, 2,  2, 'Parent', f);

					hold(ax9_2, 'on')
					plot(ax9, cell2mat(obj.valid_blocks.m.data{1, i_div}{1,1}.'), 'LineWidth', 1, 'DisplayName', 'stim')
					title(ax9, ['Xconcat Stim Trials #', num2str(obj.valid_blocks.all_trials_this_block(1)), '-', num2str(obj.valid_blocks.all_trials_this_block(end))])				

					xx = get(ax9, 'xlim');
					plot(ax9_2, cell2mat(obj.valid_blocks.m.data{1, i_div}{1,2}.'), 'LineWidth', 1, 'DisplayName', 'nostim')
					title(ax9_2, ['Xconcat No-Stim Trials #', num2str(obj.valid_blocks.all_trials_this_block(1)), '-', num2str(obj.valid_blocks.all_trials_this_block(end))])				
					xlim(ax9_2, xx)

					linkaxes([ax9, ax9_2], 'y')
				else
					error('Invalid Mode');
				end
			end
		end		
		function ax = plotNTJ(obj, Mode, ref, rxnwin_s)
			% 
			% 	>> plot >> plotNTJ
			% 	can also run alone as
			% 	>> plotNTJ
			% 
			% 	Mode options:
			% 		raster4 -- plots a raster of stimulated and unstimulated licks. top row is flick only, bottom is all licks
			% 		ecdf 	-- plots the full ecdf
			% 		hxg 	-- plots histogram and cdf of flicks and all licks (normalized by total # of licks in each category)
			% 		stimHxg -- plots only the stim hxg aligned to stim-on, excluding trials with juice at -0.5s
			% 		hxg-counts plots histogram without normalizing for # of trials
			% 		X 		-- plots concatenated X for cue-to-lick for all trials where lick occurs before EOT
			% 
			% 	Blocking: select either valid or even
			% 
			% 	division: the block # to plot for.
			% 
			if nargin < 2
				Mode = 'raster4';
			end
			if nargin < 3
				ref = 'stimOn';
			end
			if nargin < 4
				rxnwin_s = 0;
			end
			if ~strcmpi(obj.iv.stim_protocol, 'ntj')
				obj.plot(Mode, 'valid',1,ref, rxnwin_s);
			else
				obj.getBinnedLicks(ref, 30000, 30000, rxnwin_s);
				% 
				% 	Collect var pntrs
				% 
				fs_ms = obj.iv.samples_per_ms_ChR2;
				fs_X_ms = obj.iv.samples_per_ms_X;
				fs_EMG_ms = obj.iv.samples_per_ms_EMG;
				fs_CamO_ms = obj.iv.samples_per_ms_CamO;
				sb_s = obj.time.Plot.stim_on_time_s;
				se_s = obj.time.Plot.stim_off_time_s;
				rfx_s = obj.time.Plot.refrx_time_s;
				max_s = obj.time.Plot.max_time_s;
				considbound = obj.bootStat.considerationBoundary_ms/1000;
				
				f = figure;
				set(f, 'color', 'white')
				if strcmpi(Mode, 'raster4')
					obj.getBinnedLicks(ref, 5000, 5000, rxnwin_s);
					ax = subplot(1, 1,  1, 'Parent', f);
					hold(ax, 'on')
					nEvents = numel(obj.time.s.binnedLicks.f_lick_ex_s_wrtref);
					% 
					% 	Plot raster of all licks with first licks overlaid
					% 
					plot(ax, [0,0], [1,nEvents],'k-', 'DisplayName', obj.time.s.binnedLicks.ref)

					set(ax,  'YDir','reverse')
					title(ax, ['Lick Raster Aligned to ' obj.time.s.binnedLicks.ref])
					xlim(ax, [-10, max_s+1])
			        if nEvents > 1
			    		ylim(ax, [1, nEvents])
			        end
					ax.XLabel.String = ['Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
					ax.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];


				    lick_times = obj.time.s.binnedLicks.lick_ex_s;				    
				    % 
					plot(ax, obj.time.s.binnedLicks.f_lick_ex_s_wrtref, 1:numel(obj.time.s.binnedLicks.f_lick_ex_s_wrtref), 'mo', 'DisplayName', 'First Lick', 'markerfacecolor', 'm')
					%
					for itrial = 1:nEvents
						plotpnts_lick = lick_times{itrial};
						if ~isempty(plotpnts_lick)
							plot(ax, plotpnts_lick, itrial.*ones(numel(plotpnts_lick), 1),'k.')
						end	
					end
				elseif strcmpi(Mode, 'ecdf')
					error('not implemented')
				elseif strcmpi(Mode, 'hxg') || strcmpi(Mode, 'hxg-counts')
					% obj.getBinnedLicks(ref, 5000, 5000, rxnwin_s);
					axNOPreStimJuiceHXG = subplot(3, 2,  1, 'Parent', f);
					hold(axNOPreStimJuiceHXG, 'on')
					axNOPreStimJuiceFLICK = subplot(3, 2,  2, 'Parent', f);
					hold(axNOPreStimJuiceFLICK, 'on')
					axPreStimJuiceHXG = subplot(3, 2,  3, 'Parent', f);
					hold(axPreStimJuiceHXG, 'on')
					axPreStimJuiceFLICK = subplot(3, 2,  4, 'Parent', f);
					hold(axPreStimJuiceFLICK, 'on')
					obj.NTJcategories;
					all_licks_swrtSO_ex_s.all = obj.lick.swrtStimOn.lick_ex_s;
					all_licks_swrtSO_ex_s.preStimJuice = cell2mat(all_licks_swrtSO_ex_s.all(obj.valid_blocks.preStimJuiceTrials));
					all_licks_swrtSO_ex_s.NOpreStimJuice = cell2mat(all_licks_swrtSO_ex_s.all(obj.valid_blocks.NOpreStimJuiceTrials));
					all_flicks_swrtSO_ex_s.all = obj.lick.swrtStimOn.f_lick_ex_s;
					all_flicks_swrtSO_ex_s.preStimJuice = all_flicks_swrtSO_ex_s.all(obj.valid_blocks.preStimJuiceTrials);
					all_flicks_swrtSO_ex_s.NOpreStimJuice = all_flicks_swrtSO_ex_s.all(obj.valid_blocks.NOpreStimJuiceTrials);
					% 	axNOPreStimJuiceHXG
					histogram(axNOPreStimJuiceHXG, all_licks_swrtSO_ex_s.NOpreStimJuice(~isnan(all_licks_swrtSO_ex_s.NOpreStimJuice)), 30000,'displayName','licks');
					axNOPreStimJuiceHXG.YLabel.String = ['# of Licks wrt ' ref];					
					obj.overlayEventsNTJ(axNOPreStimJuiceHXG,ref);
					xlim(axNOPreStimJuiceHXG, [-15, 15])
					title(axNOPreStimJuiceHXG, 'Aligned to StimOn, NO pre-Juice @-0.5s')
					axNOPreStimJuiceHXG.XLabel.String = ['Lick Time (s wrt ' ref ')'];
					% 	axNOPreStimJuiceHXG
					histogram(axNOPreStimJuiceFLICK, all_flicks_swrtSO_ex_s.NOpreStimJuice(~isnan(all_flicks_swrtSO_ex_s.NOpreStimJuice)), 100,'displayName','first licks');
					axNOPreStimJuiceFLICK.YLabel.String = ['# of First Licks wrt ' ref];					
					obj.overlayEventsNTJ(axNOPreStimJuiceFLICK,ref);
					legend(axNOPreStimJuiceFLICK, 'show')
					xlim(axNOPreStimJuiceFLICK, [0, 20])
					title(axNOPreStimJuiceFLICK, 'Aligned to StimOn, NO pre-Juice @-0.5s')
					disp('*****NB: the histogram currently includes juice delivered at random in ROW 1, the No Pre-stim juice panels')
					axNOPreStimJuiceFLICK.XLabel.String = ['First Lick Time (s wrt ' ref ')'];
					% 	axPreStimJuiceHXG
					histogram(axPreStimJuiceHXG, all_licks_swrtSO_ex_s.preStimJuice(~isnan(all_licks_swrtSO_ex_s.preStimJuice)), 30000,'displayName','licks');
					axPreStimJuiceHXG.YLabel.String = ['# of Licks wrt ' ref];					
					obj.overlayEventsNTJ(axPreStimJuiceHXG,ref);
					xlim(axPreStimJuiceHXG, [-15, 15])
					title(axPreStimJuiceHXG, 'Aligned to StimOn, + pre-Juice @-0.5s')
					axPreStimJuiceHXG.XLabel.String = ['Lick Time (s wrt ' ref ')'];
					% 	axPreStimJuiceFLICK
					histogram(axPreStimJuiceFLICK, all_flicks_swrtSO_ex_s.preStimJuice(~isnan(all_flicks_swrtSO_ex_s.preStimJuice)), 100,'displayName','first licks');
					axPreStimJuiceFLICK.YLabel.String = ['# of First Licks wrt ' ref];					
					obj.overlayEventsNTJ(axPreStimJuiceFLICK,ref);
					xlim(axPreStimJuiceFLICK, [0, 20])
					title(axPreStimJuiceFLICK, 'Aligned to StimOn, + pre-Juice @-0.5s')
					axPreStimJuiceFLICK.XLabel.String = ['First Lick Time (s wrt ' ref ')'];
					% 
					% 	Now align to rewards...
					% 
					axRefJuiceAtStimOn = subplot(3, 2,  5, 'Parent', f);
					hold(axRefJuiceAtStimOn, 'on')
					axRefJuiceNotAtStimOnORDuringStim = subplot(3, 2,  6, 'Parent', f);
					hold(axRefJuiceNotAtStimOnORDuringStim, 'on')
					ref = 'juice';
					all_licks_swrtjuice_s.all = obj.lick.swrtJuice.lick_s;
					all_licks_swrtjuice_s.alignedToPreStimOn = cell2mat(all_licks_swrtjuice_s.all(obj.valid_blocks.Juice_Idxs_Pre_Stim));
					all_licks_swrtjuice_s.notInStimOrDuringPreStim = cell2mat(all_licks_swrtjuice_s.all(obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim));
					% 	axRefJuiceAtStimOn
					histogram(axRefJuiceAtStimOn, all_licks_swrtjuice_s.alignedToPreStimOn(~isnan(all_licks_swrtjuice_s.alignedToPreStimOn)), 30000,'displayName','licks');
					axRefJuiceAtStimOn.YLabel.String = ['# of Licks wrt ' ref];					
					plot(axRefJuiceAtStimOn, [0, 0], [0,numel(obj.time.s.juice_s)],'r-', 'DisplayName', 'Juice 0.5 s Before Stimilation')
					xlim(axRefJuiceAtStimOn, [-15, 15])
					title(axRefJuiceAtStimOn,'Licks aligned to Pre-Stim Juice @-0.5s')
					axRefJuiceAtStimOn.XLabel.String = ['Lick Time (s wrt ' ref ')'];
					% 	axRefJuiceNotAtStimOnORDuringStim
					histogram(axRefJuiceNotAtStimOnORDuringStim, all_licks_swrtjuice_s.notInStimOrDuringPreStim(~isnan(all_licks_swrtjuice_s.notInStimOrDuringPreStim)), 30000,'displayName','licks');
					axRefJuiceNotAtStimOnORDuringStim.YLabel.String = ['# of Licks wrt ' ref];					
					plot(axRefJuiceNotAtStimOnORDuringStim, [0, 0], [0,numel(obj.time.s.juice_s)],'r-', 'DisplayName', 'Juice not during prestim OR during stim')
					xlim(axRefJuiceNotAtStimOnORDuringStim, [-15, 15])
					title(axRefJuiceNotAtStimOnORDuringStim,'Licks aligned to Juice Not PreStim or During Stim')
					axRefJuiceNotAtStimOnORDuringStim.XLabel.String = ['Lick Time (s wrt ' ref ')'];
					ax = f;


				elseif strcmpi(Mode, 'X')
					error('not implemented')
					warning('Possible Problem? -- X data for trials with no-lick is entered as nans to the concat, and the whole cue-to-lick interval is included')
					ax9 = subplot(1, 2,  1, 'Parent', f);
					hold(ax9, 'on')
					ax9_2 = subplot(1, 2,  2, 'Parent', f);

					hold(ax9_2, 'on')
					plot(ax9, cell2mat(obj.valid_blocks.m.data{1, i_div}{1,1}.'), 'LineWidth', 1, 'DisplayName', 'stim')
					title(ax9, ['Xconcat Stim Trials #', num2str(obj.valid_blocks.all_trials_this_block(1)), '-', num2str(obj.valid_blocks.all_trials_this_block(end))])				

					xx = get(ax9, 'xlim');
					plot(ax9_2, cell2mat(obj.valid_blocks.m.data{1, i_div}{1,2}.'), 'LineWidth', 1, 'DisplayName', 'nostim')
					title(ax9_2, ['Xconcat No-Stim Trials #', num2str(obj.valid_blocks.all_trials_this_block(1)), '-', num2str(obj.valid_blocks.all_trials_this_block(end))])				
					xlim(ax9_2, xx)

					linkaxes([ax9, ax9_2], 'y')
				elseif strcmpi(Mode,'stimHxg')
					axNOPreStimJuiceHXG = subplot(1, 1,  1, 'Parent', f);
					hold(axNOPreStimJuiceHXG, 'on')
					set(axNOPreStimJuiceHXG,'fontsize',20);    
					obj.NTJcategories;

					all_licks_swrtSO_ex_s.all = obj.lick.swrtStimOn.lick_ex_s;
					ntrialsNoPreJuice = numel(all_licks_swrtSO_ex_s.all(obj.valid_blocks.NOpreStimJuiceTrials));
					yy = [0, round(0.5*ntrialsNoPreJuice)];
					tt = ['Trials #1-', num2str(numel(all_licks_swrtSO_ex_s.all)) ' | nNoPreJuice: ' num2str(numel(ntrialsNoPreJuice))];

					all_licks_swrtSO_ex_s.NOpreStimJuice = cell2mat(all_licks_swrtSO_ex_s.all(obj.valid_blocks.NOpreStimJuiceTrials));
					% 	axNOPreStimJuiceHXG
					histogram(axNOPreStimJuiceHXG, all_licks_swrtSO_ex_s.NOpreStimJuice(~isnan(all_licks_swrtSO_ex_s.NOpreStimJuice)), 'binwidth', 0.25,'displayName', 'No Pre-Juice','displaystyle','stairs', 'LineWidth',6);
					plot(axNOPreStimJuiceHXG, [0, 0], yy,'r-', 'DisplayName', 'Stim On')
					plot(axNOPreStimJuiceHXG, [3, 3], yy,'r-', 'DisplayName', 'Stim Off')
					axNOPreStimJuiceHXG.YLabel.String = 'Lick Rate (Licks/Trial)';	
					axNOPreStimJuiceHXG.XLabel.String = 'time (s)';
					xlim(axNOPreStimJuiceHXG, [-10, 18])
					ylim(axNOPreStimJuiceHXG,yy);
					title(axNOPreStimJuiceHXG, tt);
					
					ax = axNOPreStimJuiceHXG;

					yticks(axNOPreStimJuiceHXG, [yy(1),yy(2)])
					obj.bootStat.nTrialsInBlockStim = ntrialsNoPreJuice;
					obj.bootStat.ylimStim = yy;
					obj.bootStat.nTrialsInBlockUnstim = 0;
					obj.rescaleYAxisStim(axNOPreStimJuiceHXG);
				else
					error('Invalid Mode');
				end
			end
		end
		function overlayEventsNTJ(obj,ax,ref)
			fs_ms = obj.iv.samples_per_ms_ChR2;
			fs_X_ms = obj.iv.samples_per_ms_X;
			fs_EMG_ms = obj.iv.samples_per_ms_EMG;
			fs_CamO_ms = obj.iv.samples_per_ms_CamO;
			sb_s = obj.time.Plot.stim_on_time_s;
			se_s = obj.time.Plot.stim_off_time_s;
			rfx_s = obj.time.Plot.refrx_time_s;
			max_s = obj.time.Plot.max_time_s;
			considbound = obj.bootStat.considerationBoundary_ms/1000;
			yy = get(ax, 'ylim');
			plot(ax, [0, 0], [0,obj.iv.num_trials],'r-', 'DisplayName', ref)
			plot(ax, [se_s, se_s], [0,obj.iv.num_trials],'k-', 'DisplayName', 'Stimulation End')
			plot(ax, [rfx_s, rfx_s], [0,obj.iv.num_trials],'k-', 'DisplayName', 'Refractory Period')
			plot(ax, [max_s, max_s], [0, numel(obj.time.s.binnedLicks.refevents)],'k-', 'DisplayName', 'Max Time')
			plot(ax, [considbound, considbound], [0, numel(obj.time.s.binnedLicks.refevents)],'k-', 'DisplayName', 'Consideration Boundary')
			ylim(ax, yy)
		end
		function addEventToRasterNTJ(obj, ax, event)
			% 
			% 	>> plot >> plotNTJ >> addEventToRasterNTJ
			% 
			% 	Helper tool for plotNTJ -- allows you to add more events to the plot.
			% 
			if nargin < 3
				event = 'juice';
			end
			nEvents = numel(obj.time.s.binnedLicks.refevents);
			if strcmpi(event, 'stimOn')
				event_times = obj.time.s.binnedLicks.stimOn_s_bt_wrtref;
				cc = {'co','c'};
			elseif strcmpi(event, 'stimOff')
				event_times = obj.time.s.binnedLicks.stimOff_s_bt_wrtref;
				cc = {'ko','c'};
			elseif strcmpi(event, 'refrxOff')
				event_times = obj.time.s.binnedLicks.refrxOff_s_bt_wrtref;
				cc = {'ro','r'};
			elseif strcmpi(event, 'juice')
				event_times = obj.time.s.binnedLicks.juice_s_bt_wrtref;
				cc = {'go','g'};
			end 
			for itrial = 1:nEvents			
				plotpnts = event_times{itrial};
				if ~isempty(plotpnts)
					plot(ax, plotpnts, itrial.*ones(numel(plotpnts), 1),cc{1}, 'markerfacecolor', cc{2})
				end	
			end
		end
		function data = ntjRasterCollateHelper(obj)
			% 
			% 	NTJ + juice
			% 
			fs_ms = obj.iv.samples_per_ms_ChR2;
			sb_s = obj.time.Plot.stim_on_time_s;
			se_s = obj.time.Plot.stim_off_time_s;
			rfx_s = obj.time.Plot.refrx_time_s;
			max_s = obj.time.Plot.max_time_s;
			considbound = obj.bootStat.considerationBoundary_ms/1000;
				
		    data.fs_ms = fs_ms;
		    data.sb_s = sb_s;
		    data.se_s = se_s;
		    data.rfx_s = rfx_s;
		    data.max_s = max_s;
		    data.considbound = considbound;
		    data.time = obj.time;

		    obj.NTJcategories;
			all_licks_swrtSO_ex_s.all = obj.lick.swrtStimOn.lick_ex_s;
			all_licks_swrtSO_ex_s.preStimJuice = cell2mat(all_licks_swrtSO_ex_s.all(obj.valid_blocks.preStimJuiceTrials));
			all_licks_swrtSO_ex_s.NOpreStimJuice = cell2mat(all_licks_swrtSO_ex_s.all(obj.valid_blocks.NOpreStimJuiceTrials));
			all_flicks_swrtSO_ex_s.all = obj.lick.swrtStimOn.f_lick_ex_s;
			all_flicks_swrtSO_ex_s.preStimJuice = all_flicks_swrtSO_ex_s.all(obj.valid_blocks.preStimJuiceTrials);
			all_flicks_swrtSO_ex_s.NOpreStimJuice = all_flicks_swrtSO_ex_s.all(obj.valid_blocks.NOpreStimJuiceTrials);
			all_licks_swrtjuice_s.all = obj.lick.swrtJuice.lick_s;
			all_licks_swrtjuice_s.alignedToPreStimOn = cell2mat(all_licks_swrtjuice_s.all(obj.valid_blocks.Juice_Idxs_Pre_Stim));
			all_licks_swrtjuice_s.notInStimOrDuringPreStim = cell2mat(all_licks_swrtjuice_s.all(obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim));
            data.all_licks_swrtSO_ex_s = all_licks_swrtSO_ex_s;
            data.all_flicks_swrtSO_ex_s = all_flicks_swrtSO_ex_s;
            data.all_licks_swrtjuice_s = all_licks_swrtjuice_s;
            data.valid_blocks = obj.valid_blocks;
            obj.getBinnedLicks('stimOn', 30000, 30000, 0);
            data.juice_s_bt_wrtref = obj.time.s.binnedLicks.juice_s_bt_wrtref;
		end
		function data = ntRasterCollateHelper(obj)
			ref = 'cue';
			rxnwin_s = 0;
			obj.getBinnedLicks(ref, 30000, 30000, rxnwin_s);
			i_div = 1;
			if i_div ~=1 || ~(strcmpi(obj.valid_blocks.Mode,'fullday') || strcmpi(obj.valid_blocks.Mode,'full-day'))
				error('Not implemented at obj.valid_blocks.all_trials_this_block @getvalidblocks')
			end
			% 
			% 	Collect var pntrs
			% 
			fs_ms = obj.iv.init_variables_ChR2.time_parameters.samples_per_ms;
			rb_ms = obj.time.ms.op_rew_open_ms;
			eot = obj.time.s.ITI_time_s;
			cue = obj.time.ms.cue_on_time_ms;
			considbound = obj.bootStat.considerationBoundary_ms/1000;
			all_trials_this_block = obj.valid_blocks.all_trials_this_block;
			all_lick_times_ex_swrtc = obj.lick.swrtc.lick_ex_times_by_trial_swrtc;
			for iexc = obj.iv.ExcludedTrials
		    	all_lick_times_ex_swrtc{iexc} = [];
	    	end
			f_lick_s_wrtcue.stim = obj.valid_blocks.f_lick_s_wrtcue(i_div).stim;
			f_lick_s_wrtcue.unstim = obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim;
			stimTrials = obj.h.stimTrials;
			nostim_trials = obj.h.nostim_trials;

		    data.ref = ref;
		    data.rxnwin_s = rxnwin_s;
		    data.fs_ms = fs_ms;
		    data.rb_ms = rb_ms;
		    data.eot = eot;
		    data.cue = cue;
		    data.considbound = considbound;
		    data.all_trials_this_block = all_trials_this_block;
		    data.all_lick_times_ex_swrtc = all_lick_times_ex_swrtc;
		    data.f_lick_s_wrtcue = f_lick_s_wrtcue;
		    data.stimTrials = stimTrials;
		    data.nostim_trials = nostim_trials;
	    end
	    function ax = plotCollatedNT(obj, Mode, ax, IDs)
	    	% 
	    	% 	Use on collated naive test - juice datasets to plot the raster and histogram
	    	% 
	    	% 	Mode: 	'r&h' -- plots a raster and prob hxg together
	    	% 			'raster4' 	  -- plots composite rasters
	    	% 			'hxg' 		  -- plots composite hxg
	    	% 			'stimHxg'     -- plots only the stimulated hxg in units of licks per trial
	    	% 
	    	if nargin < 2
	    		Mode = 'r&h';
    		end
    		if nargin < 4
    			IDs = 1:numel({obj.bootStat.collateResult.ID});
			end

	    	
			rxnwin_s = 0; 
			fs_ms = obj.iv.init_variables_ChR2.time_parameters.samples_per_ms;
			rb_ms = obj.time.ms.op_rew_open_ms;
			eot = obj.time.s.ITI_time_s;
			cue = obj.time.ms.cue_on_time_ms;
			considbound = obj.bootStat.considerationBoundary_ms/1000;
			% 
			%	Get numbers of valid trials across all sets 
			% 
			all_trials_this_block = [];
			stimTrials = [];
			nostim_trials = [];
			all_lick_times_ex_swrtc = [];
			alllicktimes_stimwrtc = [];
			alllicktimes_unstimwrtc = [];
			f_lick_s_wrtcue.stim = [];
			f_lick_s_wrtcue.unstim = [];
			for iID = IDs
				if numel(obj.bootStat.collateResult(iID).all_lick_times_ex_swrtc) ~= numel(obj.bootStat.collateResult(iID).all_trials_this_block)
					thisID.all_lick_times_ex_swrtc = obj.bootStat.collateResult(iID).all_lick_times_ex_swrtc(1:end-1);
				else
					thisID.all_lick_times_ex_swrtc = obj.bootStat.collateResult(iID).all_lick_times_ex_swrtc;
                end
                if size(thisID.all_lick_times_ex_swrtc,1) ~= 1
                    thisID.all_lick_times_ex_swrtc = thisID.all_lick_times_ex_swrtc';%'
                end

				if isempty(all_trials_this_block)
					stimTrials(1:numel(obj.bootStat.collateResult(iID).stimTrials)) = obj.bootStat.collateResult(iID).stimTrials;
					nostim_trials(1:numel(obj.bootStat.collateResult(iID).nostim_trials)) = obj.bootStat.collateResult(iID).nostim_trials;
					all_trials_this_block(1:numel(obj.bootStat.collateResult(iID).all_trials_this_block)) = obj.bootStat.collateResult(iID).all_trials_this_block;
					all_lick_times_ex_swrtc = thisID.all_lick_times_ex_swrtc;
					alllicktimes_stimwrtc = thisID.all_lick_times_ex_swrtc(obj.bootStat.collateResult(iID).stimTrials);
					alllicktimes_unstimwrtc = thisID.all_lick_times_ex_swrtc(obj.bootStat.collateResult(iID).nostim_trials);
					f_lick_s_wrtcue.stim = obj.bootStat.collateResult(iID).f_lick_s_wrtcue.stim;
					f_lick_s_wrtcue.unstim = obj.bootStat.collateResult(iID).f_lick_s_wrtcue.unstim;
				else
					stimTrials(end+1:end+numel(obj.bootStat.collateResult(iID).stimTrials)) = all_trials_this_block(end)+obj.bootStat.collateResult(iID).stimTrials;
					nostim_trials(end+1:end+numel(obj.bootStat.collateResult(iID).nostim_trials)) = all_trials_this_block(end)+obj.bootStat.collateResult(iID).nostim_trials;
					all_trials_this_block(end+1:end+numel(obj.bootStat.collateResult(iID).all_trials_this_block)) = all_trials_this_block(end)+obj.bootStat.collateResult(iID).all_trials_this_block;
					all_lick_times_ex_swrtc = [all_lick_times_ex_swrtc, thisID.all_lick_times_ex_swrtc];
					alllicktimes_stimwrtc = [alllicktimes_stimwrtc, thisID.all_lick_times_ex_swrtc(obj.bootStat.collateResult(iID).stimTrials)];
					alllicktimes_unstimwrtc = [alllicktimes_unstimwrtc, thisID.all_lick_times_ex_swrtc(obj.bootStat.collateResult(iID).nostim_trials)];
					f_lick_s_wrtcue.stim = [f_lick_s_wrtcue.stim; obj.bootStat.collateResult(iID).f_lick_s_wrtcue.stim];
					f_lick_s_wrtcue.unstim = [f_lick_s_wrtcue.unstim; obj.bootStat.collateResult(iID).f_lick_s_wrtcue.unstim];
				end
			end
			% 
			if nargin < 3 || isempty(ax)
				f = figure;
				set(f, 'color', 'white');	
				if numel(IDs) == 1
					set(f, 'name', obj.bootStat.collateResult(IDs).ID)
				end
			end
			if strcmpi(Mode,'raster4')
				ax11nostim = subplot(2, 2,  1, 'Parent', f);
				hold(ax11nostim, 'on')
				ax11stim = subplot(2, 2,  2, 'Parent', f);
				hold(ax11stim, 'on')
				ax11ONLYnostim = subplot(2, 2,  3, 'Parent', f);
				hold(ax11ONLYnostim, 'on')
				ax11ONLYstim = subplot(2, 2,  4, 'Parent', f);
				hold(ax11ONLYstim, 'on')
				% 
				% 	Plot raster of all licks with first licks overlaid
				% 
				plot(ax11nostim, zeros(size(all_trials_this_block)), all_trials_this_block,'r-', 'DisplayName', 'Cue')
			    plot(ax11nostim, [considbound, considbound], [1, all_trials_this_block(end)],'b-', 'DisplayName', 'Consideration Limit')
				plot(ax11nostim, [rb_ms/1000, rb_ms/1000], [1, all_trials_this_block(end)],'k-', 'DisplayName', 'Reward Boundary')
				plot(ax11nostim, [eot, eot], [1, all_trials_this_block(end)],'k-', 'DisplayName', 'ITI Start')
				plot(ax11nostim, [rxnwin_s, rxnwin_s], [1, all_trials_this_block(end)],'g-', 'DisplayName', 'Permitted Reaction Window')
				set(ax11nostim,  'YDir','reverse')
				title(ax11nostim, ['No Stim Lick Raster Aligned to cue'])
				xlim(ax11nostim, [0, obj.time.ms.total_time_ms/1000])
		        if numel(all_trials_this_block(end)) > 1
		    		ylim(ax11nostim, [1, all_trials_this_block(end)])
		        end
				ax11nostim.XLabel.String = ['Time (s wrt cue)'];
				ax11nostim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];
				plot(ax11stim, zeros(size(all_trials_this_block)), all_trials_this_block,'r-', 'DisplayName', 'Cue')
			    plot(ax11stim, [considbound, considbound], [all_trials_this_block(1), all_trials_this_block(end)],'b-', 'DisplayName', 'Consideration Limit')
				plot(ax11stim, [rb_ms/1000, rb_ms/1000], [all_trials_this_block(1), all_trials_this_block(end)],'k-', 'DisplayName', 'Reward Boundary')
				plot(ax11stim, [eot, eot], [all_trials_this_block(1), all_trials_this_block(end)],'k-', 'DisplayName', 'ITI Start')
				set(ax11stim,  'YDir','reverse')
				title(ax11stim, '+ stimulation')
				xlim(ax11stim, [0, obj.time.ms.total_time_ms/1000])
		        if numel(all_trials_this_block(end)) > 1
		    		ylim(ax11stim, [1, all_trials_this_block(end)])
		        end
				ax11stim.XLabel.String = ['Time (s wrt cue)'];
				ax11stim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];
				plot(ax11ONLYnostim, zeros(size(all_trials_this_block)), numel(nostim_trials),'r-', 'DisplayName', 'Cue')
			    plot(ax11ONLYnostim, [considbound, considbound], [1, numel(nostim_trials)],'b-', 'DisplayName', 'Consideration Limit')
				plot(ax11ONLYnostim, [rb_ms/1000, rb_ms/1000], [1, numel(nostim_trials)],'k-', 'DisplayName', 'Reward Boundary')
				plot(ax11ONLYnostim, [eot, eot], [1, numel(nostim_trials)],'k-', 'DisplayName', 'ITI Start')	
				plot(ax11ONLYstim, [0,obj.time.ms.total_time_ms/1000], [numel(nostim_trials), numel(nostim_trials)], 'k-', 'DisplayName', '# no stim trials')				
				set(ax11ONLYnostim,  'YDir','reverse')
				title(ax11ONLYnostim, 'ONLY Consec Trials, No stim')
				xlim(ax11ONLYnostim, [0, obj.time.ms.total_time_ms/1000])
		        if numel(nostim_trials) > 1
		    		ylim(ax11ONLYnostim, [1, numel(nostim_trials)])
		        end
				ax11ONLYnostim.XLabel.String = ['Time (s wrt cue)'];
				ax11ONLYnostim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];
				plot(ax11ONLYstim, zeros(size(all_trials_this_block)), numel(stimTrials),'r-', 'DisplayName', 'Cue')
			    plot(ax11ONLYstim, [considbound, considbound], [1, numel(stimTrials)],'b-', 'DisplayName', 'Consideration Limit')
				plot(ax11ONLYstim, [rb_ms/1000, rb_ms/1000], [1, numel(stimTrials)],'k-', 'DisplayName', 'Reward Boundary')
				plot(ax11ONLYstim, [eot, eot], [1, numel(stimTrials)],'k-', 'DisplayName', 'ITI Start')
				plot(ax11ONLYstim, [0,obj.time.ms.total_time_ms/1000], [numel(stimTrials), numel(stimTrials)], 'k-', 'DisplayName', '# stim trials')
						
				set(ax11ONLYstim,  'YDir','reverse')
				title(ax11ONLYstim, 'ONLY Consec Trials, + stim')
				xlim(ax11ONLYstim, [0, obj.time.ms.total_time_ms/1000])
		        if numel(stimTrials) > 1
		    		% ylim(ax11ONLYstim, [1, numel(stimTrials)])
		    		ylim(ax11ONLYstim, [1, numel(nostim_trials)])
		        end
				ax11ONLYstim.XLabel.String = ['Time (s wrt cue)'];
				ax11ONLYstim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];

			    unstimTrialTracker = 1;
			    stimTrialTracker = 1;
			    
				plot(ax11stim, f_lick_s_wrtcue.stim, all_trials_this_block, 'mo', 'DisplayName', 'First Lick', 'markerfacecolor', 'm')
				plot(ax11nostim, f_lick_s_wrtcue.unstim, all_trials_this_block, 'mo', 'DisplayName', 'First Lick', 'markerfacecolor', 'm')
				for itrial = all_trials_this_block(1):all_trials_this_block(end)
					obj.valid_blocks.plotpnts = all_lick_times_ex_swrtc{itrial};
					if ismember(itrial, stimTrials)
						if ~isempty(obj.valid_blocks.plotpnts) && numel(f_lick_s_wrtcue.stim)>=itrial
							plot(ax11stim, obj.valid_blocks.plotpnts, itrial.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
							plot(ax11ONLYstim, obj.valid_blocks.plotpnts, stimTrialTracker.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
							plot(ax11ONLYstim, f_lick_s_wrtcue.stim(itrial), stimTrialTracker, 'mo', 'DisplayName', 'First Lick')
						end
						stimTrialTracker = stimTrialTracker +1;
					elseif ismember(itrial, nostim_trials) && numel(f_lick_s_wrtcue.unstim)>=itrial
						if ~isempty(obj.valid_blocks.plotpnts)
							plot(ax11nostim, obj.valid_blocks.plotpnts, itrial.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
							plot(ax11ONLYnostim, obj.valid_blocks.plotpnts, unstimTrialTracker.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
							plot(ax11ONLYnostim, f_lick_s_wrtcue.unstim(itrial), unstimTrialTracker, 'mo', 'DisplayName', 'First Lick')
						end
						unstimTrialTracker = unstimTrialTracker +1;
					else
						error('shouldn''t get here');
					end					
				end
				yy = get(ax11ONLYnostim, 'ylim');
				ylim(ax11ONLYstim, yy)
	            xx = get(ax11nostim, 'xlim');
	            plot(ax11nostim, xx, [numel(f_lick_s_wrtcue.unstim)+0.5, numel(f_lick_s_wrtcue.unstim)+0.5], 'r-', 'displayName', 'Consecutive Trials limit for AUC')
	            plot(ax11stim, xx, [numel(f_lick_s_wrtcue.stim)+0.5, numel(f_lick_s_wrtcue.stim)+0.5], 'r-', 'displayName', 'Consecutive Trials limit for AUC')
	            ax = ax11stim;
            elseif strcmpi(Mode, 'hxg')
	            ax10 = subplot(1, 3,  2, 'Parent', f);
				hold(ax10, 'on')
				ax103 = subplot(1, 3,  1, 'Parent', f);
				hold(ax103, 'on')
				ax102 = subplot(1, 3,  3, 'Parent', f);
				hold(ax102, 'on')
				histogram(ax10, f_lick_s_wrtcue.unstim(~isnan(f_lick_s_wrtcue.unstim)), 30, 'DisplayName', 'No Stimulation', 'Normalization', 'probability')
				histogram(ax10, f_lick_s_wrtcue.stim(~isnan(f_lick_s_wrtcue.stim)), 30,'DisplayName', '+ Stimulation', 'Normalization', 'probability')
				ax10.YLabel.String = 'Percentage of First Licks Across Category';
				yy = get(ax10, 'ylim');
				plot(ax10, [0, 0], yy,'r-', 'DisplayName', 'Cue')
			    plot(ax10, [considbound, considbound], yy,'b-', 'DisplayName', 'Consideration Limit')
				plot(ax10, [rb_ms/1000, rb_ms/1000], yy,'k-', 'DisplayName', 'Reward Boundary')
				plot(ax10, [eot, eot], yy,'k-', 'DisplayName', 'ITI Start')
				plot(ax10, [rxnwin_s, rxnwin_s], yy,'g-', 'DisplayName', 'Permitted Reaction Window')
				ylim(ax10, yy)
				title(ax10, ['Trials #', num2str(all_trials_this_block(1)), '-', num2str(all_trials_this_block(end))])
				legend(ax10, 'show')
				xlim(ax10, [-.4, obj.time.ms.total_time_ms/1000])
				ax10.XLabel.String = 'First Lick Time (s wrt cue)';
				histogram(ax103, f_lick_s_wrtcue.unstim(~isnan(f_lick_s_wrtcue.unstim)), 30, 'DisplayName', 'No Stimulation', 'Normalization', 'cdf')
				histogram(ax103, f_lick_s_wrtcue.stim(~isnan(f_lick_s_wrtcue.stim)), 30,'DisplayName', '+ Stimulation', 'Normalization', 'cdf')
				yy = get(ax103, 'ylim');
				plot(ax103, [0, 0], yy,'r-', 'DisplayName', 'Cue')
			    plot(ax103, [considbound, considbound], yy,'b-', 'DisplayName', 'Consideration Limit')
				plot(ax103, [rb_ms/1000, rb_ms/1000], yy,'k-', 'DisplayName', 'Reward Boundary')
				plot(ax103, [obj.time.s.ITI_time_s, obj.time.s.ITI_time_s], yy,'k-', 'DisplayName', 'ITI Start')
				plot(ax103, [rxnwin_s, rxnwin_s], yy,'g-', 'DisplayName', 'Permitted Reaction Window')
				ylim(ax103, yy)
				title(ax103, ['Histogram of first licks Aligned to ' obj.time.s.binnedLicks.ref])
				xlim(ax103, [-.4, obj.time.ms.total_time_ms/1000])
				ax103.XLabel.String = ['First Lick Time (s wrt ' obj.time.s.binnedLicks.ref ')'];
				ax103.YLabel.String = 'CDF of First Licks Across Category';

				warning('Hxg of all licks might not incorporate exclusions...')
				histogram(ax102, cell2mat(alllicktimes_unstimwrtc'), 30000, 'DisplayName', 'No Stimulation', 'Normalization', 'probability')
			    hold(ax102, 'on');
				histogram(ax102, cell2mat(alllicktimes_stimwrtc'), 30000, 'DisplayName', '+ Stimulation', 'Normalization', 'probability')
				yy = get(ax102, 'ylim');
				plot(ax102, [0, 0], yy,'r-', 'DisplayName', 'Cue')
			    plot(ax102, [considbound, considbound], yy,'b-', 'DisplayName', 'Consideration Limit')
				plot(ax102, [rb_ms/1000, rb_ms/1000], yy,'k-', 'DisplayName', 'Reward Boundary')
				plot(ax102, [obj.time.s.ITI_time_s, obj.time.s.ITI_time_s], yy,'k-', 'DisplayName', 'ITI Start')
				plot(ax102, [rxnwin_s, rxnwin_s], yy,'g-', 'DisplayName', 'Permitted Reaction Window')
				ylim(ax102, yy)
				title(ax102, ['Trials #', num2str(all_trials_this_block(1)), '-', num2str(all_trials_this_block(end))])
				xlim(ax102, [-5,18])
				ax102.XLabel.String = ['All Lick Times (s wrt ' obj.time.s.binnedLicks.ref ')'];
				ax102.YLabel.String = 'Percentage of Licks Across Session';
				ax = ax102;
			elseif strcmpi(Mode, 'r&h')
				ax11ONLYnostim = subplot(2, 2,  1, 'Parent', f);
				hold(ax11ONLYnostim, 'on')
				ax11ONLYstim = subplot(2, 2,  2, 'Parent', f);
				hold(ax11ONLYstim, 'on')
				ax102 = subplot(2, 2,  3, 'Parent', f);
				hold(ax102, 'on')
				ax1022 = subplot(2, 2,  4, 'Parent', f);
				hold(ax1022, 'on')
				set(ax11ONLYnostim,'fontsize',20);
				set(ax11ONLYstim,'fontsize',20);
				set(ax102,'fontsize',20);
				set(ax1022,'fontsize',20);
				% 
				% 	Final figure
				% 
				plot(ax11ONLYnostim, zeros(size(all_trials_this_block)), numel(nostim_trials),'r-', 'DisplayName', 'Cue')
			    plot(ax11ONLYnostim, [considbound, considbound], [1, numel(nostim_trials)],'b-', 'DisplayName', 'Consideration Limit')
				plot(ax11ONLYnostim, [rb_ms/1000, rb_ms/1000], [1, numel(nostim_trials)],'k-', 'DisplayName', 'Reward Boundary')
				plot(ax11ONLYnostim, [eot, eot], [1, numel(nostim_trials)],'k-', 'DisplayName', 'ITI Start')	
				plot(ax11ONLYstim, [-50,obj.time.ms.total_time_ms/1000], [numel(nostim_trials), numel(nostim_trials)], 'k-', 'DisplayName', '# no stim trials')				
				set(ax11ONLYnostim,  'YDir','reverse')
				title(ax11ONLYnostim, 'ONLY Consec Trials, No stim')
				xlim(ax11ONLYnostim, [0, obj.time.ms.total_time_ms/1000])
		        if numel(nostim_trials) > 1
		    		ylim(ax11ONLYnostim, [1, numel(nostim_trials)])
		        end
				ax11ONLYnostim.XLabel.String = ['Time (s wrt cue)'];
				ax11ONLYnostim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];
				plot(ax11ONLYstim, zeros(size(all_trials_this_block)), numel(stimTrials),'r-', 'DisplayName', 'Cue')
			    plot(ax11ONLYstim, [considbound, considbound], [1, numel(stimTrials)],'b-', 'DisplayName', 'Consideration Limit')
				plot(ax11ONLYstim, [rb_ms/1000, rb_ms/1000], [1, numel(stimTrials)],'k-', 'DisplayName', 'Reward Boundary')
				plot(ax11ONLYstim, [eot, eot], [1, numel(stimTrials)],'k-', 'DisplayName', 'ITI Start')
				plot(ax11ONLYstim, [-50,obj.time.ms.total_time_ms/1000], [numel(stimTrials), numel(stimTrials)], 'k-', 'DisplayName', '# stim trials')
						
				set(ax11ONLYstim,  'YDir','reverse')
				title(ax11ONLYstim, 'ONLY Consec Trials, + stim')
				xlim(ax11ONLYstim, [0, obj.time.ms.total_time_ms/1000])
		        if numel(stimTrials) > 1
		    		% ylim(ax11ONLYstim, [1, numel(stimTrials)])
		    		ylim(ax11ONLYstim, [1, numel(nostim_trials)])
		        end
				ax11ONLYstim.XLabel.String = ['Time (s wrt cue)'];
				ax11ONLYstim.YLabel.String = [obj.time.s.binnedLicks.ref 'Event #'];
			    unstimTrialTracker = 1;
			    stimTrialTracker = 1;
				for itrial = all_trials_this_block(1):all_trials_this_block(end)
					obj.valid_blocks.plotpnts = all_lick_times_ex_swrtc{itrial};
					if ismember(itrial, stimTrials)
						if ~isempty(obj.valid_blocks.plotpnts) && numel(f_lick_s_wrtcue.stim)>=itrial
							plot(ax11ONLYstim, obj.valid_blocks.plotpnts, stimTrialTracker.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
							plot(ax11ONLYstim, f_lick_s_wrtcue.stim(itrial), stimTrialTracker, 'mo', 'DisplayName', 'First Lick')
						end
						stimTrialTracker = stimTrialTracker +1;
					elseif ismember(itrial, nostim_trials) && numel(f_lick_s_wrtcue.unstim)>=itrial
						if ~isempty(obj.valid_blocks.plotpnts)
							plot(ax11ONLYnostim, obj.valid_blocks.plotpnts, unstimTrialTracker.*ones(numel(obj.valid_blocks.plotpnts), 1),'k.')
							plot(ax11ONLYnostim, f_lick_s_wrtcue.unstim(itrial), unstimTrialTracker, 'mo', 'DisplayName', 'First Lick')
						end
						unstimTrialTracker = unstimTrialTracker +1;
					else
						error('shouldn''t get here');
					end					
				end
				yy = get(ax11ONLYnostim, 'ylim');
				ylim(ax11ONLYstim, yy)
	            histogram(ax102, cell2mat(alllicktimes_unstimwrtc'), 30000, 'DisplayName', 'No Stimulation', 'displayStyle', 'stairs', 'LineWidth', 6)		    
				% histogram(ax102, cell2mat(alllicktimes_stimwrtc'), 30000, 'DisplayName', '+ Stimulation', 'Normalization', 'probability')
				yy = get(ax102, 'ylim');
				plot(ax102, [0, 0], yy,'r-', 'DisplayName', 'Cue')
			    plot(ax102, [considbound, considbound], yy,'b-', 'DisplayName', 'Consideration Limit')
				plot(ax102, [rb_ms/1000, rb_ms/1000], yy,'k-', 'DisplayName', 'Reward Boundary')
				plot(ax102, [obj.time.s.ITI_time_s, obj.time.s.ITI_time_s], yy,'k-', 'DisplayName', 'ITI Start')
				plot(ax102, [rxnwin_s, rxnwin_s], yy,'g-', 'DisplayName', 'Permitted Reaction Window')
				ylim(ax102, yy)
				ylim(ax102,[0, numel(alllicktimes_unstimwrtc)])
				title(ax102, ['Trials #', num2str(all_trials_this_block(1)), '-', num2str(all_trials_this_block(end))])
				xlim(ax102, [-5,18])
				ax102.XLabel.String = ['All Lick Times (s wrt ' obj.time.s.binnedLicks.ref ')'];
				ax102.YLabel.String = 'Licks/Trial';
				yticks(ax102, [0,round(0.5*numel(alllicktimes_unstimwrtc)),numel(alllicktimes_unstimwrtc)])
	            % histogram(ax1022, cell2mat(alllicktimes_unstimwrtc'), 30000, 'DisplayName', 'No Stimulation', 'Normalization', 'probability')		    
				histogram(ax1022, cell2mat(alllicktimes_stimwrtc'), 30000, 'DisplayName', '+ Stimulation', 'displayStyle', 'stairs', 'LineWidth', 6)
				% yy = get(ax1022, 'ylim');
				plot(ax1022, [0, 0], yy,'r-', 'DisplayName', 'Cue')
			    plot(ax1022, [considbound, considbound], yy,'b-', 'DisplayName', 'Consideration Limit')
				plot(ax1022, [rb_ms/1000, rb_ms/1000], yy,'k-', 'DisplayName', 'Reward Boundary')
				plot(ax1022, [obj.time.s.ITI_time_s, obj.time.s.ITI_time_s], yy,'k-', 'DisplayName', 'ITI Start')
				plot(ax1022, [rxnwin_s, rxnwin_s], yy,'g-', 'DisplayName', 'Permitted Reaction Window')
				% ylim(ax1022, yy)
				ylim(ax1022,[0, numel(alllicktimes_stimwrtc)])
				title(ax1022, ['Trials #', num2str(all_trials_this_block(1)), '-', num2str(all_trials_this_block(end))])
				xlim(ax1022, [-10,18])
				ax1022.XLabel.String = ['All Lick Times (s wrt ' obj.time.s.binnedLicks.ref ')'];
				ax1022.YLabel.String = 'Licks/Trial';
				yticks(ax1022, [0,round(0.5*numel(alllicktimes_stimwrtc)),numel(alllicktimes_stimwrtc)])
				linkaxes([ax11ONLYnostim,ax11ONLYstim,ax102,ax1022],'x')

				% linkaxes([ax102,ax1022],'y')
				obj.bootStat.nTrialsInBlockStim = numel(alllicktimes_stimwrtc);
				obj.bootStat.nTrialsInBlockUnstim = numel(alllicktimes_unstimwrtc);
				obj.rescaleYAxisStim(ax1022);
				obj.rescaleYAxisUnstim(ax102);
				ax = ax1022;
			elseif strcmpi(Mode, 'stimHxg')
				if nargin < 3 || isempty(ax)
					ax1022 = subplot(1, 1,  1, 'Parent', f);
					hold(ax1022, 'on')
					tt = ['Trials #', num2str(all_trials_this_block(1)), '-', num2str(all_trials_this_block(end)) ' | nstim: ' num2str(numel(alllicktimes_stimwrtc))];
				else
					error('You can''t actually overlay these!')
					ax1022 = ax;
					hold(ax1022, 'on')
					tt = get(get(ax, 'title'), 'String');
					tt = {tt, ['Trials #', num2str(all_trials_this_block(1)), '-', num2str(all_trials_this_block(end)) ' | nstim: ' num2str(numel(alllicktimes_stimwrtc))]};
				end
				yy = [0, round(0.15*numel(alllicktimes_stimwrtc))];
				
				
				set(ax1022,'fontsize',20);    
				histogram(ax1022, cell2mat(alllicktimes_stimwrtc'), 'binwidth', 0.25, 'DisplayName', '+ Stimulation', 'displayStyle', 'stairs', 'LineWidth', 6)
				plot(ax1022, [0, 0], yy,'r-', 'DisplayName', 'Stim On')
				plot(ax1022, [3, 3], yy,'r-', 'DisplayName', 'Stim Off')
				ylim(ax1022,yy);
				title(ax1022, tt)
				xlim(ax1022, [-10,18])
				ax1022.XLabel.String = ['time (s)'];
				ax1022.YLabel.String = 'Lick Rate (Licks/Trial)';
				yticks(ax1022, [yy(1),yy(2)])
				obj.bootStat.nTrialsInBlockStim = numel(alllicktimes_stimwrtc);
				obj.bootStat.ylimStim = yy;
				obj.bootStat.nTrialsInBlockUnstim = numel(alllicktimes_unstimwrtc);
				obj.rescaleYAxisStim(ax1022);
				ax = ax1022;
			end
    	end
    	function [nstim, nunstim] = getNTrialsInCollatedPlot(obj)
    		nstim = obj.bootStat.nTrialsInBlockStim;
    		nunstim = obj.bootStat.nTrialsInBlockUnstim;
		end
		function rescaleYAxisHistogram(obj, c, ax)
			if nargin < 3
				ax = gca;
			end
			% 
			% 	Adjusts Y-tick labels on histogram by a constant
			% 
			actualY = get(ax, 'ytick');
			scaledY = round(actualY/c,1);
			set(ax, 'yticklabel', num2cell(scaledY))
		end
		function rescaleYAxisStim(obj, ax)
			if nargin < 2
				ax = gca;
			end
			yticks(ax, obj.bootStat.ylimStim)
			[c,~] = obj.getNTrialsInCollatedPlot;
			obj.rescaleYAxisHistogram(c, ax);
		end
		function rescaleYAxisUnstim(obj, ax)
			if nargin < 2
				ax = gca;
			end
			[~,c] = obj.getNTrialsInCollatedPlot;
			obj.rescaleYAxisHistogram(c, ax);
		end
		function ax = plotCollatedNTJ(obj, Mode, ax, IDs)
			% 
	    	% 	Use on collated naive test + juice datasets to plot the histogram
	    	% 
	    	% 	Mode: 	'stimHxg'     -- plots only the stimulated hxg in units of licks per trial NO PRE JUICE!
	    	% 			'raster4' 	  -- plots the raster for all trials as per usual --- NB, only plots one ID, so if you put in more than one ID, will only use the 1st one!
	    	% 
	    	if nargin < 2 || isempty(Mode)
	    		Mode = 'stimHxg';
    		end
    		if nargin < 4 || isempty(IDs)
    			IDs = 1:numel({obj.bootStat.collateResult.ID});
			end
			
			
			ntrialsNoPreJuice = 0;
			ntrials = 0;
			% ALL = [];
			all_licks_swrtSO_ex_s.NOpreStimJuice = [];

			for iID = IDs
				ALL = obj.bootStat.collateResult(iID).all_licks_swrtSO_ex_s.all;
				ntrials = ntrials + numel(ALL);
				ntrialsNoPreJuice = ntrialsNoPreJuice + numel(ALL(obj.bootStat.collateResult(iID).valid_blocks.NOpreStimJuiceTrials));
				NOpreStimJuice = obj.bootStat.collateResult(iID).all_licks_swrtSO_ex_s.NOpreStimJuice;
				% all_licks_swrtSO_ex_s.all = [all_licks_swrtSO_ex_s(iID).all;ALL];
				all_licks_swrtSO_ex_s.NOpreStimJuice = [all_licks_swrtSO_ex_s.NOpreStimJuice; NOpreStimJuice];
			end
			% 
			if nargin < 3 || isempty(ax)
				f = figure;
				set(gcf, 'color', 'white');	
				if numel(IDs) == 1
					set(f, 'name', obj.bootStat.collateResult(IDs).ID)
				end
			end


			if strcmpi(Mode,'stimHxg')
				axNOPreStimJuiceHXG = subplot(1, 1,  1, 'Parent', f);
				hold(axNOPreStimJuiceHXG, 'on')
				set(axNOPreStimJuiceHXG,'fontsize',20);    				
				yy = [0, round(0.5*ntrialsNoPreJuice)];
				tt = ['Trials #1-', num2str(ntrials) ' | nNoPreJuice: ' num2str(ntrialsNoPreJuice)];
				% 	axNOPreStimJuiceHXG
				histogram(axNOPreStimJuiceHXG, all_licks_swrtSO_ex_s.NOpreStimJuice(~isnan(all_licks_swrtSO_ex_s.NOpreStimJuice)), 'binwidth', 0.25,'displayName', 'No Pre-Juice','displaystyle','stairs', 'LineWidth',6);
				plot(axNOPreStimJuiceHXG, [0, 0], yy,'r-', 'DisplayName', 'Stim On')
				plot(axNOPreStimJuiceHXG, [3, 3], yy,'r-', 'DisplayName', 'Stim Off')
				axNOPreStimJuiceHXG.YLabel.String = 'Lick Rate (Licks/Trial)';	
				axNOPreStimJuiceHXG.XLabel.String = 'time (s)';
				xlim(axNOPreStimJuiceHXG, [-10, 18])
				ylim(axNOPreStimJuiceHXG,yy);
				title(axNOPreStimJuiceHXG, tt);
				
				ax = axNOPreStimJuiceHXG;

				yticks(axNOPreStimJuiceHXG, [yy(1),yy(2)])
				obj.bootStat.nTrialsInBlockStim = ntrialsNoPreJuice;
				obj.bootStat.ylimStim = yy;
				obj.bootStat.nTrialsInBlockUnstim = 0;
				obj.rescaleYAxisStim(axNOPreStimJuiceHXG);
			elseif strcmpi(Mode, 'raster4')
				if numel(IDs) > 1
					error('Only input 1 ID at a time for use with raster4')
				end
				ax = subplot(1, 1,  1, 'Parent', f);
				hold(ax, 'on')
				nEvents = numel(obj.bootStat.collateResult(IDs).all_flicks_swrtSO_ex_s.all);
				% 
				% 	Plot raster of all licks with first licks overlaid
				% 
				plot(ax, [0,0], [1,nEvents],'k-', 'DisplayName', 'Stim On')
				plot(ax, [3, 3], [1,nEvents],'k-', 'DisplayName', 'Stim Off')

				set(ax,  'YDir','reverse')
				title(ax, ['Lick Raster Aligned to Stim On'])
				xlim(ax, [-10, 18])
		        if nEvents > 1
		    		ylim(ax, [1, nEvents])
		        end
				ax.XLabel.String = ['time (s)'];
				ax.YLabel.String = ['Trial #'];


			    lick_times = obj.bootStat.collateResult(IDs).all_licks_swrtSO_ex_s.all;				    
			    % 
				plot(ax, obj.bootStat.collateResult(IDs).all_flicks_swrtSO_ex_s.all, 1:numel(obj.bootStat.collateResult(IDs).all_flicks_swrtSO_ex_s.all), 'mo', 'DisplayName', 'First Lick', 'markerfacecolor', 'm')
				%
				
				for itrial = 1:nEvents
					plotpnts_lick = lick_times{itrial};
					if ismember(itrial, obj.bootStat.collateResult(IDs).valid_blocks.preStimJuiceTrials)
						cc = 'r.';
					else
						cc = 'k.';
					end
					if ~isempty(plotpnts_lick)
						plot(ax, plotpnts_lick, itrial.*ones(numel(plotpnts_lick), 1),cc)
					end	
				end
				% nEvents = numel(obj.time.s.binnedLicks.refevents);
				event_times = obj.bootStat.collateResult(IDs).juice_s_bt_wrtref;
				for itrial = 1:nEvents			
					plotpnts = event_times{itrial};
					if ~isempty(plotpnts)
						plot(ax, plotpnts, itrial.*ones(numel(plotpnts), 1),'g.', 'markersize', 20)
					end	
				end
				% event_times = obj.bootStat.collateResult(IDs).juice_s_bt_wrtref;
				% for itrial = 1:numel(event_times)
				% 	plotpnts = event_times{itrial};
				% 	if ~isempty(plotpnts)
				% 		plot(ax, plotpnts, obj.bootStat.collateResult(IDs).valid_blocks.preStimJuiceTrials(itrial).*ones(numel(plotpnts), 1),'g.', 'markersize', 20)
				% 	end	
				% end
			end
		end
		%------------------------------------------------------------------------------------------------------------------- 
		% 	SINGLE SESSION STAT METHODS
		%-------------------------------------------------------------------------------------------------------------------
		% 
		% 	Legacy stats pre-processing
		% 
		function getvalidblocks(obj,Mode)
			% 
			% 	Mode = either full-day or number of consecutive valid trials per block
			% 
			% 	n.b.! The legacy version only used consecutive trials even for full-day processing. I'll keep this for now.
			% 
			% 	flicks always ref'd to cue
			% 
			considerationBoundary_s = obj.bootStat.considerationBoundary_ms/1000;
			if nargin < 2
				Mode = 'full-day';
			end
			obj.valid_blocks.Mode = Mode;
			if strcmpi(Mode, 'full-day') || strcmpi(Mode, 'fullday')
				obj.valid_blocks.num_trials_per_block = obj.iv.init_variables_ChR2.num_trials;
			else
				obj.valid_blocks.num_trials_per_block = Mode;
			end
			obj.valid_blocks.flstim_s = nan(size(obj.lick.swrtc.stim_first_licks_swrtc));
			obj.valid_blocks.flunstim_s = nan(size(obj.lick.swrtc.unstim_first_licks_swrtc));
			obj.valid_blocks.flstim_s(obj.lick.swrtc.stim_first_licks_swrtc>considerationBoundary_s) = obj.lick.swrtc.stim_first_licks_swrtc(obj.lick.swrtc.stim_first_licks_swrtc>considerationBoundary_s);
			obj.valid_blocks.flunstim_s(obj.lick.swrtc.unstim_first_licks_swrtc>considerationBoundary_s) = obj.lick.swrtc.unstim_first_licks_swrtc(obj.lick.swrtc.unstim_first_licks_swrtc>considerationBoundary_s);

			obj.valid_blocks.allfl_s = nansum([obj.valid_blocks.flstim_s, obj.valid_blocks.flunstim_s],2); 
			obj.valid_blocks.allflconsec_idx = find(obj.valid_blocks.allfl_s > 0);
			obj.valid_blocks.allflconsec_s = obj.valid_blocks.allfl_s(obj.valid_blocks.allflconsec_idx);

			obj.valid_blocks.nconsec = length(obj.valid_blocks.allflconsec_s);
			obj.valid_blocks.ndivs = floor(obj.valid_blocks.nconsec/obj.valid_blocks.num_trials_per_block);
			if obj.valid_blocks.ndivs ~= 0
			    obj.valid_blocks.trimindicies = [1, 1 + obj.valid_blocks.num_trials_per_block*(1:obj.valid_blocks.ndivs)];
			elseif ~strcmp(obj.iv.stimtype, 'NaiveTest-No Juice') && ~strcmp(obj.iv.stimtype, 'NaiveTest+Juice')
			    obj.valid_blocks.ndivs = 1;
			    obj.valid_blocks.trimindicies = [1, obj.valid_blocks.nconsec];
			else % for naive test, include everything.
				obj.valid_blocks.ndivs = 1;
			    obj.valid_blocks.trimindicies = 1; %obj.valid_blocks.num_trials_per_block]
			end
			% 
			% 	Break into valid consecutive blocks
			% 
			for i_div = 1:length(obj.valid_blocks.trimindicies)
				% 
				% 	Find all the stim indicies and unstim indicies in this block
				% 
				if ~(length(obj.valid_blocks.trimindicies)==1) && obj.valid_blocks.trimindicies(i_div) == obj.valid_blocks.trimindicies(end)
					break
			    elseif obj.valid_blocks.trimindicies(end) == 0 || strcmp(obj.iv.stimtype, 'NaiveTest-No Juice') || strcmp(obj.iv.stimtype, 'NaiveTest+Juice')% case with no licks at all
			        obj.valid_blocks.f_lick_s_wrtcue(i_div).ntrials_this_block = numel(obj.valid_blocks.num_trials_per_block);
			    elseif obj.valid_blocks.trimindicies(i_div+1) == length(obj.valid_blocks.allflconsec_idx) + 1 % added for the unlikely but possible case the trim indicies are off by 1 for last block
			        obj.valid_blocks.f_lick_s_wrtcue(i_div).ntrials_this_block = obj.valid_blocks.allflconsec_idx(obj.valid_blocks.trimindicies(i_div+1)-1) - obj.valid_blocks.allflconsec_idx(obj.valid_blocks.trimindicies(i_div));
				else
					obj.valid_blocks.f_lick_s_wrtcue(i_div).ntrials_this_block = obj.valid_blocks.allflconsec_idx(obj.valid_blocks.trimindicies(i_div+1)) - obj.valid_blocks.allflconsec_idx(obj.valid_blocks.trimindicies(i_div));
				end
				% 
				% 
				if ~strcmp(obj.iv.stimtype, 'NaiveTest-No Juice') && ~strcmp(obj.iv.stimtype, 'NaiveTest+Juice')
			 		obj.valid_blocks.all_trials_this_block = obj.valid_blocks.allflconsec_idx(obj.valid_blocks.trimindicies(i_div)) : obj.valid_blocks.allflconsec_idx(obj.valid_blocks.trimindicies(i_div+1) - 1)';
				else
					obj.valid_blocks.all_trials_this_block = 1:obj.valid_blocks.num_trials_per_block;
				end
			    obj.valid_blocks.f_lick_s_wrtcue(i_div).stimTrials = find(obj.valid_blocks.flstim_s(obj.valid_blocks.all_trials_this_block)>0);
				obj.valid_blocks.f_lick_s_wrtcue(i_div).unstimTrials = find(obj.valid_blocks.flunstim_s(obj.valid_blocks.all_trials_this_block)>0);
				% 
				% 	Get the lick times for each trial in block
				% 
				obj.valid_blocks.f_lick_s_wrtcue(i_div).stim = obj.valid_blocks.flstim_s(obj.valid_blocks.all_trials_this_block);
				obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim = obj.valid_blocks.flunstim_s(obj.valid_blocks.all_trials_this_block);
				% 
				% 	Get the frequency of licks per trial
				% 
				obj.valid_blocks.trial_nums(i_div).stim = find(obj.lick.swrtc.stim_first_licks_swrtc > 0);
				obj.valid_blocks.trial_nums(i_div).stim = obj.valid_blocks.trial_nums(i_div).stim(obj.valid_blocks.trial_nums(i_div).stim >= obj.valid_blocks.trimindicies(i_div) & obj.valid_blocks.trial_nums(i_div).stim < obj.valid_blocks.trimindicies(i_div) + obj.valid_blocks.num_trials_per_block - 1);
				obj.valid_blocks.trial_nums(i_div).unstim = find(obj.lick.swrtc.unstim_first_licks_swrtc > 0);
				obj.valid_blocks.trial_nums(i_div).unstim = obj.valid_blocks.trial_nums(i_div).unstim(obj.valid_blocks.trial_nums(i_div).unstim >= obj.valid_blocks.trimindicies(i_div) & obj.valid_blocks.trial_nums(i_div).unstim < obj.valid_blocks.trimindicies(i_div) + obj.valid_blocks.num_trials_per_block - 1);
				% 
				% 	n.b., lick.legacy.lick_times_by_trial is okay for this. It doesn't remove any licks between 700-1000ms
				% 
				obj.valid_blocks.all_licks(i_div).nlicks_stim = sum((obj.lick.legacy.s.lick_ex_times_by_trial_s(obj.valid_blocks.trial_nums(i_div).stim, :)>0), 2);
				obj.valid_blocks.all_licks(i_div).nlicks_unstim = sum((obj.lick.legacy.s.lick_ex_times_by_trial_s(obj.valid_blocks.trial_nums(i_div).unstim, :)>0), 2);
				obj.valid_blocks.all_licks(i_div).mean_nlicks_per_trial_stim = sum(obj.valid_blocks.all_licks(i_div).nlicks_stim)/length(obj.valid_blocks.all_licks(i_div).nlicks_stim);
				obj.valid_blocks.all_licks(i_div).mean_nlicks_per_trial_unstim = sum(obj.valid_blocks.all_licks(i_div).nlicks_unstim)/length(obj.valid_blocks.all_licks(i_div).nlicks_unstim);
				% 
				%	Calculate the median, cdf, and del-cdf 
				% 
				% obj.valid_blocks.endoftrial_s = eot * fs_ms / 1000; % changed 10/18/19
				obj.valid_blocks.endoftrial_s = obj.time.s.ITI_time_s;
				% 
				% 	Get medians
				% 
				obj.valid_blocks.Median(i_div).stim = nanmedian(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim > 0));
				obj.valid_blocks.Median(i_div).unstim = nanmedian(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim > 0));
				% 
				% 	Get eCDFs for stim and unstim
				% 
			    unstim_ecdf_pntr = obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim(obj.valid_blocks.f_lick_s_wrtcue(i_div).unstim > 0);
			    stim_ecdf_pntr = obj.valid_blocks.f_lick_s_wrtcue(i_div).stim(obj.valid_blocks.f_lick_s_wrtcue(i_div).stim > 0);
			    if ~isempty(unstim_ecdf_pntr)
			    	[obj.valid_blocks.eCDF(i_div).unstim_f, obj.valid_blocks.eCDF(i_div).unstim_x] = ecdf(unstim_ecdf_pntr);
			    else
			        obj.valid_blocks.eCDF(i_div).unstim_f = [];
			        obj.valid_blocks.eCDF(i_div).unstim_x = [];
			    end
				if ~isempty(stim_ecdf_pntr)
			        [obj.valid_blocks.eCDF(i_div).stim_f, obj.valid_blocks.eCDF(i_div).stim_x] = ecdf(stim_ecdf_pntr);
			    else
			        obj.valid_blocks.eCDF(i_div).stim_f = [];
			        obj.valid_blocks.eCDF(i_div).stim_x = [];
			    end
				% 
				% 	Make sure the length of the two cdfs is the same by appending the peak cdf (eg, 1) to the position of the longer x array
				% 
			    if ~isempty(unstim_ecdf_pntr) && ~isempty(stim_ecdf_pntr)
			        if obj.valid_blocks.eCDF(i_div).unstim_x(end) < obj.valid_blocks.eCDF(i_div).stim_x(end)
			            obj.valid_blocks.trimunstim = true;
			            obj.valid_blocks.eCDF(i_div).unstim_x(end+1) = obj.valid_blocks.eCDF(i_div).stim_x(end);
			            obj.valid_blocks.eCDF(i_div).unstim_f(end+1) = max(obj.valid_blocks.eCDF(i_div).unstim_f);
			        else
			            obj.valid_blocks.trimunstim = false;
			            obj.valid_blocks.eCDF(i_div).stim_x(end+1) = obj.valid_blocks.eCDF(i_div).unstim_x(end);
			            obj.valid_blocks.eCDF(i_div).stim_f(end+1) = max(obj.valid_blocks.eCDF(i_div).stim_f);
			        end
			        % 
					%  Standard AUC
					% 
					obj.valid_blocks.eCDF(i_div).stim_AUC = trapz(obj.valid_blocks.eCDF(i_div).stim_x, obj.valid_blocks.eCDF(i_div).stim_f);
					obj.valid_blocks.eCDF(i_div).unstim_AUC = trapz(obj.valid_blocks.eCDF(i_div).unstim_x, obj.valid_blocks.eCDF(i_div).unstim_f);
					% 
					%  AUC cut off at end of trial (7 s) -- add a point with max cdf at 7 to make sure area calc is good.
					% 
					obj.valid_blocks.eCDF(i_div).stim_x_EOT = obj.valid_blocks.eCDF(i_div).stim_x(obj.valid_blocks.eCDF(i_div).stim_x <= obj.valid_blocks.endoftrial_s);
					obj.valid_blocks.eCDF(i_div).stim_f_EOT = obj.valid_blocks.eCDF(i_div).stim_f(obj.valid_blocks.eCDF(i_div).stim_x <= obj.valid_blocks.endoftrial_s);
					obj.valid_blocks.eCDF(i_div).stim_x_EOT(end+1) = obj.valid_blocks.endoftrial_s;
                    if ~isempty(max(obj.valid_blocks.eCDF(i_div).stim_f_EOT))
    					obj.valid_blocks.eCDF(i_div).stim_f_EOT(end+1) = max(obj.valid_blocks.eCDF(i_div).stim_f_EOT);
                    else
                        obj.valid_blocks.eCDF(i_div).stim_f_EOT(end+1) = 1;
                    end

					obj.valid_blocks.eCDF(i_div).unstim_x_EOT = obj.valid_blocks.eCDF(i_div).unstim_x(obj.valid_blocks.eCDF(i_div).unstim_x <= obj.valid_blocks.endoftrial_s);
					obj.valid_blocks.eCDF(i_div).unstim_f_EOT = obj.valid_blocks.eCDF(i_div).unstim_f(obj.valid_blocks.eCDF(i_div).unstim_x <= obj.valid_blocks.endoftrial_s);
					obj.valid_blocks.eCDF(i_div).unstim_x_EOT(end+1) = obj.valid_blocks.endoftrial_s;
                    if ~isempty(max(obj.valid_blocks.eCDF(i_div).unstim_f_EOT))
    					obj.valid_blocks.eCDF(i_div).unstim_f_EOT(end+1) = max(obj.valid_blocks.eCDF(i_div).unstim_f_EOT);
                    else
                        obj.valid_blocks.eCDF(i_div).unstim_f_EOT(end+1) = 1;
                    end

					obj.valid_blocks.eCDF(i_div).stim_AUC_EOT = trapz(obj.valid_blocks.eCDF(i_div).stim_x_EOT, obj.valid_blocks.eCDF(i_div).stim_f_EOT);
					obj.valid_blocks.eCDF(i_div).unstim_AUC_EOT = trapz(obj.valid_blocks.eCDF(i_div).unstim_x_EOT, obj.valid_blocks.eCDF(i_div).unstim_f_EOT);
					% 
					%  AUC, include the trials past 7 s, cdf=cdfmax
					%  
					% 		trim off the extra value
					% 
					if obj.valid_blocks.trimunstim == true
						obj.valid_blocks.eCDF(i_div).stim_x_EOTinclusive = obj.valid_blocks.eCDF(i_div).stim_x;
						obj.valid_blocks.eCDF(i_div).stim_f_EOTinclusive = obj.valid_blocks.eCDF(i_div).stim_f;		
						obj.valid_blocks.eCDF(i_div).stim_x_EOTinclusive(obj.valid_blocks.eCDF(i_div).stim_x >= obj.valid_blocks.endoftrial_s) = obj.valid_blocks.endoftrial_s;
						obj.valid_blocks.eCDF(i_div).stim_f_EOTinclusive(obj.valid_blocks.eCDF(i_div).stim_x >= obj.valid_blocks.endoftrial_s) = max(obj.valid_blocks.eCDF(i_div).stim_f);
						% 
						obj.valid_blocks.eCDF(i_div).stim_AUC_EOTinclusive = trapz(obj.valid_blocks.eCDF(i_div).stim_x_EOTinclusive, obj.valid_blocks.eCDF(i_div).stim_f_EOTinclusive);
						% 
						% 
						obj.valid_blocks.eCDF(i_div).unstim_x_EOTinclusive = obj.valid_blocks.eCDF(i_div).unstim_x(1:end-1);
						obj.valid_blocks.eCDF(i_div).unstim_f_EOTinclusive = obj.valid_blocks.eCDF(i_div).unstim_f(1:end-1);		
						obj.valid_blocks.eCDF(i_div).unstim_x_EOTinclusive(obj.valid_blocks.eCDF(i_div).unstim_x >= obj.valid_blocks.endoftrial_s) = obj.valid_blocks.endoftrial_s;
						obj.valid_blocks.eCDF(i_div).unstim_f_EOTinclusive(obj.valid_blocks.eCDF(i_div).unstim_x >= obj.valid_blocks.endoftrial_s) = max(obj.valid_blocks.eCDF(i_div).unstim_f(1:end-1));
						% 
						obj.valid_blocks.eCDF(i_div).unstim_AUC_EOTinclusive = trapz(obj.valid_blocks.eCDF(i_div).unstim_x_EOTinclusive, obj.valid_blocks.eCDF(i_div).unstim_f_EOTinclusive);
					else
						obj.valid_blocks.eCDF(i_div).stim_x_EOTinclusive = obj.valid_blocks.eCDF(i_div).stim_x(1:end-1);
						obj.valid_blocks.eCDF(i_div).stim_f_EOTinclusive = obj.valid_blocks.eCDF(i_div).stim_f(1:end-1);		
						obj.valid_blocks.eCDF(i_div).stim_x_EOTinclusive(obj.valid_blocks.eCDF(i_div).stim_x >= obj.valid_blocks.endoftrial_s) = obj.valid_blocks.endoftrial_s;
						obj.valid_blocks.eCDF(i_div).stim_f_EOTinclusive(obj.valid_blocks.eCDF(i_div).stim_x >= obj.valid_blocks.endoftrial_s) = max(obj.valid_blocks.eCDF(i_div).stim_f(1:end-1));
						% 
						obj.valid_blocks.eCDF(i_div).stim_AUC_EOTinclusive = trapz(obj.valid_blocks.eCDF(i_div).stim_x_EOTinclusive, obj.valid_blocks.eCDF(i_div).stim_f_EOTinclusive);
						% 
						% 
						obj.valid_blocks.eCDF(i_div).unstim_x_EOTinclusive = obj.valid_blocks.eCDF(i_div).unstim_x;
						obj.valid_blocks.eCDF(i_div).unstim_f_EOTinclusive = obj.valid_blocks.eCDF(i_div).unstim_f;		
						obj.valid_blocks.eCDF(i_div).unstim_x_EOTinclusive(obj.valid_blocks.eCDF(i_div).unstim_x >= obj.valid_blocks.endoftrial_s) = obj.valid_blocks.endoftrial_s;
						obj.valid_blocks.eCDF(i_div).unstim_f_EOTinclusive(obj.valid_blocks.eCDF(i_div).unstim_x >= obj.valid_blocks.endoftrial_s) = max(obj.valid_blocks.eCDF(i_div).unstim_f);
						% 
						obj.valid_blocks.eCDF(i_div).unstim_AUC_EOTinclusive = trapz(obj.valid_blocks.eCDF(i_div).unstim_x_EOTinclusive, obj.valid_blocks.eCDF(i_div).unstim_f_EOTinclusive);
					end
					% 
					% CDF difference measurements:
					% 
					obj.valid_blocks.eCDF(i_div).delAUCun_minus_stim = obj.valid_blocks.eCDF(i_div).unstim_AUC - obj.valid_blocks.eCDF(i_div).stim_AUC;
					obj.valid_blocks.eCDF(i_div).delAUC_EOTun_minus_stim = obj.valid_blocks.eCDF(i_div).unstim_AUC_EOT - obj.valid_blocks.eCDF(i_div).stim_AUC_EOT;
					obj.valid_blocks.eCDF(i_div).delAUC_EOTinclusive_un_minus_stim = obj.valid_blocks.eCDF(i_div).unstim_AUC_EOTinclusive - obj.valid_blocks.eCDF(i_div).stim_AUC_EOTinclusive;
					% 
					% 	Standard AUC plot
					% 
			    else
% 			        warning('One of the ecdf''s is empty, so we are not bothering permuting')
			        obj.valid_blocks.trimunstim = false;
			        obj.valid_blocks.eCDF(i_div).unstim_f = nan;
			        obj.valid_blocks.eCDF(i_div).unstim_x = nan;
			        obj.valid_blocks.eCDF(i_div).stim_f = nan;
			        obj.valid_blocks.eCDF(i_div).stim_x = nan;
			        obj.valid_blocks.eCDF(i_div).stim_AUC = nan;
					obj.valid_blocks.eCDF(i_div).unstim_AUC = nan;
					obj.valid_blocks.eCDF(i_div).stim_x_EOT = nan;
					obj.valid_blocks.eCDF(i_div).stim_f_EOT = nan;
					obj.valid_blocks.eCDF(i_div).stim_x_EOT(end+1) = nan;
					obj.valid_blocks.eCDF(i_div).stim_f_EOT(end+1) = nan;
					obj.valid_blocks.eCDF(i_div).unstim_x_EOT = nan;
					obj.valid_blocks.eCDF(i_div).unstim_f_EOT = nan;
					obj.valid_blocks.eCDF(i_div).unstim_x_EOT(end+1) = nan;
					obj.valid_blocks.eCDF(i_div).unstim_f_EOT(end+1) = nan;
					obj.valid_blocks.eCDF(i_div).stim_AUC_EOT = nan;
					obj.valid_blocks.eCDF(i_div).unstim_AUC_EOT = nan;
					obj.valid_blocks.eCDF(i_div).stim_x_EOTinclusive = nan;
					obj.valid_blocks.eCDF(i_div).stim_f_EOTinclusive = nan;
					obj.valid_blocks.eCDF(i_div).stim_x_EOTinclusive(obj.valid_blocks.eCDF(i_div).stim_x >= obj.valid_blocks.endoftrial_s) = nan;
					obj.valid_blocks.eCDF(i_div).stim_f_EOTinclusive(obj.valid_blocks.eCDF(i_div).stim_x >= obj.valid_blocks.endoftrial_s) = nan;
					obj.valid_blocks.eCDF(i_div).stim_AUC_EOTinclusive = nan;
					obj.valid_blocks.eCDF(i_div).unstim_x_EOTinclusive = nan;
					obj.valid_blocks.eCDF(i_div).unstim_f_EOTinclusive = nan;
					obj.valid_blocks.eCDF(i_div).unstim_x_EOTinclusive(obj.valid_blocks.eCDF(i_div).unstim_x >= obj.valid_blocks.endoftrial_s) = nan;
					obj.valid_blocks.eCDF(i_div).unstim_f_EOTinclusive(obj.valid_blocks.eCDF(i_div).unstim_x >= obj.valid_blocks.endoftrial_s) = nan;
					obj.valid_blocks.eCDF(i_div).unstim_AUC_EOTinclusive = nan;
					% 
					% CDF difference measurements:
					% 
					obj.valid_blocks.eCDF(i_div).delAUCun_minus_stim = nan;
					obj.valid_blocks.eCDF(i_div).delAUC_EOTun_minus_stim = nan;
					obj.valid_blocks.eCDF(i_div).delAUC_EOTinclusive_un_minus_stim = nan;
				end
		    
				% ---------------------------------------------- IF MOVEMENT SIGNALS ALSO PRESENT, CALCULATE VARIANCES ---------------------------------
				% 
				% 
				if ~strcmp(obj.iv.stimtype, 'NaiveTest-No Juice') && ~strcmp(obj.iv.stimtype, 'NaiveTest+Juice')
					% 
					% 	Gather the trial # indicies for stim and unstim case (i.e., the actual trial number FOR VALID TRIALS (not reactions! 6/29/18))
					% 
					obj.valid_blocks.m.trial_nums(i_div).stim = find(obj.lick.swrtc.stim_first_licks_swrtc > 0);
					obj.valid_blocks.m.trial_nums(i_div).stim = obj.valid_blocks.m.trial_nums(i_div).stim(obj.valid_blocks.m.trial_nums(i_div).stim >= obj.valid_blocks.all_trials_this_block(1) & obj.valid_blocks.m.trial_nums(i_div).stim <= obj.valid_blocks.all_trials_this_block(end));
					obj.valid_blocks.m.trial_nums(i_div).unstim = find(obj.lick.swrtc.unstim_first_licks_swrtc > 0);
					obj.valid_blocks.m.trial_nums(i_div).unstim = obj.valid_blocks.m.trial_nums(i_div).unstim(obj.valid_blocks.m.trial_nums(i_div).unstim >= obj.valid_blocks.all_trials_this_block(1) & obj.valid_blocks.m.trial_nums(i_div).unstim <= obj.valid_blocks.all_trials_this_block(end));
		        else
					disp('We are using all trials for naive test')
					obj.valid_blocks.m.trial_nums(i_div).stim = find(~isnan(obj.lick.swrtc.stim_first_licks_swrtc));
					obj.valid_blocks.m.trial_nums(i_div).unstim = find(~isnan(obj.lick.swrtc.unstim_first_licks_swrtc));			
				end
				%--------------------------------------------------
				%		METHOD B: ROBIN BLISS METHOD: permute trials, not timepoints for variance comparison (6/25/18)
				% 		B. Take all timepoints up to the first movement - from cue to the lick
				% 
				%	data = {1 x nblocks} cell
				% 	data{1, iblock} = 1x2 cell (stim, nostim)
				% 	data{1, iblock}{1, 1=stim or 2=nostim} = ntrials x 1 vector OR CELL
				%
				obj.valid_blocks.m.cue_pos = obj.iv.init_variables_X.time_parameters.first_post_cue_position;
				if 	~strcmp(obj.iv.stimtype, 'NaiveTest-No Juice') && ~strcmp(obj.iv.stimtype, 'NaiveTest+Juice')
	                obj.valid_blocks.m.stim_flick_mpos = obj.h.stim_first_licks(obj.valid_blocks.m.trial_nums(i_div).stim)/obj.iv.samples_per_ms_ChR2 * obj.iv.samples_per_ms_X;
					obj.valid_blocks.m.unstim_flick_mpos = obj.h.unstim_first_licks(obj.valid_blocks.m.trial_nums(i_div).unstim)/obj.iv.samples_per_ms_ChR2 * obj.iv.samples_per_ms_X;
			% 				obj.valid_blocks.m.stim_flick_mpos = obj.lick.pos_ChR2.stim_first_licks_posChR2(obj.valid_blocks.m.trial_nums(i_div).stim)/samples_per_ms_ChR2 * samples_per_ms_X;
			% 				obj.valid_blocks.m.unstim_flick_mpos = obj.lick.pos_ChR2.unstim_first_licks_posChR2_posChR2(obj.valid_blocks.m.trial_nums(i_div).unstim)/samples_per_ms_ChR2 * samples_per_ms_X;
				% 
				% 	If we are running naive test, we want all the trial duration included in the measurement...
				% 
				elseif strcmp(obj.iv.stimtype, 'NaiveTest-No Juice')
					disp('	...and we are looking at the entire trial duration (cue to end of ITI)');
					obj.valid_blocks.m.unstim_flick_mpos = size(obj.m.smooth_vbt, 2).*ones(numel(obj.valid_blocks.m.trial_nums(i_div).unstim), 1);
					obj.valid_blocks.m.stim_flick_mpos = size(obj.m.smooth_vbt, 2).*ones(numel(obj.valid_blocks.m.trial_nums(i_div).stim), 1);
				elseif strcmp(obj.iv.stimtype, 'NaiveTest+Juice')
					disp('	...and we are looking at the entire trial duration (cue to end of ITI)');
					error('not implemented')
					obj.valid_blocks.m.unstim_flick_mpos = size(obj.m.smooth_vbt, 2).*ones(numel(obj.valid_blocks.m.trial_nums(i_div).unstim), 1);
					obj.valid_blocks.m.stim_flick_mpos = size(obj.m.smooth_vbt, 2).*ones(numel(obj.valid_blocks.m.trial_nums(i_div).stim), 1);
				end
				% 
				% 	Initialize the data cell array {1, i_div}{1,1} = stim, {1, i_div}{1,2} = unstim
				% 
				obj.valid_blocks.m.data{1, i_div} = cell(1,2); 
				% 
				obj.valid_blocks.m.data{1, i_div}{1,1} = cell(length(obj.valid_blocks.m.trial_nums(i_div).stim),1); % for the stimulated case
				obj.valid_blocks.m.data{1, i_div}{1,2} = cell(length(obj.valid_blocks.m.trial_nums(i_div).unstim),1); % for the stimulated case
				% 
				% 	For each stimulated trial, add its first sec to the cell array of data
				% 		
			% 			warning('Possible Problem? -- X data for trials with no-lick is entered as nans to the concat, and the whole cue-to-lick interval is included')
				for i_stim_trial = 1:length(obj.valid_blocks.m.trial_nums(i_div).stim)
					% 
					% 	If the first lick occurred within the post-stim window, ignore and leave timestamp vector as []
					% 		nanvar([]) = Nan, so we can use this to just ignore these trials
					% 
					obj.valid_blocks.m.data{1, i_div}{1,1}{i_stim_trial, 1} = obj.m.bandPass_vbt(obj.valid_blocks.m.trial_nums(i_div).stim(i_stim_trial), obj.valid_blocks.m.cue_pos:round(obj.valid_blocks.m.stim_flick_mpos(i_stim_trial)));
				end
				% 
				% 	For each unstimulated trial, add its first sec to the cell array of data
				% 
				for i_unstim_trial = 1:length(obj.valid_blocks.m.trial_nums(i_div).unstim)
					% 
					% 	If the first lick occurred within the post-stim window, ignore and leave timestamp vector as []
					% 		nanvar([]) = Nan, so we can use this to just ignore these trials
					% 
					obj.valid_blocks.m.data{1, i_div}{1,2}{i_unstim_trial, 1} = obj.m.bandPass_vbt(obj.valid_blocks.m.trial_nums(i_div).unstim(i_unstim_trial), obj.valid_blocks.m.cue_pos:round(obj.valid_blocks.m.unstim_flick_mpos(i_unstim_trial)));
				end
				% 
				% 	Store the variance ratio for statistics (for now will ignore fact that the number of datapoints is not equal between sets)
				% 
				obj.valid_blocks.m.Std.cue2lick_win(i_div).Robin_variance_ratio = nanvar(cell2mat(obj.valid_blocks.m.data{1, i_div}{1,2}.'))/nanvar(cell2mat(obj.valid_blocks.m.data{1, i_div}{1,1}.')); % var unstim / var stim '
				% 
				% 
				% 	Gather the variance data for each trial in the movement signal --- moving to end now so using cue2lick window (6/28/18)
				% 		Variance within CUE TO MOVEMENT WINDOW (25ms smooth)
				%
				%		Indicies of trials
			    % 
				obj.valid_blocks.m.Std.cue2lick_win(i_div).stim_idx = obj.valid_blocks.m.trial_nums(i_div).stim;
				obj.valid_blocks.m.Std.cue2lick_win(i_div).unstim_idx = obj.valid_blocks.m.trial_nums(i_div).unstim;
				%
				%		Standard deviation by trial
				%
				obj.valid_blocks.m.Std.cue2lick_win(i_div).unstim_Std_bt = cellfun(@nanvar,obj.valid_blocks.m.data{1,1}{1,2});
				% 
				obj.valid_blocks.m.Std.cue2lick_win(i_div).stim_Std_bt = cellfun(@nanvar,obj.valid_blocks.m.data{1,1}{1,1});
				% 
				% 		Average std by trial for each case ------ MOVED TO LATER SO WE CAN USE CONCAT DATASET INSTEAD OF AVE!
				% 
				obj.valid_blocks.m.Std.cue2lick_win(i_div).unstim_Std_concat = nanstd(cell2mat(obj.valid_blocks.m.data{1, i_div}{1,2}.')); %'
				obj.valid_blocks.m.Std.cue2lick_win(i_div).stim_Std_concat = nanstd(cell2mat(obj.valid_blocks.m.data{1, i_div}{1,1}.')); %'
			end
		end
		% 
		% 	NTJ divvying by categories
		% 
		function NTJcategories(obj, preStimWindow)
			if nargin < 2
				preStimWindow = 0.5;
			end
			% 
			% 	We should divide the NTJ data into categories and then have them stored for plotting
			obj.valid_blocks = [];
			obj.valid_blocks.preStimWindow = preStimWindow;
			% 
			% 	1. Identify all trials with juice dispensed in pretrial
			% 
			% obj.valid_blocks.preStimJuiceTrials = find(ismember(round(obj.time.s.stimOn_s,3), round(obj.time.s.juice_s,3)+0.5));
			% obj.valid_blocks.NOpreStimJuiceTrials = find(~ismember(round(obj.time.s.stimOn_s,3), round(obj.time.s.juice_s,3)+0.5));
			% obj.valid_blocks.Juice_Idxs_Pre_Stim = find(ismember(round(obj.time.s.juice_s,3)+0.5, round(obj.time.s.stimOn_s,3)));
			obj.valid_blocks.preStimJuiceTrials = find(ismember(round(obj.time.s.stimOn_s,1), round(obj.time.s.juice_s,1)+0.5));
			obj.valid_blocks.NOpreStimJuiceTrials = find(~ismember(round(obj.time.s.stimOn_s,1), round(obj.time.s.juice_s,1)+0.5));
			obj.valid_blocks.Juice_Idxs_Pre_Stim = find(ismember(round(obj.time.s.juice_s,1)+0.5, round(obj.time.s.stimOn_s,1)));
			% 
			% 	2. Identify which juice events are not in the trial period
			% 		Start by finding to which trial each juice event belongs
			% 
			[~, ~, juice_events] = histcounts(obj.time.s.juice_s,obj.time.s.stimOn_s);
			[obj.valid_blocks.First_Juice_Trial_Number, obj.valid_blocks.juiceIdxIsFirstInTrial, ~] = unique(juice_events);
			if obj.valid_blocks.juiceIdxIsFirstInTrial(1) == 0
				obj.valid_blocks.juiceIdxIsFirstInTrial = obj.valid_blocks.juiceIdxIsFirstInTrial(2:end);
				obj.valid_blocks.First_Juice_Trial_Number = obj.valid_blocks.First_Juice_Trial_Number(2:end);
			end
			obj.valid_blocks.Juice_Trial_Number = nan(numel(obj.time.s.juice_s),1);
			obj.valid_blocks.Juice_Trial_Number(obj.valid_blocks.juiceIdxIsFirstInTrial) = obj.valid_blocks.First_Juice_Trial_Number;
			obj.valid_blocks.Juice_Trial_Number = fillmissing(obj.valid_blocks.Juice_Trial_Number,'previous');
			obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim = [];
			for iTrial = 1:numel(obj.time.s.stimOn_s)-1
				obj.time.s.juice_idxs_by_trial{iTrial} = find(obj.valid_blocks.Juice_Trial_Number == iTrial);
				obj.time.s.juice_times_by_trial{iTrial} = obj.time.s.juice_s(obj.time.s.juice_idxs_by_trial{iTrial});
				obj.time.s.juice_times_by_trial_wrtSO{iTrial} = obj.time.s.juice_times_by_trial{iTrial} - obj.time.s.stimOn_s(iTrial);
				% 
				% 		Now find juice idxs not in the stimulation period. Also ignores those within preStimWindow of the stimON
				% 
				idxsToAppend = obj.time.s.juice_idxs_by_trial{iTrial}(find(obj.time.s.juice_times_by_trial_wrtSO{iTrial} > obj.time.Plot.stim_off_time_s & ~(obj.time.s.juice_times_by_trial{iTrial} > obj.time.s.stimOn_s(iTrial+1))));
				obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim = [obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim;idxsToAppend];
            end
            obj.time.s.juice_idxs_by_trial{iTrial+1} = find(obj.valid_blocks.Juice_Trial_Number == iTrial+1);
            obj.time.s.juice_times_by_trial{iTrial+1} = obj.time.s.juice_s(obj.time.s.juice_idxs_by_trial{iTrial+1});
            obj.time.s.juice_times_by_trial_wrtSO{iTrial+1} = obj.time.s.juice_times_by_trial{iTrial+1} - obj.time.s.stimOn_s(iTrial+1);
            % 
            % 		Now find juice idxs not in the stimulation period. Also ignores those within preStimWindow of the stimON
            % 
            idxsToAppend = obj.time.s.juice_idxs_by_trial{iTrial+1}(find(obj.time.s.juice_times_by_trial_wrtSO{iTrial+1} > obj.time.Plot.stim_off_time_s));
            obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim = [obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim;idxsToAppend];
			% 
			% 	Remove those juice idxs that are in the pre-stim -0.5
			% 
			obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim(find(ismember(obj.valid_blocks.Juice_Idxs_Not_In_Stim_or_Pre_Stim, obj.valid_blocks.Juice_Idxs_Pre_Stim))) = [];
		end
		% 
		% 	Statistical tests runner
		% 
		function analyze(obj, collateKey, shellObj)
			% 
			% 	collateKey: 	A keyword indicating the analysis to complete
			% 				bootAUC
			%				emeryBoot 
			% 				changeTypeA 	changes obj type to stimA--to fix an error with the processing
			% 				ntRaster 		gathers raster data for all the naive test - juice datasets
			% 				ntjRaster 		gathers raster data for all the naive test + juice datasets
			%
			%	Default Protocol A delAUC_EOT bootstrap: 
			% 		obj.boot('lick', obj.bootStat.considerationBoundary_ms, 'delCDF_AUC_EOT', 10000, 'valid', 0)
			% 	Default Emery bootstrap
			% 		obj.boot('lick', obj.bootStat.considerationBoundary_ms, 'Emery', 1000000, 'EOT', 0)
			% ----------
			% 
			% 	Set default consideration boundary
			% 
			if nargin < 2
				collateKey = 'Emery';
				disp('-------------default analysis: Emery from consideration boundary to EOT=7s, post EOT set to 7s---------------')
			end
			obj.setConsiderationBound(700);
			if strcmpi(collateKey, 'bootAUC')
				% 
				% 	This only makes sense really to do for stim-typeA objs, so update that here
				% 
				obj.changeObjStimType('A');
				obj.bootStat.note = 'Changed stim type to A';
				obj.getvalidblocks;
				obj.boot('lick', obj.bootStat.considerationBoundary_ms, 'delCDF_AUC_EOT', 10000, 'valid', 0);
				legend('hide');
				tstamp = datestr(now, 'yyyy_mm_dd_HHMM');
				figurename = ['bootAUC_optolegacy_' obj.iv.init_variables_ChR2.mousename_ '_d' obj.iv.init_variables_ChR2.daynum_ '_' tstamp '_' num2str(shellObj.iv.runID)];
				set(gcf, 'name', figurename);
				switchdir = pwd;
				cd(shellObj.iv.figuresDir);
				printFigure(figurename, gcf);
				cd(switchdir);
			elseif strcmpi(collateKey, 'Emery')
				% 
				% 	This only makes sense really to do for stim-typeA objs, so update that here
				% 
				obj.changeObjStimType('A');
				obj.bootStat.note = 'Changed stim type to A';
				obj.getvalidblocks;
				obj.boot('lick', obj.bootStat.considerationBoundary_ms, 'Emery', 1000000, 'EOT', 0);
				legend('hide');
				tstamp = datestr(now, 'yyyy_mm_dd_HHMM');
				figurename = ['Emery_optolegacy_' obj.iv.init_variables_ChR2.mousename_ '_d' obj.iv.init_variables_ChR2.daynum_ '_' tstamp '_' num2str(shellObj.iv.runID)];
				set(gcf, 'name', figurename);
				switchdir = pwd;
				cd(shellObj.iv.figuresDir);
				printFigure(figurename, gcf);
				cd(switchdir);
			elseif strcmpi(collateKey, 'changeTypeA')
				obj.changeObjStimType('A');
				obj.bootStat.note = 'Changed stim type to A';
			elseif strcmpi(collateKey, 'ntRaster')
				shellObj.bootStat.note = 'concatenated raster data across sessions for naive test - juice';
				obj.getvalidblocks;
				data = obj.ntRasterCollateHelper;
				if ~isfield(shellObj.iv, 'init_variables_ChR2')
					shellObj.time.ms.total_time_ms = obj.time.ms.total_time_ms;
					shellObj.time.s.binnedLicks.ref = data.ref;
				    shellObj.time.s.rxnwin_s = data.rxnwin_s;
				    shellObj.iv.init_variables_ChR2.time_parameters.samples_per_ms = data.fs_ms;
				    shellObj.time.ms.op_rew_open_ms = data.rb_ms;
				    shellObj.time.s.ITI_time_s = data.eot;
				    shellObj.time.ms.cue_on_time_ms = data.cue;
				    shellObj.bootStat.considerationBoundary_ms = data.considbound;					
				end
				if ~isfield(shellObj.bootStat,'collateResult')
					idx = 1;
				else
					idx = numel({shellObj.bootStat.collateResult.ID})+1;
				end
				shellObj.bootStat.collateResult(idx).ID = shellObj.iv.files{idx};
				shellObj.bootStat.collateResult(idx).all_trials_this_block = data.all_trials_this_block;
				shellObj.bootStat.collateResult(idx).all_lick_times_ex_swrtc = data.all_lick_times_ex_swrtc;
				shellObj.bootStat.collateResult(idx).f_lick_s_wrtcue = data.f_lick_s_wrtcue;
				shellObj.bootStat.collateResult(idx).stimTrials = data.stimTrials;
				shellObj.bootStat.collateResult(idx).nostim_trials = data.nostim_trials;	
				shellObj.bootStat.collateResult(idx).fs_ms = data.fs_ms;
				if shellObj.iv.init_variables_ChR2.time_parameters.samples_per_ms ~= data.fs_ms
					shellObj.bootStat.collateResult(idx).FLAG = 'FLAG - this dataset had unexpected sampling rate for ChR2';
					warning('FLAG - this dataset had unexpected sampling rate for ChR2')
				end			
			elseif strcmpi(collateKey, 'ntjRaster')
				shellObj.bootStat.note = 'concatenated raster data across sessions for naive test + juice';
				data = obj.ntjRasterCollateHelper;
				if ~isfield(shellObj.iv, 'fs_ms')
                    shellObj.iv.fs_ms = data.fs_ms;
                    shellObj.iv.sb_s = data.sb_s;
                    shellObj.iv.se_s = data.se_s;
                    shellObj.iv.rfx_s = data.rfx_s;
                    shellObj.iv.max_s = data.max_s;
                    shellObj.iv.considbound = data.considbound;
				end
				if ~isfield(shellObj.bootStat,'collateResult')
					idx = 1;
				else
					idx = numel({shellObj.bootStat.collateResult.ID})+1;
				end
				shellObj.bootStat.collateResult(idx).ID = shellObj.iv.files{idx};
				shellObj.bootStat.collateResult(idx).time = data.time;
				shellObj.bootStat.collateResult(idx).all_licks_swrtSO_ex_s = data.all_licks_swrtSO_ex_s;
				shellObj.bootStat.collateResult(idx).all_flicks_swrtSO_ex_s = data.all_flicks_swrtSO_ex_s;
				shellObj.bootStat.collateResult(idx).all_licks_swrtjuice_s = data.all_licks_swrtjuice_s;
				shellObj.bootStat.collateResult(idx).juice_s_bt_wrtref = data.juice_s_bt_wrtref;
				shellObj.bootStat.collateResult(idx).valid_blocks = data.valid_blocks;
				if shellObj.iv.fs_ms ~= data.fs_ms
					shellObj.bootStat.collateResult(idx).FLAG = 'FLAG - this dataset had unexpected sampling rate for ChR2';
					warning('FLAG - this dataset had unexpected sampling rate for ChR2')
				end
			elseif strcmpi(collateKey, 'ecdf')
				obj.plot('ecdf')
                if ~isfield(shellObj.bootStat,'collateResult')
					idx = 1;
				else
					idx = numel({shellObj.bootStat.collateResult.ID})+1;
				end
				[h,p,ks2stat] = kstest2(obj.valid_blocks.eCDF.unstim_x_EOTinclusive, obj.valid_blocks.eCDF.stim_x_EOTinclusive);
				shellObj.bootStat.collateResult(idx).ID = shellObj.iv.files{idx};
				shellObj.bootStat.collateResult(idx).note = 'ecdf EOT in KS test (0.7 to 7s, consecutive)';
				shellObj.bootStat.collateResult(idx).h = h;
				shellObj.bootStat.collateResult(idx).p = p;
				shellObj.bootStat.collateResult(idx).ks2stat = ks2stat;
				disp(['		p: ' num2str(p)])	
                title(['KStest2 p: ' num2str(p)]);
				xlim([0.7,7])
				xticks([0.7,1,2,3,4,5,6,7])
				tstamp = datestr(now, 'yyyy_mm_dd_HHMM');
				figurename = ['ecdf_optolegacy_' obj.iv.init_variables_ChR2.mousename_ '_d' obj.iv.init_variables_ChR2.daynum_ '_' tstamp '_' num2str(shellObj.iv.runID)];
				set(gcf, 'name', figurename);
				switchdir = pwd;
				cd(shellObj.iv.figuresDir);
				printFigure(figurename, gcf);
				cd(switchdir);

							
			else
				error('undefined canned analysis')
			end
		end
		function resetBootResults(obj)
			% 
			% 	Reset bootStat field for single-session obj
			% 
			obj.bootStat.resultLick = {};
			obj.bootStat.resultX = {};
			obj.bootStat.resultLick.nboot = {};
			obj.bootStat.resultX.nboot = {};
		end
		function boot(obj, Signal, considerationBoundary_ms, Mode, nboot, Blocking, saveData)
			% 
			% 	>> boot
			%  OR call a canned analysis with
			% 	>> analyze >> boot
			% -------------------------------------------------------------------------------------
			% 	Blocking: 	'valid' for valid blocks
			% 				'even' 	for even blocks
			% 		"EDGE CASE HANDLING" ---- HAS DIFFERENT MEANING FOR Mode = Emery:
			% 				'EOT': we will set all lick times > EOT = EOT in s
			% 				'omit': we remove any trials with lick times > EOT from the pool
			% 
			% 	Mode:		The Bootstrap style
			% 				'median'			for median
			% 				'delCDF_AUC'		compare entire ecdf up to end of ITI
			% 				'delCDF_AUC_EOT'	compare ecdf up to end of trial
			% 				'Emery' 			Draws a stim and unstim trial from a pool of trials and compares the difference in lick times
			% 									NB! This mode is handled COMPLETELY differently from legacy modes
			% 
			% 	Signal 		Whether to use behavior or X variance
			% 
			% 				'lick'
			% 				'X'
			%	Default Protocol A delAUC_EOT bootstrap: 
			% 		obj.boot('lick', obj.bootStat.considerationBoundary_ms, 'delCDF_AUC_EOT', 10000, 'valid', 0)
			% 	Default Emery bootstrap
			% 		% obj.boot('lick', obj.bootStat.considerationBoundary_ms, 'Emery', 1000000, 'EOT', 0)
			% 
			if nargin < 2
				Signal = 'lick';
			end
			if nargin < 3
				considerationBoundary_ms = obj.bootStat.considerationBoundary_ms;
			end
			if nargin < 4 || isempty(Mode)
				if strcmpi(Signal, 'X')
					Mode = 'ts_variance_ratio';
				elseif strcmpi(Signal, 'lick')
					Mode = 'Emery';
					% Mode = 'delCDF_AUC_EOT';
				end
			end
			if nargin < 5 && strcmpi(Signal, 'lick')
				if ~strcmpi(Mode, 'emery')
					nboot = 10000;
				else
					nboot = 1000000;
				end
			elseif nargin < 5 && strcmpi(Signal, 'X')
				nboot = 1000;
			end
			if nargin < 6 || isempty(Blocking)
				if strcmpi(Mode, 'emery')
					% 
					% 	Options:
					% 		EOT, omit, EndITI, 'none'
					% 
					Blocking = 'EndITI';
				else
					Blocking = 'valid';
				end
			end
			if nargin < 7
				saveData = false;
			end

			
			% ptrs
			EOT = obj.time.s.ITI_time_s;
			Even = false; %run_.boot_matched_numbers_stim_unstim = false from alpha testing
			
			if strcmpi(Mode, 'emery')
				EdgeCaseHandling = Blocking;
				Signal = 'lick';
		    	bootTestNo = length([obj.bootStat.resultLick.nboot]) + 1;
		    	Even = 'n/a';

		    	obj.getBinnedLicks('cue', 30000, 30000, 0, 0);
				
				ref = 'n/a';
			    %
			    obj.bootStat.resultLick(bootTestNo).considerationBoundary_ms = considerationBoundary_ms;     
			    obj.bootStat.resultLick(bootTestNo).nboot = nboot;
			    obj.bootStat.resultLick(bootTestNo).Mode = Mode;
			    obj.bootStat.resultLick(bootTestNo).EdgeCaseHandling = EdgeCaseHandling;
			    obj.bootStat.resultLick(bootTestNo).EOT = EOT;

			    obj.emeryBoot(nboot, EOT, EdgeCaseHandling);
			    
			elseif strcmpi(Blocking, 'even')
				error('Not implemented')
			elseif strcmpi(Blocking, 'valid')
				% 
				% 	Recalculate the valid blocks in case we changed the consideration boundary
				% 
				obj.getvalidblocks
				% 
				obj.bootStat.fl.stim = obj.valid_blocks.f_lick_s_wrtcue.stim;
				obj.bootStat.fl.unstim = obj.valid_blocks.f_lick_s_wrtcue.unstim;
				% 
				% 	Handle the consideration boundary for flicks
				% 	
				obj.bootStat.fl.stim(obj.bootStat.fl.stim < considerationBoundary_ms/1000) = nan;
				obj.bootStat.fl.unstim(obj.bootStat.fl.unstim < considerationBoundary_ms/1000) = nan;
				% 
				% 
				obj.bootStat.fl.stim = {obj.bootStat.fl.stim}; 
				obj.bootStat.fl.unstim = {obj.bootStat.fl.unstim};


                for iblock = 1:sum(arrayfun(@(x) numel(x), obj.valid_blocks.Median))
                    obj.bootStat.median_unstim_minus_stim(iblock) = obj.valid_blocks.Median(iblock).unstim - obj.valid_blocks.Median(iblock).stim;
                    obj.bootStat.fl.data{iblock} = {obj.bootStat.fl.stim{iblock}, obj.bootStat.fl.unstim{iblock}};
                end
				data = obj.bootStat.fl.data;
				if strcmpi(Signal, 'lick')
					bootTestNo = length([obj.bootStat.resultLick.nboot]) + 1;
					if strcmpi(Mode, 'median')
						ref = obj.valid_blocks.median_unstim_minus_stim;
					elseif strcmpi(Mode, 'delCDF_AUC')
						ref = [obj.valid_blocks.eCDF.delAUCun_minus_stim];
					elseif strcmpi(Mode, 'delCDF_AUC_EOT')
						ref = [obj.valid_blocks.eCDF.delAUC_EOTun_minus_stim];
					end
					valid_blocks_bootobj = CLASS_bootstrap_roadmapv1_4(data, 'stim', 'no stim', ref, EOT);
					% --------------------------------------------------
					% 				PERMUTATION
					% --------------------------------------------------
					% 
					valid_blocks_bootobj.permute(nboot, Even, Mode, obj.iv.runID);
				    %
				    obj.bootStat.resultLick(bootTestNo).considerationBoundary_ms = considerationBoundary_ms;     
				    obj.bootStat.resultLick(bootTestNo).nboot = nboot;
				    obj.bootStat.resultLick(bootTestNo).Mode = Mode;
				    obj.bootStat.resultLick(bootTestNo).p = valid_blocks_bootobj.PermutedBlocks.Block.p;
				    obj.bootStat.resultLick(bootTestNo).ref = ref;
				    obj.bootStat.resultLick(bootTestNo).EOT = EOT;
				    obj.bootStat.resultLick(bootTestNo).Even = Even;
				elseif strcmpi(Signal, 'X')
					saveData = true;
					bootTestNo = length([obj.bootStat.resultX.nboot]) + 1;
					ref = [obj.valid_blocks.m.Std.cue2lick_win.Robin_variance_ratio];
			        valid_blocks_bootobj = CLASS_bootstrap_roadmapv1_4(obj.valid_blocks.m.data, 'stim', 'no stim', ref, EOT);
					% 
					valid_blocks_bootobj.permute(nboot, false, Mode, obj.iv.runID);
					% 
					obj.bootStat.resultX(bootTestNo).considerationBoundary_ms = considerationBoundary_ms;     
				    obj.bootStat.resultX(bootTestNo).nboot = nboot;
				    obj.bootStat.resultX(bootTestNo).Mode = Mode;
				    obj.bootStat.resultX(bootTestNo).p = valid_blocks_bootobj.PermutedBlocks.Block.p;
				    obj.bootStat.resultLick(bootTestNo).unstimvar_div_stimvar_stat = valid_blocks_bootobj.PermutedBlocks.Block.Statistic.C2_div_C1;
				    obj.bootStat.resultX(bootTestNo).ref = ref;
				    obj.bootStat.resultX(bootTestNo).EOT = EOT;
				    obj.bootStat.resultX(bootTestNo).Even = Even;
				end
				% 
				%% --- Save data ----
				if saveData
					savefilename = ['BootstrapResults_' Signal, '_' obj.iv.init_variables_ChR2.mousename_, '_day', obj.iv.init_variables_ChR2.daynum_, '_' datestr(now, 'YYYYmmDD_HH_MM') '_runIDno' num2str(obj.iv.runID)];
					save([savefilename, '.mat'], 'valid_blocks_bootobj', '-v7.3');	
					obj.save;	
				end
			elseif strcmpi(Blocking, 'composite')
				error('Not implemented')
			end
		end
		function emeryBoot(obj, nboot, EOT, Mode) % seems to be working correctly -- stepped thru it all on 11/4/19 13:21
			% 
			% 	>> analyze >> boot > emeryBoot
			% 
			Notes{1,1} = ['Emery Boot Run ' datestr(now)];
			bootTestNo = length([obj.bootStat.resultLick.nboot]);
			% 
			% 	EOT: time in s after which licks should be omitted as edge case (default=7s)
			% 	
			% 	Mode: how to handle the edge case
			% 			'EOT': 	sets all lick times > EOT = EOT in sec and keeps these in the dataset 	-- this is like plotting the cdf AUC_EOT -- keeps those later trials in the mix
			% 			'omit':	omits all trials with lick time > EOT or nan 		-- this is like plotting the cdf up till 7s, expect effect to wash out
			% 			'EndITI': sets all trials with lick time > end of ITI = time of EndITI in s wrt cue
			% 			'none': keeps all > EOT licks
			% 			
			% 	N.B. that trials with lick < considerationBoundary_ms will be omitted from pool, as will those with NEVER-LICK for EOT:
			% 
			% 	NB! If EOT, really we should not be including trials when mouse never licked at all - it's just totally undefined behavior and will tend to pull the median toward zero. 
			% 	Essentially we are including a lot of trials we should have excluded before analysis. So I will leave these out entirely to better emulate the original AUC method
			% 	So we only want to include trials where the mouse actually did something. If does nothing, we have no idea if was trying at all
			% 	ALSO NB that the AUC method ignored any non-consecutive trials, which also kinda cleaned things up. Here, we didn't do that, we just ignored any trials 
			% 	with NEVER lick, so they are slightly different, but I don't believe the analysis gives any different result.
			% 	And in this way is a bit more conservative I suppose.
			% 
			% 	0. Handle excluded trials
			%
			flicks = obj.time.s.binnedLicks.f_lick_s_wrtref; 
			flicks(obj.iv.ExcludedTrials) = nan;
			% 
			stim_pool = flicks(obj.h.stimTrials);
			unstim_pool = flicks(obj.h.nostim_trials);
			nStimTrials = numel(obj.h.stimTrials);
			nUnstimTrials = numel(obj.h.nostim_trials);
			% 
			% 	1. Handle the consideration boundary
			% 
			considbound = obj.bootStat.considerationBoundary_ms/1000;
			nStimRxns = sum(stim_pool < considbound);
			nUnstimRxns = sum(unstim_pool < considbound);
			pcStimRxns = nStimRxns/numel(stim_pool);
			pcUnstimRxns = nUnstimRxns/numel(unstim_pool);
			nStimNoLick = sum(stim_pool > EOT);
			nUnstimNoLick = sum(unstim_pool > EOT);
			pcStimNoLick = nStimNoLick/numel(stim_pool);
			pcUnstimNoLick = nUnstimNoLick/numel(unstim_pool);
			% 
			% 	2. Remove excluded trials
			% 
			nStimExcluded = sum(isnan(stim_pool));
			nUnstimExcluded = sum(isnan(unstim_pool));
			stim_pool(stim_pool < considbound) = [];
			unstim_pool(unstim_pool < considbound) = [];
			Notes{2,1} = ['Deleted Trials with Lick Time < Consideration Boundary = ' num2str(considbound) 's.'];
			stim_pool(isnan(stim_pool)) = [];
			unstim_pool(isnan(unstim_pool)) = [];
			Notes{3,1} = ['Deleted Excluded Trials:' mat2str(obj.iv.ExcludedTrials) '.'];
			% 
			% 	3. Handle Edge Cases
			% 
			EndITI = obj.time.s.total_time_s;
			if strcmpi(Mode, 'EOT')
				stim_pool(stim_pool > EndITI) = [];
				unstim_pool(unstim_pool > EndITI) = [];
				stim_pool(stim_pool > EOT) = EOT;
				unstim_pool(unstim_pool > EOT) = EOT;
				Notes{2,1} = ['Set Trials with Lick Time > EOT to EOT = ' num2str(EOT) 's and remove never-lick trials.'];
			elseif strcmpi(Mode, 'omit')
				stim_pool(stim_pool > EOT) = [];
				unstim_pool(unstim_pool > EOT) = [];
				Notes{2,1} = ['Deleted Trials with Lick Time > EOT, EOT = ' num2str(EOT) 's.'];
			elseif strcmpi(Mode, 'EndITI')
				stim_pool(stim_pool > EndITI) = EndITI;
				unstim_pool(unstim_pool > EndITI) = EndITI;
				Notes{2,1} = ['Set Trials with Lick Time > EndITI to EndITI = ' num2str(EndITI) 's.'];
			else
				Notes{2,1} = ['Included Trials with Lick Time > EOT and NEVER licks even if flick in the subsequent trial (never lick case), EOT = ' num2str(EOT) 's.'];
			end
			% 
			% 	4. Run the bootstrap
			% 
			nUnstimPool = numel(unstim_pool);
			nStimPool = numel(stim_pool);

			Notes{4,1} = ['Bootstrapped STIM - UNSTIM lick times in s. Results reported as lick time (stim) - lick time(unstim) in s'];
			bootLickTimeDifferences_StimMinusUnstim = nan(nboot,1);
			for iboot = 1:nboot
				stimDrawIdx = randi(nStimPool);
				unstimDrawIdx = randi(nUnstimPool);
				bootLickTimeDifferences_StimMinusUnstim(iboot) = stim_pool(stimDrawIdx) - unstim_pool(unstimDrawIdx);
			end
			result = sort(bootLickTimeDifferences_StimMinusUnstim);
			ll = round(0.025*nboot);
			ul = round(0.975*nboot);
			CIl = result(ll);
			CIu = result(ul);

			H = figure;
			axhxg = subplot(1,2,1);
			axbox = subplot(1,2,2);
			hold(axhxg, 'on');
			hold(axbox, 'on');
			histogram(axhxg, result, 'BinMethod', 'integer',  'Normalization', 'probability', 'displayStyle', 'stairs', 'DisplayName', 'Stim Lick Time - Unstim Lick Time');
			yy = get(axhxg, 'ylim');
			plot(axhxg, [CIl, CIl], [0,1], 'r-', 'DisplayName', '95% CI bound')
			plot(axhxg, [CIu, CIu], [0,1], 'r-', 'DisplayName', '95% CI bound')
			plot(axhxg, [mean(result), mean(result)], [0,1], 'r-', 'DisplayName', 'Mean')
			plot(axhxg, [median(result), median(result)], [0,1], 'r-', 'DisplayName', 'Mean')
			ylim(axhxg, yy);
			title(axhxg, ['Distribution of Bootstrapped Differences in Lick Time, nboot: ' num2str(nboot)]);
			axhxg.XLabel.String = 'Lick Time Difference (s)';
			axhxg.YLabel.String = 'Probability';

			boxplot(axbox, result);
			title(axbox, ['5 Number Summary of Bootstrapped Differences in Lick Time, nboot: ' num2str(nboot)]);
			axhxg.YLabel.String = 'Lick Time Difference (s)';
			FILENAME = ['EmeryBootResult_Mode' Mode '_n' num2str(nboot) '_ID' num2str(obj.iv.runID), '_' datestr(now, 'YYYYmmDD_HHMM')];
			savefig(H, FILENAME)

			obj.bootStat.resultLick(bootTestNo).median = median(bootLickTimeDifferences_StimMinusUnstim);
			obj.bootStat.resultLick(bootTestNo).mean = mean(bootLickTimeDifferences_StimMinusUnstim);
			obj.bootStat.resultLick(bootTestNo).CI95 = [CIl, CIu];
			obj.bootStat.resultLick(bootTestNo).Notes = Notes;
			obj.bootStat.resultLick(bootTestNo).nTrials.nStimTrials = nStimTrials;
			obj.bootStat.resultLick(bootTestNo).nTrials.nUnstimTrials = nUnstimTrials;
			obj.bootStat.resultLick(bootTestNo).nTrials.nUnstimPool = nUnstimPool;
			obj.bootStat.resultLick(bootTestNo).nTrials.nStimPool = nStimPool;
			obj.bootStat.resultLick(bootTestNo).nTrials.nUnstimPool = numel(unstim_pool);
			obj.bootStat.resultLick(bootTestNo).nTrials.nStimRxns = nStimRxns;
			obj.bootStat.resultLick(bootTestNo).nTrials.nUnstimRxns = nUnstimRxns;
			obj.bootStat.resultLick(bootTestNo).nTrials.pcStimRxns = pcStimRxns;
			obj.bootStat.resultLick(bootTestNo).nTrials.pcUnstimRxns = pcUnstimRxns;
			obj.bootStat.resultLick(bootTestNo).nTrials.nStimNoLick = nStimNoLick;
			obj.bootStat.resultLick(bootTestNo).nTrials.nUnstimNoLick = nUnstimNoLick;
			obj.bootStat.resultLick(bootTestNo).nTrials.pcStimNoLick = pcStimNoLick;
			obj.bootStat.resultLick(bootTestNo).nTrials.pcUnstimNoLick = pcUnstimNoLick;
			obj.bootStat.resultLick(bootTestNo).nTrials.nStimExcluded = nStimExcluded;
			obj.bootStat.resultLick(bootTestNo).nTrials.nUnstimExcluded = nUnstimExcluded;
			obj.bootStat.resultLick(bootTestNo).Pools.stim_pool = stim_pool;
			obj.bootStat.resultLick(bootTestNo).Pools.unstim_pool = unstim_pool;
			obj.bootStat.resultLick(bootTestNo).Pools.bootLickTimeDifferences_StimMinusUnstim = result;
		end
		% 
		% 	Compares lick time differences without shuffling -- e.g., paired trials, etc
		% 
		function results = lickTimeDiffDistribution(obj, Mode, Pooling)
			% 	
			% 	Mode: how to handle the edge case
			% 			'EOT': 	sets all lick times > EOT = EOT in sec and keeps these in the dataset 	-- this is like plotting the cdf AUC_EOT -- keeps those later trials in the mix
			% 			'omit':	omits all trials with lick time > EOT or nan 		-- this is like plotting the cdf up till 7s, expect effect to wash out
			% 			'EndITI': sets all trials with lick time > end of ITI = time of EndITI in s wrt cue
			% 			'none': keeps all > EOT licks
			% 
			% 	Pooling: how we determine neighbors to compare to, for example, we'd like to be in the same neighborhood to deal with the non-stationarity, but we also might want to do more than paired
			% 			'Paired' 	-- takes adjacent trials
			% 			'Sequence' 	-- for every stim trial, compares to all prior unstim trials
			% 			'Web' 		-- each stim strial compared to nearest 10 neighboring unstim trials on either side, allowing overlap
			% 			
			% 	N.B. that trials with lick < considerationBoundary_ms will be omitted from pool.
			% 
			if nargin < 2
				Mode = 'none';
			end
			if nargin < 3
				Mode = 'paired';
			end
			obj.getBinnedLicks('cue', 30000, 30000, 0, 0);
            EOT = obj.time.s.ITI_time_s;
			% 
			% 	Handle excluded trials, reactions, and things outside of consideration bound
			% 
			flicks = obj.time.s.binnedLicks.f_lick_s_wrtref; 
			flicks(obj.iv.ExcludedTrials) = nan;
			flicks(flicks < obj.bootStat.considerationBoundary_ms/1000) = nan;
			if strcmpi(Mode, 'EOT')
				flicks(flicks > EOT) = EOT;
				Notes{1} = ['Set Trials with Lick Time > EOT to EOT = ' num2str(EOT) 's.'];
			elseif strcmpi(Mode, 'omit')
				flicks(flicks > EOT) = nan;
				Notes{1} = ['Deleted Trials with Lick Time > EOT, EOT = ' num2str(EOT) 's.'];
			elseif strcmpi(Mode, 'EndITI')
				EndITI = obj.time.s.total_time_s;
				flicks(flicks > EndITI) = EndITI;
				Notes{1} = ['Set Trials with Lick Time > EndITI to EndITI = ' num2str(EndITI) 's.'];
			else
				Notes{1} = ['Included Trials with Lick Time > EOT even if flick in the subsequent trial (never lick case), EOT = ' num2str(EOT) 's.'];
			end
			if strcmpi(Pooling, 'paired')
				% 
				% 	Get the paired stim/unstim trials across the session
				% 
				stim_then_unstim = obj.h.stimTrials(ismember(obj.h.stimTrials, obj.h.nostim_trials-1));
				unstim_then_stim = obj.h.nostim_trials(ismember(obj.h.nostim_trials, obj.h.stimTrials-1));
				assert(sum(ismember(stim_then_unstim, unstim_then_stim))==0);
				% 
				% 	Now get the difference in lick time between the pairs
				% 
				lickTimeDifference_stim_minus_unstim = nan(numel(stim_then_unstim)+numel(unstim_then_stim), 1);
				lickTimeDifference_stim_minus_unstim(1:numel(stim_then_unstim)) = flicks(stim_then_unstim) - flicks(stim_then_unstim+1);
				lickTimeDifference_stim_minus_unstim(numel(stim_then_unstim)+1:end) = flicks(unstim_then_stim+1) - flicks(unstim_then_stim);
			elseif strcmpi(Pooling, 'sequence')
				lickTimeDifference_stim_minus_unstim = [];
				for iStimtrial = 1:numel(obj.h.stimTrials)
					if iStimtrial == 1 && obj.h.stimTrials(iStimtrial) ~= 1
						ntrialsthischunk = obj.h.stimTrials(iStimtrial) - 1;
						chunk = flicks(obj.h.stimTrials(iStimtrial)) - flicks(1:obj.h.stimTrials(iStimtrial)-1); 
					elseif obj.h.stimTrials(iStimtrial) == 1
						ntrialsthischunk = 0;
						chunk = [];
					else
						ntrialsthischunk = obj.h.stimTrials(iStimtrial) - obj.h.stimTrials(iStimtrial-1) - 1;
						chunk = flicks(obj.h.stimTrials(iStimtrial)) - flicks(obj.h.stimTrials(iStimtrial-1)+1:obj.h.stimTrials(iStimtrial)-1); 
					end
					lickTimeDifference_stim_minus_unstim = [lickTimeDifference_stim_minus_unstim; chunk];
				end

			elseif strcmpi(Pooling, 'web')
				lickTimeDifference_stim_minus_unstim = [];
				for iStimtrial = 1:numel(obj.h.stimTrials)
					preUnstim = obj.h.nostim_trials(find(obj.h.nostim_trials < obj.h.stimTrials(iStimtrial), 10, 'last'));
					postUnstim = obj.h.nostim_trials(find(obj.h.nostim_trials > obj.h.stimTrials(iStimtrial), 10, 'first'));
					unstimidx =  [preUnstim;postUnstim];
					chunk = flicks(obj.h.stimTrials(iStimtrial)) - flicks(unstimidx); 
					lickTimeDifference_stim_minus_unstim = [lickTimeDifference_stim_minus_unstim; chunk];
				end
			end
			% 
			% 	Ignore any cases where we excluded the trial for some reason
			% 			
			lickTimeDifference_stim_minus_unstim(isnan(lickTimeDifference_stim_minus_unstim)) = [];
			% 
			% 
			% 
			result = sort(lickTimeDifference_stim_minus_unstim);
			ll = round(0.025*numel(lickTimeDifference_stim_minus_unstim));
			ul = round(0.975*numel(lickTimeDifference_stim_minus_unstim));
			CIl = result(ll);
			CIu = result(ul);
			
			figure
			axhxg = subplot(1,2,1);
			axbox = subplot(1,2,2);
			hold(axhxg, 'on');
			hold(axbox, 'on');
			histogram(axhxg, result, 'BinMethod', 'integer', 'Normalization', 'probability', 'displayStyle', 'stairs', 'DisplayName', 'Stim Lick Time - Unstim Lick Time');
			yy = get(axhxg, 'ylim');
			plot(axhxg, [CIl, CIl], [0,1], 'r-', 'DisplayName', '95% CI bound')
			plot(axhxg, [CIu, CIu], [0,1], 'r-', 'DisplayName', '95% CI bound')
			ylim(axhxg, yy);
			title(axhxg, 'Distribution of Paired Differences in Lick Time (stim-unstim)');
			axhxg.XLabel.String = 'Lick Time Difference (s)';
			axhxg.YLabel.String = 'Probability';

			boxplot(axbox, result);
			title(axbox, '5 Number Summary of Differences in Lick Time (stim-unstim)');
			axhxg.YLabel.String = 'Lick Time Difference (s)';


			results.CI95 = [CIl, CIu];
			results.Notes = Notes;
			% 			results.nTrials.stim_then_unstim = numel(stim_then_unstim);
			% 			results.nTrials.unstim_then_stim = numel(unstim_then_stim);
			% 			results.Pools.stim_then_unstim = stim_then_unstim;
			% 			results.Pools.unstim_then_stim = unstim_then_stim;
			results.Pools.lickTimeDifference_stim_minus_unstim = lickTimeDifference_stim_minus_unstim;
		end



		%-----------------------------------------------------------------------------------------------------------
		% 	Helper methods
		%-----------------------------------------------------------------------------------------------------------
		% 
		% 	Helper tools for single-session data
		% 
		function parceExclusions(obj)
			% 
			%  Parse Exclusions
			% 
            Excluded_Trials = [];
			ichar = 1;
            while ichar <= length(obj.iv.excludedtrials_)
                if strcmp(obj.iv.excludedtrials_(ichar),' ')
                    ichar = ichar + 1;
                elseif ismember(obj.iv.excludedtrials_(ichar), '0123456789')
                    jchar = ichar;
                    next_number = '';
                    while jchar <= length(obj.iv.excludedtrials_) && ismember(obj.iv.excludedtrials_(jchar), '0123456789')% get the single numbers eg 495
                        next_number(end+1) = obj.iv.excludedtrials_(jchar);
                        jchar = jchar + 1;
                    end
                    next_number = str2double(next_number);
                    if next_number <= obj.iv.num_trials % otherwise ignore bc is not in range
                        Excluded_Trials(end + 1) = next_number;
                    end
                    ichar = jchar;
                elseif strcmp(obj.iv.excludedtrials_(ichar),'-')
                    while ichar <= length(obj.iv.excludedtrials_) && ~ismember(obj.iv.excludedtrials_(ichar), '0123456789')
                        ichar = ichar + 1;
                    end
                    jchar = ichar;
                    next_number = '';
                    while jchar <= length(obj.iv.excludedtrials_) && ismember(obj.iv.excludedtrials_(jchar), '0123456789')% get the single numbers eg 495
                        next_number(end+1) = obj.iv.excludedtrials_(jchar);
                        jchar = jchar + 1;
                    end
                    next_number = str2double(next_number);
                    if next_number <= obj.iv.num_trials && ~isempty(Excluded_Trials)
                        trials_to_append = (Excluded_Trials(end)+1:next_number);
                        Excluded_Trials = horzcat(Excluded_Trials,trials_to_append);	
                    elseif ~isempty(Excluded_Trials)
                        trials_to_append = (Excluded_Trials(end)+1:obj.iv.num_trials);
                        Excluded_Trials = horzcat(Excluded_Trials,trials_to_append);
                    else
                        warning('Should only reach this line if there''s a dash between two non-numbers')
                    end
                    ichar = ichar;
                else
                    % disp(['parse error: only use numbers, spaces and dashes. you entered: ', obj.iv.excludedtrials_(ichar)])
                    ichar = ichar + 1;
                end
            end
            obj.iv.ExcludedTrials = Excluded_Trials;
		end
		function setConsiderationBound(obj, time_ms)
			obj.bootStat.considerationBoundary_ms = time_ms;
		end
		function changeObjStimType(obj, newProtocol)
			obj.iv.stimtype = newProtocol;
		end
		function save(obj)
			ID = obj.iv.runID;
			if strcmpi(obj.iv.Mode, 'collate')
				savefilename = ['OptoLegacyAnalysis_' obj.iv.collateKey, '_'  datestr(now, 'YYYYmmDD_HH_MM') '_runIDno' num2str(ID)];
			else
				savefilename = ['OptoLegacyObj_' datestr(now, 'YYYYmmDD_HH_MM') '_runIDno' num2str(ID)];
			end
			save([savefilename, '.mat'], 'obj', '-v7.3');
		end
		% 
		% 	Binning licks to reference events
		% 
		function getParsibleLicks(obj)
			if strcmpi(obj.iv.stim_protocol, 'ntj')
				obj.time.s.trialEnd_s = obj.time.s.stimOn_s-0.001;
			else
				obj.time.s.lick_s = obj.iv.init_variables_ChR2.lick_times;
				obj.time.s.cue_s = obj.iv.init_variables_ChR2.cue_on_times;
				obj.time.s.lampOff_s = obj.iv.init_variables_ChR2.trial_start_times;
			end
		end
		function getBinnedLicks(obj, ref, s_b4, s_post, rxnwin_s, nanNeverLick)
			% 
			% 	ref = the reference event for alignment
			% 			'cue'
			% 			'lick'
			% 			'lampOff'
			% 
			% 	s_b4 	= time before ref to consider in SEC
			% 	s_post 	= time after ref to consider in SEC
			% 
			% 	rxnwin)s = permitted rxn window during exp (0 for op0, 0.5 for hyb500/op500)
			% 
			% 	nanNeverLick = logical. Normally if a first lick occurs after the next trial has started, we don't want to call this a "first lick" - set to nan
			% 					NB: in Emery method, we handle these NEVERLICKS differently, so set to zero when call binning from boot > emeryBoot
			% 
			if nargin < 3
				s_b4 = 30000;
			end
			if nargin < 4
				s_post = 30000;
			end
			if nargin < 2
				ref = 'cue';
			end
			if nargin < 5
				rxnwin_s = 0;
			end
			if nargin < 6
				nanNeverLick = 1;
			end
			obj.time.s.binnedLicks.ref = ref;
			obj.time.s.binnedLicks.s_b4 = s_b4;
			obj.time.s.binnedLicks.s_post = s_post;

			% pntrs
			if strcmpi(obj.iv.stim_protocol, 'ntj')
				stimOn 	= obj.time.s.stimOn_s;
				stimOff = obj.time.s.stimOff_s;
				lick = obj.time.s.lick_s;
				refrxOff = obj.time.s.refrxOff_s;
				juice 	= obj.time.s.juice_s;
				if strcmpi(ref, 'stimOn')
					refevents = stimOn;
				elseif strcmpi(ref, 'stimOff')
					refevents = stimOff;
				elseif strcmpi(ref, 'lick')
					refevents = lick;
				elseif strcmpi(ref, 'refrxOff')
					refevents = refrxOff;
				elseif strcmpi(ref, 'juice')
					refevents = juice;
				end
				lick_s_bt_wrtref = cell(numel(refevents), 1);
				f_lick_s_wrtref = nan(numel(refevents), 1);
				stimOn_s_bt_wrtref = cell(numel(refevents), 1);
				stimOff_s_bt_wrtref = cell(numel(refevents), 1);
				refrxOff_s_bt_wrtref = cell(numel(refevents), 1);
				juice_s_bt_wrtref = cell(numel(refevents), 1);
			else
				cue = obj.time.s.cue_s;
				lick = obj.time.s.lick_s;
				lampOff = obj.time.s.lampOff_s;
				if strcmpi(ref, 'cue')
					refevents = cue;
				elseif strcmpi(ref, 'lick')
					refevents = lick;
				elseif strcmpi(ref, 'lampOff')
					refevents = lampOff;
				elseif strcmpi(ref, 'lampOff')
					refevents = lampOff;
				end
				lick_s_bt_wrtref = cell(numel(refevents), 1);
				f_lick_s_wrtref = nan(numel(refevents), 1);
			end
			% 
			for i_ref = 1:numel(refevents)
				if strcmpi(obj.iv.stim_protocol, 'ntj')
					% 
					% 	Find the current trial No
					% 
					trialNo = find(refevents(i_ref) >= obj.time.s.stimOn_s, 1, 'first');
	                [~, lick_s_bt_wrtref{i_ref}, f_lick_s_wrtref(i_ref)] = obj.referenceBinHelper(i_ref, refevents, lick, s_b4, s_post, rxnwin_s, nanNeverLick, trialNo);
	                [~, stimOn_s_bt_wrtref{i_ref}] = obj.referenceBinHelper(i_ref, refevents, stimOn, s_b4, s_post, rxnwin_s, nanNeverLick, trialNo);
	                [~, stimOff_s_bt_wrtref{i_ref}] = obj.referenceBinHelper(i_ref, refevents, stimOff, s_b4, s_post, rxnwin_s, nanNeverLick, trialNo);
	                [~, refrxOff_s_bt_wrtref{i_ref}] = obj.referenceBinHelper(i_ref, refevents, refrxOff, s_b4, s_post, rxnwin_s, nanNeverLick, trialNo);
	                [~, juice_s_bt_wrtref{i_ref}] = obj.referenceBinHelper(i_ref, refevents, juice, s_b4, s_post, rxnwin_s, nanNeverLick, trialNo);

                else
                	[~, lick_s_bt_wrtref{i_ref}, f_lick_s_wrtref(i_ref)] = obj.referenceBinHelper(i_ref, refevents, lick, s_b4, s_post, rxnwin_s, nanNeverLick);
            	end                
			end
			% 
			% 	Handle exclusions!
			% 
			obj.time.s.binnedLicks.lick_s = lick_s_bt_wrtref;
			obj.time.s.binnedLicks.lick_ex_s = lick_s_bt_wrtref;
			

			obj.time.s.binnedLicks.f_lick_s_wrtref = f_lick_s_wrtref;
			obj.time.s.binnedLicks.f_lick_ex_s_wrtref = f_lick_s_wrtref;
			if strcmpi(obj.iv.stim_protocol, 'ntj')
				obj.time.s.binnedLicks.stimOn_s_bt_wrtref = stimOn_s_bt_wrtref;
				obj.time.s.binnedLicks.stimOff_s_bt_wrtref = stimOff_s_bt_wrtref;
				obj.time.s.binnedLicks.refrxOff_s_bt_wrtref = refrxOff_s_bt_wrtref;
				obj.time.s.binnedLicks.juice_s_bt_wrtref = juice_s_bt_wrtref;
			end

			if strcmpi(ref, 'cue') || strcmpi(ref, 'stimOn')
				for iExcl = obj.iv.ExcludedTrials
					obj.time.s.binnedLicks.lick_ex_s{iExcl} = [];
					obj.time.s.binnedLicks.f_lick_ex_s_wrtref(iExcl) = nan;
	            end
            else
            	warning('@getBinnedLicks: exclusions not taken on non cue/stimOn event-aligned binned licks')
            end

			obj.time.s.binnedLicks.refevents = refevents;
			obj.time.s.binnedLicks.rxnwin_s = rxnwin_s;
			obj.time.s.binnedLicks.nanNeverLick = nanNeverLick;
		end		
		function [targetevent_s_bt, targetevent_s_bt_wrtref, first_targetevent_s_bt_wrtref] = referenceBinHelper(obj, i_ref, refevents, targetevents, s_b4, s_post, rxnwin_s, nanNeverLick, trialNo)
			% 
			% 	>> getBinnedLicks >> referenceBinHelper
			% 
			% 	Does the binning wrt reference for any general data set and reference
			% 
			if nargin < 9
				trialNo = [];
			end
			lbe = find(targetevents > refevents(i_ref)-s_b4);
			targetevent_s_bt = targetevents(lbe(ismember(lbe, find(targetevents < refevents(i_ref) + s_post))));
            targetevent_s_bt_wrtref = targetevent_s_bt - refevents(i_ref);
            if ~isempty(targetevent_s_bt) && ~isempty(targetevent_s_bt(find(targetevent_s_bt>refevents(i_ref)+rxnwin_s, 1, 'first')))
                first_targetevent_s_bt_wrtref = targetevent_s_bt(find(targetevent_s_bt>refevents(i_ref)+rxnwin_s, 1, 'first')) - refevents(i_ref);
            else
                first_targetevent_s_bt_wrtref = nan;
            end
            if nanNeverLick
            	if strcmpi(obj.iv.stim_protocol, 'ntj') 
            		if i_ref ~= numel(refevents) && first_targetevent_s_bt_wrtref > obj.time.s.stimOn_s(trialNo+1)
            			first_targetevent_s_bt_wrtref = nan;
        			end
        		elseif first_targetevent_s_bt_wrtref > obj.time.s.total_time_s
                	first_targetevent_s_bt_wrtref = nan;
            	end
        	end
            if ~isempty(targetevent_s_bt)
                assert(targetevent_s_bt(1) > refevents(i_ref)-s_b4 && targetevent_s_bt(end) < refevents(i_ref)+s_post);
                assert(targetevent_s_bt_wrtref(1) > 0-s_b4 && targetevent_s_bt_wrtref(end) < s_post);
            end
		end




		% ----------------------------------------------------------------------------------------------------------
		% 	Methods for COLLATED data
		% ----------------------------------------------------------------------------------------------------------
		% 
		% 	Helper tools:
		% 
		function setCategoryIndicies(obj)
			% 
			% 	Run this to alter the indicies
			% 
			activationidxs = [1:2, 5:45, 69:751];
            inhibidxs = [48:68];
            noopsinidxs = [46,47, 89:108];[46,47, 76:87];
            obj.bootStat.collateHelper.categoryIndicies.activationidxs = activationidxs;
            obj.bootStat.collateHelper.categoryIndicies.inhibidxs = inhibidxs;
            obj.bootStat.collateHelper.categoryIndicies.noopsinidxs = noopsinidxs;
            % 
            % 	Check: {obj.bootStat.collateResult(obj.bootStat.collateHelper.categoryIndicies.inhibidxs).sessionID}
            % 
		end
		function resetCollateHelper(obj)
			obj.bootStat.collateHelper = {};
		end
		% 
		% 	Plot collated stats by SESSION
		% 
		function plotCollateResult(obj, Mode)
			% 
			% 	Plots the median and means of each session from the bootStat.collateResult field. 
			% 	Use with EmeryBoot method in collation
			% 
			% 	Mode: 	box -- plots boxplot
			% 			scatter -- plots scatter
			% 
			if nargin < 2
				Mode = 'scatter';
			end
			obj.setCategoryIndicies;
			activationidxs = obj.bootStat.collateHelper.categoryIndicies.activationidxs;
            inhibidxs = obj.bootStat.collateHelper.categoryIndicies.inhibidxs;
            noopsinidxs = obj.bootStat.collateHelper.categoryIndicies.noopsinidxs;
            idxs = [activationidxs,inhibidxs];
            allidxs = [activationidxs,inhibidxs, noopsinidxs];

            f = figure;
            set(f, 'color', 'white');
            ax1 = subplot(1,2,1);
            hold(ax1, 'on');
            ax1.YLabel.String = 'Median Difference in Movement Time (s)';
            set(ax1, 'fontsize', 20)
            
            ax2 = subplot(1,2,2);
            hold(ax2, 'on');
            ax2.YLabel.String = 'Mean Difference in Movement Time (s)';
            set(ax2, 'fontsize', 20)
            

            if strcmpi(Mode, 'scatter')
				plot(ax1, [1:length(activationidxs)], [obj.bootStat.collateResult(activationidxs).median], 'b.', 'markersize',30)
	            plot(ax1, [length(activationidxs)+1:length(idxs)], [obj.bootStat.collateResult(inhibidxs).median], 'r.', 'markersize',30)
	            plot(ax1, [length(idxs)+1:numel(allidxs)], [obj.bootStat.collateResult(noopsinidxs).median], 'k.', 'markersize',30)
	            yy = get(ax1, 'ylim');
	            xx = get(ax1, 'xlim');
                plot(ax1, xx, [0,0], 'k-')
	            plot(ax1, [length(activationidxs)+.5, length(activationidxs)]+.5, yy, 'k-')
	            plot(ax1, [length(idxs)+.5, length(idxs)+.5], yy, 'k-')
	            ax1.XLabel.String = 'Session';
	            plot(ax2, [1:length(activationidxs)], [obj.bootStat.collateResult(activationidxs).mean], 'b.', 'markersize',30)
	            plot(ax2, [length(activationidxs)+1:length(idxs)], [obj.bootStat.collateResult(inhibidxs).mean], 'r.', 'markersize',30)
	            plot(ax2, [length(idxs)+1:numel(allidxs)], [obj.bootStat.collateResult(noopsinidxs).mean], 'k.', 'markersize',30)
	            yy = get(ax2, 'ylim');
	            xx = get(ax2, 'xlim');
                plot(ax2, xx, [0,0], 'k-')
	            plot(ax2, [length(activationidxs)+.5, length(activationidxs)]+.5, yy, 'k-')
	            plot(ax2, [length(idxs)+.5, length(idxs)+.5], yy, 'k-')
	            ax2.XLabel.String = 'Session';
            elseif strcmpi(Mode, 'box')
            	plot(ax1, ones(numel(activationidxs),1) + (rand(numel(activationidxs),1)-0.5)./2, [obj.bootStat.collateResult(activationidxs).median], 'b.', 'markersize',25)
	            plot(ax1, 2.*ones(numel(noopsinidxs),1) + (rand(numel(noopsinidxs),1)-0.5)./2, [obj.bootStat.collateResult(noopsinidxs).median], 'k.', 'markersize',25)
	            plot(ax1, 3.*ones(numel(inhibidxs),1) + (rand(numel(inhibidxs),1)-0.5)./2, [obj.bootStat.collateResult(inhibidxs).median], 'r.', 'markersize',25)
            	boxplot(ax1, [[obj.bootStat.collateResult(activationidxs).median],[obj.bootStat.collateResult(noopsinidxs).median],[obj.bootStat.collateResult(inhibidxs).median]],[ones(numel(activationidxs),1)', 2.*ones(numel(noopsinidxs),1)', 3.*ones(numel(inhibidxs),1)'])
	            xticks(ax1,[1,2,3]);
	            xticklabels(ax1,{'stimulation', 'no opsin', 'inhibition'});
	            xtickangle(ax1,45);
                plot(ax1, [0,4], [0,0], 'k-')
                plot(ax2, ones(numel(activationidxs),1) + (rand(numel(activationidxs),1)-0.5)./2, [obj.bootStat.collateResult(activationidxs).mean], 'b.', 'markersize',25)
	            plot(ax2, 2.*ones(numel(noopsinidxs),1) + (rand(numel(noopsinidxs),1)-0.5)./2, [obj.bootStat.collateResult(noopsinidxs).mean], 'k.', 'markersize',25)
	            plot(ax2, 3.*ones(numel(inhibidxs),1) + (rand(numel(inhibidxs),1)-0.5)./2, [obj.bootStat.collateResult(inhibidxs).mean], 'r.', 'markersize',25)
	            boxplot(ax2, [[obj.bootStat.collateResult(activationidxs).mean],[obj.bootStat.collateResult(noopsinidxs).mean],[obj.bootStat.collateResult(inhibidxs).mean]],[ones(numel(activationidxs),1)', 2.*ones(numel(noopsinidxs),1)', 3.*ones(numel(inhibidxs),1)'])
	            xticks(ax2,[1,2,3]);
	            xticklabels(ax2,{'stimulation', 'no opsin', 'inhibition'});
	            xtickangle(ax2,45);
                plot(ax2, [0,4], [0,0], 'k-')
        	else
        		error('undefined plot Mode. Use scatter or box')
            end
		end
		% 
		% 	ANOVA to compare activation, inhibition and no opsin mean and median statistics across SESSIONS
		% 
		function anovaCollatedData(obj)
			% 
			% 	>> anovaCollatedData
			% 	Runs 1-way ANOVA to compare mean and median session data between the three categories, activation, inhibition, no opsin
			% 	Tested for EMERY BOOT collation
			% 
			obj.resetCollateHelper;
			obj.setCategoryIndicies;
			% 
			% 	ANOVA1 dataset
			% 
			collatedMedians = [obj.bootStat.collateResult.median]';
			collatedMeans = [obj.bootStat.collateResult.mean]';
			% 
			% 	Category Labels
			% 
			categoryLabels = cell(numel(collatedMedians), 1);
			categoryLabels(obj.bootStat.collateHelper.categoryIndicies.activationidxs, 1) = {'activation'};
			categoryLabels(obj.bootStat.collateHelper.categoryIndicies.inhibidxs) = {'inhibition'};
			categoryLabels(obj.bootStat.collateHelper.categoryIndicies.noopsinidxs) = {'no opsin'};
            %             
            %   If we have omitted any of the collated data, we should
            %   handle that now
            % 
            killIdxs = cellfun(@(x) isempty(x), categoryLabels);
            categoryLabels(killIdxs) = [];
            collatedMedians(killIdxs) = [];
            collatedMeans(killIdxs) = [];
			% 
			% 	Run ANOVA-1
			% 		H0: mean(activation) = mean(no opsin) = mean(inhibition)
			% 		H1: means are not equal for all categories
			% 
			[obj.bootStat.collateHelper.ANOVA1.median.p,obj.bootStat.collateHelper.ANOVA1.median.tbl,obj.bootStat.collateHelper.ANOVA1.median.stats] = anova1(collatedMedians,categoryLabels);
			[obj.bootStat.collateHelper.ANOVA1.mean.p,obj.bootStat.collateHelper.ANOVA1.mean.tbl,obj.bootStat.collateHelper.ANOVA1.mean.stats] = anova1(collatedMeans,categoryLabels);
		end
		function anovaCollatedData_compareMean(obj)
			% 
			%  Run first:
			%  >> anovaCollatedData
			% 	Then run this function to compare the between-category mean of session means
			% 	>> anovaCollatedData_compareMean 
			% --------------------------------------------------------------
			% 
			% 	Handle multiple comparisons to test which mean is different
			% 
			[obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.results, obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.mu_sem] = multcompare(obj.bootStat.collateHelper.ANOVA1.mean.stats);
			obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.note{1,1} = 'group 1';
			obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.note{1,2} = 'group 2';
			obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.note{1,6} = 'p';
			obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.GroupIDs{1,1} = obj.bootStat.collateHelper.ANOVA1.mean.stats.gnames{1};
			obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.GroupIDs{2,1} = obj.bootStat.collateHelper.ANOVA1.mean.stats.gnames{2};
			obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.GroupIDs{3,1} = obj.bootStat.collateHelper.ANOVA1.mean.stats.gnames{3};
			obj.printANOVA1Summary('mean');
		end
		function anovaCollatedData_compareMedian(obj)
			% 
			%  Run first:
			%  >> anovaCollatedData
			% 	Then run this function to compare the between-category mean of session medians
			% 	>> anovaCollatedData_compareMedian 
			% --------------------------------------------------------------
			% 
			% 	Handle multiple comparisons to test which mean is different
			% 
			[obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.results, obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.mu_sem] = multcompare(obj.bootStat.collateHelper.ANOVA1.median.stats);
			obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.note{1,1} = 'group 1';
			obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.note{1,2} = 'group 2';
			obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.note{1,6} = 'p';
			obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.GroupIDs{1,1} = obj.bootStat.collateHelper.ANOVA1.median.stats.gnames{1};
			obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.GroupIDs{2,1} = obj.bootStat.collateHelper.ANOVA1.median.stats.gnames{2};
			obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.GroupIDs{3,1} = obj.bootStat.collateHelper.ANOVA1.median.stats.gnames{3};
			obj.printANOVA1Summary('median');
		end
		function printANOVA1Summary(obj, style)
			if nargin < 2
				style = 'median';
			end
			if strcmpi(style, 'median')
				p = obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.results(:,6);
				groupNames = obj.bootStat.collateHelper.ANOVA1.median.stats.gnames;
				mu_sem = obj.bootStat.collateHelper.ANOVA1.median.multipleComparisons.mu_sem;
			elseif strcmpi(style, 'mean')
				p = obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.results(:,6);
				groupNames = obj.bootStat.collateHelper.ANOVA1.mean.stats.gnames;
				mu_sem = obj.bootStat.collateHelper.ANOVA1.mean.multipleComparisons.mu_sem;
			else
				error('undefined style');
			end
			disp('-------- Multiple Comparisons ANOVA1 Result --------')
			disp(['	STYLE:' style])
			disp('	GROUP IDs, mean and SEM')
			disp(['		1. ' groupNames{1} ': ' num2str(mu_sem(1,1)) ' (' num2str(mu_sem(1,1)-mu_sem(1,2)) ', ' num2str(mu_sem(1,1)+mu_sem(1,2)) ')'])
			disp(['		2. ' groupNames{2} ': ' num2str(mu_sem(2,1)) ' (' num2str(mu_sem(2,1)-mu_sem(2,2)) ', ' num2str(mu_sem(2,1)+mu_sem(2,2)) ')'])
			disp(['		3. ' groupNames{3} ': ' num2str(mu_sem(3,1)) ' (' num2str(mu_sem(3,1)-mu_sem(3,2)) ', ' num2str(mu_sem(3,1)+mu_sem(3,2)) ')'])
			disp(' ')
			disp('	P-values:')
			disp(['		1. ' groupNames{1} ' ~= ' groupNames{2} ': ' num2str(p(1))])
			disp(['		2. ' groupNames{1} ' ~= ' groupNames{3} ': ' num2str(p(2))])
			disp(['		3. ' groupNames{2} ' ~= ' groupNames{3} ': ' num2str(p(3))])
			
		end
		% 
		% 	Obsolete method to try comparing session results with bootstrap
		%
		function bootCollatedMedians(obj, nboot)
			% 
			% 	THIS IS FOR CI ON EACH INDIVIDUAL GROUP
			% 
			% 	UPDATE: Fixed! This works great now!
			% 
			obj.setCategoryIndicies;
			activationidxs = obj.bootStat.collateHelper.categoryIndicies.activationidxs;
            inhibidxs = obj.bootStat.collateHelper.categoryIndicies.inhibidxs;
            noopsinidxs = obj.bootStat.collateHelper.categoryIndicies.noopsinidxs;

			obj.resetCollateHelper
			if nargin < 2
				nboot = 1000000;
			end
			obj.bootStat.collateHelper.Note{1,1} = ['Used bootCollatedMedians(obj, nboot=' num2str(nboot) ')'];
			% 
			% 	Boot the Collated Pooled Data
			% 
			activationSets.median = [obj.bootStat.collateResult(activationidxs).median];
			activationSets.mean = [obj.bootStat.collateResult(activationidxs).mean];
			[obj.bootStat.collateHelper.activation.median.result] = obj.nonParametricBoot(activationSets.median, nboot);
			[obj.bootStat.collateHelper.activation.mean.result] = obj.nonParametricBoot(activationSets.mean, nboot);

			inhibitionsSets.median = [obj.bootStat.collateResult(inhibidxs).median];
			inhibitionsSets.mean = [obj.bootStat.collateResult(inhibidxs).mean];
			[obj.bootStat.collateHelper.inhibition.median.result] = obj.nonParametricBoot(inhibitionsSets.median, nboot);
			[obj.bootStat.collateHelper.inhibition.mean.result] = obj.nonParametricBoot(inhibitionsSets.mean, nboot);

			noOpsinSets.median = [obj.bootStat.collateResult(noopsinidxs).median];
			noOpsinSets.mean = [obj.bootStat.collateResult(noopsinidxs).mean];
			[obj.bootStat.collateHelper.noOpsin.median.result] = obj.nonParametricBoot(noOpsinSets.median, nboot);
			[obj.bootStat.collateHelper.noOpsin.mean.result] = obj.nonParametricBoot(noOpsinSets.mean, nboot);
		end
		function result = nonParametricBoot(obj, X, nboot)
			% 
			% 	>> bootCollatedMedians >> nonParametricBoot
			% 
			% 	Samples a new dataset with replacement of same size as the original dataset
			% 	reports the boot mean, median, and CI on each
			% 
			%  row = Xb(i), col = boot
			Xb = nan(numel(X), nboot);
			for iboot = 1:nboot
				Xb(:,iboot) = X(randi(numel(X),[numel(X),1]));
			end
			result.mean.Xb = sort(mean(Xb,1));
			result.mean.CI95 = [result.mean.Xb(round(0.025*nboot)), result.mean.Xb(round(0.975*nboot))];

			result.median.Xb = sort(median(Xb,1));
			result.median.CI95 = [result.median.Xb(round(0.025*nboot)), result.median.Xb(round(0.975*nboot))];			
		end	
		function plotBootMedian(obj, Mode)
			% 
			% 	>> bootCollatedMedians >> plotBootMedian
			% 
			% 	Mode = 'mean' or 'median' as the bootstrapped statistic
			% 
			if nargin < 2
				Mode = 'Median';
			end
			obj.setCategoryIndicies;
			activationidxs = obj.bootStat.collateHelper.categoryIndicies.activationidxs;
            inhibidxs = obj.bootStat.collateHelper.categoryIndicies.inhibidxs;
            noopsinidxs = obj.bootStat.collateHelper.categoryIndicies.noopsinidxs;
			if strcmpi(Mode, 'Mean')
				result.activation = obj.bootStat.collateHelper.activation.mean.result;
				result.inhibition = obj.bootStat.collateHelper.inhibition.mean.result;
				result.noOpsin = obj.bootStat.collateHelper.noOpsin.mean.result;
				empiricalStat.activation = mean([obj.bootStat.collateResult(activationidxs).mean]);
				empiricalStat.inhibition = mean([obj.bootStat.collateResult(inhibidxs).mean]);
				empiricalStat.noOpsin = mean([obj.bootStat.collateResult(noopsinidxs).mean]);
				empiricalStat.X.activation = [obj.bootStat.collateResult(activationidxs).mean];
				empiricalStat.X.inhibition = [obj.bootStat.collateResult(inhibidxs).mean];
				empiricalStat.X.noOpsin = [obj.bootStat.collateResult(noopsinidxs).mean];
			elseif strcmpi(Mode, 'Median')
				result.activation = obj.bootStat.collateHelper.activation.median.result;
				result.inhibition = obj.bootStat.collateHelper.inhibition.median.result;
				result.noOpsin = obj.bootStat.collateHelper.noOpsin.median.result;
				empiricalStat.activation = mean([obj.bootStat.collateResult(activationidxs).median]);
				empiricalStat.inhibition = mean([obj.bootStat.collateResult(inhibidxs).median]);
				empiricalStat.noOpsin = mean([obj.bootStat.collateResult(noopsinidxs).median]);
				empiricalStat.X.activation = [obj.bootStat.collateResult(activationidxs).median];
				empiricalStat.X.inhibition = [obj.bootStat.collateResult(inhibidxs).median];
				empiricalStat.X.noOpsin = [obj.bootStat.collateResult(noopsinidxs).median];
			end
			% 
			% 	Plot Histogram --- Boot of Medians
			% 
			figure
			ax4 = subplot(1,3,1);
			ax5 = subplot(1,3,2);
			ax6 = subplot(1,3,3);
			hold(ax4, 'on');
			hold(ax5, 'on');
			hold(ax6, 'on');
			set(ax4, 'fontsize', 20);
			set(ax5, 'fontsize', 20);
			set(ax6, 'fontsize', 20);
			histogram(ax4, result.activation.mean.Xb, 'normalization', 'probability', 'displayStyle', 'stairs','LineWidth', 2,  'displayName', ['Pool n=' num2str(numel(empiricalStat.X.activation))])
			xx = get(ax4, 'xlim');
			yy = get(ax4, 'ylim');
			plot(ax4, [result.activation.mean.CI95(1), result.activation.mean.CI95(1)], yy, 'r-', 'LineWidth', 2, 'DisplayName', ['2.5pct: ' num2str(round(result.activation.mean.CI95(1),2))]);
			plot(ax4, [result.activation.mean.CI95(2), result.activation.mean.CI95(2)], yy, 'r-', 'LineWidth', 2, 'DisplayName', ['97.5pct: ' num2str(round(result.activation.mean.CI95(2),2))]);
			plot(ax4, [empiricalStat.activation,empiricalStat.activation], yy, 'b-', 'LineWidth', 2, 'DisplayName', ['Empirical ' Mode ': ' num2str(empiricalStat.activation)]);
			title(ax4, ['Boot Mean of Activation ' Mode])
			% legend(ax4, 'show')

			histogram(ax5, result.inhibition.mean.Xb, 'LineWidth', 2, 'normalization', 'probability', 'displayStyle', 'stairs', 'displayName', ['Pool n=' num2str(numel(empiricalStat.X.inhibition))])
			plot(ax5, [result.inhibition.mean.CI95(1), result.inhibition.mean.CI95(1)], yy, 'r-', 'LineWidth', 2, 'DisplayName', ['2.5pct: ' num2str(round(result.inhibition.mean.CI95(1),2))]);
			plot(ax5, [result.inhibition.mean.CI95(2), result.inhibition.mean.CI95(2)], yy, 'r-', 'LineWidth', 2, 'DisplayName', ['97.5pct: ' num2str(round(result.inhibition.mean.CI95(2),2))]);
			plot(ax5, [empiricalStat.inhibition,empiricalStat.inhibition], yy, 'b-', 'LineWidth', 2, 'DisplayName', ['Empirical ' Mode ': ' num2str(empiricalStat.inhibition)]);
			title(ax5, ['Boot Mean of inhibition ' Mode]);
			% legend(ax5, 'show')

			histogram(ax6, result.noOpsin.mean.Xb, 'LineWidth', 2, 'normalization', 'probability', 'displayStyle', 'stairs', 'displayName', ['Pool n=' num2str(numel(empiricalStat.X.noOpsin))])
			plot(ax6, [result.noOpsin.mean.CI95(1), result.noOpsin.mean.CI95(1)], yy, 'r-', 'LineWidth', 2, 'DisplayName', ['2.5pct: ' num2str(round(result.noOpsin.mean.CI95(1),2))]);
			plot(ax6, [result.noOpsin.mean.CI95(2), result.noOpsin.mean.CI95(2)], yy, 'r-', 'LineWidth', 2, 'DisplayName', ['97.5pct: ' num2str(round(result.noOpsin.mean.CI95(2),2))]);
			plot(ax6, [empiricalStat.noOpsin,empiricalStat.noOpsin], yy, 'b-', 'LineWidth', 2, 'DisplayName', ['Empirical ' Mode ': ' num2str(empiricalStat.noOpsin)]);
			title(ax6, ['Boot Mean of No Opsin ' Mode])
			% legend(ax6, 'show')
			disp('----------------------------------------')
			disp(['	Nonparametric Bootstrap of ' Mode ' Results:'])
			disp(['		Activation:'])
			disp(['			nSessions: ' num2str(numel(empiricalStat.X.activation))])
			disp(['			sample mean of session ' Mode 's: ' num2str(round(empiricalStat.activation, 3)) '	bootCI95: (' num2str(round(result.activation.mean.CI95(1),2)) ',' num2str(round(result.activation.mean.CI95(2),2)) ')'])
			disp(['		Inhibition:'])
			disp(['			nSessions: ' num2str(numel(empiricalStat.X.inhibition))])
			disp(['			sample mean of session ' Mode 's: ' num2str(round(empiricalStat.inhibition, 3)) '	bootCI95: (' num2str(round(result.inhibition.mean.CI95(1),2)) ',' num2str(round(result.inhibition.mean.CI95(2),2)) ')'])
			disp(['		No Opsin:'])
			disp(['			nSessions: ' num2str(numel(empiricalStat.X.noOpsin))])
			disp(['			sample mean of session ' Mode 's: ' num2str(round(empiricalStat.noOpsin, 3)) '	bootCI95: (' num2str(round(result.noOpsin.mean.CI95(1),2)) ',' num2str(round(result.noOpsin.mean.CI95(2),2)) ')'])
		end 
		function bootCollatedMediansComparison(obj, nboot)
			% 
			% 	THIS IS FOR CI ON between group comparisons!
			% 
			%  bootCollatedMediansComparison > nonParametricBoot2GroupComparison > plotBootComparison
			% 
			obj.setCategoryIndicies;
			activationidxs = obj.bootStat.collateHelper.categoryIndicies.activationidxs;
            inhibidxs = obj.bootStat.collateHelper.categoryIndicies.inhibidxs;
            noopsinidxs = obj.bootStat.collateHelper.categoryIndicies.noopsinidxs;

			obj.resetCollateHelper
			if nargin < 2
				nboot = 1000000;
			end
			obj.bootStat.collateHelper.Note{1,1} = ['Used bootCollatedMediansComparison(obj, nboot=' num2str(nboot) ')'];
			% 
			% 	Boot the Collated Pooled Data
			% 
			activationSets.median = [obj.bootStat.collateResult(activationidxs).median];
			activationSets.mean = [obj.bootStat.collateResult(activationidxs).mean];
			inhibitionsSets.median = [obj.bootStat.collateResult(inhibidxs).median];
			inhibitionsSets.mean = [obj.bootStat.collateResult(inhibidxs).mean];
			noOpsinSets.median = [obj.bootStat.collateResult(noopsinidxs).median];
			noOpsinSets.mean = [obj.bootStat.collateResult(noopsinidxs).mean];
			% 
			% 	Compare activation to no-opsin
			% 
			disp('bootCollatedMediansComparison > stim vs no opsin median...')
			[obj.bootStat.collateHelper.activationVnoOpsin.median.result] = obj.nonParametricBoot2GroupComparison(activationSets.median, noOpsinSets.median, nboot);
			disp('bootCollatedMediansComparison > stim vs no opsin mean...')
			[obj.bootStat.collateHelper.activationVnoOpsin.mean.result] = obj.nonParametricBoot2GroupComparison(activationSets.mean, noOpsinSets.mean, nboot);
			
			disp('bootCollatedMediansComparison > inhib vs no opsin median...')
			[obj.bootStat.collateHelper.inhibitionVnoOpsin.median.result] = obj.nonParametricBoot2GroupComparison(inhibitionsSets.median, noOpsinSets.median, nboot);
			disp('bootCollatedMediansComparison > inhib vs no opsin mean...')
			[obj.bootStat.collateHelper.inhibitionVnoOpsin.mean.result] = obj.nonParametricBoot2GroupComparison(inhibitionsSets.mean, noOpsinSets.mean, nboot);

			disp('bootCollatedMediansComparison > inhib vs stim median...')
			[obj.bootStat.collateHelper.inhibitionVactivation.median.result] = obj.nonParametricBoot2GroupComparison(inhibitionsSets.median, activationSets.median, nboot);
			disp('bootCollatedMediansComparison > inhib vs stim mean...')
			[obj.bootStat.collateHelper.inhibitionVactivation.mean.result] = obj.nonParametricBoot2GroupComparison(inhibitionsSets.mean, activationSets.mean, nboot);


		end 
		function bootCollatedAUCComparison(obj, nboot)
			% 
			% 	THIS IS FOR CI ON between group comparisons!
			% 
			%  BOOTSTRAPS THE dAUCS IN EACH CATEGORY
			% 
			activationidxs = [find(strcmpi({obj.bootStat.collateResult.type},'ChR2')), find(strcmpi({obj.bootStat.collateResult.type},'tapers')), find(strcmpi({obj.bootStat.collateResult.type},'ChrimsonR'))];
            inhibidxs = find(strcmpi({obj.bootStat.collateResult.type},'gtacr2'));
            noopsinidxs = find(strcmpi({obj.bootStat.collateResult.type},'no opsin'));

			obj.resetCollateHelper
			if nargin < 2
				nboot = 1000000;
			end
			obj.bootStat.collateHelper.Note{1,1} = ['Used bootCollatedAUCComparison(obj, nboot=' num2str(nboot) ')'];
			% 
			% 	Boot the Collated Pooled Data
			% 
			activationSets.dAUC = [obj.bootStat.collateResult(activationidxs).delAUC_EOT];
			inhibitionsSets.dAUC = [obj.bootStat.collateResult(inhibidxs).delAUC_EOT];
			noOpsinSets.dAUC = [obj.bootStat.collateResult(noopsinidxs).delAUC_EOT];
			% 
			% 	Compare activation to no-opsin
			% 
			disp('bootCollatedAUCComparison > stim vs no opsin dAUC...')
			[obj.bootStat.collateHelper.activationVnoOpsin.delAUC_EOT.result] = obj.nonParametricBoot2GroupComparison(activationSets.dAUC, noOpsinSets.dAUC, nboot);
			
			disp('bootCollatedAUCComparison > inhib vs no opsin dAUC...')
			[obj.bootStat.collateHelper.inhibitionVnoOpsin.delAUC_EOT.result] = obj.nonParametricBoot2GroupComparison(inhibitionsSets.dAUC, noOpsinSets.dAUC, nboot);
			
			disp('bootCollatedAUCComparison > inhib vs stim dAUC...')
			[obj.bootStat.collateHelper.inhibitionVactivation.delAUC_EOT.result] = obj.nonParametricBoot2GroupComparison(inhibitionsSets.dAUC, activationSets.dAUC, nboot);
		end 
		function result = nonParametricBoot2GroupComparison(obj, X1, X2, nboot)
			% 
			% 	>> nonParametricBoot2GroupComparison >> nonParametricBoot
			% 
			% 	Samples a new dataset with replacement of same size as the original dataset
			% 	reports the boot mean, median, and CI on each
			% 
			%  row = Xb(i), col = boot
			delX = nan(nboot,1);
			for iboot = 1:nboot
				delX(iboot) = mean(X1(randi(numel(X1),[numel(X1),1]))) - mean(X2(randi(numel(X2),[numel(X2),1])));
			end
			result.meanX1_minus_meanX2 = mean(X1) - mean(X2);
			result.delX = sort(delX);
			result.CI95 = [result.delX(round(0.025*nboot)), result.delX(round(0.975*nboot))];
		end	
		function plotBootComparison(obj, Mode)
			% 
			% 		Set Mode = 'emery' or 'bootAUC'
			% 
			if nargin < 2
				Mode = obj.iv.collateKey;
			end
			

			f = figure;
			set(f,'color', 'w');
			
			if strcmpi(Mode, 'emery')
				if ~isfield(obj.bootStat.collateHelper, 'activationVnoOpsin')
					warning(' Need to calculate the bootCollatedMediansComparison(obj, nboot). Running now')
				end
				ax1 = subplot(1,2,1);
				hold(ax1, 'on');
				title(ax1, 'Median 2-Group Comparison')
				obj.plotMeanWithCI(obj.bootStat.collateHelper.activationVnoOpsin.median.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.activationVnoOpsin.median.result.CI95, 1, 'b', ax1);
				obj.plotMeanWithCI(obj.bootStat.collateHelper.inhibitionVnoOpsin.median.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.inhibitionVnoOpsin.median.result.CI95, 2, 'r', ax1);
				obj.plotMeanWithCI(obj.bootStat.collateHelper.inhibitionVactivation.median.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.inhibitionVactivation.median.result.CI95, 3, 'k', ax1);
				plot(ax1, [0,4],[0,0], 'k-')
				xticks(ax1,[1,2,3])
				xticklabels(ax1,{'stim v no opsin', 'inhib v no opsin', 'inhib v stim'})
				xtickangle(ax1,45)
				set(ax1,'fontsize',30)

				ax2 = subplot(1,2,2);
				hold(ax2, 'on');
				title(ax2, 'Mean 2-Group Comparison')
				obj.plotMeanWithCI(obj.bootStat.collateHelper.activationVnoOpsin.mean.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.activationVnoOpsin.mean.result.CI95, 1, 'b', ax2);
				obj.plotMeanWithCI(obj.bootStat.collateHelper.inhibitionVnoOpsin.mean.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.inhibitionVnoOpsin.mean.result.CI95, 2, 'r', ax2);
				obj.plotMeanWithCI(obj.bootStat.collateHelper.inhibitionVactivation.mean.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.inhibitionVactivation.mean.result.CI95, 3, 'k', ax2);
				plot(ax2, [0,4],[0,0], 'k-')
				xticks(ax2,[1,2,3])
				xticklabels(ax2,{'stim v no opsin', 'inhib v no opsin', 'inhib v stim'})
				xtickangle(ax2,45)
				set(ax2,'fontsize',30)
			elseif strcmpi(Mode, 'bootAUC')
				ax1 = subplot(1,1,1);
				hold(ax1, 'on');
				title(ax1, 'dAUC 2-Group Comparison')
				obj.plotMeanWithCI(obj.bootStat.collateHelper.activationVnoOpsin.delAUC_EOT.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.activationVnoOpsin.delAUC_EOT.result.CI95, 1, 'b', ax1);
				obj.plotMeanWithCI(obj.bootStat.collateHelper.inhibitionVnoOpsin.delAUC_EOT.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.inhibitionVnoOpsin.delAUC_EOT.result.CI95, 2, 'r', ax1);
				obj.plotMeanWithCI(obj.bootStat.collateHelper.inhibitionVactivation.delAUC_EOT.result.meanX1_minus_meanX2, obj.bootStat.collateHelper.inhibitionVactivation.delAUC_EOT.result.CI95, 3, 'k', ax1);
				plot(ax1, [0,4],[0,0], 'k-')
				xticks(ax1,[1,2,3])
				xticklabels(ax1,{'stim v no opsin', 'inhib v no opsin', 'inhib v stim'})
				xtickangle(ax1,45)
				set(ax1,'fontsize',30)
			end
		end
		function plotMeanWithCI(obj, Mean, CI, xpos, Color,ax)
			plot(ax,xpos, Mean, [Color 'o'], 'markerfacecolor', Color, 'markersize', 20)
			plot(ax,[xpos,xpos], [CI(1),Mean], [Color '-'], 'LineWidth', 5)
			plot(ax,[xpos,xpos], [Mean,CI(2)], [Color '-'], 'LineWidth', 5)
			plot(ax,[xpos-0.25,xpos+0.25], [CI(1),CI(1)], [Color '-'], 'LineWidth', 5)
			plot(ax,[xpos-0.25,xpos+0.25], [CI(2),CI(2)], [Color '-'], 'LineWidth', 5)
		end
		% 
		% 	Methods to do stat analysis on COLLATED data -- i.e., considering ALL session data as one session
		% 
		function poolCollatedData(obj)
			% 
			% 	>> bootCollatedData >> poolCollatedData
			% 
			% 	Also can run poolCollatedData on its own to get summary plots
			% 
			% 	Pools the stim and unstim flicks within the category across sessions of collated data. 
			% 	Runs for each cateory and stores data in obj.bootStat.collateHelper for easy access to bootstrap functions
			% 	Then it plots the summary of all data.
			% 	Use for a COMPOSITE FIGURE of opto results for each category
			%
			activationidxs = obj.bootStat.collateHelper.categoryIndicies.activationidxs;
            inhibidxs = obj.bootStat.collateHelper.categoryIndicies.inhibidxs;
            noopsinidxs = obj.bootStat.collateHelper.categoryIndicies.noopsinidxs;			
			% 
			% 	Collect the arrays
			% 
			obj.getPoolsAcrossSesh(activationidxs, inhibidxs, noopsinidxs);
			% 
			% 	Plot some hxg:
			% 
			f = figure;
			ax1 = subplot(2,3,1);
			ax2 = subplot(2,3,2);
			ax3 = subplot(2,3,3);
			hold(ax1, 'on');
			hold(ax2, 'on');
			hold(ax3, 'on');
			histogram(ax1, obj.bootStat.collateHelper.activation.unstimpool, 'normalization', 'probability', 'displayName', ['Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.activation.unstimpool))])
			histogram(ax1, obj.bootStat.collateHelper.activation.stimpool, 'normalization', 'probability', 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.activation.stimpool))])
			title(ax1, ['Activation'])
			legend(ax1, 'show')

			histogram(ax2, obj.bootStat.collateHelper.inhibition.unstimpool, 'normalization', 'probability', 'displayName', ['Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.inhibition.unstimpool))])
			histogram(ax2, obj.bootStat.collateHelper.inhibition.stimpool, 'normalization', 'probability', 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.inhibition.stimpool))])
			title(ax2, ['inhibition'])
			legend(ax2, 'show')

			histogram(ax3, obj.bootStat.collateHelper.noOpsin.unstimpool, 'normalization', 'probability', 'displayName', ['Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.noOpsin.unstimpool))])
			histogram(ax3, obj.bootStat.collateHelper.noOpsin.stimpool, 'normalization', 'probability', 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.noOpsin.stimpool))])
			title(ax3, ['No Opsin'])
			legend(ax3, 'show')
			% 
			% 	Calculate some ecdfs
			% 
			[a_f_s, a_x_s] = ecdf(obj.bootStat.collateHelper.activation.stimpool);
			[a_f_us, a_x_us] = ecdf(obj.bootStat.collateHelper.activation.unstimpool);
			[i_f_s, i_x_s] = ecdf(obj.bootStat.collateHelper.inhibition.stimpool);
			[i_f_us, i_x_us] = ecdf(obj.bootStat.collateHelper.inhibition.unstimpool);
			[no_f_s, no_x_s] = ecdf(obj.bootStat.collateHelper.noOpsin.stimpool);
			[no_f_us, no_x_us] = ecdf(obj.bootStat.collateHelper.noOpsin.unstimpool);
			ax4 = subplot(2,3,4);
			ax5 = subplot(2,3,5);
			ax6 = subplot(2,3,6);
			hold(ax4, 'on');
			hold(ax5, 'on');
			hold(ax6, 'on');
			plot(ax4, a_x_us, a_f_us, 'displayName', ['Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.activation.unstimpool))])
			plot(ax4, a_x_s, a_f_s, 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.activation.stimpool))])
			title(ax4, ['Activation'])
			legend(ax4, 'show')

			
			plot(ax5, i_x_us, i_f_us, 'displayName', ['Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.inhibition.unstimpool))])
			plot(ax5, i_x_s, i_f_s, 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.inhibition.stimpool))])
			title(ax5, ['inhibition'])
			legend(ax5, 'show')

			plot(ax6, no_x_us, no_f_us, 'displayName', ['Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.noOpsin.unstimpool))])
			plot(ax6, no_x_s, no_f_s, 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.noOpsin.stimpool))])
			title(ax6, ['No Opsin'])
			legend(ax6, 'show')
			% 
			% 	Make sure there's no missing data!
			% 
			assert(sum(isnan(obj.bootStat.collateHelper.activation.unstimpool))==0 && sum(isnan(obj.bootStat.collateHelper.activation.stimpool))==0 && sum(isnan(obj.bootStat.collateHelper.inhibition.stimpool))==0 && sum(isnan(obj.bootStat.collateHelper.inhibition.unstimpool))==0 && sum(isnan(obj.bootStat.collateHelper.noOpsin.stimpool))==0 && sum(isnan(obj.bootStat.collateHelper.noOpsin.unstimpool))==0);			
		end
		function kstestCombinedDatasets(obj)
			if ~isfield(obj.bootStat.collateHelper.activation, 'stimpool')
				obj.poolCollatedData;
			end
			[obj.bootStat.KScomposite.stim.h,obj.bootStat.KScomposite.stim.p,obj.bootStat.KScomposite.stim.ks2stat] = kstest2(obj.bootStat.collateHelper.activation.stimpool,obj.bootStat.collateHelper.activation.unstimpool);
			[obj.bootStat.KScomposite.unstim.h,obj.bootStat.KScomposite.unstim.p,obj.bootStat.KScomposite.unstim.ks2stat] = kstest2(obj.bootStat.collateHelper.inhibition.stimpool,obj.bootStat.collateHelper.inhibition.unstimpool);
			[obj.bootStat.KScomposite.noOpsin.h,obj.bootStat.KScomposite.noOpsin.p,obj.bootStat.KScomposite.noOpsin.ks2stat] = kstest2(obj.bootStat.collateHelper.noOpsin.stimpool,obj.bootStat.collateHelper.noOpsin.unstimpool);
		end
		function getPoolsAcrossSesh(obj, activationidxs, inhibidxs, noopsinidxs)
			% 
			% 	>> poolCollatedData >>
			% 
			% 	Pools the stim and unstim flicks within the category across sessions of collated data. 
			% 	Runs for each cateory and stores data in obj.bootStat.collateHelper for easy access to bootstrap functions
			% 
			if isfield(obj.bootStat.collateHelper, 'Note')
				obj.bootStat.collateHelper.Note{end+1, 1} = {'Used obj.getPoolsAcrossSesh to get pooled collated flicks for each stimulation case and stim/unstim categories'};
			else
				obj.bootStat.collateHelper.Note = {'Used obj.getPoolsAcrossSesh to get pooled collated flicks for each stimulation case and stim/unstim categories'};
			end
			[activationPoolNumel.nstim, activationPoolNumel.nunstim] = obj.getNumelOfPoolAcrossSesh(activationidxs, obj.bootStat.collateResult(activationidxs));
			[inhibitionPoolNumel.nstim, inhibitionPoolNumel.nunstim] = obj.getNumelOfPoolAcrossSesh(inhibidxs, obj.bootStat.collateResult(inhibidxs));
			[noOpsinPoolNumel.nstim, noOpsinPoolNumel.nunstim] = obj.getNumelOfPoolAcrossSesh(noopsinidxs, obj.bootStat.collateResult(noopsinidxs));
			obj.bootStat.collateHelper.numelPool.activation.nstim = activationPoolNumel.nstim;
			obj.bootStat.collateHelper.numelPool.activation.nunstim = activationPoolNumel.nunstim;
			obj.bootStat.collateHelper.numelPool.inhibition.nstim = inhibitionPoolNumel.nstim;
			obj.bootStat.collateHelper.numelPool.inhibition.nunstim = inhibitionPoolNumel.nunstim;
			obj.bootStat.collateHelper.numelPool.noOpsin.nstim = noOpsinPoolNumel.nstim;
			obj.bootStat.collateHelper.numelPool.noOpsin.nunstim = noOpsinPoolNumel.nunstim;
			% 
			% 	Fill in
			% 
			[activation.stimpool, activation.unstimpool] = obj.getPool(obj.bootStat.collateResult(activationidxs), activationidxs, obj.bootStat.collateHelper.numelPool.activation.nstim, obj.bootStat.collateHelper.numelPool.activation.nunstim)
			[inhibition.stimpool, inhibition.unstimpool] = obj.getPool(obj.bootStat.collateResult(inhibidxs), inhibidxs, obj.bootStat.collateHelper.numelPool.inhibition.nstim, obj.bootStat.collateHelper.numelPool.inhibition.nunstim)
			[noOpsin.stimpool, noOpsin.unstimpool] = obj.getPool(obj.bootStat.collateResult(noopsinidxs), noopsinidxs, obj.bootStat.collateHelper.numelPool.noOpsin.nstim, obj.bootStat.collateHelper.numelPool.noOpsin.nunstim)

			obj.bootStat.collateHelper.activation.stimpool = activation.stimpool;
			obj.bootStat.collateHelper.activation.unstimpool = activation.unstimpool;
			obj.bootStat.collateHelper.inhibition.stimpool = inhibition.stimpool;
			obj.bootStat.collateHelper.inhibition.unstimpool = inhibition.unstimpool;
			obj.bootStat.collateHelper.noOpsin.stimpool = noOpsin.stimpool;
			obj.bootStat.collateHelper.noOpsin.unstimpool = noOpsin.unstimpool;
		end
		function [nstim, nunstim] = getNumelOfPoolAcrossSesh(obj, idxs, seshpool)
			% 
			% 	>> poolCollatedData >> getPoolsAcrossSesh >>
			% 	Gets numel of flicks in each category across sessions
			%
			nstim = 0;
			nunstim = 0;
			
			for iSession = 1:numel(idxs)
				nstim = nstim + numel(seshpool(iSession).bootStat.resultLick.Pools.stim_pool);
				nunstim = nunstim + numel(seshpool(iSession).bootStat.resultLick.Pools.unstim_pool);
			end			
		end
		function [stimpool, unstimpool] = getPool(obj, seshpool, idxs, nstim, nunstim)
			% 
			% 	>> poolCollatedData >> getPoolsAcrossSesh >>
			% 	Pools the stim and unstim flicks within the category across sessions of collated data. 
			% 
			stimpool = nan(nstim,1);
			unstimpool = nan(nunstim,1);

			stimidx = 1;
			unstimidx = 1;
			
			for iSession = 1:numel(idxs)
				currentpool = seshpool(iSession).bootStat.resultLick.Pools.stim_pool
				idx2 = stimidx + numel(currentpool) - 1;
				stimpool(stimidx:idx2) = currentpool;
				stimidx = idx2+1;
				
				currentpool = seshpool(iSession).bootStat.resultLick.Pools.unstim_pool;
				idx2 = unstimidx + numel(currentpool) - 1;
				unstimpool(unstimidx:idx2) = currentpool;
				unstimidx = idx2+1;
			end
		end
		function bootCollatedData(obj, nboot)
			% 
			% 	>> bootCollatedData >> poolCollatedData
			% 
			% 	Pools collated session flicks by category and then runs emery boot difference test on POOLED flick data across sessions
			% 	Turns out to give similar results to within session bootstraps, so was not useful (11/6/19)
			% 
			obj.resetCollateHelper
			obj.setCategoryIndicies;
			obj.poolCollatedData
			if nargin < 2
				nboot = 1000000;
			end
			obj.bootStat.collateHelper.Note{1,1} = ['Used bootCollatedData(obj, nboot=' num2str(nboot) ')'];
			% 
			% 	Boot the Collated Pooled Data
			% 
			[obj.bootStat.collateHelper.activation.bootdist] = obj.bootDifference(obj.bootStat.collateHelper.activation.stimpool, obj.bootStat.collateHelper.activation.unstimpool, nboot);
			[obj.bootStat.collateHelper.inhibition.bootdist] = obj.bootDifference(obj.bootStat.collateHelper.inhibition.stimpool, obj.bootStat.collateHelper.inhibition.unstimpool, nboot);
			[obj.bootStat.collateHelper.noOpsin.bootdist] = obj.bootDifference(obj.bootStat.collateHelper.noOpsin.stimpool, obj.bootStat.collateHelper.noOpsin.unstimpool, nboot);
			% 
			% 	Get CI
			% 
			obj.bootStat.collateHelper.activation.CI95 = [obj.bootStat.collateHelper.activation.bootdist(round(0.025*nboot)), obj.bootStat.collateHelper.activation.bootdist(round(0.975*nboot))];
			obj.bootStat.collateHelper.inhibition.CI95 = [obj.bootStat.collateHelper.inhibition.bootdist(round(0.025*nboot)), obj.bootStat.collateHelper.inhibition.bootdist(round(0.975*nboot))];
			obj.bootStat.collateHelper.noOpsin.CI95 = [obj.bootStat.collateHelper.noOpsin.bootdist(round(0.025*nboot)), obj.bootStat.collateHelper.noOpsin.bootdist(round(0.975*nboot))];
			% 
			% 	Median and Mean
			% 
			obj.bootStat.collateHelper.activation.mean = mean(obj.bootStat.collateHelper.activation.bootdist);
			obj.bootStat.collateHelper.activation.median = median(obj.bootStat.collateHelper.activation.bootdist);
			obj.bootStat.collateHelper.inhibition.mean = mean(obj.bootStat.collateHelper.inhibition.bootdist);
			obj.bootStat.collateHelper.inhibition.median = median(obj.bootStat.collateHelper.inhibition.bootdist);
			obj.bootStat.collateHelper.noOpsin.mean = mean(obj.bootStat.collateHelper.noOpsin.bootdist);
			obj.bootStat.collateHelper.noOpsin.median = median(obj.bootStat.collateHelper.noOpsin.bootdist);
			% 
			% 	Plot Histogram
			% 
			figure
			ax1 = subplot(1,3,1);
			ax2 = subplot(1,3,2);
			ax3 = subplot(1,3,3);
			hold(ax1, 'on');
			hold(ax2, 'on');
			hold(ax3, 'on');
			histogram(ax1, obj.bootStat.collateHelper.activation.bootdist, 'normalization', 'probability', 'displayStyle', 'stairs', 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.activation.stimpool)) ' Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.activation.unstimpool))])
			xx = get(ax1, 'xlim');
			yy = get(ax1, 'ylim');
			plot(ax1, [obj.bootStat.collateHelper.activation.CI95(1), obj.bootStat.collateHelper.activation.CI95(1)], yy, 'r-', 'DisplayName', '2.5pct');
			plot(ax1, [obj.bootStat.collateHelper.activation.CI95(2), obj.bootStat.collateHelper.activation.CI95(2)], yy, 'r-', 'DisplayName', '97.5pct');
			plot(ax1, [obj.bootStat.collateHelper.activation.median, obj.bootStat.collateHelper.activation.median], yy, 'k-', 'DisplayName', 'median');
			plot(ax1, [obj.bootStat.collateHelper.activation.mean, obj.bootStat.collateHelper.activation.mean], yy, 'b-', 'DisplayName', 'mean');
			title(ax1, ['Activation'])
			legend(ax1, 'show')

			histogram(ax2, obj.bootStat.collateHelper.inhibition.bootdist, 'normalization', 'probability', 'displayStyle', 'stairs', 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.inhibition.stimpool)) ' Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.inhibition.unstimpool))])
			plot(ax2, [obj.bootStat.collateHelper.inhibition.CI95(1), obj.bootStat.collateHelper.inhibition.CI95(1)], yy, 'r-', 'DisplayName', '2.5pct');
			plot(ax2, [obj.bootStat.collateHelper.inhibition.CI95(2), obj.bootStat.collateHelper.inhibition.CI95(2)], yy, 'r-', 'DisplayName', '97.5pct');
			plot(ax2, [obj.bootStat.collateHelper.inhibition.median, obj.bootStat.collateHelper.inhibition.median], yy, 'k-', 'DisplayName', 'median');
			plot(ax2, [obj.bootStat.collateHelper.inhibition.mean, obj.bootStat.collateHelper.inhibition.mean], yy, 'b-', 'DisplayName', 'mean');
			title(ax2, ['inhibition'])
			legend(ax2, 'show')

			histogram(ax3, obj.bootStat.collateHelper.noOpsin.bootdist, 'normalization', 'probability', 'displayStyle', 'stairs', 'displayName', ['Stim Pool n=' num2str(numel(obj.bootStat.collateHelper.noOpsin.stimpool)) ' Unstim Pool n=' num2str(numel(obj.bootStat.collateHelper.noOpsin.unstimpool))])
			plot(ax3, [obj.bootStat.collateHelper.noOpsin.CI95(1), obj.bootStat.collateHelper.noOpsin.CI95(1)], yy, 'r-', 'DisplayName', '2.5pct');
			plot(ax3, [obj.bootStat.collateHelper.noOpsin.CI95(2), obj.bootStat.collateHelper.noOpsin.CI95(2)], yy, 'r-', 'DisplayName', '97.5pct');
			plot(ax3, [obj.bootStat.collateHelper.noOpsin.median, obj.bootStat.collateHelper.noOpsin.median], yy, 'k-', 'DisplayName', 'median');
			plot(ax3, [obj.bootStat.collateHelper.noOpsin.mean, obj.bootStat.collateHelper.noOpsin.mean], yy, 'b-', 'DisplayName', 'mean');
			title(ax3, ['No Opsin'])
			legend(ax3, 'show')
		end
		function [bootdist] = bootDifference(obj, pool1, pool2, nboot)
			% 
			% 	>> bootCollatedData >> bootDifference
			% 
			% 	Draws a sample from pool 1 and a sample from pool 2
			% 	Then, Bi = pool1i - pool2i
			% 
			% 	Emery boot algorithm
			% 
			bootdist = nan(nboot, 1);
			if nargin < 4
				nboot = 1000000;
			end

			for iboot = 1:nboot
				pool1Idx = randi(numel(pool1));
				pool2Idx = randi(numel(pool2));
				bootdist(iboot) = pool1(pool1Idx) - pool2(pool2Idx);
			end
			bootdist = sort(bootdist);
		end




		% 
		% 	K-S ecdf plotting
		% 
		function f = plotSessionKSResult(obj, pSignificant)
			if nargin < 2
				pSignificant = 0.05;
			end

			if pSignificant == 0.05
				h = [obj.bootStat.collateResult.h];
			else
				h = [obj.bootStat.collateResult.p] < pSignificant;
			end

			[f, ax] = makeStandardFigure(1,[1,1]);

			obj.bootStat.ks2ByGroup.activation = [];
			obj.bootStat.ks2ByGroup.noOpsin = [];
			obj.bootStat.ks2ByGroup.inhibition = [];

			plot(ax, [0,4], [0,0], 'k-', 'LineWidth', 5)
			for iSesh = 1:numel(obj.bootStat.collateResult)
				if h(iSesh) == 1
					markers = '.';
					markersSize = 30;
				else
					markers = 'o';
					markersSize = 10;
				end
				if strcmpi(obj.bootStat.collateResult(iSesh).type, 'ChR2')
					x = 1 + (2 * rand-1)/4;
					markers = ['b',markers];
					obj.bootStat.ks2ByGroup.activation(end+1) = obj.bootStat.collateResult(iSesh).ks2stat;
				elseif strcmpi(obj.bootStat.collateResult(iSesh).type, 'tapers')
					x = 1 + (2 * rand-1)/4;
					markers = ['c',markers];
					obj.bootStat.ks2ByGroup.activation(end+1) = obj.bootStat.collateResult(iSesh).ks2stat;
				elseif strcmpi(obj.bootStat.collateResult(iSesh).type, 'ChrimsonR')
					x = 1 + (2 * rand-1)/4;
					markers = ['r',markers];
					obj.bootStat.ks2ByGroup.activation(end+1) = obj.bootStat.collateResult(iSesh).ks2stat;
				elseif strcmpi(obj.bootStat.collateResult(iSesh).type, 'gtacr2')
					x = 3 + (2 * rand-1)/4;
					markers = ['y',markers];
					obj.bootStat.ks2ByGroup.inhibition(end+1) = obj.bootStat.collateResult(iSesh).ks2stat;
				elseif strcmpi(obj.bootStat.collateResult(iSesh).type, 'no opsin')
					x = 2 + (2 * rand-1)/4;
					markers = ['k',markers];
					obj.bootStat.ks2ByGroup.noOpsin(end+1) = obj.bootStat.collateResult(iSesh).ks2stat;
				end
				plot(ax, x, obj.bootStat.collateResult(iSesh).ks2stat, markers, 'markersize', markersSize);
			end
			boxpnts = [obj.bootStat.ks2ByGroup.activation,obj.bootStat.ks2ByGroup.noOpsin,obj.bootStat.ks2ByGroup.inhibition];
			g = [ones(numel(obj.bootStat.ks2ByGroup.activation),1);2*ones(numel(obj.bootStat.ks2ByGroup.noOpsin),1);3*ones(numel(obj.bootStat.ks2ByGroup.inhibition),1)];
			boxplot(boxpnts, g);
			title(['KS2 results, p*=' num2str(pSignificant)])
			xlim([0,4])
		end



		% 
		% 	AUC permutations
		% 
		function f = plotSessionbootAUCResult(obj, pSignificant)
			if nargin < 2
				pSignificant = 0.025;
			end

			
			h = [obj.bootStat.collateResult.p] < pSignificant | [obj.bootStat.collateResult.p] > 1-pSignificant;
			

			[f, ax] = makeStandardFigure(1,[1,1]);

			obj.bootStat.AUCByGroup.activation = [];
			obj.bootStat.AUCByGroup.noOpsin = [];
			obj.bootStat.AUCByGroup.inhibition = [];

			plot(ax, [0,4], [0,0], 'k-', 'LineWidth', 5)
			for iSesh = 1:numel(obj.bootStat.collateResult)
				if h(iSesh) == 1
					markers = '.';
					markersSize = 60;
				else
					markers = 'o';
					markersSize = 20;
				end
				if strcmpi(obj.bootStat.collateResult(iSesh).type, 'ChR2')
					x = 1 + (2 * rand-1)/4;
					markers = ['b',markers];
					obj.bootStat.AUCByGroup.activation(end+1) = obj.bootStat.collateResult(iSesh).delAUC_EOT;
				elseif strcmpi(obj.bootStat.collateResult(iSesh).type, 'tapers')
					x = 1 + (2 * rand-1)/4;
					markers = ['c',markers];
					obj.bootStat.AUCByGroup.activation(end+1) = obj.bootStat.collateResult(iSesh).delAUC_EOT;
				elseif strcmpi(obj.bootStat.collateResult(iSesh).type, 'ChrimsonR')
					x = 1 + (2 * rand-1)/4;
					markers = ['r',markers];
					obj.bootStat.AUCByGroup.activation(end+1) = obj.bootStat.collateResult(iSesh).delAUC_EOT;
				elseif strcmpi(obj.bootStat.collateResult(iSesh).type, 'gtacr2')
					x = 3 + (2 * rand-1)/4;
					markers = ['y',markers];
					obj.bootStat.AUCByGroup.inhibition(end+1) = obj.bootStat.collateResult(iSesh).delAUC_EOT;
				elseif strcmpi(obj.bootStat.collateResult(iSesh).type, 'no opsin')
					x = 2 + (2 * rand-1)/4;
					markers = ['k',markers];
					obj.bootStat.AUCByGroup.noOpsin(end+1) = obj.bootStat.collateResult(iSesh).delAUC_EOT;
				end
				plot(ax, x, obj.bootStat.collateResult(iSesh).delAUC_EOT, markers, 'markersize', markersSize);
			end
			boxpnts = [obj.bootStat.AUCByGroup.activation,obj.bootStat.AUCByGroup.noOpsin,obj.bootStat.AUCByGroup.inhibition];
			g = [ones(numel(obj.bootStat.AUCByGroup.activation),1);2*ones(numel(obj.bootStat.AUCByGroup.noOpsin),1);3*ones(numel(obj.bootStat.AUCByGroup.inhibition),1)];
			boxplot(boxpnts, g);
			title(['delAUC_EOT results, p*=' num2str(pSignificant)])
			xlim([0,4])
		end








		% 
		% 	end of methods
		% 
		function pathstr = correctPathOS(obj,pathstr)
			if ispc
    			pathstr = strjoin(strsplit(pathstr, '/'), '\');
			else
				pathstr = [strjoin(strsplit(pathstr, '\'), '/')];
			end
		end
		function printFigure(obj,name, f)
			if nargin < 3
				f = gcf;
			end
			print(f,'-depsc','-painters', [name, '.eps'])
		end






		

	end
end



