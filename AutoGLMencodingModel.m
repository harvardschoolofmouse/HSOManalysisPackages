%  AutoGLMencodingModel.m
jobID = randi(10000);
try
	STDmultiplier = 2*ones(4,1);%4.*ones(8,1);
	redmode = true;
	noControlMode = false;
    noRedButHasMove = false;
    noControlButHasRed = false;
	lamda = 0; %0.2;
	smoothing = 0; % added 12/5/19 -- shouldn't be smoothing before modeling. (was 100 before)
	% 
	% 	Created 	3/12/19	Allison Hamilos		ahamilos{at}g.harvard.edu
	% 	Modified 	12/5/19 Allison Hamilos		ahamilos{at}g.harvard.edu
	% 
	% 	Uses v3x Host Structure
	% 		HOSTObj -> SNcHost -> B5_SNc_13 -- extract the Obj automatically 
	% 
	% 	VERSION IDs:
	% 	VERSION 3.x
% 	1	- 5/19/20:	Version used for figures in Hamilos et al., 2020
	% 		v3x8 -- for v3x8 version of CLASS_photometry_roadmap1v_4 (3/12/19) - prototype version
	% 				Uses nestedGLM v2.6
	% 	Update Log
	% 		12/5/19:	Added no movement control cases to here and photometry statObj code and canned styles
	% 					SHOULDN'T BE SMOOTHING, since model will do better with noise without smooth. so set smoothing to 0 (was 100 before)
	% 		3/16/19:	Added x_style = {'MOVEdelta', STDmultiplier} to nestedGLM -- 
	% 					now we can input the STD multiplier for each dataset.
	% 						For now, we will just use 3, but we could use any to check robustness, 
	% 						and we will apply for both X data and EMG data
	% 						i.e., I'm gonna let model do its default thing
	% 						will need to make a field obj.GLM.MOVE instead of .EMG to use this method
	% 		3/15/19:	Verified that X and EMG generally can each use 2*STD to 
	% 					produce similar kernel for gfit, can use this method to check kernel
	% 					before fitting each animal.
	% 
	% 
	% -----------------------------------------------------------------------------
	% 
	%	Initialize parameters for model fit
	% 
	% 	v3x8 nesting keys (uses model v2.6, 11/26/19)
	% 
	disp('==================================================')
	disp('==             AutoGLMdraft v3x8                ==')
	disp('==================================================')
	suppressPlot = true;
	autoVersionID = 'v3x8';
	runID = ['lam0 ' datestr(now)];
	if ~redmode && ~noControlMode && ~noRedButHasMove
		nesting_keys = {'tdt',...									% tdt
						'BEST_cue_flick_tdt',...					% events + tdt
						'BEST_cue_flick_MOVEdelta_tdt',...			% events + tdt + EMG
						'BEST_cue_flick_boxes_MOVEdelta_tdt',... 	% events + tdt + EMG + t-dep baseline
						'cue_flick_MOVEdelta_tdt_stretchONES'};
						%FOR PERCENTILES: 'cue_flick_rampdelta_MOVEdelta_tdt'};		% events + tdt + EMG + stretch-time + t-dep baseline     'BEST_cue_flick_boxes_vt_MOVEdelta_tdt',...	% events + tdt + EMG + vt + t-dep baseline
		eventIdxs = {1,...
				1:3,...
				1:4,...
				1:5,...
				1:6};  %				1:6 
    elseif noRedButHasMove
        nesting_keys = {'BEST_cue_flick',...s						% events 
                        'BEST_cue_flick_MOVEdelta',...						% events
						'BEST_cue_flick_boxes_MOVEdelta',... 					% events + t-dep baseline
						'BEST_cue_flick_MOVEdelta_stretchONES'};
						%FOR PERCENTILES:'BEST_cue_flick_rampdelta_MOVEdelta'};				% events + stretch-time + t-dep baseline
		eventIdxs = {1:2,...
				1:3,...
				1:4,...
                1:5}; 
    elseif noControlButHasRed
        nesting_keys = {'tdt',...									% tdt
        				'BEST_cue_flick_tdt',...
        				'BEST_cue_flick_boxes_tdt',...
        				'cue_flick_tdt_stretchONES'};
                        % 'BEST_cue_flick',...s						% events 
                        % 'BEST_cue_flick_boxes',...					% events
						
						%FOR PERCENTILES: 'BEST_cue_flick_rampdelta'};				% events + stretch-time + t-dep baseline
		eventIdxs = {1,...
				1:3,...
				1:4,...
                1:5}; 
    elseif noControlMode
		nesting_keys = {'BEST_cue_flick',...						% events 
						'BEST_cue_flick_boxes',... 					% events + t-dep baseline
						'BEST_cue_flick_stretchONES'};
						%FOR PERCENTILES: 'BEST_cue_flick_rampdelta'};				% events + stretch-time + t-dep baseline
		eventIdxs = {1:2,...
				1:3,...
				1:4}; 
	else
		warning('THIS MODE SHOULD ONLY BE USED FOR RED CHANNEL DATA')
		nesting_keys = {'BEST_cue_flick',...
					'BEST_cue_flick_MOVEdelta',... %,... % + EMG
					'BEST_cue_flick_boxes_MOVEdelta',... % + baseline(t)%'BEST_cue_flick_boxes_vt_MOVEdelta',... % + veridical time
					'BEST_cue_flick_MOVEdelta_stretchONES'};
					%FOR PERCENTILES: 'cue_flick_rampdelta_MOVEdelta'}; % + stretch time
		eventIdxs = {1:2,...
				1:3,...
				1:4,...%1:5,...
				1:5}; 
	end




	 % {'BEST_cue_flick',... % events
	% 				'BEST_cue_flick_MOVEdelta',... %,... % + MOVEdelta
	% 				'BEST_cue_flick_boxes_MOVEdelta',... % + baseline(t)
	% 				'BEST_cue_flick_boxes_vt_MOVEdelta',... % + veridical time
	% 				'cue_flick_rampdelta_MOVEdelta'}; % + stretch time
	% nesting_keys = {'BEST_cue_flick',... % events
	% eventIdxs = {1:2,...
	% 			1:3,...%,...
	% 			1:4,...
	% 			1:5,...
	% 			1:5}; % only one baseline in the full model
	gfitStyle = 'box200000';
	yMode = 'trial2lick';
	% 
	% 
	% ==============================================================================
	% 
	% 	Outer layer: 
	% 		Autoload datasets to fit with model.
	% 			Run model on photometry data first + EMG
	% 				DEBUG: Run and compare to equivalent model for X
	% 			Then run model on control data if present (red channel)
	% 		Then save the model fit results with the obj in its folder
	% 		Then add the parameters for cross-day comparison to some larger tabulation (I think just theta values)
	%
	%	Initialize Result container for overall findings 
	RESULTS = {};
	RESULTS.id.autoVersionID = autoVersionID;
	RESULTS.id.runID = runID;
	RESULTS.id.signalname = '';
	RESULTS.id.signaltype_ = '';
	RESULTS.id.datasetMap = {}; % fill in as we go: use results.id field for each session, which has all identifiers

	RESULTS.p.smoothing = smoothing;
	RESULTS.p.lamda = lamda;
	RESULTS.p.yMode = yMode;

	RESULTS.nesting_keys = nesting_keys;

	RESULTS.nest = {}; % fill in as we go for each animal/day
	% 
	% 	nest # = the number in the nesting scheme
	% 
	% 	RESULTS.nest.FIELD - the mouse name/day/signal (to be defined by eval command)
	% 
	% --------------------------------------------------------------------------
	%
	%	1. Request user to select directory with subfolders. 
	% 
	RESULTS.id.signalname = {};
	fn = {'SNc', 'DLS', 'VTA', 'DLSright', 'DLSleftD', 'SNcred', 'VTAred', 'DLSred', 'EMG', 'X', 'Y', 'Z', 'CamO', 'ChR2', 'SNcnovir', 'VTAnovir', 'SNcgreen', 'VTAgreen', 'DLSgreen'};			
	idxrg_phot = 1:5;
	idxrg_ctrlphit = [6:8, 15, 16, 17, 18, 19];
	idxrg_move = 9:13;
	idxrg_stim = 14;
	[indx,~] = listdlg('PromptString','Select data type(s) to include in model...',...
	                           'ListString',fn);
	if isempty(indx)
		disp('cancelled.')
		return
	elseif numel(indx) == 1
		RESULTS.id.signalname = fn(indx);
		if ismember(indx, idxrg_phot) || ismember(indx, idxrg_ctrlphit)
			RESULTS.id.signalname = 'photometry';
		elseif ismember(indx, idxrg_move)
			if indx == idxrg_move(1)
				RESULTS.id.signaltype_ = 'EMG';
				disp('Switching to rect-only EMG gfit style')
				gfitStyle = {'EMG', []};
			elseif ismember(indx, idxrg_move(2:4))
				RESULTS.id.signaltype_ = 'accelerometer';
				disp('Using Abs-bandpass-X gfit style')
				gfitStyle = {'Abs-X', []};RESULTS.id
			elseif indx == idxrg_move(5)
				RESULTS.id.signaltype_ = 'camera';
				disp('Using Abs-CamOderivative gfit style')
				gfitStyle = {'Abs-CamOderivative', []};
			end
		end				
	elseif numel(indx) > 1 
		% 
		% 	The first signal will be the photometry signal of interest. If there's more than one photom signal of interest, we will include both
		% 	If a control signal is also specified, we ignore it at this point.
		% 
		RESULTS.id.signalname = fn(indx(1));
		RESULTS.id.signaltype_ = 'photometry';
		for idx = 2:numel(indx)
			if ismember(indx(idx), idxrg_phot)
				RESULTS.id.signalname(idx) = fn(indx(idx));
			elseif ismember(indx(idx), idxrg_ctrlphit)
				RESULTS.id.signalname(idx) = fn(indx(idx));
				RESULTS.id.signaltype_(idx) = 'photometry';	
			elseif ismember(indx(idx), idxrg_move)
				error('NOT IMPLEMENTED')
			elseif ismember(indx(idx), idxrg_stim)
				error('NOT IMPLEMENTED')
			end
		end
	end
	hostFolder = uigetdir('','Select host folder');
	cd(hostFolder)
	hostFiles = dir(hostFolder);
	dirFlags = [hostFiles.isdir];
	subFolders = hostFiles(dirFlags);
	folderNames = {subFolders(3:end).name};
	folderPaths = {subFolders(3:end).folder};
	obj.iv.files = folderNames;
	disp([strjoin(['The following datasets will be loaded: ' folderNames]) '\n']);



	% 
	% Inner layer: For each object in set:
	% 	Complete all the nested fits based on parameterization
	% 		Save results in {results} structure
	% 
	%%
	for iset = 1:numel(folderNames)
		try	
		    cd(folderNames{iset})
			% 
			% 	Indicate that we are initializing processing for the current subfolder
			% 	
			disp(['=====>>> Processing GLM Models for statObj in folder ' folderNames{iset} ' (' num2str(iset) '/' num2str(numel(folderNames)) ' ' datestr(now,'HH:MM AM') ') ===================']);
			% 
			% 	Check what info is available to us in the subfolder. If we want a box200 gfit, we need to load the gfit. If exclusions are present we will add them
			% 
			dirFiles = dir;
			% 
			% 	Find and load the statObj
			% 
			sObjpos = find(contains({dirFiles.name},'sObj'));
			if isempty(sObjpos)
				sObjpos = find(contains({dirFiles.name},'snpObj'));
				if isempty(sObjpos)
					sObjpos = find(contains({dirFiles.name},'statObj'));
				end
			end
			if ~isempty(sObjpos) && numel(sObjpos) < 2
				disp('		Detected statObj in folder')
				obj = load([dirFiles(sObjpos).folder, '\' dirFiles(sObjpos).name]);
		        sObjfield = fieldnames(obj);
		        eval(['obj = obj.' sObjfield{1} ';']);
			else
				error('No statObj in folder or too many! Processing data from scratch to create a sObj')
			end
			if suppressPlot
				close; % closes the statObj log
			end
			
			% -----------------------------------------------------------------------
			% 
			% 	Obj QC:
			% 
			% 	1. Sampling rate:
			% 
			if obj.Plot.samples_per_ms ~= 1/(1000*mode(obj.GLM.gtimes(2:end) - obj.GLM.gtimes(1:end-1)))
			    warning('Correcting sampling rate for this old obj')
			    obj.Plot.samples_per_ms = 1/(1000*mode(obj.GLM.gtimes(2:end) - obj.GLM.gtimes(1:end-1)));
			end
			% 
			% 	2. Check gfit:
			% 
			if ~isfield(obj.GLM, 'gfitMode')
				if isfield(obj.iv, 'gfit_box_win_')
					warning('gfit Mode not specified in obj. Setting to box based on iv')
					obj.GLM.gfitMode = ['box' num2str(obj.iv.gfit_box_win_)];
				else
					warning('gfit Mode not specified in obj. Setting to box 200000')
					obj.GLM.gfitMode = 'box200000';
				end
			end
			if isfield(obj.GLM, 'gfitMode') && ~strcmp(gfitStyle, obj.GLM.gfitMode)
				error('CLASS_photometry_roadmapv1_4:Incompatible_gfitMode', 'Incompatible gfit style')
			end
			% 
			% 	3. Align smoothing to call:
			% 
			obj.Plot.smooth_kernel = smoothing;
			% 
			%	4. Check that is GLM obj 
			% 
			if ~isfield(obj.GLM, 'Mode') || ~obj.GLM.Mode
				if strcmp(obj.iv.setStyle, 'v3x-single session')
					obj.GLM.Mode = true;
				end
			end
			% 
			% 	5. Check for appropriate EMG field. v3x objects have 'EMG', not 'EMGfit'
			% 
			resave = false;
			if ~isfield(obj.GLM, 'emgFit') && ~noControlMode
				if isfield(obj.GLM, 'EMG')
					obj.GLM.MOVE = abs(obj.GLM.EMG);
					obj.GLM.MOVEtimes = obj.GLM.EMGtimes;
					obj.iv.ctrl_signalname = 'EMG';
					% obj.GLM.emgFit = obj.GLM.EMG;
					% obj.GLM.emgTimes = obj.GLM.EMGtimes;
				else
					resave = true;
					disp('Need to collect movement data!')
					s7spos = find(~contains({dirFiles(3:end).name},'gfit') & ~contains({dirFiles(3:end).name},'excl') & ~contains({dirFiles(3:end).name},'roadmap') & ~contains({dirFiles(3:end).name},'Obj') & ~contains({dirFiles(3:end).name},'GLM') & ~contains({dirFiles(3:end).name},'obj')) + 2;
					s7s = load([dirFiles(s7spos).folder, '\', dirFiles(s7spos).name]);
					signals = fieldnames(s7s);
					if sum(contains(signals, 'EMG'))
	                    fieldIdx = find(contains(signals, 'EMG'));
						eval(cell2mat(['obj.GLM.EMG = s7s.' signals(fieldIdx) '.values;']));
						eval(cell2mat(['obj.GLM.EMGtimes = s7s.' signals(fieldIdx) '.times;']));
						obj.GLM.MOVE = abs(obj.GLM.EMG);
						obj.GLM.MOVEtimes = obj.GLM.EMGtimes;
						obj.iv.ctrl_signalname = 'EMG';
					elseif isfield(obj.GLM, 'X')
						if mean(obj.GLM.X) > 0.2
							warning('Correcting X for Abs-X gfit method...')
							obj.GLM.X = abs(obj.bandPass(obj.GLM.X));
						end
						obj.GLM.MOVE = obj.GLM.X;
						obj.GLM.MOVEtimes = obj.GLM.Xtimes;
						obj.iv.ctrl_signalname = 'X';
						warning('USING X AS CONTROL SIG! RBF!')
					elseif sum(contains(signals, 'X'))
	                    fieldIdx = find(contains(signals, 'X'));
	                    eval(cell2mat(['obj.GLM.X = s7s.' signals(fieldIdx) '.values;']));
	                    eval(cell2mat(['obj.GLM.Xtimes = s7s.' signals(fieldIdx) '.times;']));
						obj.GLM.X = abs(obj.bandPass(obj.GLM.X));
						obj.GLM.MOVE = obj.GLM.X;
						obj.GLM.MOVEtimes = obj.GLM.Xtimes;
						obj.iv.ctrl_signalname = 'X';
						warning('USING X AS CONTROL SIG! RBF!')
						% 
						% 	Since we went to the trouble to update the object, let's rewrite it
						% 
					else
						error('No movement signals for this file')
					end
				end
			elseif noControlMode
				warning('This object is for dataset with no tdt or movement control!!')
				resave = false;
				obj.GLM.MOVE = [];
				obj.GLM.MOVEtimes = [];
				obj.iv.ctrl_signalname = 'none';
			else
				error('THIS IS AN OLD OBJ, CHECK METHODS!')
				resave = true;
				obj.GLM.MOVE = abs(obj.GLM.emgFit);
				obj.GLM.MOVEtimes = obj.GLM.emgTimes;
				obj.iv.ctrl_signalname = 'EMG';
			end

			if ~isfield(obj.GLM, 'tdt') && ~redmode && ~noControlMode && ~noRedButHasMove
				resave = true;
				alert = ['ACTION NEEDED: AutoGLMdraft Job' num2str(jobID) ' requires UI Input']; 
				mailAlert(alert);
				obj.addtdt;
			end
			% 
			% 	Add the STDMultiplier for this obj
			% 
			if numel(STDmultiplier) < iset && ~noControlMode
				STDmultiplier(iset) = STDmultiplier(end);
			elseif noControlMode
				STDmultiplier(iset) = nan;
			end
			obj.GLM.STDmultiplier = STDmultiplier(iset);
			if resave
				disp(['Saving the corrected sObj to folder (' datestr(now,'HH:MM AM') ') '])
				save('sObj_Corrected.mat', 'obj', '-v7.3');
			end

			% ----------------------------------------------------------
			% 
			% 	If we pass QC, initialize the results
			% 
	        results = {};
			results.id.autoVersionID = autoVersionID;
			results.id.runID = runID;
			results.id.mousename_ = obj.iv.mousename_;
			results.id.daynum_ = obj.iv.daynum_;
			results.id.hostObjID = obj.iv.date;
	        results.id.hostObjURL = [dirFiles(sObjpos).folder, '\' dirFiles(sObjpos).name];
			results.id.signalname = obj.iv.signalname{1,1};
			results.id.ctrl_signalname = obj.iv.ctrl_signalname;
			results.id.sessionCode = [results.id.mousename_ '_' results.id.signalname '_' results.id.ctrl_signalname '_' results.id.daynum_ ];
			results.id.STDmultiplier = STDmultiplier(iset);

			results.p.smoothing = smoothing;
			results.p.lamda = lamda;
			results.p.yMode = yMode;
			results.p.STDmultiplier = STDmultiplier(iset);

			results.nest.key = {};
			results.nest.th = {};
			results.nest.X = {};
			results.nest.a = {};
			results.nest.Stat = {}; % set obj.Sta to this when plotting stuff relevant to this fit
			% 
			% 	Compile...
			% 
			RESULTS.id.datasetMap{iset} = results.id; % fill in as we go: use results.id field for each session, which has all identifiers
			% 
			% 	Now run all the models and store data...
			% 
			for ikey = 1:numel(nesting_keys)
				if ikey ~= 1
					obj.GLM.flush.recycleUniformSn = 1; % this will let us recycle the Uniform trial set from earlier nests - keeps models fair across one another...
	                % 
					% 	Delete X, a to save memory
					% 
					results.nest(ikey-1).X = [];
					results.nest(ikey-1).a = [];
				end
				% 
				% 	For each nest, run the model, extract the parameters
				% 
				[results.nest(ikey).th, results.nest(ikey).X, results.nest(ikey).a, yFit, ~, ~, results.nest(ikey).CVmat] = obj.nestedGLM(nesting_keys{ikey}, yMode, lamda); 
				results.nest(ikey).key = nesting_keys{ikey};
				results.nest(ikey).Stat = obj.Stat.GLM;
	            if isfield(obj.GLM.flush,'SnTrialsUniform')
	                results.nest(ikey).Stat.SnTrialsUniform = obj.GLM.flush.SnTrialsUniform;
	            else
	                results.nest(ikey).Stat.SnTrials_sub7 = obj.GLM.flush.SnTrials_sub7;
	            end		
				results.nest(ikey).Stat.binEdges = obj.GLM.flush.binEdges;
				results.nest(ikey).Stat.SnBinIdx = obj.GLM.flush.SnBinIdx;
				results.nest(ikey).Stat.se_model = obj.Stat.GLM.se_model;
				results.nest(ikey).Stat.se_th = obj.Stat.GLM.se_th;
				results.nest(ikey).Stat.GLM.signifCoeff = obj.Stat.GLM.signifCoeff;
				results.nest(ikey).Stat.Resid = obj.Stat.GLM.Resid;
				results.nest(ikey).Stat.std_Resid = obj.Stat.GLM.std_Resid;
				results.nest(ikey).Stat.explainedVarianceR2 = obj.Stat.GLM.explainedVarianceR2;

				if isfield(obj.Stat.GLM, 'tsScalingFactor')
					results.nest(ikey).Stat.obj.Stat.GLM.tsScalingFactor = obj.Stat.GLM.tsScalingFactor;
				end
				[results.nest(ikey).Stat.AIC, results.nest(ikey).Stat.AICc, results.nest(ikey).Stat.nAIC, results.nest(ikey).Stat.BIC] = testAIC(results.nest(ikey).a, results.nest(ikey).th, yFit);
				% 
				% 	Display the simulated CTA
				% 
				obj.simulateCTA(eventIdxs{ikey}, [], suppressPlot)
				title([results.id.mousename_ ' ' results.id.daynum_ ' ' results.id.signalname ' ' results.id.ctrl_signalname])
				% 
				% 	Compile results
				% 
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).th = results.nest(ikey).th;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).se_model = results.nest(ikey).Stat.se_model;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).se_th = results.nest(ikey).Stat.se_th;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).explainedVarianceR2 = results.nest(ikey).Stat.explainedVarianceR2;']); % fill in as we go for each animal/day			
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).meanAloss = results.nest(ikey).Stat.meanAloss;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).modelSquaredLoss = results.nest(ikey).Stat.modelSquaredLoss;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).AIC = results.nest(ikey).Stat.AIC;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).AICc = results.nest(ikey).Stat.AICc;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).nAIC = results.nest(ikey).Stat.nAIC;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).BIC = results.nest(ikey).Stat.BIC;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).eventNames = obj.Stat.GLM.eventNames;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).basisMap = obj.Stat.GLM.basisMap;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).eventMap = obj.Stat.GLM.eventMap;']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).n = numel(results.nest(ikey).a);']); % fill in as we go for each animal/day
				eval(['RESULTS.nest.' results.id.sessionCode '(ikey).CVmat = results.nest(ikey).CVmat;']); % fill in as we go for each animal/day
				% 
				% 	Now save the results to the current folder so we can save memory...:
				% 
				tRunTD = regexprep(runID, {' ', '-', ':'}, '_');
				filename = ['Nest' num2str(ikey) '_GLMresults_' results.id.sessionCode '_' tRunTD];
				save(filename, 'results', '-v7.3');

				alert = ['AutoGLMdraft Job' num2str(jobID) ' in Progress. d=' num2str(iset) '/' num2str(numel(folderNames)) ' n=' num2str(ikey) '/' num2str(numel(nesting_keys)) ' Complete']; 
			    mailAlert(alert);
			end
			% 
			% plot the final results for last nest of model
			% 
			obj.findFeature(results.nest(ikey).th, results.nest(ikey).X);
			% 
			% 	You can load results back into object to replot features if you'd like.
			% 
			cd('..')
		catch ME
			warning(['ERROR! ' ME.message])
			if strcmp(ME.identifier, 'CLASS_photometry_roadmapv1_4:Incompatible_gfitMode')
				disp(['The gfitMode for this obj didn''t match the call. Requested ' gfitStyle ' and received ' obj.GLM.gfitMode '. Passing...'])
			else
			    rethrow(ME)
		    end
		end
	end

	% 
	% 	If completes successfully, save the results struct


	tRunTD = regexprep(runID, {' ', '-', ':'}, '_');
	filename = [autoVersionID '_GLM_RESULTS_' tRunTD];
	save(filename, 'RESULTS', '-v7.3');
	alert = ['AutoGLMdraft Job' num2str(jobID) ' COMPLETE without errors!']; 
	mailAlert(alert);
catch EX
	alert = ['ERROR in AutoGLMdraft Job' num2str(jobID) ' in Progress. n=' num2str(iset) '/' num2str(numel(folderNames))]; 
    msg = ['Exception Thrown: ' EX.identifier ' | ' EX.message '\n\n' jobID];
    mailAlert(alert, msg);
    rethrow(EX);
end


disp('Completed for all files in HOST.')
disp('')
disp('============================================')
 
