% cannedStyleSheet_nestedGLMv3x8.m
% 
% 	Created 	3/27/19		Allison Hamilos 	ahamilos{at}g.harvard.edu
% 	Modified 	3/27/19		Allison Hamilos 	ahamilos{at}g.harvard.edu
% 
% 	VERSION 3.8
% 	1	- 5/19/20:	Version used for figures in Hamilos et al., 2020
% 

% 	11/26/19 BEST MODEL v2.6: To make figures for ppt and paper using nested model:
		%		load('STATOBJ.mat')
		% 		[th, X, a, yFit, x, xValidationStruct] = obj.nestedGLM('BEST_cue_flick', 'trial2lick', .2); 
		% 		obj.simulateCTA([1:2])
		% 		print(3,'-depsc','-painters', 'b5d13_clGLM_c.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clGLM_l.eps')
		% 		print(4,'-depsc','-painters', 'b5d13_clGLM_simCTA.eps')
		% 		[th, X, a, yFit, x, xValidationStruct] = obj.nestedGLM('BEST_cue_flick_EMG', 'trial2lick', .2);
		% 		obj.simulateCTA([1:3])
		%		print(3,'-depsc','-painters', 'b5d13_cleGLM_e.eps')
		%		print(3,'-depsc','-painters', 'b5d13_cleGLM_l.eps')
		%		print(3,'-depsc','-painters', 'b5d13_cleGLM_c.eps')
		%		[th, X, a, yFit, x, xValidationStruct] = obj.nestedGLM('BEST_cue_flick_boxes_EMG', 'trial2lick', .2);
		% 		print(1,'-depsc','-painters', 'b5d13_clbxseGLM_c.eps')
		% 		print(1,'-depsc','-painters', 'b5d13_clbxseGLM_l.eps')
		% 		print(1,'-depsc','-painters', 'b5d13_clbxseGLM_b1.eps')
		% 		print(1,'-depsc','-painters', 'b5d13_clbxseGLM_b2.eps')
		% 		print(1,'-depsc','-painters', 'b5d13_clbxseGLM_e.eps') 
		% 		obj.simulateCTA([1:5])
		% 		print(2,'-depsc','-painters', 'b5d13_clbxseGLM_sim.eps')
		% 		[th, X, a, yFit, x, xValidationStruct] = obj.nestedGLM('BEST_cue_flick_boxes_vt_EMG', 'trial2lick', .2); % verifical time
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsVtimeeGLM_e.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsVtimeeGLM_b2.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsVtimeeGLM_b1.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsVtimeeGLM_l.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsVtimeeGLM_c.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsVtimeeGLM_vt.eps')
		% 		obj.simulateCTA([1:6])
		% 		print(1,'-depsc','-painters', 'b5d13_clbxsVtimeeGLM_sim.eps')
		% 		[th, X, a, yFit, x, xValidationStruct] = obj.nestedGLM('cue_flick_rampdelta_EMG', 'trial2lick', .2);
		% 		obj.simulateCTA([1:6])
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsstretchtimeeGLM_stretch.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsstretchtimeeGLM_c.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsstretchtimeeGLM_l.eps')
		% 		print(3,'-depsc','-painters', 'b5d13_clbxsstretchtimeeGLM_e.eps')
		% 		print(5,'-depsc','-painters', 'b5d13_clbxsstretchtimeeGLM_sim.eps')
		% 		obj.getBinnedTimeseries(obj.GLM.gfit, 'times', 34, 10000);
		% 		obj.plot('CTA2l', [3:14], 0, 150, 'last-to-first', 1)
		% 		xlim([-5,7])
		% 		ylim([-0.04, 0.05])
		% 		print(6,'-depsc','-painters', 'b5d13_times34_CTA2l.eps')
		% 

disp(['	Using cannedStyle: ', cannedStyle])

if strcmp(cannedStyle, 'cue_firstlicktype')
	obj.Stat.GLM.eventNames = {'cue', 'flick_rxn', 'flick_early', 'flick_reward'};%, 'flick_iti'};
	basis{1} = {'cue'};
	basis{2} = {'lick'};	% rxn
	basis{3} = {'lick'};	% early
	basis{4} = {'lick'};	% reward
	% basis{5} = {'lick'};	% iti
	Events = {obj.GLM.cue_s,...
			 obj.GLM.fLick_s.rxn,...
			 obj.GLM.fLick_s.early,...
			 obj.GLM.fLick_s.reward};
			 % obj.GLM.fLick_s.iti};
    x_style = {{'delta', 1},...
    			 {'delta', 1},...
    			 {'delta', 1},...
    			 {'delta', 1}};
    			 % {'delta', 1}};
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'test_ssStretch')
	obj.Stat.GLM.eventNames = {'ssStretch'};
	basis{1} = {'timing'};		% no basis, will edit the feature map after make x-rep
	
	Events = {{obj.GLM.cue_s, obj.GLM.firstLick_s,[1000,2000,7000]}};
			 
    x_style = {{'ssStretch', true},... % we are constraining stretch features
    			};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'cue_flick_timing')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-ramp', 'timing-box'};%, 'flick_iti'};
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing
	basis{4} = {'timing'};	% timing
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s}};
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'ramp', 'none'},...
                {'boxcar', 'none'}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'timing_ramps')
	palette	% calls the style map from palette.m

    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'cue_flick_EMG')
	% obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-ramp', 'timing-box', 'EMG'};
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'EMGdelta'};
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'EMGdelta'};	% timing
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.emgFit}};
			 % {obj.GLM.emgFit, obj.GLM.emgTimes}};
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'EMGdelta', 1}};
                % {'blur', 200, 2}};
                % {'none', 2}};
                % {'ramp', 'none'},...
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end    
elseif strcmp(cannedStyle, 'cue_flick_timing_EMG')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'timing-ramp-conv', 'EMGdelta'};
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing
	basis{4} = {'ramp-conv'};	% timing
	basis{5} = {'EMGdelta'};	% timing
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.emgFit}};
			 							 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},...
                {'ramp', 'none'},...
                {'EMGdelta', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
% 
% 	11/26/19 BEST MODEL STYLES for v2.6.......................................
% 
elseif strcmp(cannedStyle, 'cue_flick_rampdelta_MOVEdelta')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time', 'MOVEdelta'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{5} = {'MOVEdelta'};	% timing
	% basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE}};
			 % {obj.GLM.emgFit}};
			 % 							 {obj.GLM.lampOff_s, obj.GLM.cue_s},...
			 % {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 							 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'MOVEdelta', STDmultiplier}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'cue_flick_rampdelta')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time'};%, 'EMGdelta'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
% 					basis{5} = {'EMGdelta'};	% timing
	% basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s}};%,...
% 							 {obj.GLM.emgFit}};
			 % 							 {obj.GLM.lampOff_s, obj.GLM.cue_s},...
			 % {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 							 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end

elseif strcmp(cannedStyle, 'BEST_cue_flick')           	

	obj.Stat.GLM.eventNames = {'cue',...				%	1
							 'flick'} %,...				%	2

	count = 1;
	basis{1} = {'cue'};		% cue
	count = count + 1;
	basis{count} = {'lick'};	% flick

    Events = {obj.GLM.cue_s,...								% 	1
			 obj.GLM.firstLick_s} %,...						%	2
			 
    x_style = {{'delta', 1},...					%	1
                {'delta', 1}}; %,...				%	2

    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'BEST_cue_flick_EMG') || strcmp(cannedStyle, 'BEST_cue_flick_boxes_EMG') || strcmp(cannedStyle, 'BEST_cue_flick_boxes_vt_EMG') || strcmp(cannedStyle, 'cue_flick_rampdelta_EMG')
	error('BEST_cue_flick_EMG, BEST_cue_flick_boxes_EMG, BEST_cue_flick_boxes_vt_EMG, and cue_flick_rampdelta_EMG obsolete. Use BEST_(x)_(x)_(x)_MOVEdelta version with UI-defined STDmultiplier for threshold');
elseif strcmp(cannedStyle, 'BEST_cue_flick_MOVEdelta')           	

	obj.Stat.GLM.eventNames = {'cue',...				%	1
							 'flick',...				%	2
							 'MOVEdelta'};				%	7
	count = 1;
	basis{1} = {'cue'};		% cue
	count = count + 1;
	basis{count} = {'lick'};	% flick
	count = count + 1;
	basis{count} = {'MOVEdelta'};	% EMG

    Events = {obj.GLM.cue_s,...								% 	1
			 obj.GLM.firstLick_s,...						%	2
			 {obj.GLM.MOVE}};
			 % {obj.GLM.emgFit}};								%	7
			 
    x_style = {{'delta', 1},...					%	1
                {'delta', 1},...				%	2
                {'MOVEdelta', STDmultiplier}};				%	7

    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'BEST_cue_flick_boxes')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- no convolution, a dummy basis
	
	
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s}};
			 
			 
			 							 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'BEST_cue_flick_boxes_tdt')
	obj.Stat.GLM.eventNames = {'tdt', 'cue', 'flick', 'timing-box'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{4} = {'tdt'}; %'ramp-delta' 'ramp-delta-norm'
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.smooth(obj.GLM.tdt, 150)}};				 		 
			 							 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},...
                {'none', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'BEST_cue_flick_boxes_MOVEdelta')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'MOVEdelta'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- no convolution, a dummy basis
	% basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{4} = {'MOVEdelta'};	% timing
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE}};
			 % {obj.GLM.emgFit}};
			 % {obj.GLM.lampOff_s, obj.GLM.cue_s},...
			 % {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 
			 							 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'MOVEdelta', STDmultiplier}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end

elseif strcmp(cannedStyle, 'BEST_cue_flick_boxes_vt_MOVEdelta')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'ramp-delta-norm','MOVEdelta'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	% basis{5} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{5} = {'MOVEdelta'};	% timing
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE}};
			 % {obj.GLM.emgFit}};
			 % 							 {obj.GLM.lampOff_s, obj.GLM.cue_s},...
			 % {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 							 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},...%{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'MOVEdelta', STDmultiplier}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
% 
% 	Basic sets (pre 2.6)
% 
elseif strcmp(cannedStyle, 'EMG')
	obj.Stat.GLM.eventNames = {'EMGdelta'};%, 'flick_iti'};
	basis{1} = {'EMGdelta'};	% cue
	
	Events = {{obj.GLM.emgFit}};
			 
    x_style = {{'EMGdelta', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'EMGdelta')
	obj.Stat.GLM.eventNames = {'EMGdelta'};
	basis{1} = {'EMGdelta'};	% cue
	
	Events = {{obj.GLM.emgFit}};
			 
    x_style = {{'EMGdelta', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'cue_flick_ramps_EMG')           	

	obj.Stat.GLM.eventNames = {'cue',...				%	1
							 'flick',...				%	2
							 'timing-ramp',...			%	3
							 'timing-box',...			%	4
							 'timing-expx',...			%	5
							 'timing-rad0.5x',...		%	6
							 'EMGdelta'};				%	7
	count = 1;
	basis{1} = {'cue'};		% cue
	count = count + 1;
	basis{count} = {'lick'};	% flick
	count = count + 1;
	basis{count} = {'timing'};	% timing ramp
	count = count + 1;
	basis{count} = {'timing'};	% timing box
	count = count + 1;
	basis{count} = {'timing'};	% timing expx
	count = count + 1;
	basis{count} = {'timing'};	% timing exp0.5x
	count = count + 1;
	basis{count} = {'EMGdelta'};	% EMG

    Events = {obj.GLM.cue_s,...								% 	1
			 obj.GLM.firstLick_s,...						%	2
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...		%	3
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...	%	4
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...		%	5
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...		%	6
			 {obj.GLM.emgFit}};								%	7
			 
    x_style = {{'delta', 1},...					%	1
                {'delta', 1},...				%	2
                {'ramp', 'none'},...			%	3
                {'boxcar', 'none'},...			%	4
                {'exp-set', 1},...				%	5
                {'rad-set', 0.5},...			%	6
                {'EMGdelta', 1}};				%	7

    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'timing')
	obj.Stat.GLM.eventNames = {'timing-box', 'stretch-time', 'EMGdelta'};%'timing-ramp-conv'
	basis{1} = {'timing'};	% timing
	basis{2} = {'timing'};	% 'ramp-conv'
	basis{3} = {'EMGdelta'};	% timing
	
	Events = {{obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.emgFit}};
			 							 
			 
    x_style = {{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'EMGdelta', 1}};
	
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'ramp-convolution')
	% error('IN PROGRESS!!!!!!!!!!')
	obj.Stat.GLM.eventNames = {'timing-ramp-conv'};
	basis{1} = {'ramp-conv'};	% ramp
	
	Events = {{obj.GLM.cue_s, obj.GLM.firstLick_s}};
			 
    x_style = {{'ramp', 'none'}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
% 
% 	Autoregression
% 
elseif strcmp(cannedStyle, 'self')
	obj.Stat.GLM.eventNames = {'self'};
	basis{1} = {'EMG'};	% just reproduces the thing	
	
	Events = {{obj.smooth(obj.GLM.gfit)}};
			 
    x_style = {{'autoregression', 0}}; % makes sure sampling rate is correct and aligns stuff as needed
    %  500 is the shift back in time in ms
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
	
elseif strcmp(cannedStyle, 'autoregression')
	obj.Stat.GLM.eventNames = {'autoregression'};
	basis{1} = {'EMG'};	% just reproduces the thing	
	
	Events = {{obj.smooth(obj.GLM.gfit)}};
			 
    x_style = {{'autoregression', 500}}; % makes sure sampling rate is correct and aligns stuff as needed
    %  500 is the shift back in time in ms
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end

% 
% 	tdt and autoregression methods added 4/4/19--------------------------------------------------------------------------------------------
% 

elseif strcmp(cannedStyle, 'tdt')
	obj.Stat.GLM.eventNames = {'tdt'};
	basis{1} = {'EMG'};	% just reproduces the thing	
	
	Events = {{obj.smooth(obj.GLM.tdt, 150)}};
			 
    x_style = {{'none', 1}}; % makes sure sampling rate is correct and aligns stuff as needed
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end


elseif strcmp(cannedStyle, 'BEST_cue_flick_tdt')           	

	obj.Stat.GLM.eventNames = {'cue',...				%	1
							 'flick',...
							 'tdt'} %,...				%	2

	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'EMG'};		% blank basis, just reproduces the thing

    Events = {obj.GLM.cue_s,...				
			 obj.GLM.firstLick_s,...
			 {obj.smooth(obj.GLM.tdt, 150)}};
			 
    x_style = {{'delta', 1},...					
                {'delta', 1},...
                {'none', 1}}; %,...				

    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end



elseif strcmp(cannedStyle, 'BEST_cue_flick_MOVEdelta_tdt')           	

	obj.Stat.GLM.eventNames = {'cue',...				%	1
							 'flick',...				%	2
							 'MOVEdelta',...
							 'tdt'};				
	basis{1} = {'cue'};			% cue
	basis{2} = {'lick'};		% flick
	basis{3} = {'MOVEdelta'};	% EMG delta basis set
	basis{4} = {'EMG'};			% blank basis, just reproduces the thing, for tdt

    Events = {obj.GLM.cue_s,...								% 	1
			 obj.GLM.firstLick_s,...						%	2
			 {obj.GLM.MOVE},...
			 {obj.smooth(obj.GLM.tdt, 150)}};
			 
			 
    x_style = {{'delta', 1},...					%	1
                {'delta', 1},...				%	2
                {'MOVEdelta', STDmultiplier},...
                {'none', 1}}; %,... for tdt	

    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end


elseif strcmp(cannedStyle, 'BEST_cue_flick_boxes_MOVEdelta_tdt')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'MOVEdelta', 'tdt'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{4} = {'MOVEdelta'};	% movement
	basis{5} = {'EMG'};		% tdt
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE},...
			 {obj.smooth(obj.GLM.tdt, 150)}};			 							 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'MOVEdelta', STDmultiplier},...
                {'none', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end

elseif strcmp(cannedStyle, 'BEST_cue_flick_boxes_vt_MOVEdelta_tdt')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'ramp-delta-norm','MOVEdelta', 'tdt'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{5} = {'MOVEdelta'};	% movement
	basis{6} = {'EMG'};		% tdt
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE},...
			 {obj.smooth(obj.GLM.tdt, 150)}};
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},...%{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'MOVEdelta', STDmultiplier},...
                {'none', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end

elseif strcmp(cannedStyle, 'cue_flick_rampdelta_MOVEdelta_tdt')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time', 'MOVEdelta', 'tdt'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- BOX - baseline offset
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{5} = {'MOVEdelta'};	% timing
	basis{6} = {'EMG'};		% tdt
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE},...
			 {obj.smooth(obj.GLM.tdt, 150)}};				 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'MOVEdelta', STDmultiplier},...
                {'none', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'cue_flick_MOVEdelta_tdt_stretchONES')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time-ONES', 'MOVEdelta', 'tdt'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- BOX - baseline offset
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{5} = {'MOVEdelta'};	% timing
	basis{6} = {'EMG'};		% tdt
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE},...
			 {obj.smooth(obj.GLM.tdt, 150)}};				 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'MOVEdelta', STDmultiplier},...
                {'none', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'BEST_cue_flick_rampdelta_MOVEdelta')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time', 'MOVEdelta'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- BOX - baseline offset
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{5} = {'MOVEdelta'};	% timing
	
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE}};				 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'MOVEdelta', STDmultiplier}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'BEST_cue_flick_MOVEdelta_stretchONES')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time-ONES', 'MOVEdelta'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- BOX - baseline offset
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{5} = {'MOVEdelta'};	% timing
	
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.GLM.MOVE}};				 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'MOVEdelta', STDmultiplier}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end

elseif strcmp(cannedStyle, 'BEST_cue_flick_rampdelta')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- BOX - baseline offset
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s}};				 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'BEST_cue_flick_stretchONES')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time-ONES'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- BOX - baseline offset
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s}};				 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end
elseif strcmp(cannedStyle, 'cue_flick_tdt_stretchONES')
	obj.Stat.GLM.eventNames = {'cue', 'flick', 'timing-box', 'stretch-time-ONES', 'tdt'}; %'ramp-delta' 'ramp-delta-norm'
	basis{1} = {'cue'};		% cue
	basis{2} = {'lick'};	% flick
	basis{3} = {'timing'};	% timing -- BOX - baseline offset
	basis{4} = {'timing'};	% timing -- no convolution, a dummy basis
	basis{5} = {'EMG'};		% tdt
	
	Events = {obj.GLM.cue_s,...
			 obj.GLM.firstLick_s,...
			 {obj.GLM.lampOff_s, obj.GLM.firstLick_s},...
			 {obj.GLM.cue_s, obj.GLM.firstLick_s},...
			 {obj.smooth(obj.GLM.tdt, 150)}};				 
			 
    x_style = {{'delta', 1},...
                {'delta', 1},...
                {'boxcar', 'none'},... %{'boxcar', 'none'},...
                {'ramp-untrimmed', 'none'},...
                {'none', 1}};
    			 
    for istyle = 1:numel(x_style)
		disp(['Event Name: ', obj.Stat.GLM.eventNames{istyle}, ' #' num2str(istyle), x_style{istyle}, 'basis: ' basis{istyle}])
    end





end