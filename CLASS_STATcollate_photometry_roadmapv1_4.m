classdef CLASS_STATcollate_photometry_roadmapv1_4 < handle
    %
    %   Github HSOM Version from Hamilos et al., 2020
	% 
	% 	Made for CLASS_photometry_roadmapv1_4.m to combine stat analyses across sessions into a plottable thing.
	% 
	% 	Created 	11/6/19	ahamilos
	% 	Modified 	11/6/19	VERSION CODE: ['CLASS_STATcollate_photometry_roadmapv1_4 v1.0 Modified 11-6-19 13:05 | obj created: ' datestr(now)];
	% 
	% 	VERSION CODE HISTORY
	% 
	properties
		iv
		collatedResults
		analysis
	end

	%-------------------------------------------------------
	%		Methods: Initialization
	%-------------------------------------------------------
	methods
		function obj = CLASS_STATcollate_photometry_roadmapv1_4(collateKey)
			% 
			% 	collateKey:			specifies the kind of analysis to capture in each object
			% 						ht 		-- computes the threshold crossing regression for all files in HOST folder
			% 						-cdf 	-- computes the cdf for all files in HOST
			% 						vt		-- computes the threshold crossing regression for all files in HOST folder
			% 						rim		-- runs RimSummary on a Mary-type statObj -- caution using on AH machines! Only works on Mary Machines
			% 									NB: renaming CLASS_MARY_STATcollate_photometry_roadmapv1_4 for MARY use until ready to combine the files.
			% 									Get version updates from here. Last update on this 12/19/19 17:05
			% 						-PCAdecoding 	-- runs the pca methods AND the decoding models
			% 						-htStiff -- ht with stiff threshold and decoding model on ht (gfit only, no tdt channel)
			% 						-ht_1_Stiff -- ht with stiff threshold and decoding model on ht (one threshold and tdt channel)
			% 						-DecodingPaperFinal -- has all the final versions of the decoding model -- 10s baseline, stiff threshs on tdt, PC methods... but too hard to run
			% 						-	1ht_stiff_ea
			% 						-	multiht_stiff_ea
			% 						-	PC1_3_1httdtstiff
			% 						=	1htPCAstiff_1tdtstiff
			% 						=	multihtPCAstiff_multitdtstiff
			% 						-	multihtPCAstiff_multiPCAtdtstiff
			% 						-PCAcollateSummary --  runs obj.plotPCA, then stores the PCs and the summary
			% 
			% 						baselineANOVAidx -- gets items for selectivity index plots
			% 						baselineANOVAwithLick -- does the 3 way ANOVA datasets
			% 						divergenceIndex -- uses new dF/F baseline bootstrap method
			% 						-singleTrialFigures 	-- bins gfit signal by single trials and then saves the figure with 200ms smoothing to check for outliers
			% 						-singleTrialOutliers -- redoes exclusions on all the datasets, then it generates a composite binned stat obj of your choosing.
			% 						-movectrlcustom -- does the special binning for the LTA2l figure and plots and saves all the figures
			% 
			%  
			% 
			obj.iv.runID = randi(10000);
			obj.iv.versionCode = ['CLASS_STATcollate_photometry_roadmapv1_4 v1.0 Modified 11-6-19 13:05 | obj created: ' datestr(now)];

			if nargin < 1 || isempty(collateKey)
				obj.iv.collateKey = 'cdf';
			else
				obj.iv.collateKey = collateKey;
			end

			obj.getDataset;
	    	% 
	    	% 	Save the collated results
	    	% 
	    	obj.save;
	    	disp('Collation of analyses is complete! The shell object is saved to Host Folder.')
			alert = ['Photometry Collation Analysis Obj Complete: ' num2str(obj.iv.runID)]; 
			reportErrors(obj);
		    mailAlert(alert);
		    if strcmpi(obj.iv.collateKey,'PCAdecoding')
		    	cd('..');
	    	elseif strcmpi(obj.iv.collateKey,'singleTrialOutliers')
	    		alert = ['ACTION NEEDED--@singleTrialOutliers Need to select folders now. Photometry Collation Analysis Obj In Prog...: ' num2str(obj.iv.runID)]; 
	    		mailAlert(alert);
	    		disp('	** We are in the process of reExcluding based on outliers, now creating the binned statObj... Select the host folder!')
	    		obj = CLASS_photometry_roadmapv1_4('v3x', 'times', 34, {'box', 200000}, 30000, [], [], 'off');
    		end
		end
		function getDataset(obj, correctionsMode, addDataMode)
			% 
			% 	To correct an existing dataset:
			% 		correctionsMode = true
			% 		Use the original host folder so it can read the original folder names that are missing data
			% 	
			% 	To add data to an existing dataset:
			% 		addDataMode = true (will overide correctionsMode)
			% 		Use a new host with folders for only the data you want to add.
			% 		Remember to move over the new saved file
			% 
			if nargin <3
				addDataMode = false;
			end
			if nargin <2
				correctionsMode = false;
			end
			if ~correctionsMode && ~addDataMode
				if strcmpi(obj.iv.collateKey, 'rim')
					instructions = sprintf('Collate RIM Analyses across Session Object Initialization -- For use with MARY stat objs only as of 12/15/19 \n\n 1. Set Up HOSTObj folder > rimSummary > then folders within this for each day named MOUSENAME_DAY_STYLE (e.g., M1R_1_hyb500 or M1R_7_2_op0). Process Photometry StatObjs as usual and have ready in folders for collation. IS OKAY TO PUT FOLDERS FROM DIFFERENT MICE \n 2. Select the Host Folder For Collation \n 3. Each folder will be processed and CollatedAnalysisObj saved to the HOST folder for the Collation.')
				else
		    		instructions = sprintf('Collate Analyses across Session Object Initialization \n\n 1. Set Up HOSTObj folder > Signal > then folders within this for each day. Each folder is named MouseID_Signal_Day#. Process Photometry StatObjs as usual and have ready in folders for collation. IS OKAY TO PUT FOLDERS FROM DIFFERENT MICE AND SIGNALS DEPENDING ON ANALYSIS \n 2. Select the Host Folder For Collation \n 3. Each folder will be processed and CollatedAnalysisObj saved to the HOST folder for the Collation.')
	    		end
    		elseif correctionsMode
    			instructions = sprintf('FIXING ERRORS \n\n 1. After fixing problems with files, Go to the original HOSTObj folder > Signal  \n 2. Select the original Host Folder For Collation \n 3. Each folder will be processed and appended to CollatedAnalysisObj saved to the HOST folder for the Collation.')
			elseif addDataMode
				instructions = sprintf('ADDING DATA \n\n 1. Put HOSTObj day folders in a new addHOST folder. Only put new files to append here.  \n 2. Select the new addHost Folder For Addition \n 3. Each folder will be processed and appended to CollatedAnalysisObj saved to the addHOST folder for the Collation. \n\n -- ** REMEMBER TO MOVE ADDED AND APPENDED FILES TO MAIN HOST FOR STORAGE!')
    		end
    		hhh = msgbox(instructions);
	    	
	    	disp('Select the HOST Folder for Collation')
	    	disp('	NB: for PCA, outer HOST contains sessionHOST, PCAfigureHOST, and decodeFigureHOST**')
            hostFolder = uigetdir('', 'Select the HOST Folder for Collation');
            obj.iv.hostFolder = hostFolder;
            if ~strcmpi(obj.iv.collateKey,'baselineANOVAidx') || ~strcmpi(obj.iv.collateKey,'ht') || ~strcmpi(obj.iv.collateKey,'vt') || ~strcmpi(obj.iv.collateKey,'rim')
                cd(hostFolder)
            	filesax = dir;
            	dirFlagsax = [filesax.isdir];
            	Foldersax = filesax(dirFlagsax);
            	name = {Foldersax(3:end).name};
				folder = {Foldersax(3:end).folder};
            	
            	if strcmpi(obj.iv.collateKey, 'cdf')
            		obj.iv.suppressNsave.cdf = hostFolder;
        		elseif strcmpi(obj.iv.collateKey,'divergenceIndex'), strcmpi(obj.iv.collateKey, 'singleTrialFigures') || strcmpi(obj.iv.collateKey,'singleTrialOutliers') || strcmpi(obj.iv.collateKey,'movectrlcustom') || strcmpi(obj.iv.collateKey, 'PCAcollateSummary')
        			cd ..
					disp(	'Select the figureHOST Folder for Collation')
					obj.iv.suppressNsave.singleTrialFigures = uigetdir('', 'Select the figureHOST Folder for Collation');
					obj.iv.suppressNsave.figuresHOST = obj.iv.suppressNsave.singleTrialFigures;
        		else
	            	if isempty(find(contains(name,'sessionHOST'),1))
						disp(	'Select the sessionHOST Folder for Collation')
						sessionHostFolder = uigetdir('', 'Select the HOST Folder for Collation');
	        		else
	        			sessionHostFolder = obj.correctPathOS([hostFolder '/' name{find(contains(name,'sessionHOST'))}]);
	        		end
	        		obj.iv.outerHostFolder = hostFolder;
	                outerHostFolder = obj.iv.outerHostFolder;
	        		obj.iv.hostFolder = sessionHostFolder;
	        		hostFolder = sessionHostFolder;
	        		obj.iv.sessionHostFolder = sessionHostFolder;
	        		
	        		if strcmpi(obj.iv.collateKey,'PCAdecoding')
						if isempty(find(contains(name,'PCAfigureHOST')))
							mkdir('PCAfigureHOST');
							PCAfigureHOSTfolder = obj.correctPathOS([outerHostFolder '/PCAfigureHOST']);
						else					
							PCAfigureHOSTfolder = obj.correctPathOS([outerHostFolder '/' name{find(contains(name,'PCAfigureHOST'))}]);
						end
						obj.iv.PCAfigureHOSTfolder = PCAfigureHOSTfolder;
						cd(PCAfigureHOSTfolder)
		            	filesax = dir;
		            	if size(filesax,1) < 8+2
			            	mkdir('PCAsummary')
			            	mkdir('PCAtestfit')
			            	mkdir('PCAwtVsLickTime')
			            	mkdir('PCAXfitAll')
			            	mkdir('PCAXfitSelected')
			            	mkdir('Xfitbinned')
			            	mkdir('PCAmeanSlope')
			            	mkdir('HOSTdecode')
			            	filesax = dir;   
		            	end         
		            	dirFlagsax = [filesax.isdir];
		            	Foldersax = filesax(dirFlagsax);
		            	name = {Foldersax(3:end).name};
						folder = {Foldersax(3:end).folder};
		            	obj.iv.suppressNsave.PCAsummary = [folder{find(contains(name,'PCAsummary'))} '/' name{find(contains(name,'PCAsummary'))}];
		            	obj.iv.suppressNsave.PCAtestfit = [folder{find(contains(name,'PCAtestfit'))} '/' name{find(contains(name,'PCAtestfit'))}];
		            	obj.iv.suppressNsave.PCAwtVsLickTime = [folder{find(contains(name,'PCAwtVsLickTime'))} '/' name{find(contains(name,'PCAwtVsLickTime'))}];
		            	obj.iv.suppressNsave.PCAXfitAll = [folder{find(contains(name,'PCAXfitAll'))} '/' name{find(contains(name,'PCAXfitAll'))}];
		            	obj.iv.suppressNsave.PCAXfitSelected = [folder{find(contains(name,'PCAXfitSelected'))} '/' name{find(contains(name,'PCAXfitSelected'))}];
		            	obj.iv.suppressNsave.Xfitbinned = [folder{find(contains(name,'Xfitbinned'))} '/' name{find(contains(name,'Xfitbinned'))}];
		            	obj.iv.suppressNsave.PCAmeanSlope = [folder{find(contains(name,'PCAmeanSlope'))} '/' name{find(contains(name,'PCAmeanSlope'))}];	

		                cd(outerHostFolder)
	                end

	            	filesax = dir;
	            	dirFlagsax = [filesax.isdir];
	            	Foldersax = filesax(dirFlagsax);
	            	name = {Foldersax(3:end).name};
					folder = {Foldersax(3:end).folder};
					if isempty(find(contains(name,'decodeFigureHOST')))
						mkdir('decodeFigureHOST');
						decodeFigureHOSTfolder = obj.correctPathOS([outerHostFolder '/decodeFigureHOST']);
					else					
						decodeFigureHOSTfolder = obj.correctPathOS([outerHostFolder '/' name{find(contains(name,'decodeFigureHOST'))}]);
					end
					obj.iv.decodeFigureHOSTfolder = decodeFigureHOSTfolder;
%                     obj.iv.suppressNsave.decodeFigureHOSTfolder = decodeFigureHOSTfolder;
					cd(decodeFigureHOSTfolder)
	            	filesax = dir;
	            	if strcmpi(obj.iv.collateKey,'PCAdecoding')
		            	if size(filesax,1) < 5+4
			            	mkdir('ht')
			            	mkdir('htPCA')
			            	mkdir('PC1_3')
			            	mkdir('PC1_3htPCA')
			            	mkdir('mislope')
			            	mkdir('htPCA-stiff')
			            	mkdir('PCA_ht_1_Stiff')
			            	filesax = dir;   
		            	end
	            	elseif strcmpi(obj.iv.collateKey,'DecodingPaperFinal') || strcmpi(obj.iv.collateKey,'1ht_stiff_ea') || strcmpi(obj.iv.collateKey,'multiht_stiff_ea') || strcmpi(obj.iv.collateKey,'PC1_3_1httdtstiff') || strcmpi(obj.iv.collateKey,'1htPCAstiff_1tdtstiff') || strcmpi(obj.iv.collateKey,'multihtPCAstiff_multitdtstiff') || strcmpi(obj.iv.collateKey,'multihtPCAstiff_multiPCAtdtstiff') 
	            		if size(filesax,1) < 6+2
			            	mkdir('1ht_stiff_ea')
			            	mkdir('multiht_stiff_ea')
			            	mkdir('PC1_3_1httdtstiff')
			            	mkdir('1htPCAstiff_1tdtstiff')
			            	mkdir('multihtPCAstiff_multitdtstiff')
			            	mkdir('multihtPCAstiff_multiPCAtdtstiff')
			            	filesax = dir;   
		            	end
	            	elseif strcmpi(obj.iv.collateKey,'htStiff') || strcmpi(obj.iv.collateKey,'ht_1_Stiff')
	            		dirFlagsax = [filesax.isdir];
		            	Foldersax = filesax(dirFlagsax);
		            	name = {Foldersax(3:end).name};
	            		if ~contains(name, 'ht-stiff')
			            	mkdir('ht-stiff')
			            	filesax = dir;   
			            	dirFlagsax = [filesax.isdir];
			            	Foldersax = filesax(dirFlagsax);
			            	name = {Foldersax(3:end).name};
		            	end
		            	folder = {Foldersax(3:end).folder};
	            	end
	            	
					% idx = strcmpi(name, 'ht');
                    
					if strcmpi(obj.iv.collateKey,'PCAdecoding')
                        filesax = dir;
                        dirFlagsax = [filesax.isdir];
                        Foldersax = filesax(dirFlagsax);
                        name = {Foldersax(3:end).name};
                        folder = {Foldersax(3:end).folder};
		            	obj.iv.suppressNsave.ht = [folder{strcmpi(name,'ht')} '/' name{strcmpi(name,'ht')}];
		            	obj.iv.suppressNsave.htPCA = [folder{strcmpi(name,'htPCA')} '/' name{strcmpi(name,'htPCA')}];
		            	obj.iv.suppressNsave.PC1_3 = [folder{strcmpi(name,'PC1_3')} '/' name{strcmpi(name,'PC1_3')}];
		            	obj.iv.suppressNsave.PC1_3htPCA = [folder{strcmpi(name,'PC1_3htPCA')} '/' name{strcmpi(name,'PC1_3htPCA')}];
		            	obj.iv.suppressNsave.mislope = [folder{strcmpi(name,'mislope')} '/' name{strcmpi(name,'mislope')}];
		            	obj.iv.suppressNsave.htPCA_stiff = [folder{strcmpi(name,'htPCA-stiff')} '/' name{strcmpi(name,'htPCA-stiff')}];
		            	obj.iv.suppressNsave.PCA_ht_1_Stiff  = [folder{strcmpi(name,'PCA_ht_1_Stiff')} '/' name{strcmpi(name,'PCA_ht_1_Stiff')}];
	            	elseif strcmpi(obj.iv.collateKey,'DecodingPaperFinal') || strcmpi(obj.iv.collateKey,'1ht_stiff_ea') || strcmpi(obj.iv.collateKey,'multiht_stiff_ea') || strcmpi(obj.iv.collateKey,'PC1_3_1httdtstiff') || strcmpi(obj.iv.collateKey,'1htPCAstiff_1tdtstiff') || strcmpi(obj.iv.collateKey,'multihtPCAstiff_multitdtstiff') || strcmpi(obj.iv.collateKey,'multihtPCAstiff_multiPCAtdtstiff') 
	            		filesax = dir;
                        dirFlagsax = [filesax.isdir];
                        Foldersax = filesax(dirFlagsax);
                        name = {Foldersax(3:end).name};
                        folder = {Foldersax(3:end).folder};
		            	obj.iv.suppressNsave.ht_stiff_ea = [folder{strcmpi(name,'1ht_stiff_ea')} '/' name{strcmpi(name,'1ht_stiff_ea')}];
		            	obj.iv.suppressNsave.multiht_stiff_ea = [folder{strcmpi(name,'multiht_stiff_ea')} '/' name{strcmpi(name,'multiht_stiff_ea')}];
		            	obj.iv.suppressNsave.PC1_3_1httdtstiff = [folder{strcmpi(name,'PC1_3_1httdtstiff')} '/' name{strcmpi(name,'PC1_3_1httdtstiff')}];
		            	obj.iv.suppressNsave.htPCAstiff_1tdtstiff = [folder{strcmpi(name,'1htPCAstiff_1tdtstiff')} '/' name{strcmpi(name,'1htPCAstiff_1tdtstiff')}];
		            	obj.iv.suppressNsave.multihtPCAstiff_multitdtstiff = [folder{strcmpi(name,'multihtPCAstiff_multitdtstiff')} '/' name{strcmpi(name,'multihtPCAstiff_multitdtstiff')}];
		            	obj.iv.suppressNsave.multihtPCAstiff_multiPCAtdtstiff = [folder{strcmpi(name,'multihtPCAstiff_multiPCAtdtstiff')} '/' name{strcmpi(name,'multihtPCAstiff_multiPCAtdtstiff')}];
	            	elseif strcmpi(obj.iv.collateKey,'htStiff') || strcmpi(obj.iv.collateKey,'ht_1_Stiff')
	            		obj.iv.suppressNsave.ht_stiff = [folder{strcmpi(name,'ht-stiff')} '/' name{strcmpi(name,'ht-stiff')}];
	        		end
        		end
        	end
        	cd(hostFolder)

            if exist('hhh', 'var')
                close(hhh);
            end
        	disp('====================================================')
		    disp('			Collated Photometry Analysis Processing 	 	  ')
		    disp('====================================================')
		    disp(' ')
		    disp(['Started: ' datestr(now)])
		    disp(' ')
            hostFiles = dir(hostFolder);
			dirFlags = [hostFiles.isdir];
			subFolders = hostFiles(dirFlags);
			folderNames = {subFolders(3:end).name};
			% folderPaths = {subFolders(3:end).folder};
			if correctionsMode
				% 
				% 	Only load new folders and folders with error
				% 
				EEidxs = obj.reportErrors(false);
                EEidxsliteral = find(EEidxs);
				folderNames = {obj.collatedResults(EEidxs).sessionID};
			elseif addDataMode
				obj.iv.files(end+1:end+numel(folderNames)) = folderNames;
			else 
				obj.iv.files = folderNames;
			end

			disp(char(['Loading the following files...' folderNames]))
			disp(' ')
	    	disp('-----------------------------------------------------------------')
	
            
			for ipos = 1:numel(folderNames)
				if correctionsMode
					iset = EEidxsliteral(ipos);
				elseif addDataMode
					iset = numel({obj.collatedResults.sessionID})+1;
				else 
					iset = ipos;
				end
				result = [];
				fprintf(['Working on file #' num2str(ipos) ': ' folderNames{ipos} '(' num2str(ipos) '/' num2str(numel(folderNames)) ' ' datestr(now,'HH:MM AM') ') \n'])
        		cd(folderNames{ipos})
        		% 
				% 	Check what info is available to us in the subfolder. If we want a box200 gfit, we need to load the gfit. If exclusions are present we will add them
				% 
				dirFiles = dir;
				% 
				% 	First, ensure statObj is already present:
				% 
				sObjpos = find(contains({dirFiles.name},'sObj'));
				if isempty(sObjpos)
					sObjpos = find(contains({dirFiles.name},'snpObj'));
					if isempty(sObjpos)
						sObjpos = find(contains({dirFiles.name},'statObj'));
					end
				end
				
				try
					if ~isempty(sObjpos)
						% 
						% 	Find the newest version of the obj
						% 
						idxdates = [dirFiles(sObjpos).datenum];
						newestObj = find(idxdates == max([dirFiles(sObjpos).datenum]));
						sObjpos = sObjpos(newestObj);

						pathstr = obj.correctPathOS([dirFiles(sObjpos).folder, '\' dirFiles(sObjpos).name]);
						sObj = load(pathstr);
	                    sObjfield = fieldnames(sObj);
	                    eval(['sObj = sObj.' sObjfield{1} ';']);
	                    if strcmpi(obj.iv.collateKey, 'rim')
	                    	istyle = strsplit(folderNames{ipos}, '_');
	                    	sObj.iv.Style = istyle{end};
                    	end

	                    result = obj.analyze(sObj, obj.iv.collateKey);
	                    obj.collatedResults(iset).analysisType = obj.iv.collateKey;
	                    obj.collatedResults(iset).sessionID = folderNames{ipos};
	                    if strcmpi(obj.iv.collateKey, 'cdf')
	                    	obj.collatedResults(iset).ecdf_f = result.ecdf_f;
	                    	obj.collatedResults(iset).ecdf_x = result.ecdf_x;
	                    	obj.collatedResults(iset).f_lick_ex_s_ecdf = result.result.f_lick_ex_s_wrtref_cdf;
	                    	obj.collatedResults(iset).lick_s = result.result.lick_s;
	                    	obj.collatedResults(iset).f_lick_ex_s_wrtref = result.result.f_lick_ex_s_wrtref;
	                    	obj.collatedResults(iset).rb_s = result.result.rb_ms/1000;
	                    	obj.collatedResults(iset).note = result.result.note;
                    	elseif strcmpi(obj.iv.collateKey, 'PCAdecoding')
                    		obj.collatedResults(iset).PCA = result.PCA;
                    		obj.collatedResults(iset).decoding.ht = result.decoding.ht;
                    		obj.collatedResults(iset).decoding.htPCA = result.decoding.htPCA;
                    		obj.collatedResults(iset).decoding.PC1_3 = result.decoding.PC1_3;
                    		obj.collatedResults(iset).decoding.PC1_3htPCA = result.decoding.PC1_3htPCA;
                    		obj.collatedResults(iset).decoding.mislope = result.decoding.mislope;
                    		obj.collatedResults(iset).decoding.htPCA_stiff = result.decoding.htPCA_stiff;
                    		obj.collatedResults(iset).decoding.PCA_ht_1_Stiff = result.decoding.PCA_ht_1_Stiff;
                    		
                    		obj.collatedResults(iset).ht_raw = result.decoding.ht_raw;
							obj.collatedResults(iset).ht_PCA = result.decoding.ht_PCA;
                    		obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
                		elseif strcmpi(obj.iv.collateKey, 'DecodingPaperFinal')
                			obj.collatedResults(iset).decoding.ht_stiff_ea = result.decoding.ht_stiff_ea;
                    		obj.collatedResults(iset).decoding.multiht_stiff_ea = result.decoding.multiht_stiff_ea;
                    		obj.collatedResults(iset).decoding.PC1_3_1httdtstiff = result.decoding.PC1_3_1httdtstiff;
                    		obj.collatedResults(iset).decoding.htPCAstiff_1tdtstiff = result.decoding.htPCAstiff_1tdtstiff;
                    		obj.collatedResults(iset).decoding.multihtPCAstiff_multitdtstiff = result.decoding.multihtPCAstiff_multitdtstiff;
                    		obj.collatedResults(iset).decoding.multihtPCAstiff_multiPCAtdtstiff = result.decoding.multihtPCAstiff_multiPCAtdtstiff;
                    		obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
                    		obj.collatedResults(iset).flagNoRed = result.flagNoRed;
                		elseif strcmpi(obj.iv.collateKey,'1ht_stiff_ea')
                    		obj.collatedResults(iset).decoding.ht_stiff_ea = result.decoding.ht_stiff_ea;
                    		obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
                    		obj.collatedResults(iset).flagNoRed = result.flagNoRed;
                		elseif strcmpi(obj.iv.collateKey,'multiht_stiff_ea')
                			obj.collatedResults(iset).decoding.multiht_stiff_ea = result.decoding.multiht_stiff_ea;
                    		obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
                    		obj.collatedResults(iset).flagNoRed = result.flagNoRed;
                		elseif strcmpi(obj.iv.collateKey,'PC1_3_1httdtstiff')
                			obj.collatedResults(iset).decoding.PC1_3_1httdtstiff = result.decoding.PC1_3_1httdtstiff;
                    		obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
                    		obj.collatedResults(iset).flagNoRed = result.flagNoRed;
                		elseif strcmpi(obj.iv.collateKey,'1htPCAstiff_1tdtstiff')
                			obj.collatedResults(iset).decoding.htPCAstiff_1tdtstiff = result.decoding.htPCAstiff_1tdtstiff;
                    		obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
                    		obj.collatedResults(iset).flagNoRed = result.flagNoRed;
                		elseif strcmpi(obj.iv.collateKey,'multihtPCAstiff_multitdtstiff')
                			obj.collatedResults(iset).decoding.multihtPCAstiff_multitdtstiff = result.decoding.multihtPCAstiff_multitdtstiff;
                    		obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
                    		obj.collatedResults(iset).flagNoRed = result.flagNoRed;
                		elseif strcmpi(obj.iv.collateKey,'multihtPCAstiff_multiPCAtdtstiff') 
                			obj.collatedResults(iset).decoding.multihtPCAstiff_multiPCAtdtstiff = result.decoding.multihtPCAstiff_multiPCAtdtstiff;
                    		obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
                    		obj.collatedResults(iset).flagNoRed = result.flagNoRed;
                		elseif strcmpi(obj.iv.collateKey, 'htStiff') || strcmpi(obj.iv.collateKey, 'ht_1_Stiff')
                			obj.collatedResults(iset).decoding.ht_stiff = result.decoding.ht_stiff;
                			obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;

                    	elseif strcmpi(obj.iv.collateKey, 'ht')
                    		obj.collatedResults(iset).note = result.note;
                    		obj.collatedResults(iset).nthresh = result.nthresh;
                    		obj.collatedResults(iset).nTrials_InRange = result.nTrialsTotal;
                    		obj.collatedResults(iset).nbins_inRange = result.nbins_inRange;

                    		obj.collatedResults(iset).binned_ntrialsperbin = result.binned_ntrialsperbin;
                    		obj.collatedResults(iset).delay_ms = result.delay_ms;
                    		obj.collatedResults(iset).smoothing_samples = result.smoothing_samples;

                    		obj.collatedResults(iset).singletrial_b = result.singletrial_b;
                    		obj.collatedResults(iset).singletrial_dev = result.singletrial_dev;
                    		obj.collatedResults(iset).singletrial_stats = result.singletrial_stats;
                    		obj.collatedResults(iset).singletrial_rsq = result.singletrial_rsq;
                    		obj.collatedResults(iset).singletrial_nbinsXing = result.singletrial_nbinsXing;
                    		obj.collatedResults(iset).singletrial_pcTrialsXing = result.singletrial_pcTrialsXing;
                    		obj.collatedResults(iset).singletrial_time2lickFromThreshXing = result.singletrial_time2lickFromThreshXing;
                    		obj.collatedResults(iset).binned_b = result.binned_b;
                    		obj.collatedResults(iset).binned_dev = result.binned_dev;
                    		obj.collatedResults(iset).binned_stats = result.binned_stats;
                    		obj.collatedResults(iset).binned_rsq = result.binned_rsq;
                    		obj.collatedResults(iset).binned_nbinsXing = result.binned_nbinsXing;
                    		obj.collatedResults(iset).binned_pcBinsXing = result.binned_pcBinsXing;
                    		obj.collatedResults(iset).binned_time2lickFromThreshXing = result.binned_time2lickFromThreshXing;
                		elseif strcmpi(obj.iv.collateKey, 'vt')
                    		obj.collatedResults(iset).note = result.note;
                    		obj.collatedResults(iset).nthresh = result.nthresh;
                    		obj.collatedResults(iset).nTrialsTotal = result.nTrialsTotal;
                    		obj.collatedResults(iset).nBinsTotal = result.nBinsTotal;
                    		obj.collatedResults(iset).thresholds = result.thresholds;

                    		obj.collatedResults(iset).binned_ntrialsperbin = result.binned_ntrialsperbin;
                    		obj.collatedResults(iset).delay_ms = result.delay_ms;
                    		obj.collatedResults(iset).smoothing_samples = result.smoothing_samples;

                    		obj.collatedResults(iset).singletrial_b = result.singletrial_b;
                    		obj.collatedResults(iset).singletrial_dev = result.singletrial_dev;
                    		obj.collatedResults(iset).singletrial_stats = result.singletrial_stats;
                    		obj.collatedResults(iset).singletrial_rsq = result.singletrial_rsq;
                    		obj.collatedResults(iset).singletrial_r = result.singletrial_r;
                    		
                    		
                    		obj.collatedResults(iset).binned_b = result.binned_b;
                    		obj.collatedResults(iset).binned_dev = result.binned_dev;
                    		obj.collatedResults(iset).binned_stats = result.binned_stats;
                    		obj.collatedResults(iset).binned_rsq = result.binned_rsq;
                    		obj.collatedResults(iset).binned_r = result.binned_r;
                		elseif strcmpi(obj.iv.collateKey,'baselineANOVAidx') || strcmpi(obj.iv.collateKey,'baselineANOVAwithLick')
                			obj.collatedResults(iset).results = result.results;
                			obj.collatedResults(iset).F_nm1 = result.F_nm1;
                			obj.collatedResults(iset).F_n = result.F_n;
                			obj.collatedResults(iset).nm1Score = result.nm1Score;
                			obj.collatedResults(iset).nScore = result.nScore;
                			obj.collatedResults(iset).sig_nm1 = result.sig_nm1;
                			obj.collatedResults(iset).sig_n = result.sig_n;
                			obj.collatedResults(iset).centers = result.centers;
                			obj.collatedResults(iset).baselineWindow = result.baselineWindow;
            			elseif strcmpi(obj.iv.collateKey, 'PCAcollateSummary')
            				obj.collatedResults(iset).PCA = result.PCA;
            				obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
        				elseif strcmpi(obj.iv.collateKey, 'divergenceIndex')
        					result.Stat.convergenceIndex.XE = [];
        					result.Stat.convergenceIndex.XR = [];
        					result.Stat.divergenceIndex.EX = [];
							result.Stat.divergenceIndex.RX = [];
            				obj.collatedResults(iset).Stat = result.Stat;
            				obj.collatedResults(iset).flick_s_wrtc = result.flick_s_wrtc;
	                	end
	                    
	                    
					else
						error('No statObj in the folder!')
					end	
					disp(' ')
					disp('	File analyzed and added to Collated Obj.')
				catch ex
					EE = getReport(ex)
					disp(EE)
					warning(['Error while processing this file. The message was:' EE])
					disp('	Skipping this file. Add it to the collated analysis obj later')

                    obj.collatedResults(iset).analysisType = obj.iv.collateKey;
                    obj.collatedResults(iset).sessionID = folderNames{ipos};
                    obj.collatedResults(iset).error = ['Error Encountered:' EE];

				end
	    		disp('-----------------------------------------------------------------')
                cd(hostFolder)
			end
			if correctionsMode || addDataMode
				obj.save;
			end
		end
		function fixErrors(obj)
			% 
			% 	Go to the original host folder after stuff fixed in folders or code and run again
			% 
			obj.getDataset(true, false);
		end
		function addData(obj)
			% 
			% 	Go to a distinct host folder with the new data host days to add and run this
			% 
			obj.getDataset(false,true);
		end


		function [idxs,ok] = reportErrors(obj, verbose)
			if nargin < 2
				verbose = true;
			end
			if verbose
				disp('----------------	Error Report -----------------')
				disp('There were errors processing the following datasets:')
				disp(' ')
			end
			if ~isfield(obj.collatedResults, 'error')
				idxs = zeros(length({obj.collatedResults.sessionID}),1);
				ok = ones(length({obj.collatedResults.sessionID}),1);
				return
			end
			if strcmpi(obj.iv.collateKey, 'cdf')
				if verbose, disp(char({obj.collatedResults(cellfun(@(x) isempty(x), {obj.collatedResults.rb_s})).sessionID}')), end
				idxs = cellfun(@(x) isempty(x), {obj.collatedResults.rb_s});
				ok = cellfun(@(x) ~isempty(x), {obj.collatedResults.rb_s});
			elseif strcmpi(obj.iv.collateKey, 'PCAdecoding') || strcmpi(obj.iv.collateKey, 'htStiff') || strcmpi(obj.iv.collateKey, 'ht_1_Stiff') || strcmpi(obj.iv.collateKey, 'baselineANOVAidx') || strcmpi(obj.iv.collateKey, 'baselineANOVAwithLick')
                if verbose, disp(char({obj.collatedResults(cellfun(@(x) ~isempty(x), {obj.collatedResults.error})).sessionID}')), end
				idxs = cellfun(@(x) ~isempty(x), {obj.collatedResults.error});
				ok = cellfun(@(x) isempty(x), {obj.collatedResults.error});
			elseif strcmpi(obj.iv.collateKey, 'ht')
				if verbose, disp(char({obj.collatedResults(cellfun(@(x) isempty(x), {obj.collatedResults.binned_time2lickFromThreshXing})).sessionID}')), end
				idxs = cellfun(@(x) isempty(x), {obj.collatedResults.binned_time2lickFromThreshXing});
				ok = cellfun(@(x) ~isempty(x), {obj.collatedResults.binned_time2lickFromThreshXing});
			end
		end

		function replaceErrors(obj)

		end

		

		function result = analyze(obj, sObj, collateKey)
			% 
			% 	collateKey: 	A keyword indicating the analysis to complete
			% 				cdf
			% 				ht
			%
			% ----------
			if nargin < 2
				collateKey = 'cdf';
				disp('-------------default analysis: cdf of first licks with exclusions, excluding trials in first 700ms -- saving ALL first lick times (including rxns) and cdfs---------------')
			end

			if strcmpi(collateKey, 'cdf')
				% 
				% 	Will work for an object with any kind of signal
				% 
				[result.ecdf_f, result.ecdf_x, result.result] = sObj.getCDFAnalysis(0, obj.iv.runID, obj.iv.suppressNsave.cdf);	
			elseif strcmpi(collateKey, 'singleTrialFigures')	
				sObj.singleTrialFigures([2,7],200, obj.iv.suppressNsave.singleTrialFigures);
				result = [];
			elseif strcmpi(collateKey, 'movectrlcustom')
				sObj.getBinnedTimeseries(sObj.GLM.gfit, 'custom', [0,2000,3333,3334,7000,17000], 30000);
				if strcmpi(sObj.iv.signaltype_, 'camera')
					smoothing = 3;
				elseif strcmpi(sObj.iv.signaltype_, 'EMG')
					smoothing = 30;
				elseif strcmpi(sObj.iv.signaltype_, 'Photometry')
					smoothing = 100;
				else
					smoothing = 30;
				end
				sObj.plot('LTA2l', [2,4], false, smoothing, 'last-to-first', 1)
				result = [];
				xticks([-5:1:5])
				title([sObj.iv.mousename_ ' ' sObj.iv.daynum_ ' ' num2str(smoothing) 'smsamp'])
				xlim([-5,5])
				figureName = ['MoveCtrl_LTA2l'];
				sObj.suppressNsaveFigure(obj.iv.suppressNsave.figuresHOST, figureName, gcf);
			elseif strcmpi(obj.iv.collateKey,'singleTrialOutliers') 
				lickTimeRange = [0,7];
				redoExclusions = true;
				stdmultiplier = 2;
				outlierTimeRange = [-1.5,7];
				sObj.singleTrialOutliers(lickTimeRange, redoExclusions,stdmultiplier,outlierTimeRange,obj.iv.suppressNsave.singleTrialFigures)
				save('sObj_Corrected.mat', 'sObj', '-v7.3');
				result = [];
			elseif strcmpi(collateKey, 'divergenceIndex')
				sObj.bootOutcomeDivergenceIndex(1000000, 100, 'all')
				sObj.plotDivergenceIndicies(obj.iv.suppressNsave.figuresHOST)
				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				result.Stat = sObj.Stat;
			elseif strcmpi(collateKey, 'PCAcollateSummary')
				sObj.interpolateForPCA;
				sObj.modelDatasetWithPCA(1:3);
				sObj.plotPCA('summary',[],[],[],obj.iv.suppressNsave.figuresHOST)
				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				result.PCA = sObj.GLM.PCA;
				
			elseif strcmpi(collateKey, 'PCAdecoding')
				warning('off','stats:glmfit:IllConditioned');
				% 
				% 	Run PCA analysis
				% 
				sObj.interpolateForPCA;
				sObj.modelDatasetWithPCA(1:3);
				sObj.plotPCA('summary',[],[],[],obj.iv.suppressNsave.PCAsummary)
				if ~isempty(sObj.GLM.PCA.trialidx)
					sObj.plotPCA('testfit', sObj.GLM.PCA.trialidx(1),[],[],obj.iv.suppressNsave.PCAtestfit)
					if numel(sObj.GLM.PCA.trialidx) >= 10
						sObj.plotPCA('testfit', sObj.GLM.PCA.trialidx(10),[],[],obj.iv.suppressNsave.PCAtestfit)
					else
						warning('Found fewer than 10 trials passing PCA critera...')
						nmax = numel(sObj.GLM.PCA.trialidx);
						sObj.plotPCA('testfit', sObj.GLM.PCA.trialidx(nmax),[],[],obj.iv.suppressNsave.PCAtestfit)
					end
					if numel(sObj.GLM.PCA.trialidx) >= 100
						sObj.plotPCA('testfit', sObj.GLM.PCA.trialidx(100),[],[],obj.iv.suppressNsave.PCAtestfit)
					else
						warning('Found fewer than 100 trials passing PCA critera...')
						nmax = numel(sObj.GLM.PCA.trialidx);
						sObj.plotPCA('testfit', sObj.GLM.PCA.trialidx(nmax),[],[],obj.iv.suppressNsave.PCAtestfit)
					end
				else
					warning('Found NO trials passing criteria for PCA analysis. sObj.GLM.PCA.trialidx is empty')
				end
				sObj.plotPCA('wtVsLickTime',[],[],[],obj.iv.suppressNsave.PCAwtVsLickTime)
				sObj.plotPCA('Xfit', 'all',[],[],obj.iv.suppressNsave.PCAXfitAll)
				sObj.plotPCA('Xfit', [1:10:numel(sObj.GLM.cue_s)],[],[],obj.iv.suppressNsave.PCAXfitSelected)
				sObj.binPCAfit('times', 7)
                sObj.plotPCA('Xfitbinned',[1:7],[1:3],[],obj.iv.suppressNsave.Xfitbinned)
% 				sObj.plotPCA('Xfitbinned',[1:7],[],[],obj.iv.suppressNsave.Xfitbinned)
				sObj.PCAmeanSlope(1:3, true, obj.iv.suppressNsave.PCAmeanSlope);
				result.PCA = sObj.GLM.PCA;
				% 
				% 	Run decoding models
				% 
				sObj.Nested_GLM_predictLickTime(1:8, 'ht', false)
				result.decoding.ht = sObj.GLM.decoding;
				result.decoding.ht_raw = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:8, obj.iv.suppressNsave.ht);
				sObj.plotDecodingModelResults('fit',1:8, obj.iv.suppressNsave.ht);

				sObj.Nested_GLM_predictLickTime(1:8, 'htPCA', false)
				result.decoding.htPCA = sObj.GLM.decoding;
				result.decoding.ht_PCA = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:8, obj.iv.suppressNsave.htPCA);
				sObj.plotDecodingModelResults('fit',1:8, obj.iv.suppressNsave.htPCA);

				sObj.Nested_GLM_predictLickTime(1:8, 'PC1-3', false)
				result.decoding.PC1_3 = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:8, obj.iv.suppressNsave.PC1_3);
				sObj.plotDecodingModelResults('fit',1:8, obj.iv.suppressNsave.PC1_3);

				sObj.Nested_GLM_predictLickTime(1:11, 'PC1-3htPCA', false)
				result.decoding.PC1_3htPCA = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:11, obj.iv.suppressNsave.PC1_3htPCA);
				sObj.plotDecodingModelResults('fit',1:11, obj.iv.suppressNsave.PC1_3htPCA);

				sObj.Nested_GLM_predictLickTime(1:6, 'mislope', false)
				result.decoding.mislope = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:6, obj.iv.suppressNsave.mislope);
				sObj.plotDecodingModelResults('fit',1:6, obj.iv.suppressNsave.mislope);

				sObj.Nested_GLM_predictLickTime(1:8, 'htPCA-stiff', false)
				result.decoding.htPCA_stiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:8, obj.iv.suppressNsave.htPCA_stiff);
				sObj.plotDecodingModelResults('fit',1:8, obj.iv.suppressNsave.htPCA_stiff);

				sObj.Nested_GLM_predictLickTime(1:10, 'pretrial_tdtht_1htstiff_PCAversion', false)
				result.decoding.PCA_ht_1_Stiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:10, obj.iv.suppressNsave.PCA_ht_1_Stiff);
				sObj.plotDecodingModelResults('fit',1:10, obj.iv.suppressNsave.PCA_ht_1_Stiff);

				

				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;

				warning('on','stats:glmfit:IllConditioned');	

			elseif strcmpi(collateKey, 'DecodingPaperFinal')
				warning('off','stats:glmfit:IllConditioned');
				disp('-----pretrial_tdtht_1htstiff')
				% 1 stiff tdt, 1 stiff gcamp
				n = 10;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_1htstiff', false)
				result.decoding.ht_stiff_ea = sObj.GLM.decoding;
				if isfield(sObj.GLM.decoding,'flagNoRed')
					result.flagNoRed = sObj.GLM.decoding.flagNoRed;
				else
					result.flagNoRed = false;
				end
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.ht_stiff_ea);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.ht_stiff_ea);
				disp('-----pretrial_tdtht_multiHTstiff')
				% 3 stiff tdt, 3 stiff gcamp
				n = 14;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_multiHTstiff', false)
				result.decoding.multiht_stiff_ea = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.multiht_stiff_ea);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.multiht_stiff_ea);
				disp('-----pretrial_tdtht_PC1-3')
				% 1 stiff tdt, 3PC |weights|
				n = 12;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_PC1-3', false)
				result.decoding.PC1_3_1httdtstiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.PC1_3_1httdtstiff);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.PC1_3_1httdtstiff);
				disp('-----pretrial_tdtht_1htstiff_PCAversion')
				% 1 stiff tdt(not pca), 1stiffPCgcamp
				n = 10;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_1htstiff_PCAversion', false)
				result.decoding.htPCAstiff_1tdtstiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.htPCAstiff_1tdtstiff);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.htPCAstiff_1tdtstiff);
				disp('-----pretrial_tdtht_multiHTstiffPCA')
				% 3 stiff tdt(not pca), 3stiffPCgcamp
				n = 14;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_multiHTstiffPCA', false)
				result.decoding.multihtPCAstiff_multitdtstiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.multihtPCAstiff_multitdtstiff);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.multihtPCAstiff_multitdtstiff);

				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				disp('-----pretrial_PCAtdtht_multiHTstiffPCA')
				% 3 stiff tdt-PCA, 3stiffPCgcamp
				if ~result.flagNoRed
					n = 14;
					sObj.Nested_GLM_predictLickTime(1:14, 'pretrial_PCAtdtht_multiHTstiffPCA', false)
					result.decoding.multihtPCAstiff_multiPCAtdtstiff = sObj.GLM.decoding;
					sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.multihtPCAstiff_multiPCAtdtstiff);
					sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.multihtPCAstiff_multiPCAtdtstiff);
					disp('-----Done.')
				else
					result.decoding.multihtPCAstiff_multiPCAtdtstiff = [];
				end

			elseif strcmpi(obj.iv.collateKey,'1ht_stiff_ea')
				warning('off','stats:glmfit:IllConditioned');
				disp('-----pretrial_tdtht_1htstiff')
				% 1 stiff tdt, 1 stiff gcamp
				n = 10;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_1htstiff', false)
				result.decoding.ht_stiff_ea = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.ht_stiff_ea);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.ht_stiff_ea);
				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				if isfield(sObj.GLM.decoding,'flagNoRed')
					result.flagNoRed = sObj.GLM.decoding.flagNoRed;
				else
					result.flagNoRed = false;
				end
			elseif strcmpi(obj.iv.collateKey,'multiht_stiff_ea')
				warning('off','stats:glmfit:IllConditioned');
				disp('-----pretrial_tdtht_multiHTstiff')
				% 3 stiff tdt, 3 stiff gcamp
				n = 14;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_multiHTstiff', false)
				result.decoding.multiht_stiff_ea = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.multiht_stiff_ea);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.multiht_stiff_ea);
				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				if isfield(sObj.GLM.decoding,'flagNoRed')
					result.flagNoRed = sObj.GLM.decoding.flagNoRed;
				else
					result.flagNoRed = false;
				end
			elseif strcmpi(obj.iv.collateKey,'PC1_3_1httdtstiff')
				warning('off','stats:glmfit:IllConditioned');
				disp('-----pretrial_tdtht_PC1-3')
				% 1 stiff tdt, 3PC |weights|
				n = 12;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_PC1-3', false)
				result.decoding.PC1_3_1httdtstiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.PC1_3_1httdtstiff);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.PC1_3_1httdtstiff);
				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				if isfield(sObj.GLM.decoding,'flagNoRed')
					result.flagNoRed = sObj.GLM.decoding.flagNoRed;
				else
					result.flagNoRed = false;
				end
			elseif strcmpi(obj.iv.collateKey,'1htPCAstiff_1tdtstiff')
				warning('off','stats:glmfit:IllConditioned');
				disp('-----pretrial_tdtht_1htstiff_PCAversion')
				% 1 stiff tdt(not pca), 1stiffPCgcamp
				n = 10;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_1htstiff_PCAversion', false)
				result.decoding.htPCAstiff_1tdtstiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.htPCAstiff_1tdtstiff);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.htPCAstiff_1tdtstiff);
				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				if isfield(sObj.GLM.decoding,'flagNoRed')
					result.flagNoRed = sObj.GLM.decoding.flagNoRed;
				else
					result.flagNoRed = false;
				end
			elseif strcmpi(obj.iv.collateKey,'multihtPCAstiff_multitdtstiff')
				warning('off','stats:glmfit:IllConditioned');
				disp('-----pretrial_tdtht_multiHTstiffPCA')
				% 3 stiff tdt(not pca), 3stiffPCgcamp
				n = 14;
				sObj.Nested_GLM_predictLickTime(1:n, 'pretrial_tdtht_multiHTstiffPCA', false)
				result.decoding.multihtPCAstiff_multitdtstiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.multihtPCAstiff_multitdtstiff);
				sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.multihtPCAstiff_multitdtstiff);
				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				if isfield(sObj.GLM.decoding,'flagNoRed')
					result.flagNoRed = sObj.GLM.decoding.flagNoRed;
				else
					result.flagNoRed = false;
				end
			elseif strcmpi(obj.iv.collateKey,'multihtPCAstiff_multiPCAtdtstiff') 
				warning('off','stats:glmfit:IllConditioned');
				if ~isfield(sObj.GLM, 'tdt')
					errorflag = sObj.addtdt(false);
                    if ~errorflag
                        sObj.GLM.decoding.flagNoRed = false;
                        result.flagNoRed = false;
                    else
                        sObj.GLM.decoding.flagNoRed = true;
                        result.flagNoRed = true;
                    end
                else
                	result.flagNoRed = false;
                end
				disp('-----pretrial_PCAtdtht_multiHTstiffPCA')
				% 3 stiff tdt-PCA, 3stiffPCgcamp
				if ~result.flagNoRed
					n = 14;
					sObj.Nested_GLM_predictLickTime(1:14, 'pretrial_PCAtdtht_multiHTstiffPCA', false)
					result.decoding.multihtPCAstiff_multiPCAtdtstiff = sObj.GLM.decoding;
					sObj.plotDecodingModelResults('summary',1:n, obj.iv.suppressNsave.multihtPCAstiff_multiPCAtdtstiff);
					sObj.plotDecodingModelResults('fit',1:n, obj.iv.suppressNsave.multihtPCAstiff_multiPCAtdtstiff);
					disp('-----Done.')
				else
					result.decoding.multihtPCAstiff_multiPCAtdtstiff = [];
				end
				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				

			elseif strcmpi(collateKey, 'htStiff')
				warning('off','stats:glmfit:IllConditioned');
				sObj.Nested_GLM_predictLickTime(1:8, 'ht-stiff', false)
				result.decoding.ht_stiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:8, obj.iv.suppressNsave.ht_stiff);
				sObj.plotDecodingModelResults('fit',1:8, obj.iv.suppressNsave.ht_stiff);

				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				warning('on','stats:glmfit:IllConditioned');	
			elseif strcmpi(collateKey, 'ht_1_Stiff')
				warning('off','stats:glmfit:IllConditioned');
				sObj.Nested_GLM_predictLickTime(1:10, 'pretrial_tdtht_1htstiff', false)
				result.decoding.ht_stiff = sObj.GLM.decoding;
				sObj.plotDecodingModelResults('summary',1:10, obj.iv.suppressNsave.ht_stiff);
				sObj.plotDecodingModelResults('fit',1:10, obj.iv.suppressNsave.ht_stiff);

				if ~isfield(sObj.GLM,'flick_s_wrtc')
		            sObj.GLM.flick_s_wrtc = nan(size(sObj.GLM.cue_s));
					sObj.GLM.flick_s_wrtc(sObj.GLM.fLick_trial_num) = sObj.GLM.firstLick_s - sObj.GLM.cue_s(sObj.GLM.fLick_trial_num);
				end
				result.flick_s_wrtc = sObj.GLM.flick_s_wrtc;
				warning('on','stats:glmfit:IllConditioned');

			elseif strcmpi(collateKey, 'ht')
				warning('off','stats:glmfit:IllConditioned');
				% 
				%	 obj.horizontalThreshold(Mode = LTA2l or CTA2l, bins = ‘all’, nthresh, delay, direction = ‘+’ (upward crossing), Plot)
				% 
				result = sObj.runHorizontalThresholdCollate(20, obj.iv.runID);
				warning('on','stats:glmfit:IllConditioned');
			elseif strcmpi(collateKey, 'vt')
				warning('off','stats:glmfit:IllConditioned');
				% 
				%	 
				% 
				result = sObj.runVerticalThresholdCollate(70, [-10000:100:7000], obj.iv.runID);
				result.thresholds =[-10000:100:7000];
				warning('on','stats:glmfit:IllConditioned');
			elseif strcmpi(collateKey, 'baselineANOVAidx')
				if ~isfield(sObj.gFitLP, 'nMultibDFF')
					if ~isfield(sObj.GLM, 'rawF')
						sObj.loadRawF;
					end
					sObj.normalizedMultiBaselineDFF(5000, 10, sObj.GLM.rawF);
					save('sObj_Corrected.mat', 'sObj', '-v7.3');
				end
				[result.results, result.F_nm1, result.F_n, result.nm1Score, result.nScore, result.sig_nm1, result.sig_n, result.centers,result.baselineWindow] = sObj.slidingBaselineANOVA('off',false);
			elseif strcmpi(collateKey, 'baselineANOVAwithLick')
				if ~isfield(sObj.gFitLP, 'nMultibDFF')
					if ~isfield(sObj.GLM, 'rawF')
						sObj.loadRawF;
					end
					sObj.normalizedMultiBaselineDFF(5000, 10, sObj.GLM.rawF);
					save('sObj_Corrected.mat', 'sObj', '-v7.3');
				end
				[result.results, result.F_nm1, result.F_n, result.nm1Score, result.nScore, result.sig_nm1, result.sig_n, result.centers,result.baselineWindow] = sObj.slidingBaselineANOVA('include',false);
			elseif strcmpi(collateKey, 'rim')
				result.f = sObj.RIMSummary(sObj.iv.mousename_, sObj.iv.daynum_, sObj.iv.Style);
			else
				error('undefined canned analysis')
			end
			warning('on','stats:glmfit:IllConditioned');
		end



		function ax = plot(obj, Mode, MouseID)
			% 
			% 	Mode options:
			%		Collate Key: cdf
			% 			cdf (3,5 distinguishment)
			% 			cdf-NHl
			% 			hxg (3,5 distinguishment)
			% 			hxg-xsesh (3,5 distinguishment)
            % 			hxg-xsesh-NHL (3,5 distinguishment)
			% 			cdf-xsesh
			% 			hxg-NHL 
			% 
			if nargin < 2
				if strcmpi(obj.iv.collateKey, 'cdf')
					Mode = 'cdf';
				else
					Mode = [];
				end
			end
			if nargin < 3
				MouseID =[];
			end
			if ~isempty(MouseID)
				idxs = contains({obj.collatedResults.sessionID}, MouseID);
			else
				idxs = 1:numel({obj.collatedResults.sessionID});
			end
					
			f = figure;
			if strcmpi(obj.iv.collateKey, 'cdf')
				if strcmpi(Mode, 'cdf') || strcmpi(Mode, 'cdf-NHL')
					rb_s = nan(numel({obj.collatedResults.sessionID}),1);
					rb_s(cellfun(@(x) ~isempty(x), {obj.collatedResults.rb_s})) = [obj.collatedResults.rb_s];
					if strcmpi(Mode, 'cdf-NHL')
						nhlidx = contains({obj.collatedResults.sessionID}, 'NHL');
						rb_s(nhlidx) = -3;
					end
					% 
					% 	Plot cdfs overlaid
					% 
					ax = subplot(1,1,1, 'parent', f);
					hold(ax, 'on');
					xlim(ax, [0, 17]);
					plot(ax, [7, 7], [0, 1],'r-', 'DisplayName', 'ITI Start-3.33 s task')
					title(ax, ['eCDF of First Licks wrt cue ' MouseID])
					for iSet = 1:numel({obj.collatedResults.sessionID})
						if idxs(iSet) == 0
							continue
						end
						if round(rb_s(iSet)) == 3
							pc = 'r-';
						elseif round(rb_s(iSet)) == -3
							pc = 'k-';
							xlim(ax, [0, 17]);
							title(ax, ['eCDF fLicks red=normal, black=NHL ' MouseID])
						elseif round(rb_s(iSet)) == 5
							pc = 'b-';
							xlim(ax, [0, 20]);
							plot(ax, [10, 10], [0, 1],'b-', 'DisplayName', 'ITI Start-5 s task')
						else
							pc = 'm-';
							warning(['undefined reward bound at ' obj.collatedResults(iSet).sessionID])
						end
						plot(ax, obj.collatedResults(iSet).ecdf_x, obj.collatedResults(iSet).ecdf_f, pc, 'displayname', obj.collatedResults(iSet).sessionID)
					end
					
					ylim(ax, [0,1])
					ax.XLabel.String = 'First Lick Time (s wrt cue)';
					ax.YLabel.String = '% of responses in session';
				elseif strcmpi(Mode, 'hxg') || strcmpi(Mode, 'hxg-NHL') || strcmpi(Mode, 'hxg-xsesh') || strcmpi(Mode, 'cdf-xsesh') || strcmpi(Mode, 'hxg-xsesh-NHL')|| strcmpi(Mode, 'cdf-xsesh-NHL') 
					rb_s = nan(numel({obj.collatedResults.sessionID}),1);
					rb_s(cellfun(@(x) ~isempty(x), {obj.collatedResults.rb_s})) = [obj.collatedResults.rb_s];
					if strcmpi(Mode, 'hxg-NHL') || strcmpi(Mode, 'hxg-xsesh-NHL')|| strcmpi(Mode, 'cdf-xsesh-NHL')
						nhlidx = contains({obj.collatedResults.sessionID}, 'NHL');
						rb_s(nhlidx) = -3;
					end
					% 
					% 	Plot hxgs overlaid
					% 
					yy = [0,0];
					ax = subplot(1,1,1, 'parent', f);
					hold(ax, 'on');
					xlim(ax, [0, 7.01]);
					title(ax, ['Normalized First Licks, omitting exclusions' MouseID])
					allTimes_33 = {};
					allTimes_5 = {};
					for iSet = 1:numel({obj.collatedResults.sessionID})
						if idxs(iSet) == 0
							continue
						end
						if round(rb_s(iSet)) == 3
							pc = 'r';
						elseif round(rb_s(iSet)) == -3
							pc = 'k';
							title(ax, ['Normalized First Licks, omitting exclusions, red=normal, black=NHL' MouseID])
						elseif round(rb_s(iSet)) == 5
							pc = 'b';
							xlim(ax, [0, 10.01]);
						else
							pc = 'm';
							warning(['undefined reward bound at ' obj.collatedResults(iSet).sessionID])
						end
						if ~strcmpi(Mode, 'hxg-xsesh') && ~strcmpi(Mode, 'cdf-xsesh') && ~strcmpi(Mode, 'hxg-xsesh-NHL')&& ~strcmpi(Mode, 'cdf-xsesh-NHL') 
							histogram(ax, obj.collatedResults(iSet).f_lick_ex_s_ecdf(obj.collatedResults(iSet).f_lick_ex_s_ecdf>0), [0:0.25:20], 'displaystyle', 'stairs', 'edgecolor', pc, 'displayname', obj.collatedResults(iSet).sessionID, 'normalization', 'probability')
							y = get(ax, 'ylim');
							yy(2) = max(yy(2), y(2));
						else
							if round(rb_s(iSet)) == 3
								allTimes_33{iSet} = obj.collatedResults(iSet).f_lick_ex_s_ecdf(obj.collatedResults(iSet).f_lick_ex_s_ecdf>0);
                            elseif round(rb_s(iSet)) == -3
								allTimes_NHL{iSet} = obj.collatedResults(iSet).f_lick_ex_s_ecdf(obj.collatedResults(iSet).f_lick_ex_s_ecdf>0);
							else
								allTimes_5{iSet} = obj.collatedResults(iSet).f_lick_ex_s_ecdf(obj.collatedResults(iSet).f_lick_ex_s_ecdf>0);
							end
						end
					end
					if strcmpi(Mode, 'cdf-xsesh')
						allTimes_33 = cell2mat(allTimes_33');
						[ecdf_f33,ecdf_x33] = ecdf(allTimes_33);
						allTimes_5 = cell2mat(allTimes_5');
						[ecdf_f5,ecdf_x5] = ecdf(allTimes_5);
						plot(ax, ecdf_x33,ecdf_f33, 'r-', 'displayname', '3.3s total', 'linewidth', 5)
						plot(ax, ecdf_x5,ecdf_f5, 'b-', 'displayname', '5s total', 'linewidth', 5)
						y = get(ax, 'ylim');
						yy(2) = max(yy(2), y(2));
					elseif strcmpi(Mode, 'hxg-xsesh')
						allTimes_33 = cell2mat(allTimes_33');
						allTimes_5 = cell2mat(allTimes_5');
						histogram(ax, allTimes_33, [0:0.25:20], 'displaystyle', 'stairs', 'edgecolor', 'r', 'displayname', '3.3s total', 'normalization', 'probability', 'linewidth', 5)
						histogram(ax, allTimes_5, [0:0.25:20], 'displaystyle', 'stairs', 'edgecolor', 'b', 'displayname', '5s total', 'normalization', 'probability', 'linewidth', 5)
						y = get(ax, 'ylim');
						yy(2) = max(yy(2), y(2));
                    elseif strcmpi(Mode, 'hxg-xsesh-NHL')
                        allTimes_33 = cell2mat(allTimes_33');
						allTimes_NHL = cell2mat(allTimes_NHL');
						histogram(ax, allTimes_33, [0:0.25:20], 'displaystyle', 'stairs', 'edgecolor', 'r', 'displayname', '3.3s total', 'normalization', 'probability', 'linewidth', 5)
						histogram(ax, allTimes_NHL, [0:0.25:20], 'displaystyle', 'stairs', 'edgecolor', 'k', 'displayname', '5s total', 'normalization', 'probability', 'linewidth', 5)
						y = get(ax, 'ylim');
						yy(2) = max(yy(2), y(2));
                    elseif strcmpi(Mode, 'cdf-xsesh-NHL')
						allTimes_33 = cell2mat(allTimes_33');
						allTimes_NHL = cell2mat(allTimes_NHL');
						[ecdf_f33,ecdf_x33] = ecdf(allTimes_33);
						[ecdf_f5,ecdf_x5] = ecdf(allTimes_NHL);
						plot(ax, ecdf_x33,ecdf_f33, 'r-', 'displayname', '3.3s total', 'linewidth', 5)
						plot(ax, ecdf_x5,ecdf_f5, 'k-', 'displayname', '5s total', 'linewidth', 5)
						y = get(ax, 'ylim');
						yy(2) = max(yy(2), y(2));
                    else
                        plot(ax, [7, 7], [yy],'r-', 'DisplayName', 'ITI Start-3.33 s task')
                        plot(ax, [10, 10], [yy],'b-', 'DisplayName', 'ITI Start-5 s task')
					end
					
					
					ax.XLabel.String = 'first lick time (s relative to cue)';
					ax.YLabel.String = 'fraction of responses in session';
				end
				set(ax, 'fontsize', 20)
			else
				error('not Implemented')	
			end
		end

		function inValidIdx = probeBestHT(obj, idx, plotInvalid)
			if nargin < 3
				plotInvalid = false;
			end
			if nargin < 2
				idx = 'all';
			end
			if strcmpi(idx, 'all')
				idx = 1:numel({obj.collatedResults.sessionID});
			end
			collatedResults = obj.collatedResults(idx);
			
			disp('----------------------------------------------------------')
			disp('	Threshold Crossing Analysis Results')
			disp('		')

			f = figure;
			ax_bin = subplot(3,2,1, 'parent', f);
			ax_st = subplot(3,2,2, 'parent', f);
			ax_bin_rsq = subplot(3,2,3, 'parent', f);
			ax_st_rsq = subplot(3,2,4, 'parent', f);
			ax_bin_th = subplot(3,2,5, 'parent', f);
			ax_st_th = subplot(3,2,6, 'parent', f);
			hold(ax_bin, 'on');
			hold(ax_st, 'on');
			hold(ax_bin_rsq, 'on');
			hold(ax_st_rsq, 'on');
			hold(ax_bin_th, 'on');
			hold(ax_st_th, 'on');
			set(ax_bin, 'fontsize', 20);
			set(ax_st, 'fontsize', 20);
			set(ax_bin_rsq, 'fontsize', 20);
			set(ax_st_rsq, 'fontsize', 20);
			set(ax_bin_th, 'fontsize', 20);
			set(ax_st_th, 'fontsize', 20);

			f2 = figure;
			ax_bin_xt2l = subplot(1,2,1, 'parent', f2);
			ax_st_xt2l = subplot(1,2,2, 'parent', f2);
			hold(ax_bin_xt2l, 'on');
			hold(ax_st_xt2l, 'on');
			set(ax_bin_xt2l, 'fontsize', 20);
			set(ax_st_xt2l, 'fontsize', 20);
			
			% legend(ax_bin_rsq,'show', 'interpreter', 'none');
			% legend(ax_st_rsq,'show', 'interpreter', 'none');
			% legend(ax_bin_th,'show', 'interpreter', 'none');
			% legend(ax_st_th,'show', 'interpreter', 'none');

			title(ax_bin, 'Binned hThreshold Crossing-% bins')
			title(ax_st, 'Single-Trial hThreshold Crossing-%trials')
			title(ax_bin_rsq, 'Binned R^2')
			title(ax_st_rsq, 'Single-Trial R^2')
			title(ax_bin_th, 'Binned Slope: xtime vs flick')
			title(ax_st_th, 'Single-Trial Slope: xtime vs flick')
			ax_bin.YLabel.String = '% of bins crossing';
			ax_st.YLabel.String = '% of trials crossing';
			ax_bin_rsq.YLabel.String = 'Rsq';
			ax_st_rsq.YLabel.String = 'Rsq';
			ax_bin_th.YLabel.String = 'Slope';
			ax_st_th.YLabel.String = 'Slope';
			% ax_bin.XLabel.String = 'threshold #';
			% ax_st.XLabel.String = 'threshold #';
			ax_bin_th.XLabel.String = 'threshold #';
			ax_st_th.XLabel.String = 'threshold #';

			title(ax_bin_xt2l, 'Binned xtime - lick time')
			title(ax_st_xt2l, 'Single-Trial xtime - lick time')
			
			ax_bin_xt2l.YLabel.String = 'time (s)';
			ax_bin_xt2l.XLabel.String = 'threshold #';
			ax_st_xt2l.XLabel.String = 'threshold #';

			binColors = 0.1:(.9-0.1)/numel({collatedResults.sessionID}):0.9;
			bin_nValidSesh = 0;
			st_nValidSesh = 0;
			bin_invalidSessions = {};
			st_invalidSessions = {};
			inValidIdx.bin = [];
			inValidIdx.st = [];
			for iSession = 1:numel({collatedResults.sessionID})
				% try
					nthresh = collatedResults(iSession).nthresh;
					if isempty(nthresh)
						warning(['It appears this dataset is missing: ' collatedResults(iSession).sessionID])
						bin_invalidSessions{end+1,1} = collatedResults(iSession).sessionID;
						st_invalidSessions{end+1,1} = collatedResults(iSession).sessionID;
						continue
					end
					% 
					% 	Find trials meeting some critera - let's say at least 100 trials crossing or 10 bins crossing AND median time2xing > 250ms
					% 
					bin_median_xt2l = cellfun(@(x) median(x), collatedResults(iSession).binned_time2lickFromThreshXing);
					st_median_xt2l = cellfun(@(x) median(x), collatedResults(iSession).singletrial_time2lickFromThreshXing);

					if ~plotInvalid
						bin_thresh_not_in_range = unique([find(collatedResults(iSession).binned_pcBinsXing.*collatedResults(iSession).nbins_inRange < 10), find(bin_median_xt2l > -0.25)]); %& bin_median_xt2l < 0.25
						% bin_thresh_not_in_range = bin_thresh_not_in_range(ismember(bin_thresh_not_in_range, find(bin_median_xt2l < 0.25)));
						st_thresh_not_in_range = unique([find(collatedResults(iSession).singletrial_pcTrialsXing.*collatedResults(iSession).nTrials_InRange < 100), find(st_median_xt2l > -0.25)]); % & st_median_xt2l < 0.25
					else
						bin_thresh_not_in_range = [];
						st_thresh_not_in_range = [];
					end
					% st_thresh_not_in_range = st_thresh_not_in_range(ismember(st_thresh_not_in_range, find(st_median_xt2l < 0.25))); % & st_median_xt2l < 0.25

					plot(ax_bin_xt2l, 1:nthresh, bin_median_xt2l, '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['median xtime - lick time: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nbins_inRange)], 'linewidth', 2);
					plot(ax_st_xt2l, 1:nthresh, st_median_xt2l, '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['median xtime - lick time: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nTrials_InRange)], 'linewidth', 2);

					
					bin_maxthreshXtimeLimit250 = find(bin_median_xt2l > -0.025, 1, 'first');
					st_maxthreshXtimeLimit250 = find(st_median_xt2l > -0.025, 1, 'first');
					if isempty(bin_maxthreshXtimeLimit250), bin_maxthreshXtimeLimit250 = nan;, end
					if isempty(st_maxthreshXtimeLimit250), st_maxthreshXtimeLimit250 = nan;, end

					plot(ax_bin, 1:nthresh, collatedResults(iSession).binned_pcBinsXing, '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['%bins: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nbins_inRange)], 'linewidth', 2);
					plot(ax_bin, [bin_maxthreshXtimeLimit250,bin_maxthreshXtimeLimit250], [0,1], '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['250ms latency bound'], 'linewidth', 2);
					plot(ax_st, 1:nthresh, collatedResults(iSession).singletrial_pcTrialsXing, '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['%bins: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nTrials_InRange)], 'linewidth', 2);
					plot(ax_st, [st_maxthreshXtimeLimit250,st_maxthreshXtimeLimit250], [0,1], '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['250ms latency bound'], 'linewidth', 2);
					

					bin_rsq_in_range = collatedResults(iSession).binned_rsq;
					if ~isempty(bin_thresh_not_in_range)
						bin_rsq_in_range(bin_thresh_not_in_range) = nan;
					end
					st_rsq_in_range = collatedResults(iSession).singletrial_rsq;
					if ~isempty(st_thresh_not_in_range)
						st_rsq_in_range(st_thresh_not_in_range) = nan;
					end
					plot(ax_bin_rsq, 1:nthresh, bin_rsq_in_range, '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['Rsq: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nbins_inRange)], 'linewidth', 2);
					plot(ax_st_rsq, 1:nthresh, st_rsq_in_range, '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['Rsq: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nTrials_InRange)], 'linewidth', 2);
					% plot(ax_bin_rsq, [bin_maxthreshXtimeLimit250,bin_maxthreshXtimeLimit250], [0,1], '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['250ms latency bound'], 'linewidth', 2);
					% plot(ax_st_rsq, [st_maxthreshXtimeLimit250,st_maxthreshXtimeLimit250], [0,1], '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['250ms latency bound'], 'linewidth', 2);

					
					th_bin = cell2mat(collatedResults(iSession).binned_b)';
					th_bin = th_bin(:,2);
					th_bin(bin_thresh_not_in_range) = nan;
					th_st = cell2mat(collatedResults(iSession).singletrial_b)';
					th_st = th_st(:,2);
					th_st(st_thresh_not_in_range) = nan;
					plot(ax_bin_th, 1:nthresh, th_bin, '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['r: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nbins_inRange)], 'linewidth', 2);
					plot(ax_st_th, 1:nthresh, th_st, '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['r: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nTrials_InRange)], 'linewidth', 2);
					% plot(ax_bin_th, [bin_maxthreshXtimeLimit250,bin_maxthreshXtimeLimit250], [0,1], '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['250ms latency bound'], 'linewidth', 2);
					% plot(ax_st_th, [st_maxthreshXtimeLimit250,st_maxthreshXtimeLimit250], [0,1], '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['250ms latency bound'], 'linewidth', 2);

				% catch
					
				% end
				if numel(bin_thresh_not_in_range) ~= nthresh
					bin_nValidSesh = bin_nValidSesh + 1;
				else
					bin_invalidSessions{end+1,1} = collatedResults(iSession).sessionID;
					inValidIdx.bin(end+1) = iSession;
				end
				if numel(st_thresh_not_in_range) ~= nthresh
					st_nValidSesh = st_nValidSesh + 1;
				else
					st_invalidSessions{end+1,1} = collatedResults(iSession).sessionID;
					inValidIdx.st(end+1) = iSession;
				end
			end
			ylim(ax_bin, [0,1]);
			ylim(ax_st, [0,1]);
			ylim(ax_bin_rsq, [0,1]);
			ylim(ax_st_rsq, [0,1]);
			ylim(ax_bin_th, [0,1]);
			ylim(ax_st_th, [0,1]);
			xlim(ax_bin, [0,nthresh])
			xlim(ax_st, [0,nthresh])
			xlim(ax_bin_rsq, [0,nthresh])
			xlim(ax_st_rsq, [0,nthresh])
			xlim(ax_bin_th, [0,nthresh])
			xlim(ax_st_th, [0,nthresh])
			% legend(ax_bin,'show', 'interpreter', 'none');
			% legend(ax_st,'show', 'interpreter', 'none');
			
			disp(['	Total Sessions: ' num2str(numel({collatedResults.sessionID}))])
			disp(' ');
			disp(['	Binned Data:'])
			disp(['		# Valid Sessions: ' num2str(bin_nValidSesh)])
			disp(['		Invalid Sessions: '])
            disp(char(bin_invalidSessions))
			disp(' ');
			disp(['	Single Trial Data:'])
			disp(['		# Valid Sessions: ' num2str(st_nValidSesh)])
			disp(['		Invalid Sessions: ' ])
            disp(char(st_invalidSessions))

		end

		function inValidIdx = probeBestVT(obj, idx, plotInvalid)
			if nargin < 3
				plotInvalid = false;
			end
			if nargin < 2
				idx = 'all';
			end
			if strcmpi(idx, 'all')
				idx = 1:numel({obj.collatedResults.sessionID});
			end
			collatedResults = obj.collatedResults(idx);
			
			disp('----------------------------------------------------------')
			disp('	Threshold Crossing Analysis Results')
			disp('		')

			f = figure;
			ax_bin = subplot(3,2,1, 'parent', f);
			ax_st = subplot(3,2,2, 'parent', f);
			ax_bin_rsq = subplot(3,2,3, 'parent', f);
			ax_st_rsq = subplot(3,2,4, 'parent', f);
			ax_bin_r = subplot(3,2,5, 'parent', f);
			ax_st_r = subplot(3,2,6, 'parent', f);
			hold(ax_bin, 'on');
			hold(ax_st, 'on');
			hold(ax_bin_rsq, 'on');
			hold(ax_st_rsq, 'on');
			hold(ax_bin_r, 'on');
			hold(ax_st_r, 'on');
			set(ax_bin, 'fontsize', 20);
			set(ax_st, 'fontsize', 20);
			set(ax_bin_rsq, 'fontsize', 20);
			set(ax_st_rsq, 'fontsize', 20);
			set(ax_bin_r, 'fontsize', 20);
			set(ax_st_r, 'fontsize', 20);

			
			title(ax_bin, 'Binned vThreshold Crossing-% bins')
			title(ax_st, 'Single-Trial vThreshold Crossing-%trials')
			title(ax_bin_rsq, 'Binned R^2')
			title(ax_st_rsq, 'Single-Trial R^2')
			title(ax_bin_r, 'Binned Correlation: xpos vs flick')
			title(ax_st_r, 'Single-Trial Correlation: xpos vs flick')
			ax_bin.YLabel.String = '% of bins crossing';
			ax_st.YLabel.String = '% of trials crossing';
			ax_bin_rsq.YLabel.String = 'Rsq';
			ax_st_rsq.YLabel.String = 'Rsq';
			ax_bin_r.YLabel.String = 'r';
			ax_st_r.YLabel.String = 'r';
			ax_bin_r.XLabel.String = 'threshold #';
			ax_st_r.XLabel.String = 'threshold #';


			binColors = 0.1:(.9-0.1)/numel({collatedResults.sessionID}):0.9;
			bin_nValidSesh = 0;
			st_nValidSesh = 0;
			bin_invalidSessions = {};
			st_invalidSessions = {};
			inValidIdx.bin = [];
			inValidIdx.st = [];
			for iSession = 1:numel({collatedResults.sessionID})
				% try
					nthresh = collatedResults(iSession).nthresh;
					if isempty(nthresh)
						warning(['It appears this dataset is missing: ' collatedResults(iSession).sessionID])
						bin_invalidSessions{end+1,1} = collatedResults(iSession).sessionID;
						st_invalidSessions{end+1,1} = collatedResults(iSession).sessionID;
						continue
					end
					

					

					

					bar(ax_bin, iSession, collatedResults(iSession).nBinsTotal, 'facecolor', [0,0,binColors(iSession)], 'linewidth', 2);
					bar(ax_st, iSession, collatedResults(iSession).nTrialsTotal, 'facecolor', [binColors(iSession),0,0], 'linewidth', 2);
					
					

					bin_rsq_in_range = collatedResults(iSession).binned_rsq;
% 					if ~isempty(bin_thresh_not_in_range)
% 						bin_rsq_in_range(bin_thresh_not_in_range) = nan;
% 					end
					st_rsq_in_range = collatedResults(iSession).singletrial_rsq;
% 					if ~isempty(st_thresh_not_in_range)
% 						st_rsq_in_range(st_thresh_not_in_range) = nan;
% 					end
					plot(ax_bin_rsq, collatedResults(iSession).thresholds, bin_rsq_in_range, '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['Rsq: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nBinsTotal)], 'linewidth', 2);
					plot(ax_st_rsq, collatedResults(iSession).thresholds, st_rsq_in_range, '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['Rsq: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nTrialsTotal)], 'linewidth', 2);
					
					
					plot(ax_bin_r, collatedResults(iSession).thresholds, collatedResults(iSession).binned_r, '-', 'color', [0,0,binColors(iSession)], 'DisplayName', ['r: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nBinsTotal)], 'linewidth', 2);
					plot(ax_st_r, collatedResults(iSession).thresholds, collatedResults(iSession).singletrial_r, '-', 'color', [binColors(iSession),0,0], 'DisplayName', ['r: ' collatedResults(iSession).sessionID, ' nbins: ' num2str(collatedResults(iSession).nTrialsTotal)], 'linewidth', 2);
					
				% if numel(bin_thresh_not_in_range) ~= nthresh
				% 	bin_nValidSesh = bin_nValidSesh + 1;
				% else
				% 	bin_invalidSessions{end+1,1} = collatedResults(iSession).sessionID;
				% 	inValidIdx.bin(end+1) = iSession;
				% end
				% if numel(st_thresh_not_in_range) ~= nthresh
				% 	st_nValidSesh = st_nValidSesh + 1;
				% else
				% 	st_invalidSessions{end+1,1} = collatedResults(iSession).sessionID;
				% 	inValidIdx.st(end+1) = iSession;
				% end
			end
			% ylim(ax_bin, [0,1]);
			% ylim(ax_st, [0,1]);
			ylim(ax_bin_rsq, [0,1]);
			ylim(ax_st_rsq, [0,1]);
			ylim(ax_bin_r, [-1,1]);
			ylim(ax_st_r, [-1,1]);
			
			xlim(ax_bin_rsq, [collatedResults(1).thresholds(1),collatedResults(1).thresholds(end)])
			xlim(ax_st_rsq, [collatedResults(1).thresholds(1),collatedResults(1).thresholds(end)])
			xlim(ax_bin_r, [collatedResults(1).thresholds(1),collatedResults(1).thresholds(end)])
			xlim(ax_st_r, [collatedResults(1).thresholds(1),collatedResults(1).thresholds(end)])
			% legend(ax_bin,'show', 'interpreter', 'none');
			% legend(ax_st,'show', 'interpreter', 'none');
			
			disp(['	Total Sessions: ' num2str(numel({collatedResults.sessionID}))])
			disp(' ');
			disp(['	Binned Data:'])
			disp(['		# Valid Sessions: ' num2str(bin_nValidSesh)])
			disp(['		Invalid Sessions: '])
            disp(char(bin_invalidSessions))
			disp(' ');
			disp(['	Single Trial Data:'])
			disp(['		# Valid Sessions: ' num2str(st_nValidSesh)])
			disp(['		Invalid Sessions: ' ])
            disp(char(st_invalidSessions))
		end


% ------------------------------------------------------------------------------------------------------------------------------------
% 	PCA Collated Results
% ------------------------------------------------------------------------------------------------------------------------------------

		function PCAvarianceExplainedOverlay(obj, d, pcs)
			% 
			% 	d: the index of the datasets to include
			% 	pcs: default is 1:10 for top 10 pcs
			% 
			if nargin < 2 || isempty(d)
				d = find(~obj.reportErrors(false));
			end
			if nargin < 3 || isempty(pcs)
				pcs = 1:10;
			end
			disp('-------------------------------------------------------')
			disp('	Plotting PCA Variance Explained Overlay for the following Sets:')
			disp(char({obj.collatedResults(d).sessionID}'))
			disp(['	pcs: ' num2str(pcs)])

			for id = d

			end
		end







		function plotDecodingModelResults(obj,idx,Mode,n,suppressNsave,Modelnum)
			% 
			% 	Mode: 	'summary'
			% 			'fit'
			% 			'fit-stem'
			% 			'fit-final'
			% 			'fit-final-compressed'
			% 	Same idea, but use refit- to use on refit datasets
			% 			'refit-summary'
			% 			'refit-final-compressed'
			% 
			% 	n: nests to plot
			% 	suppressNsave: a folder to save the plot to. will close the plot and not show (for collation)
			% 
			if nargin < 6
				Modelnum = 1;
            end
            if nargin < 3  || isempty(Mode)
				Mode = 'refit-final-compressed';
			end
			% 
			% 	Extract data
			% 
			if ~contains(Mode, 'refit')
				Dset = obj.collatedResults(idx).decoding;
				names = fieldnames(Dset);
				eval(['decoding = Dset.' names{Modelnum} ';'])
			else
				decoding = obj.collatedResults(idx).refit;
			end
			
			if nargin < 5
				suppressNsave = [];
			end
			if nargin < 4
				n = decoding.n;
			end
			
			colors = ['k','r','g','c','b','m','k','r','g','c','b','m','k','r','g'];
			y = decoding.y;
			yfit = decoding.yfit;

			if contains(Mode, 'summary')
				predictorNames = decoding.predictorNames;
				b = decoding.b;
				CImin = decoding.CImin;
				CImax = decoding.CImax;
				stats = decoding.stats;
				LossImprovement = decoding.LossImprovement;
				BIC = decoding.BIC;
				Rsq = decoding.Rsq;
				ESS = decoding.ESS;
				RSS = decoding.RSS;
				f = figure;
				ax = subplot(1,4,1, 'parent',f);
	            hold(ax, 'on');
	            title(ax,'Log-Coefficients with 95% CI')
	            plot(ax, [0,numel(predictorNames)], [0,0], 'k-')
				ax3 = subplot(1,4,3, 'parent',f);
				title(ax3,'BIC')
				hold(ax3, 'on')
	            ax4 = subplot(1,4,4, 'parent',f);
	            hold(ax4, 'on');
	            title(ax4,'Rsq')
	            for in=n
					plot(ax, [1:in] + (in-1)/(10*max(n)), b{in}, 'o', 'MarkerFaceColor', colors(in),'MarkerEdgeColor', colors(in))
					for inn = 1:in
						plot(ax, [inn + (in-1)/(10*max(n)),inn + (in-1)/(10*max(n))], [CImin{in}(inn),b{in}(inn)], [colors(in) '-'])	%ci min
						plot(ax, [inn + (in-1)/(10*max(n)),inn + (in-1)/(10*max(n))], [b{in}(inn),CImax{in}(inn)], [colors(in) '-'])	%ci max
					end
					plot(ax3, in + (in-1)/(10*max(n)), BIC(in), 'o', 'MarkerFaceColor', colors(in),'MarkerEdgeColor', colors(in))
					plot(ax4, in + (in-1)/(10*max(n)), Rsq(in), 'o', 'MarkerFaceColor', colors(in),'MarkerEdgeColor', colors(in))
				end
	            xticks(ax, 1:max(n))
	            xticklabels(ax, predictorNames(1:max(n)))
	            xticks(ax4, n)
	            xticklabels(ax4, predictorNames(n))
				% ylabel(ax4, 'Log-y fit weights')
				ax2 = subplot(1,4,2, 'parent',f);
				plot(ax2,LossImprovement(1:max(n)), 'ko-')
	            title(ax2,'Loss Improvement')
				xlabel(ax2,'Nest #')
				ylabel(ax2,'Training MSE Loss / Null Loss')
				ylabel(ax, 'Log-y fit weights')
				xticks(ax2, n)
	            xticklabels(ax2, predictorNames(n))
				
				xticks(ax3, n)
	            xticklabels(ax3, predictorNames(n))
	            xtickangle(ax,45)
	            xtickangle(ax2,45)
	            xtickangle(ax3,45)
	            xtickangle(ax4,45)

	            if isempty(suppressNsave)
					disp('		-------------------------------------------')
		 			disp(['		GLM Fit Nested Results for ' decoding.Conditioning ' | ' decoding.Link])
		 			disp(' ')
		 			disp(['		nests shown: ' mat2str(n)])
% 		 			disp(['		ntrials: ' num2str(numel(decoding.y{1})), '/' num2str(numel(obj.GLM.cue_s)) ' (' num2str(100*numel(decoding.y{1})/numel(obj.GLM.cue_s)) '%)'])

		 			disp('	')
		 			disp(cell2table({decoding.predictorNames{1:max(n)};decoding.predictorSet{1:max(n)}}))
	 			end
	 			if ~isempty(suppressNsave), figureName = ['DecodingFitResults_' decoding.predictorKey];, obj.suppressNsaveFigure(suppressNsave, figureName, f), close(f), end
			elseif strcmpi(Mode, 'fit')
	 			for in = n
                    trials_in_range_first = decoding.trials_in_range_first;
	 				f = figure;
		 			subplot(1,2,1);
		 			title(['log-all x,y | Nest ' num2str(in)])
		 			hold on
                    
		 			%  Handle case of PCA where the yfit is not same length across all nests...
		 			if numel(trials_in_range_first) ~= numel(y{in})
		 				trials_in_range_first = decoding.firstnest.trials_in_range_first;
	 				end

                    plot(trials_in_range_first+1,y{in}, 'o-', 'displayname', 'logy')
                    plot(trials_in_range_first+1,yfit{in}, 'o-', 'displayname', 'logyfit')                    
		 			xlabel('trial n')
		 			ylabel('log lick time trial n, (log-s)')
		 			subplot(1,2,2)
		 			hold on
		 			title(['original scale | Nest ' Mode])
                    plot(trials_in_range_first+1,exp(y{in}), 'o-', 'displayname', 'y')
                    plot(trials_in_range_first+1,exp(yfit{in}), 'o-', 'displayname', 'yfit')
		 			xlabel('trial n')
		 			ylabel('lick time trial n, (s)')
		 			legend('show')
		 			if ~isempty(suppressNsave), figureName = ['DecodingFit' decoding.predictorKey '_nest' num2str(in)];, obj.suppressNsaveFigure(suppressNsave, figureName, f), close(f), end
	 			end
 			elseif contains(Mode, 'fit-stem') || contains(Mode, 'fit-final')  || contains(Mode, 'fit-final-compressed')
                for in = n
                    trials_in_range_first = decoding.trials_in_range_first;
                    f = figure;
                    ax = subplot(1,1,1);
                    hold(ax, 'on');
                    title(ax,['original scale | Nest ' Mode])

                    %  Handle case of PCA where the yfit is not same length across all nests...
                    if numel(trials_in_range_first) ~= numel(y{in})
                        trials_in_range_first = decoding.firstnest.trials_in_range_first;
                    end

                    if contains(Mode, 'fit-stem')
	                    stem(trials_in_range_first+1,exp(y{in}), 'displayname', 'y')
	                    stem(trials_in_range_first+1,exp(yfit{in}), 'displayname', 'yfit')
%                     elseif contains(Mode, 'fit-final')
%                     	plot(trials_in_range_first+1,exp(y{in}), 'o-', 'displayname', 'y')
% 	                    plot(trials_in_range_first+1,exp(yfit{in}), 'o-', 'displayname', 'yfit')
                    else
                    	plot(exp(y{in}), 'o-', 'displayname', 'y')
	                    plot(exp(yfit{in}), 'o-', 'displayname', 'yfit')
	                end
                    xlabel('trial n')
                    ylabel('lick time trial n, (s)')
                    legend('show')
                    if ~isempty(suppressNsave), figureName = ['DecodingFit' decoding.predictorKey '_nest' num2str(in)];, obj.suppressNsaveFigure(suppressNsave, figureName, f), close(f), end
		
                end
            end
		end
		function flagCustom(obj)
			% 
			% 	This will specifically exclude datasets where QC was not met in the v3x QC pptx
			% 
			obj.analysis.note{2,1} = 'custom exclusions by QC taken'; 
			exclusions = {'H5_SNc_11',...
			'B5_SNc_17',...
			'B5_SNc_19',...
			'B6_SNc_13',...
			'H3_SNc_17',...
			'H3_SNc_18',...
			'H3_SNc_19',...
			'H3_SNc_20',...
			'H7_SNc_13',...
			'H7_SNc_12',...
			};
			for ii = 1:numel(obj.collatedResults), obj.collatedResults(ii).flagCustom = false; end
			for ii = 1:numel(exclusions)
				obj.collatedResults(find(strcmpi({obj.collatedResults.sessionID}, exclusions{ii}))).flagCustom = true;
			end
		end
		function extractThetas(obj, modelNum, Mode)
			% 
			% 	Pulls out all the theta info for the largest nest from all datasets
			% 
			% 	obj.analysis.thetas = [dataset# x nest]
			% 	obj.analysis.se_th = [dataset# x nest]
			% 	
			% 	Mode = 'fit' or 'refit'
			% 
			obj.flagCustom;	
			if nargin < 3
				Mode = 'fit';
			end
			if nargin < 2
				modelNum = 1; % since most datasets now have one model per dataset...
			end
			fields = fieldnames(obj.collatedResults(1).decoding);
			name = fields{modelNum};
			if strcmpi(Mode, 'fit')
				disp('** Extracting thetas on original fit')
			elseif strcmpi(Mode, 'refit')
                if ~strcmpi(name, 'multihtPCAstiff_multiPCAtdtstiff')
    				disp(['** Extracting thetas on ' num2str(obj.collatedResults(1).refit.k(1)) '-fold xval refit'])
                else
                    disp(['** Extracting thetas on ' num2str(obj.collatedResults(find(~[obj.collatedResults.flagNoRed],1,'first')).refit.k(1)) '-fold xval refit'])
                end
			else
				error('undefined mode')
			end
			try
    			eval(['nMax = numel(obj.collatedResults(1).decoding.' name '.predictorNames);'])
                eval(['numel(obj.collatedResults(1).decoding.' name '.predictorSet);'])
                hasRedFix = false;
            catch
                warning('running for hasred only')
                nMax = 14;
                hasRedFix = true;
                eval(['obj.collatedResults(1).decoding.' name '.predictorNames = obj.collatedResults(find(~[obj.collatedResults.flagNoRed], 1, ''first'')).decoding.' name '.predictorNames;'])
            end
			obj.analysis.note{1,1} = 'sets with < 20 df flagged for poor fit'; 
			obj.analysis.note{3,1} = ['Mode: ' Mode]; 
			obj.analysis.thMode = Mode;
			obj.analysis.setID = cell(numel(obj.collatedResults), 1);
			obj.analysis.thetas = nan(numel(obj.collatedResults), nMax);
			obj.analysis.se_ths = nan(numel(obj.collatedResults), nMax);
			obj.analysis.dfs = nan(numel(obj.collatedResults),1);
			for ii = 1:numel(obj.collatedResults)
				obj.analysis.setID{ii, 1} = obj.collatedResults(ii).sessionID;
                if ~hasRedFix || ~obj.collatedResults(ii).flagNoRed
                    if strcmpi(Mode, 'fit')
                        eval(['obj.analysis.thetas(ii, :) = obj.collatedResults(ii).decoding.' name '.stats{1, nMax}.beta;']) 
                        eval(['obj.analysis.se_ths(ii, :) = obj.collatedResults(ii).decoding.' name '.stats{1, nMax}.se;'])
                        eval(['obj.analysis.dfs(ii,1) = obj.collatedResults(ii).decoding.' name  '.stats{1, nMax}.dfe;'])
                    elseif strcmpi(Mode, 'refit')
                        obj.analysis.thetas(ii, :) = obj.collatedResults(ii).refit.stats{nMax,1}.beta;
                        obj.analysis.se_ths(ii, :) = obj.collatedResults(ii).refit.stats{nMax,1}.se;
                        % if obj.collatedResults(ii).flagNoRed
                        % 	% 
                        % 	% 	Because the tdts are set to be 1's across the board, we need to combine with the offset term to deal with them splitting everything
                        % 	% 
                        % 	tdtIdxs = find(cell2mat(cellfun(@(x) contains(x{1},'tdt'), obj.collatedResults(ii).refit.predictorSet, 'uniformoutput',false)));
                        % 	obj.analysis.thetas(ii, 1) = sum(obj.analysis.thetas(ii, [1, tdtIdxs]));
                        % 	obj.analysis.se_ths(ii, 1) = sum(obj.analysis.se_ths(ii, [1, tdtIdxs]));
                        % 	obj.analysis.thetas(ii, tdtIdxs) = 0;
                        % 	obj.analysis.se_ths(ii, tdtIdxs) = 0;
                        % end
                        % if sum(~isreal(obj.analysis.se_ths(ii, :))) > 0
                        % 	iiii = find(~isreal(obj.analysis.se_ths(ii, :)));
                        % 	disp(['	#' num2str(ii) ' has imaginary error on ths, setting to nan: ' mat2str(iiii)])
                        % 	obj.analysis.se_ths(ii, iiii) = nan;
                        % 	obj.analysis.ths(ii, iiii) = nan;
         %                    obj.analysis.thetas(ii, iiii) = nan;
                        % end
                        obj.analysis.dfs(ii,1) = obj.collatedResults(ii).refit.stats{nMax,1}.dfe;
                    else
                        obj.analysis.thetas(ii, :)
                        obj.analysis.se_ths(ii, :) = nan;
                        obj.analysis.dfs(ii,1) = nan;
                    end
				end
					
				if obj.analysis.dfs(ii,1) < 20
					obj.collatedResults(ii).flagPoorFit = true;
				else
					obj.collatedResults(ii).flagPoorFit = false;
				end
			end		
			if ~hasRedFix
    			eval(['obj.analysis.thNames = obj.collatedResults(1).decoding.' name '.predictorNames;'])
            else
                eval(['obj.analysis.thNames = obj.collatedResults(find(~[obj.collatedResults.flagNoRed], 1, ''first'')).decoding.' name '.predictorNames;'])
            end
		end
		function [meanTh, propagated_se_th, mdf] = getCompositeTheta(obj,idxs, Mode)
			% 
			% 	Called by decodingFigures
			% 
			% 	idxs = the datasets to use
			% 
			if nargin < 3
				Mode = 'fit';
            end
            if numel(fieldnames(obj.collatedResults(1).decoding)) > 1, warning('IS THIS THE RIGHT MODEL VERSION? USING THE FIRST ONE!'), end
			obj.extractThetas(1, Mode)
			if strcmpi(obj.analysis.thMode, 'fit')
				warning('composite theta on ORIGINAL FIT, NOT xval!!!!')
			end
			if nargin < 2
				idxs = 1:numel(obj.collatedResults);
			end
			
			
			
			ths = obj.analysis.thetas(idxs, :);
			se_ths = obj.analysis.se_ths(idxs, :);
			N = numel(idxs);
			NN = N.*ones(1, size(ths, 2));
			
			meanTh = 1/N .* nansum(ths, 1);
			propagated_se_th = 1/N .* sqrt(nansum(se_ths.^2, 1));
			mdf = sum(obj.analysis.dfs(idxs)).*ones(1, size(meanTh,2));
			% 
			% 	Now, handle the thetas with tdt separately
			% 
			thsWithtdt = find(contains(obj.analysis.thNames, 'tdt'));

			tdtIdxs = find(~[obj.collatedResults.flagNoRed]);
			tdtIdxs = tdtIdxs(ismember(tdtIdxs, idxs));
			disp(['tdt betas using only sets with tdt+. Found ' num2str(numel(tdtIdxs)), ' tdt+ sets in range.'])
 			ths = obj.analysis.thetas(tdtIdxs, thsWithtdt);
			se_ths = obj.analysis.se_ths(tdtIdxs, thsWithtdt);
			N = numel(tdtIdxs);
			NN(thsWithtdt) = N;
			meanTh(thsWithtdt) = 1/N .* nansum(ths, 1);
			propagated_se_th(thsWithtdt) = 1/N .* sqrt(nansum(se_ths.^2, 1));
			mdf(thsWithtdt) = sum(obj.analysis.dfs(tdtIdxs)).*ones(1, size(thsWithtdt,2));
			% 
			% 	Now, calculate the CI = b +/- t(0.025, n(m-1))*se
			% 
			for nn = 1:size(meanTh, 2)
				CImin(nn) = meanTh(nn) - abs(tinv(.025,numel(NN(nn))*(mdf(nn) - 1))).*propagated_se_th(nn);
				CImax(nn) = meanTh(nn) + abs(tinv(.025,numel(NN(nn))*(mdf(nn) - 1))).*propagated_se_th(nn);
%                 Tried below, too, but yields same result. Not different
%                 and I think above is correct version
%                 CIminA(nn) = meanTh(nn) - abs(tinv(.025,numel(mdf(nn))*(NN(nn) - 1))).*propagated_se_th(nn);
% 				CImaxA(nn) = meanTh(nn) + abs(tinv(.025,numel(mdf(nn))*(NN(nn) - 1))).*propagated_se_th(nn);
			end

			obj.analysis.flush.meanTh = meanTh;
			obj.analysis.flush.propagated_se_th = propagated_se_th;
			obj.analysis.flush.mdf = mdf;
			obj.analysis.flush.N = NN;
			obj.analysis.flush.CImin = CImin;
			obj.analysis.flush.CImax = CImax;
		end
		function decodingFigures(obj, Mode, Flag,Datasets, ModelNum)
			% 
			% 	Plots composite decoding model figures of all types
			% 
			% 	Mode: 	loss -- plots overlay of all loss plots
			% 			Rsq
			% 			BIC
			% 			b
			% 
			% 			loss-CI -- plots just the 2*STD bars on the loss
			% 			loss-noCI
			% 			Rsq-CI
			% 			Rsq-noCI
			% 			BIC-CI
			% 			b-CI
			% 			loss-condensed (combines the outcome stuff into one feature)
			% 			b-condensed
			% 			b-nob0
			% 			Rsq-overlay
			% 			BIC-overlay
			% 			BIC-unnormalized
			% 
			% 	Add refit- tag to any model to use the REFIT data instead
			% 
			% 	Flag: 	'hasRed' -- gets all sets with red data from ~obj.collatedResults.flagNoRed
			% 
			markers = 20;
			linewidths = 5;
			if nargin < 5 || isempty(ModelNum)
				ModelNum  = 1;
			end
			if nargin < 4 || isempty(Datasets)
				error_idxs = obj.reportErrors(false);
				goodidxs = find(~error_idxs);
			else
				goodidxs = Datasets;
			end
			if nargin < 3 || isempty(Flag)
				Flag = 'none';
			end
			if nargin < 2
				Mode = 'loss-condensed';
			end
			if contains(Mode, 'refit')
				fitMode = 'refit';
				warning('using refit xval''d data')
% 				Mode = erase(Mode, 'refit');
			else
				fitMode = 'fit';
                warning('on')
				warning('USING ORIGINAL FIT, NOT XVAL!')
			end

			if strcmpi(Flag, 'none')
				goodidxs = goodidxs;
			elseif strcmpi(Flag, 'hasRed')
				disp('		* Only using Red+ datasets')
				hasRed = find(~[obj.collatedResults.flagNoRed]);
				goodidxs = goodidxs(ismember(goodidxs, hasRed));
			elseif strcmpi(Flag, 'noRed')
				disp('		* Only using NO Red datasets')
				noRed = find([obj.collatedResults.flagNoRed]);
				goodidxs = goodidxs(ismember(goodidxs, noRed));
			else
				error('Unrecognized Dataset Flag. Options: noRed, hasRed, none');
			end

			% Deal with poor fits
			obj.extractThetas(1,fitMode)
			goodFits = find(~[obj.collatedResults.flagPoorFit]);
			goodidxs = goodidxs(ismember(goodidxs, goodFits));
			goodFits = find(~[obj.collatedResults.flagCustom]);
			goodidxs = goodidxs(ismember(goodidxs, goodFits));
            % 
			% 	Select model
			% 
			ModelType = fieldnames(obj.collatedResults(goodidxs(1)).decoding);
			ModelType = ModelType{ModelNum};
			disp(['Plotting for ' ModelType '.......'])
            if strcmpi(ModelType, 'multihtPCAstiff_multiPCAtdtstiff')
                goodidxs(ismember(goodidxs, find([obj.collatedResults.flagNoRed]))) = [];
            end
			disp(['--------- Plotting ' num2str(numel(goodidxs)) ' Successful Fits --------'])
			[meanTh, propagated_se_th] = obj.getCompositeTheta(goodidxs,fitMode);
			
			eval(['predictorNames = obj.collatedResults(goodidxs(1)).decoding.' ModelType '.predictorNames;'])
			% 
			% 	Decide whether we are refitting or no.
			% 
			if ~contains(Mode, 'refit')
				warning('Using original fit, NOT xval.')
				decoding1 = [obj.collatedResults(goodidxs).decoding];
                for ii = 1:numel(goodidxs)
                   eval(['decoding{ii,1} = decoding1(ii).' ModelType ';'])
                end
            else
				warning('Using XVALIDATED refit.')
				Mode = erase(Mode, 'refit');
				decoding = {obj.collatedResults(goodidxs).refit};
			end
				

			condensedPredictors = find(contains(predictorNames, 'n-1') & ~contains(predictorNames, 'lick time'));
            if strcmpi(Mode, 'b-condensed')
				li = nan(1,numel(predictorNames)-numel(condensedPredictors)+1);
                bi = nan(1,numel(predictorNames)-numel(condensedPredictors)+1);
                predictorNames = predictorNames(2:end);
                condensedPredictors = find(contains(predictorNames, 'n-1') & ~contains(predictorNames, 'lick time'));
            else
               li = nan(1,numel(predictorNames)-numel(condensedPredictors)+1); 
			end
			tdtPredictors = find(contains(predictorNames, 'tdt')); % don't want to include the non-tdt datasets in the averages for the loss or rsq on that term, they should be unchanged
			
			figure
			ax = subplot(1,1,1);
			hold(ax, 'on');
			for idx = 1:numel(goodidxs)
				if strcmpi(Mode, 'loss') || strcmpi(Mode, 'loss-CI') || strcmpi(Mode, 'loss-noCI')
					plot(ax, [0,numel(predictorNames)+1],[1,1], 'k-', 'linewidth', linewidths)
					% eval(['li(idx, 1:numel(predictorNames)) = obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.LossImprovement;'])
					li(idx, 1:numel(predictorNames)) = decoding{idx}.LossImprovement;
					% 
					% 	Check for no red and adjust
					% 
					if obj.collatedResults(idx).flagNoRed
						li(idx, tdtPredictors) = li(idx, tdtPredictors(1)-1);
					end
					ylabel(['Loss ' ModelType],'interpreter', 'none')
					if strcmpi(Mode, 'loss')
						plot(ax, li(idx, 1:numel(predictorNames)), '-', 'color', [0.2,0.2,0.2],'linewidth', 1, 'displayname', obj.collatedResults(goodidxs(idx)).sessionID)
					end
					ylim(ax, [0.5,1])
				elseif strcmpi(Mode, 'loss-condensed')
					% eval(['lipre(idx, 1:numel(predictorNames)) = obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.LossImprovement;'])
					warning('rbf')
					lipre(idx, 1:numel(predictorNames)) = decoding{idx}.LossImprovement;
					if obj.collatedResults(idx).flagNoRed
						lipre(idx, tdtPredictors) = lipre(idx, tdtPredictors(1)-1);
					end
					li(idx, condensedPredictors(1)) = lipre(idx, condensedPredictors(end));
					li(idx, 1:condensedPredictors(1)-1) = lipre(idx, 1:condensedPredictors(1)-1);
					li(idx, condensedPredictors(1)+1:end) = lipre(idx, condensedPredictors(end)+1:end);
					ylabel(['Loss ' ModelType],'interpreter', 'none')
                    plot(ax, li(idx, 1:numel(li(1,:))), '-', 'color', [0.2,0.2,0.2],'linewidth', 1, 'displayname', obj.collatedResults(goodidxs(idx)).sessionID)
                elseif strcmpi(Mode, 'b-condensed')
                	warning('rbf')
					% eval(['li1pre(idx, 1:numel(predictorNames)) = obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.b(numel(predictorNames)+1);'])
					li1pre(idx, 1:numel(predictorNames)) = decoding{idx}.b(numel(predictorNames)+1);
					if obj.collatedResults(idx).flagNoRed
						lipre(idx, tdtPredictors) = lipre(idx, tdtPredictors(1)-1);
					end
                    li(idx, 1:length(li1pre{end})) = li1pre{end};
					% eval(['CImin_i(idx, 1:numel(predictorNames)) = cell2mat(obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.CImin(2:end);'])
     %                eval(['CImax_i(idx, 1:numel(predictorNames)) = cell2mat(obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.CImax(2:end);'])
					bi(idx, condensedPredictors(1)) = nanmean(li(idx, condensedPredictors(end)));
					bi(idx, 1:condensedPredictors(1)-1) = li(idx, 1:condensedPredictors(1)-1);
					bi(idx, condensedPredictors(1)+1:end) = li(idx, condensedPredictors(end)+1:end);
					ylabel(['wt ' ModelType],'interpreter', 'none')
                    plot(ax, bi(idx, 1:numel(bi(1,:))), '-', 'color', [0.2,0.2,0.2],'linewidth', 1, 'displayname', obj.collatedResults(goodidxs(idx)).sessionID)
				elseif strcmpi(Mode, 'Rsq') || strcmpi(Mode, 'Rsq-CI') || strcmpi(Mode, 'Rsq-noCI')
					% eval(['li(idx, 1:numel(predictorNames)) = obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.Rsq;'])
					li(idx, 1:numel(predictorNames)) = decoding{idx}.Rsq;
					if obj.collatedResults(idx).flagNoRed
						li(idx, tdtPredictors) = li(idx, tdtPredictors(1)-1);
					end
					title(['Rsq ' ModelType],'interpreter', 'none')
					if strcmpi(Mode, 'Rsq')
						plot(ax, li(idx, 1:numel(predictorNames)), '-', 'color', [0.2,0.2,0.2],'linewidth', 1,'displayname', obj.collatedResults(goodidxs(idx)).sessionID)
					end
					ylim(ax, [0,0.5])
				elseif strcmpi(Mode, 'Rsq-overlay')
					error('depricated method, use Rsq for this')
					% eval(['lipre(idx, 1:numel(predictorNames)) = obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.Rsq;'])
					lipre(idx, 1:numel(predictorNames)) = decoding{idx}.Rsq;
					if obj.collatedResults(idx).flagNoRed
						lipre(idx, tdtPredictors) = lipre(idx, tdtPredictors(1)-1);
					end
					title(['Rsq ' ModelType],'interpreter', 'none')
					li(idx, condensedPredictors(1)) = lipre(idx, condensedPredictors(end));
					li(idx, 1:condensedPredictors(1)-1) = lipre(idx, 1:condensedPredictors(1)-1);
					li(idx, condensedPredictors(1)+1:end) = lipre(idx, condensedPredictors(end)+1:end);
					plot(ax, li(idx, 1:numel(li(1,:))), '-', 'color', [0.2,0.2,0.2],'linewidth', 1, 'displayname', obj.collatedResults(goodidxs(idx)).sessionID)
					ylim(ax, [0,0.5])
				elseif strcmpi(Mode, 'BIC') || strcmpi(Mode, 'BIC-CI') || strcmpi(Mode, 'BIC-overlay') || strcmpi(Mode, 'BIC-unnormalized')
					% eval(['li(idx, 1:numel(predictorNames)) = obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.BIC;'])
					li(idx, 1:numel(predictorNames)) = decoding{idx}.BIC;
					warning('no noredflag taken into account, showing fit BIC')
				elseif strcmpi(Mode,'b') || strcmpi(Mode, 'b-CI') || strcmpi(Mode, 'b-nob0')
					jitter = (rand*2-1)*0.2;					
					% eval(['bi(idx, 1:numel(predictorNames)) = cell2mat(obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.b(numel(predictorNames)));'])
					% eval(['CImin_i(idx, 1:numel(predictorNames)) = cell2mat(obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.CImin(numel(predictorNames)));'])
     %                eval(['CImax_i(idx, 1:numel(predictorNames)) = cell2mat(obj.collatedResults(' num2str(goodidxs(idx)) ').decoding.' ModelType '.CImax(numel(predictorNames)));'])
     				bi(idx, 1:numel(predictorNames)) = cell2mat(decoding{idx}.b(numel(predictorNames)));
					CImin_i(idx, 1:numel(predictorNames)) = cell2mat(decoding{idx}.CImin(numel(predictorNames)));
                    CImax_i(idx, 1:numel(predictorNames)) = cell2mat(decoding{idx}.CImax(numel(predictorNames)));
					if strcmpi(Mode, 'b')
						for b = 1:numel(predictorNames)
							plot(ax, b+jitter, bi(idx,b), 'k.', 'markersize', markers)
							plot(ax, [b+jitter,b+jitter], [CImin_i(idx,b),CImax_i(idx,b)], 'k-', 'linewidth', linewidths)
						end
					end
				end
			end
			% 
			%	Normalize BIC 
			% 
			if strcmpi(Mode,'BIC') || strcmpi(Mode,'BIC-CI') || strcmpi(Mode,'BIC-overlay')
				for idx = 1:numel(goodidxs)
					row = li(idx, :);
					row = row./max(row);
					% plot(ax, row, 'k-o', 'linewidth', 1)
					li(idx, :) = row;
					plot(ax, li(idx, :),'-', 'color', [0.2,0.2,0.2],'linewidth', 1,'displayname', obj.collatedResults(goodidxs(idx)).sessionID)
				end
				if strcmpi(Mode,'BIC')
					boxplot(ax,li,'extrememode','compress');
				elseif strcmpi(Mode,'BIC-overlay')
					plot(ax, mean(li,1), 'r-o', 'linewidth',linewidths, 'markersize', markers)
				else
					CI = 2*std(li,1);
					plot(ax, mean(li,1), 'r-o', 'linewidth',linewidths, 'markersize', markers)
					for ii = 1:numel(li(1,:))
						obj.plotCIbar(ax,ii,mean(li(:,ii)),CI(ii));
						obj.plotCIbar(ax,ii,mean(li(:,ii)),-CI(ii));
					end
				end
				title(['normalized BIC ' ModelType ' | ' obj.analysis.thMode],'interpreter', 'none')
			elseif strcmpi(Mode,'BIC-unnormalized')
				warning('rbf')
				for idx = 1:numel(goodidxs)
					plot(ax, li(idx, :),'-', 'color', [0.2,0.2,0.2],'linewidth', 1,'displayname', obj.collatedResults(goodidxs(idx)).sessionID)
				end
				plot(ax, mean(li,1), 'r-o', 'linewidth',linewidths, 'markersize', markers)
				title(['unnormalized BIC ' ModelType ' | ' obj.analysis.thMode],'interpreter', 'none')
			elseif strcmpi(Mode,'b') || strcmpi(Mode,'b-CI') || strcmpi(Mode,'b-nob0')
				% 
				% 	Extract all thetas and all SEMs on each theta
				% 


				if strcmpi(Mode,'b-nob0')
                	plot(ax, [0,numel(predictorNames)],[0,0], 'k-', 'linewidth', linewidths)
            	else
            		plot(ax, [0,numel(predictorNames)+1],[0,0], 'k-', 'linewidth', linewidths)
        		end
        		% 	Old version -- incorrect error propagation
				% meanb = nanmean(bi, 1);
				% meanCIn = nanmean(CImin_i, 1);
				% meanCIx = nanmean(CImax_i, 1);
				% 	Below: correct error propagation: ( I propagated sem, then applied correct dof to get the 95% t-test CI)
				meanb = obj.analysis.flush.meanTh;
				meanCIn = obj.analysis.flush.CImin;
				meanCIx = obj.analysis.flush.CImax;
				for b = 1:numel(predictorNames)
					if strcmpi(Mode,'b-nob0')
						if b == 1
							continue
						end
						obj.plotCIbar(ax,b-1,meanb(b),meanCIx(b), true, false,obj.collatedResults(goodidxs(idx)).sessionID);
						obj.plotCIbar(ax,b-1,meanb(b),meanCIn(b), true, false);
						obj.plotCIbar(ax,b-1,meanb(b),[meanCIn(b),meanCIx(b)], true, true);
					else
						obj.plotCIbar(ax,b,meanb(b),meanCIx(b), true, false,obj.collatedResults(goodidxs(idx)).sessionID);
						obj.plotCIbar(ax,b,meanb(b),meanCIn(b), true, false);
						obj.plotCIbar(ax,b,meanb(b),[meanCIn(b),meanCIx(b)], true, true);
					end
					% plot(ax, b, meanb(b), 'r.', 'markersize', 30)
					% plot(ax, [b,b], [meanCIn(b),meanCIx(b)], 'r-', 'linewidth', 3)
				end
				if strcmpi(Mode,'b-CI')
					plot(ax, nanmean(meanb,1), 'r-o', 'linewidth',linewidths, 'markersize', markers,'markerfacecolor','r')
				elseif strcmpi(Mode,'b-nob0')
					plot(ax, nanmean(meanb(:,2:end),1), 'r-o', 'linewidth',linewidths, 'markersize', markers,'markerfacecolor','r')
				end
				ylabel('predictor weight log(s)')
			elseif strcmpi(Mode,'loss-CI') || strcmpi(Mode,'Rsq-CI')
				CI = 2*std(li,1);
				for idx = 1:numel(obj.collatedResults)
					if obj.collatedResults(idx).flagNoRed
						li(idx, tdtPredictors) = nan;
					end
				end
				plot(ax, nanmean(li,1), 'r-o', 'linewidth',linewidths)
				for ii = 1:numel(li(1,:))
					obj.plotCIbar(ax,ii,mean(li(:,ii)),CI(ii));
					obj.plotCIbar(ax,ii,mean(li(:,ii)),-CI(ii));
				end
            elseif contains(Mode,'condensed')
                warning('has errors')
                xxx = nan(1, size(li,2) - numel(condensedPredictors)+1); 
                xxx([1:condensedPredictors(1)-1,condensedPredictors(1)+1:end]) = mean(li(:,[1:condensedPredictors(1)-1,condensedPredictors(end)+1:end]),1);
                xxx(condensedPredictors) = mean(mean(li(condensedPredictors),1));
                plot(ax, xxx, 'r-o', 'linewidth',linewidths, 'markersize', markers, 'markerfacecolor','r')
            else
				plot(ax, mean(li,1), 'r-o', 'linewidth',linewidths, 'markersize', markers, 'markerfacecolor','r')
			end
				
			if strcmpi(Mode,'loss-condensed') || strcmpi(Mode,'b-condensed')
				error('has errors')
                if contains(Mpde, 'loss')
    				plot(ax, [0,numel(predictorNames)+1],[1,1], 'k-', 'linewidth', linewidths)
                else
                    plot(ax, [0,numel(predictorNames)+1],[0,0], 'k-', 'linewidth', linewidths)
                end
                predictorNamesCondensed = cell(1, numel(predictorNames) - numel(condensedPredictors)+1);
				predictorNamesCondensed(1:condensedPredictors(1)-1) = predictorNames(1:condensedPredictors(1)-1);
				predictorNamesCondensed{condensedPredictors(1)} = 'outcome';
				predictorNamesCondensed(condensedPredictors(1)+1:end) = predictorNames(condensedPredictors(end)+1:end);
				xticks(1:numel(predictorNamesCondensed));
	            xticklabels(predictorNamesCondensed); 
				xlim([0,numel(predictorNamesCondensed)+1])
			elseif strcmpi(Mode,'b-nob0')
				xticks(1:numel(predictorNames)-1);
	            xticklabels(predictorNames(2:end)); 
				xlim([0,numel(predictorNames)])
            else
                xticks(1:numel(predictorNames));
	            xticklabels(predictorNames); 
				xlim([0,numel(predictorNames)+1])
			end
			xtickangle(ax,45)
            set(ax, 'fontsize', 50)
            title([Mode, ' ' Flag ' | ' obj.analysis.thMode])
        end
		function plotCIbar(obj,ax, x,y,CIspan, absMode, printstar,handlename)
			% 
			% 	CIspan is the position of the bar end relative to mean. put in absolute mode to go to the position of CI
			% 	printstar takes the full CI span. If it doesn't cross zero, it plots a star 0.1 units above the top bar
			% 		** must use in absMode
			% 
			% 	Plots one side of CI
			% 
			if nargin < 6
				absMode = false;
			end
			if nargin < 7
				printStar = false;
            end
            if nargin < 8
                handlename = [];
            end
            if ~isreal(CIspan)
                warning([handlename ' has imaginary error! This happens because tdt and b0 are identical in the model. Matlab''s glmfit accounts for this, but my xval doesn''t'])
            end
            if absMode && printstar
				if CIspan(1)>0 && CIspan(2)>0 || CIspan(1)<0 && CIspan(2)<0
					plot(ax, x, max(CIspan)+0.2, 'k*', 'markersize',20)
                end
                return
			end
			if absMode && ~printstar
				if ~isempty(handlename)
					plot(ax, [x,x], [y,CIspan], 'k-', 'linewidth',3, 'displayname', handlename)
				else
					plot(ax, [x,x], [y,CIspan], 'k-', 'linewidth',3, 'HandleVisibility', 'off')
				end
				plot(ax, [x-0.15,x+0.15], [CIspan,CIspan], 'k-', 'linewidth',3)
			else
				if ~isempty(handlename)
					plot(ax, [x,x], [y,y + CIspan], 'k-', 'linewidth',3, 'displayname', handlename)
				else
					plot(ax, [x,x], [y,y + CIspan], 'k-', 'linewidth',3, 'HandleVisibility', 'off')
				end
				plot(ax, [x-0.15,x+0.15], [y + CIspan,y + CIspan], 'k-', 'linewidth',3)
            end
		end

		function plotBaselineANOVAidx(obj,sessionIdx,Mode)
			obj.flagPoorANOVA
			if nargin < 3
				Mode = 'mean';
			end
			if nargin < 2 || isempty(sessionIdx)
				% [~,sessionIdx] = obj.reportErrors(false);
    %             sessionIdx = find(sessionIdx);
    			sessionIdx = find(~[obj.collatedResults.FlagQC]);
			elseif strcmpi(sessionIdx, 'mdt')
				sessionIdx = find([obj.collatedResults.MDT]);
                sessionIdx(ismember(find([obj.collatedResults.MDT]),find([obj.collatedResults.FlagQC]))') = [];
			end
			figure
			ax = subplot(1,1,1);
			hold(ax, 'on');
			centers = obj.collatedResults(sessionIdx(1)).centers;
			jitter = rand/10;
			if strcmpi(Mode, 'overlay')
				for ii = 1:numel(sessionIdx)
					
					plot((centers+0.5*0)./1000, obj.collatedResults(sessionIdx(ii)).nm1Score, 'k-')
					plot((centers+0.5*0)./1000, obj.collatedResults(sessionIdx(ii)).sig_nm1+jitter, 'k-', 'linewidth',5)
					plot((centers+0.5*0)./1000, obj.collatedResults(sessionIdx(ii)).nScore, 'r-')
					plot((centers+0.5*0)./1000, obj.collatedResults(sessionIdx(ii)).sig_n-2+jitter, 'r-','linewidth',5)
				end
			elseif strcmpi(Mode, 'mean')
				obj.analysis.seshIdx = [];
				obj.analysis.nm1Scores = {};
				obj.analysis.nScores = {};
				for ii = 1:numel(sessionIdx)
					obj.analysis.nm1Scores{end+1} = obj.collatedResults(sessionIdx(ii)).nm1Score;
					obj.analysis.nScores{end+1} = obj.collatedResults(sessionIdx(ii)).nScore;
					obj.analysis.seshIdx(end+1) = sessionIdx(ii);
                end
                nm1Scoresmean = nanmean(cell2mat(obj.analysis.nm1Scores'));
                nScoresmean = nanmean(cell2mat(obj.analysis.nScores'));
				plot((centers+0.5*0)./1000, nm1Scoresmean, 'k-')
				plot((centers+0.5*0)./1000, nScoresmean, 'r-')
			end
			legend('show')
			xlabel('time (s relative to lamp-off)')
			ylabel('Selectivity Index')
			set(gca, 'fontsize',30)
			xlim([(centers(1)+0.5*0)./1000, (centers(end)+0.5*0)./1000])
			set(gcf,'color','w');
		end
		function [p,tbl,stats,terms] = xMouseBaselineANOVA(obj, sessionIdx, showtable)
			if nargin < 3
				showtable = 'on'
			end
			if showtable == 1
				showtable = 'on';
			elseif showtable == 0
				showtable = 'off';
			end

			% level is 1xn list of labels
			% data is a nx1 list of medians
			A_level = obj.collatedResults(sessionIdx).results{1, 1}.cellData.A_level;
			B_level = obj.collatedResults(sessionIdx).results{1, 1}.cellData.B_level;
			data = obj.collatedResults(sessionIdx).results{1, 1}.cellData.data;
			[p,tbl,stats,terms] = anovan(data, {A_level, B_level}, 'model','interaction', 'display',showtable);
		end
		function flagPoorANOVA(obj, overwrite)
            if nargin < 2
                overwrite = false;
            end
			% 
			% 	Flags sessions based on there being no information in the ANOVA
			% 
			if ~isfield(obj.collatedResults, 'MDT') || overwrite
				obj.collatedResults(1).MDT = [];
				obj.collatedResults(1).FlagQC = [];
				obj.collatedResults(1).Nm1sig = [];
				obj.collatedResults(1).Nsig = [];
				obj.collatedResults(1).NxMsig = [];
				for ii = 1:numel(obj.collatedResults)
					if isempty(obj.collatedResults(ii).results)
						obj.collatedResults(ii).Nm1sig = false;
						obj.collatedResults(ii).Nsig = false;
						obj.collatedResults(ii).NxMsig = false;
						obj.collatedResults(ii).MDT = false;
						obj.collatedResults(ii).FlagQC = true;
					else
						[p,~,~,~] = obj.xMouseBaselineANOVA(ii, false);
						if p(1) < 0.05
							obj.collatedResults(ii).Nm1sig = true;
						else
							obj.collatedResults(ii).Nm1sig = false;
						end
						if p(2) < 0.05
							obj.collatedResults(ii).Nsig = true;
						else
							obj.collatedResults(ii).Nsig = false;
						end
						if p(3) < 0.05
							obj.collatedResults(ii).NxMsig = true;
						else
							obj.collatedResults(ii).NxMsig = false;
						end
						if contains(obj.collatedResults(ii).sessionID, {'H6','H7','B5','B6'})
							obj.collatedResults(ii).MDT = true;
						else
							obj.collatedResults(ii).MDT = false;
						end
						if ~sum([obj.collatedResults(ii).Nm1sig,obj.collatedResults(ii).Nsig,obj.collatedResults(ii).NxMsig])
							obj.collatedResults(ii).FlagQC = true;
						else
							obj.collatedResults(ii).FlagQC = false;
						end
					end
				end
			end
		end
		function flagMDT(obj)
			if ~isfield(obj.collatedResults, 'MDT')
				for ii = 1:numel(obj.collatedResults)
					if contains(obj.collatedResults(ii).sessionID, {'H6','H7','B5','B6'})
						obj.collatedResults(ii).MDT = true;
					else
						obj.collatedResults(ii).MDT = false;
					end
				end
			end
		end
		






















		%%% DECODING XVAL
		function kCheck(obj, d, n, k, modelName)
			if ~isfield(obj.collatedResults, 'kfoldSets')
				obj.collatedResults(1).kfoldSets = [];
			end
			if isempty(obj.collatedResults(d).kfoldSets) || ~isfield(obj.collatedResults(d).kfoldSets, 'k') || obj.collatedResults(d).kfoldSets.k ~= k
				obj.getKfoldSet(d, k, modelName);
			end
		end
		function getKfoldSet(obj, d, k, modelName, startOver)
			disp('Fetching k-fold set...')
			if nargin < 5
				startOver = false;
				% 
				% 	This means that we will recycle an old k-fold set with same k from the obj
				% 
			end
			if startOver || ~isfield(obj.collatedResults(d).kfoldSets, 'k') || obj.collatedResults(d).kfoldSets.k ~= k
				% 
				% 	Clear the previous results
				% 
				obj.collatedResults(d).XvalResults = [];
				if startOver
					disp(['	' num2str(k) '-fold Set is being started again from scratch.'])
				else
					disp(['	' num2str(k) '-fold Set NOT stored in obj. Fetching now.'])
				end
				fetchID = randi(1000000);
				%					
				% 	Randomly divide the set into sets with about equal numbers of trials
				% 
				eval(['nTrials = numel(obj.collatedResults(d).decoding.' modelName '.y{end});'])
				trials_per_set = floor(nTrials/k);
				shuffleSet = randperm(nTrials);
				for iSet = 1:k
					if iSet ~= k
						trialSets{iSet} = shuffleSet((iSet-1)*trials_per_set + 1:iSet*trials_per_set);
					else
						trialSets{iSet} = shuffleSet((iSet-1)*trials_per_set + 1:end);
					end
				end
				obj.collatedResults(d).kfoldSets.ID = fetchID;
				obj.collatedResults(d).kfoldSets.k = k;
				obj.collatedResults(d).kfoldSets.trialSets = trialSets;
				obj.collatedResults(d).kfoldSets.nTrials = nTrials;
			else
				disp(['	' num2str(k) '-fold Set already stored in obj. Recycling.'])
			end
		end
		function selectBestLam(obj, d, n)
			% 
			% 	We will get the best lamda for each cross-validation type (k-fold, code)
			% 		It will get updated for every iteration of xvalidate
			%
		 	%	obj.XvalResults(d).nest(n)
		 	% 
		 	% 		. kfoldID
		 	% 		. k
		 	% 		meanTestLoss
		 	% 		testLoss (as fx of lamda) -- so a matrix of k x #lam
		 	% 		bestLam
		 	% 		lamRange
		 	% ----------------------------
		 	% 
		 	% 	Get optimal lam
		 	%	obj.collatedResults(d).XvalResults(n).
		 	% 
		 	[lams, sortIdx] = sort([obj.collatedResults(d).XvalResults(n).models.lam]);
		 	nullLossTraining = obj.collatedResults(d).XvalResults(n).models(sortIdx).nullLossTraining;
			nullLossTest = obj.collatedResults(d).XvalResults(n).models(sortIdx).nullLossTest;
			trainingLoss = [obj.collatedResults(d).XvalResults(n).models(sortIdx).trainingLoss];
			testLoss = [obj.collatedResults(d).XvalResults(n).models(sortIdx).testLoss];
			% 
			% 	Deal with ill-specified models, ie if matrix was singular, we should exclude it from consideration
			% 	by making loss effectively infinite
			% 
			trainingLoss(isnan(trainingLoss)) = 10^10;
			testLoss(isnan(testLoss)) = 10^10;
			
			meanTrainingLoss = mean(trainingLoss, 1);
			meanTestLoss = mean(testLoss, 1);


		 	bestLamIdx = find(meanTestLoss == min(meanTestLoss));
		 	bestLam = lams(bestLamIdx);
            bestTestLoss = meanTestLoss(bestLamIdx);
            bestTrainingLoss = meanTrainingLoss(bestLamIdx);

			k = obj.collatedResults(d).kfoldSets.k;
	 		obj.collatedResults(d).XvalResults(n).kfoldID = obj.collatedResults(d).XvalResults(n).models(1).kfoldID;
            obj.collatedResults(d).XvalResults(n).bestLam = bestLam;
            obj.collatedResults(d).XvalResults(n).bestTestLoss = bestTestLoss;
		 	obj.collatedResults(d).XvalResults(n).k = k;
		 	obj.collatedResults(d).XvalResults(n).meanTestLoss = meanTestLoss;
		 	obj.collatedResults(d).XvalResults(n).testLoss = testLoss;
		 	obj.collatedResults(d).XvalResults(n).meanTrainingLoss = meanTrainingLoss;
		 	obj.collatedResults(d).XvalResults(n).trainingLoss = trainingLoss;
		 	obj.collatedResults(d).XvalResults(n).bestTrainingLoss = bestTrainingLoss;
		 	obj.collatedResults(d).XvalResults(n).lamRange = lams;
		 	obj.collatedResults(d).XvalResults(n).nullLossTraining = nullLossTraining;
		 	obj.collatedResults(d).XvalResults(n).nullLossTest = nullLossTest;
		 	obj.collatedResults(d).XvalResults(n).bestLamIdx = bestLamIdx;
		end
		function xValLossPlotVsLam(obj, d, n)
			figure
			ax2 = axes;
			hold(ax2, 'on');
			lamLabel = cell(1, length(obj.collatedResults(d).XvalResults(n).models));
			for ilam = 1:length(obj.collatedResults(d).XvalResults(n).models)
				lamLabel{ilam} = num2str(obj.collatedResults(d).XvalResults(n).models(ilam).lam);
				nullLossTraining(ilam) =  mean([obj.collatedResults(d).XvalResults(n).models(ilam).nullLossTraining]);
				nullLossTest(ilam) =  mean([obj.collatedResults(d).XvalResults(n).models(ilam).nullLossTest]);
			end


			% Have to sort the lams before plotting it all
			[lams, sortIdx] = sort([obj.collatedResults(d).XvalResults(n).models.lam]);
			lamLabel = lamLabel(sortIdx);
			nullLossTraining = nullLossTraining(sortIdx);
			nullLossTest = nullLossTest(sortIdx);
			trainingLoss = [obj.collatedResults(d).XvalResults(n).models(sortIdx).trainingLoss];
			testLoss = [obj.collatedResults(d).XvalResults(n).models(sortIdx).testLoss];
			meanTrainingLoss = mean(trainingLoss, 1);
			meanTestLoss = mean(testLoss, 1);
			plot(ax2, 1:numel(lams), nullLossTraining, 'b--', 'DisplayName', 'Mean Null Loss - Training')
			plot(ax2, 1:numel(lams), nullLossTest, 'r--', 'DisplayName', 'Mean Null Loss - Test')
			for ilam = 1:numel(lams)
				l = scatter(ax2, ilam.*ones(size(trainingLoss(:, ilam))), trainingLoss(:, ilam), 'b', 'filled', 'HandleVisibility', 'off');
				alpha(l, 0.2)
				l = scatter(ax2, ilam.*ones(size(testLoss(:, ilam))), testLoss(:, ilam), 'r', 'filled', 'HandleVisibility', 'off');
				alpha(l, 0.2)
			end

			plot(ax2, 1:numel(lams), meanTrainingLoss, 'b-o', 'DisplayName', 'Training Loss')
			plot(ax2, 1:numel(lams), meanTestLoss, 'r-o', 'DisplayName', 'Test Loss')
			xticks(ax2, 1:length(obj.collatedResults(d).XvalResults(n).models));
			xticklabels(ax2, lamLabel);
			title(ax2, 'Loss vs Regularization')
			xlabel(ax2, 'lambda')
			ylabel(ax2, 'loss')
			legend(ax2, 'show')
			xtickangle(ax2,90);
		end
		function plotXval(obj, d, n)
			models = obj.collatedResults(d).XvalResults(n).models;
			figure,
			ax_trainingLoss = subplot(2,2,1);
			ax_testLoss = subplot(2,2,2);
			ax_EnImprovement = subplot(2,2,3);
			ax_EImprovement = subplot(2,2,4);
			% ax_se_model = subplot(2,2,3);
			% ax_R2 = subplot(2,2,4);
			hold(ax_trainingLoss, 'on')
			hold(ax_testLoss, 'on')
			hold(ax_EnImprovement, 'on')
			hold(ax_EImprovement, 'on')
			% hold(ax_se_model, 'on')
			% hold(ax_R2, 'on')

			nXv = numel(models(1).se_model);
			jitter = 0.5*rand(nXv,1) -.25;

			lamLabel = cell(1, length(models));

			for ilam = 1:length(models)
				lamLabel{ilam} = num2str(models(ilam).lam);
				nullLossTraining(ilam) =  mean([models(ilam).nullLossTraining]);
				nullLossTest(ilam) =  mean([models(ilam).nullLossTest]);

				scatter(ax_trainingLoss, jitter + ilam*ones(nXv, 1), models(ilam).trainingLoss)
				alpha(.1);
				scatter(ax_trainingLoss, ilam, models(ilam).meanTrainingLoss, 'filled')

				scatter(ax_testLoss, jitter + ilam*ones(nXv, 1), models(ilam).testLoss)
				alpha(.1);
				scatter(ax_testLoss, ilam, models(ilam).meanTestLoss, 'filled')

				EnImprovement{ilam} = 100*(nullLossTraining(ilam) - models(ilam).trainingLoss)./nullLossTraining(ilam);
				meanEnImprovement(ilam) = mean(EnImprovement{ilam});
				EImprovement{ilam} = 100*(nullLossTest(ilam) - models(ilam).testLoss)./nullLossTest(ilam);
				meanEImprovement(ilam) = mean(EImprovement{ilam});

				scatter(ax_EnImprovement, jitter + ilam*ones(nXv, 1), EnImprovement{ilam})
				alpha(.1);

				scatter(ax_EImprovement, jitter + ilam*ones(nXv, 1), EImprovement{ilam})
				alpha(.1);

				% scatter(ax_se_model, jitter + ilam*ones(nXv, 1), models(ilam).se_model)
				% alpha(.3);

				% scatter(ax_R2, jitter + ilam*ones(nXv, 1), models(ilam).R2)
				% alpha(.3);
			end
			plot(ax_trainingLoss, nullLossTraining, 'k-o')
			plot(ax_testLoss, nullLossTest, 'k-o')
			plot(ax_EnImprovement, meanEnImprovement, 'k-o')
			plot(ax_EImprovement, meanEImprovement, 'k-o')


			xticks(ax_trainingLoss, 1:length(models));
			xticklabels(ax_trainingLoss, lamLabel);
			title(ax_trainingLoss, 'Training Loss')
			xlabel(ax_trainingLoss,'lamda')
			ylabel(ax_trainingLoss,'Training Set MSELoss')

			
			xticks(ax_testLoss, 1:length(models));
			xticklabels(ax_testLoss, lamLabel);
			title(ax_testLoss, 'Test Loss')
			xlabel(ax_testLoss,'lamda')
			ylabel(ax_testLoss,'Test Set MSELoss')

			xticks(ax_EnImprovement, 1:length(models));
			xticklabels(ax_EnImprovement, lamLabel);
			title(ax_EnImprovement, 'Training Loss % Improvement vs Sn Null')
			xlabel(ax_EnImprovement,'lamda')
			ylabel(ax_EnImprovement,'% Improvement (+ is better)')

			xticks(ax_EImprovement, 1:length(models));
			xticklabels(ax_EImprovement, lamLabel);
			title(ax_EImprovement, 'Test Loss % Improvement vs Sn Null')
			xlabel(ax_EImprovement,'lamda')
			ylabel(ax_EImprovement,'% Improvement (+ is better)')

			% xticks(ax_se_model, [models.lam]);
			% xticklabels(ax_se_model, lamLabel);
			% title(ax_se_model, 'se of model')

			% xticks(ax_R2, [models.lam]);
			% xticklabels(ax_R2, lamLabel);
			% title(ax_R2, 'R^2')
			models = models;
			figure
			ax2 = axes;
			hold(ax2, 'on');
			plot(ax2, nullLossTraining, 'b--', 'DisplayName', 'Mean Null Loss - Training')
			plot(ax2, nullLossTest, 'r--', 'DisplayName', 'Mean Null Loss - Test')

			plot(ax2, [models.meanTrainingLoss], 'b-o', 'DisplayName', 'Training Loss')
			plot(ax2, [models.meanTestLoss], 'r-o', 'DisplayName', 'Test Loss')
			xticks(ax2, 1:length(models));
			xticklabels(ax2, lamLabel);
			title(ax2, 'Loss vs Regularization')
			xlabel(ax2, 'lambda')
			ylabel(ax2, 'loss')
			legend(ax2, 'show')
			xtickangle(ax2,90);

		end
		function ax = plotBestLoss(obj, d, kfoldID, ax)
			% 
			% 	Will make training and test loss plots for the dataset in question
			% 
			if nargin < 2
				error('Must specify dataset #')
			end
			if nargin < 4
				figure, 
				ax{1} = subplot(2,2,1);
				title(ax{1}, 'Training Loss by Nest')
				xlabel(ax{1}, 'Nest #')
				ylabel(ax{1}, 'MSE Loss')
				hold(ax{1}, 'on');
				ax{2} = subplot(2,2,2);
				title(ax{2}, 'Test Loss by Nest')
				xlabel(ax{2}, 'Nest #')
				ylabel(ax{2}, 'MSE Loss')
				hold(ax{2}, 'on');
				ax{3} = subplot(2,2,3);
				title(ax{3}, 'Training Loss Improvement')
				xlabel(ax{3}, 'Nest #')
				ylabel(ax{3}, '% Improvement vs Training Null')
				hold(ax{3}, 'on');
				ax{4} = subplot(2,2,4);
				title(ax{4}, 'Test Loss Improvement')
				xlabel(ax{4}, 'Nest #')
				ylabel(ax{4}, '% Improvement vs Test Null')
				hold(ax{4}, 'on');
			end

			
			
			XvalResults = obj.collatedResults(d).XvalResults;
			nests = 1:numel(XvalResults(end).nests);
			nNests = numel(XvalResults);
			bestTestLoss = nan(nNests, 1);
			bestTrainingLoss = nan(nNests, 1);
			testLoss = cell(nNests, 1);
			trainingLoss = cell(nNests, 1);
			meanPercImprovementTrain = nan(nNests, 1);
			meanPercImprovementTest = nan(nNests, 1);
			percImprovementTrain = cell(nNests, 1);
			percImprovementTest = cell(nNests, 1);

			for n = 1:nNests
				k = XvalResults(n).k;
				bestLamIdx = XvalResults(n).bestLamIdx;
				bestTestLoss(n) = XvalResults(n).bestTestLoss;
				bestTrainingLoss(n) = XvalResults(n).bestTrainingLoss;
				trainingLoss{n} = XvalResults(n).trainingLoss(:,bestLamIdx);
				testLoss{n} = XvalResults(n).testLoss(:,bestLamIdx);
				meanPercImprovementTrain(n) = 100*(mean(XvalResults(n).nullLossTraining) - bestTrainingLoss(n))/mean(XvalResults(n).nullLossTraining);
				meanPercImprovementTest(n) = 100*(mean(XvalResults(n).nullLossTest) - bestTestLoss(n))/mean(XvalResults(n).nullLossTest);
				for il = 1:numel(trainingLoss{n})
					percImprovementTrain{n}(end+1) = 100*(XvalResults(n).nullLossTraining(il) - trainingLoss{n}(il))/XvalResults(n).nullLossTraining(il);
					percImprovementTest{n}(end+1) = 100*(XvalResults(n).nullLossTest(il) - testLoss{n}(il))/XvalResults(n).nullLossTest(il);
				end

				l = scatter(ax{1}, nests(n).*ones(size(trainingLoss{n})), trainingLoss{n}, 'b', 'filled', 'HandleVisibility', 'off');
				alpha(l, 0.1);
				l = scatter(ax{2}, nests(n).*ones(size(testLoss{n})), testLoss{n}, 'r', 'filled', 'HandleVisibility', 'off');
				alpha(l, 0.1);
				l = scatter(ax{3}, nests(n).*ones(size(percImprovementTrain{n})), percImprovementTrain{n}, 'b', 'filled', 'HandleVisibility', 'off');
				alpha(l, 0.1);
				l = scatter(ax{4}, nests(n).*ones(size(percImprovementTest{n})), percImprovementTest{n}, 'r', 'filled', 'HandleVisibility', 'off');
				alpha(l, 0.1);
			end
			plot(ax{1}, nests, mean(XvalResults(n).nullLossTraining).*ones(size(nests)), 'b--', 'DisplayName', 'Null Training Loss')
			plot(ax{1}, nests, bestTrainingLoss, 'b-o', 'DisplayName', 'Best Training Loss')
			legend(ax{1}, 'show')

			plot(ax{2}, nests, mean(XvalResults(n).nullLossTest).*ones(size(nests)), 'r--', 'DisplayName', 'Null Test Loss')
			plot(ax{2}, nests, bestTestLoss, 'r-o', 'DisplayName', 'Best Test Loss')
			legend(ax{2}, 'show')

			plot(ax{3}, nests, meanPercImprovementTrain, 'b-o')
			plot(ax{4}, nests, meanPercImprovementTest, 'r-o')

			xticks(ax{1}, 1:max(nests));
			xticks(ax{2}, 1:max(nests));
			xticks(ax{3}, 1:max(nests));
			xticks(ax{4}, 1:max(nests));

		end
		function autoXval(obj, d, n, modelName, k, nLamsPerRound, useSingle)
            nRounds = 1;
			% 
			% 	If you set nRounds to 0, this will just return the plot summary! nice
			% 
			% 	d = the collatedResults index for the dataset
			% 	n = the nest of the model to xval
			% 	
			suppressPlot = true;
			if nargin < 7
				useSingle = false;
			end
			if nargin<6
				nLamsPerRound = 20;
			end
			if nargin < 5;
				k = 5;
			end
			
			disp('~~~~~~~~~~~~~~~~~~')
			disp('~	Auto Xvalidate ~')
			disp('~~~~~~~~~~~~~~~~~~')
			% 
			% 	The function will work to choose a lam that is best 
			% 
			% 	while min is 0 or the largest lam
			% 	1. Test new range of lam and find best lam
			% 	for nRounds after finding a middle-lam that is minimum test loss,
			% 	2. Test range of lams between min lam and surrounding lams
			% 
			bestLamIdx = 1;
			% 
			% 	Check for existing k-set for this dataset
			% 	
			obj.kCheck(d,n,k,modelName);
			kfoldID = obj.collatedResults(d).kfoldSets.ID;
			% 
			% 	Recycle any lams already in our saved results for this kID
			% 
			if ~isfield(obj.collatedResults, 'XvalResults')
				obj.collatedResults.XvalResults = [];
			end
			if ~isempty(obj.collatedResults(d).XvalResults) && isfield(obj.collatedResults(d).XvalResults, 'n') && numel(obj.collatedResults(d).XvalResults.n) >= n && isfield(obj.collatedResults(d).XvalResults.n, 'kfoldID') && ~isempty(obj.collatedResults(d).XvalResults(n).kfoldID) && obj.collatedResults(d).XvalResults(n).kfoldID == kfoldID
				lams = obj.collatedResults(d).XvalResults(n).lamRange;
				obj.selectBestLam(d, n);
				bestLamIdx = obj.collatedResults(d).XvalResults(n).bestLamIdx;
			else
				lams = [0,0.001, 0.1,1,10,100,1000,10000];
				obj.xValidate(d, n, lams, 'k-fold', k, false, modelName, suppressPlot, useSingle);
				obj.selectBestLam(d, n);
				bestLamIdx = obj.collatedResults(d).XvalResults(n).bestLamIdx;
            end
            if numel(bestLamIdx)>1
                bestLamIdx = bestLamIdx(2);
            end
            iter = 1;
			while bestLamIdx==1 || bestLamIdx==2 || bestLamIdx==numel(lams) && iter <= 10								
				if bestLamIdx==1
					disp('lam=0 is still best.')
                    sortLams = sort(lams);
					lams = [lams, sortLams(2)/10];
				elseif bestLamIdx==2
					disp('lam=0 is still flanking on left.')
                    sortLams = sort(lams);
					lams = [lams, sortLams(2)/10];
				elseif bestLamIdx==numel(lams)
					disp(['lam=' num2str(lams(end)) ' is best.'])
					lams = [lams, lams(end)*10];
				end
				disp(['	Attempting with lam=' num2str(lams(end)) '... ' datestr(now)])
				obj.xValidate(d, n, lams, 'k-fold', k, false, modelName, suppressPlot, useSingle);
				obj.selectBestLam(d, n);
				lams = obj.collatedResults(d).XvalResults(n).lamRange;
				bestLamIdx = obj.collatedResults(d).XvalResults(n).bestLamIdx;
                iter = iter + 1;
                if numel(bestLamIdx)>1
                    bestLamIdx = bestLamIdx(2);
                end
				if iter > 10
					warning('Timeout - ran 10 iterations of lams and found nothing...')
                    
                    if bestLamIdx==1 || bestLamIdx==2
                        disp('lam=0 or effectively 0 is still best. Finish this round')
                        nRounds = 0;
                        break
                    else
                        disp('lam is exploding, clearly the fit is not good. proceed with largest lam')
                        nRounds = 0;
                        break
                    end
				end
			end
			disp(['----- Found a minimum lam=' num2str(lams(bestLamIdx)) ' | ' datestr(now)])
			disp([' '])
			% 
			% 	Now zoom in and find even better lam
			% 
			for l = 1:nRounds
				disp(['----- Testing subdiv lams, round ' num2str(l) '/' num2str(nRounds) ' | ' datestr(now)])
				[lamsSorted, idxs] = sort(lams);
				bestLam = lams(bestLamIdx);
				bestLamIdxSorted = find(lamsSorted == bestLam);
				if bestLamIdxSorted == 2
					% 
					% 	The range is between [0, bestlam, nextbest]
					% 
					lams = unique([lams, lamsSorted(bestLamIdxSorted)/10:lamsSorted(bestLamIdxSorted)/(nLamsPerRound/2):lamsSorted(bestLamIdxSorted), lamsSorted(bestLamIdxSorted):lamsSorted(bestLamIdxSorted+1)/(nLamsPerRound/2):lamsSorted(bestLamIdxSorted+1)]);
					obj.xValidate(d, n, lams, 'k-fold', k, false, modelName, suppressPlot, useSingle);
					obj.selectBestLam(d, n);
					bestLamIdx = obj.collatedResults(d).XvalResults(n).bestLamIdx;
					lams = obj.collatedResults(d).XvalResults(n).lamRange;
                else
					lams = unique([lams, lamsSorted(bestLamIdxSorted-1):(lamsSorted(bestLamIdxSorted)-lamsSorted(bestLamIdxSorted-1))/(nLamsPerRound/2):lamsSorted(bestLamIdxSorted), lamsSorted(bestLamIdxSorted):(lamsSorted(bestLamIdxSorted+1)-lamsSorted(bestLamIdxSorted))/(nLamsPerRound/2):lamsSorted(bestLamIdxSorted+1)]);
					obj.xValidate(d, n, lams, 'k-fold', k, false, modelName, suppressPlot, useSingle);
					obj.selectBestLam(d, n);
					bestLamIdx = obj.collatedResults(d).XvalResults(n).bestLamIdx;
					lams = obj.collatedResults(d).XvalResults(n).lamRange;
				end
			end
			% 
			% The results struct is always sorted by lam, so we can access it without fear here;
			% 
			[lamsSorted, idxs] = sort(lams);
            nullLossTest = obj.collatedResults(d).XvalResults(n).nullLossTest;
			bestLam = obj.collatedResults(d).XvalResults(n).bestLam;
            if numel(bestLam) > 1
            	obj.collatedResults(d).XvalResults(n).flag = {'more than one bestLam, removed extras', bestLam};
                bestLam = bestLam(1);
                obj.collatedResults(d).XvalResults(n).bestLam = obj.collatedResults(d).XvalResults(n).bestLam(1);
                obj.collatedResults(d).XvalResults(n).bestLamIdx = obj.collatedResults(d).XvalResults(n).bestLamIdx(1);
                obj.collatedResults(d).XvalResults(n).bestTestLoss = obj.collatedResults(d).XvalResults(n).bestTestLoss(1);
                obj.collatedResults(d).XvalResults(n).bestTrainingLoss = obj.collatedResults(d).XvalResults(n).bestTrainingLoss(1);
            end
			bestLamIdxSorted = find(lamsSorted == bestLam);
            assert(obj.collatedResults(d).XvalResults(n).bestLamIdx == bestLamIdxSorted);
			bestTestLoss = obj.collatedResults(d).XvalResults(n).bestTestLoss;
			testLoss = [obj.collatedResults(d).XvalResults(n).testLoss(:,bestLamIdxSorted)];
			meanPercImprovementTest = 100*(mean(nullLossTest) - bestTestLoss)/mean(nullLossTest);
			percImprovementTest = nan(1, k);
			for il = 1:k
				percImprovementTest(il) = 100*(nullLossTest(il) - testLoss(il))/nullLossTest(il);
			end
			% 
			% 	Report the results
			% 
			disp('Complete. --------------')
			disp('')
			disp('Report: ')
			disp(['Lams tested: ' mat2str(sort(obj.collatedResults(d).XvalResults(n).lamRange))])
			disp(['Best lam: ' mat2str(bestLam)])
			disp(['Test Loss Improvement over Null: ' mat2str(round(meanPercImprovementTest,8)) ' % Improvement'])
			disp(['Range Test Loss Improvement over Null Range: ' mat2str(round(sort(percImprovementTest),8)) ' % Improvement'])
			disp('')
			
            if bestLamIdxSorted ~= 1
                llIdx = idxs(bestLamIdxSorted - 1);
            else
                llIdx = idxs(bestLamIdxSorted);
            end
            if bestLamIdxSorted ~= numel(lamsSorted)
    			ulIdx = idxs(bestLamIdxSorted + 1);
            else
                ulIdx = idxs(bestLamIdxSorted);
            end
			
			bestLam_ll = lams(llIdx);
		 	testLoss_ll = [obj.collatedResults(d).XvalResults(n).testLoss(:,llIdx)];
		 	meanTestLoss_ll = mean(testLoss_ll, 1);
		 	meanPercImprovementTest_ll = 100*(mean(nullLossTest) - meanTestLoss_ll)/mean(nullLossTest);
			percImprovementTest_ll = nan(1, k);

			bestLam_ul = lams(ulIdx);
		 	testLoss_ul = [obj.collatedResults(d).XvalResults(n).testLoss(:,ulIdx)];
		 	meanTestLoss_ul = mean(testLoss_ul, 1);
		 	meanPercImprovementTest_ul = 100*(mean(nullLossTest) - meanTestLoss_ul)/mean(nullLossTest);			
			percImprovementTest_ul = nan(1, k);
			
			for il = 1:k
				percImprovementTest_ll(il) = 100*(nullLossTest(il) - testLoss_ll(il))/nullLossTest(il);
				percImprovementTest_ul(il) = 100*(nullLossTest(il) - testLoss_ul(il))/nullLossTest(il);
			end

			disp(['Next best lams=' mat2str([bestLam_ll, bestLam_ul])])
			disp(['Test Loss Improvement over Null: ' mat2str([round(meanPercImprovementTest_ll,8), round(meanPercImprovementTest_ul,8)]) ' % Improvement'])
			disp(['Lower Level Range Test Loss Improvement over Null Range: ' mat2str(round(sort(percImprovementTest_ll),8)) ' % Improvement'])
			disp(['Upper Level Range Test Loss Improvement over Null Range: ' mat2str(round(sort(percImprovementTest_ul),8)) ' % Improvement'])
		end
		function yFit = calcYfit(obj, th, X)
			yFit = th.'*X;
		end 
		function [se_model, se_th, CVmat, signifCoeff] = standardErrorOfModelAndTh(obj, XtX, th, yActual, yFit, lambda)
			se_model = sqrt(sum((yActual-yFit).^2./numel(yFit)));
			CVmat = (XtX+lambda*eye(size(XtX)))^-1*XtX*(XtX+lambda*eye(size(XtX)))^-1;
			se_th = se_model.*diag(CVmat).^.5; 
			distFromZero = abs(th) - 2*abs(se_th);
			signifCoeff = distFromZero > 0;
		end
		function E = MSELoss(obj,a,yFit)
			E = 1/numel(a)*sum((a - yFit).^2);
		end
		function [Resid, std_Resid, Rsq] = getModelResidualsAndR2(obj, yActual, yFit, th)
			Resid = yActual - yFit;
			std_Resid = sqrt(sum(Resid.^2)./(numel(yActual) - numel(th)));
			std_yActual = std(yActual);
% 			explainedVarianceR2 = 1 - std_Resid^2/std_yActual^2;

			% 
			%  Check consistent
			% 	
			ESS = sum((yFit - mean(yActual)).^2);
 			RSS = sum((yFit - yActual).^2);
 			Rsq = ESS/(RSS+ESS);
%  			assert(explainedVarianceR2 == Rsq)
		end

		function xValidate(obj, d, n, lam, Mode, k, recycle_k, modelName, suppressPlot, useSingle)
			disp('=========================================')
			disp('=			 crossvalidation	 		  =')
			disp('=========================================')
			disp(' ')
			% 
			% 	We will reuse the saved X and a vectors to run Xvalidations with an array of lambdas to try
			% 
			% 	The default will be to do leave-1-out Xvalidation, since we have the ability to do this!
			% 
			% obj.xValidate(d, n, lams, 'k-fold', k, false, true, true, useSingle);
			% 
			%---------------------------------
			if nargin < 9
				suppressPlot = false;
			end
			if nargin < 4
				lam = [];
			end
			if nargin < 2
				d = 1;
			end
			if nargin < 3
				n = 1;
			end
			if nargin < 7
				recycle_k = false;
			end
			
			
			if nargin < 5
				Mode = 'leave-1-out';
				nX = obj.collatedResults(d).kfoldSets.nTrials;
			elseif strcmpi(Mode, 'leave-1-out')
				nX = obj.collatedResults(d).kfoldSets.nTrials;
			elseif strcmpi(Mode, 'k-fold')
				% 
				% 	Once a k-fold set is drawn, use it for ALL nests. Stored in obj.kfoldSets
				% 
				if ~recycle_k
					if nargin < 6
						k = 5;
					end
					nX = k;
					obj.getKfoldSet(d, k)
					
					trialSets = obj.collatedResults(d).kfoldSets.trialSets;

					obj.collatedResults(d).XvalResults(n).k = k;
					obj.collatedResults(d).XvalResults(n).trialSets = trialSets;
				else
					nX = obj.collatedResults(d).kfoldSets.k;
					trialSets = obj.collatedResults(d).kfoldSets.trialSets;
				end
			elseif strcmpi(Mode, 'check-lam')
				warning('We are checking out lam effect before running the cross-validation.')
				nX = 1;
				trial_check = randi(obj.collatedResults(d).kfoldSets.nTrials);
			else
				warning('Running in debug mode - will stop after 5x leave-1-out cross validations')
			end
			disp(['Dataset: ' obj.collatedResults(d).sessionID ' Nest #' num2str(n)])
			disp(['Calculating ' num2str(nX) 'x XVs for lamda=' mat2str(lam)])
			disp([' '])
			% 
			% 	We will have a set of models for each lamda
			%
			if strcmpi(Mode, 'k-fold')
				% 
				% 	Check in lam already in our set - if so, ignore it
				% 
				if ~isempty(obj.collatedResults(d).XvalResults) && isfield(obj.collatedResults(d).XvalResults, 'n') && numel(obj.collatedResults(d).XvalResults) >= n && isfield(obj.collatedResults(d).XvalResults, 'models') && ~isempty(obj.collatedResults(d).XvalResults(n).models)
					disp(['	Detected existing models for ' obj.collatedResults(d).sessionID ' nest#' num2str(n)])
					igLamIdxs = [];
					for ilam = 1:numel(lam)
						if sum(ismember([obj.collatedResults(d).XvalResults(n).models.lam], lam(ilam)))
							igLamIdxs(end+1) = ilam;
						end
					end
					lam(igLamIdxs) = [];
                end
            end
            if isempty(lam)
                disp('All input lams were already tested and in the results set')
            else
				for l = 1:numel(lam)
					models(l).lam = lam(l);
					if strcmpi(Mode, 'k-fold')
						models(l).kfoldID = obj.collatedResults(d).kfoldSets.ID;
						models(l).k = obj.collatedResults(d).kfoldSets.k;
					end
				end
				% 
				% 	models have stats as XV x lam
				% 			lam1	lam2	lam3...
				% 	XV 1
				% 	XV 2
				% 	...
				
				[models.nullLossTraining] = deal(nan(nX, 1));
				[models.meanTrainingLoss] = deal(nan(nX, 1));
				[models.nullLossTest] = deal(nan(nX, 1));
				[models.meanTestLoss] = deal(nan);
				[models.trainingLoss] = deal(nan(nX, 1));
				[models.testLoss] = deal(nan(nX, 1));
				[models.th] = deal(cell(nX, 1));
				[models.se_model] = deal(nan(nX, 1));
				[models.se_th] = deal(cell(nX, 1));
				[models.R2] = deal(nan(nX, 1));
				%
				%	Extract the relevant X and a for Xvalidation 
				% 	
				eval(['X = obj.collatedResults(d).decoding.' modelName '.X(1:n,:);'])
				eval(['nests = obj.collatedResults(d).decoding.' modelName '.predictorSet(1:n);'])
				eval(['a = obj.collatedResults(d).decoding.' modelName '.y{end};']) % using the final nest because this is fairest -- the dataset shrinks when we add the final predictors, so fair is fair here...
				
				dd = size(X, 1);
				% 
				% 	Now, for each XV, let's get the Xxv and aXV and calculate
				% 	the loss for each lam
				% 
				if strcmpi(Mode, 'check-lam') || ~suppressPlot
					figure,
					ax = subplot(1,2,1);
					ax2 = subplot(1,2,2);
					hold(ax, 'on');
					hold(ax2, 'on');
					C = linspecer(numel(lam));
				end
				disp(['Initiating cross-validation fitting ' datestr(now)])
				for xv = 1:nX
					% if ~rem(xv, 25)
					% 	disp(['On cross-validation ' num2str(xv) ' of ' num2str(nX) ' ' datestr(now)])
					% end
					obj.progressBar(xv, nX, false, 5)
					% 
					% 	Get the X and a
					% 
					if strcmpi(Mode, 'leave-1-out')
						trial = xv;
						a_xv = a;
						X_xv = X;
						a_test = a;
						X_test = X;	
					elseif strcmpi(Mode, 'check-lam')
						trial = trial_check;
						a_xv = a;
						X_xv = X;
						a_test = a;
						X_test = X;
					elseif strcmpi(Mode, 'k-fold')
						trials_test = trialSets{xv};	
						trials_training = trialSets;
						trials_training{xv} = [];
						trials_training = cell2mat(trials_training);
						
	                    Sn_killpos = trials_test;
	                    a_xv = a;
						X_xv = X;
	                    a_xv(Sn_killpos) = [];
	                    X_xv(:,Sn_killpos) = [];
	                    S_killpos = trials_training;
						a_test = a;
						X_test = X;
	                    a_test(S_killpos) = [];
	                    X_test(:, S_killpos) = [];
					end
					
					
					if strcmpi(Mode, 'check-lam') && ~suppressPlot
						xv = 1;
						plot(ax, a_test, 'DisplayName', 'dF/F')
						title('Lambda Check')
						xlabel('samples')
						ylabel('model dF/F')
					end


					for ilam = 1:numel(lam)
						XtX = X_xv*X_xv.';
						th = (XtX+lam(ilam).*eye(dd))\X_xv*a_xv.';
						models(ilam).th = th;
						% 
						% 	Return the yFit
						% 
						yFit = obj.calcYfit(th, X_xv);
						[models(ilam).se_model(xv, 1), models(ilam).se_th{xv, 1}, ~, ~] = obj.standardErrorOfModelAndTh(XtX, th, a_xv, yFit, lam(ilam));
						[~, ~, models(ilam).R2(xv, 1)] = obj.getModelResidualsAndR2(a_xv, yFit, th);
						models(ilam).trainingLoss(xv, 1) = obj.MSELoss(a_xv,yFit);
						models(ilam).nullLossTraining(xv, 1) = obj.MSELoss(a_xv,mean(a_xv).*ones(size(a_xv)));
						models(ilam).nullLossTest(xv, 1) = obj.MSELoss(a_test,mean(a_test).*ones(size(a_test)));
						models(ilam).testLoss(xv, 1) = obj.MSELoss(a_test,obj.calcYfit(th, X_test));

						
						if strcmpi(Mode, 'check-lam') && ~suppressPlot
							disp(['	lam: ' num2str(lam(ilam)) ' | TestLoss: ' num2str(models(ilam).testLoss(xv, 1))])
							if lam(ilam) == 0
								plot(ax, obj.calcYfit(th, X_test), 'k-', 'linewidth',3, 'DisplayName', num2str(lam(ilam)))
							else							
								plot(ax, obj.calcYfit(th, X_test), 'color', C(ilam, :), 'DisplayName', num2str(lam(ilam)))
							end
						end
					end
				end
				for ilam = 1:numel(lam)
					models(ilam).meanTestLoss = mean([models(ilam).testLoss(:, 1)]);
					models(ilam).meanTrainingLoss = mean([models(ilam).trainingLoss(:, 1)]);
				end

				% 
				% save everything to the xval structure. If we have a matching set of k-codes, then we should append, else we will replace
				% 
                obj.collatedResults(d).XvalResults(n).models = models;
                obj.collatedResults(d).XvalResults(n).nests = nests;
% 				if strcmpi(Mode, 'k-fold')
% 					for iModel = 1:numel(lam)
% 						if ~sum(ismember([obj.collatedResults(d).XvalResults(n).models.lam], lam(iModel)))
% 							obj.collatedResults(d).XvalResults(n).models = [obj.collatedResults(d).XvalResults(n).models, models(iModel)];
% 						end
% 					end
% 				else
% 					error('rbf, not tested')
% 					obj.collatedResults(d).XvalResults(n).models = models;
% 				end
				obj.analysis.flush.n = n;
	            obj.analysis.flush.d = d;
	            obj.analysis.flush.X = X;
				obj.analysis.flush.a = a;
			end
			
			if strcmpi(Mode, 'check-lam') && ~suppressPlot
				legend(ax, 'show');
				plot(ax2, [models.nullLossTraining], 'b--', 'DisplayName', 'Null Loss - Training')
				plot(ax2, [models.nullLossTest], 'r--', 'DisplayName', 'Null Loss - Test')

				plot(ax2, [models.meanTrainingLoss], 'b-o', 'DisplayName', 'Training Loss')
				plot(ax2, [models.meanTestLoss], 'r-o', 'DisplayName', 'Test Loss')
				xticks(ax2, 1:length(obj.collatedResults(d).XvalResults(n).models));
				xticklabels(ax2, lamLabel);
				title(ax2, 'Loss vs Regularization')
				xlabel(ax2, 'lambda')
				ylabel(ax2, 'loss')
				legend(ax2, 'show')
				xtickangle(ax2,90);
			elseif ~suppressPlot
				obj.xValLossPlotVsLam(d, n);
			end

			% 
			% 	Update selectBestLam(obj, d, n)
			% 
			obj.selectBestLam(d, n)
            disp('Complete. ~~')
		end
		

		function mouseLevelXval(obj, d, n, modelNum, k, nLamsPerRound, useSingle)
            try
                disp('==================================================================================')
                disp('==		Running AutoXval for multiple datasets...							  ==') 
                disp('==================================================================================') 
                jobID = randi(10000000);
                % 
                % 	Autoruns autoXval(obj, d, n, nRounds, k, nLamsPerRound) for all datasets, d
                % 
                % 	d = vector of dataset idx
                % 	n = vector of nest idx or 'all'
                % 	If only one type of model saved in the set, then just leave modelNum empty, or set to 1
                % 
                if nargin < 7
                	useSingle = false;
            	elseif useSingle
            		warning('Using single point precision.')
            	end
                if nargin < 2 || isempty(d)
                    d = 1:numel(obj.collatedResults);
                end
                if nargin < 3 || isempty(n)
                    n = 'all';
                end
                if nargin < 4 || isempty(modelNum)
                    modelNum = 1;
                end
                modelName = fieldnames(obj.collatedResults(1).decoding);
                modelName = modelName{modelNum};
                if nargin < 5 || isempty(k)
                    k = 20;
                end
                if nargin < 6 || isempty(nLamsPerRound)
                    nLamsPerRound = 20;
                end
                for id = d
	                seshCode = obj.collatedResults(id).sessionID;
                	roundID = ['JobID: ' num2str(jobID), ' | Current Signal ' seshCode ' | d=' num2str(id)];
                	obj.collatedResults(id).kfoldSets.modelName = modelName;
                    disp('==================================================================================') 
                    disp(['	Initiating for dataset ' num2str(find(d == id)) ' of ' num2str(numel(d))])
                    disp('    ')
                    obj.progressBar(find(d == id), numel(d), false, 1)
                    if strcmpi(n, 'all')
                        eval(['n_d = 1:numel(obj.collatedResults(id).decoding.' modelName '.predictorSet);'])
                    else
                        n_d = n;
                    end
                    obj.collatedResults(id).jobID = jobID;
                    for in_d = n_d
                        disp(['	Working on dataset ' num2str(find(d == id)) ' of ' num2str(numel(d)) ' | Nest #'])
                        obj.progressBar(find(n_d == in_d), numel(n_d), true, 1)
                        obj.autoXval(id, in_d, modelName, k, nLamsPerRound, useSingle);
                        
                        obj.analysis.flush = []; % free up memory immediately before saving
                        % obj.save;
                        % mailAlert(['mouseLevelXval Decoding Job' num2str(jobID) ' in Progress. Now complete: ' seshCode ' n=' num2str(in_d) '/' num2str(numel(n_d)) ' d=' num2str(find(d == id)) '/' num2str(numel(d))], roundID);
                    end
                    % 
                    % 	Remove spent Stats fields to save space
                    % 
                end
            catch EX
                msg = ['Exception Thrown: ' EX.identifier ' | ' EX.message '\n\n' roundID];
                alert = ['ERROR in mouseLevelXval Decoding Job in Progress.']; 
                mailAlert(alert, msg);
                rethrow(EX)
            end
			disp('==================================================================================') 
			disp('Complete.')
			mailAlert(['mouseLevelXval Decoding Job' num2str(jobID) ' COMPLETE without errors!']);
		end
		function refitXval(obj, dd, nn, lam)
            
			if nargin < 2
				dd = 1:numel(obj.collatedResults);
			end
			if isfield(obj.collatedResults(1).kfoldSets,'modelName')
				modelName = obj.collatedResults(1).kfoldSets.modelName;
			else
				modelName = fieldnames(obj.collatedResults(1).decoding);
				modelName = modelName{1};
			end
			if nargin < 3
				eval(['nn = 1:numel(obj.collatedResults(dd(1)).decoding.' modelName '.n);'])
			end
			disp('==================================================================================')
            disp('==		Refitting Xval fits...												  ==') 
            disp('==================================================================================') 
			try
				for d = dd		
					disp(['Dataset:' obj.collatedResults(d).sessionID])
                    obj.collatedResults(d).refit = [];
					obj.progressBar(find(dd==d), numel(dd), false, 1)			
                    

					for n = nn
% 						obj.progressBar(find(nn==n), numel(nn), true, 1)
						if nargin < 4 || isempty(kfoldID)
							kfoldID = obj.collatedResults(d).kfoldSets.ID;
							k = obj.collatedResults(d).kfoldSets.k;
						end
						if nargin < 5
							lam = obj.collatedResults(d).XvalResults(n).bestLam;
% 							disp([num2str(n) '		Using best lam=' num2str(lam)])
						end						
						% 
						% 	Collect X and a from file
						% 
						
						eval(['X = obj.collatedResults(d).decoding.' modelName '.X(1:n,:);'])
						eval(['a = obj.collatedResults(d).decoding.' modelName '.y{end};'])
						eval(['dfe = obj.collatedResults(d).decoding.' modelName '.stats{1, n}.dfe;'])
						meanAloss = []; 
						modelSquaredLoss = [];
						th = []; 
						se_model = []; 
						se_th = [];
						signifCoeff= [];
						Resid = [];
						std_Resid = [];
						explainedVarianceR2 = [];
						AIC = [];
						AICc = [];
						nAIC = [];
						BIC = [];
						% 
						% 	Now refit the model with lam
						% 
						% disp(['XX.T is pseudo-invertible! Using analytical RIDGE solution, lambda = ', num2str(lam)])
						% 
						% if matrix is invertible, we will use analytical solution
						% 
						% 
						% 	Combine singular variables!
						% 
						singularIndicies = [];
						for ir = 2:size(X,1)
							if sum(X(ir,:)) == size(X,2)
								singularIndicies(end+1) = ir;
							end
                        end
                        Xfix = X;
                        if ~isempty(singularIndicies)
                            Xfix(singularIndicies,:) = []; 
                        end
						% 
						% 
						% 
						XtX = Xfix*Xfix.';
						th = (XtX+lam.*eye(size(XtX, 1)))\Xfix*a.';
						% 
						% 	Return the yFit
						% 
						yFit = th.'*Xfix;
						[se_model, se_th, CVmat, signifCoeff] = obj.standardErrorOfModelAndTh(XtX, th, a, yFit, lam);
						[Resid, std_Resid, explainedVarianceR2] = obj.getModelResidualsAndR2(a, yFit, th);	
                        [AIC, AICc, nAIC, BIC] = testAIC(a, th, yFit);
                        meanAloss = 1/numel(a)*sum((a - mean(a)).^2);
						modelSquaredLoss = 1/numel(a)*sum((a - yFit).^2);	
						% 
						% 	Add back the singular components
						% 
                        if ~isempty(singularIndicies)
                            if n == max(nn)
                                disp(['   Singular correction on nest: ' mat2str(singularIndicies)])
                            end
                            th = [th(1:singularIndicies(1)-1); nan(numel(singularIndicies),1); th(singularIndicies(1):end)];
                            se_th = [se_th(1:singularIndicies(1)-1); nan(numel(singularIndicies),1); se_th(singularIndicies(1):end)];
                        end
                            

						CImin = th - abs(tinv(.025,dfe)).*se_th;
						CImax = th + abs(tinv(.025,dfe)).*se_th;
						close


						% 
						% 	Compile results
						%
						eval(['obj.collatedResults(d).refit.predictorKey = obj.collatedResults(d).decoding.' modelName '.predictorKey;'])
						eval(['obj.collatedResults(d).refit.trials_in_range_first = obj.collatedResults(d).decoding.' modelName '.trials_in_range_first;'])
						eval(['obj.collatedResults(d).refit.Conditioning = obj.collatedResults(d).decoding.' modelName '.Conditioning;'])
						eval(['obj.collatedResults(d).refit.Link = obj.collatedResults(d).decoding.' modelName '.Link;'])
						eval(['obj.collatedResults(d).refit.predictorNames = obj.collatedResults(d).decoding.' modelName '.predictorNames;'])
						eval(['obj.collatedResults(d).refit.n = obj.collatedResults(d).decoding.' modelName '.n;'])
						obj.collatedResults(d).refit.b{n,1} = th;
						obj.collatedResults(d).refit.CImin{n,1} = CImin;
						obj.collatedResults(d).refit.CImax{n,1} = CImax;
						obj.collatedResults(d).refit.stats{n,1}.dfe = dfe;
						obj.collatedResults(d).refit.stats{n,1}.lam = lam;
						obj.collatedResults(d).refit.stats{n,1}.beta = th;
						obj.collatedResults(d).refit.stats{n,1}.s = se_model;
						obj.collatedResults(d).refit.stats{n,1}.se = se_th;
						obj.collatedResults(d).refit.stats{n,1}.covb = CVmat;
						obj.collatedResults(d).refit.y{n,1} = a;
						obj.collatedResults(d).refit.yfit{n,1} = yFit;
						obj.collatedResults(d).refit.LossImprovement(n,1) = modelSquaredLoss/meanAloss;
						obj.collatedResults(d).refit.BIC(n,1) = BIC;
						obj.collatedResults(d).refit.AIC(n,1) = AIC;
						obj.collatedResults(d).refit.AICc(n,1) = AICc;
						obj.collatedResults(d).refit.nAIC(n,1) = nAIC;
						obj.collatedResults(d).refit.Rsq(n,1) = explainedVarianceR2;
						eval(['obj.collatedResults(d).refit.ESS(n,1) = obj.collatedResults(d).decoding.' modelName '.ESS(n);'])
						eval(['obj.collatedResults(d).refit.RSS(n,1) = obj.collatedResults(d).decoding.' modelName '.ESS(n);'])
						eval(['obj.collatedResults(d).refit.predictorSet = obj.collatedResults(d).decoding.' modelName '.predictorSet;'])
						obj.collatedResults(d).refit.X = X;
						obj.collatedResults(d).refit.Xn{n,1} = X;
                        obj.collatedResults(d).refit.singularIndicies{n,1} = singularIndicies;
						obj.collatedResults(d).refit.kfoldID(n,1) = kfoldID;
						obj.collatedResults(d).refit.k(n,1) = k;						 
					end
				end
				% 
				%	Create obj with all the results together	 
				% 
				alert = ['ALL DONE: Decoding Refit XvalResults COMPLETE. d=' num2str(max(d)) ' n=' mat2str(max(nn))]; 
			    mailAlert(alert);
		    catch EX
                alert = ['ERROR in refitXval']; 
                msg = ['Exception Thrown: ' EX.identifier ' | ' EX.message ];
                mailAlert(alert, msg);
                rethrow(EX)
            end
			disp('==================================================================================') 
			disp('Complete. Ready to get stats on refit data')
		end





		function conditionPCA(obj, minPCs, killSessions)
			% 
			% 	Assign good indicies based on these criteria (default)
			% 
			% 	minPCs = 30
			% 
			if nargin < 3
				killSessions = false;
			end
			if nargin < 2
				minPCs = 30;
			end
			sessionsToKill = {'DLSred',...
			'H6_SNc_d23_1',...
			'H6_SNc_d23_2',...
			'H6_SNc_d23_3',...
			'H6_SNc_d22_1',...
			'H6_SNc_d22_2',...
			'H6_SNc_d22_3',...
			'H6_SNc_d22_4',...
			'H6_SNc_d22_5',...
			'H6_SNc_d21',...
			'H6_SNc_d20',...
			'B6_SNc_12',...
			'B6_SNc_13',...
			'B6_SNc_14_1',...
			'B6_SNc_15_1',...
			'B6_SNc_15_2',...
			'B3_SNc_21',...
			'B3_SNc_20',...
			'B5_SNc_19',...
			'B5_SNc_20_1',...
			'B5_SNc_20_2',...
			'B5_SNc_20_3',...
			'B5_SNc_20_4',...
			'B5_SNc_21_1',...
			'B5_SNc_21_2',...
			'B5_SNc_17',...
			'B5_SNc_19',...
			'B6_SNc_13',...
			'H3_SNc_17',...
			'H3_SNc_18',...
			'H3_SNc_19',...
			'H3_SNc_20',...
			'H5_SNc_11',...
			'H5_SNc_17',...
			'H5_SNc_18',...
			'H5_SNc_19',...
			'H7_SNc_12',...
			'H7_SNc_13',...
			'H7_SNc_14_2',...
			'H7_SNc_15_1',...
			'H7_SNc_15_2',...
			'H7_SNc_15_3',...
			'H7_SNc_15_4',...
			'H7_SNc_15_5',...
			'H7_SNc_16_1',...
			'H7_SNc_16_2',...
			'H7_SNc_16_3',...
			};
			for d = 1:numel(obj.collatedResults)
				if numel(obj.collatedResults(d).PCA.mu) < minPCs
					obj.collatedResults(d).QCflag = true;
				elseif killSessions && sum(contains(obj.collatedResults(d).sessionID,sessionsToKill))>0
					obj.collatedResults(d).QCflag = true;
				else
					obj.collatedResults(d).QCflag = false;
				end
			end
		end
		function plotPCA(obj,killSessions,Mode,idxs)
			obj.conditionPCA(30, killSessions);
        	% 
        	% 	Mode: 'summary': plots a summary of explained variance of PCs and top 3 PCs
        	% 		  
    		excludetime = obj.collatedResults(1).PCA.excludetime;
    		maxWindow_trimRHS = obj.collatedResults(1).PCA.maxWindow_trimRHS;
            enforcedlatency = obj.collatedResults(1).PCA.enforcedlatency;

            if nargin < 4
            	idxs = find(~[obj.collatedResults.QCflag]);
        	end
    		if nargin < 3
    			Mode = 'summary';
			end

			f = figure;
			if strcmpi(Mode, 'summary')
				ax_rsq = subplot(1,2,1);
				hold(ax_rsq, 'on')
				ylim(ax_rsq, [0, 100])
				ylabel(ax_rsq, '% Variance Explained')
				title(ax_rsq, ['n = ' num2str(numel(idxs)) '/' num2str(numel(obj.collatedResults))])
				xlabel(ax_rsq, 'PC#')
				xlim(ax_rsq, [1,10])
				ax_pc1 = subplot(3,2,2);
				ax_pc2 = subplot(3,2,4);
				ax_pc3 = subplot(3,2,6);
				title(ax_pc1,'PC #1')
				title(ax_pc2,'PC #2')
				title(ax_pc3,'PC #3')
				xlim(ax_pc1, [0,maxWindow_trimRHS/1000])
				xlim(ax_pc2, [0,maxWindow_trimRHS/1000])
				xlim(ax_pc3, [0,maxWindow_trimRHS/1000])
				hold(ax_pc1, 'on')
				hold(ax_pc2, 'on')
				hold(ax_pc3, 'on')

				plot(ax_pc1, [0,7], [0,0], 'k-', 'HandleVisibility', 'off', 'linewidth',3)
				plot(ax_pc2, [0,7], [0,0], 'k-','HandleVisibility', 'off', 'linewidth',3)
				plot(ax_pc3, [0,7], [0,0], 'k-','HandleVisibility', 'off', 'linewidth',3)
				
				for dd = 1:numel(idxs)
                    d = idxs(dd);
					rsq{dd,1} = obj.collatedResults(d).PCA.explained;
					plot(ax_rsq, rsq{dd}, '-', 'color', [0.2,0.2,0.2], 'DisplayName', obj.collatedResults(d).sessionID)
					s1{dd,1} = obj.collatedResults(d).PCA.score(:,1);
					s2{dd,1} = obj.collatedResults(d).PCA.score(:,2);
					s3{dd,1} = obj.collatedResults(d).PCA.score(:,3);
					if mean(s1{dd}) < 0
						s1{dd} = -1.*s1{dd};
					end
					if mean(s2{dd}(end-100:end)) < 0
						s2{dd} = -1.*s2{dd};
					end
					if mean(s3{dd}(end-100:end)) < 0
						s3{dd} = -1.*s3{dd};
					end
					plot(ax_pc1, linspace(excludetime,maxWindow_trimRHS/1000,maxWindow_trimRHS-excludetime*1000), s1{dd}, '-', 'color', [0.2,0.2,0.2],'HandleVisibility', 'off')
					plot(ax_pc2, linspace(excludetime,maxWindow_trimRHS/1000,maxWindow_trimRHS-excludetime*1000), s2{dd}, '-', 'color', [0.2,0.2,0.2],'HandleVisibility', 'off')
					plot(ax_pc3, linspace(excludetime,maxWindow_trimRHS/1000,maxWindow_trimRHS-excludetime*1000), s3{dd}, '-', 'color', [0.2,0.2,0.2],'HandleVisibility', 'off')
				end
				rsq10 = cellfun(@(x) x(1:10), rsq,'uniformoutput', false);
				mrsq = mean(cell2mat(rsq10'), 2);
				ms1 = mean(cell2mat(s1'),2);
				ms2 = mean(cell2mat(s2'),2);
				ms3 = mean(cell2mat(s3'),2);

				plot(ax_rsq, mrsq, 'r-o', 'linewidth', 4)
				plot(ax_pc1, linspace(excludetime,maxWindow_trimRHS/1000,maxWindow_trimRHS-excludetime*1000), ms1, 'r-', 'linewidth', 4)
				plot(ax_pc2, linspace(excludetime,maxWindow_trimRHS/1000,maxWindow_trimRHS-excludetime*1000), ms2, 'r-', 'linewidth', 4)
				plot(ax_pc3, linspace(excludetime,maxWindow_trimRHS/1000,maxWindow_trimRHS-excludetime*1000), ms3, 'r-', 'linewidth', 4)
				x1 = get(ax_pc1, 'ylim');
				x2 = get(ax_pc2, 'ylim');
				x3 = get(ax_pc3, 'ylim');
				xl = min([x1(1), x2(1), x3(1)]);
				xu = max([x1(2), x2(2), x3(2)]);
				ylim(ax_pc1,[xl, xu])
				ylim(ax_pc2,[xl, xu])
				ylim(ax_pc3,[xl, xu])
				linkaxes([ax_pc1,ax_pc2,ax_pc3],'xy')
			end
		end





		function gatherDivergenceIndicies(obj)
			nEE = 0;
			nER = 0;
			nRE = 0;
			nRR = 0;
			% obj.analysis.divergenceIndex.meanEX = {};
			% obj.analysis.divergenceIndex.EX_CIlow = {};
			% obj.analysis.divergenceIndex.EX_CIhi = {};

			% obj.analysis.divergenceIndex.meanRX = {};
			% obj.analysis.divergenceIndex.RX_CIlow = {};
			% obj.analysis.divergenceIndex.RX_CIhi = {};

			% obj.analysis.convergenceIndex.meanXE = {};
			% obj.analysis.convergenceIndex.XE_CIlow = {};
			% obj.analysis.convergenceIndex.XE_CIhi = {};

			% obj.analysis.convergenceIndex.meanXR = {};
			% obj.analysis.convergenceIndex.XR_CIlow = {};
			% obj.analysis.convergenceIndex.XR_CIhi = {};

			for d = 1:numel(obj.collatedResults)
				nEE = nEE + obj.collatedResults(d).Stat.n.nEE;
				nER = nER + obj.collatedResults(d).Stat.n.nER;
				nRE = nRE + obj.collatedResults(d).Stat.n.nRE;
				nRR = nRR + obj.collatedResults(d).Stat.n.nRR;
				if d == 1
					obj.analysis.composite.meanEX = obj.collatedResults(d).Stat.divergenceIndex.meanEX;
					obj.analysis.composite.meanRX = obj.collatedResults(d).Stat.divergenceIndex.meanRX;
					obj.analysis.composite.meanXE = obj.collatedResults(d).Stat.convergenceIndex.meanXE;
					obj.analysis.composite.meanXR = obj.collatedResults(d).Stat.convergenceIndex.meanXR;
				else
					obj.analysis.composite.meanEX = nansum([obj.analysis.composite.meanEX .* ((d-1)./d); obj.collatedResults(d).Stat.divergenceIndex.meanEX]);
					obj.analysis.composite.meanRX = nansum([obj.analysis.composite.meanRX .* ((d-1)./d); obj.collatedResults(d).Stat.divergenceIndex.meanRX]);
					obj.analysis.composite.meanXE = nansum([obj.analysis.composite.meanXE .* ((d-1)./d); obj.collatedResults(d).Stat.convergenceIndex.meanXE]);
					obj.analysis.composite.meanXR = nansum([obj.analysis.composite.meanXR .* ((d-1)./d); obj.collatedResults(d).Stat.convergenceIndex.meanXR]);
				end
				% obj.analysis.divergenceIndex(d).meanEX = obj.collatedResults(d).Stat.divergenceIndex.meanEX;
				% obj.analysis.divergenceIndex(d).meanRX = obj.collatedResults(d).Stat.divergenceIndex.meanRX;
				% obj.analysis.convergenceIndex(d).meanXE = obj.collatedResults(d).Stat.divergenceIndex.meanXE;
				% obj.analysis.convergenceIndex(d).meanXR = obj.collatedResults(d).Stat.divergenceIndex.meanXR;
				obj.analysis.n.nEE = nEE;
				obj.analysis.n.nER = nER;
				obj.analysis.n.nRE = nRE;
				obj.analysis.n.nRR = nRR;
			end


		end
		function plotDivergenceIndicies(obj, Mode, suppressNsave)
			% 
			% 	Mode: 'composite', 'overlay'
			% 
			obj.gatherDivergenceIndicies
			if nargin < 2
				Mode = 'composite';
			end
			if nargin < 3
				suppressNsave = [];
			end
			
			c = obj.collatedResults(1).Stat.centers/1000;

			[f, ax] = makeStandardFigure(4,[2,2]);
			hold(ax(1), 'on');
            hold(ax(2), 'on');
            hold(ax(3), 'on');
            hold(ax(4), 'on');
            plot(ax(1), [c(1), c(end)], [0,0], 'k-', 'LineWidth',1)
            plot(ax(2), [c(1), c(end)], [0,0], 'k-', 'LineWidth',1)
            plot(ax(3), [c(1), c(end)], [0,0], 'k-', 'LineWidth',1)
            plot(ax(4), [c(1), c(end)], [0,0], 'k-', 'LineWidth',1)
            title(ax(1), ['diverg: EE(' num2str(obj.analysis.n.nEE) ')-ER(' num2str(obj.analysis.n.nER) ')'])
            title(ax(3), ['convergence: RE-EE'])
            
            ylabel(ax(1),'nth trial selectivity')
            title(ax(2), ['divergence: RE(' num2str(obj.analysis.n.nRE) ')-RR(' num2str(obj.analysis.n.nRR) ')'])
            ylabel(ax(3),'n-1th trial selectivity')
            title(ax(4), ['convergence: RR-ER'])
            xlim(ax(1), [-10,3])
            xlim(ax(2), [-10,3])
            xlim(ax(3), [-10,3])
            xlim(ax(4), [-10,3])

            if strcmpi(Mode, 'composite')
            	plot(ax(1), c, obj.analysis.composite.meanEX, 'r-', 'LineWidth',3)				
				plot(ax(2), c, obj.analysis.composite.meanRX, 'r-', 'LineWidth',3)
				plot(ax(3), c, -1.*obj.analysis.composite.meanXE, 'r-', 'LineWidth',3)
				plot(ax(4), c, -1.*obj.analysis.composite.meanXR, 'r-', 'LineWidth',3)
        	elseif strcmpi(Mode, 'overlay')
        		for d = 1:numel(obj.collatedResults)
	        		plot(ax(1), c, obj.collatedResults(d).Stat.divergenceIndex.meanEX, 'r-', 'LineWidth',3)
					% plot(ax(1), c, obj.Stat.divergenceIndex.EX_CIlow, 'k-', 'LineWidth',2)
					% plot(ax(1), c, obj.Stat.divergenceIndex.EX_CIhi, 'k-', 'LineWidth',2)
					plot(ax(2), c, obj.collatedResults(d).Stat.divergenceIndex.meanRX, 'r-', 'LineWidth',3)
					% plot(ax(2), c, obj.Stat.divergenceIndex.RX_CIlow, 'k-', 'LineWidth',2)
					% plot(ax(2), c, obj.Stat.divergenceIndex.RX_CIhi, 'k-', 'LineWidth',2)					
					plot(ax(3), c, -1.*obj.collatedResults(d).Stat.convergenceIndex.meanXE, 'r-', 'LineWidth',3)
					% plot(ax(3), c, -1.*obj.Stat.convergenceIndex.XE_CIlow, 'k-', 'LineWidth',2)
					% plot(ax(3), c, -1.*obj.Stat.convergenceIndex.XE_CIhi, 'k-', 'LineWidth',2)
					plot(ax(4), c, -1.*obj.collatedResults(d).Stat.convergenceIndex.meanXR, 'r-', 'LineWidth',3)
					% plot(ax(4), c, -1.*obj.Stat.convergenceIndex.XR_CIlow, 'k-', 'LineWidth',2)
					% plot(ax(4), c, -1.*obj.Stat.convergenceIndex.XR_CIhi, 'k-', 'LineWidth',2)
				end
    		else 
    			error('undefined Mode')
    		end        		
			if ~isempty(suppressNsave), figureName = ['DivergenceIndex_' obj.Stat.outcomeMode.Mode];, obj.suppressNsaveFigure(suppressNsave, figureName, f), close(f), end
		end
		function compositeANOVAdataset(obj, idx)
			% 
			% 	for baselineANOVAwithLick or baselineANOVAidx objects, called by A2ElevelANOVA(obj,nLevels,idx, overwrite)
			% 	NB: switch C and D levels for 3-way on baselineANOVAidx objs with no lick stuff
			% 
			% A_level = n-1th outcome 
			% B_level = nth outcome
			% C_level = # of licks in baseline
			% D_level = session ID
			% E_level = mouse ID (redundant with session so not really used)
			% 	
			% 
			if nargin < 2
				idx = 1:numel(obj.collatedResults);
			end
			% 
			% 	Gather 1 result slot by collating all the A_level, B_level and data across the sessionIdxs
			%
			if ~isfield(obj.analysis,'mouseName')
				for ii = 1:numel(obj.collatedResults) 
					obj.analysis.mouseName{ii} = strsplit(obj.collatedResults(ii).sessionID,'_');
					obj.analysis.mouseName{ii} = obj.analysis.mouseName{ii}{1};
				end
			end
			mice = unique(obj.analysis.mouseName);
			% 
			% 	Start by gathering the datasets
			% 
			% for each position in baseline gather across sessions and also track which animal and session
			data = cell(1,numel(obj.collatedResults(1).results));
			A_level = cell(1,numel(obj.collatedResults(1).results));
			B_level = cell(1,numel(obj.collatedResults(1).results));
			C_level = cell(1,numel(obj.collatedResults(1).results));
			D_level = cell(1,numel(obj.collatedResults(1).results));
			E_level = cell(1,numel(obj.collatedResults(1).results));
			for M = 1:numel(mice) 
				midx = find(contains(obj.analysis.mouseName,mice{M}));
				for dd = 1:numel(midx)
					d = midx(dd);
					if ismember(d, idx)
						for pos = 1:numel(obj.collatedResults(d).results)
							data{1,pos} = [data{1,pos};obj.collatedResults(d).results{1, pos}.cellData.data];  	% df/f
							A_level{1,pos} = [A_level{1,pos};obj.collatedResults(d).results{1, pos}.cellData.A_level'];   % n-1 category
							B_level{1,pos} = [B_level{1,pos};obj.collatedResults(d).results{1, pos}.cellData.B_level'];  % nth category
							% C_level{1,pos} = [C_level{1,pos};obj.collatedResults(d).results{1, pos}.cellData.C_level']; % not imp yet but could be EMG spikes
							% D_level{1,pos} = [D_level{1,pos};obj.collatedResults(d).results{1, pos}.cellData.D_level']; % session number, to be assigned
							E_level{1,pos} = [E_level{1,pos};obj.collatedResults(d).results{1, pos}.cellData.lickLevel']; % mouse, to be assigned
							MM = cell(size(obj.collatedResults(d).results{1, pos}.cellData.data));
							MM(:) = mice(M);
							D_level{1,pos} = [D_level{1,pos}; MM];
							SS = cell(size(obj.collatedResults(d).results{1, pos}.cellData.data));
							SS(:) = {obj.collatedResults(d).sessionID};
							C_level{1,pos} = [C_level{1,pos}; SS];
						end
					end
				end
			end
			obj.analysis.A_level = A_level;
			obj.analysis.B_level = B_level;
			obj.analysis.C_level = C_level;
			obj.analysis.D_level = D_level;
			obj.analysis.E_level = E_level;
			obj.analysis.data = data;
			obj.analysis.nMice = numel(unique(D_level{1}));
			obj.analysis.nSesh = numel(unique(C_level{1}));
		end
		function A2ElevelANOVA(obj,nLevels,idx, overwrite)
			% 
			% 	for baselineANOVAwithLick or baselineANOVAidx objects
			% 
			% 	To just run 17 MDT seshs, use: A2ElevelANOVA(obj,4, [27:29,38,39,40,47,48,97:101,117,118,131,132],true) 
			% 
            if nargin < 4
                overwrite = false;
            end
			if nargin < 3 || isempty(idx)
				idx = 1:numel(obj.collatedResults);
			elseif strcmpi(idx, 'MDT')
				obj.flagMDT;
				idx = find([obj.collatedResults.MDT]);
				overwrite = true;
			end
			if nargin < 2
				nLevels = 4;
			end
			if ~isfield(obj.analysis, 'A_level') || overwrite
				compositeANOVAdataset(obj, idx);
			end
			A_level = obj.analysis.A_level;
			B_level = obj.analysis.B_level;
			C_level = obj.analysis.C_level;
			D_level = obj.analysis.D_level;
			E_level = obj.analysis.E_level;
			data = obj.analysis.data;
			for pos = 1:numel(obj.collatedResults(1).results)
				if nLevels == 2
					levels = {A_level{pos},B_level{pos}};
				elseif nLevels == 3
					levels = {A_level{pos},B_level{pos},E_level{pos}};%C_level{pos}};
				elseif nLevels == 4
					levels = {A_level{pos},B_level{pos},C_level{pos},E_level{pos}};
				elseif nLevels == 5
					levels = {A_level{pos},B_level{pos},C_level{pos},E_level{pos},D_level{pos}};
				end
				[p{pos},tbl{pos},stats{pos}] = anovan(data{pos},levels);
			end
			obj.analysis.p = p;
			obj.analysis.tbl = tbl;
			obj.analysis.stats = stats;

			obj.compositeANOVAFindex
		end
		function compositeANOVAFindex(obj)
			% 
			% 	Use after A2ElevelANOVA(obj,nLevels,idx, overwrite) to plot composite F-index
			% 
			c = obj.collatedResults(1).centers/1000;

			for ic = 1:numel(c)
				F_nm1(ic) = obj.analysis.tbl{1, ic}{2,6};
				F_n(ic) = obj.analysis.tbl{1, ic}{3,6};
				F_l(ic) = obj.analysis.tbl{1, ic}{4,6};
				F_s(ic) = obj.analysis.tbl{1, ic}{5,6};
				nm1Score(ic) = (F_nm1(ic) - F_n(ic))/(F_nm1(ic) + F_n(ic));
				sig_nm1(ic) = 1*(obj.analysis.tbl{1, ic}{2,7} < 0.05);
				nScore(ic) = (F_n(ic) - F_nm1(ic))/(F_n(ic) + F_nm1(ic));
				sig_n(ic) = 1*(obj.analysis.tbl{1, ic}{3,7} < 0.05);
				sig_l(ic) = 1*(obj.analysis.tbl{1, ic}{4,7} < 0.05);
				sig_s(ic) = 1*(obj.analysis.tbl{1, ic}{5,7} < 0.05);

				nInfluence(ic) = (F_n(ic))/(F_n(ic) + F_nm1(ic));
				nm1Influence(ic) = (F_nm1(ic))/(F_n(ic) + F_nm1(ic));

				n3Influence(ic) = (F_n(ic))/(F_n(ic) + F_nm1(ic) + F_l(ic));
				nm13Influence(ic) = (F_nm1(ic))/(F_n(ic) + F_nm1(ic) + F_l(ic));
				l3Influence(ic) = (F_l(ic))/(F_n(ic) + F_nm1(ic) + F_l(ic));


				n4Influence(ic) = (F_n(ic))/(F_n(ic) + F_nm1(ic) + F_l(ic) + F_s(ic));
				nm14Influence(ic) = (F_nm1(ic))/(F_n(ic) + F_nm1(ic) + F_l(ic) + F_s(ic));
				l4Influence(ic) = (F_l(ic))/(F_n(ic) + F_nm1(ic) + F_l(ic) + F_s(ic));
				s4Influence(ic) = (F_s(ic))/(F_n(ic) + F_nm1(ic) + F_l(ic) + F_s(ic));
			end
			sig_nm1(find(~sig_nm1))=nan;
			sig_n(~sig_n)=nan;
			% [f,ax] = makeStandardFigure(1,[1,1]);
			% plot(ax,(c+0.5*0), nm1Score, 'k-', 'displayname', 'n-1th selectivity')
			% hold(ax, 'on')
			% plot(ax,(c+0.5*0), sig_nm1, 'k-', 'displayname', 'n-1th p<0.05', 'linewidth',5)
			% plot(ax,(c+0.5*0), nScore, 'r-', 'displayname', 'nth selectivity')
			% plot(ax,(c+0.5*0), sig_n-2, 'r-', 'displayname', 'nth p<0.05', 'linewidth',5)
			% legend(ax,'show')
			% xlabel(ax,'time (s relative to lamp-off)')
			% ylabel(ax,'Selectivity Index')
			% xlim(ax,[(c(1)+0.5*0), (c(end)+0.5*0)])

			[f,ax] = makeStandardFigure(1,[1,1]);
			hold(ax, 'on')
			plot(ax,(c+0.5*0), nm1Influence, 'k-', 'displayname', 'n-1th selectivity')
			plot(ax,(c+0.5*0), sig_nm1, 'k-', 'displayname', 'n-1th p<0.05', 'linewidth',5)
			plot(ax,(c+0.5*0), nInfluence, 'r-', 'displayname', 'nth selectivity')
			plot(ax,(c+0.5*0), sig_n+0.1, 'r-', 'displayname', 'nth p<0.05', 'linewidth',5)
			xlabel(ax,'Time (s relative to lamp-off)')
			ylabel(ax,'Relative Influence')
			xlim(ax,[(c(1)+0.5*0), (c(end)+0.5*0)])

			[f,ax] = makeStandardFigure(1,[1,1]);
			hold(ax, 'on')
			plot(ax,(c+0.5*0), nm13Influence, 'k-', 'displayname', 'n-1th selectivity')
			plot(ax,(c+0.5*0), sig_nm1, 'k-', 'displayname', 'n-1th p<0.05', 'linewidth',5)
			plot(ax,(c+0.5*0), n3Influence, 'r-', 'displayname', 'nth selectivity')
			plot(ax,(c+0.5*0), sig_n+0.01, 'r-', 'displayname', 'nth p<0.05', 'linewidth',5)
			plot(ax,(c+0.5*0), l3Influence, 'b-', 'displayname', 'nth selectivity')
			plot(ax,(c+0.5*0), sig_l+0.02, 'b-', 'displayname', 'nth p<0.05', 'linewidth',5)
			xlabel(ax,'Time (s relative to lamp-off)')
			ylabel(ax,'Relative Influence')
			xlim(ax,[(c(1)+0.5*0), (c(end)+0.5*0)])

			[f,ax] = makeStandardFigure(1,[1,1]);
			hold(ax, 'on')
			plot(ax,(c+0.5*0), nm14Influence, 'k-', 'displayname', 'n-1th selectivity')
			plot(ax,(c+0.5*0), sig_nm1, 'k-', 'displayname', 'n-1th p<0.05', 'linewidth',5)
			plot(ax,(c+0.5*0), n4Influence, 'r-', 'displayname', 'nth selectivity')
			plot(ax,(c+0.5*0), sig_n+0.01, 'r-', 'displayname', 'nth p<0.05', 'linewidth',5)
			plot(ax,(c+0.5*0), l4Influence, 'b-', 'displayname', 'nth selectivity')
			plot(ax,(c+0.5*0), sig_l+0.02, 'b-', 'displayname', 'nth p<0.05', 'linewidth',5)
			plot(ax,(c+0.5*0), s4Influence, 'g-', 'displayname', 'nth selectivity')
			plot(ax,(c+0.5*0), sig_s+0.03, 'g-', 'displayname', 'nth p<0.05', 'linewidth',5)
			xlabel(ax,'Time (s relative to lamp-off)')
			ylabel(ax,'Relative Influence')
			xlim(ax,[(c(1)+0.5*0), (c(end)+0.5*0)])
		end
		function [results, F_nm1, F_n, nm1Score, nScore, sig_nm1, sig_n, centers,baselineWindow, dataScore] = slidingBaselineANOVA(obj,baselineLickMode,verbose)
			error('not Implemented')
			% 
			% 	For use with divergenceIndex obj!
			% 
			% 	baselineLickMode: 'off' (don't worry about lick, 2-way ANOVA)
			% 	baselineLickMode: 'exclude' (ignore trials with lick in baseline window, 2-way ANOVA)
			% 	baselineLickMode: 'include' (3-way anova including lick presence as predictor)
			% 
			% 	set datanotFstat to true if you want to use the data itself in the selectivity index. If want to use the Fstat, leave this false (default and original version)
			% 		this was added 2/23/2020, see Lab Notebook ppt on baseline analysis figure s4 finalization
			% 		(actually we get it for free, so just updating the methods. We want the dataScore output)
			% 
			if nargin < 2
				baselineLickMode = 'off';
			end
			if nargin < 3
				verbose = true;
			end
			centers = obj.collatedResults(1).Stat.centers;
			for ic = 1:numel(centers)
				results{ic} = obj.buildStaticBaselineANOVADataset(baselineWindow, obj.gFitLP.nMultibDFF.dFF, baselineLickMode, 'lampOff', centers(ic), false);
				% results{ic} = obj.buildStaticBaselineANOVADataset(baselineWindow, obj.gFitLP.nMultibDFF.dFF, 'exclude', 'lampOff', centers(ic), false);
				F_nm1(ic) = results{ic}.results.tbl{2,6};
				F_n(ic) = results{ic}.results.tbl{3,6};
				nm1Score(ic) = (F_nm1(ic) - F_n(ic))/(F_nm1(ic) + F_n(ic));
				sig_nm1(ic) = 1*(results{ic}.results.tbl{2,7} < 0.025);
				nScore(ic) = (F_n(ic) - F_nm1(ic))/(F_n(ic) + F_nm1(ic));
				sig_n(ic) = 1*(results{ic}.results.tbl{3,7} < 0.025);
				% 
				% 	Have to think here about good selectivity index from the data itself... 
				% 
			end
			sig_nm1(find(~sig_nm1))=nan;
			sig_n(~sig_n)=nan;
			if verbose
				figure
				plot((centers+0.5*0)./1000, nm1Score, 'ko-', 'displayname', 'n-1th selectivity')
				hold on
				plot((centers+0.5*0)./1000, sig_nm1, 'k-', 'displayname', 'n-1th p<0.025', 'linewidth',5)
				plot((centers+0.5*0)./1000, nScore, 'ro-', 'displayname', 'nth selectivity')
				plot((centers+0.5*0)./1000, sig_n-2, 'r-', 'displayname', 'nth p<0.025', 'linewidth',5)
				legend('show')
				xlabel('time (s relative to lamp-off)')
				ylabel('Selectivity Index')
				set(gca, 'fontsize',30)
				xlim([(centers(1)+0.5*0)./1000, (centers(end)+0.5*0)./1000])
				set(gcf,'color','w');
			end
		end
		function ANOVAparams = buildStaticBaselineANOVADataset(obj, baselineWindow, ts, baselineLickMode, refEvent, centerOffsetFromRefEvent,verbose)
			error('not Implemented')
			% 
			% 	Baseline window is relative to Lamp-Off here, unless specified
			% 		refEvent = 'lampOff' or 'cue'
			% 
			%  OLD VERSION:
			% 	Baseline window - : consider BEFORE event as baseline
			% 	Baseline window + : consider AFTER event as baseline
			%  NOW BASELINE WINDOW IS ALWAYS +, it's the center offset from ref that decides if we are before or after
			% 
			%  NEW VERSION:
			% 	centerOffsetFromRefEvent: This number of samples gets subtracted from the reference event to determine the centering of the sliding window
			% 		- : before ref event
			% 		+ : after ref event
			% 
			if nargin < 7
				verbose = true;
			end
			if nargin < 6
				centerOffsetFromRefEvent = -2500;
			end
			if nargin < 5
				refEvent = 'lampOff';
			end
			if nargin < 4
				baselineLickMode = 'exclude'; %'off', 'include'
			end
			if nargin < 3
				ts = obj.GLM.gfit;
				warning('Using 200-boxcar gfit from GLM struct')
			end
			baselineWindow = abs(baselineWindow);
			% 
			% 	We will find all trials in each 'cell' and record the mean of the baseline to the set
			% 
			% 	Factor A = (n-1 trial outcome)
			% 		level 1 = early
			% 		level 2 = rewarded
			% 
			% 	Factor B = (n trial outcome)
			% 		level 1 = early
			% 		level 2	= rewarded
			% -------------------------------------------------
        	% 
        	% 	Early vs Rew ranges in ms
        	% 
        	earlyRange = [700, 3330];
        	rewRange = [3334, 7000];
        	% earlyRange = [700, 2000];
        	% rewRange = [4000, 7000];
        	% 
            all_fl_wrtc_ms = zeros(numel(obj.GLM.lampOff_s), 1);
			all_fl_wrtc_ms(obj.GLM.fLick_trial_num) = obj.GLM.firstLick_s - obj.GLM.cue_s(obj.GLM.fLick_trial_num);
			all_fl_wrtc_ms = all_fl_wrtc_ms*1000/obj.Plot.samples_per_ms; % convert to ms
			allTrialIdx = 1:numel(all_fl_wrtc_ms);
        	% 
        	% 	We will find all trials fitting each Factor-level, then find intersections. 
        	% 	Then we will take the appropriate data from each set to make the cell-dataset
        	% 
        	ll = allTrialIdx(all_fl_wrtc_ms(1:end-1) >= earlyRange(1));
			ul = allTrialIdx(all_fl_wrtc_ms(1:end-1) <= earlyRange(2));
			A1 = ll(ismember(ll, ul));

        	ll = allTrialIdx(all_fl_wrtc_ms(1:end-1) >= rewRange(1));
			ul = allTrialIdx(all_fl_wrtc_ms(1:end-1) <= rewRange(2));
			A2 = ll(ismember(ll, ul));

        	ll = allTrialIdx(all_fl_wrtc_ms(2:end) >= earlyRange(1));
			ul = allTrialIdx(all_fl_wrtc_ms(2:end) <= earlyRange(2));
			B1 = ll(ismember(ll, ul));

        	ll = allTrialIdx(all_fl_wrtc_ms(2:end) >= rewRange(1));
			ul = allTrialIdx(all_fl_wrtc_ms(2:end) <= rewRange(2));
			B2 = ll(ismember(ll, ul));
			% NB! This is indexed as (trial n) - 1. Just for the intersections. Need to add 1 to get final trial n!
			% 
			% 	Find trial indicies for each cell.
			% 
			A1B1_idx = intersect(A1, B1) + 1;
			A2B1_idx = intersect(A2, B1) + 1;
			A1B2_idx = intersect(A1, B2) + 1;
			A2B2_idx = intersect(A2, B2) + 1;
			if isempty(A1B1_idx), error('(n-1) early (n) early factor doesn''t exist in data'),
			elseif isempty(A2B1_idx), error('(n-1) rew (n) early factor doesn''t exist in data'),
			elseif isempty(A1B2_idx), error('(n-1) early (n) rew factor doesn''t exist in data'),
			elseif isempty(A2B2_idx), error('(n-1) rew (n) rew factor doesn''t exist in data'), end
			% 
			% 	Get data for each cell:
			% 
			A_level = cell(1, numel(A1B1_idx)+numel(A2B1_idx));
			B_level = cell(1, numel(A1B1_idx)+numel(A1B2_idx));

			if strcmpi(refEvent, 'lampOff')
				if baselineWindow < 0
					error('this is obsolete. use for old version without the sliding window')
					A1B1 = nan(numel(A1B1_idx), 1);
					for iTrial = 1:numel(A1B1_idx)
						A1B1(iTrial) = mean(ts(obj.GLM.pos.lampOff(A1B1_idx(iTrial)) - abs(baselineWindow):obj.GLM.pos.lampOff(A1B1_idx(iTrial))));
						A_level{iTrial} = 'early n-1';
						B_level{iTrial} = 'early n';
					end
					A2B1 = nan(numel(A2B1_idx), 1);
					for iTrial = 1:numel(A2B1_idx)
						A2B1(iTrial) = mean(ts(obj.GLM.pos.lampOff(A2B1_idx(iTrial)) - abs(baselineWindow):obj.GLM.pos.lampOff(A2B1_idx(iTrial))));
						A_level{numel(A1B1_idx)+iTrial} = 'rewarded n-1';
						B_level{numel(A1B1_idx)+iTrial} = 'early n';
					end
					A1B2 = nan(numel(A1B2_idx), 1);
					for iTrial = 1:numel(A1B2_idx)
						A1B2(iTrial) = mean(ts(obj.GLM.pos.lampOff(A1B2_idx(iTrial)) - abs(baselineWindow):obj.GLM.pos.lampOff(A1B2_idx(iTrial))));
						A_level{numel(A1B1_idx)+numel(A2B1_idx)+iTrial} = 'early n-1';
						B_level{numel(A1B1_idx)+numel(A2B1_idx)+iTrial} = 'rewarded n';
					end
					A2B2 = nan(numel(A2B2_idx), 1);
					for iTrial = 1:numel(A2B2_idx)
						A2B2(iTrial) = mean(ts(obj.GLM.pos.lampOff(A2B2_idx(iTrial)) - abs(baselineWindow):obj.GLM.pos.lampOff(A2B2_idx(iTrial))));
						A_level{numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+iTrial} = 'rewarded n-1';
						B_level{numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+iTrial} = 'rewarded n';	
					end
				else
					A1B1 = nan(numel(A1B1_idx), 1);
					for iTrial = 1:numel(A1B1_idx)
						idxs = round(obj.GLM.pos.lampOff(A1B1_idx(iTrial))+centerOffsetFromRefEvent - abs(baselineWindow)/2 +1):round(obj.GLM.pos.lampOff(A1B1_idx(iTrial))+centerOffsetFromRefEvent + abs(baselineWindow)/2);
						A1B1(iTrial) = mean(ts(idxs));
						A_level{iTrial} = 'early n-1';
						B_level{iTrial} = 'early n';
					end
					A2B1 = nan(numel(A2B1_idx), 1);
					for iTrial = 1:numel(A2B1_idx)
						idxs = round(obj.GLM.pos.lampOff(A2B1_idx(iTrial))+centerOffsetFromRefEvent - abs(baselineWindow)/2 +1):round(obj.GLM.pos.lampOff(A2B1_idx(iTrial))+centerOffsetFromRefEvent + abs(baselineWindow)/2);
						A2B1(iTrial) = mean(ts(idxs));
						A_level{numel(A1B1_idx)+iTrial} = 'rewarded n-1';
						B_level{numel(A1B1_idx)+iTrial} = 'early n';
					end
					A1B2 = nan(numel(A1B2_idx), 1);
					for iTrial = 1:numel(A1B2_idx)
						idxs = round(obj.GLM.pos.lampOff(A1B2_idx(iTrial))+centerOffsetFromRefEvent - abs(baselineWindow)/2 +1):round(obj.GLM.pos.lampOff(A1B2_idx(iTrial))+centerOffsetFromRefEvent + abs(baselineWindow)/2);
						A1B2(iTrial) = mean(ts(idxs));
						A_level{numel(A1B1_idx)+numel(A2B1_idx)+iTrial} = 'early n-1';
						B_level{numel(A1B1_idx)+numel(A2B1_idx)+iTrial} = 'rewarded n';
					end
					A2B2 = nan(numel(A2B2_idx), 1);
					for iTrial = 1:numel(A2B2_idx)
						idxs = round(obj.GLM.pos.lampOff(A2B2_idx(iTrial))+centerOffsetFromRefEvent - abs(baselineWindow)/2 +1):round(obj.GLM.pos.lampOff(A2B2_idx(iTrial))+centerOffsetFromRefEvent + abs(baselineWindow)/2);
						A2B2(iTrial) = mean(ts(idxs));
						A_level{numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+iTrial} = 'rewarded n-1';
						B_level{numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+iTrial} = 'rewarded n';	
					end
				end
			elseif strcmpi(refEvent, 'cue')
				error('Not Implemented')
			else
				error('Not Implemented')
			end

			if strcmpi(baselineLickMode, 'include')
				% 
				% 	If using nLicks in baseline...
				% 
				if ~isfield(obj.GLM, 'pos') || ~isfield(obj.GLM.pos, 'lick'), obj.GLM.pos.lick = obj.getXPositionsWRTgfit(obj.GLM.lick_s);, end
				% 
				% 	Determine the indicies of the beginnings of each baseline Period...
				% 
				if strcmpi(refEvent, 'lampOff')
					if baselineWindow < 0
						error('obsolete')
						obj.GLM.pos.baselineStart = obj.GLM.pos.lampOff - abs(baselineWindow) + 1;
					else
% 						obj.GLM.pos.baselineStart = obj.GLM.pos.lampOff; 
                        obj.GLM.pos.baselineStart = round(obj.GLM.pos.lampOff + centerOffsetFromRefEvent - abs(baselineWindow)/2 + 1);
					end
				elseif strcmpi(refEvent, 'cue')
					error('Not Implemented')
				else
					error('Not Implemented')
				end
				obj.GLM.pos.baselineStart(end+1) = numel(obj.GLM.rawF); % tack on the full length so that we can correct the entire signal...
				
				nBaselineLicks = zeros(size(obj.GLM.pos.cue));

				if strcmpi(refEvent, 'lampOff')
					if baselineWindow < 0
						for iTrial = 1:numel(obj.GLM.pos.cue)
							nBaselineLicks(iTrial) = sum(ismember(find(obj.GLM.pos.lick > obj.GLM.pos.baselineStart(iTrial)), find(obj.GLM.pos.lick < obj.GLM.pos.baselineStart(iTrial)+abs(baselineWindow))));
	            		end
					else
						for iTrial = 1:numel(obj.GLM.pos.cue)
							nBaselineLicks(iTrial) = sum(ismember(find(obj.GLM.pos.lick > obj.GLM.pos.baselineStart(iTrial)), find(obj.GLM.pos.lick < obj.GLM.pos.baselineStart(iTrial)+baselineWindow)));
	            		end
					end
				elseif strcmpi(refEvent, 'cue')
					error('Not Implemented')
				else
					error('Not Implemented')
				end
				if verbose
		            figure 
		            bar(nBaselineLicks)
		            xlabel('Trial #')
		            ylabel('# licks in baseline')
	            end

	            trialOrder = [A1B1_idx,A2B1_idx,A1B2_idx,A2B2_idx];
	            lickLevel = zeros(size(trialOrder));
	            lickLevel(find(nBaselineLicks(trialOrder) ~= 0)) = 1;
                
	            if verbose
		            figure;
					hold on
	% 				C = linspecer(max(nBaselineLicks)*2);
	                C = colormap('hsv');
	                C = C(1:32, :);
	                colormap(C)
	                maxNLicks = max(nBaselineLicks(all_fl_wrtc_ms > 700 & all_fl_wrtc_ms < 7000));
	                caxis([1,maxNLicks])
	%                 st = max(nBaselineLicks);
	                allIdx = 1:numel(A1B1_idx);
					for iTrial = 1:numel(A1B1_idx)
						if lickLevel(allIdx(iTrial)) == ~(nBaselineLicks(A1B1_idx(iTrial)))
							error('Lick level is not correct!')
						end
						if nBaselineLicks(A1B1_idx(iTrial)) > 0
							plot(1+rand(1)/3, A1B1(iTrial), 'o', 'MarkerFaceColor', C(round(nBaselineLicks(A1B1_idx(iTrial))/maxNLicks*32), :) -0, 'MarkerEdgeColor', C(round(nBaselineLicks(A1B1_idx(iTrial))/maxNLicks*32), :), 'Markersize', 10);
						else
							plot(1+rand(1)/3, A1B1(iTrial), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end
					plot(1, mean(A1B1), 'ko', 'MarkerSize', 30);
	                allIdx = numel(A1B1_idx)+1:numel(A1B1_idx)+numel(A2B1_idx);
					for iTrial = 1:numel(A2B1_idx)
						if lickLevel(allIdx(iTrial)) ~= true(nBaselineLicks(A2B1_idx(iTrial)))
							error('Lick level is not correct!')
						end
						if nBaselineLicks(A2B1_idx(iTrial)) > 0
							plot(2+rand(1)/3, A2B1(iTrial), 'o', 'MarkerFaceColor', C(round(nBaselineLicks(A2B1_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A2B1_idx(iTrial))/maxNLicks*32),:), 'Markersize', 10);
						else
							plot(2+rand(1)/3, A2B1(iTrial), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end	
					plot(2, mean(A2B1), 'ko', 'MarkerSize', 30);
	                allIdx = numel(A1B1_idx)+numel(A2B1_idx)+1:numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx);
					for iTrial = 1:numel(A1B2_idx)
						if lickLevel(allIdx(iTrial)) ~= true(nBaselineLicks(A1B2_idx(iTrial)))
							error('Lick level is not correct!')
						end
						if nBaselineLicks(A1B2_idx(iTrial)) > 0
							plot(3+rand(1)/3, A1B2(iTrial), 'o', 'MarkerFaceColor', C(round(nBaselineLicks(A1B2_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A1B2_idx(iTrial))/maxNLicks*32),:), 'Markersize', 10);
						else
							plot(3+rand(1)/3, A1B2(iTrial), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end
					plot(3, mean(A1B2), 'ko', 'MarkerSize', 30);
					allIdx = numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+1:numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+numel(A2B2_idx);
					for iTrial = 1:numel(A2B2_idx)
						if lickLevel(allIdx(iTrial)) ~= true(nBaselineLicks(A2B2_idx(iTrial)))
							error('Lick level is not correct!')
						end
						if nBaselineLicks(A2B2_idx(iTrial)) > 0
							plot(4+rand(1)/3, A2B2(iTrial), 'o', 'MarkerFaceColor', C(round(nBaselineLicks(A2B2_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A2B2_idx(iTrial))/maxNLicks*32),:), 'Markersize', 10);
						else
							plot(4+rand(1)/3, A2B2(iTrial), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end						
					plot(4, mean(A2B2), 'ko', 'MarkerSize', 30);
					ax = gca;
					set(ax, 'fontsize', 20);
					ax.XAxis.TickValues = [1,2,3,4];
					ax.XAxis.TickLabels = {'(n-1) early|n early','(n-1) rew|n early','(n-1) early|n rew','(n-1) rew|n rew'};
				end
				% 
				% 	Run ANOVA
				% 
				ANOVAparams.earlyRange = earlyRange;
				ANOVAparams.rewRange = rewRange;
				ANOVAparams.A1 = '(n-1) early';
				ANOVAparams.A2 = '(n-1) rewarded';
				ANOVAparams.B1 = '(n) early';
				ANOVAparams.B2 = '(n) rewarded';
				ANOVAparams.L1 = 'No Baseline Licks';
				ANOVAparams.L2 = '+ Baseline Licks';
				ANOVAparams.factorIdx.A1 = A1;
				ANOVAparams.factorIdx.A2 = A2;
				ANOVAparams.factorIdx.B1 = B1;
				ANOVAparams.factorIdx.B2 = B2;
				ANOVAparams.factorIdx.nBaselineLicks = nBaselineLicks;
				ANOVAparams.cellIdx.A1B1 = A1B1_idx;
				ANOVAparams.cellIdx.A2B1 = A2B1_idx;
				ANOVAparams.cellIdx.A1B2 = A1B2_idx;
				ANOVAparams.cellIdx.A2B2 = A2B2_idx;
				ANOVAparams.cellData.A1B1 = A1B1;
				ANOVAparams.cellData.A2B1 = A2B1;
				ANOVAparams.cellData.A1B2 = A1B2;
				ANOVAparams.cellData.A2B2 = A2B2;
				ANOVAparams.cellData.A_level = A_level;
				ANOVAparams.cellData.B_level = B_level;
				ANOVAparams.cellData.lickLevel = lickLevel;
				ANOVAparams.cellData.data = [A1B1;A2B1;A1B2;A2B2];
				if verbose
					[ANOVAparams.results.p,ANOVAparams.results.tbl,ANOVAparams.results.stats,ANOVAparams.results.terms] = anovan(ANOVAparams.cellData.data, {A_level, B_level, lickLevel}, 'model','interaction');
				else
					[ANOVAparams.results.p,ANOVAparams.results.tbl,ANOVAparams.results.stats,ANOVAparams.results.terms] = anovan(ANOVAparams.cellData.data, {A_level, B_level, lickLevel}, 'model','interaction', 'display', 'off');
				end
			elseif strcmpi(baselineLickMode, 'exclude')
				% 
				% 	If using nLicks in baseline...
				% 
				if ~isfield(obj.GLM, 'pos') || ~isfield(obj.GLM.pos, 'lick'), obj.GLM.pos.lick = obj.getXPositionsWRTgfit(obj.GLM.lick_s);, end
				% 
				% 	Determine the indicies of the beginnings of each baseline Period...
				% 
				if strcmpi(refEvent, 'lampOff')
					if baselineWindow < 0
                        error('obsolete')
						obj.GLM.pos.baselineStart = obj.GLM.pos.lampOff - abs(baselineWindow) + 1;
					else
% 						obj.GLM.pos.baselineStart = obj.GLM.pos.lampOff;
                        obj.GLM.pos.baselineStart = round(obj.GLM.pos.lampOff + centerOffsetFromRefEvent - abs(baselineWindow)/2 + 1);
					end
				elseif strcmpi(refEvent, 'cue')
					error('Not Implemented')
				else
					error('Not Implemented')
				end
				
				obj.GLM.pos.baselineStart(end+1) = numel(obj.GLM.rawF); % tack on the full length so that we can correct the entire signal...
				
				nBaselineLicks = zeros(size(obj.GLM.pos.cue));
				for iTrial = 1:numel(obj.GLM.pos.cue)
					nBaselineLicks(iTrial) = sum(ismember(find(obj.GLM.pos.lick > obj.GLM.pos.baselineStart(iTrial)), find(obj.GLM.pos.lick < obj.GLM.pos.baselineStart(iTrial)+abs(baselineWindow))));
	            end
	            if verbose            	
		            figure 
		            bar(nBaselineLicks)
		            xlabel('Trial #')
		            ylabel('# licks in baseline')
	            end
	            % 
	            %	Look at all trials
	            %
	            trialOrder = [A1B1_idx,A2B1_idx,A1B2_idx,A2B2_idx];
	            lickLevel = zeros(size(trialOrder));
	            lickLevel(find(nBaselineLicks(trialOrder) ~= 0)) = 1;
                if verbose
		            figure;
		            ax1 = subplot(1,2,1);
		            ax2 = subplot(1,2,2);
					hold(ax1, 'on');
					hold(ax2, 'on');
	% 				C = linspecer(max(nBaselineLicks)*2);
	                C = colormap('hsv');
	                C = C(1:32, :);
	                colormap(C)
	                maxNLicks = max(nBaselineLicks(all_fl_wrtc_ms > 700 & all_fl_wrtc_ms < 7000));
	                caxis(ax1, [1,maxNLicks])
	%                 st = max(nBaselineLicks);
	                allIdx = 1:numel(A1B1_idx);
	                % rng(1);
					for iTrial = 1:numel(A1B1_idx)
						if lickLevel(allIdx(iTrial)) == ~(nBaselineLicks(A1B1_idx(iTrial)))
							error('Lick level is not correct!')
						end
						if nBaselineLicks(A1B1_idx(iTrial)) > 0
							scatter3(ax1, 1+rand(1)/3, A1B1(iTrial), A1B1_idx(iTrial), 30, 'o', 'filled', 'MarkerFaceColor', C(round(nBaselineLicks(A1B1_idx(iTrial))/maxNLicks*32), :) -0, 'MarkerEdgeColor', C(round(nBaselineLicks(A1B1_idx(iTrial))/maxNLicks*32), :));
						else
							scatter3(ax1, 1+rand(1)/3, A1B1(iTrial), A1B1_idx(iTrial), 15, 'o', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end
					scatter3(ax1,1, mean(A1B1), 1, 1000, 'ko', 'LineWidth', 3);
	                allIdx = numel(A1B1_idx)+1:numel(A1B1_idx)+numel(A2B1_idx);
					for iTrial = 1:numel(A2B1_idx)
						if lickLevel(allIdx(iTrial)) ~= true(nBaselineLicks(A2B1_idx(iTrial)))
							error('Lick level is not correct!')
						end
						if nBaselineLicks(A2B1_idx(iTrial)) > 0
							scatter3(ax1, 2+rand(1)/3, A2B1(iTrial), A2B1_idx(iTrial), 30, 'o', 'filled', 'MarkerFaceColor', C(round(nBaselineLicks(A2B1_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A2B1_idx(iTrial))/maxNLicks*32),:));
						else
							scatter3(ax1, 2+rand(1)/3, A2B1(iTrial), A2B1_idx(iTrial), 15, 'o', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end	
					scatter3(ax1,2, mean(A2B1), 1, 1000, 'ko', 'LineWidth', 3);
	                allIdx = numel(A1B1_idx)+numel(A2B1_idx)+1:numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx);
					for iTrial = 1:numel(A1B2_idx)
						if lickLevel(allIdx(iTrial)) ~= true(nBaselineLicks(A1B2_idx(iTrial)))
							error('Lick level is not correct!')
						end
						if nBaselineLicks(A1B2_idx(iTrial)) > 0
							scatter3(ax1, 3+rand(1)/3, A1B2(iTrial), A1B2_idx(iTrial), 30, 'o', 'MarkerFaceColor', C(round(nBaselineLicks(A1B2_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A1B2_idx(iTrial))/maxNLicks*32),:));
						else
							scatter3(ax1, 3+rand(1)/3, A1B2(iTrial), A1B2_idx(iTrial), 15, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end
					scatter3(ax1,3, mean(A1B2), 1, 1000, 'ko', 'LineWidth', 3);
					allIdx = numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+1:numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+numel(A2B2_idx);
					for iTrial = 1:numel(A2B2_idx)
						if lickLevel(allIdx(iTrial)) ~= true(nBaselineLicks(A2B2_idx(iTrial)))
							error('Lick level is not correct!')
						end
						if nBaselineLicks(A2B2_idx(iTrial)) > 0
							scatter3(ax1, 4+rand(1)/3, A2B2(iTrial), A2B2_idx(iTrial), 30, 'o', 'MarkerFaceColor', C(round(nBaselineLicks(A2B2_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A2B2_idx(iTrial))/maxNLicks*32),:));
						else
							scatter3(ax1, 4+rand(1)/3, A2B2(iTrial), A2B2_idx(iTrial), 15, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end						
					scatter3(ax1,4, mean(A2B2), 1, 1000, 'ko', 'LineWidth', 3);
					set(ax1, 'fontsize', 16);
					ax1.XAxis.TickValues = [1,2,3,4];
					ax1.XAxis.TickLabels = {'(n-1) early|n early','(n-1) rew|n early','(n-1) early|n rew','(n-1) rew|n rew'};
				end
				% 
				% 	Remove trials with licks in ITI
				% 
	            A1B1(nBaselineLicks(A1B1_idx) > 0) = [];
				A2B1(nBaselineLicks(A2B1_idx) > 0) = [];
				A1B2(nBaselineLicks(A1B2_idx) > 0) = [];
				A2B2(nBaselineLicks(A2B2_idx) > 0) = [];
				A1B1_idx(nBaselineLicks(A1B1_idx) > 0) = [];
				A2B1_idx(nBaselineLicks(A2B1_idx) > 0) = [];
				A1B2_idx(nBaselineLicks(A1B2_idx) > 0) = [];
				A2B2_idx(nBaselineLicks(A2B2_idx) > 0) = [];
				A_level(nBaselineLicks(trialOrder) > 0) = [];
				B_level(nBaselineLicks(trialOrder) > 0) = [];
				trialOrder = [A1B1_idx,A2B1_idx,A1B2_idx,A2B2_idx];
				% 
				allIdx = 1:numel(A1B1_idx);
                % rng(1);
                if verbose
					for iTrial = 1:numel(A1B1_idx)
						if nBaselineLicks(A1B1_idx(iTrial)) > 0
							scatter3(ax2, 1+rand(1)/3, A1B1(iTrial), A1B1_idx(iTrial), 30, 'o', 'filled', 'MarkerFaceColor', C(round(nBaselineLicks(A1B1_idx(iTrial))/maxNLicks*32), :) -0, 'MarkerEdgeColor', C(round(nBaselineLicks(A1B1_idx(iTrial))/maxNLicks*32), :));
						else
							scatter3(ax2, 1+rand(1)/3, A1B1(iTrial), A1B1_idx(iTrial), 15, 'o', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end
					scatter3(ax2,1, mean(A1B1), 1, 1000, 'ko', 'LineWidth', 3);
	                allIdx = numel(A1B1_idx)+1:numel(A1B1_idx)+numel(A2B1_idx);
					for iTrial = 1:numel(A2B1_idx)
						if nBaselineLicks(A2B1_idx(iTrial)) > 0
							scatter3(ax2, 2+rand(1)/3, A2B1(iTrial), A2B1_idx(iTrial), 30, 'o', 'filled', 'MarkerFaceColor', C(round(nBaselineLicks(A2B1_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A2B1_idx(iTrial))/maxNLicks*32),:));
						else
							scatter3(ax2, 2+rand(1)/3, A2B1(iTrial), A2B1_idx(iTrial), 15, 'o', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end	
					scatter3(ax2,2, mean(A2B1), 1, 1000, 'ko', 'LineWidth', 3);
	                allIdx = numel(A1B1_idx)+numel(A2B1_idx)+1:numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx);
					for iTrial = 1:numel(A1B2_idx)
						if nBaselineLicks(A1B2_idx(iTrial)) > 0
							scatter3(ax2, 3+rand(1)/3, A1B2(iTrial), A1B2_idx(iTrial), 30, 'o', 'MarkerFaceColor', C(round(nBaselineLicks(A1B2_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A1B2_idx(iTrial))/maxNLicks*32),:));
						else
							scatter3(ax2, 3+rand(1)/3, A1B2(iTrial), A1B2_idx(iTrial), 15, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end
					scatter3(ax2,3, mean(A1B2), 1, 1000, 'ko', 'LineWidth', 3);
					allIdx = numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+1:numel(A1B1_idx)+numel(A2B1_idx)+numel(A1B2_idx)+numel(A2B2_idx);
					for iTrial = 1:numel(A2B2_idx)
						if nBaselineLicks(A2B2_idx(iTrial)) > 0
							scatter3(ax2, 4+rand(1)/3, A2B2(iTrial), A2B2_idx(iTrial), 30, 'o', 'MarkerFaceColor', C(round(nBaselineLicks(A2B2_idx(iTrial))/maxNLicks*32), :), 'MarkerEdgeColor', C(round(nBaselineLicks(A2B2_idx(iTrial))/maxNLicks*32),:));
						else
							scatter3(ax2, 4+rand(1)/3, A2B2(iTrial), A2B2_idx(iTrial), 15, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
						end
					end						
					scatter3(ax2,4, mean(A2B2), 1, 1000, 'ko', 'LineWidth', 3);
					set(ax2, 'fontsize', 16);
					ax2.XAxis.TickValues = [1,2,3,4];
					ax2.XAxis.TickLabels = {'(n-1) early|n early','(n-1) rew|n early','(n-1) early|n rew','(n-1) rew|n rew'};
				end
				% 
				% 	Run ANOVA
				% 
				ANOVAparams.earlyRange = earlyRange;
				ANOVAparams.rewRange = rewRange;
				ANOVAparams.A1 = '(n-1) early';
				ANOVAparams.A2 = '(n-1) rewarded';
				ANOVAparams.B1 = '(n) early';
				ANOVAparams.B2 = '(n) rewarded';
				ANOVAparams.L1 = 'No Baseline Licks';
				ANOVAparams.L2 = '+ Baseline Licks';
				ANOVAparams.factorIdx.A1 = A1;
				ANOVAparams.factorIdx.A2 = A2;
				ANOVAparams.factorIdx.B1 = B1;
				ANOVAparams.factorIdx.B2 = B2;
				ANOVAparams.factorIdx.nBaselineLicks = nBaselineLicks;
				ANOVAparams.cellIdx.A1B1 = A1B1_idx;
				ANOVAparams.cellIdx.A2B1 = A2B1_idx;
				ANOVAparams.cellIdx.A1B2 = A1B2_idx;
				ANOVAparams.cellIdx.A2B2 = A2B2_idx;
				ANOVAparams.cellData.A1B1 = A1B1;
				ANOVAparams.cellData.A2B1 = A2B1;
				ANOVAparams.cellData.A1B2 = A1B2;
				ANOVAparams.cellData.A2B2 = A2B2;
				ANOVAparams.cellData.A_level = A_level;
				ANOVAparams.cellData.B_level = B_level;
				ANOVAparams.cellData.data = [A1B1;A2B1;A1B2;A2B2];
				if verbose
					[ANOVAparams.results.p,ANOVAparams.results.tbl,ANOVAparams.results.stats,ANOVAparams.results.terms] = anovan(ANOVAparams.cellData.data, {A_level, B_level}, 'model','interaction');
				else
					[ANOVAparams.results.p,ANOVAparams.results.tbl,ANOVAparams.results.stats,ANOVAparams.results.terms] = anovan(ANOVAparams.cellData.data, {A_level, B_level}, 'model','interaction','display','off');
				end
					
				% rng('default');
            else
				% 
				% 	Plot the means for each cell
				% 
				if verbose
					figure;
					hold on
					plot(ones(numel(A1B1_idx)), A1B1, 'o');
					plot(1, mean(A1B1), 'ko', 'MarkerSize', 30);
					plot(2*ones(numel(A2B1_idx)), A2B1, 'o');
					plot(2, mean(A2B1), 'ko', 'MarkerSize', 30);
					plot(3*ones(numel(A1B2_idx)), A1B2, 'o');
					plot(3, mean(A1B2), 'ko', 'MarkerSize', 30);
					plot(4*ones(numel(A2B2)), A2B2, 'o');
					plot(4, mean(A2B2), 'ko', 'MarkerSize', 30);
					ax = gca;
					set(ax, 'fontsize', 20);
					ax.XAxis.TickValues = [1,2,3,4];
					ax.XAxis.TickLabels = {'(n-1) early|n early','(n-1) rew|n early','(n-1) early|n rew','(n-1) rew|n rew'};
				end
				% 
				% 	Run ANOVA
				% 
				ANOVAparams.earlyRange = earlyRange;
				ANOVAparams.rewRange = rewRange;
				ANOVAparams.A1 = '(n-1) early';
				ANOVAparams.A2 = '(n-1) rewarded';
				ANOVAparams.B1 = '(n) early';
				ANOVAparams.B2 = '(n) rewarded';
				ANOVAparams.factorIdx.A1 = A1;
				ANOVAparams.factorIdx.A2 = A2;
				ANOVAparams.factorIdx.B1 = B1;
				ANOVAparams.factorIdx.B2 = B2;
				ANOVAparams.cellIdx.A1B1 = A1B1_idx;
				ANOVAparams.cellIdx.A2B1 = A2B1_idx;
				ANOVAparams.cellIdx.A1B2 = A1B2_idx;
				ANOVAparams.cellIdx.A2B2 = A2B2_idx;
				ANOVAparams.cellData.A1B1 = A1B1;
				ANOVAparams.cellData.A2B1 = A2B1;
				ANOVAparams.cellData.A1B2 = A1B2;
				ANOVAparams.cellData.A2B2 = A2B2;
				ANOVAparams.cellData.A_level = A_level;
				ANOVAparams.cellData.B_level = B_level;
				ANOVAparams.cellData.data = [A1B1;A2B1;A1B2;A2B2];
				if verbose
					[ANOVAparams.results.p,ANOVAparams.results.tbl,ANOVAparams.results.stats,ANOVAparams.results.terms] = anovan(ANOVAparams.cellData.data, {A_level, B_level}, 'model','interaction');
				else
					[ANOVAparams.results.p,ANOVAparams.results.tbl,ANOVAparams.results.stats,ANOVAparams.results.terms] = anovan(ANOVAparams.cellData.data, {A_level, B_level}, 'model','interaction', 'display', 'off');
				end
			end
			% 
			% 	Calculate expected power of test for dataset
			% 
			% k_prime = 2; % number of levels of each factor
			% n_primeA = numel(A1B1_idx) + numel(A2B1_idx);
			% ssA = 
			% phiA = sqrt()
		end










		function progressBar(obj, iter, total, nested, cutter)
			if nargin < 5
				cutter = 1000;
			end
			if nargin < 4
				nested = false;
			end
			if nested
				prefix = '		';
			else
				prefix = '';
			end
			if rem(iter,total*.1) == 0 || rem(iter, cutter) == 0
				done = {'=', '=', '=', '=', '=', '=', '=', '=', '=', '='};
				incomplete = {'-', '-', '-', '-', '-', '-', '-', '-', '-', '-'};
				ndone = round(iter/total * 10);
				nincomp = round((1 - iter/total) * 10);
				disp([prefix '	*' horzcat(done{1:ndone}) horzcat(incomplete{1:nincomp}) '	(' num2str(iter) '/' num2str(total) ') ' datestr(now)]);
			end
		end
		function pathstr = correctPathOS(obj,pathstr)
			if ispc
    			pathstr = strjoin(strsplit(pathstr, '/'), '\');
			else
				pathstr = [strjoin(strsplit(pathstr, '\'), '/')];
			end
		end

		function save(obj)
			ID = obj.iv.runID;
			savefilename = ['CollatedStatAnalysisObj_' obj.iv.collateKey, '_' datestr(now, 'YYYYmmDD_HH_MM') '_runIDno' num2str(ID)];
			save([savefilename, '.mat'], 'obj', '-v7.3');
		end


		
	end
end



