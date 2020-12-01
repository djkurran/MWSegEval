classdef GUI
    % Class gui
    % Methods gather user input parameters, validate input data, initialize
    % program parameters, monitor and report status of tasks.
    % Data structure holds command and threshold values that are input to the
    % toolbox, and holds status of tasks that are output from the toolbox to
    % be displayed by the GUI. 
    %======================================================================    
    properties        
        segmenationIniValues                % initial threshold values used to segment forward and reconstruction images using the thresholding techhnique
        task                                % MWSegEval image segmentation and evaluation toolbox task requested by GUI
        taskStatus                          % Status of task passed to GUI
        finc                                % Incident field frequency read from configuration file used to collect forward and inverse data
        contrastIteration                   % Iteration requested by user that reconstructed image formed by inversion algorithm of reconstructed image
        contoursThrhld                      % level of the contour extraction threshold
        contoursMaxNm                       % Maximum number of contours to extract from mask
        dataPath                            % Points to configuration and data folders related to reconstruction data   
        autoSegmentation                    % Read from GUI - true - use ML autosegmentation technique to sesgment images, otherwise use thresholding technique
        finiteElement_flnm                  % File name to store finite element object in the results
        input_csi_Img_Rec_flnm              % Configuration file related to FEM-CSI reconstructed image.
        input_datacollect_flnm              % Configuration file related to FEM forward model.        
        input_forwardModel_flnm             % Complex permittivity profile of forward model for rectangular mesh input models that are in a matlab format (i.e., 2D matrix of complex permittivity values).
        input_reconstructedModel_flnm       % Complex permittivity profile of reconstructed model for rectangular mesh input models that are in a matlab format (i.e., 2D matrix of complex permittivity values).               
        output_flnm                         % CSV file name of metric values for geometric and dielectric analysis
        refSeg_flnm                         % File name storing segmentation object reference model data 
        recSeg_flnm                         % File name storing segmentation object reconstructed model data 
        performanceMetrics_flnm             % File name storing the performanceMetrics object data            
        errorFlag                           % Error flag set by a GUI method or one of the image processing tasks to signal that an error has occurred
        figFormat                           % Format to store figures (i.e., -fig, -png, -eps)
        figColorMap                         % Color map to use to display images (i.e., jet, parula, gray)             
    end
    %======================================================================    
    methods
        %------------------------------------------------------------------        
        function obj = GUI(task)
            % Constructor method to create instance of object. When instance
            % object is created, the method initialize some of the properties
            % of the GUI object. It also adds paths to folders holding image
            % analysis code and code to export Matlab figures to -png, and -eps formats.
            % Paths are on the top of the search path for the current matlab session.
            % Inputs:
            %   task - (string) task that the task manager is to carry out.            
            %--------------------------------------------------------------
            currentPath = pwd;
            addpath(currentPath);
            addpath([currentPath '\export_fig']);
            obj.task = task;            
            obj.input_csi_Img_Rec_flnm = 'input_csi_Img_Rec.xml';
            obj.input_datacollect_flnm = 'input_datacollect.xml';
            obj.input_forwardModel_flnm = 'input_forwardModel.mat';
            obj.input_reconstructedModel_flnm = 'input_reconstructedModel.mat';
            obj.output_flnm = 'output_performanceEvaluatoin.csv';
            obj.finiteElement_flnm = 'inputModels.mat';
            obj.refSeg_flnm = 'referenceImagesSegmentationResults.mat';
            obj.recSeg_flnm = 'recructedImagesSegmentationResults.mat';
            obj.performanceMetrics_flnm = 'performanceMetrics.mat';
            obj.errorFlag = false;
        end % Constructor GUI
        %------------------------------------------------------------------
        function obj = checkGUIinput(obj, femCsiIn, fileAndFigProp, contourAnalysisIn, segmentationThresholdIn )            
            % This function checks and interprets input arguments from GUI. 
            % Carries out error checking of configuration and data files
            % associated with the forward and reconstructed models.
            % Inputs:
            %   femCsiIn - structure array (plane organization) - femCsiIn.contrastIteration - iteration of the contrast file
            %   fileAndFigProp - structure array (plane organization) -
            %                    fileAndFigProp.datapath - points to directory containing reconstruction and forward data and configuration files
            %                    fileAndFigProp.inputModelFormat - 
            %                    fileAndFigProp.epsFileSave - boolean indicating whether the figure should be save in eps format
            %                    fileAndFigProp.pngFileSave - boolean indicating whether the figure should be save in png format
            %                    fileAndFigProp.figFileSave - boolean indicating whether the figure should be save in fig (matlab) format
            %                    fileAndFigProp.jetColorMap - boolean indicating whether the jet color map is used for figure.
            %                    fileAndFigProp.parulaColorMap - boolean indicating whether the parula color map is used for figure.
            %                    fileAndFigProp.grayColorMap - boolean indicating whether the gray color map is used for figure.
            %                    fileAndFigProp.doNotSave - boolean indicating whether the the figure is saved
            %                    
            %   contourAnalysisIn - structure array (plane organization) -
            %                    contourAnalysisIn.MaxNmCntr - string selected by user setting maximum numbers of contours to be extracted from reconstructed masks indicating the maximum number of contours to extract from the reconstructed mask.
            %                    contourAnalysisIn.CntrTh - integer contour threshold expressed in percent. All contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis.
            %   segmentationThresholdIn - structure array (plane organization)
            %                    segmentationThresholdIn.SegmentationuseAutomatic - string - 'yes' or 'No' - yes means that an unsupervised machine learning technique that is reinforced with hypothesis testing and statistical inference to iteratively segment the ROI into fat, transition, fibroglandular, and malignant tissues.  
            %
            %--------------------------------------------------------------
            % First set object properties with user values collected by the
            % GUI.
            obj.dataPath = fileAndFigProp.datapath;
            obj.contoursThrhld = contourAnalysisIn.CntrTh;
            obj.contoursMaxNm = str2double(contourAnalysisIn.MaxNmCntr);
            
            if strcmp(segmentationThresholdIn.SegmentationuseAutomatic,'Yes')
                obj.autoSegmentation = true;
            else
                obj.autoSegmentation = false;
            end
            
            if ischar(obj.dataPath) && length(obj.dataPath)>=3
                if ~(strcmp(obj.dataPath(length(obj.dataPath)),'/') || strcmp(obj.dataPath(length(obj.dataPath)),'\'))
                    obj.dataPath = [obj.dataPath '\'];
                end
            end            
            if fileAndFigProp.epsFileSave
                obj.figFormat = '-eps';
            end
            if fileAndFigProp.pngFileSave
                obj.figFormat = '-png';
            end
            if fileAndFigProp.figFileSave
                obj.figFormat = '-fig';
            end
            if fileAndFigProp.doNotSave
                obj.figFormat = '-DoNotSave';
            end            
            if fileAndFigProp.jetColorMap
                obj.figColorMap = 'jet';
            end
            if fileAndFigProp.parulaColorMap
                obj.figColorMap = 'parula';
            end
            if fileAndFigProp.grayColorMap
                obj.figColorMap = 'gray';
            end
            %--------------------------------------------------------------
            % Error checking
            % check if valid data path has been entered.
            if ~isdir(obj.dataPath)
                obj.taskStatus = 'Error: Invalid forward and reconstruction models data path';
                obj.errorFlag = true;  
                return
            else
                % Create figures and results folders to store results
                if ~exist([obj.dataPath 'figures'],'dir')
                    status = mkdir([obj.dataPath '\figures']);                
                end
                if ~exist([obj.dataPath 'results'],'dir')
                    status = mkdir([obj.dataPath '\results']);                
                end
                % check if data folder exists
                if ~exist([obj.dataPath 'data'],'dir')
                    obj.errorFlag = true;
                    obj.taskStatus = 'Error: inversion data does not exist';
                    return
                end
                % the forward and inverse models are to be converted to a
                % rectangular mesh. Check that the x and y axis are
                % specified.
                if ~exist([obj.dataPath 'Configuration/xNodes.txt'],'file')
                    obj.errorFlag = true;
                    errorMsg = ['Error: Transformed imaging domain file xNodes.txt does not exist'];
                    obj.taskStatus = errorMsg; 
                    return
                end
                if ~exist([obj.dataPath 'Configuration/yNodes.txt'],'file')
                    obj.errorFlag = true;
                    errorMsg = ['Error: Transformed imaging domain file yNodes.txt does not exist'];
                    obj.taskStatus = errorMsg; 
                    return
                end
                % Check that the incident field frequency has been
                % specified
                if ~exist([obj.dataPath 'Configuration/frequency.txt'],'file')
                    obj.errorFlag = true;
                    errorMsg = ['Error: Incident field frequency file frequency.txt does not exist'];
                    obj.taskStatus = errorMsg; 
                    return
                end                 
                
                % if the triangular mesh input format is selected, then get
                % the iteration of the contruast profile to be reconstructed.
                if strcmp(fileAndFigProp.inputModelFormat, 'Triangular mesh')
                    
                    obj.contrastIteration = femCsiIn.contrastIteration;
                    % Check that configuration files exist for the forward
                    % and inverse models.                    
                    if ~exist([obj.dataPath 'Configuration\' obj.input_csi_Img_Rec_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: FEM-CSI inversion file ',obj.input_csi_Img_Rec_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg; 
                        return
                    end
                    if ~exist([obj.dataPath 'Configuration\' obj.input_datacollect_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Foward model file ',obj.input_datacollect_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg; 
                        return
                    end                   
                else
                    % if the rectangular mesh format has been seleccted,
                    % then check that forward and inverse data files exist.
                    obj.contrastIteration = 0;                    
                    if ~exist([obj.dataPath 'data\' obj.input_forwardModel_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Forward model data file ',obj.input_forwardModel_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg; 
                        return
                    end
                    
                    if ~exist([obj.dataPath 'data\' obj.input_reconstructedModel_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Reconstructed model data file ',obj.input_reconstructedModel_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg; 
                        return
                    end
                    
                end                
            end
            
            switch obj.task
                % If the user has selected triangular format input, then
                % check that the contrast file that the user has selected
                % exists.
                case 'Input forward and reconstructed models for processing'
                    if strcmp(fileAndFigProp.inputModelFormat, 'Triangular mesh')
                        contrast_flnm = ['contrast' num2str(obj.contrastIteration) '.txt'];
                        if ~exist([obj.dataPath 'data\' contrast_flnm],'file')
                            obj.errorFlag = true;
                            errorMsg = ['Error: Contrast file ', contrast_flnm, ' does not exist'];
                            obj.taskStatus = errorMsg;                        
                        end
                    end                                    
                % if the user has selected the segmentation task, ensure that the finite element object constructed during the input task exists.    
                case 'Segment forward and reconstructed models into tissue masks'
                    if ~exist([obj.dataPath 'results\' obj.finiteElement_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.finiteElement_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end
                % if the user selects the analysze tissue interface task, then ensure that the finite element and segmentation objects exist.     
                case 'Analyze tissue interfaces'
                    if ~exist([obj.dataPath 'results\' obj.finiteElement_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.finiteElement_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end                    
                    if ~exist([obj.dataPath 'results\' obj.refSeg_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.refSeg_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end
                    if ~exist([obj.dataPath 'results\' obj.recSeg_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.recSeg_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end
                % if the user selects the analysze geometric properties task, then ensure that the finite element, segmentation, and performanceMetrics objects exist.         
                case 'Analyze geometric properties of reconstructed tissues'
                    if ~exist([obj.dataPath 'results\' obj.finiteElement_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.finiteElement_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end
                    if ~exist([obj.dataPath 'results\' obj.refSeg_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.refSeg_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                       
                    end
                    if ~exist([obj.dataPath 'results\' obj.recSeg_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.recSeg_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;
                       
                    end
                    if ~exist([obj.dataPath 'results\' obj.performanceMetrics_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.performanceMetrics_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end                    
                % if the user selects the analysze dielectric properties task, then ensure that the finite element, segmentation, and performanceMetrics objects exist.             
                case 'Analyze dielectric properties of tissues within reconstructed masks'
                    if ~exist([obj.dataPath 'results\' obj.refSeg_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.refSeg_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end
                    if ~exist([obj.dataPath 'results\' obj.recSeg_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.recSeg_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end
                    if ~exist([obj.dataPath 'results\' obj.performanceMetrics_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.performanceMetrics_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end
                % if the user selects the output csv file task, then ensure that the performanceMetrics object exist.             
                case 'Output csv file of metric values for each region'
                    if ~exist([obj.dataPath 'results\' obj.performanceMetrics_flnm],'file')
                        obj.errorFlag = true;
                        errorMsg = ['Error: Results file ', obj.performanceMetrics_flnm, ' does not exist'];
                        obj.taskStatus = errorMsg;                        
                    end
            end                                
            
        end % function checkGUIinput
        %------------------------------------------------------------------
        function obj = checkTissueInterfacesFwdData(obj, fldList, refModel)
            % This function checks that masks and interfaces have been
            % extracted from the glandular and tumor regions of the forward model.            
            % Inputs:
            %   fldList - structure array of strings - 'glandular' and 'tumor'
            %   refModel - refModel.masks.glandular, and refModel.masks.tumor - holds tisue mask structure of glandular and malignant tissues, respectively.
            %            - refModel.interfaces.glandular, and refModel.interfaces.tumor - holds tissue interfaces of glandular and malignant tissues, respectively.   
            % Outputs:
            %   obj.errorFlag - boolean - true - error has occurred; false error has not occurred. 
            %   obj.errorMsg - string - error message to be displayed by GUI
            %--------------------------------------------------------------
            for i=1:length(fldList)
                fldnm = fldList{i};
                try
                    idref = find(refModel.masks.(fldnm),1);
                catch % exception
                    errorMsg = ['Error: No ', fldnm, ' region present in forward model'];
                    display(errorMsg);                
                    obj.errorFlag = true;                
                    obj.taskStatus = errorMsg;                    
                    return
                end
                if isempty(idref)
                    errorMsg = ['Error: No ', fldnm, ' region present in forward model'];
                    display(errorMsg);                
                    obj.errorFlag = true;                
                    obj.taskStatus = errorMsg;                   
                    return
                end
                try
                    idref = find(refModel.interfaces.(fldnm){1, 1}.pts,1);
                catch % exception
                    errorMsg = ['Error: No ', fldnm, ' region present in forward model'];
                    display(errorMsg);                
                    obj.errorFlag = true;                
                    obj.taskStatus = errorMsg;                    
                    return
                end
                if isempty(idref)
                    errorMsg = ['Error: No ', fldnm, ' region present in forward model'];
                    display(errorMsg);                
                    obj.errorFlag = true;                
                    obj.taskStatus = errorMsg;                    
                    return
                end
            end
        end % function checkTissueInterfacesFwdData
        %------------------------------------------------------------------
        function obj = checkTissueInterfacesInvData(obj, fldList, recModel)
            % This function checks to ensure that masks and interfaces have been
            % extracted from the glandular and tumor regions of the inverses model.            
            % Inputs:
            %   fldList - structure array of strings - 'glandular' and 'tumor'
            %   refModel - refModel.masks.glandular, and refModel.masks.tumor - holds tisue mask structure of glandular and malignant tissues, respectively.
            %            - refModel.interfaces.glandular, and refModel.interfaces.tumor - holds tissue interfaces of glandular and malignant tissues, respectively.   
            % Outputs:
            %   obj.errorFlag - boolean - true - error has occurred; false error has not occurred. 
            %   obj.errorMsg - string - error message to be displayed by GUI
            %--------------------------------------------------------------
            for i=1:length(fldList)
                fldnm = fldList{i};
                try
                    idref = find(recModel.masks.(fldnm),1);
                catch % exception
                    errorMsg = ['Error: No ', fldnm, ' region present in forward model'];
                    display(errorMsg);                
                    obj.errorFlag = true;                
                    obj.taskStatus = errorMsg;
                    return
                end
                if isempty(idref)
                    errorMsg = ['Error: No ', fldnm, ' region present in forward model'];
                    display(errorMsg);                
                    obj.errorFlag = true;                
                    obj.taskStatus = errorMsg;
                    return
                end
                try
                    idref = find(recModel.interfaces.(fldnm){1, 1}.pts,1);
                catch % exception
                    errorMsg = ['Error: No ', fldnm, ' region present in forward model'];
                    display(errorMsg);                
                    obj.errorFlag = true;                
                    obj.taskStatus = errorMsg;
                    return
                end
                if isempty(idref)
                    errorMsg = ['Error: No ', fldnm, ' region present in forward model'];
                    display(errorMsg);                
                    obj.errorFlag = true;                
                    obj.taskStatus = errorMsg;
                    return
                end
            end
        end % function checkTissueInterfacesInvData
        %------------------------------------------------------------------
        function obj = errorHandling(obj, errorMsg)
            % This function sets the errorFlag and passes the error message
            % to the GUI to display the error message to the user.
            % Inputs:
            %   errorMsg - string - error message received from task to be
            %   passed to GUI.
            % Outputs:
            %   obj.errorFlag - boolean - true - an error has occurred
            %   obj.taskStatus - string - message to be displayed by GUI in
            %   the GUI status field.
            %--------------------------------------------------------------
            display(errorMsg);                
            obj.errorFlag = true;                
            obj.taskStatus = errorMsg;
        end % function checkGUIinput
        %------------------------------------------------------------------
        function obj = iniSegmentationThresholds( obj, finc )
            % This function is called by the task manager when the 'Input forward
            % and reconstructed models for processing' task is requested.
            % For a given tissue type, this function calls the
            % getComplexPermittivity() function which returns the complex permittivity
            % for the tissue type. The returned value is used for the segmenation threshold
            % for the tissue type.
            % Inputs:
            %   finc - float - incident field frequency (Hz)            
            % Outputs:
            %   obj.segmenationIniValues.fwd_skin_epr_Hi_Limit  - skin region threshold limit of forward model               
            %   obj.segmenationIniValues.femcsi_fat_epr_Hi_Limit - Real component of complex permittivity (or equivalently the relative permittivity) Hi limit of adipose tissue of inverse model
            %   obj.segmenationIniValues.fwd_fat_epr_Hi_Limit - Real component of complex permittivity (or equivalently the relative permittivity) Hi limit of adipose tissue of forward model
            %   obj.segmenationIniValues.femcsi_fat_imag_Hi_Limit - Imaginary component of complex permittivity Hi limit of adipose tissue of inverse model
            %   obj.segmenationIniValues.femcsi_fat_mag_Hi_Limit - Magnitude of complex permittivity Hi limit of adipose tissue of inverse model
            %   obj.segmenationIniValues.fwd_gland_contours_threshold - contour threshold limit of glandular tissue of forward model
            %   obj.segmenationIniValues.fwd_gland_contours_maxNm  - maximum number of contours to extract from glandular region of forward model
            %   obj.segmenationIniValues.femcsi_gland_epr_Lo_Limit - Real component of complex permittivity (or equivalently the relative permittivity) Lo limit of glandular tissue (i.e., 25 percentile) of inverse model
            %   obj.segmenationIniValues.fwd_gland_epr_Lo_Limit - - Real component of complex permittivity (or equivalently the relative permittivity) Lo limit of glandular tissue (i.e., 25 percentile) of forward model
            %   obj.segmenationIniValues.femcsi_gland_imag_Lo_Limit - Imaginary component of complex permittivity Lo limit of glandular tissue (i.e., 25 percentile) of inverse model
            %   obj.segmenationIniValues.femcsi_gland_mag_Lo_Limit - Magnitude component of complex permittivity Lo limit of glandular tissue (i.e., 25 percentile) of inverse model
            %   obj.segmenationIniValues.fwd_tumor_contours_threshold - contour threshold limit of malignant tissue of forward model
            %   obj.segmenationIniValues.fwd_tumor_contours_maxNm - maximum number of contours to extract from tumor region of forward model
            %   obj.segmenationIniValues.fwd_tumor_epr_Lo_Limit  - Real component of complex permittivity Lo limit of malignant tissue of forward model
            %   obj.segmenationIniValues.femcsi_tumor_epr_Lo_Limit - percentage of maximum value (real, imaginary, or magnitude) of ROI of inverse model
            %--------------------------------------------------------------
            
            % get the incident field frequency read by the configuration
            % file by the finite element object
            obj.finc = finc;
            
            % skin regoin threshold limits    
            [skin_epr_Hi_Limit, skin_sigma_Hi_Limit, skin_imag_Hi_Limit] = obj.getComplexPermittivity( 'dry_skin', obj.finc );             
            obj.segmenationIniValues.fwd_skin_epr_Hi_Limit = skin_epr_Hi_Limit;
        
            % fat regoin threshold limits    
            [fat_epr_Hi_Limit, fat_sigma_Hi_Limit, fat_imag_Hi_Limit] = obj.getComplexPermittivity( 'group_3_75%', obj.finc ); 

            % Real component of complex permittivity (or equivalently the relative
            % permittivity) Hi limit of adipose tissue
            obj.segmenationIniValues.femcsi_fat_epr_Hi_Limit = fat_epr_Hi_Limit + 3; 
            obj.segmenationIniValues.fwd_fat_epr_Hi_Limit = obj.segmenationIniValues.femcsi_fat_epr_Hi_Limit;

            % Imaginary component of complex permittivity Hi limit of adipose tissue          
            obj.segmenationIniValues.femcsi_fat_imag_Hi_Limit = -fat_imag_Hi_Limit + 0.72;   

            % Magnitude of complex permittivity Hi limit of adipose tissue           
            obj.segmenationIniValues.femcsi_fat_mag_Hi_Limit = sqrt(obj.segmenationIniValues.femcsi_fat_epr_Hi_Limit^2 +  obj.segmenationIniValues.femcsi_fat_imag_Hi_Limit^2);

            % transition regoin threshold limits
            % contour threshold limit             
            obj.segmenationIniValues.fwd_gland_contours_threshold = 10;
            % maximum number of countours 
            obj.segmenationIniValues.fwd_gland_contours_maxNm = 10;            
            [femcsi_gland_epr_Lo_Limit, femcsi_gland_sigma_Lo_Limit, femcsi_gland_imag_Lo_Limit] = obj.getComplexPermittivity( 'group_1_25%', obj.finc ); 

            % Real component of complex permittivity (or equivalently the relative
            % permittivity) Lo limit of group 1 (25 percentile) tissue
            obj.segmenationIniValues.femcsi_gland_epr_Lo_Limit = femcsi_gland_epr_Lo_Limit - 3; 
            obj.segmenationIniValues.fwd_gland_epr_Lo_Limit = obj.segmenationIniValues.femcsi_gland_epr_Lo_Limit;

            % Imaginary component of complex permittivity Lo limit of group 1 (25 percentile) tissue            
            obj.segmenationIniValues.femcsi_gland_imag_Lo_Limit = -femcsi_gland_imag_Lo_Limit + 0.72;

            % Magnitude of complex permittivity Lo limit of group 1 (25 percentile) tissue          
            obj.segmenationIniValues.femcsi_gland_mag_Lo_Limit = sqrt(obj.segmenationIniValues.femcsi_gland_epr_Lo_Limit^2 +  obj.segmenationIniValues.femcsi_gland_imag_Lo_Limit^2);

            % tumor regoin threshold limits
            % contour threshold limit                       
            obj.segmenationIniValues.fwd_tumor_contours_threshold = 10;
            % maximum number of countours 
            obj.segmenationIniValues.fwd_tumor_contours_maxNm = 10;
            
            [tumor_epr_Lo_Limit, tumor_sigma_Lo_Limit, tumor_imag_Lo_Limit] = obj.getComplexPermittivity( 'tumor_25%', obj.finc );
            % Real component of complex permittivity (or equivalently the relative
            % permittivity)
            obj.segmenationIniValues.fwd_tumor_epr_Lo_Limit = tumor_epr_Lo_Limit - 1;
            obj.segmenationIniValues.femcsi_tumor_epr_Lo_Limit = 90;
        end
        %------------------------------------------------------------------       
        function saveResults(obj, results)
            % This function is called by task manager to save objects at the completion
            % of a task.
            % Inputs:
            %   results - object - fe, [refModel recModel], or pa - depending on the task being run
            %--------------------------------------------------------------
            switch obj.task
                case 'Input forward and reconstructed models for processing'                    
                    fe = results;                    
                    save([obj.dataPath 'results/' obj.finiteElement_flnm], 'fe');                    
                case 'Segment forward and reconstructed models into tissue masks'
                    refModel = results(1);
                    recModel = results(2);
                    save([obj.dataPath 'results/' obj.refSeg_flnm], 'refModel');
                    save([obj.dataPath 'results/' obj.recSeg_flnm], 'recModel');
                otherwise
                    pa = results;
                    save([obj.dataPath 'results/' obj.performanceMetrics_flnm], 'pa');
            end
        end
        %------------------------------------------------------------------        
        function [er_real, sigma, imagSigma] = getComplexPermittivity(obj, tissueType, freq, debyeParam)
        % For a given tissue type, this function returns the Debye parameters to the calling function. The
        % parameters of the breast are based on Lazebnik et. al. "A large-scale study of the ultrawideband microwave
        % dielectric properties of normal, benign and malignant breast tissues obtained from cancer surgeries",
        % Phys. Med. Biol. 52 (2007) 6093–6115. The dielctric properties of the skin are based on Winters D, Bond E, Van Veen B, Hagness S. Estimation of the frequency-
        % dependent average dielectric properties of breast tissue using a time-domain inverse scattering technique. IEEE Trans Antennas Propag. 2006;54:3517–3528.
        %
        % Inputs:
        %   tissueType - (string) tissue type to look-up the Debye parameters.
        %   freq - (scaler) frequency of incident field (Hz) (e.g., input 1e9 if the incident field is 1GHz).        
        %   DebyeParam - (arry of floats) in the form [epr_inf ,epr_stat , sigma_stat, tau]
        %
        % Outputs:
        %   er_real - real part of the complex permittivity corresponding to the relative permittivity.
        %   sigma - effective conductivity (S/m)
        %   imagSigma - imaginary part of the complex permittivity        
        %
        % Note that if freq is a vector then er_real, sigma, and imagSigma are both vectors
        %
        % Written on January 17, 2014 by D. Kurrant (djkurran@ucalgary.ca)
        %   $Revision:   $Date: 2020/11/11 7:06 a.m. $ - Revised for inclusion as method for GUI object.
        %===================================================================================================
            % Variable initialization
            number_of_input_variables = nargin;
            % incident field frequency in rad/s
            w = 2*pi*freq;
            % permittivity of free space (F/m)
            e_o = 8.85419e-12;
            % permeability of free space
            mu_0 = 0.0000012566370614359;
            if number_of_input_variables == 4
                tissueType = 'custom';    
            end

            %----------------------------------------------------------------------------------------------------
            switch tissueType
                case 'custom'
                    % DebyeParam = [eps_inf, delta_eps, static_cond, tau]
                    DebyeParam.epr_inf = debyeParam(1);
                    DebyeParam.epr_stat = debyeParam(2);
                    DebyeParam.sigma_stat = debyeParam(3);
                    DebyeParam.tau = debyeParam(4);  
                    DebyeParam.alpha= 0.0;

                case 'group_1_25%'        
                    % Group1 parameters 0-30 (i.e., 'dense' fibroglandular tissues 25 percentile)
                    DebyeParam.epr_inf = 9.941;
                    DebyeParam.epr_stat = 26.60 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.462;
                    DebyeParam.tau = 10.90e-12;
                    DebyeParam.alpha= 0.0;

                case 'group_1_50%'        
                    % Group1 parameters 0-30 (i.e., 'dense' fibroglandular tissues - 50 percentile)
                    DebyeParam.epr_inf = 7.821;
                    DebyeParam.epr_stat = 41.48 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.713;
                    DebyeParam.tau = 10.66e-12;
                    DebyeParam.alpha = 0.000;        

                case 'group_1_75%'        
                    % Group1 parameters 0-30 (i.e., 'dense' fibroglandular tissues 75 percentile)
                    DebyeParam.epr_inf = 6.151;
                    DebyeParam.epr_stat = 48.26 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.809;
                    DebyeParam.tau = 10.26e-12;
                    DebyeParam.alpha = 0.0;

                case 'group_2_25%'         
                    % Group2 parameters 31-85 (i.e., 'medium dense' fibroglandular tissues 25 percentile)
                    DebyeParam.epr_inf = 8.718;
                    DebyeParam.epr_stat = 17.51 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.293;
                    DebyeParam.tau = 13.17e-12;
                    DebyeParam.alpha =  0.0;

                case 'group_2_50%'        
                    % Group2 parameters 31-85 (i.e., 'medium dense' fibroglandular tissues - 50 percentile)
                    DebyeParam.epr_inf = 5.573;
                    DebyeParam.epr_stat = 34.57 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.524;
                    DebyeParam.tau = 9.149e-12;
                    DebyeParam.alpha =  0;  

                case 'group_2_75%'        
                    % Group2 parameters 31-85 (i.e., 'medium dense' fibroglandular tissues 75 percentile)
                    DebyeParam.epr_inf = 5.157;
                    DebyeParam.epr_stat = 45.81 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.766;
                    DebyeParam.tau = 8.731e-12;
                    DebyeParam.alpha =  0;

                case 'group_3_min'        
                    % adipose tissue - minimum       
                    DebyeParam.epr_inf = 2.293;
                    DebyeParam.epr_stat = 0.141 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.0020;
                    DebyeParam.tau = 16.4e-12;
                    DebyeParam.alpha = 0.251;       

                case 'group_3_25%'        
                    % Group3 parameters 85-100 (i.e., adipose tissues 25 percentile)
                    DebyeParam.epr_inf = 2.908;
                    DebyeParam.epr_stat = 1.200 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.020;
                    DebyeParam.tau = 16.88e-12;
                    DebyeParam.alpha = 0.0;

                case 'group_3_50%'        
                    % Group3 parameters 85-100 (i.e., adipose tissues 50 percentile)
                    DebyeParam.epr_inf = 3.14;
                    DebyeParam.epr_stat = 1.708 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.036;
                    DebyeParam.tau = 14.65e-12;
                    DebyeParam.alpha = 0;          

                case 'group_3_75%'        
                    % Group3 parameters 85-100 (i.e., adipose tissues 75 percentile)
                    DebyeParam.epr_inf = 4.031;
                    DebyeParam.epr_stat = 3.654 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.083;
                    DebyeParam.tau = 14.12e-12;
                    DebyeParam.alpha = 0.0;

                case 'tumor_25%'        
                    % tumor - 25 percentile
                    DebyeParam.epr_inf = 7.670;
                    DebyeParam.epr_stat = 43.92 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.748;
                    DebyeParam.tau = 10.70e-12;
                    DebyeParam.alpha = 0.028;               

                case 'tumor_50%'                
                    % tumor - 50 percentile
                    DebyeParam.epr_inf = 6.749;
                    DebyeParam.epr_stat = 50.09 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.794;
                    DebyeParam.tau = 10.50e-12;
                    DebyeParam.alpha = 0.00;     

                case 'tumor_75%'                
                    % tumor - 75 percentile
                    DebyeParam.epr_inf = 9.058;
                    DebyeParam.epr_stat = 51.31 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.899;
                    DebyeParam.tau = 10.84e-12;
                    DebyeParam.alpha = 0.022;

                case 'dry_skin' % winters high (1.0 - 11.1 GHz)
                    DebyeParam.epr_inf = 4.0;
                    DebyeParam.epr_stat = 33.0 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 1.1;
                    DebyeParam.tau = 7.23e-12;
                    DebyeParam.alpha = 0.00;

                case 'dry_skin_low' % winters low (500 MHz - 2.5 GHz)
                    DebyeParam.epr_inf = 32.9;
                    DebyeParam.epr_stat = 10.7 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.668;
                    DebyeParam.tau = 77.4e-12;
                    DebyeParam.alpha = 0.00;

                case 'wet_skin' % winters low (500 MHz - 2.5 GHz)
                    DebyeParam.epr_inf = 33.9;
                    DebyeParam.epr_stat = 13.5 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.658;
                    DebyeParam.tau = 51.0e-12;
                    DebyeParam.alpha = 0.00;   

                case 'damp_skin' % winters low (500 MHz - 2.5 GHz)
                    DebyeParam.epr_inf = 33.4;
                    DebyeParam.epr_stat = 12.1 + DebyeParam.epr_inf;
                    DebyeParam.sigma_stat = 0.663;
                    DebyeParam.tau = 64.2e-12;
                    DebyeParam.alpha = 0.00;  

                otherwise
                    er_real = 1;
                    sigma = 0;    
                    return        
            end
            
            er_Debye = DebyeParam.epr_inf + (DebyeParam.epr_stat - DebyeParam.epr_inf)./(1+(1i*w*DebyeParam.tau).^(1-DebyeParam.alpha))+ DebyeParam.sigma_stat./(1i*w*e_o);            
            er_real=real(er_Debye);
            imagSigma = imag(er_Debye);
            sigma=-w.*e_o.*imag(er_Debye);
            
        end
        %------------------------------------------------------------------        
    end % methods
    %======================================================================
end %Class GUI