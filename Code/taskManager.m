function gui = taskManager(task, femCsiIn, contourAnalysisIn, segmentationThresholdIn, fileAndFigProp)
% This function reads user parameters; validates information; implements a
% state machine to run requested task. It also coordinates all error
% handing and image processing, segmentation, and evaluation tasks.
% The tasks performed are as follows:
% (1) Input forward and reconstructed models for processing.
% (2) Methods segment forward and reconstructed images to form  tissue masks,
%     map masks to tissue types, and sample boundaries of masks to form
%     tissue interfaces.
% (3) Analyze interfaces extracted from inverse model, compare with interfaces
%     extracted from forward model.
% (4) Analyze geometric properties of inverse tissue mask.
% (5) Aanalyze dielectric properties within reconstructed tissues; and
% (6) Write metric values to CSV file. 
close all;

% Check and interpret input arguments from GUI; carry out error checking of
% configuration and data files associated with the forward and reconstructed
% models.
gui = GUI(task);
gui = gui.checkGUIinput( femCsiIn, fileAndFigProp, contourAnalysisIn, segmentationThresholdIn);
if gui.errorFlag
    return
end
%-----------------------------------------------------------------------------------------------------------------------------------
switch task
%==============================================================================================================                  
    case 'Input forward and reconstructed models for processing'             
        disp('Forward and reconstructed models are being input for processing');
        % Construct finite element object used to carry out the task and
        % store the results.
        fe = finiteElement( gui.contrastIteration, fileAndFigProp.inputModelFormat);
        % If the forward and inverse models have triangular mesh elements,
        % then transform to square mesh elements. Store information in
        % finiete element object structure.
        if strcmp(fileAndFigProp.inputModelFormat, 'Triangular mesh')
            fe = fe.tri2square( gui );
        else
            % Otherwise the forward and inverse models have
            % rectangular mesh elements, so read the model information and
            % save in finite element object structure.
            fe = fe.rectangularMesh( gui );
        end        
        if fe.errorFlag
            % Error handling if any errors occur reading model information,
            % transforming informatoin, or storing information.
            gui = gui.errorHandling(fe.status);
            return
        end
        % Automatically initialize threshold values used to segment forward 
        % and inverse models using thresholding technique. The threshold
        % values are evaluated with the Debye single pole model at the incident
        % field frequency. Debye model parameters for various tissue types
        % (e.g., adipose tissue, malignant tissue, etc.) are used for the
        % calculations.
        gui = gui.iniSegmentationThresholds(fe.finc);
        % Store the calculated threshold values in the finite element
        % object
        fe = fe.iniSegmentationFwdThresholds(gui);
        
        % Display the forward and inverse models for each component (real,
        % imaginary, and magnitude). Save figure if requested by user.
        fe.imageFwdAndInvModels(gui);
        
        % Save finite element object in results folder. The data stored in
        % the object structure is used in the tissue segmentation step.
        gui.saveResults(fe);        
        gui.taskStatus = 'Forward and reconstructed models have been input';        
 % end 'Input forward and reconstructed models for processing'                     
%==============================================================================================================                  
    case 'Segment forward and reconstructed models into tissue masks'
        disp('Segmenting forward and reconstructed models into tissue masks');
        % Load finite element object into workspace. The forward and
        % inverse models are stored in its structure.
        load([gui.dataPath 'results/' gui.finiteElement_flnm]);
        
        % Segment forward model into tissue types. 
        % The thresholding technique is implemented. It uses the threshold
        % values calculated in the previous step with the single pole Debye
        % model at the incident field frequency.
        refModel = segmentation(fe.n, fe.m, gui.autoSegmentation, fe.xNodes, fe.yNodes);        
        refModel = refModel.segmentImmersionMedium(fe.fwdMdl.complexPermittivity);        
        refModel = refModel.segmentSkinRegion(fe.fwdMdl.eprArray, fe.segmenationFwdValues.fwd_skin_epr_Hi_Limit, segmentationThresholdIn.skinThkOffset);        
        refModel = refModel.segmentfwdFatRegion(fe.fwdMdl.eprArray, fe.segmenationFwdValues.fwd_fat_epr_Hi_Limit);
        refModel = refModel.segmentfwdTumorRegion(fe.fwdMdl.eprArray, fe.segmenationFwdValues.fwd_tumor_epr_Lo_Limit, fe.segmenationFwdValues.fwd_tumor_contours_threshold, fe.segmenationFwdValues.fwd_tumor_contours_maxNm);
        refModel = refModel.segmentfwdGlandRegion(fe.fwdMdl.eprArray, fe.segmenationFwdValues.fwd_gland_epr_Lo_Limit);                        
        refModel = refModel.segmentfwdTransitionRegion(fe.segmenationFwdValues.fwd_gland_contours_threshold, fe.segmenationFwdValues.fwd_gland_contours_maxNm);
        refModel = refModel.mapFwdMasksToTissueTypes();        
        
        % Segment the real, imaginary and magnitude parts of the 
        % inverse model to form tissue masks. Map masks to tissue type maps.
        % Fist segment immersion medium and skin from breast and construct
        % a region-of-interest within the interior of the breast that is slightly
        % offset from the skin/fat interface.
        recModel = segmentation(fe.n, fe.m, gui.autoSegmentation, fe.xNodes, fe.yNodes);
        recModel = recModel.segmentInvImmersionMedium(refModel.masks.immersionMedium, refModel.masks.interior );
        recModel = recModel.segmentInvSkinRegion(refModel.masks.skin);
        % Segment the inversion model using either a thresholding
        % method, or an automatic method
        if ~gui.autoSegmentation
            % If the thresholding method is chosen, then regions within the ROI
            % of the reconstructed model are segmented using the same
            % threshold values used to segment the forward model. The user
            % may refine these values through the GUI.
            recModel = recModel.segmentinvFatRegion_Threshold(fe.invMdl.complexPermittivity, segmentationThresholdIn);            
            recModel = recModel.segmentinvTumorRegion_Threshold(fe.invMdl.complexPermittivity, segmentationThresholdIn, gui);
            recModel = recModel.segmentinvFibroGlandRegion_Threshold(fe.invMdl.complexPermittivity, segmentationThresholdIn);                    
            recModel = recModel.segmentinvTransitionRegion_Threshold(gui.contoursThrhld, gui.contoursMaxNm);                
            recModel = recModel.mapInvMasksToTissueTypes();
        else
            % Otherwise, use the iterative unsupervised machine learning
            % technique that is reinforced with hypothesis testing and
            % statistical inference.
            
            % The method starts by setting the immersion medium
            % and skin as background, while the interior is set to the 
            % region of interest (ROI).
            ROI = recModel.masks.interior - recModel.masks.skin;
            recModel = recModel.getROI(ROI);
            % Set everything outside the ROI as background.
            recModel = recModel.setBackbackground();                 
            %==============================================================================================================    
            for lx = 1:length(recModel.cmpList)            
                % Select the image related to item in the image processing list = {(real(eps), img(eps),
                % mag(eps), img(Chi), mag(eps), mag(Chi)}. The segmentation
                % method is applied to the region of interest.
                recModel = recModel.getImageFromList(lx, fe.invMdl.complexPermittivity);                       
                recModel = recModel.iniNmCl(); % initialize number of clusters to 3.                
                recModel = recModel.getOutsideTumorMask(); % Initialize outside tumor mask to interior and null tumor mask
                recModel = recModel.iniHypothesis(); % initialize by rejecting Null hypothesis. 
                recModel = recModel.getROITissueDist();  % get ROI - comprised of the breast interior minus the skin layer                        
                while(recModel.h)
                    % Apply k-means clustering to ROI to partition the ROI
                    % into k clusters. 
                    recModel = recModel.getInitialClusterCentroids();
                    recModel = recModel.partitionROI();                                
                    
                    % Construct an estimate of the tumor region by assigning
                    % the largest valued cluster (i.e., max(k)) to the tumor mask.
                    % Use mask to extract tissues from ROI and evaluate the 
                    % distrubution of this tissue.
                    recModel = recModel.getTumorMask();                    
                    recModel = recModel.getTumorTissue();
                    recModel = recModel.getTumorTissueDistribution();                    
                    
                    % Construct mask from region within the ROI but outside
                    % the estimate of the tumor region which is the union of 
                    % cluster 2 to max(k)-1. Use mask to extract tissue
                    % from ROI and evaluate the distribution of these
                    % tissues.
                    recModel = recModel.getOutsideTumorMask();
                    recModel = recModel.getTissueOutsideTumorRegion(); 
                    recModel = recModel.getOutsideTumorTissueDistribution();
                    
                    % If k > 3, then test the hypothesis that the dielectric properties
                    % (complex permittivity or contrast) within the ROI but outside the tumor region
                    % of the present iteration and the dielectric properties within the 
                    % ROI but outside the tumor region of the previous iteration are
                    % from populations that follow the same distribution at
                    % the 5 percent significance level.
                    % Likewise, test a second hypothesis that the dielectric properties
                    % within tumor region of present iteration and the dielectric properties within  
                    % within tumor region of previous iteration are from populations that follow
                    % same distribution at 5 percent significance level.
                    % If both hypothesis are supported, then the ROI has
                    % been partitioned into k clusters.                  
                    recModel = recModel.testHypothesis();
                    % If the Null hypothesis (Ho) is rejected, then increment the number of
                    % clusters (k = k + 1), and repeat the process by partitioning the ROI
                    % with k = k + 1 clusters.                        
                    recModel = recModel.incrementNmCl(); 
                    display([recModel.cmpList{lx} ' auto segmentation ' num2str(recModel.nmCl) ' clusters']);
                    % Contintue while null hypothsis is rejected.
                end % clustering of ROI complete
                
                % Map clusters to the tissue masks that are used to
                % analyze the interfaces, geometric properties, and
                % dielectric properties of each region of tissue.
                recModel = recModel.mapClustersToMasks(recModel.cmpList{lx});  
                
                % Extract interfaces from malignant and glandular tissue
                % regions (i.e., masks) and remove all countours that are
                % less than a pre-defined threshold of area of the largest
                % region, allow a pre-selected maximum number contours to be extracted.                 
                recModel = recModel.rejectTumorRegionArtifacts(recModel.cmpList{lx}, gui);
                recModel = recModel.rejectGlandularRegionArtifacts(recModel.cmpList{lx}, gui);
                
                % Finally, map the tissue masks to tissue types used to
                % analyze the interfaces, geometric properties, and
                % dielectric properties of each region of tissue.
                recModel = recModel.mapAutoSegmentMasksToTissueTypes(recModel.cmpList{lx});                
            end %  get the next image from image list (cmpList)
        end % Automatic segmentation done when all images from list have been segmented. 
        
        % Display the forward and inverse models and the corresponding
        % tissue types extracted from the images.
        refModel.imageFwdTissueTypes(fe.xNodes, fe.yNodes, gui.dataPath, gui.figFormat);
        recModel.imageInvTissueTypes(fe.xNodes, fe.yNodes, gui.dataPath, gui.figFormat);        
        
        % Display the tissue masks for the forward and inverse models.
        recModel.imageMasks( fe.xNodes, fe.yNodes, refModel, gui.dataPath, gui.figFormat);
        if gui.autoSegmentation
            recModel.imageAutoTissueMaps( refModel, fe, gui)
        else
            recModel.imageThreshTissueMaps( refModel, fe, gui);
        end
        
        % Display the interfaces extracted from the segmented masks
        recModel.imgInterfaceOnFwdMdl( fe, gui, refModel);        
        
        gui.saveResults([refModel recModel]);        
        gui.taskStatus = 'Forward and reconstructed models have been segmented into tissue masks';
 % end 'Segment forward and reconstructed models into tissue masks' task                    
%==============================================================================================================            
    case 'Analyze tissue interfaces' 
        disp('Analysing tisue interfaces');
        % Load forward and inverse models stored in finite element object into
        % workspace.
        load([gui.dataPath 'results/' gui.finiteElement_flnm]);
        % Load the tissue masks and interfaces of the reference masks stored
        % in segmentation objects into workspace.
        load([gui.dataPath 'results/' gui.refSeg_flnm]);
        load([gui.dataPath 'results/' gui.recSeg_flnm]);        
        
        % Construct the performanceMetrics object
        pa = performanceMetrics(fe.xNodes, fe.yNodes, fe.n, fe.m);
        % Check to ensure there are mask and interface data extracted from 
        % forward model.
        fldListfwd = {'glandular', 'tumor'};
        gui = gui.checkTissueInterfacesFwdData(fldListfwd, refModel);
        if gui.errorFlag
            return
        end
        
        % Check ito ensure there are masks and interface data extracted from
        % the inverse model
        fldListInv = {'glandularRe', 'glandularIm', 'glandularMag', 'tumorRe', 'tumorIm', 'tumorMag'};
        gui = gui.checkTissueInterfacesInvData(fldListInv, recModel);
        if gui.errorFlag
            return
        end        
        
        % Evaluate the performance of the reconstruction algorithm to
        % extract the glandular and tumor tissue interfaces.
        % These edge points extracted from the surface of each glandular and 
        % tumor mask collectively represent estimated interfaces of the reconstructed
        % regions. The estimated points are compared with the corresponding
        % interfaces of the reference mask to evaluate how accurately the interface
        % is extracted by the reconstruction algorithm. It is evaluated by
        % measuring the distance each edge point on a reconstructed interface
        % is to the nearest point on the corresponding reference surface. 
        % Note that the reconstructed interface is compared with only those
        % interfaces from reference regions for which the estimated interface
        % touches or intersects. 
        i = 1;
        for j =1:6
            if j < 4
                pa = pa.interfaceAnalysis(fe, fldListfwd{i}, fldListInv{j},  refModel, recModel, gui );
            else
                if j == 4
                    i = i + 1;
                end
                pa = pa.interfaceAnalysis(fe, fldListfwd{i}, fldListInv{j},  refModel, recModel, gui);                            
            end
        end
        gui.saveResults(pa);        
        gui.taskStatus = 'Analysis of interfaces extracted from MW data complete.';
        % end 'Analyze tissue interfaces' task
%==============================================================================================================        
    case 'Analyze geometric properties of reconstructed tissues'
        disp('Analysing geometric properties of reconstructed tissue masks');
        % Load forward and inverse models stored in finite element object into
        % workspace.
        load([gui.dataPath 'results/' gui.finiteElement_flnm]);
        % Load the tissue masks and interfaces of the reference masks stored
        % in segmentation objects into workspace.
        load([gui.dataPath 'results/' gui.refSeg_flnm]);
        load([gui.dataPath 'results/' gui.recSeg_flnm]);       
        % Load the performance metrics object saved from the interface
        % analysis.        
        load([gui.dataPath 'results/' gui.performanceMetrics_flnm]);
        
        % Check to ensure there are mask and interface data extracted from 
        % forward model.
        fldListfwd = {'glandular', 'tumor'};
        gui = gui.checkTissueInterfacesFwdData(fldListfwd, refModel);
        if gui.errorFlag
            return
        end
        
        % Check ito ensure there are masks and interface data extracted from
        % the inverse model
        fldListInv = {'glandularRe', 'glandularIm', 'glandularMag', 'tumorRe', 'tumorIm', 'tumorMag'};
        gui = gui.checkTissueInterfacesInvData(fldListInv, recModel);
        if gui.errorFlag
            return
        end
        
        % The segmentation objects hold the reference and reconstructed masks
        % in their data structures that contain geometric property information.
        % Therefore, metrics are applied to these regions to provide measures
        % that evaluate the accuracy with which the geometric properties of
        % the underlying structures are reconstructed. When applying the 
        % functions, the masks are evaluated within a region of interest
        % defined as the box that bounds the union between the reference
        % and reconstructed masks.
        % The analysis accommodates cases in which there is a requirement
        % for the reconstructed region to be compared with many small isolated
        % reference regions scattered throughout the interior. For this
        % scenario, the reconstructed mask is only compared with a reference
        % mask if the masks overlap.
        i = 1;
        for j =1:6
            if j < 4
                pa = pa.geometricPropertyAnalysis(fe, fldListfwd{i}, fldListInv{j},  refModel, recModel, gui );
            else
                if j == 4
                    i = i + 1;
                end
                pa = pa.geometricPropertyAnalysis(fe, fldListfwd{i}, fldListInv{j},  refModel, recModel, gui );
            end
        end
        
        % Save results in pa object        
        gui.saveResults(pa);        
        gui.taskStatus = 'Analysis of geometric properties complete.';        
%==============================================================================================================
    case 'Analyze dielectric properties of tissues within reconstructed masks'
        disp('Analysing dielectric properties of tissues within reconstructed masks');
        % Load forward and inverse models stored in finite element object into
        % workspace.
        load([gui.dataPath 'results/' gui.finiteElement_flnm]);
         % Load the tissue masks and interfaces of the reference masks stored
        % in segmentation objects into workspace.
        load([gui.dataPath 'results/' gui.refSeg_flnm]);
        load([gui.dataPath 'results/' gui.recSeg_flnm]);       
        % Load the performance metrics object saved from the interface
        % analysis.       
        load([gui.dataPath 'results/' gui.performanceMetrics_flnm]);
        
        % Check to ensure there are mask and interface data extracted from 
        % forward model.
        fldListfwd = {'glandular', 'tumor'};
        gui = gui.checkTissueInterfacesFwdData(fldListfwd, refModel);
        if gui.errorFlag
            return
        end
        
        % Check ito ensure there are masks and interface data extracted from
        % the inverse model
        fldListInv = {'glandularRe', 'glandularIm', 'glandularMag', 'tumorRe', 'tumorIm', 'tumorMag'};        
        gui = gui.checkTissueInterfacesInvData(fldListInv, recModel);
        if gui.errorFlag
            return
        end
        
        % The segmentation objects hold the reference and reconstructed masks
        % in their data structures that contain geometric property information.
        % Therefore, metrics are applied to these regions to provide measures
        % that evaluate the accuracy with which the geometric properties of
        % the underlying structures are reconstructed. When applying the 
        % functions, the masks are evaluated within a region of interest
        % defined as the box that bounds the union between the reference
        % and reconstructed masks.
        % The analysis accommodates cases in which there is a requirement
        % for the reconstructed region to be compared with many small isolated
        % reference regions scattered throughout the interior. For this
        % scenario, the reconstructed mask is only compared with a reference
        % mask if the masks overlap.
        i = 1;
        for j =1:6
            if j < 4
                pa = pa. dielectricPropertyAnalysis(fe, fldListfwd{i}, fldListInv{j},  refModel, recModel, gui );
            else
                if j == 4
                    i = i + 1;
                end
                pa = pa.dielectricPropertyAnalysis(fe, fldListfwd{i}, fldListInv{j},  refModel, recModel, gui );
            end
        end 
        
        % Save results in pa object        
        gui.saveResults(pa);       
        gui.taskStatus = 'Analysis of dielectric properties of tissues within reconstructed masks complete';        
              
%==============================================================================================================
    case  'Output csv file of metric values for each region'
        disp('Creating csv file of metric values for each region');             
        % Load the performance metrics object saved from the geometric
        % and dielectric property analysis.        
        load([gui.dataPath 'results/' gui.performanceMetrics_flnm]);        
        
        fldListfwd = {'glandular', 'tumor'};        
        fldListInv = {'glandularRe', 'glandularIm', 'glandularMag', 'tumorRe', 'tumorIm', 'tumorMag'};        
        
        pa.getCSVfile(fldListfwd, fldListInv, gui);
        
        gui.taskStatus = 'CSV file of metric values for each region has been created and saved.';
end % case structure
%==============================================================================================================
end