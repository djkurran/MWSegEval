classdef segmentation 
    % Class segmentation
    % Methods segment forward and reconstructed images into tissue types
    % that leads to the decomposition of the breast interior into disjoint
    % tissue masks. The tissue masks are, in turn, mapped to tissue types that
    % collectively form an image. Edge points are uniformly sampled from the
    % boundary (or interface) of each reconstructed and reference mask.
    % Data structure holds the tissue masks, tissue types, and interfaces.
    %======================================================================    
    properties
        eprIm                 % Complex permittivity of immersion medium.
        n                     % number of rows in dielectric property map
        m                     % number of columns in dielectric property map
        xNodes                % engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
        yNodes                % engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels        
        masks                 % Tissue masks
        tissueTypes           % tissue types
        interfaces            % Tissue interfaces that segregate distinct tissue types
                              % interfaces.skinSurface - Coordinates of samples on skin surface. Used to generate geo file.
                              % interfaces.skinFat - Coordinates of samples on skin-fat interface.
                              % interfaces.fatGland Coordinates of samples on fibroglandular interface.                                      
        status                % Status of the segmentation technique
        bkgrd                 % Value assigned to the background for the clustering technique.
        nmCl                  % Number of clusters used for the iterative clustering technique
        h                     % test descision for the Null hypothesis that the data within a region over two successive iterations of the clustering algorithm are from the same continuous distrubution
        p                     % asymptotic p-value
        ks2stat               % test statistic of the Kolmogorov-Smirnov two-sample nonparametric hypothesis test
        ROIclusters           % clusters within a region-of-interest (ROI)
        clusters              % clusters formed from the k-means clustering technique that is iteratively applied to a ROI        
        cmpList               % components of reconstruction model to segment e.g. {'realEps', 'imEps', 'magEps'}
        ROIMask               % Region-of-interest mask
        notROI                % Complement of ROI (i.e., all pixels (or data) not in the ROI)
        ROIPart               % ROI partitioned into clusters by the clustering algorithm
        ROIClst               % Clusters of the ROI
        cmpItem               % It is a string variable. Component (or member) of the component list cmpList. For example realEps is a component of the cmpList,
        component             % It is an image. For example, Real(refModel.complexPermittivity) is the real component of the complex permittivity map.        
        iniClstCent           % The initial cluster centroids
        tumorMask             % tumor mask mapped after each iteration from cluster c_max(k).
        outsideTumorMask      % clusters ck with k = {2, 3, ..., max(k)-1} within the ROI are mapped to outsideTumorMask.
        outsideTumorTissue    % outsideTumorMask extracts tissues from reconstruction image to form outsideTumorTissue.
        outsideTumorDistribution % distrubution of pixels within outsideTumorTissue after the kth iteration.
        outsideTumorDistribution_prev % distrubution of pixels within outsideTumorTissue after the kth-1 iteration.
        ROItissueDistribution % distrubution of pixels within ROI after the kth iteration.
        tissueDistribution    
        tumorTissue           % tumorMask extracts tissues from reconstruction image to form tumorTissue.                 
        tumorDistribution     % distrubution of pixels within tumorTissue after the kth iteration.   
        tumorDistribution_prev % distrubution of pixels within tumorTissue after the kth-1 iteration.
        classifiedClusters    % clusters mapped to masks    
    end
%==========================================================================    
    methods
    %----------------------------------------------------------------------    
        %------------------------------------------------------------------
        function obj = segmentation(n, m, segmentationTechnique, xNodes, yNodes )
            % Constructor method to create instance of object. When instance
            % object is created, the method initializes some of the
            % properties. Additional properties are initiated if the
            % iterative unsupervised machine learning technique is used for
            % the segmentation process.
            % Inputs:
            %   n - (integer) number of rows in complex permittivity map
            %   m - (integer) number of columns in complex permittivity map
            %   segmentationTechnique - (string) - 'Yes' - iterative unsupervised machine learning technique used
            %   xNodes - (array of floats) engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   yNodes - (array of floats) engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels
            %--------------------------------------------------------------
            obj.n = n;
            obj.m = m;
            obj.xNodes = xNodes;
            obj.yNodes = yNodes;
            
            if segmentationTechnique
                obj.status = false;
                obj.bkgrd = -100;
                obj.nmCl = 3;                
                obj.h =1;
                obj.p = 0; 
                obj.ks2stat = 0;
                obj.ROIclusters = [];                
                obj.cmpList = {'realEps', 'imEps', 'magEps'};
            end
        end % Constructor segmentation
        %------------------------------------------------------------------
        function obj = segmentImmersionMedium(obj, complexPermittivity)
            % This function segments the immersion medium from the forward
            % model complex permittivity map. 
            % Inputs:
            %   complex permittivity - n by m array representing complex permittivity profile of forward model.
            % Outputs:
            %   obj.eprIm - float - complex permittivity of immersion medium
            %   obj.masks.interior - mask of breast interior including the skin layer
            %   obj.masks.immersionMedium - mask of immersion medium
            %   obj.interfaces.surface - interface of surface of breast            
            %--------------------------------------------------------------
            % Get the complex permittivity of the immersion medium
            % Since these are numerical models - assume that the immersion
            % medium is homogenous. For experimental models and patients, a
            % range of values would be used.
            % Also create mask of interior which is the region outside (i.e.,
            % the complement) of the immersion medium.
            obj.eprIm = complexPermittivity(1,1);
            obj.masks.interior = ones(obj.n, obj.m); 
            for i = 1:obj.n
                for j = 1:obj.m
                    if real(complexPermittivity(i,j)) == real(obj.eprIm)
                        obj.masks.interior(i,j) = 0;
                    end
                end
            end
            % Next get the edge points on the boundary of the interior
            % which corresponds to the interface between the surface of the
            % breast and the immersion medium.
            [Bg,Lg] = bwboundaries(obj.masks.interior);
            if ~isempty(Bg)
                iMax = 1;
                lMax = length(Bg{iMax});            
                for i = 2:length(Bg)
                    if length(Bg{i}) > lMax
                        iMax = i;
                        lMax = length(Bg{i});
                    end
                end
                boundary = Bg{iMax};
            end
            [r,c] = find(obj.masks.interior == 0);
            [in,on] = inpolygon(c,r,boundary(:,2),boundary(:,1));
            idx = find(in);
            for i = 1:length(idx)
                obj.masks.interior(r(idx(i)),c(idx(i))) = 1;
            end             
            obj.masks.immersionMedium = ones(obj.n,obj.m) - obj.masks.interior;
            obj.interfaces.surface = boundary;
        end % function segmentImmersionMedium        
        %------------------------------------------------------------------
        function obj = segmentSkinRegion(obj, eprArray, fwd_skin_epr_Hi_Limit, skinThkOffset)
            % This function segments the skin tissue from the forward model            
            % by eroding the interior mask by an amount determined by the user selected offset.
            % An option is also provided to segment only the skin region
            % using a threholding technique over a range of values. The
            % option is enables by setting segmentActualSkin to true.
            % Inputs:
            %   eprArray - n by m array of floats - real component of forward model complex permittivity map
            %   fwd_skin_epr_Hi_Limit - (float) relative permittivity of forward model skin tissue threshold value
            %   skinThkOffset - (float) - value set by user through GUI that determines offset from skin surface to position skin/fat interface.
            % Outputs:            
            %   obj.masks.skin - n by m binary image representing mask of skin
            %--------------------------------------------------------------           
            segmentActualSkin = false; % Set true if the actual skin tissue
            % is to be segmented, otherwise set to false if the skin region
            % is just offset from the outer skin surface.
            
            % initialize breast interior and skin mask
            interior = zeros(obj.n, obj.m);
            obj.masks.skin = zeros(obj.n, obj.m);
            % Compute the skin offset
            offset = (round(skinThkOffset*2)/2)*2;            
            % Erode the interior mask by the offset amount
            % The difference between the interior mask and the mask formed
            % by erroding the interior mask by the offset leads to a mask
            % that represents the skin region.
            se = strel('disk',offset);
            interior = imerode(obj.masks.interior, se);
            
            % The alternative approach to form the skin maak is to use the
            % thresholding techique. 
            if segmentActualSkin
                skinEpr = skin.*(eprArray);                
                for i = 1:obj.n
                    for j = 1:obj.m                    
                        if (skinEpr(i,j) > (fwd_skin_epr_Hi_Limit - (fwd_skin_epr_Hi_Limit*0.2))) &&   (skinEpr(i,j) < (fwd_skin_epr_Hi_Limit + (fwd_skin_epr_Hi_Limit*0.3)))
                            obj.masks.skin(i,j) = 1;
                        end
                    end
                end  
            else
                obj.masks.skin = obj.masks.interior -  interior;
            end
        end % function segmentSkinRegion
        %------------------------------------------------------------------
        function obj = segmentfwdFatRegion(obj, eprArray, fwd_fat_epr_Hi_Limit)            
            % This function segments the adipose tissue from the real component 
            % of the complex permittivity map of the forward map.            
            % Inputs:
            %   eprArray - n by m array of floats - real component of forward model complex permittivity map
            %   fwd_fat_epr_Hi_Limit (float) - threshold value used to segment the adipose tissue within the complex permittivity profile of the forward model 
            % Outputs:            
            %   obj.masks.fat - n by m binary image representing mask of adipose tissue of forward model
            %--------------------------------------------------------------   
            obj.masks.fat = zeros(obj.n, obj.m);              
            interior = obj.masks.interior -  obj.masks.skin;
            eprArray = interior.*(eprArray);
            [row, col] = find(0 < eprArray & eprArray < fwd_fat_epr_Hi_Limit);
            for i = 1:length(row)        
                obj.masks.fat(row(i),col(i)) = 1;
            end            
        end % function segmentfwdFatRegion
        %-----------------------------------------------------------------    
        function obj = segmentfwdTumorRegion(obj, eprArray , fwd_tumor_epr_Lo_Limit, fwd_tumor_contours_threshold, fwd_tumor_contours_maxNm )
            % This function segments the malignant tissue from the real component 
            % of the complex permittivity map of the forward map.            
            % Inputs:
            %   eprArray - n by m array - real component of forward model complex permittivity map
            %   fwd_tumor_epr_Lo_Limit - threshold value used to segment the malignant tissue within the complex permittivity profile of the forward model 
            %   fwd_tumor_contours_threshold - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   fwd_tumor_contours_maxNm - maximum number of contours to extract from the mask of the malignant tissue
            % Outputs:            
            %   obj.masks.tumor -  n by m binary image representing mask of malignant tissue of forward model
            %   obj.interfaces.tumor - structure - interfaces extracted from the tumor masks
            %--------------------------------------------------------------   
            obj.masks.tumor = zeros(obj.n, obj.m);                      
            interior = obj.masks.interior -   obj.masks.skin - obj.masks.fat;
            eprArray = interior.*(eprArray);
            % step 1 - identify all pixels that have same properties as
            % malignant tissue
            [row, col] = find(eprArray > fwd_tumor_epr_Lo_Limit);
            for i = 1:length(row)        
                 obj.masks.tumor(row(i),col(i)) = 1;
            end
            
            % step 2 - extract contours from malignant tissue regions            
            [Bt,Lt] = bwboundaries( obj.masks.tumor );
            if ~isempty(Bt)
                for i = 1:length(Bt)
                    bt(i) = size(Bt{i},1);
                end
                [bd, bdx] = sort(bt,'descend');
            end
            
            % step 3 - remove all countours that are less than 10 percent
            % of area of the largest region
            countourParm.contoursThrhld = fwd_tumor_contours_threshold;
            countourParm.contoursMaxNm = fwd_tumor_contours_maxNm;
            bnd = obj.processBoundary(bdx, Bt, countourParm);
            obj.interfaces.tumor = bnd;
            obj.masks.tumor = zeros(obj.n, obj.m);
            
            % step 4 - create tumor masks from the filtered contours
            for i = 1:length(bnd)
                tumor_i = roipoly(interior, bnd{i}.pts(:,2), bnd{i}.pts(:,1));
                obj.masks.tumor = obj.masks.tumor + tumor_i;
            end            
        end % function segmentfwdTumorRegion
        %------------------------------------------------------------------        
        function obj = segmentfwdGlandRegion(obj, eprArray , fwd_gland_epr_Lo_Limit)
            % This function segments the dense fibroglandular tissue from the real component 
            % of the complex permittivity map of the forward map.            
            % Inputs:
            %   eprArray - n by m array - real component of forward model complex permittivity map
            %   fwd_gland_epr_Lo_Limit - threshold value used to segment the dense fibroglandular tissue within the complex permittivity profile of the forward model             
            % Outputs:            
            %   obj.masks.fibroglandular - n by m binary image representing mask of fibroglandular tissue of forward model            
            %--------------------------------------------------------------   
            obj.masks.fibroglandular = zeros(obj.n, obj.m);             
            interior = obj.masks.interior - obj.masks.skin - obj.masks.fat - obj.masks.tumor;
            eprArray = interior.*(eprArray);
            [row, col] = find( eprArray > fwd_gland_epr_Lo_Limit);
            for i = 1:length(row)        
                obj.masks.fibroglandular(row(i),col(i)) = 1;
            end            
        end % function segmentfwdFatRegion                
        %------------------------------------------------------------------ 
        function obj = segmentfwdTransitionRegion(obj, fwd_gland_contours_threshold, fwd_gland_contours_maxNm)
            % This function segments the less dense transition tissue from the real component 
            % of the complex permittivity map of the forward map. The fibroglandular mask is 
            % combined with the transtion tissue mask to form the glandular
            % tissue mask. Edge points are extracted from the glandular
            % mask to form the interface of the glandular tissue region.
            % Inputs:
            %   eprArray - n by m array - real component of forward model complex permittivity map
            %   fwd_gland_contours_threshold - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   fwd_gland_contours_maxNm - maximum number of contours to extract from the mask of the glandular tissue
            % Outputs:            
            %   obj.masks.transition - n by m binary image representing mask of transition tissue of forward model            
            %   obj.masks.glandular - n by m binary image representing mask of glandular tissue of forward model
            %obj.interfaces.glandular - interface of glandular tissue of forward model            
            %--------------------------------------------------------------               
            
            % Step 1 - construct the transition masks
            obj.masks.transition = zeros(obj.n, obj.m);                  
            obj.masks.transition = obj.masks.interior - obj.masks.skin - obj.masks.fat - obj.masks.fibroglandular - obj.masks.tumor;
            
            % Step 2 - form the glandular mask as the union between the
            % fibroglandular and tumor masks.
            obj.masks.glandular = obj.masks.fibroglandular + obj.masks.transition;
            
            % step 3 - extract contours from glandular tissue region
            [Bt,Lt] = bwboundaries( obj.masks.glandular );
            if ~isempty(Bt)
                for i = 1:length(Bt)
                    bt(i) = size(Bt{i},1);
                end
                [bd, bdx] = sort(bt,'descend');
            end
            countourParm.contoursThrhld = fwd_gland_contours_threshold;
            countourParm.contoursMaxNm = fwd_gland_contours_maxNm;
            bnd = obj.processBoundary(bdx, Bt, countourParm);            
            obj.interfaces.glandular = bnd;
        end % function segmentfwdGlandRegion
        %------------------------------------------------------------------
        function obj = mapFwdMasksToTissueTypes(obj)
            % This function maps the tissue masks segmented from the
            % forward model to tissue types. The union of the tissue types
            % forms a tissue type map of the breast.
            % Inputs:
            %   obj.masks of forward model that include -
            %   masks.immersionMedium, masks.skin, masks.fat,
            %   masks.transition, masks.fibroglandular, masks.glandular,
            %   masks.tumor.
            % Outputs:
            %   obj.tissueTypes of forward model that include -
            %   tissueTypes.background, tissueTypes.fat,
            %   tissueTypes.transition , tissueTypes.fibroglandular ,
            %   tissueTypes.tumor.
            %--------------------------------------------------------------
            obj.tissueTypes.background = (obj.masks.immersionMedium + obj.masks.skin);  % background
            obj.tissueTypes.fat = obj.masks.fat*2;                                      % adipose tissue - 2
            obj.tissueTypes.transition = obj.masks.transition*3;                        % transition tissue - 3
            obj.tissueTypes.fibroglandular = obj.masks.fibroglandular*4;                % fibroglandular tissue - 4.
            obj.tissueTypes.tumor = obj.masks.tumor*5;                                  % fibroglandular tissue - 5.
            
        end % mapFwdMasksToTissueTypes
        %------------------------------------------------------------------        
        function obj = mapInvMasksToTissueTypes(obj)
            % This function maps the tissue masks segmented from each
            % component (real, Im, Mag) of the inverse model to tissue types. The union of the tissue types
            % forms a tissue type map of the breast.
            % Inputs:
            %   obj.masks of inverse model that include -
            %   masks.immersionMedium, masks.skin, masks.fatRe,
            %   masks.transitionRe, etc.
            % Outputs:
            %   obj.tissueTypes of each component of the inverse model that include -
            %   tissueTypes.background, tissueTypes.fatRe,
            %   tissueTypes.transitionRe , etc.
            %--------------------------------------------------------------
            
            obj.tissueTypes.background = (obj.masks.immersionMedium + obj.masks.skin);      % background
            obj.tissueTypes.fatRe = obj.masks.fatRe*2;                                      % adipose tissue - 2
            obj.tissueTypes.transitionRe = obj.masks.transitionRe*3;                        % transition tissue - 3
            obj.tissueTypes.fibroglandularRe = obj.masks.fibroglandularRe*4;                % fibroglandular tissue - 4.
            obj.tissueTypes.tumorRe = obj.masks.tumorRe*5;                                  % fibroglandular tissue - 5.
            
            obj.tissueTypes.fatIm = obj.masks.fatIm*2;                                      % adipose tissue - 2
            obj.tissueTypes.transitionIm = obj.masks.transitionIm*3;                        % transition tissue - 3
            obj.tissueTypes.fibroglandularIm = obj.masks.fibroglandularIm*4;                % fibroglandular tissue - 4.
            obj.tissueTypes.tumorIm = obj.masks.tumorIm*5;                                  % fibroglandular tissue - 5.
            
            obj.tissueTypes.fatMag = obj.masks.fatMag*2;                                    % adipose tissue - 2
            obj.tissueTypes.transitionMag = obj.masks.transitionMag*3;                      % transition tissue - 3
            obj.tissueTypes.fibroglandularMag = obj.masks.fibroglandularMag*4;              % fibroglandular tissue - 4.
            obj.tissueTypes.tumorMag = obj.masks.tumorMag*5;                                      % fibroglandular tissue - 5.
            
        end % mapInvMasksToTissueTypes
        %------------------------------------------------------------------
        function obj = segmentInvImmersionMedium(obj, fwdMdlImmersionMedium, fwdInterior)
            % This function segments the immersion medium from the inverse
            % model complex permittivity map. It is assumed that the toolbox 
            % is applied to numerical models, so the immersion medium mask
            % segmented with the forward model is also used for the inverse
            % model. For an experimental or patient study system, this
            % method to be replaced with a threholding technique applied over
            % a range of values.           
            % Inputs:
            %   fwdMdlImmersionMedium - n by m array representing mask for immersion medium of the forward model.
            %   fwdInterior - n by m array representing mask of interior of forward model.
            % Outputs:           
            %   obj.masks.interior - mask of breast interior including the skin layer
            %   obj.masks.immersionMedium - mask of immersion medium                 
            %--------------------------------------------------------------   
            obj.masks.immersionMedium = fwdMdlImmersionMedium;
            obj.masks.interior = fwdInterior;            
        end % function segmentInvImmersionMedium
        %------------------------------------------------------------------
        function obj = segmentInvSkinRegion(obj, fwdMdlSkin)
            % This function segments the skin region from the inverse
            % model complex permittivity map. It is assumed that the toolbox 
            % is applied to numerical models, so the skin region mask
            % segmented with the forward model is also used for the inverse
            % model. For an experimental or patient study system, this
            % method to be replaced with a threholding technique applied over
            % a range of values.                                           
            % Inputs:
            %   fwdMdlSkin - n by m array - mask of skin region of forward model.            
            % Outputs:            
            %   obj.masks.skin - mask of skin region of inverse model
            %--------------------------------------------------------------                   
            obj.masks.skin = fwdMdlSkin;                      
        end % function segmentInvSkinRegion                
        %-----------------------------------------------------------------                
        function obj = segmentinvFatRegion_Threshold(obj, complexPermittivity, segmentationThresholdIn)
            % This function segments the adipose tissue from the real, imaginary, and magnitude
            % of the complex permittivity map of the inverse model using the thresholding technique.            
            % Inputs:
            %   complexPermittivity - n by m array - inverse model complex permittivity map
            %   segmentationThresholdIn - structure array (plane organization)
            %   segmentationThresholdIn.FEMCSIrealAdipose - threshold value to segment the adipose tissue from the real component 
            %   segmentationThresholdIn.FEMCSIimAdipose - threshold value to segment the adipose tissue from the imaginary component 
            %   segmentationThresholdIn.FEMCSImagAdipose - threshold value to segment the adipose tissue from the magnitude of complex permittivity
            % Outputs:            
            %   obj.masks.fatRe, obj.masks.fatIm, obj.masks.fatMag - mask of adipose tissue of forward model
            %-------------------------------------------------------------- 
            % Define the ROI to extract tissue
            interior = obj.masks.interior - obj.masks.skin;
            
            % segment real part of the ROI of inverse model
            obj.masks.fatRe = zeros(obj.n, obj.m);  
            realArray = interior.*real(complexPermittivity);            
            [row, col] = find(0 < realArray & realArray < segmentationThresholdIn.FEMCSIrealAdipose);
            for i = 1:length(row)        
                obj.masks.fatRe(row(i),col(i)) = 1;
            end
            
            % segment imaginary part of inverse model
            obj.masks.fatIm = zeros(obj.n, obj.m);                       
            imagArray = interior.*-imag(complexPermittivity);            
            [row, col] = find(0 < imagArray & imagArray < segmentationThresholdIn.FEMCSIimAdipose);
            for i = 1:length(row)        
                obj.masks.fatIm(row(i),col(i)) = 1;
            end
            
            % segment magnitude of inverse model
            obj.masks.fatMag = zeros(obj.n, obj.m);  
            magArray = interior.*abs(complexPermittivity);            
            [row, col] = find(0 < magArray & magArray < segmentationThresholdIn.FEMCSImagAdipose);
            for i = 1:length(row)        
                obj.masks.fatMag(row(i),col(i)) = 1;
            end
        end % function segmentinvFatRegion_Threshold
        %------------------------------------------------------------------
        function obj = segmentinvTumorRegion_Threshold(obj, complexPermittivity, segmentationThresholdIn, gui)
            % This function segments the malignant tissue from the real, imaginary, and magnitude
            % of the complex permittivity map of the inverse model using the thresholding technique.            
            % Inputs:
            %   complexPermittivity - n by m array - inverse model complex permittivity map
            %   segmentationThresholdIn.FEMCSITumorMaxThreshold - threshold value to segment the malignant tissue from the real, imaginary, and magnitude. Expressed as percent of maximum value detected within a ROI.
            %   gui.contoursThrhld - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   gui.contoursMaxNm - maximum number of contours to extract from the mask of the glandular tissue            
            % Outputs:            
            %   obj.masks.tumorRe, obj.masks.tumorIm, obj.masks.tumorMag - mask of malignant tissue of inverse model
            %   obj.interfaces.tumorRe, obj.interfaces.tumorIm, obj.interfaces.tumorMag - interface of malignant tissue of inverse model
            %-------------------------------------------------------------- 
            
            % The dense tissue mask is comprised of transition,
            % fibroglandular, and malignant tissues
            denseTissueMask = zeros(obj.n, obj.m);
            countourParm.contoursThrhld = gui.contoursThrhld;
            countourParm.contoursMaxNm = gui.contoursMaxNm;
            
            % First lets get the maximum value of the real part of the complex permittivity  within
            % the dense tissue structure that includes both the tumor and
            % fibroglandular tissues
            obj.masks.tumorRe = zeros(obj.n, obj.m);
            denseTissueMask = ones(obj.n, obj.m) - obj.masks.immersionMedium -  obj.masks.skin - obj.masks.fatRe;
            maxReal = max(max(denseTissueMask.*real(complexPermittivity)));            
            interior =  obj.masks.interior - obj.masks.skin -  obj.masks.fatRe;
            
            % segment real component of complex permittivity into tumor
            % tissue mask            
            realArray = interior.*real(complexPermittivity);             
            [row, col] = find(realArray > maxReal*0.01*segmentationThresholdIn.FEMCSITumorMaxThreshold);
            for i = 1:length(row)        
                obj.masks.tumorRe(row(i),col(i)) = 1;
            end
            
            % Extract contours from malignant tissue regions
            % remove all countours that are less than 10 percent
            % of area of the largest region, allow a maximum of 20 contours
            % to be extracted.
            
            bnd = obj.getTissueInterfaces(obj.masks.tumorRe, countourParm);
            obj.interfaces.tumorRe = bnd;
            
            % update masks
            obj.masks.tumorRe = obj.updateMasks(interior, bnd);
            
            % segment imaginary component of complex permittivity into tumor
            % tissue mask
            obj.masks.tumorIm = zeros(obj.n, obj.m);
            
            % Lets also get the maximum value of the imaginary part of the complex permittivity  within
            % the dense tissue structure that includes both the tumor and
            % fibroglandular tissues
            denseTissueMask = ones(obj.n, obj.m) - obj.masks.immersionMedium - obj.masks.skin - obj.masks.fatIm;
            maxIm = max(max(denseTissueMask.*-imag(complexPermittivity)));
            interior =  obj.masks.interior - obj.masks.skin - obj.masks.fatIm;
            
            imagArray = interior.*-imag(complexPermittivity);
            [row, col] = find(imagArray > maxIm*0.01*segmentationThresholdIn.FEMCSITumorMaxThreshold);
            for i = 1:length(row)        
                obj.masks.tumorIm(row(i),col(i)) = 1;
            end
            
            % Extract contours from malignant tissue regions
            % remove all countours that are less than 10 percent
            % of area of the largest region, allow a maximum of 20 contours
            % to be extracted.            
            bnd = obj.getTissueInterfaces(obj.masks.tumorIm, countourParm);
            obj.interfaces.tumorIm = bnd;
            
            % update masks
            obj.masks.tumorIm = obj.updateMasks(interior, bnd); 
            
            % Finally, lets get the maximum value of the magnitude of the complex permittivity  within
            % the dense tissue structure that includes both the tumor and
            % fibroglandular tissues
            obj.masks.tumorMag = zeros(obj.n, obj.m);
            
            % Lets also get the maximum value of the magnitude of the complex permittivity  within
            % the dense tissue structure that includes both the tumor and
            % fibroglandular tissues
            denseTissueMask = ones(obj.n, obj.m) - obj.masks.immersionMedium - obj.masks.skin - obj.masks.fatMag;
            maxMag = max(max(denseTissueMask.*abs(complexPermittivity)));
            interior =  obj.masks.interior - obj.masks.skin -  obj.masks.fatMag;
            
            magArray = interior.*abs(complexPermittivity);
            [row, col] = find(magArray > maxMag*0.01*segmentationThresholdIn.FEMCSITumorMaxThreshold);
            for i = 1:length(row)        
                obj.masks.tumorMag(row(i),col(i)) = 1;
            end
            
            % Extract contours from malignant tissue regions
            % remove all countours that are less than 10 percent
            % of area of the largest region, allow a maximum of 20 contours
            % to be extracted.            
            bnd = obj.getTissueInterfaces(obj.masks.tumorMag, countourParm);
            obj.interfaces.tumorMag = bnd;
            
            % update masks
            obj.masks.tumorMag = obj.updateMasks(interior, bnd); 
            
        end % function segmentinvTumorRegion_Threshold
        %------------------------------------------------------------------                         
        function obj = segmentinvFibroGlandRegion_Threshold(obj, complexPermittivity, segmentationThresholdIn)            
            % This function segments the fibroglandular tissue from the real, imaginary, and magnitude
            % of the complex permittivity map of the inverse model using the thresholding technique.            
            % Inputs:
            %   complexPermittivity - n by m array - inverse model complex permittivity map
            %   segmentationThresholdIn.FEMCSIrealTransition - threshold value to segment the fibroglandular tissue from the real component
            %   segmentationThresholdIn.FEMCSIimTransition - threshold value to segment the fibroglandular tissue from the imaginary component
            %   segmentationThresholdIn.FEMCSImagTransition - threshold value to segment the fibroglandular tissue from the magnitude of the complex permittivity map            
            % Outputs:            
            %   obj.masks.fibroglandularRe, obj.masks.fibroglandularIm, obj.masks.fibroglandularMag - mask of fibroglandular tissue of inverse model            
            %-------------------------------------------------------------- 
            obj.masks.fibroglandularRe = zeros(obj.n, obj.m);             
            interior = obj.masks.interior -   obj.masks.skin - obj.masks.fatRe - obj.masks.tumorRe;
            realArray = interior.*real(complexPermittivity);
            [row, col] = find(realArray > segmentationThresholdIn.FEMCSIrealTransition);
            for i = 1:length(row)        
                obj.masks.fibroglandularRe(row(i),col(i)) = 1;
            end            
            
            % segment imaginary part of FEM-CSI
            obj.masks.fibroglandularIm = zeros(obj.n, obj.m); 
            interior = obj.masks.interior -   obj.masks.skin - obj.masks.fatIm - obj.masks.tumorIm;
            imagArray = interior.*-imag(complexPermittivity);            
            [row, col] = find(imagArray > segmentationThresholdIn.FEMCSIimTransition);
            for i = 1:length(row)        
                obj.masks.fibroglandularIm(row(i),col(i)) = 1;
            end
            
            % segment magnitude of FEM-CSI
            obj.masks.fibroglandularMag = zeros(obj.n, obj.m); 
            interior = obj.masks.interior - obj.masks.skin - obj.masks.fatMag - obj.masks.tumorMag;
            magArray = interior.*abs(complexPermittivity);            
            [row, col] = find(magArray > segmentationThresholdIn.FEMCSImagTransition);
            for i = 1:length(row)        
                obj.masks.fibroglandularMag(row(i),col(i)) = 1;
            end
        end % function segmentFEM_CSIFibroGlandRegion_Threshold        
        %------------------------------------------------------------------ 
        function obj = segmentinvTransitionRegion_Threshold(obj, contoursThrhld, contoursMaxNm)
            % This function segments the transition tissue from the real, imaginary, and magnitude
            % of the complex permittivity map of the inverse model.
            % The transition and fibroglandular masks are combined to form
            % a single dense tissue mask for the glandular tissue.
            % Inputs:
            %   complexPermittivity - n by m array - inverse model complex permittivity map
            %   segmentationThresholdIn.FEMCSIrealTransition - threshold value to segment the transition malignant tissue from the real component
            %   segmentationThresholdIn.FEMCSIimTransition - threshold value to segment the transition malignant tissue from the imaginary component
            %   segmentationThresholdIn.FEMCSImagTransition - threshold value to segment the transition malignant tissue from the magnitude of the complex permittivity map            
            %   contoursThrhld - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   contoursMaxNm - maximum number of contours to extract from the mask of the glandular tissue            
            % Outputs:            
            %   obj.masks.transitionRe, obj.masks.transitionIm, obj.masks.transitionMag - mask of transition tissue of inverse model            
            %   obj.masks.fibroglandularRe, obj.masks.fibroglandularIm, obj.masks.fibroglandularMag - mask of glandular tissue (i.e., transition and fibroglandular) of inverse model            
            %   obj.interfaces.glandularRe, obj.interfaces.glandularIm, obj.interfaces.glandularMag - interfaces of the glandular structure  
            %--------------------------------------------------------------            
            countourParm.contoursThrhld = contoursThrhld;
            countourParm.contoursMaxNm = contoursMaxNm;            
            % segment real part of FEM-CSI
            obj.masks.transitionRe = zeros(obj.n, obj.m);                  
            obj.masks.transitionRe = obj.masks.interior - obj.masks.skin - obj.masks.fatRe - obj.masks.fibroglandularRe - obj.masks.tumorRe;
            obj.masks.glandularRe = obj.masks.fibroglandularRe + obj.masks.transitionRe;
            % Extract interfaces from masks, filtering small regions or artifacts 
            % from original masks - filtered regions become part of fat region.
            obj.interfaces.glandularRe = obj.getTissueInterfaces(obj.masks.glandularRe, countourParm);
            
            % segment imaginary part of inverse model
            obj.masks.transitionIm = zeros(obj.n, obj.m);                  
            obj.masks.transitionIm = obj.masks.interior - obj.masks.skin - obj.masks.fatIm - obj.masks.fibroglandularIm - obj.masks.tumorIm;
            obj.masks.glandularIm = obj.masks.fibroglandularIm + obj.masks.transitionIm;
            % Extract interfaces from masks, filtering small regions or artifacts 
            % from original masks - filtered regions become part of fat region.
            obj.interfaces.glandularIm = obj.getTissueInterfaces(obj.masks.glandularIm, countourParm);
            
            % segment magnitude of inverse model
            obj.masks.transitionMag = zeros(obj.n, obj.m);                  
            obj.masks.transitionMag = obj.masks.interior - obj.masks.skin - obj.masks.fatMag - obj.masks.fibroglandularMag - obj.masks.tumorMag;
            obj.masks.glandularMag = obj.masks.fibroglandularMag + obj.masks.transitionMag;
            % Extract interfaces from masks, filtering small regions or artifacts 
            % from original masks - filtered regions become part of fat region.
            obj.interfaces.glandularMag = obj.getTissueInterfaces(obj.masks.glandularMag, countourParm);
        end % function segmentFEM_CSIGlandRegion_Threshold
        %------------------------------------------------------------------        
        function obj = getROI(obj, ROIMask)
            % This function sets the ROIMask object property to the ROIMask
            % passed to the function.
            % Inputs:
            %   ROIMask - n by m array corresponding to a ROI mask.
            % Outputs:
            %   obj.ROIMask - n by m array ROI mask used to extract tissues
            %                 witin a region of interest to segment.
            %--------------------------------------------------------------
            obj.ROIMask = ROIMask;
        end
        %------------------------------------------------------------------
        function obj = setBackbackground(obj)
            % This function sets the region outside the region of interest
            % to some integeter value determined by obj.bkgrd.
            % Inputs:
            %   obj.ROIMask - n by m array corresponding to a ROI mask.
            %   obj.bkgrd - integer - value used to identify the region
            %   outside the region of interest which is defined as
            %   background.
            % Outputs:
            %   obj.notROI - n by m array corresponding to background which
            %               is the region outside the region of interest
            %--------------------------------------------------------------
             obj.notROI = zeros(size(obj.ROIMask,1),size(obj.ROIMask,2));
             i = find(obj.ROIMask == 0);
             obj.notROI(i) = obj.bkgrd;    
        end % setBackbackground
        %------------------------------------------------------------------
        function obj = getImageFromList(obj, idx, complexPermittivity) 
            % This function gets the next reconstructed image to segment (real(eps), img(eps),
            % real(Chi), img(Chi), mag(eps), mag(Chi))            
            % Inputs:
            %   idx - integer - the index of the component in the image
            %         array
            %   complexPermittivity - n by m array - inverse complex permittivity profile
            % Outputs:            
            %   obj.cmpItem - string - describing component of inverse model to be segmented.
            %--------------------------------------------------------------   
            obj.cmpItem = obj.cmpList{idx};
            switch obj.cmpItem
                case 'realEps'
                    obj.component = real(complexPermittivity);
                case 'imEps'
                    obj.component = -imag(complexPermittivity);                                     
                case 'magEps'
                    obj.component = abs(complexPermittivity);                
            end
        end % getComponentFromList
        %------------------------------------------------------------------
        function obj = iniNmCl(obj)
            % This function initializes the number of clusters k is initialized to three,
            % and the k-means++ algorithm is used to initialize k pixel intensities
            % as cluster centroids.         
            % Inputs:
            %   none
            % Outputs:            
            %   obj.nmCl - integer - number of clusters to segment the region of interest.
            %--------------------------------------------------------------
            obj.nmCl = 3;
        end
        %------------------------------------------------------------------
        function obj = iniHypothesis(obj)
            % This function initializes the Kolmogorov-Smirnov two sample nonparametric
            % hypothesis test by rejecting the null hypothesis that the distribution over each of the
            % groups of clusters has not significantly changed (statistically) after an iteration
            % of the algorithm.
            % Inputs:
            %   none
            % Outputs:            
            %   obj.h - integer - Kolmogorov-Smirnov two sample nonparametric
            %           hypothesis test result that infers that the distribution over each of the
            %           groups of clusters has or has not significantly changed (statistically)
            %           after an iteration of the algorithm 
            %--------------------------------------------------------------
            obj.h = 1;            
        end
        %------------------------------------------------------------------        
        function obj = getInitialClusterCentroids(obj)
            % This function gets the centroids of each of the clusters
            % outside the cluster with the highest value (corresponding)
            % the tissue that has the highest density) when the clusters are 
            % initialized.
            % Inputs:
            %   none
            % Outputs:            
            %   obj.iniClstCent - centroids of the clusters
            %--------------------------------------------------------------            
            [r,c,vg] = find(obj.outsideTumorMask.*obj.component);
            mxVg = max(vg);
            mnVg = min(vg);
            if obj.nmCl == 3
                obj.iniClstCent = [obj.bkgrd mnVg mxVg]';
                return;
            end           
            delta =  (mxVg - mnVg)/(obj.nmCl - 2);
            obj.iniClstCent = mnVg:delta:mxVg;
            obj.iniClstCent = [obj.bkgrd obj.iniClstCent]';        
        end                
        %------------------------------------------------------------------
        function obj = getOutsideTumorMask(obj)
            % This function gets the mask for the region outside the tumor
            % estimated region.
            % Inputs:
            %   none
            % Outputs:            
            %   obj.outsideTumorMask - n by m array corresponding to the 
            %           mask of region outside tumor region.
            %--------------------------------------------------------------       
            if obj.nmCl == 3
                obj.tumorMask = [];                
                obj.outsideTumorMask = obj.ROIMask;
            else
                obj.outsideTumorMask = obj.ROIMask - obj.tumorMask;
            end
        end % getGlandMask
        %------------------------------------------------------------------
        function obj = getROITissueDist(obj)
            % This function gets the tissue distribution within the region
            % of interest corresponding to malignant tissue.
            % Inputs:
            %   none
            % Outputs:            
            %   obj.ROItissueDistribution - distribution of the dielectric
            %               properties within the tumor region
            %--------------------------------------------------------------
            [r,c,obj.ROItissueDistribution] = find(obj.ROIMask.*obj.component);            
        end % getTumorTissue
        %------------------------------------------------------------------
        function obj = partitionROI(obj)
            % This function uses the k-means clustering technique to
            % parittion the region of interest into k clusters.
            % Inputs:
            %   none
            % Outputs:
            %   obj.ROIclusters{obj.nmCl}.clusterMap - region of interest
            %                               that has been segmented into clusters.
            %   obj.ROIclusters{obj.nmCl}.labels - (integer) corresponds to
            %                               label identifying each cluster.      
            %   obj.ROIclusters{obj.nmCl}.centroids - centroids of each of
            %                               the clusters.                  
            %--------------------------------------------------------------
            % First extract the tissue within the region of interest. Set
            % the region outside the region of interest as background.
            % The background is the immersion medium and skin; the ROI is
            % the breast interior comprised of dense tissue (malginant,
            % fibroglandular), medium dense (transition tissue), and low
            % dense tissue (adipose). 
            ROI = (obj.ROIMask.*obj.component) + obj.notROI;            
            [m, n] = size(ROI);
            % Schedule number of replicates that kmeans uses based on the
            % number of clusters to partitions the images. The schedule was
            % formulated based on an analysis of variation of the cost
            % function values as the number of clusters increases. The
            % analysis was carried out on various numerical breast models
            % with a variety of tissue distrubtions and densities. The
            % replicate parameter specifies the number of times to repeat 
            % clustering using new initial cluster centroud positions. The
            % algorithm initializes the replications separately using k-means++.
            if obj.nmCl < 6
                repl = 10;
            elseif (obj.nmCl > 5 & obj.nmCl < 8)
                repl = 30;
            elseif (obj.nmCl > 7 & obj.nmCl < 9)
                repl = 50;
            else
                repl = 100;
            end            
            % If the cluster losses all of its member observations, a new
            % cluster of the one point furthest from its centroid is
            % chosen (i.e., 'EmptyAction' -> 'singleton'). Kmeans uses the
            % squared Euclidean (i.e., L2 norm) distance for minimizatoin.
            % That is, each centroid is the mean of the points of that
            % cluster. The kmeans performs an online update phase in addition
            % to a batch update phase. Although the online phase can be
            % time-consuming, it guarentees a solution that is a local
            % minimum of the distance criterion. That is, the algorithm
            % finds a partition of the data in which moving any single
            % point to a different cluster increases the total sum of
            % distances. 
                        
            % Apply the k-means clustering algorithm to the region of
            % interest.
            [idx,C]  = kmeans(ROI(:), obj.nmCl, 'EmptyAction', 'singleton','Distance','sqeuclidean','Replicates',repl,'OnlinePhase','on');
            obj.ROIPart = reshape(idx,[m n]);
            obj.ROIClst = obj.ROIPart;
            D = sort(C);
            for i = 1:length(D)
                m = find(D(i) == C);
                obj.ROIClst(find(obj.ROIPart == m)) = i;
            end
            obj.ROIclusters{obj.nmCl}.centroids = D;
            obj.ROIclusters{obj.nmCl}.labels =  unique(idx);
            obj.ROIclusters{obj.nmCl}.clusterMap = obj.ROIClst;
        end % function partitionROI           
        %-----------------------------------------------------------------
        function imageClusters( obj, xNodes, yNodes)
            % Helper function used when troubleshooting to image the clusters
            % Inputs:
            %   xNodes - (array of floats) engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   yNodes - (array of floats) engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels
            %--------------------------------------------------------------            
            figure;
            imagesc(imrotate(obj.ROIClst,90),'XData',xNodes,'YData',-yNodes);
            ylabel('y (m)','FontSize',14);
            xlabel('x (m)','FontSize',14);
            xlim([xNodes(1) xNodes(end)]);
            ylim([yNodes(1) yNodes(end)]);             
            grid on;
            set(gca,'YDir','normal');
            set(gca,'fontsize',14);
            cmap = jet(obj.nmCl);
            t = 1:1:obj.nmCl;
            t1 = [];
            for i = 1:obj.nmCl                
                t1{i} = num2str(t(i));
            end
            cmap(1,:) = [1 1 1];
            cmap(2,:) = [0 0 1];
            cmap(length(cmap),:) = [0 0 0];
            colormap(cmap);            
            colorbar('Ticks', t,'TickLabels', t1);   
        end % function imageClusters          
        %------------------------------------------------------------------
        function obj = getTumorMask(obj)
            % This function gets the region within the ROI segmented as malignant
            % This function gets the region within the ROI segmented as malignant
            % tissue. In this case, the malignant tissue has the highest density, so
            % is associated with the cluster with the highest valued cluster.
            % Inputs:
            %   none
            % Outputs:
            %   obj.tumorMask - highest valued cluter region.               
            %--------------------------------------------------------------
            obj.tumorMask = zeros(size(obj.component,1),size(obj.component,2));
            i = find(obj.ROIClst == obj.nmCl); 
            obj.tumorMask(i) = 1;
        end % getTumorMask
        %------------------------------------------------------------------        
        function obj = getTumorTissue(obj)
            % This function extracts the tissues from the component image
            % using the tumor mask constructed from the segmentation step.
            % Inputs:
            %   none
            % Outputs:
            %   obj.tumorTissue - tissues within the tumor mask segmented from the component image.               
            %--------------------------------------------------------------            
            obj.tumorTissue = obj.tumorMask.*obj.component; 
        end % getTumorTissue
        %------------------------------------------------------------------        
        function obj = getTumorTissueDistribution(obj)
            % This function gets the distribution of the tissues within the
            % tumor mask segmented in the segmentation step.
            % Inputs:
            %   none
            % Outputs:
            %   obj.tumorDistribution_prev - tissues within the tumor mask segmented in the previous iteration.
            %   obj.tumorTissue - tissues within the tumor mask segmented in the present iteration.
            %--------------------------------------------------------------        
            if obj.nmCl > 3
                obj.tumorDistribution_prev = obj.tumorDistribution;
            end                
            [r,c,obj.tumorDistribution] = find(obj.tumorTissue);
        end % getTumorTissueDistribution
        %------------------------------------------------------------------
        function obj = getTissueOutsideTumorRegion(obj)
            % This function extracts the tissues from the component image
            % using the mask for the region outside the tumor constructed from
            % the segmentation step.
            % Inputs:
            %   none
            % Outputs:
            %   obj.outsideTumorTissue - tissues within the mask outside the
            %               tumor region segmented from the component image.               
            %--------------------------------------------------------------     
            obj.outsideTumorTissue =  obj.outsideTumorMask.*obj.component;
        end % getTissueOutsideTumorRegion  
        %------------------------------------------------------------------                
        function obj = getOutsideTumorTissueDistribution(obj)
            % This function gets the distribution of the tissues within the
            % mask outside the tumor segmented in the segmentation step.
            % Inputs:
            %   none
            % Outputs:
            %   obj.outsideTumorDistribution_prev - tissues within mask outside tumor region segmented in the previous iteration.
            %   obj.outsideTumorDistribution - tissues within mask outside tumor region segmented in the present iteration.
            %--------------------------------------------------------------        
            if obj.nmCl > 3
                obj.outsideTumorDistribution_prev = obj.outsideTumorDistribution;
            end 
             [r,c,obj.outsideTumorDistribution] = find(obj.outsideTumorTissue);
        end % getGldTissueDist                
        %------------------------------------------------------------------        
        function obj = testHypothesis(obj)
            % This function, after each iteration, compares the probability
            % distribution of the pixels within the largest valued cluster
            % of the present iteration with the probability distribution of
            % the pixels over the largest valued cluster of the previous
            % iteration. It also compares the probability distribution of the
            % pixels outside the largest cluster with the probability distribution of
            % the pixels outside the largest cluster in the previous iteration. The 
            % algorithm terminates when a statistical test (the Kolmogorov-Smirnov two
            % sample nonparametric hypothesis test) infers that the distribution
            % over each of the groups of clusters has not significantly changed
            % (statistically) after an iteration of the algorithm. .
            % Inputs:
            %   none
            % Outputs:
            % Outputs:            
            %   obj.h - integer - Kolmogorov-Smirnov two sample nonparametric
            %           hypothesis test result that infers that the distribution over each of the
            %           groups of clusters has or has not significantly changed (statistically)
            %           after an iteration
            %--------------------------------------------------------------        
            % If k > 3, then test the hypothesis that the dielectric properties
            % (complex permittivity or contrast) within the ROI outside the tumor
            % region of the present iteration and the dielectric properties within
            % the ROI outside the tumor region of the previous iteration are from
            % populations that follow the same distribution at the 5
            % percent significance level
            if obj.nmCl > 3
                [ho, obj.p, obj.ks2stat] = kstest2(obj.outsideTumorDistribution, obj.outsideTumorDistribution_prev, 'Alpha',0.99999999);
                [ht, pt, ks2statTmr] = kstest2(obj.tumorDistribution, obj.tumorDistribution_prev, 'Alpha',0.9999999999);                 
                obj.h = ho && ht; % ho - hypothesis test result related to region outiside tumor
                                  % ht - hypothesis test result related to tumor region
            end               
        end % testHypothesis       
        %------------------------------------------------------------------
        function obj = incrementNmCl(obj)
            % This function increments the number of clusters to segment the
            % region of interest.
            % Inputs:
            %   none
            % Outputs:
            %   obj.nmCl - number of clusters to segment the region of interst            
            %--------------------------------------------------------------    
            if obj.h == 1
                obj.nmCl = obj.nmCl + 1;            
            end                
        end
        %------------------------------------------------------------------    
        function obj = rejectTumorRegionArtifacts(obj, component, gui)
            % This function used to sense and then reject artifacts
            % that are incorrectly attributed to malignant tissue.            
            % Inputs:
            %   component - which component of the inverse model is being
            %           segmented (i.e., real, imaginary, or magnitude.
            %   gui.contoursThrhld - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   gui.contoursMaxNm - maximum number of contours to extract from the mask of the malignant tissue            
            % Outputs:
            %   obj.masks.tumorRe, obj.masks.tumorIm, obj.masks.tumorMag - n by m mask of segmented malignant tissue segmented from real, imaginary, or magnitude component.
            %   obj.interfaces.tumorRe, obj.interfaces.tumorIm, obj.interfaces.tumorMag - interfaces extracted from the tumor masks of each of the components (i.e., real, imaginary, magnitude)
            %--------------------------------------------------------------  
            countourParm.contoursThrhld = gui.contoursThrhld;
            countourParm.contoursMaxNm = gui.contoursMaxNm;            
            [n,m] = size(obj.ROIMask);
            switch component
                case 'realEps'
                    % Store the tumor mask, before processing
                    tumorMasks_prev = obj.masks.tumorRe;
                    
                    % Get all of the edge points for each region of
                    % malignant tissue, and sort the interfaces in terms of
                    % largest area bound by the interface to smallest area
                    % bound by the interface.
                    [Bt,Lt] = bwboundaries(obj.masks.tumorRe);
                    if ~isempty(Bt)
                        for i = 1:length(Bt)
                            bt(i) = size(Bt{i},1);
                        end
                        [bd, bdx] = sort(bt,'descend');
                    end
                    
                    % process the interfaces by first limiting the number
                    % of boundaries to contoursMaxNm.
                    % For those contoursMaxNm interfaces, delete all those
                    % interfaces that bound an area that is less than
                    % contoursThrhld percent of the largest area bound. For
                    % example, if contoursMaxNm = 7, then only select 7
                    % interfaces and delete the remaining interfaces. Of
                    % the remaining interfaces, if contoursThrhld = 10, then
                    % delete all those interfaces that bound an area
                    % less than 10 % of the largest area bound.            
                    bnd = obj.processBoundary(bdx, Bt, countourParm);
                    % Save all of interfacess in the interfaces structure.
                    obj.interfaces.tumorRe = bnd;
                    
                    % Identify all of the interfaces that were rejected.
                    % For example, if if contoursMaxNm = 7, then identify
                    % eighth, ninth, etc. interfaces. Moreover, identify 
                    % the remaining interfaces that are less than
                    % contoursThrhld (10 percent, for example) less than
                    % the maximum area. Once identified, create a mask of
                    % these interfaces and name them 'artifact'. Remove the
                    % 'artifact' mask from the tumor mask tumorMasks_prev
                    % and save as tumor mask.
                    bndArt = obj.identifyArtifacts(bdx, Bt, countourParm);
                    
                    if ~isempty(bndArt)
                        artifact = zeros(obj.n, obj.m);
                        for i = 1:length(bndArt)
                            artifact_i = roipoly(obj.ROIMask, bndArt{i}.pts(:,2), bndArt{i}.pts(:,1));
                            artifact = artifact + artifact_i;
                        end
                        obj.masks.tumorRe = tumorMasks_prev - artifact;
                        
                        % Errors occur when creating masks from interfaces,
                        % so remove these small artifaces with
                        % 'bwareaopen'.
                        obj.masks.tumorRe = bwareaopen(obj.masks.tumorRe, 30);
                        % Add the artifacts to the fibroglandular region.
                        obj.masks.fibroglandularRe = obj.masks.fibroglandularRe + (tumorMasks_prev - obj.masks.tumorRe);
                        obj.masks.glandularRe = obj.masks.fibroglandularRe + obj.masks.transitionRe;
                    end
                    
                case 'imEps'
                    % Store the tumor mask, before processing
                    tumorMasks_prev = obj.masks.tumorIm;
                    
                    % Get all of the edge points for each region of
                    % malignant tissue, and sort the interfaces in terms of
                    % largest area bound by the interface to smallest area
                    % bound by the interface.
                    [Bt,Lt] = bwboundaries(obj.masks.tumorIm);
                    if ~isempty(Bt)
                        for i = 1:length(Bt)
                            bt(i) = size(Bt{i},1);
                        end
                        [bd, bdx] = sort(bt,'descend');
                    end
                    
                    % process the interfaces by first limiting the number
                    % of boundaries to contoursMaxNm.
                    % For those contoursMaxNm interfaces, delete all those
                    % interfaces that bound an area that is less than
                    % contoursThrhld percent of the largest area bound. For
                    % example, if contoursMaxNm = 7, then only select 7
                    % interfaces and delete the remaining interfaces. Of
                    % the remaining interfaces, if contoursThrhld = 10, then
                    % delete all those interfaces that bound an area
                    % less than 10 % of the largest area bound.                   
                    bnd = obj.processBoundary(bdx, Bt, countourParm);
                    % Save all of interfacess in the interfaces structure.
                    obj.interfaces.tumorIm = bnd;
                    
                    % Identify all of the interfaces that were rejected.
                    % For example, if if contoursMaxNm = 7, then identify
                    % eighth, ninth, etc. interfaces. Moreover, identify 
                    % the remaining interfaces that are less than
                    % contoursThrhld (10 percent, for example) less than
                    % the maximum area. Once identified, create a mask of
                    % these interfaces and name them 'artifact'. Remove the
                    % 'artifact' mask from the tumor mask tumorMasks_prev
                    % and save as tumor mask.
                    bndArt = obj.identifyArtifacts(bdx, Bt, countourParm);                    
                    if ~isempty(bndArt)
                        artifact = zeros(obj.n, obj.m);
                        for i = 1:length(bndArt)
                            artifact_i = roipoly(obj.ROIMask, bndArt{i}.pts(:,2), bndArt{i}.pts(:,1));
                            artifact = artifact + artifact_i;
                        end
                        obj.masks.tumorIm = tumorMasks_prev - artifact;
                        
                        % Errors occur when creating masks from interfaces,
                        % so remove these small artifaces with
                        % 'bwareaopen'.
                        obj.masks.tumorIm = bwareaopen(obj.masks.tumorIm, 30);
                        % Add the artifacts to the fibroglandular region.
                        obj.masks.fibroglandularIm = obj.masks.fibroglandularIm + (tumorMasks_prev - obj.masks.tumorIm);
                        obj.masks.glandularIm = obj.masks.fibroglandularIm + obj.masks.transitionIm;                        
                    end
                    
                otherwise
                    % Store the tumor mask, before processing
                    tumorMasks_prev = obj.masks.tumorMag;
                    
                    % Get all of the edge points for each region of
                    % malignant tissue, and sort the interfaces in terms of
                    % largest area bound by the interface to smallest area
                    % bound by the interface.
                    [Bt,Lt] = bwboundaries(obj.masks.tumorMag);
                    if ~isempty(Bt)
                        for i = 1:length(Bt)
                            bt(i) = size(Bt{i},1);
                        end
                        [bd, bdx] = sort(bt,'descend');
                    end
                    
                    % process the interfaces by first limiting the number
                    % of boundaries to contoursMaxNm.
                    % For those contoursMaxNm interfaces, delete all those
                    % interfaces that bound an area that is less than
                    % contoursThrhld percent of the largest area bound. For
                    % example, if contoursMaxNm = 7, then only select 7
                    % interfaces and delete the remaining interfaces. Of
                    % the remaining interfaces, if contoursThrhld = 10, then
                    % delete all those interfaces that bound an area
                    % less than 10 % of the largest area bound.                       
                    bnd = obj.processBoundary(bdx, Bt, countourParm);
                    % Save all of interfacess in the interfaces structure.
                    obj.interfaces.tumorMag = bnd;
                    
                    % Identify all of the interfaces that were rejected.
                    % For example, if if contoursMaxNm = 7, then identify
                    % eighth, ninth, etc. interfaces. Moreover, identify 
                    % the remaining interfaces that are less than
                    % contoursThrhld (10 percent, for example) less than
                    % the maximum area. Once identified, create a mask of
                    % these interfaces and name them 'artifact'. Remove the
                    % 'artifact' mask from the tumor mask tumorMasks_prev
                    % and save as tumor mask.
                    bndArt = obj.identifyArtifacts(bdx, Bt, countourParm);
                    
                    if ~isempty(bndArt)
                        artifact = zeros(obj.n, obj.m);
                        for i = 1:length(bndArt)
                            artifact_i = roipoly(obj.ROIMask, bndArt{i}.pts(:,2), bndArt{i}.pts(:,1));
                            artifact = artifact + artifact_i;
                        end
                        obj.masks.tumorMag = tumorMasks_prev - artifact; 
                        
                        % Errors occur when creating masks from interfaces,
                        % so remove these small artifaces with
                        % 'bwareaopen'.
                        obj.masks.tumorMag = bwareaopen(obj.masks.tumorMag, 30);
                        % Add the artifacts to the fibroglandular region.
                        obj.masks.fibroglandularMag = obj.masks.fibroglandularMag + (tumorMasks_prev - obj.masks.tumorMag);
                        obj.masks.glandularMag = obj.masks.fibroglandularMag + obj.masks.transitionMag;
                    end
            end
        end % function rejectTumorRegionArtifacts
        %------------------------------------------------------------------
        function obj = rejectGlandularRegionArtifacts(obj, component, gui)            
            % This function used to sense and then reject artifacts
            % that are incorrectly attributed to dense tissue.            
            % Inputs:
            %   component - which component of the inverse model is being
            %           segmented (i.e., real, imaginary, or magnitude.
            %   gui.contoursThrhld - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   gui.contoursMaxNm - maximum number of contours to extract from the mask of the malignant tissue            
            % Outputs:
            %   obj.masks.glandularRe, obj.masks.glandularIm, obj.masks.glandularMag - n by m mask of segmented dense tissue segmented from real, imaginary, or magnitude component.
            %   obj.interfaces.glandularRe, obj.interfaces.glandularIm, obj.interfaces.glandularMag - interfaces extracted from the glandular tissue masks of each of the components (i.e., real, imaginary, magnitude)
            %--------------------------------------------------------------  
            countourParm.contoursThrhld = gui.contoursThrhld;
            countourParm.contoursMaxNm = gui.contoursMaxNm;  
            switch component                
                case 'realEps'
                    % Store the glandular mask, before processing
                    glandularMasks_prev = obj.masks.glandularRe;
                    % Get all of the edge points for each region of
                    % dense tissue, and sort the interfaces in terms of
                    % largest area bound by the interface to smallest area
                    % bound by the interface.                   
                    [Bt,Lt] = bwboundaries(obj.masks.glandularRe, 'noholes');
                    if ~isempty(Bt)
                        for i = 1:length(Bt)
                            bt(i) = size(Bt{i},1);
                        end
                        [bd, bdx] = sort(bt,'descend');
                    end
                    
                    % process the interfaces by first limiting the number
                    % of boundaries to contoursMaxNm.
                    % For those contoursMaxNm interfaces, delete all those
                    % interfaces that bound an area that is less than
                    % contoursThrhld percent of the largest area bound. For
                    % example, if contoursMaxNm = 7, then only select 7
                    % interfaces and delete the remaining interfaces. Of
                    % the remaining interfaces, if contoursThrhld = 10, then
                    % delete all those interfaces that bound an area
                    % less than 10 % of the largest area bound.
                    % Save all of interfacess in the interfaces structure.
                    obj.interfaces.glandularRe = obj.processBoundary(bdx, Bt, countourParm);                                   
                    
                    % Identify all of the interfaces that were rejected.
                    % For example, if if contoursMaxNm = 7, then identify
                    % eighth, ninth, etc. interfaces. Moreover, identify 
                    % the remaining interfaces that are less than
                    % contoursThrhld (10 percent, for example) less than
                    % the maximum area. Once identified, create a mask of
                    % these interfaces and name them 'artifact'. Remove the
                    % 'artifact' mask from the glandular region mask glandularMasks_prev
                    % and save as glandular mask.
                    bndArt = obj.identifyArtifacts(bdx, Bt, countourParm);
                    
                    if ~isempty(bndArt)
                        artifact = zeros(obj.n, obj.m);
                        for i = 1:length(bndArt)
                            artifact_i = roipoly(obj.ROIMask, bndArt{i}.pts(:,2), bndArt{i}.pts(:,1));
                            artifact = artifact + artifact_i;
                        end
                        obj.masks.glandularRe = glandularMasks_prev - artifact;  
                        % Errors occur when creating masks from interfaces,
                        % so remove these small artifaces.
                        obj.masks.glandularRe(find(obj.masks.glandularRe == -1)) = 0;                        
                        obj.masks.glandularRe = bwareaopen(obj.masks.glandularRe, 20);
                        obj.masks.fatRe = obj.masks.fatRe + (glandularMasks_prev - obj.masks.glandularRe);
                        % Update fat, transition and fibroglandular masks accordingly
                        m = obj.updateGlandularMasks(glandularMasks_prev, obj.masks.glandularRe, obj.masks.transitionRe, obj.masks.fibroglandularRe, obj.masks.tumorRe);                        
                        obj.masks.fibroglandularRe = obj.masks.fibroglandularRe - m.deltaF;
                        obj.masks.transitionRe = obj.masks.transitionRe - m.deltaT;                    
                    end                    
                    
                case 'imEps'
                    % Store the glandular mask, before processing
                    glandularMasks_prev = obj.masks.glandularIm;
                    % Get all of the edge points for each region of
                    % dense tissue, and sort the interfaces in terms of
                    % largest area bound by the interface to smallest area
                    % bound by the interface.                                  
                    [Bt,Lt] = bwboundaries(obj.masks.glandularIm, 'noholes');
                    if ~isempty(Bt)
                        for i = 1:length(Bt)
                            bt(i) = size(Bt{i},1);
                        end
                        [bd, bdx] = sort(bt,'descend');
                    end
                    % process the interfaces by first limiting the number
                    % of boundaries to contoursMaxNm.
                    % For those contoursMaxNm interfaces, delete all those
                    % interfaces that bound an area that is less than
                    % contoursThrhld percent of the largest area bound. For
                    % example, if contoursMaxNm = 7, then only select 7
                    % interfaces and delete the remaining interfaces. Of
                    % the remaining interfaces, if contoursThrhld = 10, then
                    % delete all those interfaces that bound an area
                    % less than 10 % of the largest area bound.
                    % Save all of interfacess in the interfaces structure.
                    obj.interfaces.glandularIm = obj.processBoundary(bdx, Bt, countourParm);                                        
                    
                    % Identify all of the interfaces that were rejected.
                    % For example, if if contoursMaxNm = 7, then identify
                    % eighth, ninth, etc. interfaces. Moreover, identify 
                    % the remaining interfaces that are less than
                    % contoursThrhld (10 percent, for example) less than
                    % the maximum area. Once identified, create a mask of
                    % these interfaces and name them 'artifact'. Remove the
                    % 'artifact' mask from the glandular region mask glandularMasks_prev
                    % and save as glandular mask.
                    bndArt = obj.identifyArtifacts(bdx, Bt, countourParm);
                    
                    if ~isempty(bndArt)
                        artifact = zeros(obj.n, obj.m);
                        for i = 1:length(bndArt)
                            artifact_i = roipoly(obj.ROIMask, bndArt{i}.pts(:,2), bndArt{i}.pts(:,1));
                            artifact = artifact + artifact_i;
                        end
                        obj.masks.glandularIm = glandularMasks_prev - artifact;
                        % Errors occur when creating masks from interfaces,
                        % so remove these small artifaces.
                        obj.masks.glandularIm(find(obj.masks.glandularIm == -1)) = 0;
                        obj.masks.glandularIm = bwareaopen(obj.masks.glandularIm, 20);
                        obj.masks.fatIm = obj.masks.fatIm + (glandularMasks_prev - obj.masks.glandularIm);
                        % Update fat, transition and fibroglandular masks accordingly
                        m = obj.updateGlandularMasks(glandularMasks_prev, obj.masks.glandularIm, obj.masks.transitionIm, obj.masks.fibroglandularIm, obj.masks.tumorIm);                        
                        obj.masks.fibroglandularIm = obj.masks.fibroglandularIm - m.deltaF;
                        obj.masks.transitionIm = obj.masks.transitionIm - m.deltaT;                    
                    end   
                    
                otherwise
                    % Store the glandular mask, before processing
                    glandularMasks_prev = obj.masks.glandularMag;
                    % Get all of the edge points for each region of
                    % dense tissue, and sort the interfaces in terms of
                    % largest area bound by the interface to smallest area
                    % bound by the interface.                          
                    [Bt,Lt] = bwboundaries(obj.masks.glandularMag, 'noholes');
                    if ~isempty(Bt)
                        for i = 1:length(Bt)
                            bt(i) = size(Bt{i},1);
                        end
                        [bd, bdx] = sort(bt,'descend');
                    end
                    % process the interfaces by first limiting the number
                    % of boundaries to contoursMaxNm.
                    % For those contoursMaxNm interfaces, delete all those
                    % interfaces that bound an area that is less than
                    % contoursThrhld percent of the largest area bound. For
                    % example, if contoursMaxNm = 7, then only select 7
                    % interfaces and delete the remaining interfaces. Of
                    % the remaining interfaces, if contoursThrhld = 10, then
                    % delete all those interfaces that bound an area
                    % less than 10 % of the largest area bound.
                    % Save all of interfacess in the interfaces structure.
                    obj.interfaces.glandularMag = obj.processBoundary(bdx, Bt, countourParm);                                        
                    
                    % Identify all of the interfaces that were rejected.
                    % For example, if if contoursMaxNm = 7, then identify
                    % eighth, ninth, etc. interfaces. Moreover, identify 
                    % the remaining interfaces that are less than
                    % contoursThrhld (10 percent, for example) less than
                    % the maximum area. Once identified, create a mask of
                    % these interfaces and name them 'artifact'. Remove the
                    % 'artifact' mask from the glandular region mask glandularMasks_prev
                    % and save as glandular mask.
                    bndArt = obj.identifyArtifacts(bdx, Bt, countourParm);
                    
                    if ~isempty(bndArt)
                        artifact = zeros(obj.n, obj.m);
                        for i = 1:length(bndArt)
                            artifact_i = roipoly(obj.ROIMask, bndArt{i}.pts(:,2), bndArt{i}.pts(:,1));
                            artifact = artifact + artifact_i;
                        end
                        obj.masks.glandularMag = glandularMasks_prev - artifact;
                        
                        % Errors occur when creating masks from interfaces,
                        % so remove these small artifaces.
                        obj.masks.glandularMag(find(obj.masks.glandularMag == -1)) = 0;
                        obj.masks.glandularMag = bwareaopen(obj.masks.glandularMag, 20);
                        obj.masks.fatMag = obj.masks.fatMag + (glandularMasks_prev - obj.masks.glandularMag);
                        % Update fat, transition and fibroglandular masks accordingly
                        m = obj.updateGlandularMasks(glandularMasks_prev, obj.masks.glandularMag, obj.masks.transitionMag, obj.masks.fibroglandularMag, obj.masks.tumorMag);                        
                        obj.masks.fibroglandularMag = obj.masks.fibroglandularMag - m.deltaF;
                        obj.masks.transitionMag = obj.masks.transitionMag - m.deltaT;                    
                    end   
            end
        end % function 
        %------------------------------------------------------------------        
        function obj = mapAutoSegmentMasksToTissueTypes(obj,component)
            % This function maps the masks formed by segmenting the
            % interior into cluster to tissue types. That is the background
            % mask is mapped to background, the fatty tissue mask is mapped
            % to adipose tissue, and so on.            
            % Inputs:
            %   component - which component of the inverse model is being
            %           segmented (i.e., real, imaginary, or magnitude).
            % Outputs:
            %   obj.tissueTypes.background, obj.tissueTypes.fatRe, obj.tissueTypes.transitionRe, obj.tissueTypes.fibroglandularRe, etc.  - n by m tissue type 
            %--------------------------------------------------------------  
            % Map tissue mask to tissue type for each component
            switch component
                case 'realEps'
                    obj.tissueTypes.background = obj.masks.background;                  % background
                    obj.tissueTypes.fatRe = obj.masks.fatRe*2;                            % adipose tissue - 2
                    obj.tissueTypes.transitionRe = obj.masks.transitionRe*3;              % transition tissue - 3
                    obj.tissueTypes.fibroglandularRe = obj.masks.fibroglandularRe*4;      % fibroglandular tissue - 4.
                    obj.tissueTypes.tumorRe = obj.masks.tumorRe*5;                      % malignant tissue - 5.  
                    
                case 'imEps'
                    obj.tissueTypes.background = obj.masks.background;                  % background
                    obj.tissueTypes.fatIm = obj.masks.fatIm*2;                            % adipose tissue - 2
                    obj.tissueTypes.transitionIm = obj.masks.transitionIm*3;              % transition tissue - 3
                    obj.tissueTypes.fibroglandularIm = obj.masks.fibroglandularIm*4;      % fibroglandular tissue - 4.
                    obj.tissueTypes.tumorIm = obj.masks.tumorIm*5;                      % malignant tissue - 5.  
                    
                otherwise
                    obj.tissueTypes.background = obj.masks.background;                  % background
                    obj.tissueTypes.fatMag = obj.masks.fatMag*2;                            % adipose tissue - 2
                    obj.tissueTypes.transitionMag = obj.masks.transitionMag*3;              % transition tissue - 3
                    obj.tissueTypes.fibroglandularMag = obj.masks.fibroglandularMag*4;      % fibroglandular tissue - 4.
                    obj.tissueTypes.tumorMag = obj.masks.tumorMag*5;                      % malignant tissue - 5.  
            end
        end % function mapAutoSegmentMasksToTissueTypes
        %------------------------------------------------------------------ 
        function interfaces = getTissueInterfaces(obj, masks, countourParm)
            % This function extracts tissue interfaces from segmented
            % masks. The contour parameters set by the user determine the
            % maximum number interfaces extracted and which of the
            % interfaces extracted from the smaller regions are discared
            % (or filtered out).            
            % Inputs:
            %   masks - n by m tissue segmented from an image. Edge points
            %   from thed countour of the masks are extracted with
            %   algorihtm.
            %   countourParm.contoursThrhld - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   countourParm.contoursMaxNm - maximum number of contours to extract from the mask of the malignant tissue            
            % Outputs:
            %   interfaces.pts - edge points of the mask
            %   interfaces.area - area of region bound by interface
            %   interfaces.centre - centre of region bound by interface
            %   interfaces.perimeter - perimeter of the regoin bound by interface            
            %--------------------------------------------------------------  
            
            % step 1 - Sample countour of mask to extract boundary points of interface of the
            % reference and reconstructed models
            [Bd,Ld] = bwboundaries(masks, 'noholes');            
            
            % step 2 - Sort the boundaries in descending order according to number of points extracted
            % from the interface.
            if ~isempty(Bd)
                for i = 1:length(Bd)
                    bd(i) = size(Bd{i},1);
                end
                [bd, bdx] = sort(bd,'descend');
            end                        
            
            % step 3 - process contours by finding the N (user specified) largest contours,
            % and evaluate the centre and area of the region that
            % is enclosed by each of these contours. Finally, apply a user
            % defined threshold to the N largest contours. Reject all
            % contours for which the proportion of area of the region enclosed by the contour 
            % to area of largest region is less than the threhold value.
            interfaces = obj.processBoundary(bdx, Bd, countourParm);            
        end % function getGlandularInterface
        %------------------------------------------------------------------
        function bnd = processBoundary(obj, bdx, B, countourParm)
            % For each interface, this function converts the points to engineering units (m),
            % calls getContourGeom to calculates the area of the region
            % bound by the interface, determines the centre of the region
            % bound by the interface, and calculates the perimieter of the
            % interface.
            % The function first selects the first
            % countourParm.contoursMaxNm interfaces from the list of
            % interfaces. The remaining interfaces are discarded as
            % artifacts. For the selected interfaces, the ratio of area bound by
            % the interface and the largest area is evaluated. If the ratio is less
            % than countourParm.contoursThrhld, then the interface is
            % discaded as noise.            
            % Inputs:
            %   masks - n by m tissue segmented from an image. Edge points
            %   from thed countour of the masks are extracted with
            %   algorihtm.
            %   countourParm.contoursThrhld - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   countourParm.contoursMaxNm - maximum number of contours to extract from the mask of the malignant tissue            
            % Outputs:
            %   bnd{i}.pts - edge points of the ith mask
            %   bnd{i}.area - area of region bound by ith interface
            %   bnd{i}.centre - centre of region bound by ith interface
            %   bnd{i}.perimeter - perimeter of the regoin bound by ith interface            
            %--------------------------------------------------------------
            % Select the first countourParm.contoursMaxNm with the largest
            % number of sampled edge points.
            n = min(countourParm.contoursMaxNm, length(bdx));
            for i = 1:n
                % Get the edge points
                bnd{i}.pts = B{bdx(i)};
                % Convert the edge points to engineering units (m)
                X = [obj.xNodes(bnd{i}.pts(:,2))' obj.yNodes(bnd{i}.pts(:,1))'];
                % Call getContourGeo to calculate area, centre, and
                % perimeter
                geom = obj.getContourGeom(X(:,1), X(:,2) );
                bnd{i}.area = geom(1);
                bnd{i}.centre=[geom(2) geom(3)];
                bnd{i}.perimeter = geom(4);
                
                % Calculate the ratio of the area bound by interface with
                % largest area. Discard interface if ratio is less than countourParm.contoursThrhld
                ratio(i) = (bnd{i}.area/bnd{1}.area)*100;
                if ratio(i) < countourParm.contoursThrhld
                    ratio(i) = 0;
                end
            end
            % Return boundary points, area, centre, permiter to interface
            % structure for the query tissue tpe.
            idx = find(ratio~=0);            
            bnd = bnd([idx]);                        
        end % processBoundary        
        %------------------------------------------------------------------
        function bnd = identifyArtifacts(obj, bdx, B, countourParm) 
            % This function identifies all interfaces that are rejected as noise based
            % on the contour parameters (threshold and maximum number of interfaces).             
            % Inputs:
            %   bdx - index of list of interface extracted.
            %   B - edge pooints of interface samples.
            %   countourParm.contoursThrhld - all contours that enclose an area below the extraction threshold (i.e., ratio of area enclosed by contour to the largest area is less than threshold), are removed from the analysis
            %   countourParm.contoursMaxNm - maximum number of contours to extract from the mask of the malignant tissue            
            % Outputs:
            %   bnd{i}.pts - edge points of the ith mask
            %   bnd{i}.area - area of region bound by ith interface
            %   bnd{i}.centre - centre of region bound by ith interface
            %   bnd{i}.perimeter - perimeter of the regoin bound by ith interface            
            %--------------------------------------------------------------
            for i = 1:length(bdx)
                bnd{i}.pts = B{bdx(i)};                
                X = [obj.xNodes(bnd{i}.pts(:,2))' obj.yNodes(bnd{i}.pts(:,1))'];                
                geom = obj.getContourGeom(X(:,1), X(:,2) );
                bnd{i}.area = geom(1);
                bnd{i}.centre=[geom(2) geom(3)];
                bnd{i}.perimeter = geom(4);
                ratio(i) = (bnd{i}.area/bnd{1}.area)*100;
                if ratio(i) < countourParm.contoursThrhld
                    ratio(i) = 0;
                end
            end
            idx = find(ratio == 0);
            bnd = bnd([idx]);                        
        end % processBoundary
        %-----------------------------------------------------------------
        function masks = updateMasks(obj, ROI, bnd)
            % This function creates mask from regions bound by all 
            % interfaces that are rejected as noise.
            % Inputs:
            %   ROI - region of interest to construct masks from rejected
            %   interfaces.
            %   bnd - list of interfaces rejected as noise.
            % Outputs:
            %   mask - n by m mask created from rejected interfaces.
            %--------------------------------------------------------------
            masks = zeros(obj.n, obj.m);            
            % create updated masks from the filtered contours
            for i = 1:length(bnd)
                mask_i = roipoly(ROI, bnd{i}.pts(:,2), bnd{i}.pts(:,1));
                masks = masks + mask_i;
            end            
        end % function updateMasks
        %------------------------------------------------------------------
        function fatMasks = updateFatMasks(obj, glandularMasks_prev, glandularMasks_new, fatMasks )
            % This function updates the fat masks with artifacts removed
            % from the glaandular masks.
            % Inputs:
            %   glandularMasks_prev - n by m glandular mask with artifacts            
            %   glandularMasks_new - n by m glandular mask with artifacts removed
            %   fatMasks - n by m fat mask before artifacts from glandular
            %   masks are added.
            % Outputs:
            %   fatMasks - n by m mask with additional artifacts removed from
            %   glandular mask
            %--------------------------------------------------------------            
            fatMasks = fatMasks + (glandularMasks_prev - glandularMasks_new);
        end % updateFatMasks
        %------------------------------------------------------------------
        function m = updateGlandularMasks(obj, glandularMasks_prev, glandularMasks_new, transistionMask, fibroglandularMask, tumorMask )            
            % Errors are created when creating masks of glandular region from interfaces.
            % For example, there may be malignant within a region bound by an interface
            % that is not taking into accounnt when creating a mask of the glandular
            % structure from an interface. This function is applied to the glandualr
            % masks created from countours to correct for these types
            % errors.            
            % Inputs:
            %   glandularMasks_prev - n by m glandular mask with artifacts            
            %   glandularMasks_new - n by m glandular mask with artifacts removed
            %   transistionMask - n by m mask of region with transition
            %   or medium density tissue. This mask is constructed from the
            %   segmentation process, not with interfaces.
            %   fibroglandularMask - n by m mask of region with dense
            %   tissue. Like the transition mask, it is contructed from the
            %   segmentation process, not with interfaces.            
            %   tumorMask - n by m mask of malignant tissue.
            % Outputs:
            %   fatMasks - n by m mask with additional artifacts removed from
            %   glandular mask
            %--------------------------------------------------------------             
            m.deltaG = zeros(obj.n, obj.m); % Change in glandular mask due to errors created when creating mask with interface.
            m.deltaT = zeros(obj.n, obj.m); % Change in transition mask due to error
            m.deltaF = zeros(obj.n, obj.m); % CHange in fibroglandular mask due to error
            % First create tissue map from transition, fibroglandular, and
            % tumor masks created by segmentation process.
            tissueMap = 3*transistionMask + 4*fibroglandularMask + 5*tumorMask;
            % Itentify which regions of the artifacts are actually part of
            % the transition, fibroglandular, or malignant tissue.
            artifacts = tissueMap.*(glandularMasks_prev - glandularMasks_new);
            % Associate each region of the artifact identified to either
            % the glandular, transition, or fibroglandulr region
            m.deltaG(find(tissueMap.*glandularMasks_new == 5)) = 1;
            m.deltaT(find(artifacts == 3)) = 1;                    
            m.deltaF(find(artifacts == 4)) = 1;
        end % function updateGlandularMasks
        %------------------------------------------------------------------
        function fibroglandularMasks = updateFibroglandularMasks(obj, glandularMasks_prev, glandularMasks, fibroglandularMasks)
            % This function is updates the fibroglandular mask after removing 
            % artifacts from the glandular mask.
            % Inputs:
            %   glandularMasks_prev - n by m glandular mask with artifacts            
            %   glandularMasks_new - n by m glandular mask with artifacts removed            
            %   fibroglandularMask - n by m mask of region with dense
            %   tissue contructed from the segmentation process, not with interfaces.                       
            % Outputs:
            %   fibroglandularMasks - n by m mask with additional artifacts removed 
            %--------------------------------------------------------------
            Y = zeros(obj.n, obj.m);
            X = fibroglandularMasks + (glandularMasks_prev - glandularMasks);
            i = find( X == 2); 
            Y(i) = 1;
            fibroglandularMasks = fibroglandularMasks - Y;
        end % function updateTransitionMasks
        %------------------------------------------------------------------        
        function geom = getContourGeom(obj, x, y )
            % This function evaluates the area, centre and perimeter of a region
            % inputs: 
            %   x,y: x and y-coordinates of the contour
            % outputs:
            %   geom: geometery of region bound by closed contour
            %   geom(1) - area
            %   geom(2) - x_cen
            %   geom(3) - y_cen
            %   geom(4) - perimeter
            %--------------------------------------------------------------
            % check if inputs are same size
                if ~isequal( size(x), size(y) )
                    error( 'X and Y must be the same size');
                end
            % number of vertices
                [ x, ns ] = shiftdim( x );
                [ y, ns ] = shiftdim( y );
                [ n, c ] = size( x );
            % temporarily shift data to mean of vertices for improved accuracy
                xm = mean(x);
                ym = mean(y);
                x = x - xm*ones(n,1);
                y = y - ym*ones(n,1);
            % delta x and delta y
                dx = x( [ 2:n 1 ] ) - x;
                dy = y( [ 2:n 1 ] ) - y;
            % summations for CW boundary integrals
            % Area within closed contour    
                A = sum( y.*dx - x.*dy )/2;
                Axc = sum( 6*x.*y.*dx -3*x.*x.*dy +3*y.*dx.*dx +dx.*dx.*dy )/12;
                Ayc = sum( 3*y.*y.*dx -6*x.*y.*dy -3*x.*dy.*dy -dx.*dy.*dy )/12;
             % Length of contour (or perimeter)   
                P = sum( sqrt( dx.*dx +dy.*dy ) );
            % check for CCW versus CW contour
                if A < 0
                    A = -A;
                    Axc = -Axc;
                    Ayc = -Ayc;                
                end
            % centroidal moments
                xc = Axc / A;
                yc = Ayc / A;
            % replace mean of vertices
                x_cen = xc + xm;
                y_cen = yc + ym;         
            % geometery of region bound by closed contour
                geom = [ A  x_cen  y_cen  P ];
        end % function getContourGeom  
        %------------------------------------------------------------------
        function obj = getTissueDistributionOverEachCluster(obj)
            % This function gets the dielectric distribution of the tissues
            % over each cluster.
            % Inputs:
            %   none                                
            % Outputs:
            %   obj.tissueDistribution - distribution of dielectric
            %   properties over each cluster.
            %--------------------------------------------------------------
            [n,m] = size(obj.component);
            obj.tissueDistribution = zeros(n*m,obj.nmCl);
            for k = 1:obj.nmCl
                mask_k = zeros(n,m);
                tissue = [];
                i = find(obj.ROIClst == k);                
                mask_k(i) = 1;
                tissue = obj.component(find(mask_k.*obj.component));
                obj.tissueDistribution(1:length(tissue),k) = tissue;
            end
        end % getClsTissueDistribution                    
        %------------------------------------------------------------------
        function obj = mapClustersToTissueTypes(obj, component)           
            % This function maps the segmented clusters to tissue types.
            % Inputs:
            %   none                                
            % Outputs:
            %   obj.tissueDistribution - distribution of dielectric
            %   properties over each cluster.
            %--------------------------------------------------------------
            switch component
                case 'realEps'
                    obj.clusters.Re = obj.ROIClst;
                case 'imEps'
                    obj.clusters.Im = obj.ROIClst;
                otherwise
                    obj.clusters.Mag = obj.ROIClst;
            end
            obj.classifiedClusters =  obj.ROIClst;            
            obj.classifiedClusters(find(obj.ROIClst == 1)) = 1;                     % Cluster 1 - background
            obj.classifiedClusters(find(obj.ROIClst == 2)) = 2;                     % Cluster 2 - adipose tissue
            obj.classifiedClusters(find(obj.ROIClst > 2 & obj.ROIClst < 5)) = 3;    % Clusters 3 to 4 - transition tissue
            obj.classifiedClusters(find(obj.ROIClst > 4 & obj.ROIClst < max(obj.ROIclusters{obj.nmCl}.labels))) = 4; % Clusters 5 - max(k)-1 - fibroglandular tissue
            obj.classifiedClusters(find(obj.ROIClst == max(obj.ROIclusters{obj.nmCl}.labels))) = 5; % Cluster max(k) - malignant tissue
        end % mapClustersToTissueTypes        
        %------------------------------------------------------------------
        function obj = mapClustersToMasks(obj, component)
            % This function maps the segmented clusters to tissue masks.
            % Inputs:
            %   component - real, imaginary, or magnitude of inverse model being segmented.                                
            % Outputs:
            %   obj.clusters.Re, Im, Mag - segmented clusters of real, imaginary or magnitude of inverse model.
            %   obj.masks.background, fatRe, fatIm, fatMag, transitionRe, transitionIm, transitionMag, etc. - n by m mask of region of tissue segmented with clusters.
            %--------------------------------------------------------------
            switch component
                case 'realEps'
                    obj.clusters.Re = obj.ROIClst;
                    obj.masks.background = zeros(obj.n, obj.m);
                    obj.masks.background(find(obj.ROIClst == 1)) = 1; 
                    obj.masks.fatRe = zeros(obj.n, obj.m);
                    obj.masks.fatRe(find(obj.ROIClst == 2)) = 1;
                    obj.masks.transitionRe = zeros(obj.n, obj.m);
                    obj.masks.transitionRe(find(obj.ROIClst > 2 & obj.ROIClst < 5)) = 1;    
                    obj.masks.fibroglandularRe = zeros(obj.n, obj.m);
                    obj.masks.fibroglandularRe(find(obj.ROIClst > 4 & obj.ROIClst < max(obj.ROIclusters{obj.nmCl}.labels))) = 1; 
                    obj.masks.tumorRe = zeros(obj.n, obj.m);
                    obj.masks.tumorRe(find(obj.ROIClst == max(obj.ROIclusters{obj.nmCl}.labels))) = 1;                     
                    obj.masks.glandularRe = obj.masks.transitionRe + obj.masks.fibroglandularRe;                    
                    
                case 'imEps'
                    obj.clusters.Im = obj.ROIClst;
                    obj.masks.background = zeros(obj.n, obj.m);
                    obj.masks.background(find(obj.ROIClst == 1)) = 1;
                    obj.masks.fatIm = zeros(obj.n, obj.m);
                    obj.masks.fatIm(find(obj.ROIClst == 2)) = 1;
                    obj.masks.transitionIm = zeros(obj.n, obj.m);
                    obj.masks.transitionIm(find(obj.ROIClst > 2 & obj.ROIClst < 5)) = 1;    
                    obj.masks.fibroglandularIm = zeros(obj.n, obj.m);
                    obj.masks.fibroglandularIm(find(obj.ROIClst > 4 & obj.ROIClst < max(obj.ROIclusters{obj.nmCl}.labels))) = 1; 
                    obj.masks.tumorIm = zeros(obj.n, obj.m);
                    obj.masks.tumorIm(find(obj.ROIClst == max(obj.ROIclusters{obj.nmCl}.labels))) = 1; 
                    obj.masks.glandularIm = obj.masks.transitionIm + obj.masks.fibroglandularIm;                    
                otherwise
                    obj.clusters.Mag = obj.ROIClst;
                    obj.masks.background = zeros(obj.n, obj.m);
                    obj.masks.background(find(obj.ROIClst == 1)) = 1; 
                    obj.masks.fatMag = zeros(obj.n, obj.m);
                    obj.masks.fatMag(find(obj.ROIClst == 2)) = 1;                     
                    obj.masks.transitionMag = zeros(obj.n, obj.m);
                    obj.masks.transitionMag(find(obj.ROIClst > 2 & obj.ROIClst < 5)) = 1; 
                    obj.masks.fibroglandularMag = zeros(obj.n, obj.m);
                    obj.masks.fibroglandularMag(find(obj.ROIClst > 4 & obj.ROIClst < max(obj.ROIclusters{obj.nmCl}.labels))) = 1; 
                    obj.masks.tumorMag = zeros(obj.n, obj.m);
                    obj.masks.tumorMag(find(obj.ROIClst == max(obj.ROIclusters{obj.nmCl}.labels))) = 1; 
                    obj.masks.glandularMag = obj.masks.transitionMag + obj.masks.fibroglandularMag;                    
            end
            
        end % classifyGlndClusters
        %------------------------------------------------------------------
        function imageFwdTissueTypes( obj, xNodes, yNodes, dataPath, figFormat )
            % This is a helper function to display the tissue types within the ROI extracted
            % from the forward model.
            % Inputs:            
            %   xNodes - array - engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   yNodes - array - engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels
            %   dataPath - string - Points to results folder to store figures  
            %   figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------
            classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fat + obj.tissueTypes.transition + obj.tissueTypes.fibroglandular + obj.tissueTypes.tumor; 
            tissueType = 'forwardModel';
            obj.imageTissueTypes(classifiedTissues ,tissueType ,xNodes, yNodes, dataPath, figFormat );            
        end % function imageFwdTissueTypes        
        %------------------------------------------------------------------
        function imageInvTissueTypes( obj, xNodes, yNodes, dataPath, figFormat )
            % This is a helper function to display the tissue types within the ROI extracted
            % from the inverse model.
            % Inputs:            
            %   xNodes - array - engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   yNodes - array - engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels
            %   dataPath - string - Points to results folder to store figures  
            %   figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------
            for i=1:3
                switch i
                    case 1
                        classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatRe + obj.tissueTypes.transitionRe + obj.tissueTypes.fibroglandularRe + obj.tissueTypes.tumorRe; 
                        tissueType = 'ReEps';
                    case 2
                        classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatIm + obj.tissueTypes.transitionIm + obj.tissueTypes.fibroglandularIm + obj.tissueTypes.tumorIm; 
                        tissueType = 'ImEps';
                    otherwise
                        classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatMag + obj.tissueTypes.transitionMag + obj.tissueTypes.fibroglandularMag + obj.tissueTypes.tumorMag; 
                        tissueType = 'MagEps';
                end
                obj.imageTissueTypes(classifiedTissues ,tissueType , xNodes, yNodes,dataPath, figFormat );
            end
        end % function imageInvTissueTypes
        %------------------------------------------------------------------
        function imageTissueTypes(obj, classifiedTissues, tissueType, xNodes, yNodes, dataPath, figFormat )                                        
            % This function displays the tissue types within the ROI extracted
            % from the forward and all components of the inverse model.
            % Inputs:            
            %   classifiedTissues - n by m tissue map of fat, transition, fibroglandular, and malignant tissue types
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   xNodes - array - engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   yNodes - array - engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels
            %   dataPath - string - Points to results folder to store figures  
            %   figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------
            figure;
            cmap = jet(5);
            cmap = flipud(cmap(1:5,:));
            cmap(1,:) = [1 1 1]; % White background
            cmap(2,:) = [0 0 1]; % Blue fatty tissue mask
            cmap(3,:) = [0 1 1]; % Cyan transition tissue mask 
            cmap(4,:) = [1 0.8588 0]; % Orange fibro-glandular tissue mask
            cmap(5,:) = [1 0 0]; % Red malignant tissue mask
            colormap(cmap);            
            imagesc(imrotate(classifiedTissues,90),'XData',xNodes,'YData',-yNodes);
            ylabel('y (m)','FontSize',14);
            xlabel('x (m)','FontSize',14);
            xlim([xNodes(1) xNodes(end)])
            ylim([yNodes(1) yNodes(end)])            
            grid on;
            set(gca,'YDir','normal');
            set(gca,'fontsize',14);                
            labels = {'backgnd','fatty','transition','fibrogland.','malig.'};
            p = lcolorbar(labels,'fontsize',12);            
            set(p, 'Position', [.835 .143 .03 .78]);
            
            fileName = [tissueType '_segmentedTissueType'];            
            obj.saveFigure(fileName, dataPath, figFormat); 
        end % function imageTissueTypes
        %------------------------------------------------------------------            
        function obj = mapTissueTypesToMasks(obj, lx)
            % This function maps the segmented clusters to tissue masks.
            % Inputs:
            %   lx - (integer) - index of image of component list to apply method.            
            %--------------------------------------------------------------
            [n,m] = size(obj.classifiedClusters);
            switch obj.cmpList{lx}                                  
                case 'realEps'
                    obj.masks.fatRe = zeros(n,m);
                    obj.masks.fatRe(find(obj.classifiedClusters == 2)) = 1;
                    obj.masks.transitionRe = zeros(n,m);
                    obj.masks.transitionRe(find(obj.classifiedClusters == 3)) = 1;
                    obj.masks.fibroglandularRe = zeros(n,m);
                    obj.masks.fibroglandularRe(find(obj.classifiedClusters == 4)) = 1;
                    obj.masks.glandularRe = zeros(n,m);
                    obj.masks.glandularRe(find(obj.classifiedClusters > 2 & obj.classifiedClusters < 5)) = 1;                    
                    obj.masks.tumorRe = zeros(n,m);
                    obj.masks.tumorRe(find(obj.classifiedClusters == 5)) = 1;
                case 'imEps'                   
                    obj.masks.fatIm = zeros(n,m);
                    obj.masks.fatIm(find(obj.classifiedClusters == 2)) = 1;
                    obj.masks.transitionIm = zeros(n,m);
                    obj.masks.transitionIm(find(obj.classifiedClusters == 3)) = 1;
                    obj.masks.fibroglandularIm = zeros(n,m);
                    obj.masks.fibroglandularIm(find(obj.classifiedClusters == 4)) = 1;
                    obj.masks.glandularIm = zeros(n,m);
                    obj.masks.glandularIm(find(obj.classifiedClusters > 2 & obj.classifiedClusters < 5)) = 1;                    
                    obj.masks.tumorIm = zeros(n,m);
                    obj.masks.tumorIm(find(obj.classifiedClusters == 5)) = 1;                    
                otherwise                
                    obj.masks.fatMag = zeros(n,m);
                    obj.masks.fatMag(find(obj.classifiedClusters == 2)) = 1;
                    obj.masks.transitionMag = zeros(n,m);
                    obj.masks.transitionMag(find(obj.classifiedClusters == 3)) = 1;
                    obj.masks.fibroglandularMag = zeros(n,m);
                    obj.masks.fibroglandularMag(find(obj.classifiedClusters == 4)) = 1;
                    obj.masks.glandularMag = zeros(n,m);
                    obj.masks.glandularMag(find(obj.classifiedClusters > 2 & obj.classifiedClusters < 5)) = 1;                    
                    obj.masks.tumorMag = zeros(n,m);
                    obj.masks.tumorMag(find(obj.classifiedClusters == 5)) = 1;                 
            end % switch           
        end % function mapTissueTypesToMasks
        %------------------------------------------------------------------
        function imageTissue(obj, mask, xNodes, yNodes)
            % This is a helper function used to display the component of
            % the inverse model that is being segmented. It is used
            % primarilty for troubleshooting purposes.
            % Inputs:
            %   mask - n by m mask applied to the component of the inverse model being segmented to display to the user what is being input to the segmentation method.
            %   xNodes - array - engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   yNodes - array - engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels
            %--------------------------------------------------------------
            figure;             
            cmap = jet(256);
            %cmap = flipud(cmap(1:256,:));
            cmap(1,:) = [1 1 1]; % White background
            colormap(cmap);            
            I = mask.*obj.component;
            imagesc(imrotate(I,90),'XData',xNodes,'YData',-yNodes);
            ylabel('y (m)','FontSize',14);
            xlabel('x (m)','FontSize',14);
            xlim([xNodes(1) xNodes(end)])
            ylim([yNodes(1) yNodes(end)])            
            grid on;
            set(gca,'YDir','normal');
            set(gca,'fontsize',14);  
            colorbar
        end % function imageTissue
        %------------------------------------------------------------------
        function imageMasks( obj, xNodes, yNodes, refModel, dataPath, figFormat)
            % This function displays the masks for the real, imaginary,
            % and magnitude component of each tissue type (fat, transition,
            % fibroglandular, glandular, and malignant). 
            % The titles have been commented out to 'declutter' the
            % images.
            % Inputs:            
            %   xNodes - array - engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   yNodes - array - engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels
            %   refModel - object - masks of forward model tissue type masks are stored in structure.
            %   dataPath - string - Points to results folder to store figures  
            %   figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------
            figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.90, 0.90]);
            annotation('rectangle',[0.05,0.935,0.93,0.060],'Edgecolor','k');
            annotation('rectangle',[0.002,0.05,0.03,0.90],'Edgecolor','k');    
            
            % Display fat masks
            subplot(3,4,1);imagesc(imrotate(refModel.masks.fat,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');
            title({'Forward Model';'Real'},'Units', 'normalized','Position',[0.5,1.00,0] );
            ylabel('fat masks (m)','Units', 'normalized','Position',[-0.22,0.50,0]);            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.55 0.92 1.30 1.30]) % stretch its width and height                        
            drawnow            
                        
            subplot(3,4,2);imagesc(imrotate(obj.masks.fatRe,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');set(gca, 'yticklabel','');
            title('Real', 'Units', 'normalized','Position',[0.5,1.00,0] )
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.88 0.92 1.30 1.30]) % stretch its width and height 
            drawnow
            
            subplot(3,4,3);imagesc(imrotate(obj.masks.fatIm,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');           
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');set(gca, 'yticklabel','');
            title({'Inverse Model';'Imaginary'}, 'Units', 'normalized','Position',[0.5,1.00,0] )
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.96 0.92 1.30 1.30]) % stretch its width and height
            drawnow
            
            subplot(3,4,4);imagesc(imrotate(obj.masks.fatMag,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');set(gca, 'yticklabel','');
            title('Magnitude', 'Units', 'normalized','Position',[0.5,1.00,0] )
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.0 0.92 1.30 1.30]) % stretch its width and height
            drawnow
            
            % Display glandular masks
            subplot(3,4,5);imagesc(imrotate(refModel.masks.glandular,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');           
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');            
            ylabel('glandular masks (m)','Units', 'normalized','Position',[-0.22,0.50,0]);             
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.55 0.900 1.30 1.30]) % stretch its width and height                         
            drawnow            
                        
            subplot(3,4,6);imagesc(imrotate(obj.masks.glandularRe,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');
            %title([recType ' real fat'], 'FontSize',10);
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');set(gca, 'yticklabel','');
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.88 0.900 1.30 1.30]) % stretch its width and height 
            drawnow
            
            subplot(3,4,7);imagesc(imrotate(obj.masks.glandularIm,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');set(gca, 'yticklabel','');
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.96 0.900  1.30 1.30]) % stretch its width and height
            drawnow
            
            subplot(3,4,8);imagesc(imrotate(obj.masks.glandularMag,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');set(gca, 'yticklabel','');
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.00 0.900 1.30 1.30]) % stretch its width and height
            drawnow
            
            % Display tumor masks
            subplot(3,4,9);imagesc(imrotate(refModel.masks.tumor,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);                      
            ylabel('Tumor masks (m)','Units', 'normalized','Position',[-0.22,0.50,0]); 
            xlabel('(m)');
            xlabh = get(gca,'xLabel');
            xtickangle(90)
            set(xlabh,'Position',get(xlabh,'Position') + [0 .01 0]);
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position            
            set(gca,'position',sub_pos.*[0.55 0.72 1.30 1.30]) % stretch its width and height
            drawnow            
                        
            subplot(3,4,10);imagesc(imrotate(obj.masks.tumorRe,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            xlabel('(m)');
            xlabh = get(gca,'xLabel');
            set(xlabh,'Position',get(xlabh,'Position') + [0 .01 0]);
            set(gca, 'yticklabel','');
            xtickangle(90)
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.88 0.72 1.30 1.30]) % stretch its width and height
            drawnow
            
            subplot(3,4,11);imagesc(imrotate(obj.masks.tumorIm,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            xlabel('(m)');
            xlabh = get(gca,'xLabel');
            set(xlabh,'Position',get(xlabh,'Position') + [0 .01 0]);
            set(gca, 'yticklabel','');
            xtickangle(90);
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.96 0.72  1.30 1.30]) % stretch its width and height
            drawnow
            
            subplot(3,4,12);imagesc(imrotate(obj.masks.tumorMag,90),'XData',xNodes,'YData',-yNodes);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position
            colormap(gca,'parula')            
            set(gca,'fontsize',14);
            xlabel('(m)');
            xlabh = get(gca,'xLabel');
            set(xlabh,'Position',get(xlabh,'Position') + [0 .01 0]);
            set(gca, 'yticklabel','');
            xtickangle(90);            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.00 0.72 1.30 1.30]) % stretch its width and height
            drawnow
            
            obj.saveFigure('tissueMasks_all', dataPath, figFormat);
        end % function imageMasks      
        %------------------------------------------------------------------
        function imageThreshTissueMaps( obj, refModel, fe, gui)
            % This function displays the tissue maps when tissue types have been
            % segmented with the threholding technique.
            % Inputs:
            %   refModel - object - tissue types of forward model are stored in structure.
            %   fe.fwdMdl.complexPermittivity - n by m - forward model complex permittivity map.
            %   fe.invMdl.complexPermittivity - n by m - inverse model complex permittivity map.
            %   fe.xNodes - array - engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   fe.yNodes - array - engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels            
            %   gui.dataPath - string - Points to results folder to store figures  
            %   gui.figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------
            figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.95, 0.80]);
            annotation('rectangle',[0.05,0.920,0.93,0.080],'Edgecolor','k');
            annotation('rectangle',[0.002,0.05,0.03,0.90],'Edgecolor','k');  
            
            % Display permittivity maps
            % Forward model - 1, Re(epr) - 2, Im(epr) - 3, Mag(epr) - 4
            subplot(2,4,1);imagesc(imrotate(real(fe.fwdMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;               
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');             
            ylabel('Permittivity (m)','Units', 'normalized','Position',[-0.27,0.50,0]);
            title({'Forward Model';'Real'},'Units', 'normalized','Position',[0.5,1.00,0] );
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.55 0.93 1.10 1.10]) % stretch its width and height                        
            drawnow     
            
            subplot(2,4,2);imagesc(imrotate(real(fe.invMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;    
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');
            set(gca, 'yticklabel','');
            title('Real','Units', 'normalized','Position',[0.5,1.00,0] );            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.89 0.93 1.10 1.10]) % stretch its width and height                     
            drawnow
            
            subplot(2,4,3);imagesc(imrotate(-imag(fe.invMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;    
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            title({'Inverse model';'Imaginary'},'Units', 'normalized','Position',[0.5,1.00,0] ); 
            set(gca, 'xticklabel','');
            set(gca, 'yticklabel','');            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.965 0.93 1.10 1.10]) % stretch its width and height                     
            drawnow
            
            subplot(2,4,4);imagesc(imrotate(abs(fe.invMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes);
            colormap(gca, gui.figColorMap); 
            ch = colorbar;
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position  
            title('Magnitude','Units', 'normalized','Position',[0.5,1.00,0] );            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');
            set(gca, 'yticklabel','');            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.00 0.93 1.10 1.10]) % stretch its width and height                     
            drawnow
            
            % Display the tissue type images
            % Forward model - 5, Re(epr) - 6, Im(epr) - 7, Mag(epr) - 8
            
            subplot(2,4,5);            
            classifiedTissues = refModel.tissueTypes.background + refModel.tissueTypes.fat + refModel.tissueTypes.transition + refModel.tissueTypes.fibroglandular + refModel.tissueTypes.tumor;
            imagesc(imrotate(classifiedTissues,90),'XData',fe.xNodes,'YData',-fe.yNodes);                     
            cmap = obj.getTissueColorMap();            
            colormap(gca,cmap);            
            set(gca,'fontsize',14);            
            ylabel('Tissue Type (m)','Units', 'normalized','Position',[-0.27,0.50,0]);
            xlabel('(m)');
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.55 1.05 1.11 1.11]) % stretch its width and height                        
            drawnow
            
            subplot(2,4,6);            
            classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatRe + obj.tissueTypes.transitionRe + obj.tissueTypes.fibroglandularRe + obj.tissueTypes.tumorRe;
            imagesc(imrotate(classifiedTissues,90),'XData',fe.xNodes,'YData',-fe.yNodes);                     
            colormap(gca,cmap);            
            set(gca,'fontsize',14);
            set(gca, 'yticklabel','');
            xlabel('(m)');
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.885 1.05 1.11 1.11]) % stretch its width and height                        
            drawnow     
            
            subplot(2,4,7);            
            classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatIm + obj.tissueTypes.transitionIm + obj.tissueTypes.fibroglandularIm + obj.tissueTypes.tumorIm;
            imagesc(imrotate(classifiedTissues,90),'XData',fe.xNodes,'YData',-fe.yNodes);                     
            colormap(gca,cmap);            
            set(gca,'fontsize',14);
            set(gca, 'yticklabel','');
            xlabel('(m)');
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.965 1.05 1.11 1.11]) % stretch its width and height                        
            drawnow 
            
            subplot(2,4,8);            
            classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatMag + obj.tissueTypes.transitionMag + obj.tissueTypes.fibroglandularMag + obj.tissueTypes.tumorMag;
            imagesc(imrotate(classifiedTissues,90),'XData',fe.xNodes,'YData',-fe.yNodes);                     
            colormap(gca,cmap);
            cbh = colorbar;
            cbh.Ticks = [1, 2, 3, 4, 5];            
            cbh.TickLabels = {'backgrnd','fatty','transtion','fibrogld','malignant'};           
            set(gca,'fontsize',14);
            set(gca, 'yticklabel','');
            xlabel('(m)');
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.00 1.05 1.60 1.10]) % stretch its width and height                        
            cbh.Position = cbh.Position - [0.01 0 0 0];
            drawnow
            
            obj.saveFigure('tissueTypeImages_ThresholdMethod_all', gui.dataPath, gui.figFormat);            
        end % function imageThreshTissueMaps
        %------------------------------------------------------------------
        function imageAutoTissueMaps( obj, refModel, fe, gui)
            % This function images the tissue maps when tissue types have been
            % segmented with the iterative machine learning technique.
            % Therefore, an image of the final iteration of the clusters is
            % provided with each tissue map so that the user can compare the inverse 
            % model with the clusters used to segment the image and the final tissue map
            % fosr which the clusters are mapped to.
            % Inputs:
            %   refModel - object - tissue types of forward model are stored in structure.
            %   fe.fwdMdl.complexPermittivity - n by m - forward model complex permittivity map.
            %   fe.invMdl.complexPermittivity - n by m - inverse model complex permittivity map.
            %   fe.xNodes - array - engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   fe.yNodes - array - engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels            
            %   gui.dataPath - string - Points to results folder to store figures  
            %   gui.figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------        
            figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.95, 0.95]);
            annotation('rectangle',[0.05,0.935,0.93,0.065],'Edgecolor','k');
            annotation('rectangle',[0.002,0.05,0.03,0.90],'Edgecolor','k');
            annotation('textarrow', [0.013 0.013], [0.3 0.3], 'String', 'Cluster map (m)','HeadStyle' ,'none','LineStyle' ,'none', 'FontName','Helvetica' ,'FontSize', 14,'TextRotation',90);
            
            % Display permittivity maps
            % Forward model - 1, Re(epr) - 2, Im(epr) - 3, Mag(epr) - 4
            subplot(3,4,1);imagesc(imrotate(real(fe.fwdMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;             
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');                  
            ylabel('Permittivity map (m)','Units', 'normalized','Position',[-0.27,0.50,0]);
            title({'Forward Model';'Real'},'Units', 'normalized','Position',[0.5,1.02,0] );
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.55 0.92 1.15 1.30]) % stretch its width and height                        
            drawnow     
            
            subplot(3,4,2);imagesc(imrotate(real(fe.invMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;    
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');
            set(gca, 'yticklabel','');
            title('Real','Units', 'normalized','Position',[0.5,1.02,0] );
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.92 0.92 1.15 1.30]) % stretch its width and height                     
            drawnow
            
            subplot(3,4,3);imagesc(imrotate(-imag(fe.invMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;    
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            title({'Inverse model';'Imaginary'},'Units', 'normalized','Position',[0.5,1.02,0] );
            set(gca, 'xticklabel','');
            set(gca, 'yticklabel','');            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.00 0.92 1.15 1.30]) % stretch its width and height                     
            drawnow
            
            subplot(3,4,4);imagesc(imrotate(abs(fe.invMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes);
            colormap(gca, gui.figColorMap); 
            ch = colorbar;
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            title('Magnitude','Units', 'normalized','Position',[0.5,1.02,0] );
            set(gca, 'xticklabel','');
            set(gca, 'yticklabel','');            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.03 0.92 1.15 1.30]) % stretch its width and height                     
            drawnow
            
            % Display the tissue type images
            % Forward model - 5, Re(epr) - 6, Im(epr) - 7, Mag(epr) - 8            
            subplot(3,4,5);            
            classifiedTissues = refModel.tissueTypes.background + refModel.tissueTypes.fat + refModel.tissueTypes.transition + refModel.tissueTypes.fibroglandular + refModel.tissueTypes.tumor;
            imagesc(imrotate(classifiedTissues,90),'XData',fe.xNodes,'YData',-fe.yNodes);                     
            cmap = obj.getTissueColorMap();            
            colormap(gca,cmap);            
            set(gca,'fontsize',14);
            ylabel('Tissue type map (m)','Units', 'normalized','Position',[-0.27,0.50,0]);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.55 0.88 1.15 1.30]) % stretch its width and height                        
            drawnow
            
            subplot(3,4,6);            
            classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatRe + obj.tissueTypes.transitionRe + obj.tissueTypes.fibroglandularRe + obj.tissueTypes.tumorRe;
            imagesc(imrotate(classifiedTissues,90),'XData',fe.xNodes,'YData',-fe.yNodes);                     
            colormap(gca,cmap);            
            set(gca,'fontsize',14);
            set(gca, 'yticklabel','');
            set(gca, 'xticklabel','');                        
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.92 0.88 1.15 1.30]) % stretch its width and height                        
            drawnow     
            
            subplot(3,4,7);            
            classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatIm + obj.tissueTypes.transitionIm + obj.tissueTypes.fibroglandularIm + obj.tissueTypes.tumorIm;
            imagesc(imrotate(classifiedTissues,90),'XData',fe.xNodes,'YData',-fe.yNodes);                     
            colormap(gca,cmap);            
            set(gca,'fontsize',14);
            set(gca, 'yticklabel','');
            set(gca, 'xticklabel','');             
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.00 0.88 1.15 1.30]) % stretch its width and height                        
            drawnow 
            
            subplot(3,4,8);            
            classifiedTissues = obj.tissueTypes.background + obj.tissueTypes.fatMag + obj.tissueTypes.transitionMag + obj.tissueTypes.fibroglandularMag + obj.tissueTypes.tumorMag;
            imagesc(imrotate(classifiedTissues,90),'XData',fe.xNodes,'YData',-fe.yNodes);                     
            colormap(gca,cmap);
            cbh = colorbar;
            cbh.Ticks = [1, 2, 3, 4, 5];            
            cbh.TickLabels = {'bckgnd','fatty','transt.','fibrogld','malig.'};           
            set(gca,'fontsize',14);
            set(gca, 'yticklabel','');
            set(gca, 'xticklabel','');            
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.03 0.88 1.15 1.30]) % stretch its width and height                        
            cbh.Position = cbh.Position - [0.01 0 0 0];
            drawnow
            
            % Display the final iteration of the clusters that were mapped
            % to tissue types
            % No Fwd model clusters - 9 , Re(epr) - 10, Im(epr) - 11, Mag(epr) - 12
            subplot(3,4,10);            
            imagesc(imrotate(obj.clusters.Re,90),'XData',fe.xNodes,'YData',-fe.yNodes);            
            xlim([fe.xNodes(1) fe.xNodes(end)])
            ylim([fe.yNodes(1) fe.yNodes(end)]) 
            %axis('equal');            
            set(gca,'YDir','normal');
            set(gca,'fontsize',14);            
            nmCl = max(max(obj.clusters.Re));
            cmap = jet(nmCl);
            t = 1:1:nmCl;
            t1 = [];
            for i = 1:nmCl                
                t1{i} = num2str(t(i));
            end
            cmap(1,:) = [1 1 1];
            cmap(2,:) = [0 0 1];
            cmap(length(cmap),:) = [0 0 0];
            colormap(gca, cmap);            
            colorbar('Ticks', t,'TickLabels', t1);
            cbh = colorbar;            
            xlabel('(m)','Units', 'normalized','Position',[0.50,-0.09,0]);            
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            cbh.Position = cbh.Position - [0.00000 0 0 0];
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.92 0.51 1.65 1.30]) % stretch its width and height                        
            drawnow
            
            subplot(3,4,11);            
            imagesc(imrotate(obj.clusters.Im,90),'XData',fe.xNodes,'YData',-fe.yNodes);            
            xlim([fe.xNodes(1) fe.xNodes(end)])
            ylim([fe.yNodes(1) fe.yNodes(end)]) 
            %axis('equal');            
            set(gca,'YDir','normal');
            set(gca,'fontsize',14);
            nmCl = max(max(obj.clusters.Im));
            cmap = jet(nmCl);
            t = 1:1:nmCl;
            t1 = [];
            for i = 1:nmCl                
                t1{i} = num2str(t(i));
            end
            cmap(1,:) = [1 1 1];
            cmap(2,:) = [0 0 1];
            cmap(length(cmap),:) = [0 0 0];
            colormap(gca,cmap);            
            colorbar('Ticks', t,'TickLabels', t1);
            cbh = colorbar;
            set(gca, 'yticklabel','');
            xlabel('(m)','Units', 'normalized','Position',[0.50,-0.09,0]); 
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            cbh.Position = cbh.Position - [0.000 0 0 0];
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.00 0.51 1.65 1.30]) % stretch its width and height                        
            drawnow
            
            subplot(3,4,12);            
            imagesc(imrotate(obj.clusters.Mag,90),'XData',fe.xNodes,'YData',-fe.yNodes);            
            xlim([fe.xNodes(1) fe.xNodes(end)])
            ylim([fe.yNodes(1) fe.yNodes(end)]) 
            %axis('equal');            
            set(gca,'YDir','normal');
            set(gca,'fontsize',14);
            
            nmCl = max(max(obj.clusters.Mag));
            cmap = jet(nmCl);
            t = 1:1:nmCl;
            t1 = [];
            for i = 1:nmCl                
                t1{i} = num2str(t(i));
            end
            cmap(1,:) = [1 1 1];
            cmap(2,:) = [0 0 1];
            cmap(length(cmap),:) = [0 0 0];
            colormap(gca, cmap);            
            colorbar('Ticks', t,'TickLabels', t1);
            cbh = colorbar;            
            set(gca, 'yticklabel','');
            xlabel('(m)','Units', 'normalized','Position',[0.50,-0.09,0]); 
            set(gca,'YDir','normal');
            cbh.Position = cbh.Position - [0.000 0 0 0];
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);                      
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.03 0.51 1.65 1.30]) % stretch its width and height                        
            drawnow
            
            obj.saveFigure('tissueTypeImages_AutoMethod_all', gui.dataPath, gui.figFormat);

        end % function imageAutoTissueMaps
        %------------------------------------------------------------------
        function saveFigure(obj, fileName, dataPath, figFormat)
            % This is an helper function called by a figure display function
            % to save the figure in the format requested by the user.
            % Inputs:
            %   fileName - string - name of file to save figure.           
            %   dataPath - string - Points to results folder to store figures  
            %   figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------        
            switch figFormat
                case '-DoNotSave'
                    % Do not save files
                case '-fig'
                    savefig([dataPath 'figures/' fileName]);
                otherwise
                    export_fig([dataPath 'figures/' fileName], figFormat);
            end            
        end % function saveFigure
        %------------------------------------------------------------------
        function cmap = getTissueColorMap(obj)
            % This is an helper function used by figure display function
            % used to generate the color map used to display the tissue
            % type map.
            % Inputs:
            %   none
            % Outputs:
            %   cmap - array - color map used to display tissue type map.
            %--------------------------------------------------------------    
            cmap = jet(5);
            cmap = flipud(cmap(1:5,:));
            cmap(1,:) = [1 1 1]; % White background
            cmap(2,:) = [0 0 1]; % Blue fatty tissue mask
            cmap(3,:) = [0 1 1]; % Cyan transition tissue mask            
            cmap(4,:) = [1 0.8588 0]; % Orange fibro-glandular tissue mask
            cmap(5,:) = [1 0 0]; % Red malignant tissue mask
        end % function getTissueColorMap
        %------------------------------------------------------------------        
        function cmap = getTumorMaskColorMap(obj)
            % This is an helper function used by figure display function
            % used to generate the color map used to display the union of
            % masks from real, imaginary, and magnitude components of
            % inverse model.
            % Inputs:
            %   none
            % Outputs:
            %   cmap - array - color map used to display tissue type map.
            %--------------------------------------------------------------                          
            cmap = jet(4);
            cmap = flipud(cmap(1:3,:));
            cmap(1,:) = [1 1 1]; % No masks
            cmap(2,:) = [0.86 0.86 0.86]; % One mask
            cmap(3,:) = [0.55 0.55 0.55]; % Two masks intersect
            cmap(4,:) = [0 0 0]; % Three masks intersect
        end % function getTumorMaskColorMap
        %------------------------------------------------------------------
        function imgInterfaceOnFwdMdl(obj, fe, gui, refModel)
            % This function displays interfaces of masks from the forward and inverse model
            % superimposed onto the forward model, and the union of the tumor masks
            % for the real, imaginary, and magnitude components.
            % Inputs:
            %   refModel.interfaces.tumor - structure - forward model interface of the tumor
            %   fe.fwdMdl.complexPermittivity - n by m - forward model complex permittivity map.
            %   fe.invMdl.complexPermittivity - n by m - inverse model complex permittivity map.
            %   fe.xNodes - array - engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   fe.yNodes - array - engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels            
            %   gui.dataPath - string - Points to results folder to store figures  
            %   gui.figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------                 
            scTr = 0.3; % transparency scale (0 - transparent, 1 - no transparency)
            figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.90, 0.60]);
            
            % Display permittivity maps
            % Forward model - 1, Re(epr) - 2, Im(epr) - 3, Mag(epr) - 4
            subplot(1,3,1);
            imagesc(imrotate(ones(fe.n,fe.m) + obj.masks.tumorRe + obj.masks.tumorIm + obj.masks.tumorMag,90),'XData',fe.xNodes,'YData',-fe.yNodes);
            cmap = obj.getTumorMaskColorMap();
            colormap(gca, cmap);            
            cbh = colorbar;
            cbh.Ticks = [1, 2, 3, 4];                        
            cbh.TickLabels = {'0','1','2','3'};
            ylabel({'Number of intersecting masks';'(m)'},'Units', 'normalized','Position',[-0.05,0.50,0]);
            %set(gca, 'yticklabel',''); 
            title('Union of tumor masks','Units', 'normalized','Position',[0.5,1.02,0] );
            xlabel('(m)','FontSize',14);
            xlim([fe.xNodes(1) fe.xNodes(end)])
            ylim([fe.yNodes(1) fe.yNodes(end)])            
            grid on;
            set(gca,'YDir','normal');
            set(gca,'fontsize',14); 
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.450 1.7 1.73 0.80]) % stretch its width and height
            hold on;
            bdRef = refModel.interfaces.tumor;
            for i=1:length(bdRef)
                Xref = [fe.xNodes(bdRef{i}.pts(:,2))' fe.yNodes(bdRef{i}.pts(:,1))']; 
                if i == length(bdRef)
                    pRec = plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
                else
                    plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on                    
                end                    
            end
            drawnow
            
            subplot(1,3,2);
            imagesc(imrotate(real(fe.fwdMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes, 'AlphaData', scTr);hold on;
            colormap(gca, gui.figColorMap);                      
            lh = legend('Box','off');
            lh.NumColumns = 2;
            lh.FontSize = 10;
            line([50,60], [75, 85], 'linewidth',3,'Color','red');    % Forward model
            line([50,60], [75, 85], 'linewidth',3,'Color','green');  % Im reconstruction          
            line([50,60], [75, 85], 'linewidth',3,'Color','blue');   % Real reconstruction
            line([50,60], [75, 85], 'linewidth',3,'Color','black');  % Mag reconstruction
            set(gca,'fontsize',14);            
            set(gca,'YDir','normal');
            bdRef = refModel.interfaces.tumor;
            for i=1:length(bdRef)
                Xref = [fe.xNodes(bdRef{i}.pts(:,2))' fe.yNodes(bdRef{i}.pts(:,1))']; 
                if i == length(bdRef)
                    pRef = plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
                else
                    plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
                end
            end
            bdRec = obj.interfaces.tumorRe;
            for i=1:length(bdRec)
                Xrec = [fe.xNodes(bdRec{i}.pts(:,2))' fe.yNodes(bdRec{i}.pts(:,1))']; 
                if i == length(bdRec)
                    pRec = plot(Xrec(:,2), Xrec(:,1), '-b','LineWidth',4);hold on                    
                else
                    plot(Xrec(:,2), Xrec(:,1), '-b','LineWidth',4);hold on                    
                end                    
            end
            bdRec = obj.interfaces.tumorIm;
            for i=1:length(bdRec)
                Xrec = [fe.xNodes(bdRec{i}.pts(:,2))' fe.yNodes(bdRec{i}.pts(:,1))']; 
                if i == length(bdRec)
                    pRec = plot(Xrec(:,2), Xrec(:,1), '-g','LineWidth',4);hold on                    
                else
                    plot(Xrec(:,2), Xrec(:,1), '-g','LineWidth',4);hold on                    
                end                    
            end
            bdRec = obj.interfaces.tumorMag;
            for i=1:length(bdRec)
                Xrec = [fe.xNodes(bdRec{i}.pts(:,2))' fe.yNodes(bdRec{i}.pts(:,1))']; 
                if i == length(bdRec)
                    pRec = plot(Xrec(:,2), Xrec(:,1), '-k','LineWidth',4);hold on                    
                else
                    plot(Xrec(:,2), Xrec(:,1), '-k','LineWidth',4);hold on                    
                end                    
            end            
            xlim([fe.xNodes(1) fe.xNodes(end)])
            ylim([fe.yNodes(1) fe.yNodes(end)])
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            title('Tumor interfaces','Units', 'normalized','Position',[0.5,1.02,0] );
            xlabel('(m)');            
            set(gca, 'yticklabel',''); 
            set(gca,'position',p); % restore position
            legend({'Fwd','Im(\epsilon(r))', 'Re(\epsilon(r))', 'Mag(\epsilon(r))'});
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.97 1.7 1.35 0.80]) % stretch its width and height                        
            drawnow 
            
            subplot(1,3,3);
            imagesc(imrotate(real(fe.fwdMdl.complexPermittivity),90),'XData',fe.xNodes,'YData',-fe.yNodes, 'AlphaData', scTr);hold on;
            colormap(gca, gui.figColorMap);            
            set(gca,'fontsize',14);            
            set(gca,'YDir','normal');
            bdRef = refModel.interfaces.glandular;
            for i=1:length(bdRef)
                Xref = [fe.xNodes(bdRef{i}.pts(:,2))' fe.yNodes(bdRef{i}.pts(:,1))']; 
                if i == length(bdRef)
                    pRef = plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
                else
                    plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
                end
            end
            bdRec = obj.interfaces.glandularRe;
            for i=1:length(bdRec)
                Xrec = [fe.xNodes(bdRec{i}.pts(:,2))' fe.yNodes(bdRec{i}.pts(:,1))']; 
                if i == length(bdRec)
                    pRec = plot(Xrec(:,2), Xrec(:,1), '-b','LineWidth',4);hold on                    
                else
                    plot(Xrec(:,2), Xrec(:,1), '-b','LineWidth',4);hold on                    
                end                    
            end
            bdRec = obj.interfaces.glandularIm;
            for i=1:length(bdRec)
                Xrec = [fe.xNodes(bdRec{i}.pts(:,2))' fe.yNodes(bdRec{i}.pts(:,1))']; 
                if i == length(bdRec)
                    pRec = plot(Xrec(:,2), Xrec(:,1), '-g','LineWidth',4);hold on                    
                else
                    plot(Xrec(:,2), Xrec(:,1), '-g','LineWidth',4);hold on                    
                end                    
            end
            bdRec = obj.interfaces.glandularMag;
            for i=1:length(bdRec)
                Xrec = [fe.xNodes(bdRec{i}.pts(:,2))' fe.yNodes(bdRec{i}.pts(:,1))']; 
                if i == length(bdRec)
                    pRec = plot(Xrec(:,2), Xrec(:,1), '-k','LineWidth',4);hold on                    
                else
                    plot(Xrec(:,2), Xrec(:,1), '-k','LineWidth',4);hold on                    
                end                    
            end
            xlim([fe.xNodes(1) fe.xNodes(end)])
            ylim([fe.yNodes(1) fe.yNodes(end)])
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            set(0, 'DefaultAxesTitleFontWeight','normal');
            title('Glandular interfaces','Units', 'normalized','Position',[0.5,1.02,0] );
            xlabel('(m)');            
            set(gca, 'yticklabel',''); 
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position            
            set(gca,'position',sub_pos.*[1.01 1.7 1.35 0.80]) % stretch its width and height
            drawnow   
                        
            switch gui.figFormat
                case '-DoNotSave'
                    % Do not save files
                case '-fig'
                    savefig([gui.dataPath 'figures/' 'tissueInterfaceImages_all']);
                otherwise
                    export_fig([gui.dataPath 'figures/' 'tissueInterfaceImages_all'], gui.figFormat);
            end
            
        end % imgInterfaceOnFwdMdl
        %-----------------------------------------------------------------
        
    end % of methods
    %----------------------------------------------------------------------
end % of ClassDef
%==========================================================================

