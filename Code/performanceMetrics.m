classdef performanceMetrics
    % Class performanceMetrics
    % Methods applied to masks and interfaces to analyze geometric and
    % dielectric properties. For the interfaces, the edge points extracted
    % from the boundary of the reconstructed masks collectively represent
    % estimated interfaces of the reconstructed regions. The estimated points
    % are compared with the corresponding boundary of the reference mask to 
    % evaluate how accurately the interface is extracted by the reconstruction
    % algorithm.
    % Data structure holds values of the metrics.
%======================================================================    
    properties        
        n                       % Dimension of model/mask along y-axis (i.e., number of rows)
        m                       % Dimension of model/mask along x-axis (i.e., number of columns)
        xNodes                  % coordinates of centre of cell along x-axis
        yNodes                  % coordinates of centre of cell along y-axis
        interfaces              % Interface analysis metrics
        geometric               % Geometric property analysis metrics
        dielectric              % Dielectric property analysis metrics
    end
%======================================================================
    methods
        %------------------------------------------------------------------
        function obj = performanceMetrics( xNodes, yNodes, n, m )
            % Constructor method to create instance of object. When instance
            % object is created, the method initializes the xNodes, yNodes,
            % n and m properties.
            % Inputs:            
            %   xNodes - (array of floats) engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
            %   yNodes - (array of floats) engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels
            %   n - (integer) number of rows in complex permittivity map
            %   m - (integer) number of columns in complex permittivity map            
            %--------------------------------------------------------------
            obj.xNodes = xNodes;
            obj.yNodes = yNodes;
            obj.n = n;
            obj.m = m;
        end % Constructor performanceMetrics        
        %------------------------------------------------------------------                
        function obj = interfaceAnalysis(obj, fe, fldnmfwd, fldnmInv, refModel, recModel, gui)                
            % This function is applied to interfaces extracted from the tissue masks 
            % to carry out the interface analysis.            
            % Inputs:
            %   fe - finite element object - structure stores the forward and inverse models.
            %   fldnmfwd - string - field name of the reference model object
            %   fldnmInv - string - field name of the reconstructed model object
            %   refModel - segmentation object - holds the reference model masks, interfaces, and tissue types.
            %   recModel - segmentation object - holds the  reconstructed model masks, interfaces, and tissue types.
            %   gui - GUI object - holds the data path, figure format.            
            % Outputs:
            %   obj.interfaces.(fldnmInv).distances  - array - distance between each point on the reconstructed interface and the nearest point on the reference interface.            
            %--------------------------------------------------------------
            % Get the interfaces for the query reconstructed region and the
            % corresponding reference region
            bdRef = refModel.interfaces.(fldnmfwd);
            bdRec = recModel.interfaces.(fldnmInv);
            % step 1 - find the distance between each point of the reconstructed
            % interface with the corresponding referenc surface.  
            obj.interfaces.(fldnmInv).distances = obj.cmpBnd(bdRef, bdRec);            
            for n = 1:length(obj.interfaces.(fldnmInv).distances)
                obj.displayInterface( fe, fldnmfwd, fldnmInv, bdRef, bdRec, gui, n);
            end 
        end % function interfaceAnalysis
        %------------------------------------------------------------------        
        function clst = getNearestInterface(obj, bdRec, bdRef, i)
            % This function is applied to the query interface to find the
            % nearest reference interface. This is determined by comparing
            % the distance between the centre of the region bound by the
            % reconstructed interface with the centre of each region bound
            % by the set of reference interfaces.            
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   i - integer - index of interface from query reconstructed region.               
            % Outputs:
            %   clst.dmin - float - distance between the centre of region bound by the ith reconstructed interface from query region and the nearest centre of the reference region.
            %   clst.idx  - integer - index of the nearest reference region.
            %--------------------------------------------------------------
            for j = 1:length(bdRef)
                dc(j) = norm(bdRec{i}.centre - bdRef{j}.centre,2);
            end
            [clst.dmin, clst.idx] = min(dc);            
        end % function getNearestInterface
        %------------------------------------------------------------------
        function displayInterface(obj, fe, fldnmfwd, fldnmInv, bdRef, bdRec, gui, n)
            % This function plots reference and reconstructed interfaces along with
            % a subset of points on the reconstructed interface. The distance
            % between select reference and reconstructed interface are plotted
            % next to the corresponding point. General statistics related to the error
            % distribution are provided. 
            % Inputs:
            %   fe - finite element object - structure stores the forward and inverse models.
            %   fldnmfwd - string - field name of the reference model object
            %   fldnmInv - string - field name of the reconstructed model object
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   gui - GUI object - holds the data path, figure format.
            %   n - integer - index of the nth reconstructed interface of query region being analyzed.
            %--------------------------------------------------------------
            scTr = 0.3; % transparency scale (0 - transparent, 1 - no transparency)
            scFt = 1e3; % 1 - m, 1e-3 - mm
            offset = 8e-3*scFt;
            imgType = obj.getRecType(fldnmInv);
            xNodes = fe.xNodes*scFt;
            yNodes = fe.yNodes*scFt;            
            Xrec = [xNodes(bdRec{n}.pts(:,2))' yNodes(bdRec{n}.pts(:,1))'];
            
            figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.90, 0.70]);            
            % Display reference and reconstructed interfaces onto the
            % forward model            
            subplot(1,3,1);
            imagesc(imrotate(real(fe.fwdMdl.complexPermittivity),90),'XData',xNodes,'YData',-yNodes, 'AlphaData', scTr);hold on;
            colormap(gca, gui.figColorMap);            
            set(gca,'fontsize',14);            
            set(gca,'YDir','normal');
            % Go through all of the reference interfaces and plot all those
            % interfaces that points of the recontructed interface are
            % within the area bound by the reference interface.
            flag = 1; % The query interface is on or within at least one reference interface
            for i = 1:length(bdRef)
                %Only display the reference interface if the reconstructed
                %interface is on the reference interface or within the area
                %bound by the reference interface.
                Xref = [xNodes(bdRef{i}.pts(:,2))' yNodes(bdRef{i}.pts(:,1))'];                
                if obj.isRecOnOrtInRef(Xref,Xrec)
                    plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
                    flag = 0;
                end
            end
            
            % However, if the query point is not on or within a region bound
            % by a reference interface, then plot the reference interface
            % nearest the reconstructed interface.
            if flag
                clst = obj.getNearestInterface(bdRec, bdRef, n);
                Xref = [xNodes(bdRef{clst.idx}.pts(:,2))' yNodes(bdRef{clst.idx}.pts(:,1))'];
                plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
            end
            pRec = plot(Xrec(:,2), Xrec(:,1), '-b','LineWidth',4);hold on 
            xlim([xNodes(1) xNodes(end)])
            ylim([yNodes(1) yNodes(end)])
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);            
            xlabel([(fldnmfwd),' ', imgType, ' interface No. ', num2str(n), ' x (mm)']);            
            ylabel('y (mm)');
            ylabh = get(gca,'yLabel');
            set(ylabh,'Position',get(ylabh,'Position') + [7.5 0 0]);  
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[0.45 1.9 1.35 0.75]) % stretch its width and height                        
            drawnow     
            
            % Display distance between points on the reference and reconstructed 
            % along with the reference and reconstructed interfaces.
            subplot(1,3,2);
            error = obj.interfaces.(fldnmInv).distances{n}.error*scFt; 
            flag = 1;
            plot(Xrec(:,2), Xrec(:,1),'b','LineWidth',3);hold on 
            for i = 1:length(bdRef)
                Xref = [xNodes(bdRef{i}.pts(:,2))' yNodes(bdRef{i}.pts(:,1))'];
                %Only display the reference interface if the reconstructed
                %interface is on the reference interface or within the area
                %bound by the reference interface.
                if obj.isRecOnOrtInRef(Xref,Xrec)
                    plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
                    flag = 0;                    
                end
            end
            % However, if the query point is not on or within a region bound
            % by a reference interface, then plot the reference interface
            % nearest the reconstructed interface.
            if flag
                clst = obj.getNearestInterface(bdRec, bdRef, n);
                Xref = [xNodes(bdRef{clst.idx}.pts(:,2))' yNodes(bdRef{clst.idx}.pts(:,1))'];
                plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
            end
            
            centre = [(bdRec{n}.centre(:,1))' (bdRec{n}.centre(:,2))']*scFt;
                           
            [XrecOffset textAlgn theta] = obj.getOffsetPts(Xrec, centre, offset);
            for j = 1:size(Xrec,1)
                if ~mod(j,9)
                    plot(Xrec(j,2),Xrec(j,1),'bo','MarkerSize',12);hold on;
                    plot(Xrec(j,2),Xrec(j,1),'bo','MarkerSize',5);hold on;
                    plot(Xrec(j,2),Xrec(j,1),'k*','MarkerSize',5);hold on;
                    if strcmp(textAlgn(j),'R')
                        text(XrecOffset(j,2),XrecOffset(j,1),num2str(error(j),'%4.1f'),'Color','black','FontSize',11,'HorizontalAlignment','right');
                    else
                        if strcmp(textAlgn(j),'C')
                            text(XrecOffset(j,2),XrecOffset(j,1),num2str(error(j),'%4.1f'),'Color','black','FontSize',11,'HorizontalAlignment','center');
                        else
                            text(XrecOffset(j,2),XrecOffset(j,1),num2str(error(j),'%4.1f'),'Color','black','FontSize',11,'HorizontalAlignment','left');
                        end
                    end   
                end
            end
            
            set(gca, 'yticklabel','');
            xlabel('x (mm)','FontSize',14);
            %ht = text(-53, 63,'Forward model - red line','Color','red');
            %ht = text(-53, 59,'Reconstruction model - blue line','Color', 'blue');
            lh = legend('Box','off');
            lh.NumColumns = 2;
            lh.FontSize = 10;
            line([50,60], [75, 85], 'linewidth',3,'Color','red');    % Forward model
            line([50,60], [75, 85], 'linewidth',3,'Color','blue');  % Im reconstruction 
            xlim([xNodes(1) xNodes(end)])
            ylim([yNodes(1) yNodes(end)])
            set(gcf,'color','w');
            grid on;
            set(gca,'fontsize',14);
            legend({'Reconstruction model','Forward model'});
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[2.8 1.9 1.35 0.75]) % stretch its width and height                        
            drawnow
            
            % Plot histogram of error distribution (i.e., deviation between points on the 
            % estimated interface and the reference interface. Display general
            % statistics related to the error distribution.
            subplot(1,3,3);
            d = (obj.interfaces.(fldnmInv).distances{n}.error)*1e3;
            h = histogram(d,200, 'FaceColor', 'red');            
            maxh = max(h.BinCounts);
            maxd = max(d);
            mind = min(d);
            maxAbs = max([abs(maxd) abs(mind) ]);   
            ylim([0 (maxh + maxh*0.1)])
            xlim([-(maxAbs + maxAbs*0.20) (maxAbs + maxAbs*0.20)]);
%             text(mind,ceil((maxh*95)/100),['\mu =' num2str(mean(d),'%4.2f')],'Color','black','FontSize',14);
%             text(mind,ceil((maxh*87)/100),['\sigma =' num2str(std(d),'%4.2f')],'Color','black','FontSize',14);
            text((-maxAbs*0.8),((maxh + maxh*0.1)*95)/100,['\mu =' num2str(mean(d),'%4.2f')],'Color','black','FontSize',14);
            text((-maxAbs*0.8),((maxh + maxh*0.1)*87)/100,['\sigma =' num2str(std(d),'%4.2f')],'Color','black','FontSize',14);
            set(gcf,'color','w');
            grid on;
            ylabel('number of instances','FontSize',14);
            ylabh = get(gca,'yLabel');
            %set(ylabh,'Position',get(ylabh,'Position') + [0.4 0 0]);
            ylabh.Position(1) = ylabh.Position(1) + abs(ylabh.Position(1)*0.02);
            xlabel('error (mm)','FontSize',14);
            set(gca,'FontSize',14)
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[5.4 2.0 1.3 0.70]) % stretch its width and height                        
            drawnow
            
            % Save figures in format selected by user
            fileNm = ['tissueInterfaceImages_',(fldnmInv)];
            obj.saveFigure(fileNm, gui.dataPath, gui.figFormat);            
        end % function imgInterface        
        %------------------------------------------------------------------
        function obj = geometricPropertyAnalysis(obj, fe, flnmfwd, flnmInv,  refModel, recModel, gui )
            % This function evaluates individual sub-regions of the glandular and malignant tissue
            % formed from isolated groups of reconstructed tissue that are
            % bound by contours extracted during the tissue interface analysis
            % task. Each contour and associated region extracted are compared
            % with nearest reference contour and region.        
            % Inputs:
            %   fe - finite element object - structure stores the forward and inverse models.
            %   fldnmfwd - string - field name of the reference model object
            %   fldnmInv - string - field name of the reconstructed model object
            %   refModel - segmentation object - holds the reference model masks, interfaces, and tissue types.
            %   recModel - segmentation object - holds the  reconstructed model masks, interfaces, and tissue types.
            %   gui - GUI object - holds the data path, figure format.            
            % Outputs:
            %   obj.geometric.(flnmInv).xcorrGeo   - structure - normalized cross-correlation between query reconstructed mask and reference mask.
            %   obj.geometric.(flnmInv).goodnessOfFit - structure - goodness of fit between the reference and reconstructed masks 
            %   obj.geometric.(flnmInv).similarity  - structure - calculates the Jaccard similarity coefficient, Dice similarity coefficient, ratio of detection, artefact rejection ratio, specificity, and precision
            %   obj.geometric.(flnmInv).avgHd - structure - average Hausdorff distance between query reconstructed interface and corresponding reference interface.
            %--------------------------------------------------------------            
            
            % Pre-processing steps
            % First, get the interior and skin masks to form a region of
            % interest for which the glandular and tumor regions are
            % located.
            ROI = recModel.masks.interior - recModel.masks.skin;
            
            % Second, get the foward and corresponding reconstruction (or
            % inverse) models tissue maps within the ROI. These allow the
            % regional metrics to correct for malignant tissue emmbedded in
            % glandular tissue.
            
            % Get the tissue map for the forward model within the ROI
            tissueMapFwd = refModel.tissueTypes.fat + refModel.tissueTypes.transition + refModel.tissueTypes.fibroglandular + refModel.tissueTypes.tumor;
            
            % Get the reference mask over the breast interior
            refMask = refModel.masks.(flnmfwd);    
            
            % Next, get the tissue map for the inverse model within the ROI
            imgType = obj.getRecType(flnmInv);
            switch imgType
                case 'Real'
                    fwdMdl = real(fe.fwdMdl.complexPermittivity);
                    invMdl = real(fe.invMdl.complexPermittivity);
                    tissueMapInv = recModel.tissueTypes.fatRe + recModel.tissueTypes.transitionRe + recModel.tissueTypes.fibroglandularRe + recModel.tissueTypes.tumorRe;                                   
                    
                case 'Imag'
                    fwdMdl = -imag(fe.fwdMdl.complexPermittivity);
                    invMdl = -imag(fe.invMdl.complexPermittivity);
                    tissueMapInv = recModel.tissueTypes.fatIm + recModel.tissueTypes.transitionIm + recModel.tissueTypes.fibroglandularIm + recModel.tissueTypes.tumorIm;
                    
                otherwise
                    fwdMdl = abs(fe.fwdMdl.complexPermittivity);
                    invMdl = abs(fe.invMdl.complexPermittivity);
                    tissueMapInv = recModel.tissueTypes.fatMag + recModel.tissueTypes.transitionMag + recModel.tissueTypes.fibroglandularMag + recModel.tissueTypes.tumorMag;
                    
            end
            
            % Get all of the reference and reconstructed interfaces for the 
            % query region (fldnm).
            bdRef = refModel.interfaces.(flnmfwd);
            bdRec = recModel.interfaces.(flnmInv);
            
            obj.geometric.(flnmInv).xcorrGeo = obj.xCorrRegionGeo(bdRef, bdRec, ROI, flnmfwd, tissueMapFwd, tissueMapInv, fwdMdl, invMdl, refMask);
            obj.geometric.(flnmInv).goodnessOfFit = obj.goodnessOfFitGeo(bdRef, bdRec, ROI, flnmfwd, tissueMapFwd, tissueMapInv, fwdMdl, invMdl, refMask);
            obj.geometric.(flnmInv).similarity = obj.regionSimilarity(bdRef, bdRec, ROI, flnmfwd, tissueMapFwd, tissueMapInv, fwdMdl, invMdl, refMask);
            obj.geometric.(flnmInv).avgHd = obj.avgHausdorffDist(bdRef, bdRec);            
            obj.displayGeoPropAnalysis(bdRef, bdRec, ROI, flnmfwd, flnmInv, fwdMdl, invMdl, tissueMapFwd, tissueMapInv, refMask, gui.dataPath, gui.figFormat);
        end % function geometricPropertyAnalysis        
        %------------------------------------------------------------------
        function xcorrGeo = xCorrRegionGeo(obj, bdRef, bdRec, ROI, tissueType, tissueMapFwd, tissueMapInv, fwdMdl, invMdl, refMask)
            % This function evaluates the normalized cross-correlation
            % between query reconstructed mask and reference mask.            
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   tissueMapFwd - n by m tissue map of formed by segmenting the forward model.
            %   tissueMapInv - n by m tissue map formed by segmenting the inverse model.
            %   fwdMdl - n by m query component of forward model complex permittivity map
            %   invMdl - n by m query component of inverse model complex permittivity map
            %   refMask - n by m query reference mask.
            % Outputs:
            %   xcorrGeo   - structure - calculated normalized cross-correlation between query reconstructed mask and reference mask.            
            %-------------------------------------------------------------
            
            xcorrGeo = [];   
            for i = 1:length(bdRec)
                roi = obj.getROIMasks(i, bdRec, bdRef, ROI, tissueMapInv, tissueMapFwd, tissueType, fwdMdl, invMdl, refMask );                
                % Vectorize reconstructed and reference masks that have
                % been constructed within the ROI, and evaluate the
                % normalized cross correlation of the two vectors.
                n = size(roi.recMask_i,1);
                m = size(roi.recMask_i,2);
                
                recVector = reshape(roi.recMask_i,n*m,1); 
                if strcmp(class(recVector), 'logical')
                    recVector = double(recVector);
                end
                refVector = reshape(roi.refMask_j,n*m,1);
                if strcmp(class(refVector), 'logical')
                    refVector = double(refVector);
                end
                xcorrGeo{i} = (refVector'*recVector)/( norm(refVector,2)*norm(recVector,2) );
            end % repeat for all reconstructed contours related to the region
        end % function xCorrRegionGeo
        %------------------------------------------------------------------        
        function goodnessOfFit = goodnessOfFitGeo(obj, bdRef, bdRec, ROI, tissueType, tissueMapFwd, tissueMapInv, fwdMdl, invMdl, refMask)
            % This function evaluates the goodness of fit between the reference and
            % reconstructed masks is evaluated with the normalized root mean
            % square error (NRMSE) cost function.
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   tissueMapFwd - n by m tissue map of formed by segmenting the forward model.
            %   tissueMapInv - n by m tissue map formed by segmenting the inverse model.
            %   fwdMdl - n by m query component of forward model complex permittivity map
            %   invMdl - n by m query component of inverse model complex permittivity map
            %   refMask - n by m query reference mask.
            % Outputs:
            %   goodnessOfFit   - structure - calculated he goodness of fit between the reference and reconstructed masks.
            %-------------------------------------------------------------
            goodnessOfFit = [];            
            for i = 1:length(bdRec)
                roi = obj.getROIMasks(i, bdRec, bdRef, ROI, tissueMapInv, tissueMapFwd, tissueType, fwdMdl, invMdl, refMask );                
                % Vectorize reconstructed and reference masks that have
                % been constructed within the ROI, and evaluate the
                % normalized cross correlation of the two vectors.
                n = size(roi.recMask_i,1);
                m = size(roi.recMask_i,2);              
                
                refVector = reshape(roi.refMask_j,n*m,1);
                if strcmp(class(refVector), 'logical')
                    refVector = double(refVector);
                end
               goodnessOfFit{i} = 1 - (norm( roi.refMask_j - roi.recMask_i) / norm(roi.refMask_j - mean(refVector)*ones(n,m) )) ;
            end % repeat for all reconstructed contours related to the region          
        end % function goodnessOfFitGeo        
        %------------------------------------------------------------------
        function similarity = regionSimilarity(obj, bdRef, bdRec, ROI, tissueType, tissueMapFwd, tissueMapInv, fwdMdl, invMdl, refMask)
            % This function evaluates an array of metrics to measure the
            % similarity between thethe goodness of fit between the query
            % reconstructed mask and corresponding reference mask.            
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   tissueMapFwd - n by m tissue map of formed by segmenting the forward model.
            %   tissueMapInv - n by m tissue map formed by segmenting the inverse model.
            %   fwdMdl - n by m query component of forward model complex permittivity map
            %   invMdl - n by m query component of inverse model complex permittivity map
            %   refMask - n by m query reference mask.
            % Outputs:
            %   similarity{i}.Jaccard  - float - Jaccard similarity coefficient is defined as the size of the intersection divided by the size of the union of the regions 
            %   similarity{i}.DiceCoef - float - Dice similarity coefficient is defined as the size of the intersection divided by the average size of the reconstructed and reference regions
            %   similarity{i}.RD - float - ratio of detection (RD) metric measures the proportion of the reference mask that has been correctly reconstructed
            %   similarity{i}.AR - float - artefact rejection ratio (AR) metric measures the proportion of tissue incorrectly reconstructed as reference tissue outside the reference region
            %   similarity{i}.Specificity - 
            %   similarity{i}.Precision - float - the proportion of the reconstructed mask that has been correctly reconstructed as the reference tissue            
            %-------------------------------------------------------------     
            similarity = [];            
            for i = 1:length(bdRec)
                roi = obj.getROIMasks(i, bdRec, bdRef, ROI, tissueMapInv, tissueMapFwd, tissueType, fwdMdl, invMdl, refMask );
                            
                m = roi.refMask_j(:); % m is a ground truth mask
                o = roi.recMask_i(:);  % o is the same tissue in segmented image          
                common=sum(m & o); % common - represents true positive (TP) pixels
                union=sum(m | o);            
                cm=sum(m); % the number of voxels in m
                co=sum(o); % the number of voxels in o
                t = size(roi.refMask_j,1)*size(roi.refMask_j,2);  % t - total number of pixels in the region of interest               
                ts = sum(t);
                tn = ts - union; %tn - represents the true negative (TN) pixels - all pixels in the region of interest that are not in the reference or reconstructed masks
                fp = co - common; % fp - represents the false positive (FP) pixels
                similarity{i}.Jaccard=common/union;
                similarity{i}.DiceCoef=(2*common)/(cm+co);
                similarity{i}.RD = common/cm; % ratio of detection - proportion of the reference mask that has been correctly reconstructed. 
                similarity{i}.AR = 1 - (co-common)/cm; % artifact rejection ratio (AR) - proportion of tissue incorrectly reconstructed as reference tissue outside the reference region 
                similarity{i}.Specificity = tn/(tn+fp);
                similarity{i}.Precision = common/co; % precision - proportion of the reconstructed mask that has been correctly reconstructed. 
            end % repeat for all reconstructed contours related to the region fldnm
        end % function regionSimilarity
        %------------------------------------------------------------------
        function avgHd = avgHausdorffDist(obj, bdRef, bdRec)
            % This function evaluates the average Hausdorff distance between edge points
            % of reconstructed sub-glandular and sub-malignant tissue regions
            % and the corresponding reference sub-region.             
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.            
            % Outputs:
            %   avgHd - float - The Hausdorff distance measures the degree of mismatch between the two sets by measuring the distance of the point on the reconstructed interface that is furthest from any point on the reference interface and visa versa.
            %-------------------------------------------------------------                
            avgHd = [];            
            for i = 1:length(bdRec)
                % Step 1 - for each reconstructed contour, find the nearest
                % corresponding reference region.
                Xrec_i = [obj.xNodes(bdRec{i}.pts(:,2))' obj.yNodes(bdRec{i}.pts(:,1))'];                
                avgHd_j = [];
                for j = 1:length(bdRef)
                    Xref_j = [obj.xNodes(bdRef{j}.pts(:,2))' obj.yNodes(bdRef{j}.pts(:,1))'];                    
                    if obj.isRecOnOrtInRef(Xref_j,Xrec_i)
                        avgHd_j(j) = obj.getAvgHausdorffDist(Xrec_i, Xref_j );                        
                    end
                end
                avgHd{i} = max(avgHd_j);
            end
         end % function avgHausdorffDist
         %-----------------------------------------------------------------
         function hd = HausdorffDist(obj, bdRef, bdRec)
            % Function evaluates the Hausdorff distance between edge points
            % of reconstructed sub-glandular and sub-malignant tissue regions
            % and the corresponding reference sub-region.
            hd = [];            
            for i = 1:length(bdRec)
                % Step 1 - for each reconstructed contour, find the nearest
                % corresponding reference region.
                Xrec_i = [obj.xNodes(bdRec{i}.pts(:,2))' obj.yNodes(bdRec{i}.pts(:,1))'];                
                avgHd_j = [];
                for j = 1:length(bdRef)
                    Xref_j = [obj.xNodes(bdRef{j}.pts(:,2))' obj.yNodes(bdRef{j}.pts(:,1))'];                    
                    if obj.isRecOnOrtInRef(Xref_j,Xrec_i)
                        hd_j(j) = obj.getHausdorffDist(Xrec_i, Xref_j );                       
                    end
                end
                try
                    hd{i} = max(hd_j);
                catch
                    hd(i) = -1;
                end
            end            
         end % function HausdorffDist                  
         %------------------------------------------------------------------
         function avgHd = getAvgHausdorffDist(obj,  A, B )
            % This function computes the average Hausdorff distance between edge
            % points of a reconstructed region (A) and the edge points of the
            % corresponding reference region (B).
            % The average Forward Hausdorff Distance (fhd) is defined as:
            % fhd(A,B) = 1/Na sum a in A [ min b in B [ ||a-b|| ] ].            
            % Likewise, the Reverse Hausdorff Distance is defined as:
            % rhd(B,A) = 1/Nb sum b in B [ min a in A [ ||a-b|| ] ]
            % The Hausdorff Distance is defined as max{fhd(A,B),rhd(B,A)}
            % Inputs:
            %   A -array - interface from the query reconstructed region.
            %   B -array - interface from the corresponding reference region.            
            % Outputs:
            %   avgHd - float - The Hausdorff distance measures the degree of mismatch between the two sets by measuring the distance of the point on the reconstructed interface that is furthest from any point on the reference interface and visa versa.
            %--------------------------------------------------------------
            % Compute the sizes of the input point sets
            Asize = size(A);
            Bsize = size(B);

            % Check if the points have the same dimensions
            if Asize(2) ~= Bsize(2)
                error('The dimensions of points in the two sets are not equal');
            end

            % Evaluate the Hausdorff distance in the forward direction
            fhd = 0;                          % Initialize forward distance to 0
            for a = 1:Asize(1)             % Travel the set A to find avg of d(A,B)
                mindist = Inf;              % Initialize minimum distance to Inf
                for b = 1:Bsize(1)         % Travel set B to find the min(d(a,B))
                    tempdist = norm(A(a,:)-B(b,:));
                    if tempdist < mindist
                        mindist = tempdist;
                    end
                end
                fhd = fhd + mindist;        % Sum the forward distances
            end
            fhd = fhd/Asize(1);             %  average forward Hausdorff distance - fhd(A,B) = 1/Na sum a in A [ min b in B [ ||a-b|| ] ]

            % Now, evaluate the reverse Hausdorff distance

            rhd = 0;                    % Initialize reverse distance to 0
            for b = 1:Bsize(1)          % Travel the set B to find avg of d(B,A)
                mindist = Inf;          % Initialize minimum distance to Inf
                for a = 1:Asize(1)      % Travel set A to find the min(d(b,A))
                    tempdist = norm(A(a,:)-B(b,:));
                    if tempdist < mindist
                        mindist = tempdist;
                    end
                end
                rhd = rhd + mindist;    % Sum the reverse distances
            end
            rhd = rhd/Bsize(1);           % average reverse Hausdorff distance - rhd(B,A) = 1/Nb sum b in B [ min a in A [ ||a-b|| ] ]

            avgHd = max(fhd,rhd);         % The average Hausdorff Distance - max{fhd(A,B),rhd(B,A)}                                                  
         end % function avgHausdorffDist         
         %------------------------------------------------------------------
         function hd = getHausdorffDist(obj,  A, B )
            % The function computes the Hausdorff distance between edge
            % points of a reconstructed region (A) and the edge points of the
            % corresponding reference region (B).
            % The Forward Hausdorff Distance (fhd) is defined as:
            % fhd(A,B) = max a in A [ min b in B [ ||a-b|| ] ].
            % Intuitively fhd finds the point a from the set of edge points on the reconstructed region
            % A that is farthest from any point in the set of edge points
            % on the reference region B and measures the distance from a to its nearest neighbor
            % in B. Likewise, the Reverse Hausdorff Distance is defined as:
            % rhd(B,A) = max b in B [ min a in A [ ||a-b|| ] ]
            % The Hausdorff Distance is defined as max{fhd(A,B),rhd(B,A)}
            %--------------------------------------------------------------
            % Compute the sizes of the input point sets
            Asize = size(A);
            Bsize = size(B);

            % Check if the points have the same dimensions
            if Asize(2) ~= Bsize(2)
                error('The dimensions of points in the two sets are not equal');
            end

            % Evaluate the Hausdorff distance in the forward direction
            fhdMindist = zeros(Asize(1),1);  % Initialize sest of forward distance to 0
            for a = 1:Asize(1)                     % Travel the set A to find fhd of d(A,B)
                mindist = Inf;                      % Initialize minimum distance to Inf
                for b = 1:Bsize(1)                 % Travel set B to find the min(d(a,B))
                    tempdist = norm(A(a,:)-B(b,:));
                    if tempdist < mindist
                        mindist = tempdist;
                    end
                end
                fhdMindist(a) = mindist;        % forward distance for point ai
            end
            fhd = max(fhdMindist);             % fhd(A,B) = max a in A [ min b in B [ ||a-b|| ] ]

            % Now, evaluate the Reverse Hausdorff distance

            rhdMindist = zeros(Bsize(1),1);
            for b = 1:Bsize(1)                  % Travel the set B to find avg of d(B,A)
                mindist = Inf;                   % Initialize minimum distance to Inf
                for a = 1:Asize(1)              % Travel set A to find the min(d(b,A))
                    tempdist = norm(A(a,:)-B(b,:));
                    if tempdist < mindist
                        mindist = tempdist;
                    end
                end
                rhdMindist(b) = mindist;        % reverse distance for point bj
            end
            rhd = max(rhdMindist);          % Reverse Hausdorff distance - rhd(B,A) = max b in B [ min a in A [ ||a-b|| ] ]

            hd = max(fhd,rhd);         % The Hausdorff Distance - max{fhd(A,B),rhd(B,A)}
         end % function HausdorffDist
        %------------------------------------------------------------------
        function displayGeoPropAnalysis(obj, bdRef, bdRec, ROI, tissueType, flnminv, fwdMdl, invMdl, tissueMapFwd, tissueMapInv, refMask, dataPath, figFormat)
            % This function displays: (1) the query reconstructed mask, (2) relevent regions of forward
            % model that the reconstructed region is compared with, (3) interface of
            % query reconstructed masks and reference interface that the interface
            % is compared with superimposed onto the foreard model.             
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   flnminv - string - reconstructed region, e.g., glandularRe, tumorRe, glandularIm, tumorIm, glandularMag, tumorMag.           
            %   fwdMdl - n by m query component of forward model complex permittivity map
            %   invMdl - n by m query component of inverse model complex permittivity map
            %   tissueMapFwd - n by m tissue map of segmented forward model
            %   tissueMapInv - n by m tissue map of segmented inverse model
            %   refMask - n by m query reference mask.
            %   dataPath - string - points to directory that figure folder to save figure is located.
            %   figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------            
            
            for i = 1:length(bdRec)
                % Find the nearest reference contour to the
                % reconstructed contour. Form masks from the region bound by
                % these contours and form a bounding box around the masks
                % and extract the masks from the bounding box. Apply thses
                % masks to the forward and reconstructed models to extract
                % regions of interest. The metric is applied to this region
                % of interest tissue.
                roi = obj.getROIMasks( i, bdRec, bdRef, ROI, tissueMapInv, tissueMapFwd, tissueType, fwdMdl, invMdl, refMask );                                   
                imgType = obj.getRecType(flnminv);
                
                scTr = 0.3; % transparency scale (0 - transparent, 1 - no transparency)
                scFt = 1; % 1 - m, 1e-3 - mm                              
                xNodes = obj.xNodes*scFt;
                yNodes = obj.yNodes*scFt;
                Xrec_i = [xNodes(bdRec{i}.pts(:,2))' yNodes(bdRec{i}.pts(:,1))'];
                
                fh = figure;                               
                set(fh, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.80, 0.70]);
%                 set(fh, 'MenuBar','none');
                set(fh, 'ToolBar','none');
            
                subplot(2,3,1);
                imagesc(imrotate(roi.recMask_i,90),'XData',xNodes,'YData',-yNodes);hold on;                           
                set(gca,'fontsize',14);            
                set(gca,'YDir','normal');
                set(gca, 'yticklabel','');
                set(gca, 'xtickLabel','');
                xlabel(['Recon. mask ', tissueType ,' ', imgType]); 
                p = get(gca,'position'); % get subplot axis position
                set(gca,'position',p.*[0.4 0.90 1.3 1.3]) % stretch its width and height                        
                drawnow 
                
                subplot(2,3,2);
                imagesc(imrotate(roi.refMask_j,90),'XData',xNodes,'YData',-yNodes);hold on;                                                                     
                set(gca, 'yticklabel','');
                set(gca, 'xtickLabel','');
                xlabel(['Forward mask ', tissueType ,' ', imgType]); 
                set(gca,'fontsize',14);            
                set(gca,'YDir','normal');
                p = get(gca,'position'); % get subplot axis position
                set(gca,'position',p.*[0.9 0.90 1.3 1.3]) % stretch its width and height 
                drawnow                
                
                pf = [obj.geometric.(flnminv).xcorrGeo{i};
                  obj.geometric.(flnminv).goodnessOfFit{i};
                  obj.geometric.(flnminv).similarity{i}.Jaccard;
                  obj.geometric.(flnminv).similarity{i}.DiceCoef;
                  obj.geometric.(flnminv).similarity{i}.RD;
                  obj.geometric.(flnminv).similarity{i}.AR;
                  obj.geometric.(flnminv).similarity{i}.Specificity;
                  obj.geometric.(flnminv).similarity{i}.Precision;
                  obj.geometric.(flnminv).avgHd{i}*1e3];
                try
                    pfstr = {num2str(pf(1),'%8.4f');
                    num2str(pf(2),'%8.4f');
                    num2str(pf(3),'%8.4f');
                    num2str(pf(4),'%8.4f');
                    num2str(pf(5),'%8.4f');
                    num2str(pf(6),'%8.4f');
                    num2str(pf(7),'%8.4f');
                    num2str(pf(8),'%8.4f');
                    num2str(pf(9),'%8.4f')};
                
                    ha = subplot(2,3,3);
                    cnames = {'Value'};
                    rnames = {'Xcorr', 'GOF', 'Jaccard', 'DiceCoef','Ratio-of-Detection','artefact rejection ratio','Specificity','Precision','avg Hd (mm)'};
                    ht = uitable('ColumnName', cnames, 'RowName', rnames, 'Data', pfstr,'ColumnWidth', {99});                 
                    hh = uicontrol('Style','text','Position',[750 590 400 20], 'String', 'Geometric Analysis');               
                    set(hh, 'FontSize',13);
                    pos = get(ha,'Position');                
                    delete(ha);
                    set(ht, 'units', 'normalized');
                    set(ht, 'position', pos);
                    set(ht,'ColumnName',cnames);
                    set(ht,'RowName', rnames);
                    p = get(ht,'position'); % get subplot axis position
                    set(ht,'position',p.*[0.95 1.05 1.533 0.90]) % stretch its width and height 
                catch
                     ha = subplot(2,3,3);
                     hh1 = uicontrol('Style','text','Position',[750 590 400 20], 'String', 'Geometric Analysis');
                     set(hh1, 'FontSize',13);
                     hh2 = uicontrol('Style','text','Position',[740 570 400 20], 'String', 'Reconstructed and forward regions do not overlap');
                     set(hh2, 'FontSize',13);
                     delete(ha);
                end
                drawnow
                
                subplot(2,3,6);
                imagesc(imrotate(fwdMdl,90),'XData',xNodes,'YData',-yNodes, 'AlphaData', scTr);hold on;
                colormap(gca, 'jet');            
                set(gca,'fontsize',14);                 
                set(gca,'YDir','normal'); 
                lh = legend('Box','off');
                lh.NumColumns = 2;
                lh.FontSize = 10;                
                plot(Xrec_i(:,2), Xrec_i(:,1), '-b','LineWidth',4);hold on 
                for j = 1:length(bdRef)
                    Xref_j = [obj.xNodes(bdRef{j}.pts(:,2))' obj.yNodes(bdRef{j}.pts(:,1))'];                    
                    if obj.isRecOnOrtInRef(Xref_j,Xrec_i)                        
                        plot(Xref_j(:,2), Xref_j(:,1), '-r','LineWidth',4);hold on
                    end
                end
                
                flag = 1; % The query interface is on or within at least one reference interface
                for j = 1:length(bdRef)
                    %Only display the reference interface if the reconstructed
                    %interface is on the reference interface or within the area
                    %bound by the reference interface.
                    Xref_j = [xNodes(bdRef{j}.pts(:,2))' yNodes(bdRef{j}.pts(:,1))'];                
                    if obj.isRecOnOrtInRef(Xref_j, Xrec_i)
                        plot(Xref_j(:,2), Xref_j(:,1), '-r','LineWidth',4);hold on
                        flag = 0;
                    end
                end

                % However, if the query reconstructed interface is not on or
                % within a region bound by a reference interface, then plot
                % the reference interface nearest to the reconstructed interface.
                if flag
                    clst = obj.getNearestInterface(bdRec, bdRef, i);
                    Xref = [xNodes(bdRef{clst.idx}.pts(:,2))' yNodes(bdRef{clst.idx}.pts(:,1))'];
                    plot(Xref(:,2), Xref(:,1), '-r','LineWidth',4);hold on
                end                
                xlim([xNodes(1) xNodes(end)])
                ylim([yNodes(1) yNodes(end)])
                p=get(gca,'position'); % save position            
                set(gca,'fontsize',14);            
                xlabel([tissueType,' ', imgType, ' interface No. ', num2str(i), ' (m)']);            
                ylabel('y (m)');                
                ylabh = get(gca,'yLabel');
                ylabh.Position(1) = ylabh.Position(1) + abs(ylabh.Position(1)*0.15);
                set(gca,'position',p); % restore position
                legend({'Rec','Fwd'});
                sub_pos = get(gca,'position'); % get subplot axis position
                set(gca,'position',p.*[1.02 0.97 1.2 1.3]) % stretch its width and height  
                drawnow
                
                % Save figures in format selected by user
                fileNm = ['geometricAnalysisImages_' flnminv ,'_' num2str(i)];
                obj.saveFigure(fileNm, dataPath, figFormat);  
            end
        end
         %-----------------------------------------------------------------
         function obj = dielectricPropertyAnalysis(obj, fe, flnmfwd, flnmInv, refModel, recModel, gui )             
            % This function applies metrics to segmented regions to evaluate
            % the accuracy with which both the geometric and dielectric properties
            % of these underlying structures are reconstructed. 
            % Each reference mask is applied to the forward model to extract
            % the region of the model that corresponds to the tissue group
            % represented by the mask (e.g., tumor region). These segmented
            % property values are referred to as the reference tissue,
            % of the region. The process is repeated for each mask, so for each
            % reference mask there is a corresponding reference tissue. Likewise,
            % the reconstructed masks are applied to the reconstructed images.
            % These segmented property values are referred to as the reconstructed tissue,
            % of the region.
            % Inputs:
            %   fe - finite element object - structure stores the forward and inverse models.
            %   fldnmfwd - string - field name of the reference model object
            %   fldnmInv - string - field name of the reconstructed model object
            %   refModel - segmentation object - holds the reference model masks, interfaces, and tissue types.
            %   recModel - segmentation object - holds the  reconstructed model masks, interfaces, and tissue types.
            %   gui - GUI object - holds the data path, figure format.            
            % Outputs:
            %   obj.dielectric.(flnmInv).xcorrGeo   - structure - normalized cross-correlation between query reconstructed mask and reference mask.
            %   obj.dielectric.(flnmInv).goodnessOfFit - structure - goodness of fit between the reference and reconstructed masks             
            %--------------------------------------------------------------                
            
            % Pre-processing steps            
            % First get the foward and corresponding reconstruction (or
            % inverse) models
            imgType = obj.getRecType(flnmInv);  
            tissueMapFwd = refModel.tissueTypes.fat + refModel.tissueTypes.transition + refModel.tissueTypes.fibroglandular + refModel.tissueTypes.tumor;
            
            % Get the reference mask over the breast interior
            refMask = refModel.masks.(flnmfwd);   
                    
            switch imgType
                case 'Real'
                    fwdMdl = real(fe.fwdMdl.complexPermittivity);
                    invMdl = real(fe.invMdl.complexPermittivity);
                    tissueMapInv = recModel.tissueTypes.fatRe + recModel.tissueTypes.transitionRe + recModel.tissueTypes.fibroglandularRe + recModel.tissueTypes.tumorRe;
                    
                case 'Imag'
                    fwdMdl = -imag(fe.fwdMdl.complexPermittivity);
                    invMdl = -imag(fe.invMdl.complexPermittivity);
                    tissueMapInv = recModel.tissueTypes.fatIm + recModel.tissueTypes.transitionIm + recModel.tissueTypes.fibroglandularIm + recModel.tissueTypes.tumorIm;
                otherwise
                    fwdMdl = abs(fe.fwdMdl.complexPermittivity);
                    invMdl = abs(fe.invMdl.complexPermittivity);
                    tissueMapInv = recModel.tissueTypes.fatMag + recModel.tissueTypes.transitionMag + recModel.tissueTypes.fibroglandularMag + recModel.tissueTypes.tumorMag;
            end
            
            % Next, get all of the reference and reconstructed interfaces for the 
            % query region (fldnm).
            bdRef = refModel.interfaces.(flnmfwd);
            bdRec = recModel.interfaces.(flnmInv);
            
            % Finally, get the interior and skin masks to form a region of
            % interest for which the glandular and tumor regions are
            % located. This ROI will be refined by the methods.
            ROI = recModel.masks.interior - recModel.masks.skin;            
             
            obj.dielectric.(flnmInv).xcorrDiel = obj.xCorrRegionDiel( bdRef, bdRec, ROI, flnmfwd, fwdMdl, invMdl, tissueMapFwd, tissueMapInv, refMask);
            obj.dielectric.(flnmInv).GOFgeoDiel = obj.goodnessOfFitGeoDiel( bdRef, bdRec, ROI, flnmfwd, fwdMdl, invMdl, tissueMapFwd, tissueMapInv, refMask);
            obj.displayDielPropAnalysis(bdRef, bdRec, ROI, flnmfwd, flnmInv, fwdMdl, invMdl, tissueMapFwd, tissueMapInv, refMask, gui.dataPath, gui.figFormat);
            
        end % function evalRegionalGeometricDielectricPerformance        
        %------------------------------------------------------------------
        function xcorrDiel = xCorrRegionDiel( obj, bdRef, bdRec, ROI, tissueType, fwdMdl, invMdl, tissueMapFwd, tissueMapInv, refMask)
            % This function evaluates the goodness of fit between the reference and
            % reconstructed tissue is evaluated with the normalized root mean
            % square error (NRMSE) cost function.       
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   fwdMdl - n by m query component of forward model complex permittivity map
            %   invMdl - n by m query component of inverse model complex permittivity map            
            %   tissueMapFwd - n by m tissue map of formed by segmenting the forward model.
            %   tissueMapInv - n by m tissue map formed by segmenting the inverse model.           
            %   refMask - n by m query reference mask.
            % Outputs:
            %   xcorrDiel   - structure - calculated normalized cross-correlation between query reconstructed tissue and reference tissue.            
            %-------------------------------------------------------------            
            xcorrDiel = [];            
            for i = 1:length(bdRec)
                % Find the nearest reference contour to the
                % reconstructed contour. Form masks from the region bound by
                % these contours and form a bounding box around the masks
                % and extract the masks from the bounding box. Apply thses
                % masks to the forward and reconstructed models to extract
                % regions of interest. The metric is applied to this region
                % of interest tissue.
                roi = obj.getROIMasks( i, bdRec, bdRef, ROI, tissueMapInv, tissueMapFwd, tissueType, fwdMdl, invMdl, refMask );

                % Vectorize reconstructed and reference tissues that have
                % been constructed within the ROI, and evaluate the
                % normalized cross correlation of the two vectors.
                [n, m] = size(roi.invMdl_i);
                recVector = reshape(roi.invMdl_i,n*m,1);                
                refVector = reshape(roi.fwdMdl_j,n*m,1);
                xcorrDiel{i} = (refVector'*recVector)/( norm(refVector,2)*norm(recVector,2) );
            end      
        end % function xCorrRegionGeoDiel
        %------------------------------------------------------------------
        function GOFgeoDiel = goodnessOfFitGeoDiel( obj, bdRef, bdRec, ROI, tissueType, fwdMdl, invMdl, tissueMapFwd, tissueMapInv, refMask )
            % This function evaluates the normalized cross-correlation
            % between query reconstructed tissue and reference tissue.            
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   fwdMdl - n by m query component of forward model complex permittivity map
            %   invMdl - n by m query component of inverse model complex permittivity map            
            %   tissueMapFwd - n by m tissue map of formed by segmenting the forward model.
            %   tissueMapInv - n by m tissue map formed by segmenting the inverse model.           
            %   refMask - n by m query reference mask.
            % Outputs:
            %   goodnessOfFit   - structure - calculated he goodness of fit between the reference and reconstructed tissue.
            %-------------------------------------------------------------   
            GOFgeoDiel = [];            
            for i = 1:length(bdRec)
                % Find the nearest reference contour to the
                % reconstructed contour. Form masks from the region bound by
                % these contours and form a bounding box around the masks
                % and extract the masks from the bounding box. Apply thses
                % masks to the forward and reconstructed models to extract
                % regions of interest. The metric is applied to this region
                % of interest tissue.
                roi = obj.getROIMasks( i, bdRec, bdRef, ROI, tissueMapInv, tissueMapFwd, tissueType, fwdMdl, invMdl, refMask );
                
                [n, m] = size(roi.fwdMdl_j);                              
                refVector = reshape(roi.fwdMdl_j,n*m,1);
                GOFgeoDiel{i} = 1 - (norm( roi.fwdMdl_j - roi.invMdl_i) / norm(roi.fwdMdl_j - mean(refVector)*ones(n,m) )) ;                                
            end
        end % function goodnessOfFitGeoDiel
        %------------------------------------------------------------------
        function displayDielPropAnalysis(obj, bdRef, bdRec, ROI, tissueType, flnminv, fwdMdl, invMdl, tissueMapFwd, tissueMapInv, refMask, dataPath, figFormat)
            % This function displays: (1) the query reconstructed tissue, (2) relevent regions of forward
            % model that the reconstructed region is compared with, and (3)
            % a table of the metric values.                 
            % Inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRec -structure - interface from the corresponding reference region.
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   flnminv - string - reconstructed region, e.g., glandularRe, tumorRe, glandularIm, tumorIm, glandularMag, tumorMag.           
            %   fwdMdl - n by m query component of forward model complex permittivity map
            %   invMdl - n by m query component of inverse model complex permittivity map
            %   tissueMapFwd - n by m tissue map of segmented forward model
            %   tissueMapInv - n by m tissue map of segmented inverse model
            %   refMask - n by m query reference mask.
            %   dataPath - string - points to directory that figure folder to save figure is located.
            %   figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------          
            for i = 1:length(bdRec)
                % Find the nearest reference contour to the
                % reconstructed contour. Form masks from the region bound by
                % these contours and form a bounding box around the masks
                % and extract the masks from the bounding box. Apply thses
                % masks to the forward and reconstructed models to extract
                % regions of interest. The metric is applied to this region
                % of interest tissue.
                roi = obj.getROIMasks( i, bdRec, bdRef, ROI, tissueMapInv, tissueMapFwd, tissueType, fwdMdl, invMdl, refMask );
                maxValue = max(max(roi.fwdMdl_j));
                imgType = obj.getRecType(flnminv);
                
                fh = figure;                               
                set(fh, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.80, 0.50]);
                set(fh, 'MenuBar','none');
                %set(fh, 'ToolBar','none');
            
                subplot(1,3,1);
                imagesc(imrotate(roi.invMdl_i,90),'XData',obj.xNodes,'YData',-obj.yNodes);hold on;
                colormap(gca, 'jet');            
                set(gca,'fontsize',14);            
                set(gca,'YDir','normal');
                set(gca, 'yticklabel','');
                set(gca, 'xtickLabel','');                
                %ylabh.Position(1) = ylabh.Position(1) + abs(ylabh.Position(1)*0.15);
                xlabel(['Inverse model ', tissueType ,' ', imgType]); 
                p = get(gca,'position'); % get subplot axis position
                set(gca,'position',p.*[0.4 1.6 1.4 0.9]) % stretch its width and height                        
                drawnow 
                
                subplot(1,3,2);
                imagesc(imrotate(roi.fwdMdl_j,90),'XData',obj.xNodes,'YData',-obj.yNodes);hold on;
                colormap(gca, 'jet'); 
                ch = colorbar;
                caxis([0 maxValue]);
                set(gca, 'yticklabel','');
                set(gca, 'xtickLabel','');
                xlabel(['Forward model ', tissueType ,' ', imgType]); 
                set(gca,'fontsize',14);            
                set(gca,'YDir','normal');
                p = get(gca,'position'); % get subplot axis position
                set(gca,'position',p.*[0.9 1.6 1.7 0.9]) % stretch its width and height 
                drawnow
                
                pf = [obj.dielectric.(flnminv).xcorrDiel{i}; obj.dielectric.(flnminv).GOFgeoDiel{i}];
                if isnan(pf(1))
                    ha = subplot(1,3,3);
                    hh1 = uicontrol('Style','text','Position',[815 265 250 20], 'String', 'Dielectric Analysis');
                    set(hh1, 'FontSize',13);
                    hh2 = uicontrol('Style','text','Position',[815 240 250 20], 'String', 'Reconstructed and forward');
                    set(hh2, 'FontSize',13);
                    hh2 = uicontrol('Style','text','Position',[815 220 250 20], 'String', 'regions do not overlap');
                    set(hh2, 'FontSize',13);
                    delete(ha);
                else                    
                    pfstr = {num2str(pf(1),'%6.3f'); num2str(pf(2),'%6.3f')};                
                    ha = subplot(1,3,3);
                    cnames = {'Value'};
                    rnames = {'Xcorr', 'GOF'};
                    ht = uitable('ColumnName', cnames, 'RowName', rnames, 'Data', pfstr,'ColumnWidth', {69}); 
                    hh = uicontrol('Style','text','Position',[815 255 250 20], 'String', 'Dielectric Property Analysis');               
                    set(hh, 'FontSize',14);
                    pos = get(ha,'Position');
                    un = get(ha, 'Units');
                    delete(ha);
                    set(ht, 'units', 'normalized');
                    set(ht, 'position', pos);
                    set(ht,'ColumnName',cnames);
                    set(ht,'RowName', rnames);
                    p = get(ht,'position'); % get subplot axis position
                    set(ht,'position',p.*[1.1 4.0 0.57 0.16]) % stretch its width and height 
                    drawnow
                end
                % Save figures in format selected by user
                fileNm = ['dielectricAnalysisImages_' flnminv ,'_' num2str(i)];
                obj.saveFigure(fileNm, dataPath, figFormat);  
                
            end
        end
        %------------------------------------------------------------------
        function roi = getROIMasks(obj, i, bdRec, bdRef, ROI, tissueMapInv, tissueMapFwd, tissueType, fwdMdl, invMdl, refMask )              
            % This function is used to identify region of interest and use 
            % this ROI to construct reference and reconstructed masks relative
            % to the ROI. If dielectric property analysis is required, then
            % dielectric property maps extracted with the masks relative to
            % the ROI also constructed.
            % Inputs:
            %   i - integer - index of the ith reconstructed mask
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRef -structure - interface from the corresponding reference region.
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.
            %   tissueMapFwd - n by m tissue map of formed by segmenting the forward model.
            %   tissueMapInv - n by m tissue map formed by segmenting the inverse model. 
            %   tissueType - string - indicates what,i.e., forward model, or ReEps, ImEps, magEps inverse model.
            %   fwdMdl - n by m query component of forward model complex permittivity map
            %   invMdl - n by m query component of inverse model complex permittivity map                                
            %   refMask - n by m query reference mask.
            % Outputs:
            %   roi.recMask_i - reconstructed mask being querried
            %   roi.refMask_j - reference regions that overlap with the query reconstructed mask
            %   roi.invMdl_i - reconstructed tissue extracted when roi.recMask_i is applied to the inverse model
            %   roi.fwdMdl_j - reference tissue extracted when roi.fwdMdl_j is applied to the forward model
            %--------------------------------------------------------------
            % initialize roi structure
            roi = [];
            
            % construct reconstructed mask from interface
            recMask_i = zeros(obj.n, obj.m);
            recMask_i = roipoly(ROI,bdRec{i}.pts(:,2), bdRec{i}.pts(:,1)); 
            
            refMask_j = zeros(obj.n, obj.m);  
            % Get all reference masks within the ROI for which the ith 
            % reconstructed countour touches or is within the region bound
            % by a reference mask.
            Xrec_i = [obj.xNodes(bdRec{i}.pts(:,2))' obj.yNodes(bdRec{i}.pts(:,1))']; 
            refMask_j = obj.getRefMask(bdRef, Xrec_i, ROI);
            
            % X1 is the Intersection of the query reconstructed mask 
            % (recMask_j) and all reference masks within the interior.
            X1 = zeros(obj.n, obj.m);            
            X = recMask_i + refMask;            
            X1(find(X == 2)) = 1; 
            % Update refMask_j as the Union of X1 and reference masks that
            % are touching the reconstructed mask but are not inside            
            refMask_j = X1 + refMask_j;
            % Correct for mistakes made when mask is constructed from a
            % contour
            refMask_j(find(refMask_j == 2)) = 1;
                
            % The contours used to construct the masks are estimated by
            % sampling the edge of the tissue masks. Consequently, isolated
            % region of malignant tissue will not be extracted with this technique.
            % This means that the tissue maps constructed in the segmentation
            % step are used to identify malignant tissue embedded in glandular
            % tissue. The regional masks are corrected by using the tissue maps.
            mapInv = recMask_i.*(tissueMapInv);
            recMask_i = zeros(obj.n, obj.m);
            mapFwd = refMask_j.*(tissueMapFwd);
            refMask_j = zeros(obj.n, obj.m);
            if strcmp(tissueType, 'glandular')
                recMask_i(find(mapInv > 2 & mapInv < 5)) = 1;                    
                refMask_j(find(mapFwd > 2 & mapFwd < 5)) = 1;
            else                    
                recMask_i(find(mapInv == 5)) = 1;
                refMask_j(find(mapFwd == 5)) = 1;
            end
            
            % Use masks to identify region of interest and use this ROI to construct
            % reference and reconstructed masks relative to the ROI.
            ROI = obj.getROI(recMask_i, refMask_j);
            roi.recMask_i = recMask_i([ROI.rowStart:ROI.rowEnd],[ROI.colStart:ROI.colEnd],:);
            roi.refMask_j = refMask_j([ROI.rowStart:ROI.rowEnd],[ROI.colStart:ROI.colEnd],:);
            invMdl_i = invMdl([ROI.rowStart:ROI.rowEnd],[ROI.colStart:ROI.colEnd],:);
            fwdMdl_j = fwdMdl([ROI.rowStart:ROI.rowEnd],[ROI.colStart:ROI.colEnd],:);
                
            roi.invMdl_i = invMdl_i.*roi.recMask_i;                
            roi.fwdMdl_j = fwdMdl_j.*roi.refMask_j; 
        end                
        %------------------------------------------------------------------  
        function refMask = getRefMask(obj, bdRef, Xrec, ROI)
            % This function is used to construct reference mask from reference
            % interfaces if the reconstructed interface is on the reference
            % interface or within the area bound by the reference interface.
            % Inputs:            
            %   bdRef -structure - interface from the corresponding reference region.
            %   Xrec - array - interface of reconstructed mask
            %   ROI - n by m region of interest that contains the reconstructed and reference masks.           
            % Outputs:            %   
            %   refMask - n by m mask of reference regions that overlap with the query reconstructed mask            
            %--------------------------------------------------------------
            % Construct reference mask from reference interfaces if the
            % reconstructed interface is on the reference interface or 
            % within the area bound by the reference interface.
            refMask = zeros(obj.n, obj.m);
            for j = 1:length(bdRef)
                Xref = [obj.xNodes(bdRef{j}.pts(:,2))' obj.yNodes(bdRef{j}.pts(:,1))'];
                refMask_j = zeros(obj.n, obj.m);
                % Only include reference interface in the analysis if the
                % reconstructed interface is on the reference interface
                % or within the area bound by the reference interface.
                if obj.isRecOnOrtInRef(Xref,Xrec)
                    refMask_j = roipoly(ROI,bdRef{j}.pts(:,2), bdRef{j}.pts(:,1));
                end
                refMask = refMask + refMask_j;
            end
        end % function getRefMask
        %------------------------------------------------------------------                      
        function distances = cmpBnd(obj, bdRef, bdRec)
            % Function cmpBnd compares each point on the extracted
            % reconstructed tissue interface with the nearest point on the
            % reference tissue interface.
            % inputs:
            %   bdRec -structure - interface from the query reconstructed region.
            %   bdRef -structure - interface from the corresponding reference region.  
            % outputs:
            %   distances - array - distance between each point on the reconstructed interface with the nearest point on the reference interface.
            %--------------------------------------------------------------       
            for i = 1:length(bdRec)
                % For each reconstructed interface, get all of its points.                
                Xrec_i = [obj.xNodes(bdRec{i}.pts(:,2))' obj.yNodes(bdRec{i}.pts(:,1))'];
                % Initialize a distance vector of the distances between
                % points on the reconstructed interface and the nearest
                % point on a reference interface (where there may be 1 or
                % more reference interfaces to compare with).
                di = zeros(length(Xrec_i),1);
                % For each point (i.e., the query point) on the estimated
                % interface, find the nearest distance between that point
                % and all reference interfaces.
                for m = 1:length(Xrec_i)
                    dmin_j = zeros(length(bdRef),1);
                    for j = 1:length(bdRef)
                        Xref_j = [obj.xNodes(bdRef{j}.pts(:,2))' obj.yNodes(bdRef{j}.pts(:,1))'];
                        [k,dmin_j(j)] = dsearchn(Xref_j, Xrec_i(m,:) );                        
                    end
                    % Determine the smallest distance between the query
                    % point and the reference interfaces, and the reference
                    % contour for which the queary point is closest to.
                    [di(m), idx] = min(dmin_j);
                    % Determine if the queary point is inside or outside
                    % the reference interface that the query point is closest.
                    Xref_j = [obj.xNodes(bdRef{idx}.pts(:,2))' obj.yNodes(bdRef{idx}.pts(:,1))'];
                    if inpolygon(Xrec_i(m,1),Xrec_i(m,2),Xref_j(:,1),Xref_j(:,2) )
                        di(m) = -di(m);
                    end
                end
                distances{i}.error = di;
                distances{i}.errorMax = max(di);
                distances{i}.errorMin = min(di);
                distances{i}.errorMean = mean(di);
                distances{i}.errorStd = std(di);
            end            
        end % cmpBnd
        %------------------------------------------------------------------        
        function in_Or_on = isRecOnOrtInRef(obj,Xref,Xrec)
            % This function checks if the reconstructed interface is on the 
            % reference interface or within the region bound by the interface.
            % inputs:
            %   Xref - array - interface from the query reconstructed region.
            %   Xrec - array - interface from the corresponding reference region.  
            % outputs:
            %   in_Or_on - true if point is on or within the region, false otherwise.
            %--------------------------------------------------------------
            [in,on] = inpolygon(Xrec(:,1),Xrec(:,2),Xref(:,1),Xref(:,2));
            idx_in = find(in == 1);
            idx_on = find(on == 1);
            if ~isempty(idx_in) | ~isempty(idx_on)
                in_Or_on = true;
            else
                in_Or_on = false;
                % Check for case where reference interface is completely
                % within the region bound by the reconstructed interface.
                [in,on] = inpolygon(Xref(:,1),Xref(:,2),Xrec(:,1),Xrec(:,2));
                idx_in = find(in == 1);
                idx_on = find(on == 1);
                if ~isempty(idx_in) | ~isempty(idx_on)
                    in_Or_on = true;
                end                    
            end
        end % function isRecOnOrtInRef
        %------------------------------------------------------------------               
        function imgType = getRecType(obj, flnmInv)
            % This function displays: (1) the query reconstructed tissue, (2) relevent regions of forward
            % model that the reconstructed region is compared with, and (3)
            % a table of the metric values.                 
            % Inputs:            
            %   flnminv - string - reconstructed region, e.g., glandularRe, tumorRe, glandularIm, tumorIm, glandularMag, tumorMag.           
            % Inputs:            
            %   imgType - string - reconstructed region type - {Real, Imag, Mag}.
            %--------------------------------------------------------------       
            n = length(flnmInv);
            t = extractBetween(flnmInv, n-1, n);
            switch t{1}
                case 'Re'
                    imgType = 'Real';
                case 'Im'
                    imgType = 'Imag';
                otherwise
                    imgType = 'Mag';
            end                    
        end % function getRecType
        %------------------------------------------------------------------
        function v = translation2D(obj, v, translateVector )
            % This function translates a vector in 2D space by Tx in the x direction and Ty in the y direction which
            % are components in the translatateVector = [Tx Ty].                 
            % Inputs:            
            %   v - array - vector to be translated.
            %   translateVector - array - translate vector
            % Outputs:
            %   v - array - vector after it has been translated.
            %--------------------------------------------------------------          
            % construct translation matrix
            flag = 0;
            if size(v,1) == 1
                v = v';
                flag = 1;
            end
            T = [1 0 0; 0 1 0; translateVector(1) translateVector(2) 1];
            v3 = [v;1];
            v3 = v3'*T;
            v = v3(1:2)';
            if flag
                v = v';
            end
        end % translation2D         
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
        end % function
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
        function xy_translated = translateCoordinates(obj, xy, Tx, Ty )            
            % This algorithm translates each of the points in array by Tx and Ty.
            % Inputs:            
            %   xy - array - vector to be translated.
            %   Tx, Ty - float - x and y components of translate vector
            % Outputs:
            %   xy_translated - array - vector after it has been translated.
            %--------------------------------------------------------------
            T = [1 0 0;0 1 0;Tx Ty 1];
            [n, m] = size(xy);
            xy = [xy ones(n,1)];
            xy_new = xy*T;
            xy_translated = xy_new(:,1:2);
        end
        %------------------------------------------------------------------
        function [XrecOffset textAlgn theta] = getOffsetPts(obj, Xrec, centre, offset)            
            % This function evaluates the coordinates of a point that is offset from
            % a query point.
            % Inputs:            
            %   Xrec - array - points of the reconstructed interface.
            %   centre - coordinates of the centre of the region bound by the reconstructed interface.
            %   offset - distance from the point on the interface.
            % Outputs:
            %   XrecOffset - coordinates of point offset from the point on the reconstructed interface.
            %   textAlgn - character - alignment of the text with respect to a point based on the angular displacement of the point from the x-axis ('L', 'R', 'C').
            %   theta - float - angular displacement of point in degrees
            %--------------------------------------------------------------        
            v = Xrec;
            % figure;
            % plot(v(:,1),v(:,2));hold on;
            % plot(centre(1),centre(2),'r*','MarkerSize',5);hold on;
            translateVector = -centre;
            for i = 1:size(Xrec,1)
                v(i,:) = obj.translation2D(v(i,:), translateVector );               
            end
            geom = obj.getContourGeom( v(:,1), v(:,2) );
            centreTrn = [geom(2) geom(3)];            
            % plot(v(:,1),v(:,2));hold on;
            % plot(centreTrn(1), centreTrn(2));hold on;
            theta = (atan2(v(:,2),v(:,1)))*(180/pi);
            r = sqrt( (v(:,1)-centreTrn(1)).^2 + (v(:,2)-centreTrn(2)).^2 ) + offset;
            XrecOffset = [r.*cosd(theta) r.*sind(theta)];
            %  plot(XrecOffset(:,1),XrecOffset(:,2),'g');hold on;
            translateVector = centre;
            for i = 1:size(Xrec,1)
               XrecOffset(i,:) = obj.translation2D(XrecOffset(i,:), translateVector );               
            end
            % plot(XrecOffset(:,1),XrecOffset(:,2),'m');hold on;
            % geom = obj.getContourGeom( XrecOffset(:,1), XrecOffset(:,2) );
            % centreTrn = [geom(2) geom(3)]; 
            textAlgn = obj.getTextAlginment(theta);
        end % function getOffsetPts
        %---------------------------------------------------------------------------
        function textAlgn = getTextAlginment(obj, theta)
            % This function determines the alignment of the text with respect to a
            % point based on the angular displacement of the point from the
            % x-axis.
            % Inputs:            
            %   theta - float - angular displacement of point in degrees            
            % Outputs:
            %   textAlgn - character - {'L','C','R'} - means to the left, centre, or right of point.            
            %--------------------------------------------------------------                
            for i = 1:length(theta)
                if abs(theta(i)) >= 40 && abs(theta(i)) <=110
                    textAlgn(i) = 'C';
                else
                    if abs(theta(i)) < 40
                        textAlgn(i) = 'L';
                    else                        
                        textAlgn(i) = 'R';
                    end
                end
            end
        end % function getTextAlginment                
        %------------------------------------------------------------------
        function getCSVfile(obj, fldListfwd, fldListInv, gui)
            % This function constructs the comma separated values (CSV) file of the
            % interface, geometric, and dielectric properties of the
            % regions dominated by dense (glandular) or malignant tissue.
            % The file is written to the the results folder pointed to by
            % the gui.dataPath.
            % Inputs:            
            %   fldListfwd - field name of the forward model (eg., glandular, tumor, etc.)
            %   fldListInv - field name of the inverse model (eg., glandularRe, tumorRe, etc.)
            %   gui.dataPath - string - Points to results folder to store figures  
            %   gui.figFormat - string - Format to store figures (i.e., -fig, -png, -eps).                     
            %--------------------------------------------------------------  
            propertyList = {'xcorr Geometric', 'Goodness of Fit Geometric', 'Jaccard', 'DiceCoef', 'ratio-of-detection', 'artifact rejection ratio', 'Specificity', 'Precision', 'average Hausdorff Distance (mm)', 'xcorr Dielectric', 'Goodness of Fit Dielectric'};
            
            % Open csv file- write data to file            
            fid = fopen([gui.dataPath '\results/' gui.output_flnm],'w');
            
            % Write header information
            t = date;
            c = clock;
            tmp = ['/* CSV file constructed on ' num2str(c(4)) ':' num2str(c(5)) ' ' t '*/\n']; fprintf(fid,tmp);
            tmp = ['Performance measures for reconstruction. \n']; fprintf(fid,tmp);
            tmp = ['\n']; fprintf(fid,tmp);
            i = 1;
            for j =1:6
                if j < 4
                    for k = 1:length(obj.geometric.(fldListInv{j}).xcorrGeo)
                        tmp = [fldListfwd{i} ',' obj.getRecType(fldListInv{j}) ', number ' num2str(k) '\n']; fprintf(fid,tmp);                        
                        tmp = [propertyList{1} ','  num2str(obj.geometric.(fldListInv{j}).xcorrGeo{k},'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{2} ','  num2str(obj.geometric.(fldListInv{j}).goodnessOfFit{k},'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{3} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.Jaccard,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{4} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.DiceCoef,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{5} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.RD,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{6} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.AR,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{7} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.Specificity,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{8} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.Precision,'%15.6f') '\n']; fprintf(fid,tmp);                        
                        tmp = [propertyList{9} ','  num2str(obj.geometric.(fldListInv{j}).avgHd{k}*1e3,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{10} ',' num2str(obj.dielectric.(fldListInv{j}).xcorrDiel{k},'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{11} ',' num2str(obj.dielectric.(fldListInv{j}).GOFgeoDiel{k},'%15.6f') '\n']; fprintf(fid,tmp);
                    end
                else
                    if j == 4
                        i = i + 1;
                    end
                    for k = 1:length(obj.geometric.(fldListInv{j}).xcorrGeo)
                        tmp = [fldListfwd{i} ',' obj.getRecType(fldListInv{j}) ', number ' num2str(k) '\n']; fprintf(fid,tmp);                        
                        tmp = [propertyList{1} ','  num2str(obj.geometric.(fldListInv{j}).xcorrGeo{k},'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{2} ','  num2str(obj.geometric.(fldListInv{j}).goodnessOfFit{k},'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{3} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.Jaccard,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{4} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.DiceCoef,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{5} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.RD,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{6} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.AR,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{7} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.Specificity,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{8} ','  num2str(obj.geometric.(fldListInv{j}).similarity{k}.Precision,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{9} ','  num2str(obj.geometric.(fldListInv{j}).avgHd{k}*1e3,'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{10} ',' num2str(obj.dielectric.(fldListInv{j}).xcorrDiel{k},'%15.6f') '\n']; fprintf(fid,tmp);
                        tmp = [propertyList{11} ',' num2str(obj.dielectric.(fldListInv{j}).GOFgeoDiel{k},'%15.6f') '\n']; fprintf(fid,tmp);
                    end
                end
            end
            fclose(fid);
        end % function getCSVfile
        %------------------------------------------------------------------
        function roiBox = getROI(obj,recMask, refMask)
            % This function finds the box that bounds the region of interest
            % formed by the union of the reconstructed and reference masks.
            % Inputs:            
            %   recMask - n by m - query reconstructed mask
            %   refMask - n by m - reference mask            
            % Outputs:
            %   roiBox.colStart - (float) - start x coordinate
            %   roiBox.colEnd - (float) - end x coordinate
            %   roiBox.rowStart - (float) - start y coordinate
            %   roiBox.rowEnd - (float) - end y coordinate
            %--------------------------------------------------------------      
            
            % Function finds the box that bounds the region of interest
            % formed by the union of the reconstructed and reference masks.
            
            % Find region of interest. Namely, the union between the
            % reference and reconstructed masks.
            roi =recMask | refMask;
            % Find the box that bounds the region of interest
            bdb=regionprops(roi,'BoundingBox','Area');
            a = zeros(size(bdb,1),1);
            for i = 1:size(bdb,1)
                a(i) = bdb(i).Area;
            end
            [aMax, idx] = max(a);
            % Get the coordinates of the bounding box.
            roiBox.colStart = floor(bdb(idx).BoundingBox(1)); % start x
            roiBox.colEnd  = ceil(bdb(idx).BoundingBox(1)) + bdb(idx).BoundingBox(3); % end x
            roiBox.rowStart = floor(bdb(idx).BoundingBox(2)); % start y
            roiBox.rowEnd = ceil(bdb(idx).BoundingBox(2)) + bdb(idx).BoundingBox(4) + 1; % end y
        end % function getROI
        %------------------------------------------------------------------
    end % methods
    %======================================================================    
end % Class performanceMetrics