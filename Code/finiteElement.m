classdef finiteElement
    % Class finiteElement
    % Data structure holds the parameters related to the finite element
    % method contrast source inversion results (FEM-CSI), and the forward
    % model used by the FEM forward solver.
    % Methods are used to recover the complex permittivity profile from the
    % estimated contrast profile, transform the triangluar mesh
    % elements to square mesh elements, and segment the forward and inverse
    % models into reference and reconstructed tissue masks.
    %======================================================================    
    properties
        finc                    % Incident field frequency (Hz)
        status
        errorFlag
        fwdMdl                  % Forward model used by FEM forward solver to generate field data. 
        invMdl                  % Model reconstructed by FEM-CSI algorithm        
        contrastIteration       % Iteration of the FEM-CSI algorithm to recover contrast profile
        eprIm                   % Complex permittivity of immersion medium.
        xNodes                  % engineering represenation of the centre of each node along the x-axis of the imaging domain with rectangular pixels
        yNodes                  % engineering represenation of the centre of each node along the y-axis of the imaging domain with rectangular pixels        
        n                       % Dimension of the model/mask along the y-axis (i.e., number of rows)
        m                       % Dimension of the model/mask along the x-axis (i.e., number of columns)
        segmenationFwdValues    % Forward model segmentation thresholds
    end
    %======================================================================
    methods
        %------------------------------------------------------------------
        function obj = finiteElement( contrastIteration, inputModelFormat)
            % Constructor method to create instance of object. When instance
            % object is created, the method initializes the
            % contrastIteration and errorFlag properties.
            % Inputs:
            %   contrastIteration - (integer) iteration to construct complex permittivity profile related to inverse model.
            %   inputModelFormat - (string) input model format - 'Triangular mesh' or 'Rectangular mesh'            
            %--------------------------------------------------------------            
            if strcmp(inputModelFormat, 'Triangular mesh')
                obj.contrastIteration = contrastIteration;            
            end
            obj.errorFlag = false;
        end % Constructor finiteElement
        %------------------------------------------------------------------         
        function obj = rectangularMesh(obj, gui)
            % This function: (1) recovers the complex permittivity map at 
            % the requested iteration from the contrast profile, (2) locates
            % nodes and triangles within the reconstruction domain, (3) maps the
            % contrast profile to the nodes, and (4) transforms the triangular
            % mesh elements of the reconstruction model to square mesh elements
            % used by the imaging domain.
            % Inputs:
            %   gui.dataPath - string - Points to results folder to store figures  
            %   gui.figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------
            filePathLong = gui.dataPath;                           
            dataFilePath  = [filePathLong 'data/'];
            configurationPath = [filePathLong 'Configuration/'];
            
             % Load frequency file            
            try 
                obj.finc = (dlmread([configurationPath 'frequency.txt']))';                                    
            catch
                obj.status = 'Error: frequency file frequency.txt does not exist';
                obj.errorFlag = true;
                return
            end
            w = 2*pi*obj.finc;
            e_o = 8.85419e-12;            
            
            % Load imaging domain information related to the mesh elements
            % representing the forward and reconstructed model.
            try 
                obj.xNodes = (dlmread([configurationPath 'xNodes.txt']))';
                obj.yNodes = (dlmread([configurationPath 'yNodes.txt']))';
                obj.invMdl.xNodes = obj.xNodes;
                obj.invMdl.yNodes = obj.yNodes;
                obj.fwdMdl.xNodes = obj.xNodes;
                obj.fwdMdl.yNodes = obj.yNodes;                
            catch
                obj.status = 'Error: xNodes.txt or yNodes.txt does not exist';
                obj.errorFlag = true;
                return
            end
            obj.m = length(obj.xNodes);
            obj.n = length(obj.yNodes);
            obj.fwdMdl.m = obj.m;
            obj.fwdMdl.n = obj.n;
            obj.invMdl.m = obj.m;
            obj.invMdl.n = obj.n;
                        
            % Load the complex permittivity of the forward and
            % reconstructed models (or inversion models); Store this into the data structure of the 
            % finite element object.
            
            load([dataFilePath gui.input_forwardModel_flnm]);
           
            obj.fwdMdl.complexPermittivity = input_forwardModel;
            obj.fwdMdl.eprArray = real(obj.fwdMdl.complexPermittivity);
            obj.fwdMdl.sigmaArray = -w.*e_o.*imag(obj.fwdMdl.complexPermittivity);            
            
            load([dataFilePath gui.input_reconstructedModel_flnm]);            
            obj.invMdl.complexPermittivity = input_reconstructedModel;
            obj.invMdl.eprArray = real(obj.invMdl.complexPermittivity);
            obj.invMdl.sigmaArray = -w.*e_o.*imag(obj.invMdl.complexPermittivity);  
        end
        %------------------------------------------------------------------
        function obj = tri2square(obj, gui)
            % This function: (1) recovers the complex permittivity map at 
            % the requested iteration from the contrast profile, (2) locates
            % nodes and triangles within the reconstruction domain, (3) maps the
            % contrast profile to the nodes, and (4) transforms the triangular
            % mesh elements of the reconstruction model to square mesh elements
            % used by the imaging domain.
            % Inputs:
            %   gui.dataPath - string - Points to results folder to store figures  
            %   gui.figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %--------------------------------------------------------------
            FilePathLong = gui.dataPath;                           
            ContrastFilePath  = [FilePathLong 'data/'];
            ConfigurationPath = [FilePathLong 'Configuration/'];
            if exist([ContrastFilePath 'Contrast',num2str(obj.contrastIteration),'.txt'], 'file') ~= 2
                error(['Contrast',num2str(obj.contrastIteration),'.txt does not exist. Enter a different interation.' ]);
            end
            
            % Create imaging domain having square mesh elements
            try 
                obj.xNodes = (dlmread([ConfigurationPath 'xNodes.txt']))';
                obj.yNodes = (dlmread([ConfigurationPath 'yNodes.txt']))';
                obj.invMdl.xNodes = obj.xNodes;
                obj.invMdl.yNodes = obj.yNodes;
                obj.fwdMdl.xNodes = obj.xNodes;
                obj.fwdMdl.yNodes = obj.yNodes;                
            catch
                obj.status = 'Error: xNodes.txt or yNodes.txt does not exist';
                obj.errorFlag = true;
                return
            end
            obj.m = length(obj.xNodes);
            obj.n = length(obj.yNodes);
            [XU,YU]=meshgrid(obj.invMdl.xNodes,obj.invMdl.yNodes);

            % Read information related to reconstruction model (i.e.,
            % incident field frequency, boundary conditions, measurment
            % chamber radius, etc.) from the xml file.
           try
               [DomainData,InversionData] = obj.InputDataXML(ConfigurationPath, gui.input_csi_Img_Rec_flnm);                
           catch
               obj.status = 'Error: Check your input files.';
               obj.errorFlag = true; 
               return
           end           
           
            % Incident frequency used to generate the forward data
            obj.finc = DomainData.Freq;
            
            % Extract geometric information from GMSH file related to each triangular mesh element
            % used by the reconstruction model.            
            meshInvData = obj.getTriMeshGeometry(DomainData);
            
            % Next, get the background permittivity so that we can retrieve
            % the reconstructed complex permittivity profile from the reconstructed
            % contrast profile.
            if DomainData.inhomogEnable == 1
                % For a inhomogeneous background (i.e., background permittivity acquired from 
                % prior information), map the inhomogeneous values onto the
                % mesh.                
                epr_bkg = obj.getEprBackground(DomainData, meshInvData, 1);
            else
                % otherwise the immersion medium is used as the background
                % complex permittivity (i.e., an homogeneous background
                % complex permittivity is assumed).
                 epr_bkg = DomainData.RealBackground + 1i*DomainData.ImagBackground;
            end
                        
            % Get physics information related to inhomogeneous background
            if DomainData.inhomogEnable == 1 && DomainData.inhomogMapEnable == 1
                MappingTmp = dlmread(DomainData.inhomogMapFile);
                meshInvData.Triangles.TriPhysics = MappingTmp;
            end
            
            % Transform the triangular mesh elements of the reconstruction model
            % to square mesh elements.            
            % Read the dynamic contrast values extracted from the requested iteration (contrastIteration)            
            Chi = dlmread([ContrastFilePath 'Contrast',num2str(obj.contrastIteration),'.txt']);
            
            w = 2*pi*DomainData.Freq;

            % Find all triangular mesh elements within the imaging domain
            % for which the contrast profile is dynamic (i.e., the estimated contrast values 
            % change with each iteration). 
            [ImgDomainTriMatrix , ImgDomainTriID, ImgDomainNodes, BoundaryNodes, NonBoundaryNodes, BoundaryTriangles] = obj.findDynamicElements(DomainData, meshInvData);
            
            switch DomainData.MaterialFileType
                case {4,5,6}
                    % It is assumed that each contrast value is mapped to
                    % a finite element node. Hence, ensure that the number of nodes
                    % inside imaging domain matches the number of entries in the file
                    % storing the contrast values.            
                    if length(Chi(:,1)) ~= length(ImgDomainNodes)
                        % If not, then throw an error.
                        obj.status = 'Error: Mismatch in size of output file and number of unknowns inside imaging domain.';
                        obj.errorFlag = true;
                        return;
                    end
                    % Map contrast values to the nodes of the mesh elements contained in the
                    % imaging domain.
                    Chi_T = zeros(meshInvData.Nodes.nNodes, 1);
                    Chi_T(ImgDomainNodes) = Chi;
                    % Finally, retrieve the complex permittivity profile from the background permittivity 
                    % map and the reconstructed profile, and interpolate the triangular cell elements
                    % of the reconstruction model to square cell elements used by the imaging domain.
                    % This step is carried out similtaneously with the griddata function.
                    obj.invMdl.complexPermittivity=griddata(meshInvData.Nodes.xNodes,meshInvData.Nodes.yNodes,epr_bkg.*(Chi_T+1),XU,YU,'nearest')';
                    % Furthermore, map the contrast profile from the
                    % trianglular mesh elements to the square pixels.
                    obj.invMdl.contrast=griddata(meshInvData.Nodes.xNodes,meshInvData.Nodes.yNodes,Chi_T,XU,YU,'nearest')';
                case {8,9}
                    % It is assumed that each contrast value is mapped to
                    % the centroid of a finite element. Hence, ensure that the
                    % number of triangular element centroids inside the imaging domain
                    % (i.e., the number of dynamic finite elements) matches the number
                    % of entries in the file storing the contrast values.    
                    if length(Chi(:,1)) ~= length(ImgDomainTriID)
                        % If not, then throw an error.
                        obj.status = 'Error: Mismatch in size of output file and number of unknowns inside imaging domain.';
                        obj.errorFlag = true;
                        return;                        
                    end
                    % Map contrast values to the centroids of the mesh elements contained 
                    % in the imaging domain.
                    Chi_T = zeros(meshInvData.Triangles.nTriElements, 1);
                    Chi_T(ImgDomainTriID) = Chi;
                    % Finally, retrieve the complex permittivity profile from the background permittivity 
                    % map and the reconstructed profile, and interpolate the triangular cell elements
                    % of the reconstruction model to square cell elements used by the imaging domain.
                    % This step is carried out similtaneously with the griddata function.
                    obj.invMdl.complexPermittivity=griddata(meshInvData.Triangles.xCentroids,meshInvData.Triangles.yCentroids,epr_bkg.*(Chi_T+1),XU,YU,'nearest')';
                    % Furthermore, map the contrast profile from the
                    % trianglular mesh elements to the square pixels.
                    obj.invMdl.contrast=griddata(meshInvData.Triangles.xCentroids,meshInvData.Triangles.yCentroids,Chi_T,XU,YU,'nearest')';                    
                otherwise
                    % throw an error, as the material type is not
                    % recognized.
                    obj.status = 'Error: Unknown materials file type. Check the function''s help.';
                    obj.errorFlag = true;
                    return;                    
            end
            
            % Parse the real and imaginary part from the complex
            % permittivity map, and map the imaginary part to the
            % corresponding conductivity map.
            obj.invMdl.complexPermittivity = real(obj.invMdl.complexPermittivity) -j*imag(obj.invMdl.complexPermittivity);
            obj.invMdl.eprArray = real(obj.invMdl.complexPermittivity);
            obj.invMdl.sigmaArray = -w.*DomainData.e_o.*imag(obj.invMdl.complexPermittivity); 

            % Image the reconstructed profile that is now represented on
            % square pixels, rather than triangular mesh elements.
            try
               figure;imagesc(imrotate(real(obj.invMdl.complexPermittivity),90), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;            
            catch
               imagesc(real(obj.invMdl.complexPermittivity),'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;            
            end
            set(gca,'YDir','normal');title('Real(\epsilon)','FontSize',14);
            colormap(gca, gui.figColorMap);
            ch = colorbar;
            set(gca,'fontsize',14);
            
            try
                figure;imagesc(imrotate(-imag(obj.invMdl.complexPermittivity),90), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;
            catch
                imagesc(-imag(obj.invMdl.complexPermittivity), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;
            end
            set(gca,'YDir','normal');title('Imag(\epsilon)','FontSize',14)
            colormap(gca, gui.figColorMap);
            ch = colorbar;
            set(gca,'fontsize',14);
            
            try
               figure;imagesc(imrotate(abs(obj.invMdl.complexPermittivity),90), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;            
            catch
               imagesc(real(obj.invMdl.complexPermittivity),'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;            
            end
            set(gca,'YDir','normal');title('Mag(\epsilon)','FontSize',14);
            colormap(gca, gui.figColorMap);
            ch = colorbar;
            set(gca,'fontsize',14);
            
            try
                figure;imagesc(imrotate(real(obj.invMdl.contrast),90), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;
            catch
                imagesc(real(obj.invMdl.contrast), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;
            end                
            set(gca,'YDir','normal');title('Real(\chi)','FontSize',14)
            colormap(gca, gui.figColorMap);
            ch = colorbar;
            set(gca,'fontsize',14);
            
            try
                figure;imagesc(imrotate(-imag(obj.invMdl.contrast),90), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;
            catch
                imagesc(-imag(obj.invMdl.contrast), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;
            end
            set(gca,'YDir','normal');title('Imag(\chi)','FontSize',14)
            colormap(gca, gui.figColorMap);
            ch = colorbar;
            set(gca,'fontsize',14);
            
            % Next, transform the triangular mesh elements of the FEM forward model
            % to square mesh elements.
            % Start by extracting geometric information from GMSH file related to each triangular mesh element
            % used by the FEM forward model. 
           try
               [DomainData,InversionData] = obj.InputDataXML(ConfigurationPath, gui.input_datacollect_flnm);                
           catch               
               obj.status = 'Error: Check your input files.';
               obj.errorFlag = true;
               return;      
           end  
            
            meshFwdData= obj.getTriMeshGeometry(DomainData);
            FwdMat = dlmread(DomainData.MaterialFile); 
            FwdMat=FwdMat(:,2)+1i*FwdMat(:,3);
            % Make sure data array length is the same size as the number of triangle elements.
            if length(FwdMat(:,1)) ~= meshFwdData.Triangles.nTriElements
                % If not, then throw an error.                
                obj.status = 'Mismatch in length of forward model material array and number of elements used to represent the forward model.';
                obj.errorFlag = true;
                return;     
            end

            % interpolate the triangular cell elements of the FEM forward model 
            % to square cell elements. 
            obj.fwdMdl.complexPermittivity=griddata(meshFwdData.Triangles.xCentroids,meshFwdData.Triangles.yCentroids,FwdMat,XU,YU,'nearest')';
            obj.fwdMdl.complexPermittivity = real(obj.fwdMdl.complexPermittivity)-j*imag(obj.fwdMdl.complexPermittivity);
            obj.fwdMdl.eprArray = real(obj.fwdMdl.complexPermittivity);
            obj.fwdMdl.sigmaArray = -w.*DomainData.e_o.*imag(obj.fwdMdl.complexPermittivity); 
            
            try
                figure;imagesc(imrotate(real(obj.fwdMdl.complexPermittivity),90), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;            
            catch
                imagesc(real(obj.fwdMdl.complexPermittivity), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;            
            end    
            set(gca,'YDir','normal');title('Real(\epsilon)','FontSize',14)
            colormap(gca, gui.figColorMap);
            ch = colorbar;
            set(gca,'fontsize',14);
                        
            try
                figure;imagesc(imrotate(-imag(obj.fwdMdl.complexPermittivity),90), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;
            catch
                imagesc(-imag(obj.fwdMdl.complexPermittivity), 'XData', obj.invMdl.xNodes, 'YData', -obj.invMdl.yNodes);hold on;
            end
            set(gca,'YDir','normal');title('Imag(\epsilon)','FontSize',14);
            colormap(gca, gui.figColorMap);
            ch = colorbar;
            set(gca,'fontsize',14);
        end % function tri2square                
        %------------------------------------------------------------------
        function obj = iniSegmentationFwdThresholds(obj, gui)
            % This function sets the threshold values calculated by the GUI
            % when the input task is run.
            % Inputs:
            %   gui.segmenationIniValues - structure - stores various threshold values calculated with the GUI object.
            % Outputs:
            %   gui.segmenationFwdValues - structure - stores various threshold values calculated with the GUI object used to segment the forward model.            
            %--------------------------------------------------------------
            obj.segmenationFwdValues.fwd_skin_epr_Hi_Limit = gui.segmenationIniValues.fwd_skin_epr_Hi_Limit;
            obj.segmenationFwdValues.fwd_fat_epr_Hi_Limit = gui.segmenationIniValues.fwd_fat_epr_Hi_Limit;
            obj.segmenationFwdValues.fwd_gland_epr_Lo_Limit = gui.segmenationIniValues.fwd_gland_epr_Lo_Limit;
            obj.segmenationFwdValues.fwd_gland_contours_threshold = gui.segmenationIniValues.fwd_gland_contours_threshold;
            obj.segmenationFwdValues.fwd_gland_contours_maxNm = gui.segmenationIniValues.fwd_gland_contours_maxNm;            
            obj.segmenationFwdValues.fwd_tumor_epr_Lo_Limit = gui.segmenationIniValues.fwd_tumor_epr_Lo_Limit;                        
            obj.segmenationFwdValues.fwd_tumor_contours_threshold = gui.segmenationIniValues.fwd_tumor_contours_threshold;
            obj.segmenationFwdValues.fwd_tumor_contours_maxNm = gui.segmenationIniValues.fwd_tumor_contours_maxNm;
        end % function iniSegmentationFwdThresholds
        %------------------------------------------------------------------
        function MeshData = getTriMeshGeometry(obj, DomainData)
            % Function extracts geometric information from GMSH file related to each triangular
            % mesh element used in the reconstructed model. 
            %
            % Read the number of nodes in the mesh
            % Inputs:
            %   DomainData - structure - contains information about the FEM solver parameters
            % Outputs:
            %   MeshData - structure holds the mesh element information,
            %           and has the  following fields:
            %           - .Nodes     : nodes id and coordinates.
            %           - .Points     : points in the mesh with physics.
            %           - .Lines      : lines in the mesh with physics.
            %           - .Triangles : triangles in the mesh - Triangles.nTriElements, Triangles.NodesMatrix, 
            %                              Triangles.TriNodes, Triangles.TriPhysics, Triangles.TriID.
            %--------------------------------------------------------------
            file_name = DomainData.MeshLocation;                                
            Nodes.nNodes = dlmread(file_name, ' ',[4 0 4 0]');

        % Read the NodesID and the coordinates of the nodes
            tmpData = dlmread(file_name, ' ',[5 0 (Nodes.nNodes+4) 3]);  
            Nodes.NodesID = tmpData(:,1);
            Nodes.xNodes  = tmpData(:,2);
            Nodes.yNodes  = tmpData(:,3);

        % Read the number of elements in the mesh
            nElements = dlmread(file_name, ' ',[(Nodes.nNodes+7) 0 (Nodes.nNodes+7) 0]);

        % Read all the elements (nodes + line + triangular)
            try
                tmpData = dlmread(file_name,' ',[(Nodes.nNodes+8) 0 (nElements + Nodes.nNodes+7) 8]);
                tmpData(:,6) = [];
                display('Gmsh file version 2.4.2 or below detected...');
            catch %#ok<CTCH>
                tmpData = dlmread(file_name,' ',[(Nodes.nNodes+8) 0 (nElements + Nodes.nNodes+7) 7]);
                display('Gmsh file version 2.5 or higher detected...');
            end

        % Seperating line elements from triangular elements
        % Point elements
            itemp = tmpData(tmpData(:,2) == 15,7);
            PointElements = tmpData(itemp, :);
            Points.PointsID = itemp;
            Points.nPoints = length(PointElements);
            Points.PointPhysics = PointElements(:,4);

        % Line elements
            itemp = tmpData(tmpData(:,2) == 1);
            LineElements = tmpData(itemp, :);
            Lines.nLineElements = length(itemp);
            Lines.LineMatrix = [LineElements(:,6) LineElements(:,7)];
            Lines.LineNodes = unique([LineElements(:,6);LineElements(:,7)]);
            Lines.LinePhysics = LineElements(:,4);

        % Boundary Nodes
            Nodes.BdNodes = Lines.LineNodes;

        % Non-boundary nodes (aka free nodes);
            tmp = zeros(Nodes.nNodes,1);
            tmp(Nodes.BdNodes) = 1;
            Nodes.NonBdNodes = find(tmp == 0);

        % Calculating the midpoint of each line segment
            Lines.xMidpoint = (Nodes.xNodes(Lines.LineMatrix(:,1)) + Nodes.xNodes(Lines.LineMatrix(:,2)))/2;
            Lines.yMidpoint = (Nodes.yNodes(Lines.LineMatrix(:,1)) + Nodes.yNodes(Lines.LineMatrix(:,2)))/2;

        % Calculating the unit normal vector for each line element

        % (1) Calculating the slope of each line segment and the slope of its normal
            SegSlope = (Nodes.yNodes(Lines.LineMatrix(:,1)) - Nodes.yNodes(Lines.LineMatrix(:,2)))./(Nodes.xNodes(Lines.LineMatrix(:,1)) - Nodes.xNodes(Lines.LineMatrix(:,2)));
            SegNormalSlope = -1./SegSlope;

        % (2) Calculating the length of each line segment element
            Lines.Length    = sqrt((Nodes.xNodes(Lines.LineMatrix(:,1))-Nodes.xNodes(Lines.LineMatrix(:,2))).^2 ...
                + (Nodes.yNodes(Lines.LineMatrix(:,1))-Nodes.yNodes(Lines.LineMatrix(:,2))).^2);

        % (3) Calculating the normal unit vector of each segment
            xNormal = 1./sqrt(1 + SegNormalSlope.^2);
            yNormal = SegNormalSlope ./sqrt(1+SegNormalSlope.^2);

        % (4) Verifying that the normal unit vector is pointing outwards with respect to the problem domain
        % If it is inwards, it is forced outwards.
            itmpX = find(Lines.xMidpoint < 0);
            xNormal(itmpX) = -xNormal(itmpX);
            yNormal(itmpX) = -yNormal(itmpX);

            Lines.UnitNormal = [xNormal yNormal];

            xNaN = isnan(Lines.UnitNormal(:,1));
            yNaN = isnan(Lines.UnitNormal(:,2));

            xNaN = find(xNaN == 1);
            yNaN = find(yNaN == 1);

            % Assumption the origin of the problem is (0,0)
            if isempty(xNaN) == 0
                Lines.UnitNormal(xNaN,1) = Lines.xMidpoint(xNaN)./abs(Lines.xMidpoint(xNaN));
            end

            if isempty(yNaN) == 0
                Lines.UnitNormal(yNaN,2) = Lines.yMidpoint(yNaN)./abs(Lines.yMidpoint(yNaN));
            end

        % Triangular elements
            itemp = tmpData(tmpData(:,2) == 2);
            TriElements  = tmpData(itemp, :);
            Triangles.nTriElements = length(itemp);
            Triangles.NodesMatrix  = [TriElements(:,6) TriElements(:,7) TriElements(:,8)];
            Triangles.TriNodes     = unique([TriElements(:,6);TriElements(:,7);TriElements(:,8)]);
            Triangles.TriPhysics   = TriElements(:,4);
            Triangles.TriID        = 1:Triangles.nTriElements;

        % Calculate the centroid of each triangle
            x1 = Nodes.xNodes(Triangles.NodesMatrix(:,1)); x2 = Nodes.xNodes(Triangles.NodesMatrix(:,2)); x3 = Nodes.xNodes(Triangles.NodesMatrix(:,3));
            y1 = Nodes.yNodes(Triangles.NodesMatrix(:,1)); y2 = Nodes.yNodes(Triangles.NodesMatrix(:,2)); y3 = Nodes.yNodes(Triangles.NodesMatrix(:,3));

            Triangles.xCentroids = (x1 + x2 + x3)/3;
            Triangles.yCentroids = (y1 + y2 + y3)/3;
            Triangles.rCentroids = sqrt(Triangles.xCentroids.^2 + Triangles.yCentroids.^2);
            Triangles.thetaCentroids = angle(Triangles.xCentroids + 1i.*Triangles.yCentroids);

        % Calculating the area of each triangle
            Triangles.Area = abs(((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1))./ 2.0);

        % Verification that the local number scheme is counter-clockwise
            TriangleVector1 = [(x2 - x1)';(y2 - y1)'];
            TriangleVector2 = [(x3 - x1)';(y3 - y1)'];
            TriangleNormals = (TriangleVector1(1,:).*TriangleVector2(2,:) - TriangleVector2(1,:).*TriangleVector1(2,:))';

            if any(TriangleNormals < 0)
                error('For some triangular elements, the nodes are numbered in clockwise manner! Fix Mesh. Program Terminated.');
            end

        % Finding the unique edges shared by the triangles
            AllTriEdges = [Triangles.NodesMatrix(:,[1,2]); Triangles.NodesMatrix(:,[2,3]); Triangles.NodesMatrix(:,[3,1])];        % All the edges of the triangles - Not Unique
            AllTriEdgesID = (1:size(AllTriEdges))'; % Assigning ids to all the edges
            [tmp1, tmp2] = sort(AllTriEdges,2);     % Sorting the edges
            tmp2 = tmp2(:,2) - tmp2(:,1);

            Triangles.EdgesDirection = reshape(tmp2, Triangles.nTriElements, 3); % The direction of the edge in each triangle

            [TriEdges,tmp,jEdges] = unique(tmp1,'rows'); % Unique edges in the mesh
            TriEdgesID     = AllTriEdgesID(jEdges);    % Unique edge ids
            Triangles.EdgesLength   = sqrt((Nodes.xNodes(TriEdges(:,1)) - Nodes.xNodes(TriEdges(:,2))).^2 + ...
            (Nodes.yNodes(TriEdges(:,1)) - Nodes.yNodes(TriEdges(:,2))).^2); % Length of each edge

        % Global Tangential vectors for each edge
            Triangles.Global_T_x = (Nodes.xNodes(TriEdges(:,2)) - Nodes.xNodes(TriEdges(:,1)))./Triangles.EdgesLength;
            Triangles.Global_T_y = (Nodes.yNodes(TriEdges(:,2)) - Nodes.yNodes(TriEdges(:,1)))./Triangles.EdgesLength;

            Triangles.EdgesMidpoint = [(Nodes.xNodes(TriEdges(:,1)) + Nodes.xNodes(TriEdges(:,2)))./2 ...
                (Nodes.yNodes(TriEdges(:,1)) + Nodes.yNodes(TriEdges(:,2)))./2]; % Midpoint of each edge
        % Unique edges in each triangle
            Triangles.EdgesMatrix    = [TriEdgesID(Triangles.TriID), TriEdgesID(Triangles.TriID + Triangles.nTriElements), TriEdgesID(Triangles.TriID+2*Triangles.nTriElements)];

        % Triangles.TriEdgesID = TriEdgesID;
            Triangles.TriEdges = TriEdges;
            Triangles.nEdges   = length(Triangles.EdgesLength);

        % Boundary Edges

            tmp1 = sort(Lines.LineMatrix,2);
            [tmp1, tmp2, Triangles.BdEdges] = intersect(tmp1, Triangles.TriEdges, 'rows');
            Triangles.BdEdgesPhysics = Lines.LinePhysics(tmp2);

        % Non-boundary Edges (aka free edges);
            tmp = zeros(Triangles.nEdges,1);
            tmp(Triangles.BdEdges) = 1;
            Triangles.NonBdEdges = find(tmp == 0);

        % Finding the unit vectors associated with each edge for a given triangle
            x1AllTriEdges = [Nodes.xNodes(Triangles.NodesMatrix(:,1)) Nodes.xNodes(Triangles.NodesMatrix(:,2)) Nodes.xNodes(Triangles.NodesMatrix(:,3))];
            x2AllTriEdges = [Nodes.xNodes(Triangles.NodesMatrix(:,2)) Nodes.xNodes(Triangles.NodesMatrix(:,3)) Nodes.xNodes(Triangles.NodesMatrix(:,1))];

            y1AllTriEdges = [Nodes.yNodes(Triangles.NodesMatrix(:,1)) Nodes.yNodes(Triangles.NodesMatrix(:,2)) Nodes.yNodes(Triangles.NodesMatrix(:,3))];
            y2AllTriEdges = [Nodes.yNodes(Triangles.NodesMatrix(:,2)) Nodes.yNodes(Triangles.NodesMatrix(:,3)) Nodes.yNodes(Triangles.NodesMatrix(:,1))];

            dxAllTriEdges = x2AllTriEdges - x1AllTriEdges;
            dyAllTriEdges = y2AllTriEdges - y1AllTriEdges;

            Triangles.AllEdgesLength = sqrt(dxAllTriEdges.^2 + dyAllTriEdges.^2);
            Triangles.x_vector_edge = dxAllTriEdges./Triangles.AllEdgesLength;
            Triangles.y_vector_edge = dyAllTriEdges./Triangles.AllEdgesLength;

        % Edge to triangle connectivity matrix

            Triangles.Edges2Triangle     = zeros(Triangles.nEdges,2);
            Triangles.Edges2TriangleLoc  = zeros(Triangles.nEdges,2);
            eCounter  = ones(Triangles.nEdges,1);
            for iTri = 1:Triangles.nTriElements
                e1 = Triangles.EdgesMatrix(iTri,1);
                e2 = Triangles.EdgesMatrix(iTri,2);
                e3 = Triangles.EdgesMatrix(iTri,3);
                % Edge 1
                Triangles.Edges2Triangle(e1,eCounter(e1)) = iTri;
                Triangles.Edges2TriangleLoc(e1,eCounter(e1)) = 1;
                eCounter(e1) = eCounter(e1)+1;
                % Edge 2
                Triangles.Edges2Triangle(e2,eCounter(e2)) = iTri;
                Triangles.Edges2TriangleLoc(e2,eCounter(e2)) = 2;
                eCounter(e2) = eCounter(e2)+1;
                % Edge 3
                Triangles.Edges2Triangle(e3,eCounter(e3)) = iTri;
                Triangles.Edges2TriangleLoc(e3,eCounter(e3)) = 3;
                eCounter(e3) = eCounter(e3)+1;
            end
            Triangles.eCounter = eCounter - 1;

        % Calculating the normal of the problem edges
        % Calculating the slope of each segment and its normal
            BdEdgex1 = Nodes.xNodes(TriEdges(Triangles.BdEdges,1));
            BdEdgex2 = Nodes.xNodes(TriEdges(Triangles.BdEdges,2));
            BdEdgey1 = Nodes.yNodes(TriEdges(Triangles.BdEdges,1));
            BdEdgey2 = Nodes.yNodes(TriEdges(Triangles.BdEdges,2));

            xSegMidpoint = Triangles.EdgesMidpoint(Triangles.BdEdges,1);

            SegSlope = (BdEdgey1 - BdEdgey2)./(BdEdgex1 - BdEdgex2);
            SegNormalSlope = -1./SegSlope;

        % Calculating the normal for each element
            xNormal = 1./sqrt(1 + SegNormalSlope.^2);
            yNormal = SegNormalSlope ./sqrt(1+SegNormalSlope.^2);

            itmpX = find(xSegMidpoint < 0);
            xNormal(itmpX) = -xNormal(itmpX);
            yNormal(itmpX) = -yNormal(itmpX);

            Triangles.BdEgNormal = [xNormal yNormal];

        % Nodes2Triangle
            T1 = Triangles.NodesMatrix(:,1);
            T2 = Triangles.NodesMatrix(:,2);
            T3 = Triangles.NodesMatrix(:,3);

            Nodes2Tri = zeros(Nodes.nNodes, 10);
            Nodes2TriCounter = zeros(Nodes.nNodes,1);
            for iTri = 1:Triangles.nTriElements    
                N1 = T1(iTri,1); N2 = T2(iTri,1); N3 = T3(iTri,1);       
                Nodes2TriCounter(N1) = Nodes2TriCounter(N1) + 1;
                Nodes2TriCounter(N2) = Nodes2TriCounter(N2) + 1;
                Nodes2TriCounter(N3) = Nodes2TriCounter(N3) + 1;    
                Nodes2Tri(N1, Nodes2TriCounter(N1)) = iTri;
                Nodes2Tri(N2, Nodes2TriCounter(N2)) = iTri;
                Nodes2Tri(N3, Nodes2TriCounter(N3)) = iTri;
            end

            maxSharing = max(Nodes2TriCounter);
            Triangles.Nodes2Tri = Nodes2Tri(:,1:maxSharing);
            Triangles.Nodes2TriCounter = Nodes2TriCounter;

            % MeshData
            MeshData.Nodes = Nodes;
            MeshData.Points = Points;
            MeshData.Lines = Lines;
            MeshData.Triangles = Triangles;
        end % function getTriMeshGeometry
        %------------------------------------------------------------------        
        function MeshData = GmshReadFwd(obj, meshFilePath)                    
            % This function extracts information related to each triangular mesh element of the 
            % FEM forward model. The mesh element information extracted is saved in 
            % structure MeshData.                         
            % Inputs:
            %   meshFilePath - string - points to folder holding file containing mesh information.
            % Outputs:
            %   MeshData - structure holds the mesh element information,
            %           and has the  following fields:
            %           - .Nodes     : nodes id and coordinates.
            %           - .Points     : points in the mesh with physics.
            %           - .Lines      : lines in the mesh with physics.
            %           - .Triangles : triangles in the mesh - Triangles.nTriElements, Triangles.NodesMatrix, 
            %                              Triangles.TriNodes, Triangles.TriPhysics, Triangles.TriID.
            %--------------------------------------------------------------
            % Read the number of nodes in the mesh
            file_name = meshFilePath;        
            Nodes.nNodes = dlmread(file_name, ' ',[4 0 4 0]');

            % Read the NodesID and the coordinates of the nodes

            tmpData = dlmread(file_name, ' ',[5 0 (Nodes.nNodes+4) 3]);
  
            Nodes.NodesID = tmpData(:,1);
            Nodes.xNodes  = tmpData(:,2);
            Nodes.yNodes  = tmpData(:,3);

            % Read the number of elements in the mesh

            nElements = dlmread(file_name, ' ',[(Nodes.nNodes+7) 0 (Nodes.nNodes+7) 0]);

            % Read all the elements (nodes + line + triangular)

            try
                tmpData = dlmread(file_name,' ',[(Nodes.nNodes+8) 0 (nElements + Nodes.nNodes+7) 8]);
                tmpData(:,6) = [];
                display('Gmsh file version 2.4.2 or below detected...');
            catch %#ok<CTCH>
                tmpData = dlmread(file_name,' ',[(Nodes.nNodes+8) 0 (nElements + Nodes.nNodes+7) 7]);
                display('Gmsh file version 2.5 or higher detected...');
            end

            % Seperate line elements from triangular elements

            % Point elements
            itemp = tmpData(tmpData(:,2) == 15,7);
            PointElements = tmpData(itemp, :);
            Points.PointsID = itemp;
            Points.nPoints = length(PointElements);
            Points.PointPhysics = PointElements(:,4);

            % Line elements
            itemp = tmpData(tmpData(:,2) == 1);
            LineElements = tmpData(itemp, :);
            Lines.nLineElements = length(itemp);
            Lines.LineMatrix = [LineElements(:,6) LineElements(:,7)];
            Lines.LineNodes = unique([LineElements(:,6);LineElements(:,7)]);
            Lines.LinePhysics = LineElements(:,4);

            % Boundary Nodes
            Nodes.BdNodes = Lines.LineNodes;

            % Non-boundary nodes (aka free nodes);
            tmp = zeros(Nodes.nNodes,1);
            tmp(Nodes.BdNodes) = 1;
            Nodes.NonBdNodes = find(tmp == 0);

        % Calculating the midpoint of each line segment
            Lines.xMidpoint = (Nodes.xNodes(Lines.LineMatrix(:,1)) + Nodes.xNodes(Lines.LineMatrix(:,2)))/2;
            Lines.yMidpoint = (Nodes.yNodes(Lines.LineMatrix(:,1)) + Nodes.yNodes(Lines.LineMatrix(:,2)))/2;

        % Calculating the unit normal vector for each line element

        % (1) Calculating the slope of each line segment and the slope of its normal
            SegSlope = (Nodes.yNodes(Lines.LineMatrix(:,1)) - Nodes.yNodes(Lines.LineMatrix(:,2)))./(Nodes.xNodes(Lines.LineMatrix(:,1)) - Nodes.xNodes(Lines.LineMatrix(:,2)));
            SegNormalSlope = -1./SegSlope;

        % (2) Calculating the length of each line segment element
            Lines.Length    = sqrt((Nodes.xNodes(Lines.LineMatrix(:,1))-Nodes.xNodes(Lines.LineMatrix(:,2))).^2 ...
                + (Nodes.yNodes(Lines.LineMatrix(:,1))-Nodes.yNodes(Lines.LineMatrix(:,2))).^2);

        % (3) Calculating the normal unit vector of each segment
            xNormal = 1./sqrt(1 + SegNormalSlope.^2);
            yNormal = SegNormalSlope ./sqrt(1+SegNormalSlope.^2);

        % (4) Verifying that the normal unit vector is pointing outwards with respect to the problem domain
        % If it is inwards, it is forced outwards.
            itmpX = find(Lines.xMidpoint < 0);
            xNormal(itmpX) = -xNormal(itmpX);
            yNormal(itmpX) = -yNormal(itmpX);

            Lines.UnitNormal = [xNormal yNormal];

            xNaN = isnan(Lines.UnitNormal(:,1));
            yNaN = isnan(Lines.UnitNormal(:,2));

            xNaN = find(xNaN == 1);
            yNaN = find(yNaN == 1);

        % Assume the origin of the problem is (0,0)
            if isempty(xNaN) == 0
                Lines.UnitNormal(xNaN,1) = Lines.xMidpoint(xNaN)./abs(Lines.xMidpoint(xNaN));
            end

            if isempty(yNaN) == 0
                Lines.UnitNormal(yNaN,2) = Lines.yMidpoint(yNaN)./abs(Lines.yMidpoint(yNaN));
            end
        % Triangular elements
            itemp = tmpData(tmpData(:,2) == 2);
            TriElements  = tmpData(itemp, :);
            Triangles.nTriElements = length(itemp);
            Triangles.NodesMatrix = [TriElements(:,6) TriElements(:,7) TriElements(:,8)];
            Triangles.TriNodes = unique([TriElements(:,6);TriElements(:,7);TriElements(:,8)]);
            Triangles.TriPhysics = TriElements(:,4);
            Triangles.TriID = 1:Triangles.nTriElements;

        % Calculate the centroid of each triangle
            x1 = Nodes.xNodes(Triangles.NodesMatrix(:,1)); x2 = Nodes.xNodes(Triangles.NodesMatrix(:,2)); x3 = Nodes.xNodes(Triangles.NodesMatrix(:,3));
            y1 = Nodes.yNodes(Triangles.NodesMatrix(:,1)); y2 = Nodes.yNodes(Triangles.NodesMatrix(:,2)); y3 = Nodes.yNodes(Triangles.NodesMatrix(:,3));

            Triangles.xCentroids = (x1 + x2 + x3)/3;
            Triangles.yCentroids = (y1 + y2 + y3)/3;
            Triangles.rCentroids = sqrt(Triangles.xCentroids.^2 + Triangles.yCentroids.^2);
            Triangles.thetaCentroids = angle(Triangles.xCentroids + 1i.*Triangles.yCentroids);

        % Calculating the area of each triangle
            Triangles.Area = abs(((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1))./ 2.0);

        % Verification that the local number scheme is counter-clockwise
            TriangleVector1 = [(x2 - x1)';(y2 - y1)'];
            TriangleVector2 = [(x3 - x1)';(y3 - y1)'];
            TriangleNormals = (TriangleVector1(1,:).*TriangleVector2(2,:) - TriangleVector2(1,:).*TriangleVector1(2,:))';

            if any(TriangleNormals < 0)
                error('For some triangular elements, the nodes are numbered in clockwise manner! Fix Mesh. Program Terminated.');
            end

         % Finding the unique edges shared by the triangles
            AllTriEdges = [Triangles.NodesMatrix(:,[1,2]); Triangles.NodesMatrix(:,[2,3]); Triangles.NodesMatrix(:,[3,1])];        % All the edges of the triangles - Not Unique
            AllTriEdgesID = (1:size(AllTriEdges))'; % Assigning ids to all the edges
            [tmp1, tmp2] = sort(AllTriEdges,2);     % Sorting the edges
            tmp2 = tmp2(:,2) - tmp2(:,1);

            Triangles.EdgesDirection = reshape(tmp2, Triangles.nTriElements, 3); % The direction of the edge in each triangle

            [TriEdges,tmp,jEdges] = unique(tmp1,'rows'); % Unique edges in the mesh
            TriEdgesID     = AllTriEdgesID(jEdges);    % Unique edge ids
            Triangles.EdgesLength   = sqrt((Nodes.xNodes(TriEdges(:,1)) - Nodes.xNodes(TriEdges(:,2))).^2 + ...
                 (Nodes.yNodes(TriEdges(:,1)) - Nodes.yNodes(TriEdges(:,2))).^2); % Length of each edge

        % Global Tangential vectors for each edge
            Triangles.Global_T_x = (Nodes.xNodes(TriEdges(:,2)) - Nodes.xNodes(TriEdges(:,1)))./Triangles.EdgesLength;
            Triangles.Global_T_y = (Nodes.yNodes(TriEdges(:,2)) - Nodes.yNodes(TriEdges(:,1)))./Triangles.EdgesLength;

            Triangles.EdgesMidpoint = [(Nodes.xNodes(TriEdges(:,1)) + Nodes.xNodes(TriEdges(:,2)))./2 ...
                (Nodes.yNodes(TriEdges(:,1)) + Nodes.yNodes(TriEdges(:,2)))./2]; % Midpoint of each edge
        % Unique edges in each triangle
            Triangles.EdgesMatrix    = [TriEdgesID(Triangles.TriID), TriEdgesID(Triangles.TriID + Triangles.nTriElements), TriEdgesID(Triangles.TriID+2*Triangles.nTriElements)];

        % Triangles.TriEdgesID = TriEdgesID;
            Triangles.TriEdges = TriEdges;
            Triangles.nEdges   = length(Triangles.EdgesLength);

        % Boundary Edges
            tmp1 = sort(Lines.LineMatrix,2);
            [tmp1, tmp2, Triangles.BdEdges] = intersect(tmp1, Triangles.TriEdges, 'rows');
            Triangles.BdEdgesPhysics = Lines.LinePhysics(tmp2);

        % Non-boundary Edges (aka free edges);
            tmp = zeros(Triangles.nEdges,1);
            tmp(Triangles.BdEdges) = 1;
            Triangles.NonBdEdges = find(tmp == 0);

        % Finding the unit vectors associated with each edge for a given triangle
            x1AllTriEdges = [Nodes.xNodes(Triangles.NodesMatrix(:,1)) Nodes.xNodes(Triangles.NodesMatrix(:,2)) Nodes.xNodes(Triangles.NodesMatrix(:,3))];
            x2AllTriEdges = [Nodes.xNodes(Triangles.NodesMatrix(:,2)) Nodes.xNodes(Triangles.NodesMatrix(:,3)) Nodes.xNodes(Triangles.NodesMatrix(:,1))];

            y1AllTriEdges = [Nodes.yNodes(Triangles.NodesMatrix(:,1)) Nodes.yNodes(Triangles.NodesMatrix(:,2)) Nodes.yNodes(Triangles.NodesMatrix(:,3))];
            y2AllTriEdges = [Nodes.yNodes(Triangles.NodesMatrix(:,2)) Nodes.yNodes(Triangles.NodesMatrix(:,3)) Nodes.yNodes(Triangles.NodesMatrix(:,1))];

            dxAllTriEdges = x2AllTriEdges - x1AllTriEdges;
            dyAllTriEdges = y2AllTriEdges - y1AllTriEdges;

            Triangles.AllEdgesLength = sqrt(dxAllTriEdges.^2 + dyAllTriEdges.^2);
            Triangles.x_vector_edge = dxAllTriEdges./Triangles.AllEdgesLength;
            Triangles.y_vector_edge = dyAllTriEdges./Triangles.AllEdgesLength;

        % Edge to triangle connectivity matrix
            Triangles.Edges2Triangle     = zeros(Triangles.nEdges,2);
            Triangles.Edges2TriangleLoc  = zeros(Triangles.nEdges,2);

         % edge counter   
            eCounter  = ones(Triangles.nEdges,1);

            for iTri = 1:Triangles.nTriElements
                e1 = Triangles.EdgesMatrix(iTri,1);
                e2 = Triangles.EdgesMatrix(iTri,2);
                e3 = Triangles.EdgesMatrix(iTri,3);
        % Edge 1
                Triangles.Edges2Triangle(e1,eCounter(e1)) = iTri;
                Triangles.Edges2TriangleLoc(e1,eCounter(e1)) = 1;
                eCounter(e1) = eCounter(e1)+1;
        % Edge 2
                Triangles.Edges2Triangle(e2,eCounter(e2)) = iTri;
                Triangles.Edges2TriangleLoc(e2,eCounter(e2)) = 2;
                eCounter(e2) = eCounter(e2)+1;
        % Edge 3
                Triangles.Edges2Triangle(e3,eCounter(e3)) = iTri;
                Triangles.Edges2TriangleLoc(e3,eCounter(e3)) = 3;
                eCounter(e3) = eCounter(e3)+1;
            end
            Triangles.eCounter = eCounter - 1;

        % Calculating the normal of the problem edges
        % Calculating the slope of each segment and its normal
            BdEdgex1 = Nodes.xNodes(TriEdges(Triangles.BdEdges,1));
            BdEdgex2 = Nodes.xNodes(TriEdges(Triangles.BdEdges,2));
            BdEdgey1 = Nodes.yNodes(TriEdges(Triangles.BdEdges,1));
            BdEdgey2 = Nodes.yNodes(TriEdges(Triangles.BdEdges,2));

            xSegMidpoint = Triangles.EdgesMidpoint(Triangles.BdEdges,1);

            SegSlope = (BdEdgey1 - BdEdgey2)./(BdEdgex1 - BdEdgex2);
            SegNormalSlope = -1./SegSlope;

        % Calculating the normal for each element
            xNormal = 1./sqrt(1 + SegNormalSlope.^2);
            yNormal = SegNormalSlope ./sqrt(1+SegNormalSlope.^2);

            itmpX = find(xSegMidpoint < 0);
            xNormal(itmpX) = -xNormal(itmpX);
            yNormal(itmpX) = -yNormal(itmpX);

            Triangles.BdEgNormal = [xNormal yNormal];
        % Nodes2Triangle
            T1 = Triangles.NodesMatrix(:,1);
            T2 = Triangles.NodesMatrix(:,2);
            T3 = Triangles.NodesMatrix(:,3);

            Nodes2Tri = zeros(Nodes.nNodes, 10);
            Nodes2TriCounter = zeros(Nodes.nNodes,1);
            for iTri = 1:Triangles.nTriElements    
                N1 = T1(iTri,1); N2 = T2(iTri,1); N3 = T3(iTri,1);
                Nodes2TriCounter(N1) = Nodes2TriCounter(N1) + 1;
                Nodes2TriCounter(N2) = Nodes2TriCounter(N2) + 1;
                Nodes2TriCounter(N3) = Nodes2TriCounter(N3) + 1;
                Nodes2Tri(N1, Nodes2TriCounter(N1)) = iTri;
                Nodes2Tri(N2, Nodes2TriCounter(N2)) = iTri;
                Nodes2Tri(N3, Nodes2TriCounter(N3)) = iTri;
            end
            maxSharing = max(Nodes2TriCounter);
            Triangles.Nodes2Tri = Nodes2Tri(:,1:maxSharing);
            Triangles.Nodes2TriCounter = Nodes2TriCounter;

        % save function results in structure MeshData
            MeshData.Nodes = Nodes;
            MeshData.Points = Points;
            MeshData.Lines = Lines;
            MeshData.Triangles = Triangles;
        end % Function Ends
        %------------------------------------------------------------------
        function [ImgDomainTriMatrix , ImgDomainTriID, ImgDomainNodes, BoundaryNodes, NonBoundaryNodes, BoundaryTriangles] = findDynamicElements(obj, DomainData, MeshData)
            %   Some elements within the reconstruction model are static and do
            %   not change while others are updated by the FEM-CSI algorithm at each iteration.
            %   Accordingly, this function finds all triangular mesh elements of the contrast profile
            %   that are dynamic and are permitted to change within the reconstruction model. 
            %   Input:        
            %           DomainData  : contains information about the problem being solved.
            %           MeshData  :    contains information about the mesh used by reconstruction model.
            %   Output:
            %           ImgDomainTriMatrix : Matrix containing the nodes of the triangular elements
            %                       within the reconstruction model.
            %           ImgDomainTriID : Vector containing the I.D. of the triangles within the
            %                      reconstruction model.
            %--------------------------------------------------------------
            x = MeshData.Nodes.xNodes;
            y = MeshData.Nodes.yNodes;
            NodesID = MeshData.Nodes.NodesID;
            nelematrix = MeshData.Triangles.NodesMatrix;

            switch DomainData.ImgDomainType
                case 'Rectangular'
                % Initializations
                    xmin = DomainData.xMin;
                    xmax = DomainData.xMax;
                    ymin = DomainData.yMin;
                    ymax = DomainData.yMax;

               % Find the triangles inside the imaging domain
               % Find the indices of the nodes with x-coordinates within the limits of xmin and xmax
                    x_index = find( x>=xmin & x<=xmax);
        
               % Find the indices of the nodes with y-coordinates within the limits of ymin and ymax
                    y_index = find( y>=ymin & y<=ymax);
        
               % Find the actual x- and y- coordinates of the nodes within the limits
               % x-limits and y-limits respectively
                    x_inside = x(x_index);
                    y_inside = y(y_index);
        
               % Finding the i.d. of the nodes within the limits
                    nn_inside1 = NodesID(x_index);
                    nn_inside2 = NodesID(y_index);
        
               % Find the intersection between nodes within the x-limits and nodes
               % within the y-limits
                    nn_inside = intersect(nn_inside1, nn_inside2);
                case 'CentCircle'
                    itemp1 = sqrt(MeshData.Nodes.xNodes.^2 + MeshData.Nodes.yNodes.^2);
                    itemp2 = find(itemp1 <= DomainData.xMin);
                    nn_inside = itemp2;
                case 'OffCentCir'
                    itemp1 = sqrt((MeshData.Nodes.xNodes - DomainData.ImgDomainCirxCent).^2 + ...
                        (MeshData.Nodes.yNodes - DomainData.ImgDomainCiryCent).^2);
                    itemp2 = find(itemp1 <= DomainData.ImgDomainCirRad);
                    nn_inside = itemp2;
                case 'Physics-Img'
                    % Elements with Physics == ImgDomainPhysics are part of the imaging domain
                    itmp = ismember(MeshData.Triangles.TriPhysics, DomainData.ImgDomainPhysics);
                    itmp = itmp(itmp ~= 0);
                    nn_inside = unique(MeshData.Triangles.NodesMatrix(itmp,:)); 
                case 'Physics-NonImg'
                    % Elements with Physics == NonImgDomainPhysics are not part of the imaging domain
                    itmp = ismember(MeshData.Triangles.TriPhysics, DomainData.NonImgDomainPhysics);
                    itmp = find(itmp == 0);
                    nn_inside = unique(MeshData.Triangles.NodesMatrix(itmp,:)); %#ok<FNDSB>
            end

            % Find the x- and y- coordinates of the node within [xmin, xmax, ymin,
            % ymax]
            x_ninside = x(nn_inside);
            y_ninside = y(nn_inside);

            % Find the triangular elements that contain the nodes within the
            % given x/y limits
            tmp_1a = ismember(nelematrix(:,1),nn_inside);
            tmp_2a = ismember(nelematrix(:,2),nn_inside);
            tmp_3a = ismember(nelematrix(:,3),nn_inside);

            tmp1 = find(tmp_1a == 1);
            tmp2 = find(tmp_2a == 1);
            tmp3 = find(tmp_3a == 1);

            tmp4 = union(tmp1,tmp2);    
            tmp5 = union(tmp3,tmp4);    

            tmp6 = intersect(tmp1,tmp2);
            tmp7 = intersect(tmp3,tmp6);

            tmp8a = ismember(tmp5, tmp7);
            tmp8 = tmp5(tmp8a ~= 1);

            BoundaryTriangles = tmp8;

            % Return nodes of triangles inside domain
            ImgDomainTriMatrix = nelematrix(tmp5,:);

            % Return i.d. of triangles
            ImgDomainTriID = tmp5;
            ImgDomainNodes = nn_inside;            

            % Finding Boundary nodes
            bd_nodes = unique([MeshData.Triangles.NodesMatrix(tmp8,1);MeshData.Triangles.NodesMatrix(tmp8,2);MeshData.Triangles.NodesMatrix(tmp8,3)]);
            tmp9a = ismember(ImgDomainNodes, bd_nodes);
            BoundaryNodes = ImgDomainNodes(tmp9a == 1);
            NonBoundaryNodes = ImgDomainNodes(tmp9a ~= 1);            
        end  % Function findDynamicElements      
        %------------------------------------------------------------------        
        function DielectricData = getEprBackground(obj, DomainData, MeshData, iFreq)
            % This function maps the background complex permittivity values onto mesh elements. 
            %   Input:
            %           DomainData   :  contains information about the problem being solved.
            %           MeshData      : Structure containing the information about the mesh
            %           iFreq             : Frequency Index
            %   Output:
            %           DielectricData: background complex permittivity map
            %           
            %--------------------------------------------------------------
            switch DomainData.MaterialFileType
                case {4,5,6}
                    % Type 4-6 - transform physics value of region (e.g, 3001, 3002, 3003)
                    % to the corresponding complex permittivity value (e.g, 31.30-j11.23,
                    % 5.50-j0.26, 41.46-j23.81).                    
                    DielectricData = obj.phy2eprNodes(DomainData, MeshData, iFreq);                
                case {8,9}
                    % Type 8-9 - The background permittivity is already
                    % represented as a complex permittivity map (similar to
                    % forward map). Therefore, just need to read the values
                    % from the material file which is the complex
                    % permittivity at each node.
                    DielectricData = dlmread(DomainData.MaterialFile);
                    [nrows ncols] = size(DielectricData); 
                    if ncols == 3
                        DielectricData = DielectricData(:,2) + 1i.* DielectricData(:,3);
                        elseif ncols == 2
                            DielectricData = DielectricData(:,1) + 1i.* DielectricData(:,2);
                    end
                    if DomainData.MaterialFileType == 8
                        itmp = find(real(DielectricData == 1));
                        if length(DomainData.RealBackground) > 1
                            epr_bkg = DomainData.RealBackground(iFreq) + 1i*DomainData.ImagBackground(iFreq);
                        else
                            epr_bkg = DomainData.RealBackground + 1i*DomainData.ImagBackground;
                        end
                            DielectricData(itmp) = epr_bkg;
                    end
                otherwise                    
                    % otherwise, then throw an error.                
                    obj.status = 'Error: Unknown materials file type. Check the function''s help.';
                    obj.errorFlag = true;
                    return;  
            end
        end % Function mapContrastValues        
        %------------------------------------------------------------------        
        function DielectricData = phy2eprNodes(obj, DomainData, MeshData, iFreq) 
            %  This function creates the dielectric profile from the information given in the
            %  materials file specified in DomainData. Dielectrics are specified at
            %  the nodes.
            %  Supported material files types
            %   
            %  There are three types of material files:
            %  Type 4: contains the real and imaginary contrasts of the scatterers
            %  Type 5: contains the relative permittivity and the conductivity of the
            %           scatterers.
            %  Type 6: contains the real and imaginary parts of the relative
            %           permittivity of the scatterers.
            %  Input:
            %       DomainData : Structure containing information about the problem
            %                domain.
            %       MeshData   : Structure containing information about the mesh of the
            %                domain.
            %  Output:
            %       DielectricData : Relative permittivity values at all the nodes of the
            %                    mesh.
            %--------------------------------------------------------------
            e0 = DomainData.e_o; % Permittivity of free-space
            MaterialMatrix = dlmread(DomainData.MaterialFile);
            MaterialFileType = DomainData.MaterialFileType;
            freq = DomainData.Freq(iFreq);

            if length(DomainData.RealBackground) > 1
                epr_bkg = DomainData.RealBackground(iFreq) + 1i*DomainData.ImagBackground(iFreq);
            else
                epr_bkg = DomainData.RealBackground + 1i*DomainData.ImagBackground;
            end

            DielectricData = epr_bkg .* ones(MeshData.Nodes.nNodes,1);

            switch MaterialFileType
            case 4
                for ii = 1:length(MaterialMatrix(:,1))
                    itemp = find(MeshData.Triangles.TriPhysics == MaterialMatrix(ii,1));
                    itemp2 = unique([MeshData.Triangles.NodesMatrix(itemp,1);MeshData.Triangles.NodesMatrix(itemp,2);MeshData.Triangles.NodesMatrix(itemp,3)]);
                    DielectricData(itemp2) = epr_bkg * (MaterialMatrix(ii, 2) + 1i*MaterialMatrix(ii, 3) + 1);
                end
            case 5
                for ii = 1:length(MaterialMatrix(:,1))
                    itemp = find(MeshData.Triangles.TriPhysics == MaterialMatrix(ii,1));
                    itemp2 = unique([MeshData.Triangles.NodesMatrix(itemp,1);MeshData.Triangles.NodesMatrix(itemp,2);MeshData.Triangles.NodesMatrix(itemp,3)]);
                    DielectricData(itemp2) = MaterialMatrix(ii, 2) - 1i*MaterialMatrix(ii, 3)/(2*pi*freq*e0);
                end
            case 6
                for ii = 1:length(MaterialMatrix(:,1))
                    itemp = find(MeshData.Triangles.TriPhysics == MaterialMatrix(ii,1));
                    itemp2 = unique([MeshData.Triangles.NodesMatrix(itemp,1);MeshData.Triangles.NodesMatrix(itemp,2);MeshData.Triangles.NodesMatrix(itemp,3)]);
                    DielectricData(itemp2) = MaterialMatrix(ii, 2) + 1i*MaterialMatrix(ii, 3);
                end
            end
        end % Function phy2eprNodes
        %------------------------------------------------------------------
        function imageFwdAndInvModels( obj, gui)
            % This is an helper function to image the masks for the real, imaginary,
            % and magnitude component of each tissue type (fat, transition,
            % fibroglandular, glandular, and malignant).             
            % Inputs:
            %   gui.dataPath - string - Points to results folder to store figures  
            %   gui.figFormat - string - Format to store figures (i.e., -fig, -png, -eps).
            %-------------------------------------------------------------- 
            figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.85, 0.85]);
            annotation('rectangle',[0.12,0.95,0.85,0.05],'Edgecolor','k');
            annotation('rectangle',[0.002,0.10,0.08,0.85],'Edgecolor','k');  
            set(0, 'DefaultAxesTitleFontWeight','normal');
            
            % Display permittivity maps
            % Forward model: Re - 1, Im - 2, Mag - 3
            subplot(2,3,1);imagesc(imrotate(real(obj.fwdMdl.complexPermittivity),90),'XData',obj.xNodes,'YData',-obj.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;               
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');            
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');             
            ylabel({'Forward model'; '(m)'});
            set(0, 'DefaultAxesTitleFontWeight','normal');
            title('Real', 'Units', 'normalized','Position',[0.5,1.125,0] );
            ylabh = get(gca,'yLabel');
            set(ylabh,'Position',get(ylabh,'Position') + [-0.010 0 0]); 
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.02 0.924 1.10 1.10]) % stretch its width and height                        
            drawnow     
            
            subplot(2,3,2);imagesc(imrotate(-imag(obj.fwdMdl.complexPermittivity),90),'XData',obj.xNodes,'YData',-obj.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;    
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');
            set(gca, 'yticklabel','');
            set(0, 'DefaultAxesTitleFontWeight','normal');
            title('Imaginary', 'Units', 'normalized','Position',[0.5,1.125,0] );            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.03 0.924 1.10 1.10]) % stretch its width and height                     
            drawnow
            
            subplot(2,3,3);imagesc(imrotate(abs(obj.fwdMdl.complexPermittivity),90),'XData',obj.xNodes,'YData',-obj.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;    
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            set(gca, 'xticklabel','');
            set(gca, 'yticklabel','');
            set(0, 'DefaultAxesTitleFontWeight','normal');
            title('Magnitude', 'Units', 'normalized','Position',[0.5,1.125,0] );            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.03 0.924 1.10 1.10]) % stretch its width and height                     
            drawnow
            
            % Display permittivity maps
            % FEM-CSI reconstruction (or inverse) model: Re - 1, Im - 2, Mag - 3
            subplot(2,3,4);imagesc(imrotate(real(obj.invMdl.complexPermittivity),90),'XData',obj.xNodes,'YData',-obj.yNodes);
            colormap(gca, gui.figColorMap); 
            ch = colorbar;
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);  
            xlabel('(m)');
            ylabel({'Inverse model'; '(m)'});            
            ylabh = get(gca,'yLabel');
            set(ylabh,'Position',get(ylabh,'Position') + [-0.010 0 0]);            
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.02 1.10 1.10 1.10]) % stretch its width and height                     
            drawnow            
            
            subplot(2,3,5);imagesc(imrotate(-imag(obj.invMdl.complexPermittivity),90),'XData',obj.xNodes,'YData',-obj.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;    
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            xlabel('(m)');
            set(gca, 'yticklabel','');               
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.03 1.10 1.10 1.10]) % stretch its width and height                     
            drawnow
            
            subplot(2,3,6);imagesc(imrotate(abs(obj.invMdl.complexPermittivity),90),'XData',obj.xNodes,'YData',-obj.yNodes);
            colormap(gca, gui.figColorMap);
            ch = colorbar;    
            set(gca,'fontsize',14);
            set(gca,'YDir','normal');
            p=get(gca,'position'); % save position            
            set(gca,'fontsize',14);
            xlabel('(m)');
            set(gca, 'yticklabel','');                     
            set(gca,'position',p); % restore position
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1.03 1.10 1.10 1.10]) % stretch its width and height                     
            drawnow
            
            switch gui.figFormat
                case '-DoNotSave'
                    % Do not save files
                case '-fig'
                    savefig([gui.dataPath 'figures/' 'complexPermittivityImages_all'])
                otherwise
                    export_fig([gui.dataPath 'figures/' 'complexPermittivityImages_all'], gui.figFormat);
            end
              
        end % function imageFwdAndInvModels
        %------------------------------------------------------------------        
        function [DomainData, varargout] = InputDataXML(obj, InputDataFilePath, InputDataFile)
            % InputDataXML
            %   DomainData = InputDataXML(InputDataFilePath, InputDataFile)
            %
            %   This function reads the InputDataFile to setup the input parameters
            %   for the FEM solver and any other solver to be created later.
            %   Input:
            %   -----------------------------------------------------------------------
            %
            %   InputDataFilePath : contains the file path containing the files for the
            %                       solver.
            %   InputDataFile     : contains the name of the input data file containing
            %                       the parameters of the FEM solver.
            %   Output:
            %   -----------------------------------------------------------------------
            %   DomainData      : Structure containing information about the FEM solver
            %                     parameters
            %
            %   Created by: Gabriel Faucher for Prof. Joe LoVetri Research Group
            %               University of Manitoba, Winnipeg, Canada.
            %
            %   $Revision: 1.0 $  $Date: 2014/05/07 $
            %   $Revision: 1.3 $  $Date: 2014/05/26 $
            %   $Revision: 1.4 $  $Date: 2014/06/20 10:14 a.m. $ - Added 'ABC-PEC' boundary option.
            %   $Revision: 1.5 $  $Date: 2014/07/14 12:11 p.m. $ - Noise is under CSICONFIG
            %   $Revision: 1.6 $  $Date: 2014/07/29 11:00 a.m. $ - Inhomogeneous background added under CSICONFIG
            %   $Revision: 1.7 $  $Date: 2014/07/31 11:36 a.m. $ - Realbackground duplicated in-case of multiple frequencies added under CSICONFIG
            %   $Revision: 1.8 $  $Date: 2014/08/10 12:00 p.m. $ - Multifrequency Inversion
            %   $Revision: 1.9 $  $Date: 2014/08/10 12:00 p.m. $ - Adding discrete-level regularization
            %   $Revision: 2.0 $  $Date: 2014/08/22 10:18 a.m. $ - Adding options for 3D

            %% Function Starts
            %  Reads XML file and returns a Document Object Model node representing the document.
            headNode = xmlread([InputDataFilePath InputDataFile]);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Preparing Forward Data Struct (DomainData) %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Read run type
            DomainData.RunType = char(headNode.getElementsByTagName('FEMPROJECT').item(0).getAttributeNode('run').getValue());

            % Read Problem Dimension
            dimensionNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getAttributeNode('dimension');
            if(~isempty(dimensionNode))
                DomainData.ProblemDimension = char(dimensionNode.getValue());  % Problem Dimension
            else % Default is 2D
                DomainData.ProblemDimension = '2D';
            end


            %% Reading the Frequency File: This contains information about the
            % frequencies of the simulation

            dataNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('DATADOMAIN').item(0);
            dataPath = char(dataNode.getAttributeNode('path').getValue());
            dataPath = InputDataFilePath;
            FrequencyFile = [dataPath char(dataNode.getAttributeNode('freqfilename').getValue())];
            Freq = dlmread(FrequencyFile);

            %% Assigning values to DomainData
            DomainData.e_o = 8.8541878176e-12;                        % Permittivity of free space
            DomainData.u_o = 4*pi*1e-7;                               % Premeability of free space
            DomainData.c_o = 1/sqrt(DomainData.e_o * DomainData.u_o); % Speed of light in free space
            DomainData.k_o = 2*pi*Freq/DomainData.c_o;                % Wave Number in free space
            DomainData.Z_o = sqrt(DomainData.u_o/DomainData.e_o);     % Wave impedance of free space
            DomainData.Freq = Freq;                                   % Frequencies
            DomainData.nFreq = length(Freq);                          % Number of frequency points

            extNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('EXTERNALPHYSICS').item(0);
            DomainData.RealBackground = 1.0;       % Default Real Background Permittivity
            attrNode = extNode.getAttributeNode('epsr_real');
            if(~isempty(attrNode))
                DomainData.RealBackground = str2double(attrNode.getValue());       % Real Background Permittivity
            end

            DomainData.ImagBackground = 0.0;       % Default Imaginary Background Permittivity
            attrNode = extNode.getAttributeNode('epsr_imag');
            if(~isempty(attrNode))
                DomainData.ImagBackground = str2double(attrNode.getValue());       % Imaginary Background Permittivity
            end

            DomainData.CondFlag = 0;               % default Conductivity Flag
            DomainData.CondBackground = 0;         % default Conductivity
            attrNode = extNode.getAttributeNode('cond');
            if(~isempty(attrNode))
                DomainData.CondBackground = str2double(attrNode.getValue());       % Conductivity of the background
                DomainData.CondFlag = 1;             % %if a conductivity value is given then set the conductivity flag.
            end

            attrNode = extNode.getAttributeNode('DL_epr');
            if(~isempty(attrNode))
                temp = strsplit(char(attrNode.getValue()), ';');
                [ii,jj] = size(temp); %#ok<ASGLU>
                DL_epr = [];
                for l = 1:jj
                    DL_epr = [DL_epr, str2double(temp{l})];  %#ok<AGROW>
                end
                DomainData.DL_epr = DL_epr;
            end

            % Reading permittivity from file if provided.
            attrNode = extNode.getAttributeNode('permFilePath');
            if(~isempty(attrNode))
                try
                    PermFilePath = char(attrNode.getValue());
                    try
                        tmp = interpol_salt_soln_perm([dataPath PermFilePath], Freq);
                    catch
                        tmp = interpol_salt_soln_perm(PermFilePath, Freq);
                    end
                    %%% tmp = water_permittivity(Freq);
                    DomainData.RealBackground = real(tmp);
                    DomainData.ImagBackground = imag(tmp);
                catch %#ok<*CTCH>
                    display('The relative permittivity was not read from a file.');
                end
            end


            if DomainData.CondFlag ~= 1
                DomainData.k_bkg = 2*pi*Freq.*sqrt(DomainData.RealBackground + 1i*DomainData.ImagBackground)./DomainData.c_o;
                DomainData.Z_bkg = sqrt(DomainData.u_o./(DomainData.e_o*(DomainData.RealBackground + 1i*DomainData.ImagBackground)));
                if length(DomainData.RealBackground) == 1
                    DomainData.Z_bkg = DomainData.Z_bkg * ones(DomainData.nFreq,1);
                end
            else
                DomainData.ImagBackground = - DomainData.CondBackground./(2*pi*Freq*DomainData.e_o);
                DomainData.k_bkg = 2*pi*Freq.*sqrt(DomainData.RealBackground - 1i * DomainData.CondBackground./(2*pi*Freq*DomainData.e_o))/DomainData.c_o;
                DomainData.Z_bkg = sqrt(DomainData.u_o./(DomainData.e_o*(DomainData.RealBackground + 1i*DomainData.ImagBackground)));
                DomainData.RealBackground = DomainData.RealBackground * ones(length(DomainData.ImagBackground),1);
            end

            % Type of Boundary
            bndNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('BOUNDARY').item(0);
            DomainData.BdType = char(bndNode.getAttributeNode('surfacetype').getValue()); % Type of Boundary

            switch DomainData.BdType
                case 'abc'
                    DomainData.BdType = 'ABC';
                    DomainData.ABCType = str2double(bndNode.getAttributeNode('order').getValue());  % Order of ABC
                    DomainData.ABContour = char(bndNode.getAttributeNode('contour').getValue());  % Contour type
                    DomainData.ContourRadius = str2double(bndNode.getAttributeNode('radius').getValue());  % Radius of the circular contour
                    DomainData.SourceType = str2double(dataNode.getAttributeNode('sourcetype').getValue());  % Type of Source

                case 'pec'
                    DomainData.BdType = 'PEC';
                    DomainData.SourceType = str2double(dataNode.getAttributeNode('sourcetype').getValue());  % Type of Source

                    if (strcmp(DomainData.ProblemDimension, '2D') && DomainData.SourceType == 7)
                        DomainData.BoundaryShape = char(bndNode.getAttributeNode('bdShape').getValue());
                        DomainData.TXpEdge = str2double(bndNode.getAttributeNode('txPerTxedge').getValue());
                        DomainData.BounadryLength = str2double(bndNode.getAttributeNode('bdLength').getValue());
                    end

                case {'pml-pec', 'pml-abc'}
                    DomainData.BdType = upper(DomainData.BdType);
                    DomainData.PMLPara = str2double(bndNode.getAttributeNode('absorber').getValue());    % Absorption parameter of PML
                    DomainData.PMLcoord = char(bndNode.getAttributeNode('coordinateType').getValue());   % PML coordinate type
                    DomainData.PMLthk   = str2double(bndNode.getAttributeNode('thickness').getValue());   % PML thickness
                    DomainData.SourceType = str2double(dataNode.getAttributeNode('sourcetype').getValue()); % Type of Source

                case 'abc-pec'
                    DomainData.BdType = 'ABC-PEC';
                    DomainData.ABCType = str2double(bndNode.getAttributeNode('order').getValue());  % Order of ABC
                    DomainData.ABContour = char(bndNode.getAttributeNode('contour').getValue());  % Contour type
                    DomainData.ContourRadius = str2double(bndNode.getAttributeNode('radius').getValue());  % Radius of the circular contour
                    DomainData.SourceType = str2double(dataNode.getAttributeNode('sourcetype').getValue()); % Type of Source

                    if (strcmp(DomainData.ProblemDimension, '2D') && DomainData.SourceType == 7)
                        DomainData.BoundaryShape = char(bndNode.getAttributeNode('bdShape').getValue());
                        DomainData.TXpEdge = str2double(bndNode.getAttributeNode('txPerTxedge').getValue());
                        DomainData.BounadryLength = str2double(bndNode.getAttributeNode('bdLength').getValue());
                    end

                    DomainData.PhysicsPEC = str2double(bndNode.getAttributeNode('pec_physics').getValue()); % Physics number of the PEC boundary
                    DomainData.PhysicsABC = str2double(bndNode.getAttributeNode('abc_physics').getValue()); % Physics number of the ABC boundary

                otherwise
                    error('Unknown type of boundary. Please check the help file.');
            end

%             DomainData.TxPath = [dataPath char(dataNode.getAttributeNode('txfilename').getValue())];
%             DomainData.RxPath = [dataPath char(dataNode.getAttributeNode('rxfilename').getValue())];
% 
%             % Read in Tx and Rx Locations
%             DomainData.TxPositions = dlmread(DomainData.TxPath);
%             DomainData.RxPositions = dlmread(DomainData.RxPath);
% 
%             [DomainData.nTx tmp] = size(DomainData.TxPositions); %#ok<NASGU>
%             [DomainData.nRx tmp] = size(DomainData.RxPositions); %#ok<NASGU>
%             DomainData.nRx = DomainData.nRx/DomainData.nTx;

            % Location of the mesh file
            meshNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('MESH').item(0);
            DomainData.MeshLocation = [dataPath char(meshNode.getAttributeNode('filename').getValue())];

            % Information to create the dielectric profile of the mesh
            % Location of the dielectric data file
            physNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('PHYSICS').item(0);
            DomainData.MaterialFile = [dataPath char(physNode.getAttributeNode('materialsFilename').getValue())];
            DomainData.MaterialFileType = str2double(physNode.getAttributeNode('materialsFileType').getValue());

            %% Type of Observation Point
            % interpolation: Extracted using connectivity matrix Ms.
            % huygen: Calculated using Huygens Principle (surface used is outer boundary surface of mesh)
            DomainData.ObsTypeStr = char(dataNode.getAttributeNode('obsType').getValue());
            if strcmpi('interpolation', DomainData.ObsTypeStr)
                DomainData.ObsType = 1;
            elseif strcmpi('huygen', DomainData.ObsTypeStr)
                DomainData.ObsType = 2;
            else
                error('Unknown type of observation point type. Please check the help file.');
            end


            %% Receiver Arrangment Type - 2D
            % 0 - Synthetic
            % 1 - MWI system at UofM
            % 2 - UPC
            % 3 - Fresnel
            % For 3D Fresnel is (1,2)
            DomainData.RxType = str2double(dataNode.getAttributeNode('rxArrangement').getValue());
            %% 2D Specific Settings
            % Added on 25th of April 2012
            if strcmp(DomainData.ProblemDimension, '2D')

                if DomainData.SourceType == 5
                    phiRxTypeNode = dataNode.getAttributeNode('phiReceiverType');
                    if (~isempty(phiRxTypeNode))
                        DomainData.PhiReceivers = str2double(phiRxTypeNode.getValue());
                    else
                        display('Phi source and using phi receivers!');
                        DomainData.PhiReceivers = 0;
                    end
                end

                % Added on 1st of Nov 2010
                if DomainData.MaterialFileType == 13
                    DomainData.AntennaPhys = str2double(physNode.getAttributeNode('AntennaPhys').getValue());
                    DomainData.AntennaSkip = str2double(physNode.getAttributeNode('AntennaSkip').getValue());
                    DomainData.AntennaBackground = str2double(physNode.getAttributeNode('AntennaBackgroundRe').getValue()) + 1i * str2double(physNode.getAttributeNode('AntennaBackgroundIm').getValue());
                    DomainData.FilePath = dataPath;
                end

            end

            %% 3D Specific Settings
            if strcmp(DomainData.ProblemDimension, '3D')
                % Polarization of 3D Sources
                switch DomainData.SourceType
                    case {4,5}
                        DomainData.J_t = str2num(dataNode.getAttributeNode('srcPolarization3D').getValue()); %#ok<ST2NM>
                end

                % Converting transmitter and receiver coordinates from cartesian to spherical
                xTx = DomainData.TxPositions(:,1); yTx = DomainData.TxPositions(:,2); zTx = DomainData.TxPositions(:,3);
                xRx = DomainData.RxPositions(:,1); yRx = DomainData.RxPositions(:,2); zRx = DomainData.RxPositions(:,3);

                rTx = sqrt(xTx.^2 + yTx.^2 + zTx.^2);
                phiTx = angle(xTx + 1i.*yTx);
                thetaTx = angle(zTx + 1i.*sqrt(xTx.^2 + yTx.^2));

                rRx = sqrt(xRx.^2 + yRx.^2 + zRx.^2);
                phiRx = angle(xRx + 1i.*yRx);
                thetaRx = angle(zRx + 1i.*sqrt(xRx.^2 + yRx.^2));

                DomainData.TxPositionsSph = [rTx phiTx thetaTx];
                DomainData.RxPositionsSph = [rRx phiRx thetaRx];

                DomainData.TxPositionsSphD = [rTx (180/pi)*phiTx (180/pi)*thetaTx];
                DomainData.RxPositionsSphD = [rRx (180/pi)*phiRx (180/pi)*thetaRx];

                [DomainData.nTx nCols] = size(DomainData.TxPositions);

                [nRows tmp] = size(DomainData.RxPositions); %#ok<NASGU>

                if nCols == 3
                    DomainData.nRx = nRows/DomainData.nTx; % Assuming same number of receivers for each transmitter
                else
                    DomainData.nRx = DomainData.TxPositions(:,4);
                end
                DomainData.TxPositions = DomainData.TxPositions(:,1:3);

            end
            if ~strcmp(DomainData.RunType, 'csi')
                InversionData = 0;
                 varargout{1} = InversionData;
                return
            end

            %% Inversion Settings for 2D FEMCSI Problems
            if (strcmp(DomainData.RunType, 'csi') && strcmp(DomainData.ProblemDimension, '2D'))

                % Inhomogeneous Background Options

                % Check if a node in the xml input file exists to use an inhomogeneous
                % background
                inhomogNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('CSICONFIG').item(0).getElementsByTagName('INHOMOGENEOUSBKG').item(0);

                if (~isempty(inhomogNode))

                    % Enable inhomogeneous background
                    DomainData.inhomogEnable = str2double(inhomogNode.getAttributeNode('enable').getValue());

                    if DomainData.inhomogEnable == 1

                        % Check if an alternative mesh for the inhomogeneous background
                        % has been provided. This mesh will be used to calculate the
                        % numerical incident fields used by the inversion algorithm.
                        inhomogMeshNode = inhomogNode.getAttributeNode('mesh');

                        if (~isempty(inhomogMeshNode))

                            % Flag used in the inversion algorithms to use the
                            % alternative inhomogeneous background mesh.
                            DomainData.inhomogDiffMesh = 1;
                            DomainData.inhomogMesh = [dataPath char(inhomogNode.getAttributeNode('mesh').getValue())];                      % Inhomogeneou background mesh
                            DomainData.inhomogMaterialFile = [dataPath char(inhomogNode.getAttributeNode('materialsFilename').getValue())]; % Inhomogeneous background materials file name
                            DomainData.inhomogMaterialFileType = str2double(inhomogNode.getAttributeNode('materialsFileType').getValue());  % Inhomogeneous background materials file type
                        else
                            DomainData.inhomogDiffMesh = 0;
                        end

                        % Check if a node exists that will use a map to create the
                        % materials file for the inhomogeneous background. This feature
                        % is used in imaging the human forearms
                        inhomogMapNode = inhomogNode.getAttributeNode('inhomogMapFilename');

                        if(~isempty(inhomogMapNode))

                            % Enable flag
                            DomainData.inhomogMapEnable = 1;

                            % Name of mapping file
                            DomainData.inhomogMapFile = [dataPath char(inhomogNode.getAttributeNode('inhomogMapFilename').getValue())];

                        else
                            % Enable flag
                            DomainData.inhomogMapEnable = 0;
                        end

                    end

                else
                    % Enable flag for inhomogeneous flag
                    DomainData.inhomogEnable = 0;
                end

                if DomainData.inhomogEnable == 0
                    DomainData.inhomogDiffMesh = 0;
                    DomainData.inhomogMapEnable = 0;
                end


                % Imaging domain coordinates: used when needed
                imgNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('CSICONFIG').item(0).getElementsByTagName('IMAGINGDOMAIN').item(0);
                imgDomType = char(imgNode.getAttributeNode('ImgDomainType').getValue());

                switch imgDomType
                    case 'UniformGrid'
                        DomainData.ImgDomainType = imgDomType;
                        DomainData.xMin = str2double(imgNode.getAttributeNode('xMin').getValue());
                        DomainData.xMax = str2double(imgNode.getAttributeNode('xMax').getValue());
                        DomainData.yMin = str2double(imgNode.getAttributeNode('yMin').getValue());
                        DomainData.yMax = str2double(imgNode.getAttributeNode('yMax').getValue());
                        DomainData.nx = str2double(imgNode.getAttributeNode('nx').getValue());
                        DomainData.ny = str2double(imgNode.getAttributeNode('ny').getValue());

                        DomainData.dx = (DomainData.xMax - DomainData.xMin)/DomainData.nx;
                        DomainData.dy = (DomainData.yMax - DomainData.yMin)/DomainData.ny;
                        DomainData.dA = DomainData.dx * DomainData.dy;

                        % Area of the imaging domain
                        DomainData.ImgDomainArea = (DomainData.xMax - DomainData.xMin) * (DomainData.yMax - DomainData.yMin);

                    case 'Rectangular'
                        DomainData.ImgDomainType = imgDomType;
                        DomainData.xMin = str2double(imgNode.getAttributeNode('xMin').getValue());
                        DomainData.xMax = str2double(imgNode.getAttributeNode('xMax').getValue());
                        DomainData.yMin = str2double(imgNode.getAttributeNode('yMin').getValue());
                        DomainData.yMax = str2double(imgNode.getAttributeNode('yMax').getValue());

                    case 'Circular'
                        DomainData.ImgDomainCirRad = str2double(imgNode.getAttributeNode('radius').getValue()); % Circle Radius
                        DomainData.ImgDomainCirxCent = str2double(imgNode.getAttributeNode('xCenter').getValue()); % Circle x-center
                        DomainData.ImgDomainCiryCent = str2double(imgNode.getAttributeNode('yCenter').getValue()); % Circle y-center
                        DomainData.xMin = DomainData.ImgDomainCirxCent - DomainData.ImgDomainCirRad;
                        DomainData.xMax = DomainData.ImgDomainCirxCent + DomainData.ImgDomainCirRad;
                        DomainData.yMin = DomainData.ImgDomainCiryCent - DomainData.ImgDomainCirRad;
                        DomainData.yMax = DomainData.ImgDomainCiryCent + DomainData.ImgDomainCirRad;

                        if (DomainData.ImgDomainCirxCent == 0 && DomainData.ImgDomainCiryCent == 0)
                            DomainData.ImgDomainType = 'CentCircle';
                        else
                            DomainData.ImgDomainType = 'OffCentCir';
                        end

                    case 'Physics-Img'
                        DomainData.ImgDomainType = imgDomType;
                        temp = strsplit(char(imgNode.getAttributeNode('physNum').getValue()), ';');
                        [i,j] = size(temp);
                        DomainData.ImgDomainPhysics = [];
                        for l = 1:j DomainData.ImgDomainPhysics = [DomainData.ImgDomainPhysics, str2double(temp{l})]; end

                    case 'Physics-NonImg'
                        DomainData.ImgDomainType = imgDomType;
                        temp = strsplit(char(imgNode.getAttributeNode('physNum').getValue()), ';');
                        [i,j] = size(temp);
                        DomainData.NonImgDomainPhysics = [];
                        for l = 1:j DomainData.NonImgDomainPhysics = [DomainData.NonImgDomainPhysics, str2double(temp{l})]; end

                    otherwise
                        error('Unknown type of imaging domain. Please check the help file.');
                end

            end
            %% Inversion Settings for 3D FEMCSI Problems
            if (strcmp(DomainData.RunType, 'csi') && strcmp(DomainData.ProblemDimension, '3D'))
                % Imaging domain coordinates: used when needed
                imgNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('CSICONFIG').item(0).getElementsByTagName('IMAGINGDOMAIN').item(0);
                imgDomType = char(imgNode.getAttributeNode('ImgDomainType').getValue());

                switch imgDomType
                    case 'Rectangular'
                        DomainData.ImgDomainType = imgDomType;
                        DomainData.xMin = str2double(imgNode.getAttributeNode('xMin').getValue());
                        DomainData.xMax = str2double(imgNode.getAttributeNode('xMax').getValue());
                        DomainData.yMin = str2double(imgNode.getAttributeNode('yMin').getValue());
                        DomainData.yMax = str2double(imgNode.getAttributeNode('yMax').getValue());
                        DomainData.zMin = str2double(imgNode.getAttributeNode('zMin').getValue());
                        DomainData.zMax = str2double(imgNode.getAttributeNode('zMax').getValue());

                        % Area of the imaging domain
                        DomainData.ImgDomainVolume = (DomainData.xMax - DomainData.xMin) * ...
                            (DomainData.yMax - DomainData.yMin) * ...
                            (DomainData.zMax - DomainData.zMin) ;

                    case 'Spherical'
                        DomainData.ImgDomainType = imgDomType;
                        DomainData.ImgDomainRadius = str2double(imgNode.getAttributeNode('radius').getValue()); % Sphere Radius
                        DomainData.ImgDomainVolume = (4/3) * pi * DomainData.ImgDomainRadius^3; % Sphere Volume

                    otherwise

                        error('Unknown type of imaging domain. Please check the help file.');
                end
            end

            %% Inversion Algorithm Parameters
            if strcmp(DomainData.RunType, 'csi')

                % Parameters for the inversion
                paramNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('CSICONFIG').item(0).getElementsByTagName('PARAMETERS').item(0);
                InversionData.nIterations = str2double(paramNode.getElementsByTagName('ITERATIONS').item(0).getAttributeNode('iter_count').getValue());

                % % Tolerance parameters
                % InversionData.dataMisfit  = str2double(paramNode.getElementsByTagName('TOLERANCE').item(0).getAttributeNode('value').getValue());
                % InversionData.diffSuccMisfit = str2double(paramNode.getElementsByTagName('TOLERANCE').item(0).getAttributeNode('diffSuccMisfit').getValue());
                % InversionData.diffSuccContrasts = str2double(paramNode.getElementsByTagName('TOLERANCE').item(0).getAttributeNode('diffSuccContrasts').getValue());

                % % Regularization Parameters
                % InversionData.RegType = char(paramNode.getElementsByTagName('REGULARIZATION').item(0).getAttributeNode('type').getValue());
                % InversionData.RegPara = str2double(paramNode.getElementsByTagName('REGULARIZATION').item(0).getAttributeNode('regPara').getValue());

                % Multiple Frequency Inversion
                InversionData.MultFreq = str2double(paramNode.getAttributeNode('multFreq').getValue());
                % % Multifrequency inversion options
                % 0: Multifrequency mode off
                % 1: Simultenous frequencies
                % 2: Marching-on-Frequency, 3: Marching-on-Frequency-and-Background

                switch InversionData.MultFreq

                    case 0

                    case 1
                        InversionData.MeasurementFileExtension = char(paramNode.getAttributeNode('measurementFileExtension').getValue());
                    case {2,3}
                        DomainData.Freqs                   = DomainData.Freq;
                        DomainData.nFreqs                  = DomainData.nFreq;
                        DomainData.k_o_AllFreqs            = DomainData.k_o;
                        DomainData.k_bkg_AllFreqs          = DomainData.k_bkg;
                        DomainData.Z_bkg_AllFreqs          = DomainData.Z_bkg;
                        DomainData.RealBackground_AllFreqs = DomainData.RealBackground;
                        DomainData.ImagBackground_AllFreqs = DomainData.ImagBackground;

                        InversionData.MeasurementFileExtension = char(paramNode.getAttributeNode('measurementFileExtension').getValue());
                        InversionData.nIterationsPerFreq = str2double(paramNode.getElementsByTagName('ITERATIONS').item(0).getAttributeNode('iter_count_per_freq').getValue());

                end


                % Location of Measurement Data
%                 InversionData.MeasurementFile = fullfile(dataPath, char(paramNode.getAttributeNode('measurementFile').getValue()));

                % Path for saving results
%                 InversionData.FilePath = fullfile(dataPath, char(headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('CSICONFIG').item(0).getElementsByTagName('OUTPUT').item(0).getAttributeNode('directory').getValue()));
% 
%                 % Try Making Directory
%                 if ~exist(InversionData.FilePath)
%                     try
%                         mkdir(InversionData.FilePath);
%                     catch
%                         error('Cannot create output directory! Check path.');
%                     end
%                 end


                % Check whether to apply a constraint or not
                InversionData.Constraint = str2double(paramNode.getElementsByTagName('CONSTRAINTS').item(0).getAttributeNode('applyConstraint').getValue());
                if (InversionData.Constraint == 1)
                    InversionData.maxConstRe = str2double(paramNode.getElementsByTagName('CONSTRAINTS').item(0).getAttributeNode('max_eps_r').getValue());
                    InversionData.minConstRe = str2double(paramNode.getElementsByTagName('CONSTRAINTS').item(0).getAttributeNode('min_eps_r').getValue());
                    InversionData.maxConstIm = str2double(paramNode.getElementsByTagName('CONSTRAINTS').item(0).getAttributeNode('max_eps_i').getValue());
                    InversionData.minConstIm = str2double(paramNode.getElementsByTagName('CONSTRAINTS').item(0).getAttributeNode('min_eps_i').getValue());
                end

                % Initial Guess
                if (~isempty(paramNode.getElementsByTagName('INITIALGUESS').item(0)))
                    InversionData.InitialGuess = str2double(paramNode.getElementsByTagName('INITIALGUESS').item(0).getAttributeNode('value').getValue());
                else
                    InversionData.InitialGuess = 1;
                end

                % Noise argument
                noiseNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('CSICONFIG').item(0).getElementsByTagName('NOISE').item(0);
                if (~isempty(noiseNode))
                    DomainData.NoiseFlag = str2double(noiseNode.getAttributeNode('addNoise').getValue());
                    DomainData.NoiseType = char(noiseNode.getAttributeNode('noiseType').getValue());
                    DomainData.NoisePercentage = str2double(noiseNode.getAttributeNode('noisePercent').getValue());
                else
                    DomainData.NoiseFlag = 0;
                end

                varargout{1} = InversionData;

            end

            if strcmp(DomainData.ProblemDimension, '3D')
                try
                    solverNode = headNode.getElementsByTagName('FEMPROJECT').item(0).getElementsByTagName('SOLVERSETUP').item(0);
                    DomainData.SolverType = char(solverNode.getAttributeNode('type').getValue());

                    if ~strcmp(DomainData.SolverType,'LU')
                        DomainData.MaxIter = str2double(solverNode.getAttributeNode('MaxIter').getValue());
                        DomainData.tol     = str2double(solverNode.getAttributeNode('tolerance').getValue());
                        DomainData.rest    = str2double(solverNode.getAttributeNode('rest').getValue());
                        DomainData.init    = char(solverNode.getAttributeNode('init').getValue());
                        DomainData.precond = char(solverNode.getAttributeNode('precond').getValue());

                        switch DomainData.precond
                            case 'BlkDiagonal'
                                DomainData.BlkSize = str2double(solverNode.getAttributeNode('blockSize').getValue());
                            case 'ilu'
                                DomainData.iluSetup.type =  char(solverNode.getAttributeNode('type').getValue());
                                DomainData.iluSetup.milu =  char(solverNode.getAttributeNode('milu').getValue());
                                DomainData.iluSetup.droptol = str2double(solverNode.getAttributeNode('droptol').getValue());
                        end
                    end

                catch
                    display('Solver should be specified.');
                    error('Solver should be specified.');
                end

            end
        end % function InputDataXML
%------------------------------------------------------------------
    end % methods
    %======================================================================
end % Class finiteElement
%==========================================================================
