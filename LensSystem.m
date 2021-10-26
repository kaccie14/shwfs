classdef LensSystem < handle
    %OPTICALSYSTEM models a paraxial optical system
    %   After constructing an instance of the class, individual surfaces
    %   are added by the user to model an optical system.
    %   By default:
    %   - constructed optical system will have an object (OBJ) and an image
    %   (IMA) surface with the latter fully absorbant. 
    %   - origin (z = 0) is defined at surface 2
    %
    %   Frame of reference originates from the vertex of the first surface,
    %   and light is assumed to travel from left to right. only paraxial
    %   (Gaussian) properties are supported.
    
    properties (Access = private)
        varNames = {'Surface', 'Radius', 'Thickness', 'Index',...
            'SemiDiameter', 'Conic'};
        varTypes = {'string', 'double', 'double', 'double', 'double',...
            'double'};
    end
    
%     properties(Constant, Access = private)
%     end
    
    properties (SetAccess = private, GetAccess = public)
        lensData % table
        lensDataReverse % reversed optical system has no image plane
        stop % [position, size] aperture-stop position and diameter (mm)
    end
    
    methods (Access = private)
        
        function reverse(obj)
            sz = [(height(obj.lensData) - 1) width(obj.lensData)];
            c = flipud(table2cell(obj.lensData));
            obj.lensDataReverse = table('Size', sz,...
                'VariableTypes', obj.varTypes,...
                'VariableNames', obj.varNames);
            obj.lensDataReverse(:,1) = cell2table(c(2:end,1)); % surface
            obj.lensDataReverse.Radius = -cell2mat(c(1:end-1,2));
            obj.lensDataReverse.Thickness = cell2mat(c(2:end,3));
            obj.lensDataReverse.Index = cell2mat(c(2:end,4));
            obj.lensDataReverse.SemiDiameter = cell2mat(c(1:end-1,5));
            obj.lensDataReverse.Conic = cell2mat(c(1:end-1,6));
            obj.lensDataReverse(end+1,:) = {'IMA', inf, 0, 1j*inf, 0, 0};
        end % this function might be removed
        
    end
    
    methods (Access = public)
        
        function obj = LensSystem()
            %OPTICALSYSTEM Construct an instance of this class
            %   Default object space is air with 0 thickness
            
            obj.lensData = table('Size', [2 6],...
                'VariableTypes', obj.varTypes,...
                'VariableNames', obj.varNames);
            obj.lensData(1,:) = {'OBJ', inf, 0, 1, inf, 0};
            obj.lensData(2,:) = {'IMA', inf, 0, 1j*inf, 0, 0};
            obj.lensDataReverse = cell2table({'OBJ', inf, 0, 1, inf, 0},...
                'VariableNames', obj.varNames);
            obj.stop.size = inf;
            obj.stop.position = 0;
        end
        
        function addSurface(obj, name, radius, thickness, index,...
                semiDiameter, conic)
            %ADDSURFACE adds optical surface to model
            %   Surfaces are always added from the right
            
            % Insert new surface just before IMA
            row = table(name, radius, thickness, index, semiDiameter,...
                conic, 'VariableNames', obj.varNames);
            ima = obj.lensData(end,:);
            obj.lensData(end,:) = row;
            obj.lensData = [obj.lensData; ima];
            
            % peg OBJ subtense to smallest surface clear aperture
            if obj.lensData(1,:).SemiDiameter > semiDiameter
                obj.lensData(1,:).SemiDiameter = semiDiameter;
            end   
            obj.reverse(); % Update reverse lens data table
        end
        
%         function addStop(obj, position, semiDiameter)   
%         end  
        
        function updateFirstSurface(obj, name, radius, thickness, index,...
                semiDiameter, conic)
            %UPDATEFIRSTSURFACE replaces first surface's properties with
            %whatever the user inputs
            arguments
                obj
                name
                radius (1,1) double = obj.lensData.Radius(1)
                thickness (1,1) double = obj.lensData.Thickness(1)
                index (1,1) double = obj.lensData.Index(1)
                semiDiameter (1,1) double = obj.lensData.SemiDiameter(1)
                conic (1,1) double = obj.lensData.Conic(1)
            end
            
            obj.lensData(1,:) = table(name, radius, thickness, index,...
                semiDiameter, conic, 'VariableNames', obj.varNames);
            obj.reverse(); % reverse version needs updating            
        end
        
        function updateLastSurface(obj, name, radius, thickness, index,...
                semiDiameter, conic)
            %UPDATELASTSURFACE replaces last surface's properties with
            %whatever the user inputs
            %   Sometimes one desires to adjust the image surface 
            
            obj.lensData(end,:) = table(name, radius, thickness, index,...
                semiDiameter, conic, 'VariableNames', obj.varNames);
            obj.reverse(); % reverse version needs updating
        end
              
        function H = principalPlanes(obj, s0, s1)
            % Compute principal planes. Results are given with respect to
            % 1st and last optical surfaces for H1 and H2 respectively
            % (i.e. H < 0 means plane is located INSIDE optical system)
            arguments
                obj
                s0 (1,1) {mustBeNumeric}
                s1 (1,1) {mustBeNumeric} = height(obj.lensData) - 1
            end
            n = height(obj.lensData); % total number of surfaces
            v0 = [1; 0]; % parallel ray input
            
            % Find principal planes
            M = obj.rayTransferMatrix(s0, s1);
            v1 = M * v0;
            H.secondary = (1 - v1(1)) / tan(v1(2));
            M = obj.rayTransferMatrix(n-s1+1, n-s0+1, true);
            v1 = M * v0;    
            H.primary = (1 - v1(1)) / tan(v1(2));
        end
        
        function M = rayTransferMatrix(obj, s0, s1, reverse)
            %RAYTRACE praxial ray trace
            % M: ray transfer matrix from surface s
            arguments
                obj
                s0 (1,1) {mustBeNumeric}
                s1 (1,1) {mustBeNumeric} = height(obj.lensData) - 1
                reverse (1,1) logical = false
            end
            
            M = eye(2);
            if s0 < 2 % skip object space
                return
            else
                if reverse
                    l = obj.lensDataReverse;
                else
                    l = obj.lensData;
                end
                
                for s = s0:s1
                    R = l.Radius(s);
                    n0 = real(l.Index(s-1));
                    n = real(l.Index(s));
                    if s < s1
                        d = l.Thickness(s);
                    else
                        d = 0;
                    end
                    M = [1 d; 0 1]*[1 0; (1/R)*(n0/n - 1) n0/n]*M;
                end
                
            end   
        end
         
        
    end% Public methods 
    
    

    
    
end

