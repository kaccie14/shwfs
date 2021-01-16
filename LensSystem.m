classdef LensSystem < handle
    %OPTICALSYSTEM models a paraxial optical system
    %   After constructing an instance of the class, individual surfaces
    %   are added by the user to model an optical system.
    %   By default:
    %   - constructed optical system will have an object (OBJ) and an image
    %   (IMA) surface with the latter fully absorbant. 
    %   - stop is in
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
            
            % peg OBJ subtense to smallest surface clear apertur
            if obj.lensData(1,:).SemiDiameter > semiDiameter
                obj.lensData(1,:).SemiDiameter = semiDiameter;
            end 
              
            % Update reverse lens data table
            obj.reverse();
        end
        
        function addStop(obj, position, semiDiameter)
            
        end  
        
        function updateFirstSurface(obj, name, radius, thickness, index,...
                semiDiameter, conic)
            %UPDATEFIRSTSURFACE replaces first surface's properties with
            %whatever the user inputs
            
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
        
    end% Public methods    
end

