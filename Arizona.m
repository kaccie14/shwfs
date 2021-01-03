classdef Arizona < handle
    %ARIZONA Arizona eye model
    %   Paraxial eye models can estimate cardinal points, pupils,
    %   magnification, etc. (first-order effects). Arizona eye model is
    %   designed to also match clinical levels of aberration and
    %   incorporate slight variations due to accommodation.
    %
    %   Arizona model doesn't include semi-diameters of each surface
    %   Origin is defined at corneal apex
    
    properties (Constant, Access = private)
        sd = [6 6 5 5 0.5] % semi-diameters (mm)
        cornealThickness = 0.55; % central thickness of cornea (mm)
    end
    
    properties (Access = private)
        lensSystem
        aqueousThickness = 2.97; % unaccomodated aqueous thickness
    end
    
    properties (SetAccess = private, GetAccess = public)
        data
        stop % iris position and size
    end
    
    methods (Access = public)
        function obj = Arizona(A)
            %ARIZONA Construct an instance of Arizona Eye Model
            %   Initialize lens system with Arizona eye properties with
            %   user input diopters of accommodation (A).
            
            ct = obj.cornealThickness;
            at = obj.aqueousThickness - 0.04 * A;
            
            % Optical parameters are unitless or in mm
            ls = LensSystem();
            ls.addSurface("Cornea", 7.8, ct, 1.377, obj.sd(1), -0.25)
            ls.addSurface("Aqueous", 6.5, at, 1.337, obj.sd(2), -0.25)
            ls.addSurface("Lens", 12.0 - 0.4*A, 3.767 + 0.04*A,...
                1.42 + 0.00256*A - 0.00022*A^2, obj.sd(3), -7.518749 + 1.285720*A)
            ls.addSurface("Vitreous", -5.224557 + 0.2*A, 16.713, 1.336,...
                obj.sd(4), -1.353971 - 0.431762*A)
            ls.updateLastSurface("Retina", -13.4, 0, 0, obj.sd(5), 0)
            
            % Initialize eye-model object properties
            obj.lensSystem = ls;
            obj.aqueousThickness = at;
            obj.data = ls.lensData;
            obj.stop.size = 3;
            obj.stop.position = ct + at;
        end
        
        function a = anteriorChamberDepth(obj)
            a = obj.aqueousThickness + obj.cornealThickness;          
        end
        
        function l = axialLength(obj)
            l = sum(obj.data.Thickness);
        end
        
        function m = entrancePupil(obj)
            m = 1;
        end
          
        function draw(obj)
            
            % Optical axis
            p = 5; % padding on either side of plot            
            z_axis.x = [-p (sum(obj.data.Thickness)+p)];
            z_axis.y = [0 0];
            figure
            plot(z_axis.x, z_axis.y, 'k', 'LineWidth', 3)
            hold on
            
%             sd = obj.data.SemiDiameter(2);
%             xc = 0 + R; % yc = 0 always
%             theta = ceil(atand(6/R));
%             t = (-theta:theta)';
%             x = xc - R * cosd(t);
%             y = R * sind(t);        

            % Anterior cornea
            R2 = obj.data.Radius(2);
            K2 = obj.data.Conic(2);           
            acd = obj.aqueousThickness + obj.cornealThickness;
            z2 = (0:0.1:acd)'; % extends out to whatever corresponds to ACD
            y2 = sqrt(R2^2 - (z2*(K2 + 1) - R2).^2) / sqrt(K2 + 1);
            plot([z2 z2], [y2 -y2], 'k')
            
            % Posterior cornea
            t2 = obj.data.Thickness(2);
            R3 = obj.data.Radius(3);
            K3 = obj.data.Conic(3);   
            y3 = (-max(y2):0.1:max(y2))';
            z3 = t2 + sag_nontoric(y3, R3, K3);
            plot(z3, y3, 'k')
            
            % Anterior lens
            t3 = obj.data.Thickness(3);
            R4 = obj.data.Radius(4);
            K4 = obj.data.Conic(4);
            y4 = (-obj.data.SemiDiameter(4):0.1:obj.data.SemiDiameter(4))';
            z4 = t2 + t3 + sag_nontoric(y4, R4, K4);  
            plot(z4, y4, 'k')
            
            % Posterior lens
            t4 = obj.data.Thickness(4);
            R5 = obj.data.Radius(5);
            K5 = obj.data.Conic(5);
            z5 = t2 + t3 + t4 + sag_nontoric(y4, R5, K5);  
            plot(z5, y4, 'k')
            
            % Retina (fovea)
            t5 = obj.data.Thickness(5);
            R6 = obj.data.Radius(6);            
            y6 = (-obj.data.SemiDiameter(6):0.1:obj.data.SemiDiameter(6))';
            z6 = t2 + t3 + t4 + t5 + sag_nontoric(y6, R6);
            plot(z6, y6, 'k')
            
            % Outer sclera
            
            hold off
            xlim(z_axis.x)
            xlabel("z (mm)"), ylabel("y (mm)")
            grid on, axis equal
        end       
    end
end

function Z = sag_nontoric(Y, R, K)
arguments
    Y (:, 1) double % Radial distance
    R (1,1) double % Radius of curvature
    K (1,1) double = 0 % Conic constant (unitless)
end

Z = ((1/R) * Y.^2) ./ (1 + sqrt(1 - (1 + K) * (Y.^2) / R^2));
end


