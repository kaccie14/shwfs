classdef Arizona < handle
    %ARIZONA Arizona eye model
    %   Paraxial eye models can estimate cardinal points, pupils,
    %   magnification, etc. (first-order effects). Arizona eye model is
    %   designed to also match clinical levels of aberration and
    %   incorporate slight variations due to accommodation. The Arizona
    %   model doesn't include semi-diameters of each surface
    
    properties (Constant, Access = private)
        y = [6 6 5 5 0.5] % semi-diameters (mm)
    end
    
    properties (Access = private)
        lensSystem
    end
    
    properties (Access = public)
        data
    end
    
    methods
        function obj = Arizona(A)
            %ARIZONA Construct an instance of Arizona Eye Model
            %   Initialize lens system with Arizona eye properties with
            %   user input diopters of accommodation (A).
            
            % Surface parameters are either unitless or in mm
            ls = LensSystem();
            ls.addSurface("Cornea", 7.8, 0.55, 1.377, obj.y(1), -0.25)
            ls.addSurface("Aqueous", 6.5, 2.97 - 0.04*A, 1.337, obj.y(2), -0.25)
            ls.addSurface("Lens", 12.0 - 0.4*A, 3.767 + 0.04*A,...
                1.42 + 0.00256*A - 0.00022*A^2, obj.y(3), -7.518749 + 1.285720*A)
            ls.addSurface("Vitreous", -5.224557 + 0.2*A, 16.713, 1.336,...
                obj.y(4), -1.353971 - 0.431762*A)
            ls.updateLastSurface("Retina", -13.4, 0, 0, obj.y(5), 0)
            
            obj.lensSystem = ls;
            obj.data = ls.lensData;
        end
        
        function draw(obj)
            p = 5; % padding on either side of plot            
            z_axis.x = [-p (sum(obj.data.Thickness)+p)];
            z_axis.y = [0 0];
            
            figure
            
            % Sample full circle with 360 points
            R = obj.data.Radius(2);
            K = obj.data.Conic(2);
            sd = obj.data.SemiDiameter(2);
 
            xc = 0 + R; % yc = 0 always
            
     
            theta = ceil(atand(6/R));
            t = (-theta:theta)';
            x = xc - R * cosd(t);
            y = R * sind(t);
            
 
            plot(z_axis.x, z_axis.y, 'k', 'LineWidth', 3)           
            hold on, plot(x, y, 'b')
            
            y = (-6:0.5:6)';
            x = ((1/R) * y.^2) ./ (1 + sqrt(1 - (1 + K) * (y.^2) / R^2));
            
            
            plot(x, y, 'r')
            hold off
            xlim(z_axis.x)
            xlabel("z (mm)"), ylabel("y (mm)")
            grid on, axis equal
        end
        
                
    end
end

