classdef Arizona < handle
    %ARIZONA Arizona (AZ) eye model
    %   Paraxial eye models can estimate cardinal points, pupils,
    %   magnification, etc. (first-order effects). The AZ model also
    %   matches clinical levels of aberration and models accommodation.
    %
    %   Following extensions are not part of the AZ model but are included:
    %   - Semi-diameters of each optical surface
    %   - Different retinal curvatures for foveal and nonfoveal retina
    %
    %   Origin is defined at corneal apex
    %
    %   retina: result only used for sclera modeling
    
    properties (Constant, Access = private)
        sd = [6 6 5 5 0.5] % semi-diameters (mm)
        cornealThickness = 0.55; % central thickness of cornea (mm)
        vitreousThickness = 16.713; % mm
        acd = 3.52; % nominal anterior chamber depth (mm)
    end
    
    properties (Access = private)
        lensSystem
        aqueousThickness = 2.97; % nominal aqueous thickness (mm)
    end
    
    properties (SetAccess = private, GetAccess = public)
        data
        stop % iris position and size
        
        % Characteristics R, K, (Cz, Cy), (Ez, Ey)
        cornea; % anterio cornea assumed
        corneaCreated = false;
        retina; % retinal shape (non-foveal)
        retinaCreated = false;
        sclera; % (outer) scleral shape
    end
    
    methods (Access = public)
        
        function obj = Arizona(A)
            %ARIZONA Construct an instance of Arizona Eye Model
            %   Initialize lens system with Arizona eye properties with
            %   user input diopters of accommodation (A).
            
            ct = obj.cornealThickness;
            at = obj.aqueousThickness - 0.04 * A;
            vt = obj.vitreousThickness;
            
            % Optical parameters are unitless or in mm
            ls = LensSystem();
            ls.addSurface("Cornea", 7.8, ct, 1.377, obj.sd(1), -0.25)
            ls.addSurface("Aqueous", 6.5, at, 1.337, obj.sd(2), -0.25)
            ls.addSurface("Lens", 12.0 - 0.4*A, 3.767 + 0.04*A,...
                1.42 + 0.00256*A - 0.00022*A^2, obj.sd(3), -7.518749 + 1.285720*A)
            ls.addSurface("Vitreous", -5.224557 + 0.2*A, vt, 1.336,...
                obj.sd(4), -1.353971 - 0.431762*A)
            ls.updateLastSurface("Retina", -13.4, 0, 0, obj.sd(5), 0)
            
            % Initialize eye-model object properties
            obj.lensSystem = ls;
            obj.aqueousThickness = at;
            obj.data = ls.lensData;
            obj.stop.size = 3;
            obj.stop.position = ct + at;
            
            % Generate eye shape
            obj.createCornea(); % anterior cornea
            obj.createRetina(); % retina (non-foveal)
            obj.createSclera(); % outer sclera assumed
        end
        
        function createCornea(obj)
            R = obj.data.Radius(2);
            K = obj.data.Conic(2);
            obj.cornea.Cz = 0;
            obj.cornea.Cy = 0;
            obj.cornea.Ez = obj.acd;
            obj.cornea.Ey = radial(obj.acd, R, K);
            obj.cornea.R = R;
            obj.cornea.K = K;
            obj.corneaCreated = true;
        end
        
        function createRetina(obj)
            if ~obj.corneaCreated
                error("Cannot create retina before cornea is created")
            end
            y0 = obj.cornea.Ey;
            z0 = sag_nontoric(y0, obj.data.Radius(3), obj.data.Conic(3),...
                obj.cornealThickness);
            z1 = obj.axialLength; % foveal plane (wrt corneal apex)
            y1 = obj.sd(end); % Radial position of foveal edg
            R = obj.data.Radius(end); % Foveal radius of curvature
            sag = sag_nontoric(y1, R);
            s = (R - sag) / sqrt(R^2 - (R - sag)^2);
            m = -1/s; % slope of line normal to retina at (positive) edge
            b = -m * (obj.axialLength + R);
            
            % Center and radius of curvature
            obj.retina.Cz = (z1^2 + y1^2 - z0^2 - y0^2 + 2*b*(y0-y1)) / ...
                (2*(z1 - z0) + 2*m*(y1 - y0));
            obj.retina.Cy = m * obj.retina.Cz + b;
            obj.retina.R = sqrt((z0 - obj.retina.Cz)^2 +...
                (y0 - obj.retina.Cy)^2);
            obj.retina.K = 0;
            obj.retinaCreated = true;
        end
        
        function createSclera(obj)
            % Outer sclera assumed because that's usually the case when
            % referring to the "white" of the eye
            if ~obj.retinaCreated
                error("Cannot create sclera before retina is created")
            end
            
            % Anterior cornea
            za = obj.acd;
            ya = radial(za, obj.data.Radius(2), obj.data.Conic(2));
            
            % Outer sclera
            obj.sclera.Cz = obj.retina.Cz;
            obj.sclera.Cy = 0;
            obj.sclera.R = sqrt((za - obj.sclera.Cz)^2 + ya^2);
            obj.sclera.Ez = obj.sclera.Cz + obj.sclera.R;          
            obj.sclera.K = 0;
        end
        
        function c = centroid(obj)
            [zc, c] = obj.centroidCornea;
            [zs, s] = obj.centroidSclera;           
            c = (zc + zs) / (c + s);  
        end
        
        function a = anteriorChamberDepth(obj)
            a = obj.aqueousThickness + obj.cornealThickness;
        end
        
        function l = axialLength(obj)
            T = obj.data.Thickness;
            l = sum(T(2:end));
        end
        
        %% TODO
        %1/10 defining stops in lenSystem
        
        function m = entrancePupil(obj)
            m = 1;
            
        end
        
        function m = exitPupil(obj)
            %obj.lensSystem.trace(from, to)
            
            m = 1;
            
            % g
            
            
            
        end
        
        
        %% end
        
        function draw(obj, internal)
            arguments
                obj
                internal (1,1) logical = false
            end
            p = 5; % padding on either side of plot
            dz = 0.01;
            dy = 0.1;
            
            % Optical axis
            z_axis.x = [-p (sum(obj.data.Thickness)+p)];
            z_axis.y = [0 0];
            figure;
            plot(z_axis.x, z_axis.y, 'k')
            hold on
            
            % Anterior cornea
            R2 = obj.data.Radius(2);
            K2 = obj.data.Conic(2);
            z2 = (0:dz:obj.acd)'; % extends out to whatever corresponds to ACD
            y2 = radial(z2, R2, K2);
            plot([z2 z2], [y2 -y2], 'r')
            
            if internal
                % Posterior cornea
                t2 = obj.data.Thickness(2);
                R3 = obj.data.Radius(3);
                K3 = obj.data.Conic(3);
                y3 = (-max(y2):dy:max(y2))';
                z3 = sag_nontoric(y3, R3, K3, t2);
                plot(z3, y3, 'r')
                
                % Anterior lens
                t3 = obj.data.Thickness(3);
                R4 = obj.data.Radius(4);
                K4 = obj.data.Conic(4);
                y4 = (-obj.data.SemiDiameter(4):dy:obj.data.SemiDiameter(4))';
                z4 = sag_nontoric(y4, R4, K4, t2 + t3);
                plot(z4, y4, 'r')
                
                % Posterior lens
                t4 = obj.data.Thickness(4);
                R5 = obj.data.Radius(5);
                K5 = obj.data.Conic(5);
                z5 = sag_nontoric(y4, R5, K5, t2 + t3 + t4);
                plot(z5, y4, 'r')
                
                % Fovea
                R6 = obj.data.Radius(6);
                l = obj.axialLength;
                y6 = (-obj.data.SemiDiameter(6):dy:obj.data.SemiDiameter(6))';
                z6 = sag_nontoric(y6, R6, 0, l);
                plot(z6, y6, 'r')
                
                % Non-foveal retina
                or = obj.retina.Cz - obj.retina.R;
                z_end = round(l + sag_nontoric(max(y6),R6,l), 1);
                zr = (max(z3):dz:z_end)';
                yr = obj.retina.Cy + radial(zr, obj.retina.R,...
                    obj.retina.K, or);
                plot([zr zr], [yr -yr], 'k')
            end
            
            % Sclera
            os = obj.sclera.Cz - obj.sclera.R;
            zs = (obj.acd:dz:ceil(obj.sclera.Ez))';
            ys = radial(zs, obj.sclera.R, obj.sclera.K, os);
            plot([zs zs], [ys -ys], 'k')
            
            hold off, xlim(z_axis.x)
            xlabel("z (mm)"), ylabel("y (mm)")
            grid on, axis equal
        end
        
        function verification(obj)
            % Verify eyeball centroid calculation (issue #11)
            [L, Z] = obj.matte();
            c = sum(sum(Z.*L)) / sum(sum(L)); % numerical centroid
            if (abs(obj.centroid - c) < 0.1)            
                fprintf(['Pass: Eyeball centroid position (z) agrees',...
                    'with numerical estimate within 0.1 mm, %1.2f mm'], c);
            else
                fprintf("Fail: Eyeball centroid");
            end
            %imagesc(L); axis image, colormap gray, axis off
            
            % Verify pupil location

            % Scleral thickness 0.50 +/- 0.05 mm (tentative)
            t = obj.scleralThickness;
        end
    end
    
    %% Private methods
    
    methods (Access = private)
              
        function [L, Z] = matte(obj)
            % Generat a binary image of eyeball
            z = 0:0.1:obj.sclera.Ez;
            y = zeros(size(z));
            ind_acd = sum(z < obj.acd);
            y(1:ind_acd) = radial(z(1:ind_acd), obj.cornea.R, obj.cornea.K);
            y(ind_acd+1:end) = radial(z(ind_acd+1:end), obj.sclera.R,...
                obj.sclera.K, obj.sclera.Cz - obj.sclera.R);        
            Y = abs((-max(y):0.1:max(y))' * ones(1,numel(z)));
            L = false(size(Y));
            for ind = 1:size(Y,2)
                c = Y(:,ind) < y(ind);
                L(:,ind) = c;
            end
            Z = ones(size(Y,1),1) * z;
        end
        
        function [c, m] = centroidCornea(obj)
            R = obj.cornea.R;
            K = obj.cornea.K;
            
            % limits of integration
            u1 = 3.52 * (K + 1) - R;
            u0 = -R;
            
            % centroid integral formulas
            n = @(u) (1/(K+1)^(5/2)) * ((R/2)*(u*sqrt(R^2 - u^2) +...
                (R^2)*asin(u/R)) - (1/3)*(R^2 - u^2)^(3/2));
            d = @(u) (0.5/(K+1)^(3/2)) * (u*sqrt(R^2 - u^2) +...
                (R^2)*asin(u/R));
            
            % characteristic function
            c = n(u1) - n(u0);
            
            % measure
            m = d(u1) - d(u0);
        end
        
        function [c, m] = centroidSclera(obj)
            R = obj.sclera.R;
            Cz = obj.sclera.Cz;
            u0 = obj.acd - Cz; %u = z - Cz
            u1 = R;       

            % centroid integral formulas
            n = @(u)  (Cz/2) * (u*sqrt(R^2 - u^2) + (R^2)*asin(u/R)) -...
                (1/3)*(R^2 - u^2)^(3/2);
            d = @(u) 0.5 * (u*sqrt(R^2 - u^2) + (R^2)*asin(u/R));
            
            % characteristic function
            c = n(u1) - n(u0);
            
            % measure
            m = d(u1) - d(u0);      
        end
        
        function t = scleralThickness(obj)
            % Anterior cornea
            za = obj.acd;
            ya = radial(za, obj.data.Radius(2), obj.data.Conic(2));
            
            % Posterior cornea
            yb = ya; % this will most likely change in future
            zb = sag_nontoric(yb, obj.data.Radius(3), obj.data.Conic(3),...
                obj.data.Thickness(2));
            
            % Scleral thickness
            t = sqrt((za - zb)^2 + (ya - yb)^2);
        end
         
    end
    
end


%% Helper functions

function Z = sag_nontoric(Y, R, K, O)
arguments
    Y (:,1) double % Radial distance
    R (1,1) double % Radius of curvature
    K (1,1) double = 0 % Conic constant (unitless)
    O (1,1) double = 0 % z-axis offset
end

Z = O + ((1/R) * Y.^2) ./ (1 + sqrt(1 - (1 + K) * (Y.^2) / R^2));
end

function Y = radial(Z, R, K, O)
% Radial distance as from sag; essentially the inverse of "sag_nontoric"
arguments
    Z (:,1) double % Sag
    R (:,1) double % Radius of curvature
    K (1,1) double = 0 % Conic constant (unitlessw)
    O (1,1) double = 0 % z-axis offset
end

Z = Z - O; % apply offset
Y = real(sqrt(R^2 - (Z * (K + 1) - R).^2) / sqrt(K + 1));
end
