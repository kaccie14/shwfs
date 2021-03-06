classdef Roorda < handle
    %MODEL roordalab Shack-Hartmann wavefront sensor
    %   Comform to ANSI Z80.28 Method of Reporting Optical Aberrations of
    %   Eyes.
    
    properties (Constant, Access = private)
        lensletFocalLength = 7.6e-3 % meters
        lensletPitch = 0.3e-3 % meters
        pixelPitch = 7.5e-6 % meters
        magnification = 0.8
        
        scaleForProcessing = 0.25
        thresholdSensitivity = 0.3
        theta = linspace(0, 2*pi, 360)
    end
    
    properties %(Access = private)
        bkgd % background image
        meas % measurment image w/ background subtracted
        
        pupil = struct('x', 0, 'y', 0, 'r', 0) % pixels
    end
    
    methods (Access = public)
        
        function obj = Roorda(background)
            %MODEL Construct an instance of this class
            %   Require a background image to instantiate
            obj.bkgd = background;
        end
        
        function outputArg = process(obj, rawMeasurement)
            %PROCESS main entry to process spot-field image
            %   rawMeasurement is assumed to be acquired by the same device
            %   as that of the background image
            outputArg = 5;
            
            obj.meas = uint8(rawMeasurement - obj.bkgd);
            
            meas_scaled = imresize(obj.meas, obj.scaleForProcessing,...
                "AntiAliasing", true);
            
            
            
            
            meas_bin = imbinarize(meas_scaled, "adaptive",...
                "Sensitivity", obj.thresholdSensitivity);
            [y, x] = find(meas_bin);
            [c, r] = minboundcircle(x, y);
            obj.pupil.x = c(1) / obj.scaleForProcessing;
            obj.pupil.y = c(2) / obj.scaleForProcessing;
            obj.pupil.r = r / obj.scaleForProcessing;
        end
        
        function show(obj)
            x = obj.pupil.x + obj.pupil.r * cos(obj.theta);
            y = obj.pupil.y + obj.pupil.r * sin(obj.theta);
            
            imagesc(obj.meas)
            hold on
            plot(x, y, 'g')
            axis off
            colormap gray
            hold off
        end
        
    end
end

