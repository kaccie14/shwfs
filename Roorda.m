classdef Roorda < handle
    %MODEL roordalab Shack-Hartmann wavefront sensor
    %   MKS units assumed
    
    properties (Constant, Access = private)
        lensletFocalLength = 7.6e-3
    end
    
    properties (Access = private)
        bkgd % background image
        meas % measurment image w/ background subtracted
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
            
            obj.meas = imsubtract(rawMeasurement, obj.bkgd);
        end
        
        function show(obj)
            imagesc(obj.meas)
            axis off
            colormap gray
        end
        
    end
end

