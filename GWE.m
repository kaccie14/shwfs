classdef GWE
    %GWE image processing attributed to Gonzalez, Woods and Eddings
    %   Custom thresholding and morphology
    %   This helper class modularizes self-contained calculations into
    %   private methods called by user-facing public methods
    
    %% Private Methods
    
    methods(Static, Access = private)
        
        function gd = imdilate(f, sz)
            % Image (f) assumed to be 8 bit unsigned integer
            % Structuring element dimensions (sz) assumed to be odd
            % Structuring element assumed to be flat
            
            w = uint8(ones(sz));
            [rs, cs] = size(f);
            a = floor(0.5*sz);
            gd = f;
            
            for r = (1+a):(rs-a)
                for c = (1+a):(cs-a)
                    gd(r,c) = max(f(r-a:r+a, c-a:c+a) .* w, [], "all");
                end
            end
        end%imdilate
        
        function ge = imerode(f, sz)
            % same assumptions as imdilate
            
            w = uint8(ones(sz));
            [rs, cs] = size(f);
            a = floor(0.5*sz);
            ge = f;
            
            for r = (1+a):(rs-a)
                for c = (1+a):(cs-a)
                    ge(r,c) = min(f(r-a:r+a, c-a:c+a) .* w, [], "all");
                end
            end
        end  
    end
    
    %% Public methods
    
    methods(Static, Access = public)
        
        function go = imopen(f, sz)
            % morphological opening
            go = GWE.imdilate(GWE.imerode(f, sz), sz);
        end
        
        function T = graythresh(f)
            %GRAYTHRESH Gonzalez and Wood's threshold algorithm
            %   iterative thresholding algorithm that's an alternative to
            %   the more well-known Otsu method
            
            f = double(f(:));
            T = 0.5 * (min(f) + max(f));
            done = false;
            while ~done
                g = f >= T;
                Tnext = 0.5 * (mean(f(g)) + mean(f(~g)));
                done = abs(T - Tnext) < 0.5;
                T = Tnext;
            end
            T = uint8(T);
        end%graythresh
        
    end
    
end% GWE

