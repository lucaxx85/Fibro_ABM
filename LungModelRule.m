classdef LungModelRule < handle
    properties (Abstract)
        Model;
    end

    methods (Abstract, Static)
        apply_rule(this);
        finalize(this);
    end
    
    methods
        function nbhd = neighborhood_count(this, mat, state)
            % Count the number of cells that match a given state in each
            % cell's Moore neighborhood.
            [m, n] = size(mat);
            
            up = [2:m 1];
            down = [m 1:m-1];
            left = [2:n 1];
            right = [n 1:n-1];
            
            nbhd = (mat(up, :) == state) + ...
                   (mat(up, left) == state) + ...
                   (mat(up, right) == state) + ...
                   (mat(down, :) == state) + ...
                   (mat(down, left) == state) + ...
                   (mat(down, right) == state) + ...
                   (mat(:, left) == state) + ...
                   (mat(:, right) == state);
        end
        
        function diffused = diffuse(mat)
            % Given a matrix containing a fluid intensity (0 to 1),
            % diffuse the fluid across the matrix evenly.
            [m, n] = size(mat);
            
            up = [2:m 1];
            down = [m 1:m-1];
            left = [2:n 1];
            right = [n 1:n-1];
            
            diffused = 1/9 * (mat(up, left) + ...
                       mat(up, :) + ...
                       mat(up, right) + ...
                       mat(:, left) + ...
                       mat(:, :) + ...
                       mat(:, right) + ...
                       mat(down, left) + ...
                       mat(down, :) + ...
                       mat(down, right));
        end
    end
end