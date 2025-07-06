function slope_vals = get_slopes(x,y,size_vec)
    N = numel(size_vec);
    slope_vals = zeros(N,1);
    
    for j = 1:N
        % Extract segment
        if j == 1
            start_idx = 1;
            end_idx = size_vec(1);
        else
            start_idx = 1 + sum(size_vec(1:j-1));
            end_idx = sum(size_vec(1:j));
        end
        
        X_vals = x(start_idx:end_idx);
        Y_vals = y(start_idx:end_idx);
        
        % Filter values where |x| < 2
        valid_mask = abs(X_vals) < 2;
        X_filtered = X_vals(valid_mask);
        Y_filtered = Y_vals(valid_mask);
        
        % Calculate slope only if enough points remain
        if numel(X_filtered) >= 2
            c = polyfit(X_filtered, Y_filtered, 1);
            slope_vals(j) = c(1);
        else
            c = polyfit(X_vals, Y_vals, 1);
            slope_vals(j) = c(1); % Handle insufficient points
        end
    end
end

% function slope_vals = get_slopes(x,y,size_vec)
% N = size(size_vec,1);
% slope_vals = zeros(N,1);
% for j=1:N
%            if j == 1
%                X_vals = x(1:size_vec(1));
%                Y_vals = y(1:size_vec(1));
%            else
%                X_vals = x(1+sum(size_vec(1:j-1)):sum(size_vec(1:j)));
%                Y_vals = y(1+sum(size_vec(1:j-1)):sum(size_vec(1:j)));
%            end
%            c = polyfit(X_vals, Y_vals,1);
%            slope_vals(j) = c(1); %new arguments
% end
% end
