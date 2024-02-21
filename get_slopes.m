function slope_vals = get_slopes(x,y,size_vec)
N = size(size_vec,1);
slope_vals = zeros(N,1);
for j=1:N
           if j == 1
               X_vals = x(1:size_vec(1));
               Y_vals = y(1:size_vec(1));
           else
               X_vals = x(1+sum(size_vec(1:j-1)):sum(size_vec(1:j)));
               Y_vals = y(1+sum(size_vec(1:j-1)):sum(size_vec(1:j)));
           end
           c = polyfit(X_vals, Y_vals,1);
           slope_vals(j) = c(1); %new arguments
end
end
