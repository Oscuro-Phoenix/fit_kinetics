function [c_cost] = cost_func(tht,x,y,size_vec,temp_vec,conc_vec,rxn_type)
%   Model Predictions
    f_x = predict(x,tht,size_vec,temp_vec,conc_vec,rxn_type);
%   Calculate the cost function
    c_cost = norm(abs(f_x)-abs(y)).^2;
    %c_cost =  norm(log(f_x)-log(y)).^2 + norm((f_x-y)./(tht(1)*(1-conc_vec))).^2 ;
    %c_cost =  norm((abs(x)<tht(2)).*(log(f_x)-log(y))+(abs(x)>tht(2)).*((f_x-y)./(tht(1).*(1-conc_vec)))).^2;
end