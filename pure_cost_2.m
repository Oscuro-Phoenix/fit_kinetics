function [pure_cost,chi_squared] = pure_cost_2(tht,x,y,size_vec,temp_vec,conc_vec,rxn_type,activity,slopes)
%   Model Predictions
    [f_x,~] = predict_2(x,tht,size_vec,temp_vec,conc_vec,rxn_type,activity,slopes);
%   Calculate the cost function
    pure_cost = norm((f_x-y)).^2;
    chi_squared = sum(((f_x-y).^2)./abs(y),'all');
end