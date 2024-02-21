function c_cost = get_rmse_with_dof(tht,x,y,size_vec,temp_vec,conc_vec,rxn_type,activity)
%   Model Predictions
    [f_x,~] = predict_2(x,tht,size_vec,temp_vec,conc_vec,rxn_type,activity);
    c_cost = norm(f_x-y);
    c_cost = c_cost/sqrt(size(y,1)-size(tht,1));
end