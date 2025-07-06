function c_cost = cost_func_2(tht,x,y,size_vec,temp_vec,conc_vec,rxn_type,activity,slopes)
%   Model Predictions
    [f_x,I_ex_pred] = predict_2(x,tht,size_vec,temp_vec,conc_vec,rxn_type,activity,slopes);
%   Calculate the cost function
    if rxn_type == "ICET_symmetric" || rxn_type == "ICET_asymmetric" || ...
            rxn_type == "ICET_symmetric_film" || rxn_type == "ICET_asymmetric_film" ...
            || rxn_type == "ICET_symmetric_film2" || rxn_type == "ICET_asymmetric_film2" ...
            || rxn_type == "ECIT_normal" || rxn_type=="ECIT_normal2" || rxn_type == "BV_film" || ...
            rxn_type == "ECIT_film"
        slope_vals = get_slopes(x,y,size_vec);
        pure_cost = norm(f_x-y).^2;
        c_cost = pure_cost + 0.1*norm(I_ex_pred-slope_vals).^2;
    else
        c_cost = norm(f_x-y).^2;
    end
end