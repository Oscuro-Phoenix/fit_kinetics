function J = get_gradient(tht,x,size_vec,temp_vec,conc_vec,rxn_type,activity)
    func = @(t) predict_2(x,t,size_vec,temp_vec,conc_vec,rxn_type,activity);
    N = size(tht,1);
    M = size(x,1);
    J = zeros(M,N);
    h = 1e-3;
    [fval,~] = func(tht);
    for i=1:N
        t = tht;
        t(i) = t(i) + h;
    %   Model Predictions
        [fup,~] = func(t);
        J(:,i) = (fup-fval)/h; 
    end
%   Calculate the cost function
%     if rxn_type == "ICET_symmetric" || rxn_type == "ICET_asymmetric" || ...
%             rxn_type == "ICET_symmetric_film" || rxn_type == "ICET_asymmetric_film" ...
%             || rxn_type == "ICET_symmetric_film2" || rxn_type == "ICET_asymmetric_film2" ...
%             || rxn_type == "ECIT_normal" || rxn_type=="ECIT_normal2" || rxn_type == "BV_film" || ...
%             rxn_type == "ECIT_film"
%         slope_vals = get_slopes(x,y,size_vec);
%         pure_cost = norm(f_x-y).^2;
%         c_cost = pure_cost + 0.1*norm(I_ex_pred-slope_vals).^2;
%     else
%         c_cost = norm(f_x-y).^2;
%     end
end