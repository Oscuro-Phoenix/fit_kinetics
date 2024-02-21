function I_th = BV_film(etaf,k0,alpha,R_film,cc,c_fct,activity,T)  %put argument in the function
       i0 = k0*((1.0-cc).^(alpha))*(cc^alpha);
       opts = optimset('display', 'none');
       I_th = fzero(@(x) get_loss(x,etaf,i0,alpha,R_film),(exp(-alpha*etaf)-exp((1-alpha)*etaf)),opts);
       I_th = i0*I_th;
end

function err = get_loss(x,etaf,i0,alpha,R_film)
err = x - (-exp((1-alpha)*(etaf-x*i0*R_film))+exp(-alpha*(etaf-x*i0*R_film)));
end

function M = get_M(etaf,b_o,b_r,lmbdr,alpha)
       M_inf = (0.5*(erf(sqrt(alpha*(1-alpha)*lmbdr))+1)+exp(-alpha*(1-alpha)*lmbdr)/sqrt(4*alpha*(1-alpha)*pi*lmbdr));
       M_inf = 1/M_inf;
       M_0 = pi*sqrt(alpha*(1-alpha))/sin(pi*alpha);
      % M = (M_0-M_inf)/((1+exp(-(etaf)/(b_o+lmbdr)))*(1+exp((etaf)/(b_r+lmbdr))))+M_inf;
       M = (M_0-M_inf)/((1+exp(-(etaf+b_o)/(lmbdr)))*(1+exp((etaf-b_r)/(lmbdr))))+M_inf;
end