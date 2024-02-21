function I_th = ECIT_film(etaf,k0,lambda,dG,w_ads,R_film,cc,c_fct,activity,T)  %put argument in the function
        lmbdr = lambda*293/T;
        w_ads = w_ads*293/T;
        dG = dG*293/T;
        c_plus = activity*exp(-w_ads)/(1+activity*exp(-w_ads));  %new
        etaf = etaf - log(cc/c_plus);
        i0 = k0*exp(-dG)*(1-cc)^c_fct;
       opts = optimset('display', 'none');
       I_th = fzero(@(x) get_loss(x,etaf,i0,cc,c_plus,lmbdr,R_film),(c_plus/(1+exp(etaf))-cc/(1+exp(-(etaf))))*(1-erf((lmbdr-sqrt(1+sqrt(lmbdr)+etaf^2))...
           /(2*sqrt(lmbdr)))),opts);
       I_th = i0*I_th;
end

function err = get_loss(x,etaf,i0,cc,c_plus,lmbdr,R_film)
err = x - (c_plus/(1+exp(etaf-x*i0*R_film))-cc/(1+exp(-(etaf-x*i0*R_film))))*(1-erf((lmbdr-sqrt(1+sqrt(lmbdr)+(etaf-x*i0*R_film)^2))/(2*sqrt(lmbdr))));
end

function M = get_M(etaf,b_o,b_r,lmbdr,alpha)
       M_inf = (0.5*(erf(sqrt(alpha*(1-alpha)*lmbdr))+1)+exp(-alpha*(1-alpha)*lmbdr)/sqrt(4*alpha*(1-alpha)*pi*lmbdr));
       M_inf = 1/M_inf;
       M_0 = pi*sqrt(alpha*(1-alpha))/sin(pi*alpha);
      % M = (M_0-M_inf)/((1+exp(-(etaf)/(b_o+lmbdr)))*(1+exp((etaf)/(b_r+lmbdr))))+M_inf;
       M = (M_0-M_inf)/((1+exp(-(etaf+b_o)/(lmbdr)))*(1+exp((etaf-b_r)/(lmbdr))))+M_inf;
end