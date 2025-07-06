function I_th = ECIT_i0_fix(etaf,i0,lambda,w_ads,cc,c_fct,activity,T)  %put argument in the function
        lmbdr = lambda*293/T;
        w_ads = w_ads*293/T;
        c_plus = activity*exp(-w_ads)/(1+activity*exp(-w_ads));  %new
        etaf = etaf - log(cc/c_plus);
        i0_norm = (c_plus*cc/(c_plus+cc))*(1-erf((lmbdr-sqrt(1+sqrt(lmbdr)+(-log(cc/c_plus))^2))/(2*sqrt(lmbdr))));
        I_th = (-i0/i0_norm)*(c_plus/(1+exp(etaf))-cc/(1+exp(-etaf)))*(1-erf((lmbdr-sqrt(1+sqrt(lmbdr)+etaf^2))/(2*sqrt(lmbdr))));
end 