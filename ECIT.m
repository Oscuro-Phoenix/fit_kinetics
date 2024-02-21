function I_th = ECIT(etaf,k0,lambda,dG,w_ads,cc,c_fct,activity,T)  %put argument in the function
        lmbdr = lambda*293/T;
        w_ads = w_ads*293/T;
        dG = dG*293/T;
        c_plus = activity*exp(-w_ads)/(1+activity*exp(-w_ads));  %new
        etaf = etaf - log(cc/c_plus);
        I_th = k0*exp(-dG)*(1-cc)^c_fct*(c_plus/(1+exp(etaf))-cc/(1+exp(-etaf)))*(1-erf((lmbdr-sqrt(1+sqrt(lmbdr)+etaf^2))/(2*sqrt(lmbdr))));
end