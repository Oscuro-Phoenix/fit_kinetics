function I_ex = ECIT_exchange(k0,lambda,dG,w_ads,cc,c_fct,activity,T)  %put argument in the function
        lmbdr = lambda*298/T;
        w_ads = w_ads*298/T;
        dG = dG*298/T;
        c_plus = activity*exp(-w_ads)/(1+activity*exp(-w_ads));  %new
        etaf = - log(cc/c_plus);
        I_ex = k0*exp(-dG)*(1-cc)^c_fct*0.5*(cc*c_plus/(cc+c_plus))*(1-erf((lmbdr-sqrt(1+sqrt(lmbdr)+etaf^2))/(2*sqrt(lmbdr))));
end
