function I_ex = ICET_exchange(k0,lambda,b_r,b_o,w_ads,cc,c_fct,activity,T)  %put argument in the function
        lmbdr = lambda*298/T;
        alpha = b_r/(b_o+b_r);
        w_ads = w_ads*298/T;
        b_r = b_r*298/T; b_o = b_o*298/T;
        c_plus = activity*exp(-w_ads)/(1+activity*exp(-w_ads));  %new
        etaf = log(c_plus/cc);
         M_inf = (0.5*(erf(sqrt(alpha*(1-alpha)*lmbdr))+1)+exp(-alpha*(1-alpha)*lmbdr)/sqrt(4*alpha*(1-alpha)*pi*(lmbdr+eps)));
           M_inf = 1/M_inf;
         M_0 = pi*sqrt(alpha*(1-alpha))/sin(pi*alpha);
      % M = (M_0-M_inf)/((1+exp(-(etaf)/(b_o+lmbdr)))*(1+exp((etaf)/(b_r+lmbdr))))+M_inf;
       M = (M_0-M_inf)/((1+exp(-(etaf+b_o)/(lmbdr+eps)))*(1+exp((etaf-b_r)/(lmbdr+eps))))+M_inf;
       dG_icet = alpha*(1-alpha)*(b_r+b_o+lmbdr);
       if lmbdr ~= 0
       I_ex = k0*exp(-dG_icet)*c_plus^(1-alpha)*cc^(alpha)*(1-cc)^c_fct*M/sqrt(pi*alpha*(1-alpha)*lmbdr);
       else
           I_ex = k0*exp(-dG_icet)*c_plus^(1-alpha)*cc^(alpha)*(1-cc)^c_fct;
       end
       I_ex = abs(I_ex);
end
