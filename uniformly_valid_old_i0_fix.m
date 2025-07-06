function I_th = uniformly_valid_old_i0_fix(etaf,i0,lambda,b_red,b_ox,w_ads,cc,c_fct,activity,T)  %put argument in the function
     lmbdr = lambda*298/T;
        b_red = b_red*298/T;
        b_ox = b_ox*298/T;
        w_ads = w_ads*298/T;
        c_plus = activity*exp(-w_ads)/(1+activity*exp(-w_ads));  %new
        etaf = etaf - log(cc/c_plus); 
        alpha = b_red/(b_red+b_ox);
        M_inf = 1/(0.5*(erf(sqrt(alpha*(1-alpha)*lmbdr))+1)+exp(-alpha*(1-alpha)*lmbdr)/sqrt(4*alpha*(1-alpha)*pi*lmbdr));
        M_0 = pi*sqrt(alpha*(1-alpha))/sin(pi*alpha);
        %M = (M_0-M_inf)/(1+exp(-(etaf)/(b_ox+lmbdr)))/(1+exp((etaf)/(b_red+lmbdr)))+M_inf;
        M = (M_0-M_inf)/(1+exp(-(etaf+b_ox)/sqrt(lmbdr)))/(1+exp((etaf-b_red)/sqrt(lmbdr)))+M_inf;

        b_o = b_ox;
        b_r = b_red;
        if etaf > -b_ox
            kred = exp(-alpha*((1-alpha)*lmbdr+b_o))*exp(-alpha*etaf)/sqrt(4*alpha*(1-alpha)*pi*lmbdr);
        else
            kred = 0.5*(erf(sqrt(alpha*(1-alpha)*lmbdr))-erf((2*(1-alpha)*lmbdr+b_o+etaf)/(2*sqrt((1-alpha)*lmbdr/alpha))))+exp(-alpha*(1-alpha)*lmbdr)/sqrt(4*alpha*(1-alpha)*pi*lmbdr);
        end
%       kox calc
        if etaf < b_red
             kox = exp(-(1-alpha)*(alpha*lmbdr+b_r))*exp((1-alpha)*etaf)/sqrt(4*pi*alpha*(1-alpha)*lmbdr);
        else
             kox = 0.5*(erf(sqrt(alpha*(1-alpha)*lmbdr))-erf((2*alpha*lmbdr+b_r-etaf)/(2*sqrt(alpha*lmbdr/(1-alpha)))))+exp(-alpha*(1-alpha)*lmbdr)/sqrt(4*alpha*(1-alpha)*pi*lmbdr);
        end
%       Current
        %I_th = (k0*M)*((1.0-cc).^c_fct)*(c_plus*kred-cc*kox);   %new
        i0_norm = (exp(-alpha*((1-alpha)*lmbdr+b_o))*(cc^alpha)*(c_plus^(1-alpha)))/sqrt(4*alpha*(1-alpha)*pi*lmbdr);
        I_th = (-i0/i0_norm)*(c_plus*kred-cc*kox);
end
