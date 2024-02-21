function I_th = uniformly_valid(etaf,k0,lambda,b_red,b_ox,w_ads,cc,c_fct,activity,T)  %put argument in the function
     lmbdr = lambda*298/T;
        b_red = b_red*298/T;
        b_ox = b_ox*298/T;
        w_ads = w_ads*298/T;
        c_plus = activity*exp(-w_ads)/(1+activity*exp(-w_ads));  %new
        etaf = etaf - log(cc/c_plus); 
        alpha = b_red/(b_red+b_ox);
        M_inf = 1/(0.5*(erf(sqrt(alpha*(1-alpha)*lmbdr))+1)+exp(-alpha*(1-alpha)*lmbdr)/sqrt(4*alpha*(1-alpha)*pi*lmbdr));
        M_0 = pi*sqrt(alpha*(1-alpha))/sin(pi*alpha);
       % M = (M_0-M_inf)/(1+exp(-(etaf+b_ox)/(lmbdr)))/(1+exp((etaf-b_red)/(lmbdr)))+M_inf;
        M = (M_0-M_inf)/((1+exp(-(etaf)/(b_ox+lmbdr)))*(1+exp((etaf)/(b_red+lmbdr))))+M_inf;
        C = alpha/(4*(1-alpha)*lmbdr);
        D = 2*lmbdr*(1-alpha)+b_ox;
        A = (1-alpha)/(4*alpha*lmbdr);
        B = 2*alpha*lmbdr + b_red;
        if etaf > b_red
            kred = 0.5*sqrt(C/A)*exp(1/(4*A)-B)*erfc(sqrt(A)*(etaf+1/(2*A)-B));
        elseif etaf > -b_ox
            kred = 0.5*sqrt(C/A)*exp(1/(4*A)-B)*erfc(sqrt(A)*(b_red+1/(2*A)-B))+exp(-alpha*(1-alpha)*lmbdr-alpha*b_ox)*(exp(-alpha*etaf)-exp(-alpha*b_red))...
                /sqrt(4*pi*alpha*(1-alpha)*lmbdr);
        else
            kred =  0.5*sqrt(C/A)*exp(1/(4*A)-B)*erfc(sqrt(A)*(b_red+1/(2*A)-B))+exp(-alpha*(1-alpha)*lmbdr-alpha*b_ox)*(exp(alpha*b_ox)-exp(-alpha*b_red))/sqrt(4*pi*alpha*(1-alpha)*lmbdr) + 0.5*(erf(sqrt(C)*(D-b_ox))-erf(sqrt(C)*(D+etaf)));
        end
%       kox calc
        if etaf < -b_ox
            kox = (1/2)*sqrt(A/C)*exp(-D+1/(4*C))*(1+erf(sqrt(C)*(etaf-1/(2*C)+D)));
        elseif etaf < b_red
             kox = (1/2)*sqrt(A/C)*exp(-D+1/(4*C))*(1+erf(sqrt(C)*(-b_ox-1/(2*C)+D))) ...
                 + exp(-alpha*(1-alpha)*lmbdr-(1-alpha)*b_red)*(exp((1-alpha)*etaf)-exp(-(1-alpha)*b_ox))/sqrt(4*pi*alpha*(1-alpha)*lmbdr);
        else
             kox = (1/2)*sqrt(A/C)*exp(-D+1/(4*C))*(1+erf(sqrt(C)*(-b_ox-1/(2*C)+D))) ...
                 + exp(-alpha*(1-alpha)*lmbdr-(1-alpha)*b_red)*(exp((1-alpha)*b_red)-exp(-(1-alpha)*b_ox))/sqrt(4*pi*alpha*(1-alpha)*lmbdr) ...
                 + (1/2)*(erf(sqrt(A)*(B-b_red))-erf(sqrt(A)*(B-etaf)));
        end
%       Current
        I_th = (k0*M)*((1.0-cc).^c_fct)*(c_plus*kred-cc*kox);   %new
end
