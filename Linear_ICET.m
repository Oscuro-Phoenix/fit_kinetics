function I_th = Linear_ICET(etaf,k0,alpha,dG_icet,cc,T)  %put argument in the function
        %I_th = i0*sqrt(cc*(1-cc))*etaf;
        dG_icet = dG_icet*298/T;
        I_th = 1e8*k0*sqrt(T/298)*exp(-dG_icet)*cc^(alpha)*(1-cc)*etaf;
%         w_ads = w_ads*298/T;
%         a = 1.6;
%         c_plus = a*exp(-w_ads)/(1+a*exp(-w_ads));
%         I_th = 1e8*k0*exp^(-dG_icet*298/T)*(1-cc)*0.5*cc*c_plus/(c_plus+cc)*etaf;
end
