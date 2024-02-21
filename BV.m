function I_th = BV(eta,k0,alpha,w_ads,activity,cc,T)  %put argument in the function
        w_ads = w_ads*293/T;
        c_plus = activity*exp(-w_ads)/(1+exp(-w_ads)); 
        I_th = k0*(c_plus^(1-alpha)*cc^alpha)*(exp(-alpha*eta)-exp((1-alpha)*eta));   %new
end