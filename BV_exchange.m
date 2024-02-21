function I_ex = BV_exchange(k0,alpha,cc,c_fct,activity,T)  %put argument in the function
       i0 = k0*((1.0-cc).^alpha)*(cc^alpha);
       I_ex = abs(i0);
end