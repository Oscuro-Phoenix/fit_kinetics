function p  = get_params(tht,rxn_type)
N = size(tht,1);
if rxn_type == "ICET_asymmetric" 
    Numtemps = (N-2)/3;
    k0    = abs(tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3:(Numtemps+2)));
    b_o = lmbda+abs(tht((Numtemps+3):(2*Numtemps+2)));
    w_ads = tht((2*Numtemps+3:end));
    p = [k0;lmbda;b_r;b_o;w_ads];
elseif rxn_type == "ICET_asymmetric_film" 
    Numtemps = (N-3)/3;
    k0    = abs(tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3:(Numtemps+2)));
    b_o = lmbda+abs(tht((Numtemps+3):(2*Numtemps+2)));
    w_ads = tht((2*Numtemps+3:end-1));
    R_film = abs(tht(end));
    p = [k0;lmbda;b_r;b_o;w_ads;R_film];
elseif rxn_type == "uniformly_valid"
    Numtemps = (N-2)/3;
    k0    = abs(tht(1));
    lmbda = abs(tht(2));
    b_r = abs(tht(3:(Numtemps+2)));
    b_o = abs(tht((Numtemps+3):(2*Numtemps+2)));
    w_ads = tht((2*Numtemps+3:end));
    p = [k0;lmbda;b_r;b_o;w_ads];
elseif rxn_type == "ICET_symmetric"
    Numtemps = (N-2)/2;
    k0    = abs(tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3:(Numtemps+2)));
    w_ads = (tht((Numtemps+3:end)));
    p = [k0;lmbda;b_r;w_ads];
elseif rxn_type == "ICET_symmetric_film"
    Numtemps = (N-3)/2;
    k0    = abs(tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3:(Numtemps+2)));
    w_ads = (tht((Numtemps+3:end-1)));
    R_film = abs(tht(end));
    p = [k0;lmbda;b_r;w_ads;R_film];
elseif rxn_type == "ECIT_normal2"
    Numtemps  = (N-3);
    k0 = abs(tht(1));
    lmbda = abs(tht(2));
    dG = abs(tht(3:(Numtemps+2)));
    w_ads = tht(Numtemps+3);
    p = [k0;lmbda;dG;w_ads];
elseif rxn_type == "ECIT_normal"
    Numtemps  = (N-2)/2;
    k0 = abs(tht(1));
    lmbda = abs(tht(2));
    dG = lmbda/4 + abs(tht(3:(Numtemps+2)));
    w_ads = tht(Numtemps+3:end);
    p = [k0;lmbda;dG;w_ads];
elseif rxn_type == "uniformly_valid2"
    Numtemps = (N-3)/2;
    k0    = tht(1);
    lmbda = abs(tht(2));
    b_r = abs(tht((3):(Numtemps+2)));
    b_o = abs(tht((Numtemps+3):(2*Numtemps+2)));
    w_ads = tht(2*Numtemps+3);
    p = [k0;lmbda;b_r;b_o;w_ads];
end
end