function p  = get_params_2(tht,rxn_type,num_unique)
N = size(tht,1);
if rxn_type == "ICET_asymmetric" 
    Numtemps = (N-4);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3));
    b_o = lmbda+abs(tht(4));
    w_ads = tht((5:end));
    p = [k0;lmbda;b_r;b_o;w_ads];
elseif rxn_type == "ICET_asymmetric_film" 
    Numtemps = (N-4-num_unique);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3));
    b_o = lmbda+abs(tht(4));
    w_ads = tht(5:(4+Numtemps));
    R_film = abs(tht(5+Numtemps:end));
    p = [k0;lmbda;b_r;b_o;w_ads;R_film];
elseif rxn_type == "ICET_asymmetric_film2" 
    Numtemps = (N-4-num_unique);
    k0 = tht(1);
    b_r = abs(tht(2));
    b_o = abs(tht(3));
    w_ads = tht(4);
    R_film = abs(tht(5:end));
    p = [k0;b_r;b_o;w_ads;R_film];
elseif rxn_type == "ECIT_film" 
    Numtemps = (N-4-num_unique);
    k0 = (tht(1));
    lmbda = abs(tht(2));
    dG = abs(tht(3));
    w_ads = tht(4);
    R_film = abs(tht(5:end));
    p = [k0;lmbda;dG;w_ads;R_film];
elseif rxn_type == "uniformly_valid"
    Numtemps = (N-4);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = abs(tht(3));
    b_o = abs(tht(4));
    w_ads = tht((5:end));
    p = [k0;lmbda;b_r;b_o;w_ads];
elseif rxn_type == "uniformly_valid_i0_fix"
    lmbda = abs(tht(1));
    b_r = abs(tht(2));
    b_o = abs(tht(3));
    w_ads = tht((4:end));
    p = [lmbda;b_r;b_o;w_ads];
elseif rxn_type == "uniformly_valid_i0_fix_wfix"
    lmbda = abs(tht(1));
    b_r = abs(tht(2));
    b_o = abs(tht(3));
    w_ads = tht(4);
    p = [lmbda;b_r;b_o;w_ads];
elseif rxn_type == "ICET_symmetric"
    Numtemps = (N-3);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3));
    w_ads = (tht((4:end)));
    p = [k0;lmbda;b_r;w_ads];
elseif rxn_type == "ICET_symmetric_film"
    Numtemps = (N-3-num_unique);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3));
    w_ads = tht(4:(3+Numtemps));
    R_film = abs(tht((4+Numtemps):end));
    p = [k0;lmbda;b_r;w_ads;R_film];
elseif rxn_type == "ICET_symmetric_film2"
    Numtemps = (N-3-num_unique);
    k0    = (tht(1));
    b_r = 10+abs(tht(2));
    w_ads = tht(3);
    R_film = abs(tht(4:end));
    p = [k0;b_r;w_ads;R_film];
elseif rxn_type == "BV_film"
    k0    = (tht(1:num_unique));
    alpha = 1/(1+abs(tht(num_unique+1))^2);
    R_film = abs(tht(num_unique+2:end));
    p = [k0;alpha;R_film];
elseif rxn_type == "ECIT_normal2"
    Numtemps  = 4;
    k0 = (tht(1));
    lmbda = abs(tht(2));
    dG = abs(tht(3));
    w_ads = tht(4);
    p = [k0;lmbda;dG;w_ads];
elseif rxn_type == "ECIT_normal"
    Numtemps  = (N-3);
    k0 = (tht(1));
    lmbda = abs(tht(2));
    dG = abs(tht(3));
    w_ads = tht(4:end);
    p = [k0;lmbda;dG;w_ads];
elseif rxn_type == "ECIT_i0_fix"
    lmbda = abs(tht(1));
    w_ads = tht(2:end);
    p = [lmbda;w_ads];
elseif rxn_type == "ECIT_i0_fix_wfix"
    lmbda = abs(tht(1));
    w_ads = tht(2);
    p = [lmbda;w_ads];
elseif rxn_type == "uniformly_valid2"
    Numtemps = 5;
    k0    = tht(1);
    lmbda = abs(tht(2));
    b_r = abs(tht(3));
    b_o = abs(tht(4));
    w_ads = tht(5);
    p = [k0;lmbda;b_r;b_o;w_ads];
elseif rxn_type == "Linear"
    i0    = tht(1:num_unique);
    p = i0;
elseif rxn_type == "Linear_ICET"
    k0    = tht(1);
    alpha_wads = 1/(1+abs(tht(2))^2);
    %alpha_wads = tht(2);
    dG_icet = abs(tht(3));
    p = [k0;alpha_wads;dG_icet];
elseif rxn_type == "BV"
    k0    = tht(1);
    alpha = tht(2);
    w_ads = tht(3);
    p = [k0;alpha;w_ads];
end
end