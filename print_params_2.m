function print_params_2(tht,rxn_type,num_unique)
N = size(tht,1);
k0 = 0; lmbda = 0; b_r = 0; b_o = 0; w_ads =0; dG = 0;
R_film = 0;
if rxn_type == "ICET_asymmetric" 
    Numtemps = (N-4);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3));
    b_o = lmbda+abs(tht(4));
    w_ads = tht((5:end));
    params_name = {'k0','lambda','beta_red','beta_ox','w_ads'};
elseif rxn_type == "ICET_asymmetric_film" 
    Numtemps = (N-4-num_unique);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3));
    b_o = lmbda+abs(tht(4));
    w_ads = tht(5:(4+Numtemps));
    R_film = abs(tht(5+Numtemps:end));
    params_name = {'k0','lambda','beta_red','beta_ox','w_ads','R_film'};
elseif rxn_type == "ICET_asymmetric_film2" 
    Numtemps = (N-4-num_unique);
    k0    = (tht(1));
    b_r = lmbda+abs(tht(2));
    b_o = lmbda+abs(tht(3));
    w_ads = tht(4);
    R_film = abs(tht(5:end));
    params_name = {'k0','beta_red','beta_ox','w_ads','R_film'};
elseif rxn_type == "uniformly_valid"
    Numtemps = (N-4);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = abs(tht(3));
    b_o = abs(tht(4));
    w_ads = tht((5:end));
    params_name = {'k0','lambda','beta_red','beta_ox','w_ads'};
elseif rxn_type == "uniformly_valid_i0_fix"
    Numtemps = (N-4);
    lmbda = abs(tht(1));
    b_r = abs(tht(2));
    b_o = abs(tht(3));
    w_ads = tht((4:end));
    params_name = {'lambda','beta_red','beta_ox','w_ads'};
elseif rxn_type == "uniformly_valid_i0_fix_wfix"
    Numtemps = (N-4);
    lmbda = abs(tht(1));
    b_r = abs(tht(2));
    b_o = abs(tht(3));
    w_ads = tht(4);
    params_name = {'lambda','beta_red','beta_ox','w_ads'};
elseif rxn_type == "ICET_symmetric"
    Numtemps = (N-3);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3));
    w_ads = (tht((4:end)));
    params_name = {'k0','lambda','beta_red','w_ads'};
elseif rxn_type == "ICET_symmetric_film"
    Numtemps = (N-3-num_unique);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    b_r = lmbda+abs(tht(3));
    w_ads = tht(4:(3+Numtemps));
    R_film = abs(tht((4+Numtemps):end));
    params_name = {'k0','lambda','beta_red/ox','w_ads','R_film'};
elseif rxn_type == "ICET_symmetric_film2"
    Numtemps = (N-3-num_unique);
    k0    = (tht(1));
    b_r = lmbda+abs(tht(2));
    w_ads = tht(3);
    R_film = abs(tht(4:end));
    params_name = {'k0','beta_red/ox','w_ads','R_film'};
elseif rxn_type == "ECIT_i0_fix"
    Numtemps  = 4;
    lmbda = abs(tht(1));
    w_ads = tht(2:end);
    params_name = {'lambda','w_ads'};
elseif rxn_type == "ECIT_i0_fix_wfix"
    Numtemps  = 4;
    lmbda = abs(tht(1));
    w_ads = tht(2);
    params_name = {'lambda','w_ads'};
elseif rxn_type == "ECIT_normal2"
    Numtemps  = 4;
    k0 = (tht(1));
    lmbda = abs(tht(2));
    dG = abs(tht(3));
    w_ads = tht(4);
    params_name = {'k0','lambda','dG_IT','w_ads'};
elseif rxn_type == "ECIT_film" 
    Numtemps = (N-4-num_unique);
    k0    = (tht(1));
    lmbda = abs(tht(2));
    dG = abs(tht(3));
    w_ads = tht(4);
    R_film = abs(tht(5:end));
    params_name = {'k0','lambda','dG_IT','w_ads','R_film'};
elseif rxn_type == "ECIT_normal"
    Numtemps  = (N-3);
    k0 = (tht(1));
    lmbda = abs(tht(2));
    dG =  abs(tht(3));
    w_ads = tht(4:end);
    params_name = {'k0','lambda','dG_IT','w_ads'};
elseif rxn_type == "uniformly_valid2"
    Numtemps = 5;
    k0    = tht(1);
    lmbda = abs(tht(2));
    b_r = abs(tht(3));
    b_o = abs(tht(4));
    w_ads = tht(5);
    params_name = {'k0','lambda','beta_red','beta_ox','w_ads'};
elseif rxn_type == "BV_film"
    k0    = tht(1:num_unique);
    alpha = cos(abs(tht(num_unique+1)));
    R_film = abs(tht(2+num_unique:end));
    params_name = {'i0','alpha','R_film'};
elseif rxn_type == "Linear"
    k0    = tht(1:num_unique);
    params_name = {'i0'};
elseif rxn_type == "BV"
    k0    = tht(1);
    alpha = tht(2);
    w_ads = tht(3);
    params_name = {'i0','alpha','w_ads'};
elseif rxn_type == "Linear_ICET"
    k0    = tht(1);
    alpha_wads = 1/(1+abs(tht(2))^2);
%    alpha_wads = tht(2);
    dG_icet = abs(tht(3));
    params_name = {'k0','w_ads','dG_icet'};
end
formatSpace = '%s : \n %7.2f \n';
formatSp2 = '\n %7.2f \n';
if rxn_type == "Linear_ICET"
        fprintf(formatSpace,params_name{1},k0);
        fprintf(formatSpace,params_name{2},alpha_wads);
        fprintf(formatSpace,params_name{3},dG_icet);
elseif rxn_type == "Linear"
        fprintf(formatSpace,params_name{1},k0(1));
        fprintf(formatSp2,k0(2:end));   
        
elseif rxn_type == "BV"
        fprintf(formatSpace,params_name{1},k0);
        fprintf(formatSpace,params_name{2},alpha);
        fprintf(formatSpace,params_name{3},w_ads);
elseif rxn_type == "ECIT_film"
        fprintf(formatSpace,params_name{1},k0);
        fprintf(formatSpace,params_name{2},lmbda);
        fprintf(formatSpace,params_name{3},dG);
        fprintf(formatSpace,params_name{4},w_ads); 
        fprintf(formatSpace,params_name{5},R_film(1));  
        fprintf(formatSp2,R_film(2:end));
elseif rxn_type == "BV_film"
        fprintf(formatSpace,params_name{1},k0(1));
        fprintf(formatSp2,k0(2:end));   
        fprintf(formatSpace,params_name{2},alpha);
        fprintf(formatSpace,params_name{3},R_film(1));
        fprintf(formatSp2,R_film(2:end));
elseif rxn_type == "ECIT_i0_fix"
         fprintf(formatSpace,params_name{1},lmbda);
elseif rxn_type == "ECIT_i0_fix_wfix"
         fprintf(formatSpace,params_name{1},lmbda);
elseif rxn_type == "uniformly_valid_i0_fix"
        fprintf(formatSpace,params_name{1},lmbda);
        fprintf(formatSpace,params_name{2},b_r);
        fprintf(formatSpace,params_name{3},b_o);
elseif rxn_type == "uniformly_valid_i0_fix_wfix"
        fprintf(formatSpace,params_name{1},lmbda);
        fprintf(formatSpace,params_name{2},b_r);
        fprintf(formatSpace,params_name{3},b_o);
elseif size(params_name,2) == 6
     fprintf(formatSpace,params_name{1},k0);
        fprintf(formatSpace,params_name{2},lmbda);
        fprintf(formatSpace,params_name{3},b_r);
        fprintf(formatSpace,params_name{4},b_o);
        fprintf(formatSpace,params_name{6},R_film(1));
        fprintf(formatSp2,R_film(2:end));   
elseif size(params_name,2) == 5 && strcmp(params_name{5},'R_film') == 0
        fprintf(formatSpace,params_name{1},k0);
        fprintf(formatSpace,params_name{2},lmbda);
        fprintf(formatSpace,params_name{3},b_r);
        fprintf(formatSpace,params_name{4},b_o);
elseif  strcmp(params_name{4},'R_film') == 1
        fprintf(formatSpace,params_name{1},k0);
        fprintf(formatSpace,params_name{2},b_r);
        fprintf(formatSpace,params_name{4},R_film(1));
        fprintf(formatSp2,R_film(2:end));
elseif  size(params_name,2) == 5 && strcmp(params_name{5},'R_film') == 1
        fprintf(formatSpace,params_name{1},k0);
        fprintf(formatSpace,params_name{2},b_r);
        fprintf(formatSpace,params_name{3},b_o);
        fprintf(formatSpace,params_name{5},R_film(1));
        fprintf(formatSp2,R_film(2:end));
elseif strcmp(params_name{3},'dG_IT') == 1
        fprintf(formatSpace,params_name{1},k0);
        fprintf(formatSpace,params_name{2},lmbda);
        fprintf(formatSpace,params_name{3},dG);
else
     fprintf(formatSpace,params_name{1},k0);
     fprintf(formatSpace,params_name{2},lmbda);
     fprintf(formatSpace,params_name{3},b_r);
end
if rxn_type ~= "BV_film" && rxn_type ~= "Linear" && rxn_type ~= "Linear_ICET"
        if size(w_ads,1) > 1
        fprintf(formatSpace,'w_ads',w_ads(1));
        fprintf(formatSp2,w_ads(2:end))
        else
        fprintf(formatSpace,'w_ads',w_ads);
        end
end
end
