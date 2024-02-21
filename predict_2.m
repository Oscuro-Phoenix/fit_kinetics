function [f_mdl,I_ex_pred] = predict_2(x,tht,size_vec,temp_vec,conc_vec,rxn_type,activity)
    Ndata = length(x);
    f_mdl = zeros(Ndata,1);
    I_ex_pred = zeros(size(temp_vec));
    num_unique = size(unique(temp_vec),2);
    p = get_params_2(tht,rxn_type,num_unique);
    %activity = 1.6;        %1.6 for LiClO4, 1.9 for LiPF6 in EMC/EC (0.02,0.11,7.6 for clyte diff) 1 for EMC solvent                       %input5: activity in liquid
    c_fct = 1.0;
    N = size(size_vec,1);
    k = 1;
    if rxn_type == "ICET_asymmetric" 
    for j=1:N
        I_ex_pred(j) = ICET_exchange(p(1),p(2),p(3),p(4),p(4+j),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            if j > 1
                %Args to ICET: x,k0,lambda,br,bo,wads
            f_mdl(k) = ICET(x(i+sum(size_vec(1:j-1))),p(1),p(2),p(3),p(4),p(4+j),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ICET(x(i),p(1),p(2),p(3),p(4),p(4+j),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
    elseif rxn_type == "ICET_asymmetric_film" 
    for j=1:N
        I_ex_pred(j) = ICET_exchange(p(1),p(2),p(3),p(4),p(4+j),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            nval =0;
            if j > (N-num_unique+1) 
                nval = j-(N-num_unique+1);
            end
            if j > 1
                %Args to ICET: x,k0,lambda,br,bo,wads,R_film
            f_mdl(k) = ICET_film(x(i+sum(size_vec(1:j-1))),p(1),p(2),p(3),p(4),p(4+j),p(5+N+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ICET_film(x(i),p(1),p(2),p(3),p(4),p(4+j),p(5+N+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
elseif rxn_type == "ICET_asymmetric_film2" 
    for j=1:N
        I_ex_pred(j) = ICET_exchange(p(1),0,p(2),p(3),p(4),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            nval =0;
            if j > (N-num_unique+1) 
                nval = j-(N-num_unique+1);
            end
            if j > 1
                %Args to ICET: x,k0,lambda,br,bo,wads,R_film
            f_mdl(k) = ICET_film(x(i+sum(size_vec(1:j-1))),p(1),0,p(2),p(3),p(4),p(5+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ICET_film(x(i),p(1),0,p(2),p(3),p(4),p(5+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
    elseif rxn_type == "uniformly_valid"
    for j=1:N
        for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = uniformly_valid_old(x(i+sum(size_vec(1:j-1))),...
                p(1),p(2),p(3),p(4),p(4+j),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = uniformly_valid_old(x(i),...
                p(1),p(2),p(3),p(4),p(4+j),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
    elseif rxn_type == "ICET_symmetric"
     for j=1:N
        I_ex_pred(j) = ICET_exchange(p(1),p(2),p(3),p(3),p(3+j),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = ICET(x(i+sum(size_vec(1:j-1))),p(1),p(2),p(3),p(3),p(3+j),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ICET(x(i),p(1),p(2),p(3),p(3),p(3+j),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
     end
    elseif rxn_type == "ICET_symmetric_film"
     for j=1:N
         I_ex_pred(j) = ICET_exchange(p(1),p(2),p(3),p(3),p(3+j),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            nval =0;
            if j > (N-num_unique+1) 
                nval = j-(N-num_unique+1);
            end
            if j > 1
            f_mdl(k) = ICET_film(x(i+sum(size_vec(1:j-1))),p(1),p(2),p(3),p(3),p(3+j),p(4+N+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ICET_film(x(i),p(1),p(2),p(3),p(3),p(3+j),p(4+N+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
     end
    elseif rxn_type == "BV_film" %BVfilm needs k0,alpha,w_ads,R_film
     for j=1:N
        nval =0;
        if j > (N-num_unique+1) 
                nval = j-(N-num_unique+1);
        end
         I_ex_pred(j) = BV_exchange(p(1+nval),p(1+num_unique),conc_vec(j),c_fct,activity,temp_vec(j));
         for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = BV_film(x(i+sum(size_vec(1:j-1))),p(1+nval),p(1+num_unique),p(2+num_unique+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = BV_film(x(i),p(1+nval),p(1+num_unique),p(2+num_unique+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
     end
    elseif rxn_type == "ICET_symmetric_film2"
     for j=1:N
         I_ex_pred(j) = ICET_exchange(p(1),0,p(2),p(2),p(3),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            nval =0;
            if j > (N-num_unique+1) 
                nval = j-(N-num_unique+1);
            end
            if j > 1
            f_mdl(k) = ICET_film(x(i+sum(size_vec(1:j-1))),p(1),0,p(2),p(2),p(3),p(4+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ICET_film(x(i),p(1),0,p(2),p(2),p(3),p(4+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
   elseif rxn_type == "ECIT_normal"
     for j=1:N
      I_ex_pred(j) = ECIT_exchange(p(1),p(2),p(3),p(3+j),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = ECIT(x(i+sum(size_vec(1:j-1))),p(1),p(2),p(3),p(3+j),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ECIT(x(i),p(1),p(2),p(3),p(3+j),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
     end
      elseif rxn_type == "ECIT_normal2"
     for j=1:N
        I_ex_pred(j) = ECIT_exchange(p(1),p(2),p(3),p(4),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = ECIT(x(i+sum(size_vec(1:j-1))),p(1),p(2),p(3),p(4),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ECIT(x(i),p(1),p(2),p(3),p(4),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
     end
 elseif rxn_type == "ECIT_film" 
    for j=1:N
        I_ex_pred(j) = ECIT_exchange(p(1),p(2),p(3),p(4),conc_vec(j),c_fct,activity,temp_vec(j));
        for i = 1:size_vec(j)
            nval =0;
            if j > (N-num_unique+1) 
                nval = j-(N-num_unique+1);
            end
            if j > 1
                %Args to ECIT film: x,k0,dG,lambda,wads,R_film
            f_mdl(k) = ECIT_film(x(i+sum(size_vec(1:j-1))),p(1),p(2),p(3),p(4),p(5+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = ECIT_film(x(i),p(1),p(2),p(3),p(4),p(5+nval),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
    elseif rxn_type == "uniformly_valid2"     
    for j=1:N
        for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = uniformly_valid_old(x(i+sum(size_vec(1:j-1))),...
                p(1),p(2),p(3),p(4),p(5),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            else
            f_mdl(k) = uniformly_valid_old(x(i),...
                p(1),p(2),p(3),p(4),p(5),conc_vec(j),c_fct,activity,temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
    elseif rxn_type == "Linear"
    for j=1:N
        nval =0;
        if j > (N-num_unique+1) 
                nval = j-(N-num_unique+1);
        end
        for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = Linear(x(i+sum(size_vec(1:j-1))),...
                p(1+nval),conc_vec(j));   %new arguments
            else
            f_mdl(k) = Linear(x(i),...
                p(1+nval),conc_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
    elseif rxn_type == "Linear_ICET"
    for j=1:N
        for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = Linear_ICET(x(i+sum(size_vec(1:j-1))),...
                p(1),p(2),p(3),conc_vec(j),temp_vec(j));   %new arguments
            else
            f_mdl(k) = Linear_ICET(x(i),...
                p(1),p(2),p(3),conc_vec(j),temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
    elseif rxn_type == "BV"
    for j=1:N
        for i = 1:size_vec(j)
            if j > 1
            f_mdl(k) = BV(x(i+sum(size_vec(1:j-1))),...
                p(1),p(2),p(3),activity,conc_vec(j),temp_vec(j));   %new arguments
            else
            f_mdl(k) = BV(x(i),...
                p(1),p(2),p(3),activity,conc_vec(j),temp_vec(j));   %new arguments
            end
            k=k+1;
        end
    end
    end    
end
