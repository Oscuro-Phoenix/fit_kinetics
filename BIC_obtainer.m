clearvars;

rxn_types = {'ECIT_normal2','ICET_asymmetric_film2','BV'};
filenames = {'ECIT','ICET+Film','BV'};

for i=1:size(rxn_types,2)
%% BIC for ECIT
temp_vec = [293,293,293,293,308,318];
conc_vec = [0.5,0.6,0.7,0.8,0.5,0.5];
rxn_type = rxn_types{i};
activity_LiClO4 = 1.6;
activity_LPF = 1.9;

%LiClO4
[tht_fit_LCO,f_fit_LCO,Y_LCO] = get_fit_params({'LCO_Li05.csv','LCO_Li06.csv','LCO_Li07.csv','LCO_Li08.csv','LCO_LiClO4_T_dependent_35C.csv','LCO_LiClO4_T_dependent_45C.csv'}...
    ,['./Fit_Results/LCO_x/',filenames{i},'/fitted_params.csv'],temp_vec, conc_vec, rxn_type, activity_LiClO4);

temp_vec = [293,293,293,293];
conc_vec = [0.5,0.6,0.7,0.8];
[tht_fit_NMC,f_fit_NMC,Y_NMC] = get_fit_params({'NMC_Li05.csv','NMC_Li06.csv','NMC_Li07.csv','NMC_Li08.csv'}...
    ,['./Fit_Results/NMC_LiClO4/',filenames{i},'/fitted_params.csv'], temp_vec, conc_vec, rxn_type,activity_LiClO4);

temp_vec = [293,293,293,293];
conc_vec = [0.5,0.6,0.7,0.8];
%LPF 
[tht_fit_LCOPF,f_fit_LCOPF,Y_LCOPF] = get_fit_params({'LCO_LiPF6_Li05.csv','LCO_LiPF6_Li06.csv','LCO_LiPF6_Li07.csv','LCO_LiPF6_Li08.csv'}...
    ,['./Fit_Results/LCO_x_LPF/',filenames{i},'/fitted_params.csv'],temp_vec, conc_vec, rxn_type, activity_LPF);
[tht_fit_NMCPF,f_fit_NMCPF,Y_NMCPF] = get_fit_params({'NMC_LiPF6_Li05.csv','NMC_LiPF6_Li06.csv','NMC_LiPF6_Li07.csv','NMC_LiPF6_Li08.csv'},...
    ['./Fit_Results/NMC_LPF/',filenames{i},'/fitted_params.csv'], temp_vec, conc_vec, rxn_type,activity_LPF);

temp_vec = [293,293,293];
conc_vec = [0.6,0.7,0.8];
[tht_fit_NMC811, f_fit_NMC811, Y_NMC811] = get_fit_params({'811_Li06.csv','811_Li07.csv',...
    '811_Li08.csv'},...
    ['./Fit_Results/NMC811_LiClO4_x_dependent/',filenames{i},'/fitted_params.csv'],temp_vec,conc_vec, rxn_type,activity_LiClO4);

Y = [Y_LCO;Y_NMC;Y_LCOPF;Y_NMCPF;Y_NMC811];
f_fit = [f_fit_LCO;f_fit_NMC;f_fit_LCOPF;f_fit_NMCPF;f_fit_NMC811];
tht_fit = [tht_fit_LCO;tht_fit_NMC;tht_fit_LCOPF;tht_fit_NMCPF;tht_fit_NMC811];

sprintf('Reaction:%s\n',rxn_type);
calc_BIC_values(Y,f_fit,size(tht_fit,1));
end


%% Helpers
function [tht_fit,f_fit,Y] = get_fit_params(filenames,paramfile,temp_vec, conc_vec, rxn_type,activity)
[X,Y,size_vec] = get_file_data(filenames);
tht_fit = retreive_params(paramfile);
f_x = @(tht) predict_2(X,tht,size_vec,temp_vec,conc_vec,rxn_type,activity);
[f_fit,~] = f_x(tht_fit);
end

