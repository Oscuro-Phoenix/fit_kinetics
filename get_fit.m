
function [f_fit,tht_fit,cost_val,BIC,size_vec,xdata,ydata] = get_fit(temp_vec, conc_vec, filenames,rxn_type)
% temp_vec = [298,298,298,298];
% conc_vec = [0.5,0.6,0.7,0.8];
% filenames = {'LCO_LiPF6_Li05.csv','LCO_LiPF6_Li06.csv','LCO_LiPF6_Li07.csv','LCO_LiPF6_Li08.csv'};
% temp_vec = [298,298,298];
% conc_vec = [0.4,0.6,0.8];
% filenames = {'LFP_Li04.csv','LFP_Li06.csv','LFP_Li08.csv'};

% Load files
% % Change filenames here
% temp_vec = [298,298,298,298,308,318]; %Temperature values
% conc_vec = [0.5,0.6,0.7,0.8,0.5,0.5]; %Concentration values
% filenames = {'LCO_Li05.csv','LCO_Li06.csv','LCO_Li07.csv','LCO_Li08.csv',...
%     'LCO_LiClO4_T_dependent_35C.csv','LCO_LiClO4_T_dependent_45C.csv'}; %Corresponding filenames


% temp_vec = [298,298,298,298];
% conc_vec = [0.5,0.6,0.7,0.8];
% filenames = {'NMC_Li05.csv','NMC_Li06.csv','NMC_Li07.csv','NMC_Li08.csv'};

% temp_vec = [298,298,298,298];
% conc_vec = [0.5,0.6,0.7,0.8];
% filenames = {'NMC_LiPF6_Li05.csv','NMC_LiPF6_Li06.csv','NMC_LiPF6_Li07.csv','NMC_LiPF6_Li08.csv'};

Num_temps = size(temp_vec,2); %Number of temp,conc combinations

%Load data into xdata,x and ydata,y
x = []; y = []; size_vec = [];
for i=1:size(filenames,2)
    M = csvread(filenames{i},1,0);  % from 1st line, 1st row
    fin = size(M,1);
    last_row=find(M(:,1), 1, 'last'); %search last non-zero entry for experiments
    %%   Loaded Data
    x = [x;M(1:last_row,1)];
    y = [y;M(1:last_row,2)];
    size_vec = [size_vec; last_row];
end
xdata = x;
ydata = y;

Nukn = 0;

if rxn_type == "ICET_symmetric"
    Nukn = 2+2*Num_temps;
elseif rxn_type == "ICET_asymmetric"
    Nukn = 2+3*Num_temps;
elseif rxn_type == "ECIT_normal"
    Nukn = 2+2*Num_temps;
elseif rxn_type == "ECIT_normal2"
    Nukn = 3+Num_temps;
elseif rxn_type == "uniformly_valid"
    Nukn = 2+3*Num_temps;
elseif rxn_type == "uniformly_valid2"
    Nukn = 3+2*Num_temps;
elseif rxn_type == "ICET_symmetric_film"
    Nukn = 3+2*Num_temps;
elseif rxn_type == "ICET_asymmetric_film"
    Nukn = 3+3*Num_temps;
end

fun = @(tht) cost_func(tht,x,y,size_vec,temp_vec,conc_vec,rxn_type); %Cost function

%Initial Guess for the parameters (first param is k0)
if rxn_type == "ICET_asymmetric_film"
tht_i = [2.40779498306884;1.89100093191450;-1.56622880219696e-14;0.498250095667500;0.872592026366421;0.799271583582608;-7.03707027844717e-07;0.843740549075189;2.04691443944642;1.68360168887711;-0.697293789913788;0.698193805729373;0.852972375231014;1.26630200684295;9.03153692841622]; % Random initial Guess
else
tht_i = [1;rand(Nukn-1,1)];
end

options = optimoptions('fminunc','MaxFunctionEvaluations',1000000);

[tht_fit,cost_val] = fminunc(fun,tht_i,options); %Optimize

BIC = size(xdata,1)*log(cost_val/size(xdata,1))+Nukn*log(size(xdata,1));

f_x = @(tht) predict(xdata,tht,size_vec,temp_vec,conc_vec,rxn_type);
f_fit = f_x(tht_fit);
print_params(tht_fit,rxn_type);
end