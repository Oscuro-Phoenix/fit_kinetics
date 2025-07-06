function [f_fit,tht_fit,cost_val,chi_squared,BIC,size_vec,xdata,ydata,xvals] = get_fit_2(temp_vec, conc_vec, ...
    filenames,rxn_type,activity)

Num_temps = size(temp_vec,2); %Number of temp,conc combinations

%Load data into xdata,x and ydata,y
x = []; y = []; size_vec = [];
for i=1:size(filenames,2)
    disp(filenames{i});
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
num_unique = size(unique(temp_vec),2);

if rxn_type == "ICET_symmetric"
    Nukn = 3+Num_temps;
elseif rxn_type == "ICET_asymmetric"
    Nukn = 4+Num_temps;
elseif rxn_type == "ECIT_normal"
    Nukn = 3+Num_temps;
elseif rxn_type == "ECIT_normal2"
    Nukn = 4;
elseif rxn_type == "ECIT_i0_fix"
    Nukn = 1+Num_temps;
elseif rxn_type == "ECIT_i0_fix_wfix"
    Nukn = 2;
elseif rxn_type == "ECIT_film"
    Nukn = 4+num_unique;
elseif rxn_type == "uniformly_valid"
    Nukn = 4+Num_temps;
elseif rxn_type == "uniformly_valid_i0_fix"
    Nukn = 3+Num_temps;
elseif rxn_type == "uniformly_valid_i0_fix_wfix"
    Nukn = 4;
elseif rxn_type == "uniformly_valid2"
    Nukn = 5;
elseif rxn_type == "ICET_symmetric_film"
    Nukn = 3+Num_temps+num_unique;
elseif rxn_type == "ICET_symmetric_film2"
    Nukn = 3+num_unique;
elseif rxn_type == "ICET_asymmetric_film"
    Nukn = 4+Num_temps+num_unique;
elseif rxn_type == "ICET_asymmetric_film2"
    Nukn = 4+num_unique;
elseif rxn_type == "BV_film"
    Nukn = 1+2*num_unique;
elseif rxn_type == "Linear"
    Nukn = Num_temps;
elseif rxn_type == "Linear_ICET"
    Nukn = 3;
elseif rxn_type == "BV"
    Nukn = 3;
end

disp("Slope Values : \n");
slopes =get_slopes(x,y,size_vec);
disp(slopes);

fun = @(tht) cost_func_2(tht,x,y,size_vec,temp_vec,conc_vec,rxn_type,activity,slopes); %Cost function

%Initial Guess for the parameters (first param is k0)
% if rxn_type == "ICET_asymmetric_film"
% tht_i = [2.40779498306884;1.89100093191450;-1.56622880219696e-14;0.498250095667500;0.872592026366421;0.799271583582608;-7.03707027844717e-07;0.843740549075189;2.04691443944642;1.68360168887711;-0.697293789913788;0.698193805729373;0.852972375231014;1.26630200684295;9.03153692841622]; % Random initial Guess
% else
% if rxn_type == "ECIT_normal" || rxn_type=="ECIT_normal2"
% tht_i = [1e10;5;20;rand(Nukn-3,1)];
% else
if rxn_type == "ICET_asymmetric_film2" || rxn_type == "ICET_symmetric_film2"
    tht_i = [-100;0.1*ones(Nukn-num_unique-2,1);10*ones(num_unique+1,1)];
    lb = [-Inf;-1*ones(Nukn-num_unique-2,1);-5;ones(num_unique,1)];
    ub = [0;1*ones(Nukn-num_unique-2,1);5;100*ones(num_unique,1)];
elseif rxn_type == "ECIT_film"
    tht_i = [-100;ones(Nukn-num_unique-2,1);ones(num_unique+1,1)];
    lb = [-Inf;0;-90*ones(Nukn-num_unique-2,1);zeros(num_unique,1)];
    ub = [0;10;90*ones(Nukn-num_unique-2,1);10*ones(num_unique,1)];
elseif rxn_type == "BV_film"
    tht_i = rand(Nukn,1);
    lb = [-Inf*ones(num_unique,1);0;zeros(num_unique,1)];
    ub = [zeros(num_unique,1);pi/2;50*ones(num_unique,1)];
elseif rxn_type == "Linear" 
    tht_i = rand(Nukn,1);
    lb = [0]*ones(size(tht_i));
    ub = [Inf]*ones(size(tht_i));
elseif rxn_type == "Linear_ICET"
    tht_i = rand(Nukn,1);
    lb = [0;-10;0];
    ub = [Inf;10;Inf];
elseif rxn_type == "BV"
    tht_i = rand(Nukn,1);
    lb = [-Inf;0;-3];
    ub = [0;1;3];
elseif rxn_type == "uniformly_valid2"
    tht_i = rand(Nukn,1);
    lb = [-Inf;2;0;0;-10];
    ub = [0;10;10;10;10];
elseif rxn_type == "uniformly_valid_i0_fix"
    tht_i = rand(Nukn,1);
    lb = [3;0.1;0.1;3*ones(Nukn-3,1)];
    ub = [20;10;10;6.5*ones(Nukn-3,1)];
elseif rxn_type == "ECIT_i0_fix_wfix"
    tht_i = rand(Nukn,1);
    lb = [5;3];
    ub = [10;7];
elseif rxn_type == "uniformly_valid_i0_fix_wfix"
    tht_i = rand(Nukn,1);
    lb = [20;1;1;3];
    ub = [30;5;5;8];
elseif rxn_type == "uniformly_valid"
    tht_i = rand(Nukn,1);
    lb = [-Inf;5;0;0;-2*ones(Nukn-4,1)];
    ub = [0;10;10;10;2*ones(Nukn-4,1)];
elseif rxn_type == "ECIT_i0_fix"
    tht_i = rand(Nukn,1);
    lb = [5;-1*ones(Nukn-1,1)];
    ub = [15;1*ones(Nukn-1,1)]; %lambda, w_ads
else
    tht_i = rand(Nukn,1);
    lb = [-Inf;9.6;0;-5*ones(Nukn-3,1)];
    ub = [0;9.7;0;5*ones(Nukn-3,1)];
end


options = optimoptions('fmincon','MaxFunctionEvaluations',Inf,'MaxIterations',Inf);
problem = createOptimProblem('fmincon','objective',fun,...
    'x0', tht_i,'lb',lb,'ub',ub,'options',options);
gs = GlobalSearch('NumStageOnePoints',1e2,'NumTrialPoints',5e2,'PenaltyThresholdFactor',0.5,...
    'StartPointsToRun','bounds');
[tht_fit,~,~,~] = run(gs,problem);

% options = optimoptions('fminsearch','MaxFunEvals',Inf,'MaxIter',Inf);
%[tht_fit,~] = fminunc(fun,tht_i,options); %Optimize
% tht_fit = tht_fit + 0.01*tht_i;
% [tht_fit,cost_val] = fminunc(fun,tht_fit,options); %Reoptimize

[cost_val,chi_squared] = pure_cost_2(tht_fit,x,y,size_vec,temp_vec,conc_vec,rxn_type,activity,slopes);
BIC = size(xdata,1)*log(cost_val/size(xdata,1))+Nukn*log(size(xdata,1));
xvals = repmat(linspace(-5,5,200)',size(filenames,2),1);
f_x = @(tht) predict_2(xvals,tht,200*ones(size(filenames,2),1),temp_vec,conc_vec,rxn_type,activity,slopes);
[f_fit,~] = f_x(tht_fit);
Z = get_gradient(tht_fit,xdata,size_vec,temp_vec,conc_vec,rxn_type,activity,slopes);
J = get_gradient(tht_fit,xvals,200*ones(size(filenames,2),1),temp_vec,conc_vec,rxn_type,activity,slopes);
FI = pinv(Z'*Z)*cost_val/(size(xdata,1)-Nukn);
disp(sqrt(FI));
print_params_2(tht_fit,rxn_type,num_unique);
updated_data = [f_fit(:), xvals(:)];
writematrix(updated_data, 'fitted_data.csv');
% Check if file exists
if exist('fitted_params.csv', 'file')
    existing_params = readmatrix('fitted_params.csv');
    updated_params = [existing_params; tht_fit(:)];
else
    updated_params = tht_fit(:);
end
writematrix(updated_params, 'fitted_params.csv');

end