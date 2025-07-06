clearvars; clc;
temp_vec = [293,293,293,293];
conc_vec = [0.5,0.6,0.7,0.8];
activity = 1.6;
filenames = {'./LCO/LCO_Li05.csv'};
tht_fit = readmatrix('fitted_params.csv');
rxn_type = 'ICET_asymmetric_film2';
param_num = 1;
scalefacs = [0.1,0.2,0.5,1,2,5,10];

% %% Plot baseline
% betas = 0.01:0.5:3;
% x = -20:0.01:20;
% beta_n = {'\beta = 0.01','\beta = 0.5','\beta = 1','\beta = 1.5','\beta = 2','\beta = 2.5'};
figure(1);
set(findall(gcf,'-property','FontSize'),'FontSize',32);
size_vec_2 = 200*ones(size(filenames,2),1);

for i=1:size(scalefacs,2)
    tht_curr = tht_fit;
    tht_curr(param_num) = tht_curr(param_num)*scalefacs(i);
    [xvalues,f_fit,xdata,ydata,size_vec] = get_fitted_curent(tht_curr,filenames,temp_vec,conc_vec,rxn_type,activity);
    [X_vals,idx] = sort(xvalues(1:size_vec_2(1)));
    [X_vals2,idx2] = sort(xdata(1:size_vec(1)));
    Y_vals = f_fit(1:size_vec_2(1));
    plot(X_vals,(Y_vals(idx)),'LineWidth',3);
     hold on;
end
hold on;
scatter(xdata,ydata);
legend(string(scalefacs)+"x",'Location','northwest');

%% Enhance plot appearance
xlabel('Overpotential (dimless)', 'FontSize', 12, 'FontWeight', 'bold'); % X-axis label
ylabel('Current Density (A/cm^2)', 'FontSize', 12, 'FontWeight', 'bold'); % Y-axis label
title('Influence of Changing pre-factor (fixed R_{film})', 'FontSize', 14, 'FontWeight', 'bold'); % Plot title
%legend(beta_n,'Location','southeast');
%% Set axis properties
ax = gca; % Get current axis
ax.FontSize = 10; % Font size for axis labels
ax.LineWidth = 1; % Line width for axis lines
ax.XColor = 'k'; % Color of X-axis
ax.YColor = 'k'; % Color of Y-axis

%% Set grid and box properties
grid on; % Turn on the grid
box on; % Turn on the box around the plot

%% Set figure properties
set(gcf, 'PaperPositionMode', 'auto'); % Ensure figure saves in the same size as on screen
set(gcf, 'Position', [100, 100, 600, 400]); % Size of the figure (Width x Height)

function I= get_curr_arr(x,b_ox)
I = zeros(size(x));
activity = 1.6;
w_ads = 0.5;
b_red = b_ox;
b_ox = b_ox;
lambda = 5;
cc = 0.5;
T = 298;
k0 = -1*exp(b_ox/2);
c_fct =1;
for i=1:size(x,2)
    I(i) = uniformly_valid_old(x(i),k0,lambda,b_red,b_ox,w_ads,cc,c_fct,activity,T);
end
end

function [xvals,f_fit,xdata,ydata,size_vec] = get_fitted_curent(tht_fit,filenames,temp_vec,conc_vec,rxn_type,activity)
    %Load data into xdata,x and ydata,y
    x = []; y = []; size_vec = [];
    for i=1:size(filenames,2)
        M = csvread(filenames{i},1,0);  % from 1st line, 1st row
        last_row=find(M(:,1), 1, 'last'); %search last non-zero entry for experiments
        %%   Loaded Data
        x = [x;M(1:last_row,1)];
        y = [y;sign(M(1:last_row,1)).*M(1:last_row,2)];
        size_vec = [size_vec; last_row];
    end
    xdata = x;
    ydata = y;

    xvals = repmat(linspace(-8,8,200)',size(filenames,2),1);
    f_x = @(tht) predict_2(xvals,tht,200*ones(size(filenames,2),1),temp_vec,conc_vec,rxn_type,activity);
    [f_fit,~] = f_x(tht_fit);
end
