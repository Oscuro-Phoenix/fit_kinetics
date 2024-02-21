clearvars;
%Load data into xdata,x and ydata,y
filenames = {'./x_0.60/exp1.csv','./x_0.60/exp2.csv','./x_0.60/exp3.csv',...
    './x_0.60/exp4.csv','./x_0.60/exp5.csv','./x_0.60/exp6.csv',...
    './x_0.60/exp7.csv','./x_0.60/exp8.csv','./x_0.60/exp9.csv',...
    './x_0.60/exp10.csv','./x_0.60/exp11.csv','./x_0.60/exp12.csv',...
    './x_0.60/exp13.csv','./x_0.60/exp14.csv'};


x = [];
y = [];
for i=1:size(filenames,2)
    M = csvread(filenames{i},1,0);  % from 1st line, 1st row
    fin = size(M,1);
    last_row=find(M(:,1), 1, 'last'); %search last non-zero entry for experiments
    %%   Loaded Data
    x = [x;M(1:last_row,1)/0.0257];
    y = [y;M(1:last_row,2)];
end
xdata = x;
ydata = y;

M = csvread('./ICET_asymm_film_fit/fitted_data.csv',1,0);  % from 1st line, 1st row
fin = size(M,1);
last_row=find(M(:,1), 1, 'last'); %search last non-zero entry for experiments
%%   Loaded Data
x_fit = M(last_row-400:last_row-200,1);
y_fit = M(last_row-400:last_row-200,2);
F = abs(x_fit)<8;
figure
scatter(x, y,'filled')
hold on
scatter(y_fit(F),x_fit(F),'filled','r')
shadeWidth = 0.1;
patch([x; flipud(x)], [min(y,[],2)-shadeWidth; flipud(max(y,[],2))+shadeWidth], [1 1 1]*0.5, 'FaceAlpha',0.25, 'EdgeColor','none')
hold off
xlabel('Dimensionless Overpotential','FontSize',16);
ylabel('Current (mA)','FontSize',16);
legend('Data','Fit');