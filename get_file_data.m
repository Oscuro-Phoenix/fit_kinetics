function [xdata,ydata,size_vec] = get_file_data(filenames)
%Load data into xdata,x and ydata,y
x = []; y = []; size_vec = [];
for i=1:size(filenames,2)
    M = csvread(filenames{i},1,0);  % from 1st line, 1st row
    fin = size(M,1);
    last_row=find(M(:,1), 1, 'last'); %search last non-zero entry for experiments
    %%   Loaded Data
    x = [x;M(1:last_row,1)];
    y = [y;sign(M(1:last_row,1)).*M(1:last_row,2)];
    size_vec = [size_vec; last_row];
end
xdata = x;
ydata = y;
end