function calc_BIC_values(trueValues,predictedValues, numParameters)

% Number of data points
N = numel(trueValues);

% Residuals
residuals = trueValues - predictedValues;

% Mean squared error
mse = sum(residuals.^2) / N;

% Number of parameters in the model (adjust this based on your model)

% BIC calculation
bic = N * log(mse) + numParameters * log(N);

fprintf('BIC value: %.4f\n', bic);
end
