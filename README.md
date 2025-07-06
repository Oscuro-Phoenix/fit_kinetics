# Fit Kinetics

A MATLAB-based electrochemical kinetics fitting toolbox for analyzing and modeling electrochemical reaction data. This project provides tools for fitting experimental current-voltage data to various electrochemical reaction models including ICET (Ion Coupled Electron Transfer), ECIT (Electron Coupled Ion Transfer), BV (Butler-Volmer), and uniformly valid models.

## Overview

This toolbox is designed for researchers and engineers working with electrochemical systems, particularly those studying electrode kinetics, charge transfer processes, and electrochemical reactions. It provides a comprehensive framework for:

- Fitting experimental current-voltage data to theoretical models
- Comparing different reaction mechanisms
- Analyzing kinetic parameters
- Visualizing fitted results
- Statistical analysis of model performance

## Supported Reaction Models

### 1. ICET (Ion Coupled Electron Transfer) Models when Reactions are IT Limited
- **ICET_symmetric**: Symmetric ICET with equal forward/reverse barriers
- **ICET_asymmetric**: Asymmetric ICETwith different forward/reverse barriers
- **ICET_symmetric_film**: ICET with film resistance effects
- **ICET_asymmetric_film**: Asymmetric ICET with film resistance
- **ICET_symmetric_film2**: Simplified symmetric ICET with film
- **ICET_asymmetric_film2**: Simplified asymmetric ICET with film

### 2. ECIT (Electron Coupled Ion Transfer) Models when Reactions are ET Limited
- **ECIT_normal**: Standard electrochemical charge transfer
- **ECIT_normal2**: Modified ECIT model
- **ECIT_i0_fix**: ECIT with fixed exchange current
- **ECIT_i0_fix_wfix**: ECIT with fixed exchange current and adsorption
- **ECIT_film**: ECIT with film resistance effects

### 3. Butler-Volmer Models
- **BV**: Standard Butler-Volmer equation
- **BV_film**: Butler-Volmer with film resistance

### 4. Uniformly Valid Models
- **uniformly_valid**: Uniformly valid approximation
- **uniformly_valid_i0_fix**: Uniformly valid with fixed exchange current
- **uniformly_valid_i0_fix_wfix**: Uniformly valid with fixed exchange current and adsorption
- **uniformly_valid2**: Modified uniformly valid model

### 5. Linear Models
- **Linear**: Linear current-voltage relationship
- **Linear_ICET**: Linear approximation of ICET

## Key Functions

### Main Fitting Functions
- `get_fit_2.m`: Main fitting function that optimizes model parameters
- `predict_2.m`: Predicts current values using fitted parameters
- `print_params_2.m`: Displays fitted parameters in a readable format

### Model Implementation Functions
- `ICET.m`: Interfacial Charge Transfer model
- `ECIT.m`: Electrochemical Charge Transfer model
- `BV.m`: Butler-Volmer model
- `uniformly_valid.m`: Uniformly valid approximation
- `Linear.m`: Linear model

### Utility Functions
- `get_slopes.m`: Calculates experimental slopes
- `get_gradient.m`: Computes gradients for optimization
- `cost_func_2.m`: Cost function for optimization
- `pure_cost_2.m`: Pure cost calculation
- `plot_data.m`: Visualization of experimental and fitted data

## Data Format

### Input Data
Experimental data should be provided as CSV files with the following format:
- Column 1: Overpotential (V) or dimensionless overpotential
- Column 2: Current (mA) or current density

Example:
```csv
-0.5,2.3
-0.4,3.1
-0.3,4.2
...
```

### Output Files
The fitting process generates:
- `fitted_data.csv`: Predicted current values
- `fitted_params.csv`: Optimized model parameters

## Usage Example

### Basic Usage
```matlab
% Define experimental conditions
temp_vec = [298, 308, 318];  % Temperature in Kelvin
conc_vec = [0.1, 0.1, 0.1]; % Concentration
filenames = {'exp1.csv', 'exp2.csv', 'exp3.csv'};
rxn_type = "ICET_asymmetric";
activity = 1.6;  % Activity coefficient

% Perform fitting
[f_fit, tht_fit, cost_val, chi_squared, BIC, size_vec, xdata, ydata, xvals] = ...
    get_fit_2(temp_vec, conc_vec, filenames, rxn_type, activity);

% Print results
print_params_2(tht_fit, rxn_type, length(unique(temp_vec)));
```

### Plotting Results
```matlab
% Load and plot experimental data
plot_data;  % Uses predefined file paths

% Or customize plotting
figure;
scatter(xdata, ydata, 'filled', 'DisplayName', 'Experimental');
hold on;
scatter(xvals, f_fit, 'filled', 'r', 'DisplayName', 'Fitted');
xlabel('Overpotential (V)');
ylabel('Current (mA)');
legend;
```

## Model Parameters

### Common Parameters
- `k0`: Exchange current density
- `lambda`: Reorganization energy
- `beta_red/beta_ox`: Reduction/oxidation transfer coefficients
- `w_ads`: Adsorption energy
- `R_film`: Film resistance (for film models)
- `alpha`: Transfer coefficient (for BV models)
- `dG_IT`: Free energy change (for ECIT models)

### Parameter Constraints
Each model has specific parameter constraints to ensure physical meaning:
- Exchange current densities: Positive values
- Transfer coefficients: Between 0 and 1
- Reorganization energies: Positive values
- Film resistances: Positive values

## Optimization

The fitting process uses:
- **Global Search**: Multi-start optimization to find global minimum
- **Constrained Optimization**: Physical constraints on parameters
- **Cost Function**: Weighted least squares with regularization
- **Statistical Analysis**: BIC (Bayesian Information Criterion) for model comparison

## Statistical Analysis

### Goodness of Fit Metrics
- **Chi-squared**: Measure of fit quality
- **BIC**: Model comparison criterion (lower is better)
- **RMSE**: Root mean square error
- **R-squared**: Coefficient of determination

### Model Selection
Use BIC values to compare different models:
- Lower BIC indicates better model
- Consider physical meaning of parameters
- Check parameter uncertainties

## Requirements

- MATLAB R2018b or later
- Optimization Toolbox
- Global Optimization Toolbox
- Statistics and Machine Learning Toolbox

## File Structure

```
fit_kinetics/
├── README.md                 # This file
├── get_fit_2.m              # Main fitting function
├── predict_2.m              # Prediction function
├── print_params_2.m         # Parameter display
├── get_params_2.m           # Parameter extraction
├── cost_func_2.m            # Cost function
├── pure_cost_2.m            # Pure cost calculation
├── plot_data.m              # Visualization
├── get_slopes.m             # Slope calculation
├── get_gradient.m           # Gradient computation
├── fitted_data.csv          # Output: fitted data
├── fitted_params.csv        # Output: fitted parameters
├── ICET.m                   # ICET model
├── ECIT.m                   # ECIT model
├── BV.m                     # Butler-Volmer model
├── uniformly_valid.m        # Uniformly valid model
├── Linear.m                 # Linear model
└── [other model files]      # Additional model implementations
```

## Example Data

The project includes example data in the `x_0.60/` directory (referenced in `plot_data.m`):
- 14 experimental datasets (exp1.csv through exp14.csv)
- Current-voltage measurements at different conditions
- Used for testing and demonstration

## Contributing

To contribute to this project:
1. Add new reaction models in separate `.m` files
2. Update `get_fit_2.m` to include new model types
3. Add parameter extraction in `get_params_2.m`
4. Update parameter printing in `print_params_2.m`
5. Test with experimental data

## Citation

If you use this toolbox in your research, please cite the relevant electrochemical kinetics literature and include a reference to this software.

## License

This project is provided as-is for research and educational purposes. Please ensure proper attribution when using in publications.

## Contact

For questions or issues, please refer to the code comments or create an issue in the repository.
