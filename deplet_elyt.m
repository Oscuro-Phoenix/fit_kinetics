D = 1e-10; % D of elyte 
L = 5e-5; % L of trode
params.tauD = L^2/D; % Diff time
c0 = 1e3; %1 M
A = pi*(12.7/2*1e-3)^2; %m^2 Area
V = A*L;
params.ndimI = 96500*V*c0/params.tauD;
F = get_integrand(-I_arr,t_arr,params);
c = c0*ones(size(t_arr));
for i=4:size(t_arr,1)
c(i) = c0*exp(-t_arr(end)./params.tauD)*(1+trapz(t_arr(1:i)./params.tauD,F(1:i)));
end
plot(c);
function F = get_integrand(I,t,params)
F = (I*1e-3./params.ndimI).*exp(t./params.tauD);
end