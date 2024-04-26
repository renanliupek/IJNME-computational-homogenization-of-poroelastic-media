wfunction [rho_11,rho_12,rho_22]=get_biot_coeffs_densities(phi, rhof, rhos, alpha_inf)
% test call1:    get_biot_coeffs_gedanken(0.9,10.1e6,142e3, 6.2e6)

rho_11 =  rhos * (1-phi) + rhof*phi*(alpha_inf-1);
rho_22 =  rhof*phi*alpha_inf;
rho_12 = -rhof*phi*(alpha_inf-1);

end