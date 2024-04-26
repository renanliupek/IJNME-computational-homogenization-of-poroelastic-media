function alpha_inf = get_analytical_tortuosity(r1,l1,r2,l2)
  %design1:          get_analytical_tortuosity(5e-6, 5e-6, 4.75e-5,  2*4.75e-5)
  %design7:          get_analytical_tortuosity(2.5e-5, 2.5e-5, 4.75e-5, 7.5e-5)


% assume out-of-plane dimension = 1
S1 = r1*1;
S2 = r2*1;

% calculate tortuosity according to Allard (5.23), check the assumptions
alpha_inf = (l1*S2+l2*S1)*(l1*S1+l2*S2) / ( (l1+l2)^2*S1*S2 );

end