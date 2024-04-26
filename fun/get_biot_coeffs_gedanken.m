function [P,Q,R]=get_biot_coeffs_gedanken(phi,Kb,Kf,N,Ks)
% test call1:    get_biot_coeffs_gedanken(0.9,10.1e6,142e3, 6.2e6)

switch nargin
    
    case 5
        P = ( (1-phi)*(1-phi-Kb/Ks)*Ks + (phi*Ks*Kb)/Kf  ) / ...
            (  1 - phi - Kb/Ks +  phi*Ks/Kf              ) ...
            +...
            4/3*N;
        
        Q = ( (1-phi-Kb/Ks)*phi*Ks                       )/ ...
            (  1 - phi - Kb/Ks +  phi*Ks/Kf              );
        
        R = (   phi^2 * Ks                               )/ ...
            (  1 - phi - Kb/Ks +  phi*Ks/Kf              );
    
    case 4   % Ks ->  inf (the material the frame is made of 
             %             is not compressible)
        P =  (1-phi)^2/phi*Kf  +  Kb +  4/3*N;
        
        Q = ( (1-phi)*Kf                       );
        
        R = (   phi * Kf                                 );  

end

end