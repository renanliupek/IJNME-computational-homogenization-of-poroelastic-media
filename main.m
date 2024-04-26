%% DESCRIPTION
% Renan Liupekevicius Carnielli TU/e
% start in 29-09-2022
% computational homogenization of poroelastic medium including internal
% dynamics of the RVE.

% TOPOLOGY: Fluid at top-bottom-right-left and solid a must at the corners
% TRANSPARENT ELEMENTS:Here you can handle the transparent massless 
% solid elements
% PERIODICITY ONLY IN SOLID DOMAIN


clear; close all;



%% CHOOSE DESIGN, RVE SIZE, STRUT THICKNESS

disp('codeblock: SELECT DESIGN')
% select unit cell design
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
design = 'hexagon_loc' ;
disp(design);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

%% CONVERT TO QUADRATIC MESH

% read mesh text file comsol 5.4
%--------------------------------------------------------------------------
disp('codeblock: READ & CONVERT LINEAR TO QUADRATIC ELEMENT MESH')
% read text file
l   = read_mphtxt_54(design);

% convert to quadratic element
mesh=mesh2d_lin2qua_uc(l);
%--------------------------------------------------------------------------


% read mesh text file comsol 5.6
%--------------------------------------------------------------------------
% disp('codeblock: READ MESH COMSOL 5.6')
% mesh   = read_mphtxt_56(design);
%--------------------------------------------------------------------------

% copy to local variables (please double click on 'mesh' struct in 
% workspace to see the meaning of each cell)
%--------------------------------------------------------------------------
  x    = mesh.x{2,2};    % tensor form
  mx   = mesh.x{2,1};    % matrix form
  conn = mesh.elem{2,1}; % conn matrix of element
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% mesh parameters
  m  = size(  conn, 1); % number of elements in the mesh
  n  = size(  x, 1);    % number of nodes in the mesh 
  fprintf('NUMBER OF MESH ELEMENTS: %d\n', m)
  fprintf('NUMBER OF MESH NODES: %d\n', n)
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit cell size
 % unit cell size


  ax = 2*max(  mx(:,1));
  ay = 2*max(  mx(:,2));
  

if(max(mx(:,1))+min(mx(:,1))>1e-8); warning('RVE is not x centralized');end
if(max(mx(:,1))+min(mx(:,1))>1e-8); warning('RVE is not y centralized');end  

if(abs(ax-ay)<1e-8); a=ax;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%% DESIGNS: LIST OF TAGS (MANUALLY ADD FOR EACH DESIGN)
disp('codeblock: LOAD DESIGN TAGS')
% MANUALLY import tags from comsol: check acoustic-structure boundary
% list under 'boundary selection'. Same procedure for selecting tag-list of
% the solid phase(s).

switch sscanf(design, '%c')
    
   
    case 'hexagon' % 
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1 2 4:8 11:18 20:23 26:36];
    % solid1 elements tag
    list_solid1_tag = [1 2 4 5 8 15:17 20 21 23 30:33 36];                   
    % solid elements tag
    list_solid2_tag = [];  
    % solid 3 
    list_solid3_tag = [];
    % solid 4 =fluid elements that are the cell opening
    list_solid4_tag = [6 7 11:14 18 22 26:29 34 35];
    % interface edge elements tag  
    list_itr_edges_tag=[6 8 12 16 18 20:22 25 27 29 31 33 37:41 43:46 48:50 54 58 60 62 64 65 67:69 71 73 75 79:83 85 87:89 91 93 96 97];

    % dry or wet?
    dry = false;
    % multiples solid phases
    multiple_solid_phase = true;
    %----------------------------------------------------------------------



    case 'hexagon_loc' % 
    %----------------------------------------------------------------------
    % solid elements tag
 % solid elements tag
    list_solid_tag = [1 2 4:8 11:20 22:27 30:40];
    % solid1 elements tag
    list_solid1_tag = [1 2 4 5 8 15 16 18 22 23 25 34:37 40];                   
    % solid elements tag
    list_solid2_tag = [20 26];  
    % solid 3
    list_solid3_tag = [17 27];
    % solid 4 =fluid elements that are the cell opening
    list_solid4_tag = [6 7 11:14 19 24 30:33 38 39];
    % interface edge elements tag  
    list_itr_edges_tag=[6 8 12 16 18 20:22 25 27 29 31 33 37:44 46:49 51:54 60 64 66 72 74:78 80:82 84 86 88 90 93:97 99 101:103 105 107 110 111];

    % dry or wet?
    dry = false;
    % multiples solid phases
    multiple_solid_phase = true;
    %----------------------------------------------------------------------
    
end

%% INITIAL DEFINITIONS: MATERIALS, TENSOR BASIS AND CANONICAL TENSORS
disp('codeblock: LOAD MATERIAL PROPERTIES')
%--------------------------------------------------------------------------
dim   = 2;  % problem dimension
thickout = 1;  % thickness (out-of-plane direction), in [m]
%--------------------------------------------------------------------------


% fluid properties
%--------------------------------------------------------------------------
warning('fluid material parameter: air')
rhof     = 1.225 ; % [kg/m3]   
c        = 343 ;  % [m/s]
% pA       = 101325; [Pa]   % atmosferic p not used for order-1 equations.
%--------------------------------------------------------------------------


% fluid properties
%--------------------------------------------------------------------------
% rhof     = 200 ; % [kg/m3]
% % rhof     = 1.225e3 ; % [kg/m3]    TEST !!! PERHAPS BUG SOMEWHERE
% c        = 400 ;  % [m/s]
% % pA       = 101325; [Pa]   % atmosferic p not used for order-1 equations.
% warning('changed air material parameters');
%--------------------------------------------------------------------------


% IF METAFOAM DESIGN
    if multiple_solid_phase
    
    n_solid_phases = 4;
   
    %----------------------------------------------------------------------
    warning('solid material parameters: matrix=epoxy, rubber=PU, heavymass=lead, massless=epoxy')
    % define solid phase parameters
    % matrix (part 1)
    matProp{1}.E       =   3.6e9 ; % Young's modulus
    matProp{1}.nu      =   0.368; % poisson ratio
    matProp{1}.rho     =   1180; % mass density
    matProp{1}.G       =   matProp{1}.E/2/(1 +   matProp{1}.nu);
    matProp{1}.kappa   =   matProp{1}.E/3/(1 - 2*matProp{1}.nu);

    % rubber-like (part 2)
    matProp{2}.E       =   50e6; warning('PU diff than foam case')
    matProp{2}.nu      =   0.4; % poisson ratio
    matProp{2}.rho     =   1000; % mass density
    matProp{2}.G       =   matProp{2}.E/2/(1 +   matProp{2}.nu);
    matProp{2}.kappa   =   matProp{2}.E/3/(1 - 2*matProp{2}.nu);

    % heavy mass (part 3)
    matProp{3}.E       =   40.82e9 ; % original
    matProp{3}.nu      =   0.37; %
    matProp{3}.rho     =   11600;
    % matProp{3}.rho     =   3000; % test light mass 
    matProp{3}.G       =   matProp{3}.E/2/(1 +   matProp{3}.nu);
    matProp{3}.kappa   =   matProp{3}.E/3/(1 - 2*matProp{3}.nu);


    % opening massless, same elasticity as material 1 (part 4) 
    matProp{4}.E       =   matProp{1}.E; % Young's modulus
    matProp{4}.nu      =   matProp{1}.nu; % poisson ratio
    matProp{4}.rho     =   0; % mass density
%     matProp{4}.rho     =   matProp{1}.rho; warning('opening has mass equals PU');
    matProp{4}.G       =   matProp{4}.E/2/(1 +   matProp{4}.nu);
    matProp{4}.kappa   =   matProp{4}.E/3/(1 - 2*matProp{4}.nu);
    %----------------------------------------------------------------------


    % warning('TEST: solid material parameters: matrix=epoxy and openings with mass')
    % % define solid phase parameters
    % %----------------------------------------------------------------------
    % % matrix (part 1)
    % matProp{1}.E       =   3.6e9 ; % Young's modulus
    % matProp{1}.nu      =   0.368; % poisson ratio
    % matProp{1}.rho     =   1180; % mass density
    % matProp{1}.G       =   matProp{1}.E/2/(1 +   matProp{1}.nu);
    % matProp{1}.kappa   =   matProp{1}.E/3/(1 - 2*matProp{1}.nu);
    % 
    % % rubber-like (part 2)
    % matProp{2}.E       =   1e6; % Young's modulus
    % matProp{2}.nu      =   0.4; % poisson ratio
    % matProp{2}.rho     =   1000; % mass density
    % matProp{2}.G       =   matProp{2}.E/2/(1 +   matProp{2}.nu);
    % matProp{2}.kappa   =   matProp{2}.E/3/(1 - 2*matProp{2}.nu);
    % 
    % % heavy mass (part 3)
    % matProp{3}.E       =   40.82e9 ; % original
    % matProp{3}.nu      =   0.37; %
    % matProp{3}.rho     =   11600;
    % % matProp{3}.rho     =   3000; % test light mass 
    % matProp{3}.G       =   matProp{3}.E/2/(1 +   matProp{3}.nu);
    % matProp{3}.kappa   =   matProp{3}.E/3/(1 - 2*matProp{3}.nu);
    % 
    % 
    % % opening massless, same elasticity as material 1 (part 4) 
    % matProp{4}.E       =   matProp{1}.E; % Young's modulus
    % matProp{4}.nu      =   matProp{1}.nu; % poisson ratio
    % % matProp{4}.rho     =   0; % mass density
    % matProp{4}.rho     =   matProp{1}.rho; warning('opening has mass equals mat1');
    % matProp{4}.G       =   matProp{4}.E/2/(1 +   matProp{4}.nu);
    % matProp{4}.kappa   =   matProp{4}.E/3/(1 - 2*matProp{4}.nu);
    % %----------------------------------------------------------------------


    % ELSE: 1 SOLID PHASE
else
    %--------------------------------------------------------------------------
    %solid material properties PU from Mirka ch3 (IN TERMS OF E AND nu)
      rhos     =  1180;   % in [kg/m3]
      E        =  3.6e9 ;    % [Pa]
      nu       =  0.368;    % []    
      G        =   E/2/(1+    nu);
      kappa    =   E/3/(1-2*  nu);
    %--------------------------------------------------------------------------
%     C4_isotropicRVE = [kappa+4/3*G kappa-2/3*G 0 0
%                          kappa-2/3*G kappa+4/3*G 0 0
%                          0           0           G G
%                          0           0           G G];
end

   


% define tensor basis
%--------------------------------------------------------------------------
b  = {'e1'; 'e2'};
ee = cartesianbasis2d(  b{1},   b{2});
e1 = ee(1);
e2 = ee(2);
%--------------------------------------------------------------------------


% calculate fundamental tensors
%--------------------------------------------------------------------------
  I   = identity(2, ee);
  I4  = identity(4, ee);
  I4S = 1/2 * (  I4 + rtranspose( I4));
%--------------------------------------------------------------------------




%% SORTING SOLID ELEMENTS

% TOTAL SOLID ELEMENTS 
%--------------------------------------------------------------------------    
    % solid elements
      solid_elems =[];
      for i=list_solid_tag
      solid_elems = [solid_elems find( mesh.elem{2,2}==i).'];
      end
    % connectivity of the solid nodes
      conns       =   conn(  solid_elems,:);
    % solid nodes
      nodes_s     = unique(reshape(  conns,1,[]));
    % number of solid elements in the mesh
      ms          = size(  conns,1);
    % number of solid nodes
      ns          = length(  nodes_s);
%--------------------------------------------------------------------------


if multiple_solid_phase
%--------------------------------------------------------------------------

    
    % SOLID 1 - MATRIX
    %---------------------------------------------------------------------- 
        % solid elements
          solid1_elems =[];
          for i=list_solid1_tag
          solid1_elems = [solid1_elems find( mesh.elem{2,2}==i).'];
          end
        % connectivity of the solid nodes
          connss{1}       =   conn(  solid1_elems,:);
        % solid nodes
          nodes_ss{1}     = unique(reshape(  connss{1},1,[]));
        % number of solid elements in the mesh
          mss{1}          = size  (  connss{1}, 1);
        % number of solid nodes
          nss{1}          = length(  nodes_ss{1});
    %----------------------------------------------------------------------
    
    % SOLID 2 - RUBBER-LIKE
    %----------------------------------------------------------------------
        % solid elements
          solid2_elems =[];
          for i=list_solid2_tag
          solid2_elems = [solid2_elems find( mesh.elem{2,2}==i).'];
          end
        % pconnectivity of the solid nodes
          connss{2}       =   conn(  solid2_elems,:);
        % solid nodes
          nodes_ss{2}     = unique(reshape(  connss{2},1,[]));
        % number of solid elements in the mesh
          mss{2}          = size  (  connss{2}, 1);
        % number of solid nodes
          nss{2}          = length(  nodes_ss{2});
    %----------------------------------------------------------------------

    % SOLID 3 - HEAVY MASS
    %----------------------------------------------------------------------  
        % solid elements
          solid3_elems =[];
          for i=list_solid3_tag
          solid3_elems = [solid3_elems find( mesh.elem{2,2}==i).'];
          end
        % pconnectivity of the solid nodes
          connss{3}       =   conn(  solid3_elems,:);
        % solid nodes
          nodes_ss{3}     = unique(reshape(  connss{3},1,[]));
        % number of solid elements in the mesh
          mss{3}          = size  (  connss{3}, 1);
        % number of solid nodes
          nss{3}          = length(  nodes_ss{3});
    %----------------------------------------------------------------------

    % SOLID 4 - OPENING
    %----------------------------------------------------------------------  
        % solid elements
          solid4_elems =[];
          for i=list_solid4_tag
          solid4_elems = [solid4_elems find( mesh.elem{2,2}==i).'];
          end
        % pconnectivity of the solid nodes
          connss{4}       =   conn(  solid4_elems,:);
        % solid nodes
          nodes_ss{4}     = unique(reshape(  connss{4},1,[]));
        % number of solid elements in the mesh
          mss{4}          = size  (  connss{4}, 1);
        % number of solid nodes
          nss{4}          = length(  nodes_ss{4});
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------


%% SORTING FLUID ELEMENTS

% OPENING both fluid and solid (already included in conns)
%--------------------------------------------------------------------------
    % opening elements
      opening_elems =[];
      for i=list_solid4_tag
      opening_elems = [opening_elems find( mesh.elem{2,2}==i).'];
      end
    % connectivity of the opening nodes
      conno       =   conn(  opening_elems,:);
    % opening nodes
      nodes_o     = unique(reshape(  conno,1,[]));
    % number of solid elements in the mesh
      mo          = size  (  conno,1);
    % number of solid nodes
      no          = length(  nodes_o);

    % if dry, then no fluid phase at all
    if dry opening_elems =[];end
%--------------------------------------------------------------------------


% FLUID
%--------------------------------------------------------------------------
    % fluid elements
      fluid_elems = setdiff(1:  m,  solid_elems);
    % add opening elements
%       fluid_elems = [fluid_elems opening_elems];
      fluid_elems = union(fluid_elems,opening_elems);
    % conn fluid
      connf       =   conn(  fluid_elems,:);
    % fluid nodes
      nodes_f     = unique(reshape(  connf,1,[]));
    % number of fluid elements in the mesh
      mf          = size(  connf,1);
    % number of fluid nodes
      nf          = length(  nodes_f);
%--------------------------------------------------------------------------

clear conn;


%% INTERFACE NODES (GET TAGS FROM COMSOL MODEL)
%--------------------------------------------------------------------------


% display interface settings
%--------------------------------------------------------------------------
disp(['PLEASE, CHECK FLUID-STRUCTURE INTERFACE & SOLID ELEMENT TAGS ' ...
      'IN THE CORRESPONDING COMSOL MODEL']);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
 itr_edges  = []; % edge stands for 'edge element' or 'surface element'
for i= list_itr_edges_tag
    itr_edges  = [  itr_edges find(mesh.edge{2,2}==i).'];
end
clear i;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% mesh.edge('conn')( interface_edges , :)
  conn_edges =   mesh.edge{2,1};             % temporary assignment
  conni      =   conn_edges(itr_edges,:);    % conn matrix for interface
  mi         =   size(  conni,1);            % # of interface edge-elements
  clear   conn_edges;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Get interface nodes (2D)
% - reshape conn (edges) matrix on a concatenated vector
% - keep 'unique' nodes of concatenated vector, i.e. remove duplicates
% - interface_nodes is therefore the index of elements' conn matrix
  nodes_itr      = unique(reshape(   conni     ,1,[] ) );
%--------------------------------------------------------------------------


%plot fluid-structure interface
% if nf~=0
%     figure(1)
%     clf
%     daspect([1 1 1]);
%     hold on;
%     ve = x(conni(:,3)) - x(conni(:,1));
%     quiver(x(conni(:,1)), ve)
%     quiver(x(conni(:,2)), dot(ve,e2)*e1 - dot(ve,e1)*e2) %rotate 90 clockwise
% %     plot_edges(conni(:,[1 3]),x,ee)
%     title("Interface edges before swapping elements")
% end


% MAKE A FUNCTION
% swap nodes in the edge connectivity to satisfy the following edge element
%--------------------------------------------------------------------------

%              solid              
%      1--------2---------3  (master nodes)
%              fluid            
%--------------------------------------------------------------------------

solid_normal = zeros(1, ee, size(  conni,1),1);
%--------------------------------------------------------------------------
for e = 1:size(  conni,1) % loop over interface edges
  
    % edge e and its coordinates
    iie      =   conni(e, :);   % nodes of the current edge e
    xe       =   x(iie);        % coordinates of these nodes
    ve       =   xe(3)-xe(1);   % edge vector (of edge e)
    
    % guess fluid normal vector
    nor_f    =   -(dot(ve,e2)*e1 - dot(ve,e1)*e2);%rotate ve 90 counterclockwise
    nor_f    = nor_f/norm(nor_f);

    % get row of middle point within connf matrix
    [row,col]    =   find(connf==iie(2));
    if col<5; error('corner of the element was found'); end
    % get the oposite node within the fluid element
    if col==5 || col==6; col=col+2; else; col=col-2; end
    
    % solid normal points inwards to the fluid (for sure)
    nor_s       = x(connf(row,col))-xe(2);
    nor_s       = nor_s/norm(nor_s);
    
    % if nor_s is aligned to nor_f, then swap node positions
    if dot(nor_s,nor_f)>0; conni(e,[3 1]) = conni(e,[1 3]);   end
    
    % save normal for consistency check
    solid_normal(e) =  nor_s;
end

% apparently everything is reversed according to master element defined
% conni(:,[3 1]) = conni(:,[1 3]);

clear nor_s; clear nor_f; clear row; clear col;
%--------------------------------------------------------------------------










%% PLOT MESH




   
if multiple_solid_phase


    %--------------------------------------------------------------------------
    % plot each solid phase and fluid phase
    figure(2)
    clf
    daspect([1 1 1]);
    hold on;
    plotMesh(  1e3*mx,  connf, [0 0.5 1] );   %blue
    plotMesh(  1e3*mx,  connss{1}, [.7 .7 .7] ); %gray
    plotMesh(  1e3*mx,  connss{2}, [1 0.8 0] ); %yellow
    plotMesh(  1e3*mx,  connss{3}, 'r' );
    plotMesh(  1e3*mx,  connss{4}, [0.4 0.8 0] ); %green
    hold off;
%     plot_edges(conni(:,[1 3]),x,ee)
    % exportgraphics(gca,'plot.png','BackgroundColor','none')
    xlabel('mm');
    ylabel('mm');
    box on;
    %--------------------------------------------------------------------------

end




%% NODE SORTING OF BOUNDARIES 
disp('codeblock: NODE SORTING OF BOUNDARIES AND SPRING CONNECTION')

% NOTE: this codeblock assumes the top and bottom boundaries are solid
% boundaries.

%--------------------------------------------------------------------------
% mesh precision. Check in COMSOL for correct use of it. Here I set 1e-8 by
% try and error which is not ideal.
precision = 1e-8; 
%--------------------------------------------------------------------------


% sort nodes of the solid and fluid boundaries
%--------------------------------------------------------------------------
% nodes on the left boundary of the mesh: line x = -a/2
[  l,~]        = find(abs(  mx(:,1) +   ax/2 )<precision ); % 
  left_f       = intersect(  l,  nodes_f)';  % fluid boundary
  left_s       = intersect(  l,  nodes_s)';  % solid boundary

% nodes on the right boundary of the mesh: line x = +a/2                  
[  r,~]        = find(abs(  mx(:,1) -   ax/2 )<precision ); % 
  right_f      = intersect(  r,  nodes_f)';  % fluid boundary
  right_s      = intersect(  r,  nodes_s)';  % solid boundary

% nodes on the bottom boundary of the mesh: line y = -a/2
[  b,~]        = find(abs(  mx(:,2) +   ay/2 )<precision ); % 
   bottom_f    = intersect( b, nodes_f)'; % fluid boundary
   bottom_s    = intersect( b, nodes_s)'; % solid boundary

% nodes on the bottom boundary of the mesh: line y = +a/2
[  t,~]    = find(abs(  mx(:,2) -   ay/2 )<precision ); % 
   top_f    = intersect( t, nodes_f)'; % fluid boundary
   top_s    = intersect( t, nodes_s)'; % solid boundary

% corner nodes
  corners_s = [   intersect(  left_s ,  bottom_s) ...
                  intersect(  right_s,  bottom_s) ...
                  intersect(  right_s,  top_s   ) ...
                  intersect(  left_s ,  top_s   )];
% prescribed corners
  pcorners = [  corners_s(1)   corners_s(2)   corners_s(4)];


% exclude corners from solid node lists
  left_s   = setdiff(  left_s  ,   corners_s)';
  right_s  = setdiff(  right_s ,   corners_s)';
  top_s    = setdiff(  top_s   ,   corners_s)';
  bottom_s = setdiff(  bottom_s,   corners_s)';

% we assume there is no fluid at the corners  


clear l r b t;
%--------------------------------------------------------------------------



%% INDEX TABLE
% __________________________________________________________________________
% FLUID (F)
% SOLID (S)
%
%        CURRENT_INDEX      
%               |  
%               V
%__________________________________Reordered____________________Reordered__
% CONVERTION | mesh   | nodes_s |   -----      |    dofs_re   |   ---
%    RULE    | mesh   | nodes_f |   nodes_f    |    nodes_f   |  nodes_f
%____________|__(FS) _|_(S)_(F)_|__(S)_____(F)_|___(S)___(F)__|__(S)_(F)__
%            |    1   |   1 1   |  dofs_c   1  |  dofs_pc  1  |    p  pp
%            |    2   |   2 2   |  dofs_un  2  |  dofs_un  2  |    f  ff
%   ORDER    |    .   |   . .   |  dofs_in  .  |  dofs_in  .  |
%            |    .   |   . .   |  dofs_de  .  |           .  |
%            |    .   |   . .   |           nf |           nf | 
%            |    n   |  ns nf  |              |              |
%____________|___ ____|_________|______________|______________|____________


%%  4 INTEGRATION POINTS WITHIN MASTER ELEMENT
%--------------------------------------------------------------------------
% Notation from TensorLab
% quadratic quadrilateral element:
%         4----7----3
%         |         |
%         8         6
%         |         |
%         1----5----2
%      then CONN = [1 2 3 4 5 6 7 8].
%
% The coordinates of integration points x (in the coordinate system of the
% master element). Exact integration up to polynomial of order 5.
%          --------- 
%         | x  x  x |
%         | x  x  x |
%         | x  x  x |
%          ---------
% and the weight factors
%--------------------------------------------------------------------------
xi = [ -sqrt(3/5) *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 + sqrt(3/5) *e2
       -sqrt(3/5) *e1 + sqrt(3/5) *e2
                0 *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 +         0 *e2       
                0 *e1 + sqrt(3/5) *e2                
       -sqrt(3/5) *e1 +         0 *e2
                0 *e1 +         0 *e2   ];
w  = [ 25/81
       25/81
       25/81
       25/81
       40/81
       40/81
       40/81
       40/81
       64/81     ];
%--------------------------------------------------------------------------


%% VOLUME FRACTIONS AND POROSITY COMPUTATION

% VOLUME OF OPENING ELEMENTS (COMPUTE VOLUME FRACTION ASSUMING SOLID1)
%--------------------------------------------------------------------------
% % initialize the volume of the opening to zero
% Vo=0;
% 
% for e = 1:size(  conno,1)
%     % extract nodes
%     iie =   conno(e, :); % nodes of the current element e
%     xe  =   x(iie);            % coordinates of these nodes
% 
%     % get element area/volume and cumulate on variable Vo
%     Vo=Vo+get_element_area(xe,ee);
% end
%--------------------------------------------------------------------------


if multiple_solid_phase

% VOLUME OF SOLID ELEMENTS ( MULTIPLE PHASE MODEL)
%--------------------------------------------------------------------------

% initialize total volume of solid phase
Vs = 0;

%loop
for i =1:n_solid_phases

 % initialize area of each solid phase to zero
 V{i}=0;

    for e = 1:mss{i} % loop over all solid elements of phase i
     
    % extract nodes
    iie =   connss{i}(e, :); % nodes of the current element e
    xe  =   x(iie);            % coordinates of these nodes

    % get element area
    V{i}=V{i}+get_element_area(xe,ee);

   end % end of element loop
end

% compute total solid volume ( PU and Heavy mass only)
Vs = V{1} + V{2}+ V{3};
%--------------------------------------------------------------------------


else


% VOLUME OF SOLID ELEMENTS ( SINGLE PHASE MODEL)
%--------------------------------------------------------------------------
%   % initialize area of solid phase to zero (for porosity)
%   Vs=0;
% 
%   for e = 1:size(  conns,1) % loop over all solid elements
%     % extract nodes
%     iie =   conns(e, :); % nodes of the current element e
%     xe  =   x(iie);     % coordinates of these nodes
%         
%     % get element area
%     Vs=Vs+get_element_area(xe,ee);    
% 
%    end % end of element loop
%--------------------------------------------------------------------------


end



% VOLUME OF FLUID ELEMENTS 
%--------------------------------------------------------------------------
Vf = (ax*ay) - (Vs);
%--------------------------------------------------------------------------


% POROSITY
%--------------------------------------------------------------------------
phi = Vf/(ax*ay) ;
%--------------------------------------------------------------------------


%% ASSEMBLE STIFFNESS AND MASS MATRICES (SOLID ELEMENTS)
disp('codeblock: ASSEMBLING SOLID MASS AND STIFFNESS MATRICES')
%--------------------------------------------------------------------------
% Knn and Mnn are only temporary matrix and they are not the matrices K and
% M from discretized form in the notes. Note their size is the total number
% of nodes which include the fluid nodes. It's only relevant property is
% that it preserves the node indexing.
%--------------------------------------------------------------------------


% declare stiffness matrix tensor with initially zero-2nd-order tensors
%--------------------------------------------------------------------------
Knn = zeros(2, ee,   n,   n); 
Mnn = zeros(  n,   n);
%--------------------------------------------------------------------------





if multiple_solid_phase

% SOLID ELEMENTS ( MULTIPLE PHASE MODEL)
%--------------------------------------------------------------------------


%loop
for i =1:n_solid_phases

 % define material parameters from phase i   
 tC    =  matProp{i}.kappa*I*I + 2* matProp{i}.G*(  I4S - 1/3*  I*  I);
 rhos  =  matProp{i}.rho;


    for e = 1:mss{i} % loop over all solid elements of phase i
    
    % display computation percentage
    if mod(floor(e/mss{i}*100),10)==0
    fprintf('solid %d: assembled %.2f ', i, e/mss{i}*100); disp('%');
    end
    
    % extract nodes
    iie =   connss{i}(e, :); % nodes of the current element e
    xe  =   x(iie);            % coordinates of these nodes
        
    % compute element matrices in integration point loop
    Ke = zeros(2, ee, 8, 8); % zero matrix of 8x8 2nd order tensors
                             % in basis ee.
    Me = zeros(8,8);         % zero matrix of 8x8 scalars

    for k = 1:length(w) % loop over 9 integration points 
        
        xi1 = dot(xi(k), e1); % x1 coordinate of the current integr. point
        xi2 = dot(xi(k), e2); % x2 coordinate of the current integr. point
      
       
       %______________TEST_________________________________________________ 
       % % shape functions pressure
       %  Ne_p = [ 1/4*(1-xi1)*(1-xi2)
       %         1/4*(1+xi1)*(1-xi2)
       %         1/4*(1+xi1)*(1+xi2)
       %         1/4*(1-xi1)*(1+xi2) ];
       % 
       %  % pressure
       %  gradxiNe_p = [ -1/4*(1-xi2)*e1 - 1/4*(1-xi1)*e2
       %                1/4*(1-xi2)*e1 - 1/4*(1+xi1)*e2
       %                1/4*(1+xi2)*e1 + 1/4*(1+xi1)*e2
       %               -1/4*(1+xi2)*e1 + 1/4*(1-xi1)*e2 ];
       % 
       %  J_p = gradxiNe_p' * xe([1 2 3 4]); % pressure
        %______________TEST________________________________________________

       % column of the shape functions
        Ne = [ -1/4*(1-xi1  )*(1-xi2  )*(1+xi1+xi2)
               -1/4*(1+xi1  )*(1-xi2  )*(1-xi1+xi2)
               -1/4*(1+xi1  )*(1+xi2  )*(1-xi1-xi2)
               -1/4*(1-xi1  )*(1+xi2  )*(1+xi1-xi2)
                1/2*(1-xi1^2)*(1-xi2  )
                1/2*(1+xi1  )*(1-xi2^2)
                1/2*(1-xi1^2)*(1+xi2  )
                1/2*(1-xi1  )*(1-xi2^2)           ];
        
        % column of the gradient of the shape functions
        % with respect to the local coordinates of the mater element
        gradxiNe = [ - 1/4 * (-1 + xi2   )* (2*xi1 +   xi2) *e1 ...
                     - 1/4 * (-1 + xi1   )* (xi1   + 2*xi2) *e2
                     - 1/4 * (2*xi1 - xi2)* (-1    + xi2  ) *e1 ...
                     - 1/4 * (1 + xi1    )* (xi1   - 2*xi2) *e2
                       1/4 * (1 + xi2    )* (2*xi1 + xi2  ) *e1 ...
                     + 1/4 * (1 + xi1    )* (xi1   + 2*xi2) *e2
                       1/4 * (2*xi1 - xi2)* (1     + xi2  ) *e1 ...
                     + 1/4 * (-1 + xi1   )* (xi1 - 2*xi2  ) *e2
                                      xi1 * (-1 + xi2     ) *e1 ...
                                    + 1/2 * (-1 + xi1^2   ) *e2
                                      1/2 * ( 1 - xi2^2   ) *e1 ...
                                    - xi2 * ( 1 + xi1     ) *e2
                                    - xi1 * ( 1 + xi2     ) *e1 ...
                                    + 1/2 * ( 1 - xi1^2   ) *e2
                                      1/2 * (-1 + xi2^2   ) *e1 ...
                                    + xi2 * (-1 + xi1     ) *e2 ];
        % Jacobian
        J = gradxiNe' * xe;
       
        % column of the gradient of the shape functions
        % with respect to the global coordinates of the mesh
        gradNe = dot(inv(J), gradxiNe);
        % element matrix and right hand side
        Ke = Ke + w(k) * dot(gradNe,   tC, gradNe') *   thickout * det(J);
        Me = Me + w(k) *  Ne *  rhos*   Ne'  *   thickout * det(J);     
   
    end % end of interation point loop
   
   % assembly
    Knn(iie, iie) = Knn(iie, iie) + Ke;
    Mnn(iie, iie) = Mnn(iie, iie) + Me;


  

   end % end of element loop


% assignment actual mass and stiffness matrices
K=Knn(  nodes_s,  nodes_s);
M=Mnn(  nodes_s,  nodes_s);

% force symmetrization (correct numerical errors)
K = (K+K')/2;
M = (M+M')/2;
end

%--------------------------------------------------------------------------


else


% SOLID ELEMENTS ( ONE SOLID PHASE MODEL)
%--------------------------------------------------------------------------
 
  % define stiffness tensor
  tC   =   kappa *  I*  I + 2*  G * (  I4S - 1/3*  I*  I);

  for e = 1:size(  conns,1) % loop over all solid elements
    
    % display computation percentage
    if mod(floor(e/ms*100),10)==0
    fprintf('solid: assembled %.2f ', e/ms*100);disp('%');
    end


    iie =   conns(e, :); % nodes of the current element e
    xe  =   x(iie);     % coordinates of these nodes
        
    % compute element matrices in integration point loop
    Ke = zeros(2, ee, 8, 8); % zero matrix of 8x8 2nd order tensors
                             % in basis ee.
    Me = zeros(8,8);         % zero matrix of 8x8 scalars

    for k = 1:length(w) % loop over 9 integration points 
        
        xi1 = dot(xi(k), e1); % x1 coordinate of the current integr. point
        xi2 = dot(xi(k), e2); % x2 coordinate of the current integr. point

       % column of the shape functions
        Ne = [ -1/4*(1-xi1  )*(1-xi2  )*(1+xi1+xi2)
               -1/4*(1+xi1  )*(1-xi2  )*(1-xi1+xi2)
               -1/4*(1+xi1  )*(1+xi2  )*(1-xi1-xi2)
               -1/4*(1-xi1  )*(1+xi2  )*(1+xi1-xi2)
                1/2*(1-xi1^2)*(1-xi2  )
                1/2*(1+xi1  )*(1-xi2^2)
                1/2*(1-xi1^2)*(1+xi2  )
                1/2*(1-xi1  )*(1-xi2^2)           ];
        
        % column of the gradient of the shape functions
        % with respect to the local coordinates of the mater element
        gradxiNe = [ - 1/4 * (-1 + xi2   )* (2*xi1 +   xi2) *e1 ...
                     - 1/4 * (-1 + xi1   )* (xi1   + 2*xi2) *e2
                     - 1/4 * (2*xi1 - xi2)* (-1    + xi2  ) *e1 ...
                     - 1/4 * (1 + xi1    )* (xi1   - 2*xi2) *e2
                       1/4 * (1 + xi2    )* (2*xi1 + xi2  ) *e1 ...
                     + 1/4 * (1 + xi1    )* (xi1   + 2*xi2) *e2
                       1/4 * (2*xi1 - xi2)* (1     + xi2  ) *e1 ...
                     + 1/4 * (-1 + xi1   )* (xi1 - 2*xi2  ) *e2
                                      xi1 * (-1 + xi2     ) *e1 ...
                                    + 1/2 * (-1 + xi1^2   ) *e2
                                      1/2 * ( 1 - xi2^2   ) *e1 ...
                                    - xi2 * ( 1 + xi1     ) *e2
                                    - xi1 * ( 1 + xi2     ) *e1 ...
                                    + 1/2 * ( 1 - xi1^2   ) *e2
                                      1/2 * (-1 + xi2^2   ) *e1 ...
                                    + xi2 * (-1 + xi1     ) *e2 ];
        % Jacobian
        J = gradxiNe' * xe;
       
        % column of the gradient of the shape functions
        % with respect to the global coordinates of the mesh
        gradNe = dot(inv(J), gradxiNe);
        % element matrix and right hand side
        Ke = Ke + w(k) * dot(gradNe,   tC, gradNe') *   thickout * det(J);
        Me = Me + w(k) *  Ne *  rhos*   Ne'  *   thickout * det(J);     
        


    end % end of interation point loop
   

   % assembly
    Knn(iie, iie) = Knn(iie, iie) + Ke;
    Mnn(iie, iie) = Mnn(iie, iie) + Me;

   % get element area
   V=V+get_element_area(xe,ee);    
   

   end % end of element loop

% assignment actual mass and stiffness matrices
K=Knn(  nodes_s,  nodes_s);
M=Mnn(  nodes_s,  nodes_s);

% force symmetrization (correct numerical errors)
K = (K+K')/2;
M = (M+M')/2;
%--------------------------------------------------------------------------


end

%--------------------------------------------------------------------------
clear Me; clear Ke;
clear Knn; clear Mnn;
%--------------------------------------------------------------------------


%% ASSEMBLE STIFFNESS AND MASS MATRICES IN A FLUID ELEMENT LOOP 
disp('codeblock: ASSEMBLING FLUID MASS AND STIFFNESS MATRICES')
%--------------------------------------------------------------------------
% following the notation in my notes
% A - [1/Pa]      equivalent "mass matrix", term that multiplies d^2/dt^2 p 
% B - [1/(kg/m3)] equivalent "stiffness matrix", term that multiplies  p 
% d - [m/s2]      RHS term "external forces"
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
Ann = zeros(  n,   n); % zero matrix of nxn scalars 
Bnn = zeros(  n,   n); % zero matrix of nxn scalars


   for e = 1:size(  connf,1)% loop over all fluid elements
    
    % display computation percentage
    if mod(floor(e/mf*100),10)==0
    fprintf('fluid: assembled %.2f ', e/mf*100); disp('%');
    end

    iie =   connf(e, :); % nodes of the current element e
    xe  =   x(iie);     % coordinates of these nodes
        
    % compute element matrices in integration point loop       
    Ae = zeros(8, 8);        % zero matrix of 8x8 scalars     
    Be = zeros(8, 8);        % zero matrix of 8x8 scalars

    for k = 1:length(w) % loop over 9 integration points 
        
        xi1 = dot(xi(k), e1); % x1 coordinate of the current integr. point
        xi2 = dot(xi(k), e2); % x2 coordinate of the current integr. point

       % column of the shape functions
        Ne = [ -1/4*(1-xi1  )*(1-xi2  )*(1+xi1+xi2)
               -1/4*(1+xi1  )*(1-xi2  )*(1-xi1+xi2)
               -1/4*(1+xi1  )*(1+xi2  )*(1-xi1-xi2)
               -1/4*(1-xi1  )*(1+xi2  )*(1+xi1-xi2)
                1/2*(1-xi1^2)*(1-xi2  )
                1/2*(1+xi1  )*(1-xi2^2)
                1/2*(1-xi1^2)*(1+xi2  )
                1/2*(1-xi1  )*(1-xi2^2)           ];
        
        % column of the gradient of the shape functions
        % with respect to the local coordinates of the mater element
        gradxiNe = [ - 1/4 * (-1 + xi2   )* (2*xi1 +   xi2) *e1 ...
                     - 1/4 * (-1 + xi1   )* (xi1   + 2*xi2) *e2
                     - 1/4 * (2*xi1 - xi2)* (-1    + xi2  ) *e1 ...
                     - 1/4 * (1 + xi1    )* (xi1   - 2*xi2) *e2
                       1/4 * (1 + xi2    )* (2*xi1 + xi2  ) *e1 ...
                     + 1/4 * (1 + xi1    )* (xi1   + 2*xi2) *e2
                       1/4 * (2*xi1 - xi2)* (1     + xi2  ) *e1 ...
                     + 1/4 * (-1 + xi1   )* (xi1 - 2*xi2  ) *e2
                                      xi1 * (-1 + xi2     ) *e1 ...
                                    + 1/2 * (-1 + xi1^2   ) *e2
                                      1/2 * ( 1 - xi2^2   ) *e1 ...
                                    - xi2 * ( 1 + xi1     ) *e2
                                    - xi1 * ( 1 + xi2     ) *e1 ...
                                    + 1/2 * ( 1 - xi1^2   ) *e2
                                      1/2 * (-1 + xi2^2   ) *e1 ...
                                    + xi2 * (-1 + xi1     ) *e2 ];
        % Jacobian
        J = gradxiNe' * xe;
       
        % column of the gradient of the shape functions
        % with respect to the global coordinates of the mesh
        gradNe = dot(inv(J), gradxiNe);
        % element matrix and right hand side
        Ae = Ae + w(k)*(1/  rhof)*Ne*(1/  c^2)*Ne'*  thickout * det(J);
        Be = Be + w(k)*(1/  rhof)*dot(gradNe,gradNe')*  thickout*det(J);     
   
    end % end of interation point loop
   
    % assembly
    Ann(iie, iie) = Ann(iie, iie) + Ae;
    Bnn(iie, iie) = Bnn(iie, iie) + Be;
   end % end of fluid element loop

% define actual stiffness and mass matrices.
A=Ann(  nodes_f,  nodes_f);
B=Bnn(  nodes_f,  nodes_f);

% force symmetrization (correct numerical errors)
A = (A+A')/2;
B = (B+B')/2;

%--------------------------------------------------------------------------
clear Ae; clear Be;
clear Bnn; clear Ann;
%--------------------------------------------------------------------------


%% INTEGRATION POINTS WITHIN MASTER EDGE (ELEMENT)
%--------------------------------------------------------------------------
% use here the same shape functions for both fluid and solid edges
% the edges have three nodes refered as 1, 2 and 3 below.
%
%              solid              
%      1--------2---------3  (master nodes)
%       --x-----x-----x--    (integration points)

%              fluid             | (solid normal)   
%                                V 
% then CONN = [1 3 2].
% coordinates of integration points sketched as x.  Exact integration up
% to polynomial of order 5.
%
% in the coordinate system of the master element and the weight factors
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
xi = [ -sqrt(3/5)
              0
        sqrt(3/5) ];
w  = [ 5/9
       8/9
       5/9        ];
%--------------------------------------------------------------------------

%% ASSEMBLE COUPLING MATRIX, LOOP ON INTERFACE EDGES 
disp('codeblock: ASSEMBLING COUPLING MATRIX')
%--------------------------------------------------------------------------
% use here the same shape functions for both fluid and solid elements
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
Cnn = zeros(1, ee,   n,   n);        % zero matrix of nxn vectors


   for e = 1:size(  conni,1) % loop over interface edges
    % display computation percentage
    if mod(floor(e/mi*100),10)==0
    fprintf('coupling: assembled %.2f ', e/mi*100); disp('%');
    end


    % edge e and its coordinates
    iie      =   conni(e, :); % nodes of the current edge e
    xe       =   x(iie);      % coordinates of these nodes
    
    % solid normal vector
    ve       = xe(3)-xe(1); % edge vector (of edge e)
    normal_s = ( dot(ve,e2)*e1 - dot(ve,e1)*e2)/...
                sqrt( dot(ve,e1)^2+dot(ve,e2)^2 ); % solid normal

%   normal_f = (-dot(ve,e2)*e1 + dot(ve,e1)*e2)/norm(ve); % fluid normal

    % compute element matrices in integration point loop
    Ce =  zeros(1, ee, 3, 3);        % 3x3 zero matrix of vectors

    for k = 1:length(w) % loop over 3 integration points 
        
       % xi(k) is the coordinate of the current integr. point
       % column of the shape functions
        Ne = [ -1/2* (1-xi(k))*    xi(k)
                     (1-xi(k))* (1+xi(k))
                1/2*    xi(k) * (1+xi(k))  ];

       % elementary length ds=|dx|d(xi)       =  1/2*norm(ve)*d(xi) 
       % elementary area   dA=|dx|d(xi) thick =  1/2*norm(ve)*thick*d(xi)
       
    
       % element matrix 
        Ce = Ce + w(k) * (Ne * Ne') * normal_s * 1/2*norm(ve)*   thickout ;
        
         
   
    end % end of interation point loop
   
% assembly
    Cnn(iie, iie) = Cnn(iie, iie) + Ce;

   end % end of interface edge loop

C = Cnn(  nodes_s,  nodes_f);
clear Cnn;
%--------------------------------------------------------------------------

 
%% INDEX OF DISPLACEMENT/PRESSURE DOF UNDER NODES_S/NODES_F ORDERING
%--------------------------------------------------------------------------
% The matrices Mnn, Knn, Ann, Bnn, and Cnn have indices that correspond to
% the node number, i.g., line 1 are the terms related to node 1.
% The matrices M, K, A, B, and C, however, have indices that do not 
% correspond to the node number. Their indices correspond to the
% displacements degree of freedom or pressure degrees of freedom.
% Suppose, for example, you want to find the index in matrix M that
% correspond to the bottom-left corner displacement dof. This corner node
% is stored in the vector [corner], say first entry, corner(1). The index
% we are looking for, sometimes refered as 'new' index, is found using the
% following command:
% 
% find(nodes_s==corner(1)).
%
% get_ind() was implemented for this operation
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% fluid
    % bottom nodes
    doff_b = get_ind(  nodes_f,   bottom_f ); 
    % top nodes
    doff_t = get_ind(  nodes_f,   top_f    ); 
    % left nodes
    doff_l = get_ind(  nodes_f,   left_f   ); 
    % right nodes
    doff_r = get_ind(  nodes_f,   right_f  ); 
    % interface nodes
    doff_i = get_ind(  nodes_f,   nodes_itr);  % error here when triangle
                                               % elements on fluid phase
% solid 
    % bottom nodes
    dofs_b = get_ind(  nodes_s,   bottom_s ); 
    % top nodes
    dofs_t = get_ind(  nodes_s,   top_s    ); 
    % left nodes
    dofs_l = get_ind(  nodes_s,   left_s   ); 
    % right nodes
    dofs_r = get_ind(  nodes_s,   right_s  ); 
    % interface nodes
    dofs_i = get_ind(  nodes_s,   nodes_itr);
    % corner nodes
    dofs_c = get_ind(  nodes_s,   corners_s);
    % prescribed corners
    dofs_pc = get_ind(  nodes_s,   pcorners);
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------
% TODO:afterwards, check which sets (above) of dof were actually used.
%--------------------------------------------------------------------------

%__________________________________________________________________________
% INDEX TABLE
% FLUID (F)
% SOLID (S)
%
%                   CURRENT_INDEX      
%                          |  
%                          V
%__________________________________Reordered____________________Reordered__
% CONVERTION | mesh   | nodes_s |   -----      |    dofs_re   |   ---
%    RULE    | mesh   | nodes_f |   nodes_f    |    nodes_f   |  nodes_f
%____________|__(FS) _|_(S)_(F)_|__(S)_____(F)_|___(S)___(F)__|__(S)_(F)__
%            |    1   |   1 1   |  dofs_c   1  |  dofs_pc  1  |    p  pp
%            |    2   |   2 2   |  dofs_un  2  |  dofs_un  2  |    f  ff
%   ORDER    |    .   |   . .   |  dofs_in  .  |  dofs_in  .  |
%            |    .   |   . .   |  dofs_de  .  |           .  |
%            |    .   |   . .   |           nf |           nf | 
%            |    n   |  ns nf  |              |              |
%____________|___ ____|_________|______________|______________|____________

%% ELIMINATING DEPENDENT NODES OF SOLID PHASE


%--------------------------------------------------------------------------
% constrain relation of periodicity 
% u(top)    = u(bottom) +  (u(coners(4))-u(coners(1))) ones()
% u(right)  = u(left)   +  (u(coners(2))-u(coners(1))) ones()
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    % dependent dofs of periodic bc
    dofs_de = [dofs_t dofs_r];    % [top right], excluding corners
    
    % independent nodes of periodic bc
    dofs_in = [dofs_b dofs_l];    % [bottom left], excluding corners
    

    % safety check
    if(length(dofs_in)~= length(dofs_de)) 
        error("The periodic bc must connect same number of nodes");
    end

% unconstrained nodes ( nodes not at the boundary)
dofs_un = setdiff(1:  ns, [dofs_in dofs_de dofs_c]); % exclude corners                           
% include corners as first elements to follow the notes notation, ideally
% we should chenge the variable name here.  
dofs_un = [dofs_c dofs_un];           
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% reordering mass, stiffness, and coupling matrix components
% mass tensor matrix
% M=  Muu Mui Mud   
%     Miu Mii Mid  
%     Mdu Mdi Mdd  
% transform to tensor matrix
M  =   I*M;
% split
Muu=M(dofs_un, dofs_un); Mui=M(dofs_un, dofs_in); Mud=M(dofs_un, dofs_de);  
Miu=M(dofs_in, dofs_un); Mii=M(dofs_in, dofs_in); Mid=M(dofs_in, dofs_de);  
Mdu=M(dofs_de, dofs_un); Mdi=M(dofs_de, dofs_in); Mdd=M(dofs_de, dofs_de);

% stiffness tensor matrix
% K=  Kuu Kui Kud   
%     Kiu Kii Kid  
%     Kdu Kdi Kdd  
Kuu=K(dofs_un, dofs_un); Kui=K(dofs_un, dofs_in); Kud=K(dofs_un, dofs_de);  
Kiu=K(dofs_in, dofs_un); Kii=K(dofs_in, dofs_in); Kid=K(dofs_in, dofs_de);  
Kdu=K(dofs_de, dofs_un); Kdi=K(dofs_de, dofs_in); Kdd=K(dofs_de, dofs_de);  

% coupling tensor matrix
% C=  Cu    
%     Ci   
%     Cd  
Cu=C(dofs_un,:);
Ci=C(dofs_in,:);
Cd=C(dofs_de,:);
%--------------------------------------------------------------------------

%__________________________________________________________________________
% INDEX TABLE
% FLUID (F)
% SOLID (S)
%
%                               CURRENT_INDEX      
%                                      |  
%                                      V
%__________________________________Reordered____________________Reordered__
% CONVERTION | mesh   | nodes_s |   -----      |    dofs_re   |   ---
%    RULE    | mesh   | nodes_f |   nodes_f    |    nodes_f   |  nodes_f
%____________|__(FS) _|_(S)_(F)_|__(S)_____(F)_|___(S)___(F)__|__(S)_(F)__
%            |    1   |   1 1   |  dofs_c   1  |  dofs_pc  1  |    p  pp
%            |    2   |   2 2   |  dofs_un  2  |  dofs_un  2  |    f  ff
%   ORDER    |    .   |   . .   |  dofs_in  .  |  dofs_in  .  |
%            |    .   |   . .   |  dofs_de  .  |           .  |
%            |    .   |   . .   |           nf |           nf | 
%            |    n   |  ns nf  |              |              |
%____________|___ ____|_________|______________|______________|____________


%--------------------------------------------------------------------------
% exclude corners again to stick with notes notation, they will be treated
% separately by the projection matrix
dofs_un([1 2 3 4])=[];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% projection matrix submatrices 
O_1u = zeros(1,length(dofs_un));
O_u1 = O_1u';

O_1B = zeros(1,length(dofs_b));
O_B1 = O_1B';

O_1L = zeros(1,length(dofs_l));
O_L1 = O_1L';

O_uB = zeros(length(dofs_un),length(dofs_b));
O_Bu = O_uB';

O_uL = zeros(length(dofs_un),length(dofs_l));
O_Lu = O_uL';

O_BL = zeros(length(dofs_b),length(dofs_l));
O_LB = O_BL';

I_uu = eye(length(dofs_un),length(dofs_un));
I_BB = eye(length(dofs_b),length(dofs_b));
I_LL = eye(length(dofs_l),length(dofs_l));

I_B1 = ones(length(dofs_b),1);
I_L1 = ones(length(dofs_l),1);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% projection matrix
%     u_1   = T * u_1
%     u_2         u_2
%     u_3         u_4
%     u_4         u_un
%     u_un        u_in
%     u_in
%     u_de

% projection matrix
T  = [  1       0      0      O_1u    O_1B    O_1L
        0       1      0      O_1u    O_1B    O_1L
       -1       1      1      O_1u    O_1B    O_1L
        0       0      1      O_1u    O_1B    O_1L
       O_u1    O_u1   O_u1    I_uu    O_uB    O_uL
       O_B1    O_B1   O_B1    O_Bu    I_BB    O_BL
       O_L1    O_L1   O_L1    O_Lu    O_LB    I_LL
      -I_B1    O_B1   I_B1    O_Bu    I_BB    O_BL
      -I_L1    I_L1   O_L1    O_Lu    O_LB    I_LL];
% clear matrix components
clear O_1B O_1L O_1u O_B1 O_BL O_Bu  ...
      O_L1 O_LB O_Lu O_u1 O_uB O_uL I_B1 I_BB I_L1 I_LL I_uu
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% reduced mass tensor matrix (dependent nodes eliminated)
M_ =  T.'*[Muu Mui Mud;
           Miu Mii Mid;
           Mdu Mdi Mdd]*T;

% reduced stiffness tensor matrix (dependent nodes eliminated)
K_ =  T.'*[Kuu Kui Kud;
           Kiu Kii Kid;
           Kdu Kdi Kdd]*T;

% reduced mass tensor matrix (dependent nodes eliminated)
C_ =  T.'*[Cu;
          Ci;
          Cd];
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% set of remaining nodes, this is the vector that contains the rules for
% changing the indices from 'before projection' to 'after projection'.
         
dofs_re = [dofs_c(1)  dofs_c(2)  dofs_c(4)  dofs_un dofs_in];
ns_re   = length(dofs_re);
%--------------------------------------------------------------------------

%__________________________________________________________________________
% INDEX TABLE
% FLUID (F)
% SOLID (S)
%
%                                               CURRENT_INDEX      
%                                                      |  
%                                                      V
%__________________________________Reordered____________________Reordered__
% CONVERTION | mesh   | nodes_s |   -----      |    dofs_re   |   ---
%    RULE    | mesh   | nodes_f |   nodes_f    |    nodes_f   |  nodes_f
%____________|__(FS) _|_(S)_(F)_|__(S)_____(F)_|___(S)___(F)__|__(S)_(F)__
%            |    1   |   1 1   |  dofs_c   1  |  dofs_pc  1  |    p  pp
%            |    2   |   2 2   |  dofs_un  2  |  dofs_un  2  |    f  ff
%   ORDER    |    .   |   . .   |  dofs_in  .  |  dofs_in  .  |
%            |    .   |   . .   |  dofs_de  .  |           .  |
%            |    .   |   . .   |           nf |           nf | 
%            |    n   |  ns nf  |              |              |
%____________|___ ____|_________|______________|______________|____________



%% PRESCRIBED AND FREE NODES SPLIT
%--------------------------------------------------------------------------
% INDEX OF DISPLACEMENT DOF IN NODES_S (ELIMINATED DEPENDENT DOFS)
% same index change procedure must be executed but only in the solid phase.

%|__NOTES_NOTATION_|____ DESCRIPTION_______________________________________
%|                 |
%|        p        | prescribed nodes in the solid phase
%|_________________|_______________________________________________________
%|                 |
%|        f        | free nodes in the solid phase                                                      
%|                 |
%|_________________|_______________________________________________________
%|                 |
%|        p'       | prescribed nodes in the fluid phase
%|_________________|_______________________________________________________
%|                 |
%|        f'       | free nodes in the fluid phase
%|_________________|_______________________________________________________
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% solid
p  = get_ind(dofs_re,dofs_pc); % = [1 2 3] by construction
f  = setdiff(1:ns_re,p);

% fluid
pp = [doff_l doff_r doff_b doff_t];
fp = setdiff(1:nf,pp);
%--------------------------------------------------------------------------


%__________________________________________________________________________
% INDEX TABLE
% FLUID (F)
% SOLID (S)
%
%                                                            CURRENT_INDEX      
%                                                                   |  
%                                                                   V
%__________________________________Reordered____________________Reordered__
% CONVERTION | mesh   | nodes_s |   -----      |    dofs_re   |   ---
%    RULE    | mesh   | nodes_f |   nodes_f    |    nodes_f   |  nodes_f
%____________|__(FS) _|_(S)_(F)_|__(S)_____(F)_|___(S)___(F)__|__(S)_(F)__
%            |    1   |   1 1   |  dofs_c   1  |  dofs_pc  1  |    p  pp
%            |    2   |   2 2   |  dofs_un  2  |  dofs_un  2  |    f  ff
%   ORDER    |    .   |   . .   |  dofs_in  .  |  dofs_in  .  |
%            |    .   |   . .   |  dofs_de  .  |           .  |
%            |    .   |   . .   |           nf |           nf | 
%            |    n   |  ns nf  |              |              |
%____________|___ ____|_________|______________|______________|____________

%% PARTITIONING

%--------------------------------------------------------------------------
% solid
    % mass tensor matrix
    % M =  Mpp Mpf   
    %      Mfp Mff 
    M_p_p   = M_(p , p ); M_p_f   = M_(p , f );
    M_f_p   = M_(f , p ); M_f_f   = M_(f , f );
    % stiffness tensor matrix
    % K =  Kpp Kpf   
    %      Kfp Kff 
    K_p_p   = K_(p , p ); K_p_f   = K_(p , f );
    K_f_p   = K_(f , p ); K_f_f   = K_(f , f );
    
% fluid
    % 'mass' scalar matrix 
    % B =  Bp'p' Bp'f'   
    %      Bf'p' Bf'f' 
    A_pp_pp = A (pp, pp); A_pp_fp = A (pp, fp);
    A_fp_pp = A (fp, pp); A_fp_fp = A (fp, fp);

    % 'stiffness' scalar matrix 
    % B =  Bp'p' Bp'f'   
    %      Bf'p' Bf'f' 
    B_pp_pp = B (pp, pp); B_pp_fp = B (pp, fp);
    B_fp_pp = B (fp, pp); B_fp_fp = B (fp, fp);

% coupling tensor matrix ( converts fluid pressure to solid displacement)
    % C    =  Cpp'  Cpf'   
    %         Cfp'  Cff' 
    C_p_pp  = C_(p , pp);C_p_fp  = C_(p , fp);
    C_f_pp  = C_(f , pp);C_f_fp  = C_(f , fp);
% coupling tensor matrix ( converts solid displacement to fluid pressure)
    % -C^T =  Cp'p  Cp'f = -(Cpp')^T  -(Cfp')^T  
    %         Cf'p  Cf'f   -(Cpf')^T  -(Cff')^T

    %     C_pp_p = - C_p_pp.' ; C_pp_f = - C_f_pp.' ;
    %     C_fp_p = - C_p_fp.' ; C_fp_f = - C_f_fp.' ;


%--------------------------------------------------------------------------



%% TRANFORM TENSOR MATRIX TO EQUIVALENT SCALAR MATRIX
disp('codeblock: TRANFORM TENSOR MATRIX TO EQUIVALENT SCALAR MATRIX')

%--------------------------------------------------------------------------

% mass matrix M
sM_p_p = tm2sm_v( M_p_p, ee , 2 ); sM_p_f = tm2sm_v( M_p_f, ee , 2 );
sM_f_p = tm2sm_v( M_f_p, ee , 2 ); sM_f_f = tm2sm_v( M_f_f, ee , 2 );

% stiffness matrix K
sK_p_p = tm2sm_v( K_p_p, ee , 2 ); sK_p_f = tm2sm_v( K_p_f, ee , 2 );
sK_f_p = tm2sm_v( K_f_p, ee , 2 ); sK_f_f = tm2sm_v( K_f_f, ee , 2 );

% coupling matrix C
sC_p_pp = tm2sm_v(C_p_pp, ee, 1); sC_p_fp = tm2sm_v(C_p_fp, ee, 1);
sC_f_pp = tm2sm_v(C_f_pp, ee, 1); sC_f_fp = tm2sm_v(C_f_fp, ee, 1);


% zero matrices declaration directly as equivalent scalar matrix
O_pp_p = sparse(zeros( length(pp) , size(ee,1)*length(p)));
O_fp_p = sparse(zeros( length(fp) , size(ee,1)*length(p)));
O_fp_f = sparse(zeros( length(fp) , size(ee,1)*length(f)));
O_pp_f = sparse(zeros( length(pp) , size(ee,1)*length(f)));

O_p_pp = O_pp_p.';
O_f_fp = O_fp_f.';
O_p_fp = O_fp_p.';
O_f_pp = O_pp_f.';
%--------------------------------------------------------------------------


%% ASSEMBLE NONSYMMETRIC SYSTEM OF EQUATIONS 
disp('codeblock: SOLVE')

% recap on -C^T relations
%--------------------------------------------------------------------------
%     C_pp_p = - C_p_pp.' ; C_pp_f = - C_f_pp.' ;
%     C_fp_p = - C_p_fp.' ; C_fp_f = - C_f_fp.' ;
%--------------------------------------------------------------------------


% lamda submatrices
%--------------------------------------------------------------------------
lam_P_P = [  sM_p_p          O_p_pp;
            -(sC_p_pp).'      A_pp_pp ];

lam_P_F = [  sM_p_f          O_p_fp;
           -(sC_f_pp).'      A_pp_fp ];

lam_F_P = [  sM_f_p          O_f_pp;
           -(sC_p_fp).'      A_fp_pp ];

lam_F_F = [  sM_f_f          O_f_fp;
           -(sC_f_fp).'      A_fp_fp ];
%--------------------------------------------------------------------------


% mu submatrices
%--------------------------------------------------------------------------
mu_P_P  = [  sK_p_p              sC_p_pp ;
            O_pp_p              B_pp_pp    ];  

mu_P_F  = [  sK_p_f              sC_p_fp ;
            O_pp_f              B_pp_fp    ];  

mu_F_P  = [  sK_f_p              sC_f_pp ;
            O_fp_p              B_fp_pp    ];  

mu_F_F  = [  sK_f_f              sC_f_fp ;
            O_fp_f              B_fp_fp    ];  
%--------------------------------------------------------------------------

%% EIGEN: EIGEN STUDY, PRINT EIGENFREQUENCIES AND EIGENMODES

%--------------------------------------------------------------------------
% define desired number of modes to compute/display
n_modes = 5;
%-------------------------------------------------------------------------

% iterative solver
%--------------------------------------------------------------------------
disp('ITERATIVE SOLVER')  
%     warning("The iterative solver doesn't converge for RVEs which" + ...
%             "the solid phase is connected only with springs. It will" + ...
%             "likely cause the non physical eigenfrenquencies with" + ...
%             "imaginary part");
    % Subset of eigenvalues -> CHEAPER %%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % eigs is recommended by matlab when sparse nonsymmetric matrices.
    % % compute eigen modes V, and eigen values D with smallest absolute
    % value right problem.
    [phi_F_Q, Dr] = eigs(mu_F_F ,lam_F_F ,n_modes,...
                   'smallestabs',...
                   'Display',0,...
                   'Tolerance',1e-15,...
                   'FailureTreatment','replacenan');
    % left problem (psi appears always transposed along the expressions)
    [psi_F_Q, Dl] = eigs(mu_F_F',lam_F_F',n_modes,...
                   'smallestabs',...
                   'Display',0,...
                   'FailureTreatment','replacenan');
    
    % CHECK whether the frequencies in Dr and Dl are enoughly equal. For
    % the current example they are equal up to the 5th digit after comma.
    

    % print frequencies f
    fprintf('Eigen frequencies:\n');
    freqs = diag(sqrt(Dr)/(2*pi));
    fprintf('%f Hz\n', freqs);
%--------------------------------------------------------------------------



%direct solver
% %--------------------------------------------------------------------------
% %disp('DIRECT SOLVER')
% 
% % compute F eigenfrequencies
% [phi_F_Q, Dr, psi_F_Q] = eig(full(mu_F_F) ,full(lam_F_F),'qz');
% 
% % sort(smallest) and print frequencies
% [freqs,ir] = sort(diag(sqrt(Dr)/(2*pi)));
% fprintf('Eigen frequencies:\n');
% fprintf('%f Hz\n', freqs(1:n_modes));
% 
% % reorder matrices from smallest eigenfrequency
% Dr = Dr(ir,ir);
% phi_F_Q = phi_F_Q(:,ir);
% psi_F_Q = psi_F_Q(:,ir);
% 
% % disregard larger eigenfrequencies
% Dr = Dr(1:n_modes,1:n_modes);
% phi_F_Q = phi_F_Q(:,1:n_modes);
% psi_F_Q = psi_F_Q(:,1:n_modes);
% %--------------------------------------------------------------------------


% EIGEN:  PHASE CORRECTION (not necessary when direct solver)
%--------------------------------------------------------------------------
% IMPORTANT NOTE: PHASE CORRECTION BEFORE NORMALIZATION 
% Compute round(psi'*lamFF*phi). It should be positive diag matrix - note 
% that this product should be approx positive diag because the the
% frequencies must be positive-. Since the right and left eigen modes were
% computed with different function calls from the function eigs, they can
% potentially converge to a shape mode with phase difference of pi. In 
% this case, psi'*lamFF*phi ( or round(psi'*lamFF*phi) for getting rid 
% of almost zero off diag elements) may have some of the diagonal elements
% as -negative instead of +positive. To correct this effect it is necessary
% to multiply either phi or psi by the 'sign' matrix 'I'= psi'*lamFF*phi.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% get diagonal
diag_psi_lam_phi   = diag(psi_F_Q'*lam_F_F*phi_F_Q);
% make it a vector with only the sign of each element
diag_psi_lam_phi   = diag_psi_lam_phi./abs(diag_psi_lam_phi);
% make 'identity' matrix with sign correction
I_phase_correction = diag(diag_psi_lam_phi);

% correct the sign of left eigen modes with respect to the right modes
psi_F_Q = psi_F_Q * I_phase_correction;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% save version of phi for plotting
phi_plot = phi_F_Q;
%--------------------------------------------------------------------------


% EIGEN:  MODE NORMALIZATION WITH RESPESCT TO THE LAMDA MATRIX
%--------------------------------------------------------------------------
% ready to mass normalization, now, the product psi'*lamFF*phi has only
% postive diagonal compared to previous case
% get norm of each mode

vec_norms = sqrt(diag(psi_F_Q'*lam_F_F*phi_F_Q));
% figure(8);subplot(2,1,1); loglog(abs(vec_norms)); hold on;

phi_F_Q=phi_F_Q./vec_norms'; %divide each column of phi by element of vec_norms 
psi_F_Q=psi_F_Q./vec_norms'; %divide each column of phi by element of vec_norms 

% subplot(2,1,2); loglog(abs(sqrt(diag(psi_F_Q'*lam_F_F*phi_F_Q))));
% CONSISTENCY CHECK: 
% Compute again psi'*lamFF*phi. It should be approx. the identity matrix.
% when eigs is used the phase problem doen't exist. It can be seen by
% computing psi_F_Q(:,1:n_modes)'*lam_F_F*phi_F_Q(:,1:n_modes).
% The diagonal is made out of positive integers.
%--------------------------------------------------------------------------





% EIGEN: ASSIGN DISP. SHAPE MODES TO TENSOR VECTOR & PLOT
%--------------------------------------------------------------------------
% displacement 'ured' is followed is 'u' displacement reduced 'red' thanks
% to the elimitation of the dependend points due to periodic bc
% 'ured' stands for 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% prescribed (reduced) displacements are zero
    % declare & assign
ured_p = zeros(1, ee, length(p), size(phi_F_Q,2) ); 
        
% NOT TAKE REAL PART free (reduced) displacements are the eigenvectors    
ured_f = phi_plot(1:2:2*length(f)-1,:)*e1+ phi_plot(2:2:2*length(f),:)*e2;

% (use the following lines when complex displacements appear)
% TAKE REAL PART free (reduced) displacements are the eigenvectors    
% ured_f = real(phi_plot(1:2:2*length(f)-1,:))*e1 + ...
%          real(phi_plot(2:2:2*length(f),:))*e2;


% total reduced displacements (indices ordered as [dofs_pc dofs_un dofs_in])
    % declare ured(unconstrained+corners+independent nodes, modes)
    ured = zeros(1, ee, length(p)+length(f), size(phi_F_Q,2));
    % assign
    % assign prescribed nodes
    ured(p,:) = ured_p;
    % assign free nodes
    ured(f,:) = ured_f;
    % free unused variables
    clear ured_f; clear uref_p;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% declare displacement tensor vectors (total # of dofs)
vu     = zeros(1, ee,   ns, size(phi_F_Q,2) );
vu_aux = zeros(1, ee,   ns, size(phi_F_Q,2) ); 
% assign. Note vu_aux is in order [dofs_c dofs_un dofs_in dofs_de]
vu_aux = T*ured;
% assign to vu. Note vu is in nodes_s order
vu([dofs_c dofs_un dofs_in dofs_de],:)  = vu_aux;
% free unused variables
clear vu_aux;

% prepare for plot: length(vu_plot)>length(vu)
    % declare disp. vector for plot -  whole mesh
    vu_plot = zeros(1, ee,   ns, size(phi_F_Q,2) );
    % assign vu to vector for plot, only solid nodes are non-zero
    vu_plot(  nodes_s,:) = vu;
%--------------------------------------------------------------------------


% ASSIGN PRESSURE SHAPE MODES TO SCALAR VECTOR vp
%--------------------------------------------------------------------------
% extract pressure values from bc and eigen vectors
    % prescribed pressures are zero (same for every mode)    
    p_pp     = zeros(length(pp), size(phi_F_Q,2)); 
    
    % free pressures (different for each mode)
    p_fp     = phi_plot(end-length(fp)+1:end,:);

% declare pressure scalar vector
vp_eig = zeros(  nf,size(phi_F_Q,2));

% assign to scalar vector solution of size
vp_eig(pp,:)=p_pp;
vp_eig(fp,:)=p_fp;

% prepare for plot
    % declare pressure scalar vector for plot
    vp_plot = zeros(  n,size(phi_F_Q,2));
    
    % assign vp to vector for plot
    vp_plot(  nodes_f,:) = vp_eig;
%--------------------------------------------------------------------------




% PLOT EIGEN PRESSURE/DISP MODE SHAPES (RESCALED DISP.)
%--------------------------------------------------------------------------
disp('codeblock: PLOT')

warning("variable vp is overriden in in the next codeblock");

sfreqs=freqs; %sorted frequencies vector
% mode=1;


for mode=1:n_modes

if nf==0
% plot the  displacement
%--------------------------------------------------------------------------
figure
clf
daspect([1 1 1]);
sign=1; % control modes shape phase ( 0 or 180 degrees)
axis equal
hold on
vu_plot_n=sign*vu_plot(:,mode)/norm(norm(vu(:,mode)))*a; %rescale to see
% vu_plot_n=sign*vu_plot(:,mode);
femplot(  x,   connss{1}, 'Color', 'c');
femplot(  x+vu_plot_n,   conns, 'Color', 'k');
title("Mode shape of frequency "  + num2str(sfreqs(mode)) + " Hz" );
hold off
%--------------------------------------------------------------------------
else
% plot the pressure solution and displacement
%--------------------------------------------------------------------------
figure
clf
daspect([1 1 1]);
sign=1; % control modes shape phase ( 0 or 180 degrees in time)
scatter(dot( x(nodes_f),e1),dot(x(nodes_f),e2),20,sign*vp_eig(:,mode),'filled')
colorbar;
pmax = max(abs(vp_plot(:,mode)));
caxis([-pmax pmax]);
colormap('jet');
axis equal
hold on
vu_plot_n=sign*vu_plot(:,mode)/norm(norm(vu(:,mode)))*ax; %rescale to see
femplot(  x+vu_plot_n,   connss{1}, 'Color', 'k');
femplot(  x+vu_plot_n,   connss{2}, 'Color', 'k');
femplot(  x+vu_plot_n,   connss{3}, 'Color', 'k');
 % plotMesh(  [dot(x+vu_plot_n,e1) dot(x+vu_plot_n,e2)],  connss{1}, [.7 .7 .7] ); %gray
 % plotMesh(  [dot(x+vu_plot_n,e1) dot(x+vu_plot_n,e2)],  connss{2}, [1 0.8 0] ); %yellow
 % plotMesh(  [dot(x+vu_plot_n,e1) dot(x+vu_plot_n,e2)],  connss{3}, 'r' ); %gray

title("Mode shape "+ num2str(mode) + " of frequency "  + num2str(sfreqs(mode)) + " Hz" );
hold off
%--------------------------------------------------------------------------
end


end


% plot for review
%--------------------------------------------------------------------------
% for mode=1:n_modes
% figure
% plot(dot( x(top_f),e1)/ax,sign*vp_plot(top_f,mode))
% sign=1; % control modes shape phase ( 0 or 180 degrees in time)
% pmax = max(abs(vp_plot(:,mode)));
% caxis([-pmax pmax]);
% axis equal
% end
% grid on;
%--------------------------------------------------------------------------


%% STEADY: PRESCRIBED DISPLACEMENT, PRESSURE AND ITS GRADIENTS



   % reference position vector 
   warning('xR=0');
   %--------------------------------------------------------------------------
    % node to be fixed displacement
    x1u = 0*e1 ; % RVE bottom-left corner
    % node to be imposed pressure 
        x1p = 0*e1;     % some node on the bottom RVE boundary
    %--------------------------------------------------------------------------


    % ones
    %--------------------------------------------------------------------------
    ones_p = ones(length(p) ,1);
    ones_pp = ones(length(pp),1);
    %--------------------------------------------------------------------------


    % prescribed solid points are [corners_s]
    %--------------------------------------------------------------------------
        % MODEL INPUT
    %     uM     = 1e-5*(e1+e2);
        uM     = 0.00*e1;
        graduM = +0.0 *e1*e1 + 0.0   *e1*e2 +...
                  0.0 *e2*e1 + 0.00  *e2*e2;
        % first-order comp. homogenization, prescribed disp. at the corners
        % note input index of x() is always "mesh" type.
        % delta xs_p
        dxs_p = x(  pcorners)- x1u*ones_p;
        u_p    =  uM *ones_p + dot(graduM.', dxs_p);


        % transform to equivalent scalar (no need for sparsicity on a vector)
        su_p   = full(tm2sm_v(u_p, ee, 1));

    if length(p)~=length(  pcorners); error('something wrong');end
    %--------------------------------------------------------------------------


    % precribed fluid points are [left_f right_f bottom_f top_f]
    %--------------------------------------------------------------------------
        % MODEL INPUT
        pM     =  01*10^2;
        gradpM =  1*10^6*(e1);
        %delta xf_pp
        dxf_pp = x([left_f right_f bottom_f top_f])- x1p*ones_pp;
        % first-order comp. homogenization
        p_pp   =  pM *ones_pp + dot(gradpM,  dxf_pp );
    if length(pp)~=length([left_f right_f bottom_f top_f]); error('something wrong');end

    % stack both in a hybrid state vector containing disp. and pressure
    wP = [  su_p
             p_pp];

    clear su_p;
    %--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% %% STATIONARY: SOLVE RIGHT/LEFT LINEAR SYSTEM OF EQUATIONS

%--------------------------------------------------------------------------
% right constrain modes
% S_F_P = - inv(mu_F_F  ) * mu_F_P  ;
S_F_P = - mu_F_F\mu_F_P  ;
% S_F_P = - pinv(full(mu_F_F)  ) * mu_F_P  ;
% right stationary hybrid state vector
wr_f  =  S_F_P * wP;

% left constrain modes (Y appears always transposed along the expressions)
% Y_F_P = - inv(mu_F_F.') * mu_P_F.';
Y_F_P = - (mu_F_F.') \ mu_P_F.';
% Y_F_P = - pinv(full(mu_F_F.')) * mu_P_F.';
% left stationary hybrid state vector
wl_f  =  Y_F_P * wP;
%--------------------------------------------------------------------------



% %% STATIONARY: ASSIGN DISP. SOLUTION TO TENSOR VECTOR
%--------------------------------------------------------------------------
% 'ured' is  reduced displacement thanks
% to the elimitation of the dependend points due to periodic bc
%
%  length(ured) = length([dofs_pc; dofs_un; dofs_in])
%  length(u   ) = length([dofs_c ; dofs_un; dofs_in; dofs_de])
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% free (reduced) displacements is the right stationary solution
    % declare
    % nodes are between index 1 to 2*length(f)
    ured_f = wr_f(1:2:2*length(f)-1)*e1+ wr_f(2:2:2*length(f))*e2;
% total reduced displacements (indices ordered as [dofs_un dofs_in])
    % declare ured(ordered: [unconstrained independent] nodes)
    ured = zeros(1, ee, length(p)+length(f), 1);
    % assign prescribed nodes
    ured(p,:) = u_p;
    % assign free nodes
    ured(f,:) = ured_f;
    % free unused variables
    clear ured_f; clear u_p;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% declare displacement tensor vectors (total # of dofs)
vu     = zeros(1, ee,   ns, 1);
vu_aux = zeros(1, ee,   ns, 1); 
% assign. Note vu_aux is in order [dofs_c dofs_un dofs_in dofs_de]
vu_aux = T*ured;
% assign to vu. Note vu is in [nodes_s] order
vu([dofs_c dofs_un dofs_in dofs_de])  = vu_aux;
% free unused variables
clear vu_aux; clear ured;

% prepare for plot: length(vu_plot)>length(vu)
    % declare disp. vector for plot -  whole mesh
    vu_plot = zeros(1, ee,   n, 1);
    % assign vu to vector for plot, only solid nodes are non-zero
    vu_plot(  nodes_s) = vu;
    % assign macro displacement to fluid nodes
%     vu_plot(nodes_f) = uM*ones(nf,1);
%--------------------------------------------------------------------------



% %% STATIONARY: ASSIGN PRESSURE SHAPE MODES TO SCALAR VECTOR vp
%--------------------------------------------------------------------------
% extract pressure values from state vector wr_f
    % free pressures 
    p_fp = wr_f(end-length(fp)+1:end);

% declare pressure scalar vector
vp_steady = zeros(  nf,1);

% assign to scalar vector solution of size
vp_steady(pp,:)=p_pp;
vp_steady(fp,:)=p_fp;

% prepare for plot length(vp_plot)>length(vp)
    % declare pressure scalar vector for plot - whole mesh vector
    vp_plot = zeros(  n,1);
    
    % assign vp to vector for plot
    vp_plot(  nodes_f) = vp_steady;
%--------------------------------------------------------------------------

disp('codeblock: PLOT')

% %% STATIONARY: PLOT PRESSURE/DISP SHAPE (RESCALED DISP.)


% plot the pressure solution and difesplacement
%--------------------------------------------------------------------------
% if ns~=0 && nf~=0
%     figure(5)
%     clf
%     daspect([1 1 1]);
%     % scatter(dot(  x,e1),dot(  x,e2),30,vp_plot,'filled')
%     scatter(dot( x(nodes_f),e1),dot(x(nodes_f),e2),30,vp,'filled')
%     colorbar;
%     pmax = max(abs(vp_plot));
%     caxis([-pmax pmax]);
%     colormap('jet');
%     axis equal
%     hold on
% %     vu_plot_n= vu_plot/norm(norm(vu))*a; rescale='yes';% rescale
%     vu_plot_n= vu_plot; rescale='none';% dont't rescale
%     femplot(  x          ,   conns, 'Color', 'c');
%     femplot(  x+vu_plot_n,   connss{1}, 'Color', 'k');
%     femplot(  x+vu_plot_n,   connss{2}, 'Color', 'r');
%     % title("Stationary problem (check if normalized)" );
% %     title("Steady(check if normalized), x_R=(" ...
% %            +string(dot(xR,e1))+","+string(dot(xR,e2))+")" );
%     title(['Steady (rescaling: ' rescale ')' ] );
%     hold off
% end
%--------------------------------------------------------------------------


% plot the pressure solution and difesplacement
%--------------------------------------------------------------------------
if ns~=0 && nf~=0
    figure(5)
    clf
    daspect([1 1 1]);
    % scatter(dot(  x,e1),dot(  x,e2),30,vp_plot,'filled')
    scatter(dot( x(nodes_f),e1),dot(x(nodes_f),e2),40,vp_steady,'filled')
    % colorbar;
    pmax = max(abs(vp_plot));
    caxis([-pmax pmax]);
    colormap('jet');
    axis equal
    hold on
%     vu_plot_n= vu_plot/norm(norm(vu))*a; rescale='yes';% rescale
    vu_plot_n= vu_plot; rescale='none';% dont't rescale
    % femplot(  x          ,   conns, 'Color', 'c');
%     femplot(  x+vu_plot_n,   connss{1}, 'Color', 'k');
%     femplot(  x+vu_plot_n,   connss{2}, 'Color', 'r');
    mxdeformed =[dot(x+vu_plot_n,e1) dot(x+vu_plot_n,e2)];
    % %colors
    % plotMesh( mxdeformed,  connss{1}, [.7 .7 .7] ); %gray
    % plotMesh( mxdeformed,  connss{2}, [1 0.8 0] ); %yellow
    % plotMesh( mxdeformed,  connss{3}, 'r' );
   
    %white  underformed
    plotMesh( mx,  connss{1}, [1 1 1] ); %gray
    plotMesh( mx,  connss{2}, [1 1 1] ); %yellow
    plotMesh( mx,  connss{3}, [1 1 1] );

    %black
    plotMesh( mxdeformed,  connss{1}, [0 0 0] ); %gray
    plotMesh( mxdeformed,  connss{2}, [0 0 0] ); %yellow
    plotMesh( mxdeformed,  connss{3}, [0 0 0] );

    % %black underformed
    % plotMesh( mx,  connss{1}, [1 1 1] ); %gray
    % plotMesh( mx,  connss{2}, [1 1 1] ); %yellow
    % plotMesh( mx,  connss{3}, [1 1 1] );


    % title("Stationary problem (check if normalized)" );
    title(['Steady (rescaling: ' rescale ')' ] );
    hold off
end

%--------------------------------------------------------------------------


% plot the displacement
%--------------------------------------------------------------------------
if dry
    figure(5)
    clf
    daspect([1 1 1]);
    colorbar;
    axis equal
    hold on
%     vu_plot_n= vu_plot/norm(norm(vu))*a;rescale='yes'; % rescale
    vu_plot_n= vu_plot; rescale='none';% dont't rescale
    femplot(  x          ,   conns    , 'Color', 'c');
    femplot(  x+vu_plot_n,   connss{1}, 'Color', 'k');
    femplot(  x+vu_plot_n,   connss{2}, 'Color', 'r');
    % title("Stationary problem (check if normalized)" );
%     title("Steady(check if normalized), x_R=(" ...
%            +string(dot(xR,e1))+","+string(dot(xR,e2))+")" );
     title(['Steady (rescaling: ' rescale ')' ] );
    hold off
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
clear mxdeformed vu_plot vu_plot_n rescale pmax
%--------------------------------------------------------------------------


%% REDUCED COUPLED DYNAMIC MODEL
% Comments
%--------------------------------------------------------------------------
% 1. check about transposition .' or hermitian transposition '
% regular transposition is being used.
% 2. notation is, for instance, tlam_P_P means matrix with second order
% tensor components of size = (P,P).
% t - matrix with second order tensor components
% v - matrix with first order (vector) tensor components 
% s - matrix with scalar components
%--------------------------------------------------------------------------


% compute reduced matrices on P nodes with Q eigenmodes
%--------------------------------------------------------------------------
%steady part
lam_qs  =  full( ...
                 lam_P_P   +  lam_P_F * S_F_P   +   Y_F_P.' * lam_F_P  ...
                + Y_F_P.'  *  lam_F_F * S_F_P ...
                                                                         );
mu_qs   =  full( ...
                 mu_P_P   +   mu_P_F * S_F_P   +   Y_F_P.' *  mu_F_P  ...
                + Y_F_P.' *   mu_F_F * S_F_P  ...
                                                                         );
lam_Q_P =   psi_F_Q.' * lam_F_P  +   psi_F_Q.' * lam_F_F * S_F_P  ;

lam_P_Q =   lam_P_F   * phi_F_Q  +     Y_F_P.' * lam_F_F * phi_F_Q;


% dynamics part
I_Q_Q    =    eye(n_modes);
LAM_Q_Q  =    psi_F_Q(:,1:n_modes)'*mu_F_F*phi_F_Q(:,1:n_modes);
O_Q_1    =    zeros(n_modes,1);
O_Q_2    =    zeros(n_modes,2);
%--------------------------------------------------------------------------


% convert equivalent scalar matrix back to tensor matrix form
%--------------------------------------------------------------------------
% lam_qs hybrid matrix partitioning
% lam_qs = tlam_p_p   vlam_p_p'
%          vlam_p'_p  slam_p'_p' 
tlam_p_p   = sm2tm(lam_qs(1:dim*length(p),1:dim*length(p)),ee,2);
vlam_p_pp  = sm2tm(lam_qs(1:dim*length(p),dim*length(p)+1:end),ee,1);
% transpose because contraction happens on the column
vlam_pp_p  = sm2tm(lam_qs(dim*length(p)+1:end,1:dim*length(p)).',ee,1).';
slam_pp_pp =       lam_qs(dim*length(p)+1:end,dim*length(p)+1:end);

% mu_qs hybrid matrix partitioning
% mu_qs =  tmu_p_p   vmu_p_p'
%          vmu_p'_p  mu_p'_p'
% note vmu_p'_p is identically zero from the theory
tmu_p_p   = sm2tm(mu_qs(1:dim*length(p),1:dim*length(p)),ee,2);
vmu_p_pp  = sm2tm(mu_qs(1:dim*length(p),dim*length(p)+1:end),ee,1);
% transpose because contraction happens on the column
vmu_pp_p  = sm2tm(mu_qs(dim*length(p)+1:end,1:dim*length(p)).',ee,1).';
smu_pp_pp =       mu_qs(dim*length(p)+1:end,dim*length(p)+1:end);

% lam_P_Q hybrid matrix partitioning (for constitutive relations)
% lam_P_Q = vm_p_Q   
%           sm_p'_Q 
vm_p_Q   = sm2tm(lam_P_Q(1:dim*length(p)    ,:),ee,1);
sm_pp_Q  =       lam_P_Q(dim*length(p)+1:end,:);

% lam_Q_P hybrid matrix partitioning (for evolution equation)
% lam_Q_P = [vm_Q_p  sm_Q_p'] 
% transpose because contraction happens on the column
vm_Q_p   = sm2tm(lam_Q_P(: , 1:dim*length(p)).',ee,1).';
sm_Q_pp  =       lam_Q_P(: , dim*length(p)+1:end  );
%--------------------------------------------------------------------------



% volume fraction
%--------------------------------------------------------------------------
Vf = (ax*ay)*phi    ;
Vs = (ax*ay)*(1-phi);
%--------------------------------------------------------------------------

% compute homogenized macroscopic solid stress   
%--------------------------------------------------------------------------
%
% !! still need to left transpose hEs and hAs !!
% this can be done after
%
% note: h stands for 'homogenized' and s for 'stress'

hsA = 1/(ax*ay)*    l3transpose(dxs_p.' * tlam_p_p  * ones_p  ,ee,3) ;
hsB = 1/(ax*ay)*     ltranspose(dxs_p.' * tlam_p_p  * dxs_p)  ;
hsC = 1/(ax*ay)*                dxs_p.' * vlam_p_pp * ones_pp ;
hsD = 1/(ax*ay)*                dxs_p.' * vlam_p_pp * dxf_pp  ;
hsE = 1/(ax*ay)*    l3transpose(dxs_p.' * tmu_p_p   * ones_p  ,ee,3) ;
hsF = 1/(ax*ay)*     ltranspose(dxs_p.' * tmu_p_p   * dxs_p)  ;
hsG = 1/(ax*ay)*                dxs_p.' * vmu_p_pp * ones_pp  ;
hsH = 1/(ax*ay)*                dxs_p.' * vmu_p_pp * dxf_pp   ;
hsL = 1/(ax*ay)*               (vm_p_Q.'* dxs_p)             ; %corrected 20/04/2023-
% for the evolution equation
hsLs = 1/(ax*ay)*              (vm_Q_p  * dxs_p )             ;%corrected 20/04/2023-
hsl  =                         (vm_Q_p  * dxs_p )             ;


% % convert to matrix
% tC_homog = tm2sm_v(hsF,ee,4);

% full contraction of each tensor for an estimation of order
abshsA = sqrt(  dddot(hsA,hsA) );
abshsB = sqrt( ddddot(hsB,hsB) );
abshsC = sqrt(   ddot(hsC,hsC) );
abshsD = sqrt(  dddot(hsD,hsD) );
abshsE = sqrt(  dddot(hsE,hsE) );
abshsF_BIOT = sqrt( ddddot(hsF,hsF) );
abshsG_BIOT = sqrt(   ddot(hsG,hsG) );
abshsH = sqrt(  dddot(hsH,hsH) );

abshsL=[];
for i=1:size(hsL,1)
abshsL =   [abshsL  sqrt( ddot(hsL(i),hsL(i)) )];
end
abshsL = abs(abshsL); 

abshsLs=[];
for i=1:size(hsLs,1)
abshsLs =   [abshsLs  sqrt( ddot(hsLs(i),hsLs(i)) )];
end
abshsLs = abs(abshsLs); 



% compute homogenized macroscopic solid momentum rate
%--------------------------------------------------------------------------
% still need to left transpose hEm and hAm
% note: h stands for 'homogenized' and m for 'momentum'
hmA = 1/(ax*ay)*                ones_p.' * tlam_p_p  * ones_p  ;
hmB = 1/(ax*ay)*                ones_p.' * tlam_p_p  * dxs_p   ;
hmC = 1/(ax*ay)*                ones_p.' * vlam_p_pp * ones_pp ;
hmD = 1/(ax*ay)*                ones_p.' * vlam_p_pp * dxf_pp  ;
hmE = 1/(ax*ay)*                ones_p.' * tmu_p_p   * ones_p  ;
hmF = 1/(ax*ay)*                ones_p.' * tmu_p_p   * dxs_p   ;
hmG = 1/(ax*ay)*                ones_p.' * vmu_p_pp  * ones_pp ;
hmH = 1/(ax*ay)*                ones_p.' * vmu_p_pp  * dxf_pp  ;
hmL = 1/(ax*ay)*               (ones_p.' * vm_p_Q).'           ;%corrected 20/04/2023-
% for the evolution equation
hmLs= 1/(ax*ay)*               (vm_Q_p  * ones_p)             ; %corrected 20/04/2023-
hml =                          (vm_Q_p  * ones_p)             ;


% full contraction of each tensor for an estimation of order
abshmA_BIOT = sqrt(    ddot(hmA,hmA) );
abshmB = sqrt(   dddot(hmB,hmB) );
abshmC = sqrt(     dot(hmC,hmC) );
abshmD = sqrt(    ddot(hmD,hmD) );
abshmE = sqrt(    ddot(hmE,hmE) );
abshmF = sqrt(abs(dddot(hmF,hmF)) );    %check why negativev
abshmG = sqrt(     dot(hmG,hmG) );
abshmH_BIOT = sqrt(    ddot(hmH,hmH) );

abshmL = [];
for i=1:size(hmL,1)
abshmL =  [abshmL  sqrt( dot(hmL(i),hmL(i)) )];
end

abshmLs = [];
for i=1:size(hmLs,1)
abshmLs =  [abshmLs  sqrt( dot(hmLs(i),hmLs(i)) )];
end


% compute homogenized macroscopic fluid acceleration ( v dot)
%--------------------------------------------------------------------------
% note: h stands for 'homogenized' and m for 'momentum'
hvA = -1/(ax*ay)*                dxf_pp.' * vlam_pp_p  * ones_p  ;
hvB = -1/(ax*ay)*                dxf_pp.' * vlam_pp_p  * dxs_p   ;
hvC = -1/(ax*ay)*                dxf_pp.' * slam_pp_pp * ones_pp ;
hvD = -1/(ax*ay)*                dxf_pp.' * slam_pp_pp * dxf_pp  ;
hvE = -1/(ax*ay)*                dxf_pp.' * vmu_pp_p   * ones_p  ;
hvF = -1/(ax*ay)*                dxf_pp.' * vmu_pp_p   * dxs_p   ;
hvG = -1/(ax*ay)*                dxf_pp.' * smu_pp_pp  * ones_pp ;
hvH = -1/(ax*ay)*                dxf_pp.' * smu_pp_pp  * dxf_pp  ;
hvL = -1/(ax*ay)*               (sm_pp_Q.'* dxf_pp)             ; %corrected 20/04/2023-
% for the evolution equation
hvLs= -1/(ax*ay)*               (sm_Q_pp  * dxf_pp)             ; %corrected 04/05/2023
hvl =                          -(sm_Q_pp  * dxf_pp)             ; %corrected 04/05/2023

% full contraction of each tensor for an estimation of order
abshvA_BIOT = sqrt(    ddot(hvA,hvA) );
abshvB = sqrt(   dddot(hvB,hvB) );
abshvC = sqrt(     dot(hvC,hvC) );
abshvD = sqrt(    ddot(hvD,hvD) );
abshvE = sqrt(    ddot(hvE,hvE) );
abshvF = sqrt(   dddot(hvF,hvF) );
abshvG = sqrt(     dot(hvG,hvG) );
abshvH_BIOT = sqrt(    ddot(hvH,hvH) );

abshvL = [];
for i=1:size(hvL,1)
abshvL =  [abshvL  sqrt( dot(hvL(i),hvL(i)))];
end

abshvLs = [];
for i=1:size(hvLs,1)
abshvLs =  [abshvLs  sqrt( dot(hvLs(i),hvLs(i)))];
end


% compute homogenized macroscopic fluid volumetric rate rate
%--------------------------------------------------------------------------
% note: h stands for 'homogenized' and m for 'momentum'
heA = -1/(ax*ay)*                ones_pp.' * vlam_pp_p  * ones_p  ;
heB = -1/(ax*ay)*                ones_pp.' * vlam_pp_p  * dxs_p   ;
heC = -1/(ax*ay)*                ones_pp.' * slam_pp_pp * ones_pp ;
heD = -1/(ax*ay)*                ones_pp.' * slam_pp_pp * dxf_pp  ;
heE = -1/(ax*ay)*                ones_pp.' * vmu_pp_p   * ones_p  ;
heF = -1/(ax*ay)*                ones_pp.' * vmu_pp_p   * dxs_p   ;
heG = -1/(ax*ay)*                ones_pp.' * smu_pp_pp  * ones_pp ;
heH = -1/(ax*ay)*                ones_pp.' * smu_pp_pp  * dxf_pp  ;
heL = -1/(ax*ay)*               (ones_pp.' * sm_pp_Q).'           ;%corrected 20/04/2023-
% for the evolution equation
heLs = -1/(ax*ay)*              (sm_Q_pp   * ones_pp)             ;%corrected 04/05/2023
hel  =                         -(sm_Q_pp   * ones_pp)             ;%corrected 04/05/2023

% full contraction of each tensor for an estimation of order
absheA = sqrt(     dot(heA,heA) );
absheB_BIOT = sqrt(    ddot(heB,heB) );
absheC_BIOT = sqrt(         heC*heC  );
absheD = sqrt(     dot(heD,heD) );
absheE = sqrt(     dot(heE,heE) );
absheF = sqrt(    ddot(heF,heF) );
absheG = sqrt(         heG*heG  );
absheH = sqrt(     dot(heH,heH) );

absheL = [];
for i=1:size(heL,1)
absheL =  [absheL   sqrt( heL(i)*heL(i) )];
end

absheLs = [];
for i=1:size(heLs,1)
absheLs =  [absheLs   sqrt( heLs(i)*heLs(i) )];
end


     

   

      
%% CONNECT (u,p) to (u,U) coefficients 
% disp("COMPUTE STATIC DENSITIES ACCORDING TO BIOT")

    trho_11   = hmA - dot(hmH,inv(hvH),hvA)
    trho_22   = - phi^2*inv(hvH)
    trho_12   = phi* dot(hmH,inv(hvH))  
    trho_12c  = phi* dot(inv(hvH),hvA) 
    tC_11     = hsF - hsG*heB/heC
    tC_22     = - phi^2/ heC
    tC_12     = phi*hsG/heC  
    tC_12c    = phi*heB/heC             
%     
%     disp("the following should be equal according to Biot")
%     hvA
%     hmH
%     disp("the following should be equal according to Biot")
%     hsG
%     heB
%     % some important homogenized material parameters
    hrhos     = (trho_11  +trho_12)/(1-phi)
    hrhof     = (trho_12c +trho_22)/(phi)
    alpha_inf = dot(trho_22,inv(trho_12+trho_22))


%     if multiple_solid_phase     
%     rhos_weigthed =(V{1}*matProp{1}.rho+V{2}*matProp{2}.rho)/(V{1}+V{2})
%     end    
% 
% 
%    rho_foam  = phi*rhof+(1-phi)*rhos_weigthed_solid    
%    hrho_foam = trho_11  +trho_12+trho_12c +trho_22


%% PLOT COEFFICIENTS' MAGNITUDE VALUE
clear etaticklabels

% build the x tick labels
% mode labels
etaticklabels = cell(1,n_modes);
for i=1:n_modes
    etaticklabels{i} =  ['(' num2str(round(freqs(i),0)) ...
                         ' Hz)$ \eta_{' num2str(i) '} $'];
end     

% make artificial x axis for plot
xx=1:8+n_modes;

figure(6)
clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  

% SOLID MOMENTUM
%--------------------------------------------------------------------------
subplot(2,2,1)
Ym = [abshmA_BIOT abshmB abshmC abshmD abshmE abshmF ...
      abshmG abshmH_BIOT abshmL];
Xm = {'(BIOT) $ \ddot{\vec{u}}_M $', ...
          ' $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
          ' $ \ddot{\mathrm{p}}_{_{M}} $', ...
          ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
          ' $ \vec{u}_{_{M}} $', ...
          ' $\vec{\nabla} \vec{u}_{_{M}} $', ...
          ' $\mathrm{p}_{_{M}} $', ...
          ' (BIOT) $\vec{\nabla} \mathrm{p}_{_{M}} $'};
Xm = [Xm etaticklabels];
plot(xx,log10(Ym), 'kx');
hold on;
plot(xx([1 8]),log10(Ym([1 8])), 'bx');
plot(xx(9:8+n_modes),log10(Ym(9:8+n_modes)), 'gx');
mincoeff = min(Ym(1),Ym(8)); % A and H are Biot coeffs
plot(log10(mincoeff*ones(1,8)), "Color", "b");
hold off;
grid on;
set(gca, 'XTick', 1:length( Xm), 'XTickLabels',  Xm);
title( '$\dot{\vec{\pi}}_{_{M}} $ Solid Momentum coefficients', 'interpreter', 'latex')
ylabel('log10( absolute     value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
clear Xm  mincoeff;


% SOLID STRESS
%--------------------------------------------------------------------------
subplot(2,2,2)
Ys = [abshsA abshsB abshsC abshsD abshsE abshsF_BIOT abshsG_BIOT ...
      abshsH abshsL];
Xs = {' $ \ddot{\vec{u}}_M $', ...
      ' $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
      ' $ \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $ \vec{u}_{_{M}} $', ...
      ' (BIOT) $\vec{\nabla} \vec{u}_{_{M}} $', ...
      ' (BIOT) $\mathrm{p}_{_{M}} $', ...
      ' $\vec{\nabla} \mathrm{p}_{_{M}} $'};
Xs = [Xs etaticklabels];
plot(xx,log10(Ys), 'kx');
hold on;
mincoeff = min(Ys(6),Ys(7)); % F and G are Biot coeffs
plot(xx([6 7]),log10(Ys([6 7])), 'bx');
plot(xx(9:8+n_modes),log10(Ys(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)), "Color", "b");
title( '$\sigma_{_{M}} $ Solid Stress coefficients', 'interpreter', 'latex')
hold off;
grid on;
set(gca, 'XTick', 1:length( Xs), 'XTickLabels',  Xs);
ylabel('log10( absolute value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
clear Xs  mincoeff;


% FLUID VELOCITY
%--------------------------------------------------------------------------
subplot(2,2,3)
Yv = [abshvA_BIOT abshvB abshvC abshvD abshvE abshvF abshvG ...
      abshvH_BIOT abshvL];
Xv = {'(BIOT) $ \ddot{\vec{u}}_M $', ...
      ' $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
      ' $ \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $ \vec{u}_{_{M}} $', ...
      ' $\vec{\nabla} \vec{u}_{_{M}} $', ...
      ' $\mathrm{p}_{_{M}} $', ...
      ' (BIOT) $\vec{\nabla} \mathrm{p}_{_{M}} $'};
Xv = [Xv etaticklabels];
plot(xx,log10(Yv), 'kx');
hold on;
mincoeff = min(Yv(1),Yv(8)); % A and H are Biot coeffs
plot(xx([1 8]),log10(Yv([1 8])), 'bx');
plot(xx(9:8+n_modes),log10(Yv(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)), "Color", "b");
title( '$\dot{\vec{v}}_{_{M}} $ Fluid Velocity coefficients', 'interpreter', 'latex')
hold off;
grid on;
set(gca, 'XTick', 1:length( Xv), 'XTickLabels',  Xv);
ylabel('log10( absolute value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
clear Xv  mincoeff;



% FLUID VOLUMETRIC STRAIN
%--------------------------------------------------------------------------
subplot(2,2,4)
Ye = [absheA absheB_BIOT absheC_BIOT absheD absheE absheF absheG ...
      absheH absheL];
Xe = {' $ \ddot{\vec{u}}_M $', ...
      ' (BIOT) $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
      ' (BIOT) $ \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $ \vec{u}_{_{M}} $', ...
      ' $\vec{\nabla} \vec{u}_{_{M}} $', ...
      ' $\mathrm{p}_{_{M}} $', ...
      ' $\vec{\nabla} \mathrm{p}_{_{M}} $'};
Xe = [Xe etaticklabels];
plot(xx,log10(Ye), 'kx');
hold on;
mincoeff = min(Ye(2),Ye(3)); % B and C are Biot coeffs
plot(xx([2 3]),log10(Ye([2 3])), 'bx');
plot(xx(9:8+n_modes),log10(Ye(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)), "Color", "b");
% legend('coeff. abs', 'smallest biot');
title( '$\ddot{e}_{_{M}} $ Fluid Volumetric Strain coefficients', 'interpreter', 'latex')
hold off;
grid on;
set(gca, 'XTick', 1:length( Xe), 'XTickLabels',  Xe);
ylabel('log10( absolute value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
clear Xe  mincoeff;

clear etaticklabels axi;


% DIMENSIONLESS DRIVING FIELDS (CONSTITUTIVE EQUATIONS)

% characteristic values
%--------------------------------------------------------------------------
% l0 = ay;
% l0 = thi;

% p0 = 0.01;     % conversation
% p0 = 0.2;     % car (max)
% p0 = 0.02;     % car (min)
% p0 = 2;     % jackhammer

% max min
p0 = 1e4;  
% p0 = 1;    

% max min
f0 = 1000;
% f0 = 100;

l0 = 2/f0;    %WAVELENGTH from dispersion curve
u0 = 1e-6*l0;   % small strain
%--------------------------------------------------------------------------

% macroscopic solid stress absolute value
%--------------------------------------------------------------------------
dimless_abshsA      = u0*f0^2     * abshsA      ;
dimless_abshsB      = u0*f0^2/l0  * abshsB      ;
dimless_abshsC      = p0*f0^2     * abshsC      ;
dimless_abshsD      = p0*f0^2/l0  * abshsD      ;
dimless_abshsE      = u0          * abshsE      ;
dimless_abshsF_BIOT = u0/l0       * abshsF_BIOT ;
dimless_abshsG_BIOT = p0          * abshsG_BIOT ;
dimless_abshsH      = p0/l0       * abshsH      ;
dimless_abshsL      = f0^2        * abshsL      ;
%--------------------------------------------------------------------------


% macroscopic solid momentum rate absolute value
%--------------------------------------------------------------------------
dimless_abshmA_BIOT = u0*f0^2     *abshmA_BIOT  ;
dimless_abshmB      = u0*f0^2/l0  *abshmB       ;
dimless_abshmC      = p0*f0^2     *abshmC       ;
dimless_abshmD      = p0*f0^2/l0  *abshmD       ;
dimless_abshmE      = u0          *abshmE       ;
dimless_abshmF      = u0/l0       *abshmF       ;
dimless_abshmG      = p0          *abshmG       ;
dimless_abshmH_BIOT = p0/l0       *abshmH_BIOT  ;
dimless_abshmL      = f0^2        *abshmL       ;
%--------------------------------------------------------------------------


% macroscopic fluid acceleration (v dot) absolute value
%--------------------------------------------------------------------------
dimless_abshvA_BIOT = u0*f0^2     *abshvA_BIOT  ;
dimless_abshvB      = u0*f0^2/l0  *abshvB       ;
dimless_abshvC      = p0*f0^2     *abshvC       ;
dimless_abshvD      = p0*f0^2/l0  *abshvD       ;
dimless_abshvE      = u0          *abshvE       ;
dimless_abshvF      = u0/l0       *abshvF       ;
dimless_abshvG      = p0          *abshvG       ;
dimless_abshvH_BIOT = p0/l0       *abshvH_BIOT  ;
dimless_abshvL      = f0^2        *abshvL       ;
%--------------------------------------------------------------------------

% macroscopic fluid volumetric rate rate absolute value
%--------------------------------------------------------------------------
dimless_absheA      = u0*f0^2     *absheA       ;
dimless_absheB_BIOT = u0*f0^2/l0  *absheB_BIOT  ;
dimless_absheC_BIOT = p0*f0^2     *absheC_BIOT       ;
dimless_absheD      = p0*f0^2/l0  *absheD       ;
dimless_absheE      = u0          *absheE       ;
dimless_absheF      = u0/l0       *absheF       ;
dimless_absheG      = p0          *absheG  ;
dimless_absheH      = p0/l0       *absheH       ;
dimless_absheL      = f0^2        *absheL       ;
%--------------------------------------------------------------------------



% PLOT DIMENSIONLESS COEFFICIENTS' MAGNITUDE VALUE

% mode labels
etaticklabels = cell(1,n_modes);
for i=1:n_modes
    etaticklabels{i} =  ['(' num2str(round(freqs(i),0)) ...
                         ' Hz)$ \eta_{' num2str(i) '}^* $'];
end  

% make artificial x axis for plot
xx=1:8+n_modes;

% without dimensions

figure(7)
clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  
sgtitle([ 'DIMENSIONLESS :'...
         '$ f_0 $ = ' num2str(f0) ...
         '[Hz], $l_0$ =' num2str(l0*100,'%.1f') ...
         '[cm], $ p_0 $ =' num2str(p0) '[Pa]'], 'Interpreter','latex');  
% SOLID MOMENTUM
%--------------------------------------------------------------------------
subplot(2,2,1)
Ym = [dimless_abshmA_BIOT dimless_abshmB dimless_abshmC ... 
      dimless_abshmD dimless_abshmE dimless_abshmF ...
      dimless_abshmG dimless_abshmH_BIOT dimless_abshmL];
Xm = {'(BIOT) $ \ddot{\vec{u}}_M^{*}  $', ...
      ' $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
      ' $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
      ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
      ' $ \vec{u}_{_{M}}^{*}  $', ...
      ' $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*}  $', ...
      ' $\mathrm{p}_{_{M}}^{*}  $', ...
      ' (BIOT) $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
Xm = [Xm etaticklabels];
plot(xx,log10(Ym), 'kx');
hold on;
mincoeff = min(Ym(1),Ym(8)); % A and H are Biot coeffs
% plot(xx([1 8]),log10(Ym([1 8])), 'bx');
scatter(xx([1 8]),log10(Ym([1 8])), 'Marker', 'x', 'Color', [1.0000    0.4118    0.1608]);
plot(xx(9:8+n_modes),log10(Ym(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)), 'Color', [1.0000    0.4118    0.1608]);
hold off;
grid on;
set(gca, 'XTick', 1:length( Xm), 'XTickLabels',  Xm);
title( '$\dot{\vec{\pi}}_{_{M}} $ Solid Momentum coefficients', 'interpreter', 'latex')
ylabel('log10( absolute     value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
% clear Xm Ym mincoeff;


% SOLID STRESS
%--------------------------------------------------------------------------
subplot(2,2,2)
Ys = [dimless_abshsA dimless_abshsB dimless_abshsC ... 
      dimless_abshsD dimless_abshsE dimless_abshsF_BIOT ...
      dimless_abshsG_BIOT dimless_abshsH dimless_abshsL];
Xs = {' $ \ddot{\vec{u}}_M^{*}  $', ...
      ' $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
      ' $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
      ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
      ' $ \vec{u}_{_{M}}^{*}  $', ...
      ' (BIOT) $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*}  $', ...
      ' (BIOT) $\mathrm{p}_{_{M}}^{*}  $', ...
      ' $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
Xs = [Xs etaticklabels];
plot(xx,log10(Ys), 'kx');
hold on;
mincoeff = min(Ys(6),Ys(7)); % F and G are Biot coeffs
% plot(xx([6 7]),log10(Ys([6 7])), 'bx');
scatter(xx([6 7]),log10(Ys([6 7])),  'Marker', 'x', 'Color', [1.0000    0.4118    0.1608]);
plot(xx(9:8+n_modes),log10(Ys(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)),  'Color', [1.0000    0.4118    0.1608]);
title( '$\sigma_{_{M}} $ Solid Stress coefficients', 'interpreter', 'latex')
hold off;
grid on;
set(gca, 'XTick', 1:length( Xs), 'XTickLabels',  Xs);
ylabel('log10( absolute value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
% clear Xs Ys mincoeff;


% FLUID VELOCITY
%--------------------------------------------------------------------------
subplot(2,2,3)
Yv = [dimless_abshvA_BIOT dimless_abshvB dimless_abshvC ...
      dimless_abshvD dimless_abshvE dimless_abshvF dimless_abshvG ...
      dimless_abshvH_BIOT dimless_abshvL];
Xv = {'(BIOT) $ \ddot{\vec{u}}_M^{*}  $', ...
      ' $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
      ' $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
      ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
      ' $ \vec{u}_{_{M}} ^{*} $', ...
      ' $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*} $', ...
      ' $\mathrm{p}_{_{M}}^{*} $', ...
      ' (BIOT) $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
Xv = [Xv etaticklabels];
plot(xx,log10(Yv), 'kx');
hold on;
mincoeff = min(Yv(1),Yv(8)); % A and H are Biot coeffs
% plot(xx([1 8]),log10(Yv([1 8])), 'bx');
scatter(xx([1 8]),log10(Yv([1 8])),  'Marker', 'x', 'Color', [1.0000    0.4118    0.1608]);
plot(xx(9:8+n_modes),log10(Yv(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)), 'Color', [1.0000    0.4118    0.1608]);
% legend('coeff. abs', 'smallest biot');
title( '$\dot{\vec{v}}_{_{M}} $ Fluid Velocity coefficients', 'interpreter', 'latex')
hold off;
grid on;
set(gca, 'XTick', 1:length( Xv), 'XTickLabels',  Xv);
ylabel('log10( absolute value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
% clear Xv Yv mincoeff;



% FLUID VOLUMETRIC STRAIN
%--------------------------------------------------------------------------
subplot(2,2,4)
Ye = [dimless_absheA dimless_absheB_BIOT dimless_absheC_BIOT ...
      dimless_absheD dimless_absheE dimless_absheF ...
      dimless_absheG dimless_absheH dimless_absheL];
Xe = {' $ \ddot{\vec{u}}_M^{*}  $', ...
      ' (BIOT) $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
      ' (BIOT) $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
      ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
      ' $ \vec{u}_{_{M}}^{*}  $', ...
      ' $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*}  $', ...
      ' $\mathrm{p}_{_{M}}^{*}  $', ...
      ' $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
Xe = [Xe etaticklabels];
plot(xx,log10(Ye), 'kx');
hold on;
mincoeff = min(Ye(2),Ye(3)); % B and C are Biot coeffs
% plot(xx([2 3]),log10(Ye([2 3])), 'bx');
scatter(xx([2 3]),log10(Ye([2 3])), 'Marker', 'x', 'Color', [1.0000    0.4118    0.1608]);
plot(xx(9:8+n_modes),log10(Ye(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)), 'Color', [1.0000    0.4118    0.1608]);
% legend('coeff. abs', 'smallest biot');
title( '$\ddot{e}_{_{M}} $ Fluid Volumetric Strain coefficients', 'interpreter', 'latex')
hold off;
grid on;
set(gca, 'XTick', 1:length( Xe), 'XTickLabels',  Xe);
ylabel('log10( absolute value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------



% %--------------------------------------------------------------------------
% clear Xe Ye mincoeff
% clear etaticklabels axi
% clear abshmA_BIOT abshmB abshmC abshmD abshmE abshmF ...
%       abshmG abshmH_BIOT abshmL ...
%       abshsA abshsB abshsC abshsD abshsE abshsF_BIOT abshsG_BIOT ...
%       abshsH abshsL ...
%       abshvA_BIOT abshvB abshvC abshvD abshvE abshvF abshvG ...
%       abshvH_BIOT abshvL ...
%       absheA absheB_BIOT absheC_BIOT absheD absheE absheF absheG ...
%       absheH absheL 
% clear dimless_abshmA_BIOT dimless_abshmB dimless_abshmC ... 
%       dimless_abshmD dimless_abshmE dimless_abshmF ...
%       dimless_abshmG dimless_abshmH_BIOT dimless_abshmL ...
%       dimless_abshsA dimless_abshsB dimless_abshsC ... 
%       dimless_abshsD dimless_abshsE dimless_abshsF_BIOT ...
%       dimless_abshsG_BIOT dimless_abshsH dimless_abshsL ...
%       dimless_abshvA_BIOT dimless_abshvB dimless_abshvC ...
%       dimless_abshvD dimless_abshvE dimless_abshvF dimless_abshvG ...
%       dimless_abshvH_BIOT dimless_abshvL ...
%       dimless_absheA dimless_absheB_BIOT dimless_absheC_BIOT ...
%       dimless_absheD dimless_absheE dimless_absheF ...
%       dimless_absheG dimless_absheH dimless_absheL ...
% %--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% PLOT SOLID STIFFNESS TENSOR PER TOTAL(OR SOLID) VOLUME


% stiffness
C4  = tm2sm_v(hsF,ee,4);
C4v = t2voigt(C4 ,ee,4)
% C1111_over_C2121= C4v(1,1)/C4v(3,3)

% % if the following bulk modulus are different then it is not isotropic RVE
BulkModulus1111 = C4v(1,1)-4/3*C4v(3,3)
BulkModulus1122 = C4v(1,2)+2/3*C4v(3,3)
% 
% % Efoam = (C4v(1,1)-C4v(1,2))*(C4v(1,1)+2*C4v(1,2))/(C4v(1,1)+C4v(1,2))
% % 
% % 
% % C11  = tm2sm_v(tC_11,ee,4);
% % C11v = t2voigt(C11 ,ee,4)
% % BulkModulus1111 = C11v(1,1)-4/3*C11v(3,3)
% % BulkModulus1122 = C11v(1,2)+2/3*C11v(3,3)
% 
% if dry
% Kb= (BulkModulus1111+BulkModulus1122 )/2;
% N = C4v(4,4);
% [P,Q,R]=get_biot_coeffs_gedanken(phi,Kb,rhof*c^2,N,matProp{1}.kappa);
% C4_BIOT = t2voigt( tm2sm_v( (P-2*N)*I*I+ 2*N*I4S,ee,4 ) ,ee,4)
% else
% 
% % hC_11 = t2voigt(tm2sm_v(tC_11,ee,4),ee,4)
% end
% 
% if multiple_solid_phase     
% rhos_weigthed_solid =(V{1}*matProp{1}.rho+V{2}*matProp{2}.rho)/(V{1}+V{2})
% end




%% FREE SOME MEMORY
clear A_pp_pp A A_fp_fp A_fp_pp A_pp_fp ...
      B B_fp_fp B_fp_pp B_pp_fp B_pp_pp ...
      C C_ C_f_fp C_f_pp C_p_fp C_p_pp Cd Ce Ci Cu ...
      K K_ K_f_f K_f_p K_p_f K_p_p Kdd Kdi Kdu Kid ...
      Kii Kiu Kud Kui Kuu ...
      M M_ M_f_f M_f_p M_p_f M_p_p Mdd Mdi Mdu Mid ...
      Miu Mii Mud Mui Muu ...
      O_pp_p O_pp_f O_p_pp O_p_fp O_fp_p O_fp_f O_f_pp O_f_fp






%% DISPERSION BLOCH ANALYSIS (polarization)



% % pol
% figure(23)
% hold on
% grid on
% scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)),30,real(dispersion_bloch(:,3)),Marker="+")
% % scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)),30,'ko')
% 
% colormap("cool");
% ylim([0 1600])
% box on;
% 
% % one color (black)
%  % scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)),"k+"');


% black crosses
% figure(23)
hold on
grid on
scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)),30,real(dispersion_bloch(:,3)),"MarkerEdgeColor","black",Marker="+")
colormap("cool");
 ylim([0 8000])


% % pol   
% % figure(23)
% hold on
% grid on
% scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)), "MarkerEdgeColor","black", Marker="+")
% % colormap("cool");
% box on;
% ylim([0 1600])


%% DISPERSION CURVE BIOT STEADY TERMS 

% test eigenvector othogonality
isorthogonal=[];

% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_branches = 3; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
% range_of_dimless_k = -1:kstep:2;
% range_of_dimless_k = 0:kstep:0.25;
% range_of_dimless_k = 0:kstep:0.02;
% range_of_dimless_k = 0:kstep:1;
range_of_dimless_k = -1:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda    = zeros(number_of_branches,length(range_of_dimless_k)); 
wavepol_t = zeros(number_of_branches,length(range_of_dimless_k));
wavepol_s = zeros(number_of_branches,length(range_of_dimless_k)); 
wavepol_f = zeros(number_of_branches,length(range_of_dimless_k)); 
waveborn  = zeros(number_of_branches,length(range_of_dimless_k)); 
%--------------------------------------------------------------------------


i=1;
% loop in k
%--------------------------------------------------------------------------
    for dimless_k = range_of_dimless_k
    
    
    % wavenumber
    %--------------------------------------------------------------------------
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end
    %--------------------------------------------------------------------------

    
    % build matrix
    %--------------------------------------------------------------------------
    % Density matrix (k)
    % first line of M
    Muu= full(tm2sm_v(hmA,ee,2));
    Mup= [0;0];
    
    % second line of M
    Mpu= full(1i*tm2sm_v(dot(-heB,k)+dot(k,hvA),ee,1)).';
    Mpp= heC ;  
    
  
    % Elasticity matrix (k)
    % first line of K
    Kuu=full(tm2sm_v(dot(k,hsF,k),ee,2));
     % Kuu= tm2sm_v(hmE+dot(k,hsF,k),ee,2);
    Kup=full(1i*tm2sm_v(dot(-hmH,k)+dot(k,hsG),ee,1));
    
    % second line of K
    Kpu=    [0 0];
    Kpp=    dot(k,hvH,k) ;
    
    % density matrix 
    M=[Muu Mup 
       Mpu Mpp];
    % elasticity matrix 
    K=[Kuu Kup 
       Kpu Kpp];
    %--------------------------------------------------------------------------



    % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------

    % split displacement and pressure/omega2
    %--------------------------------------------------------------------------
    us               = eigenvectors([1 2],:);
    p_over_omega2    = eigenvectors(3,:) ./omega2.';
    %--------------------------------------------------------------------------            

    % compute fluid-disp/aggregate using constitutive relatio
    %--------------------------------------------------------------------------
    uf =1/phi*[tm2sm_v(hvA,ee,2) 1i*tm2sm_v(dot(hvH,k),ee,1)]*[us        ;
                                                             p_over_omega2];
    %--------------------------------------------------------------------------                                            

   
   % eigen frequencies
   %--------------------------------------------------------------------------
   lambda(:,i) = omega2;
   %--------------------------------------------------------------------------


   % total displacement 
   ut = (1-phi)*us + phi*uf;
   

   % wave born
   %--------------------------------------------------------------------------
   waveborn(:,i)   = sqrt(sum(uf.*conj(uf),1))./ sqrt(sum(us.*conj(us),1));
   %--------------------------------------------------------------------------


   % total displacement polarization
   %--------------------------------------------------------------------------
   % in k direction, not k is real so conj is useless
   ut_dot_ek = sum( ut.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);
  
   % wave polarization 
   wavepol_t(:,i)   = sqrt(sum(ut_dot_ek.*conj(ut_dot_ek),1))./ sqrt(sum(ut.*conj(ut),1));
   %--------------------------------------------------------------------------
   
   % solid displacement polarization
   %--------------------------------------------------------------------------
   % in k direction, not k is real so conj is useless
   us_dot_ek = sum( us.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);
  
   % wave polarization
   wavepol_s(:,i)   = sqrt(sum(us_dot_ek.*conj(us_dot_ek),1))./ sqrt(sum(us.*conj(us),1));
   %--------------------------------------------------------------------------

   % fluid displacement polarization
   %--------------------------------------------------------------------------
   % in k direction, not k is real so conj is useless
   uf_dot_ek = sum( uf.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);
  
   % wave polarization
   wavepol_f(:,i)   = sqrt(sum(uf_dot_ek.*conj(uf_dot_ek),1))./ sqrt(sum(uf.*conj(uf),1));
   %-------------------------------------------------------------------------- 
    

   % % are the eigenvectors orthonormal?
   % %--------------------------------------------------------------------------
   % isorthogonal = [isorthogonal ...
   %                   abs(eigenvectors(:,1)'*eigenvectors(:,2) + ...
   %                       eigenvectors(:,1)'*eigenvectors(:,3) + ...
   %                       eigenvectors(:,2)'*eigenvectors(:,3))  ];
   % %-------------------------------------------------------------------------- 
   
   % % are the eigenvectors orthogonal?
   %  %--------------------------------------------------------------------------
   %  isorthogonal = [isorthogonal ...
   %                   abs(eigenvectors(:,1)'*eigenvectors(:,2)) + ...
   %                   abs(eigenvectors(:,1)'*eigenvectors(:,3)) + ...
   %                   abs(eigenvectors(:,2)'*eigenvectors(:,3))  ];
   %  %--------------------------------------------------------------------------

   % are the eigenvectors orthonormal?
   %--------------------------------------------------------------------------
   isorthogonal = [isorthogonal ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,2))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,3))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,3)))  ];
   %-------------------------------------------------------------------------- 

   % % are the eigenvectors orthonormal?
   % %--------------------------------------------------------------------------
   % isorthogonal = [isorthogonal ...
   %                   abs(imag(eigenvectors(:,1))'*imag(eigenvectors(:,2)) + ...
   %                       imag(eigenvectors(:,1))'*imag(eigenvectors(:,3)) + ...
   %                       imag(eigenvectors(:,2))'*imag(eigenvectors(:,3)) )  ];
   % %-------------------------------------------------------------------------- 



    i=i+1;
    end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
%--------------------------------------------------------------------------



% plot log10 imag part of lambda
%--------------------------------------------------------------------------
% figure
% hold on
% grid on
% log10_im_lambda = log10(abs(imag(lambda)));
% for i=1:number_of_branches
% % scatter(k_plot,real(frequencies(i,:)),'b.')
% scatter(range_of_dimless_k,log10_im_lambda(i,:), '.')
% end
%--------------------------------------------------------------------------



% plot orthogonality ratio
%--------------------------------------------------------------------------
% % c_gx =  2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure
% hold on;
% grid on;
% for i=1:number_of_branches
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(isorthogonal),'filled' );
% colorbar;
% end
% % ylim([0 9000])
% % title("slope \Gamma-X is c = " + num2str(c_gx)+"[m/s]") ; 
% grid on;
%--------------------------------------------------------------------------



% % plot dispersion (one color)
% %--------------------------------------------------------------------------

% c_gx =  2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X

hold on;
for i=1:number_of_branches
% scatter(range_of_dimless_k,real(frequencies(i,:)),'b.') ;
plot(range_of_dimless_k,real(frequencies(i,:)),'b--') ;
% scatter(k_plot,real(frequencies(i,:)),5, wavepol_t(i,:),'filled') ;
% colormap("cool");
% colorbar('Ticks',[0.01,0.99],...
%          'TickLabels',{'u \perp k', 'u || k'})
end
% ylim([0 9000])
% title("slope \Gamma-X is c = " + num2str(c_gx)+"[m/s]") ; 
grid on;
% %--------------------------------------------------------------------------





% % plot dispersion and wave polarization
% %--------------------------------------------------------------------------
% frequencies = sqrt(lambda)/2/pi;
% c_gx =  2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% 
% hold on;
% for i=1:number_of_branches
% % plot(range_of_dimless_k,real(frequencies(i,:)),'b--') ;
% scatter(range_of_dimless_k,real(frequencies(i,:)),2, wavepol_s(i,:),'filled') ;
% 
% 
% colormap("cool");
% colorbar('Ticks',[0.01,0.99],...
%          'TickLabels',{'u \perp k', 'u || k'})
% end
% % ylim([0 9000])
% title("slope \Gamma-X is c = " + num2str(c_gx)+"[m/s]") ; 
% grid on;
% %--------------------------------------------------------------------------


% plot preferential propagation per branch
% %--------------------------------------------------------------------------
% % figure(4)
% hold on;
% frequencies = sqrt(lambda)/2/pi;
% for i=1:number_of_branches
% % scatter(range_of_dimless_k,real(frequencies(i,:)), '.')
% % scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10( phi/(1-phi)*waveborn(i,:)),'filled' );
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(             waveborn(i,:)),'filled' );
% end
% grid on;
% colorbar('Ticks',[0.001,2.99],...
%          'TickLabels',{'solid', 'fluid'})
% 
% hold off
% %--------------------------------------------------------------------------

%% DISPERSION CURVE BIOT TERMS + m-MODES

% test eigenvector othogonality
isorthogonal=[];

% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_branches = 3+n_modes; % propagating 
%------------------------------------------------   --------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
% range_of_dimless_k = -1:kstep:2;
% range_of_dimless_k = 0:kstep:0.2;
range_of_dimless_k = 0:kstep:1;
% range_of_dimless_k = -1:kstep:1;
% range_of_dimless_k = 0:kstep:0.05;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda    = zeros(number_of_branches,length(range_of_dimless_k)); 
wavepol_t = zeros(number_of_branches,length(range_of_dimless_k));
wavepol_s = zeros(number_of_branches,length(range_of_dimless_k)); 
wavepol_f = zeros(number_of_branches,length(range_of_dimless_k)); 
waveborn  = zeros(number_of_branches,length(range_of_dimless_k)); 
%--------------------------------------------------------------------------

i=1;
% loop in k
%--------------------------------------------------------------------------
    for dimless_k = range_of_dimless_k
    

    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k   *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    % k =                 0*pi/ax *e1 -  dimless_k    *pi/ay *e2;warning("Gamma-X and Gamma-Y");
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2; 
    end
    
    % Density matrix (k)
    % first line of M
    Muu=(tm2sm_v(hmA,ee,2));
    Mup=(0* tm2sm_v(hmC,ee,1));
    Mue=(tm2sm_v(hmL.',ee,1))   ;
    
    % second line of M
    Mpu= (1i*tm2sm_v(dot(-heB,k)+dot(k,hvA),ee,1)).';
    Mpp=        heC ;
    Mpe=        0*heL.'                 ;
    
    % third line of M
    Meu=(tm2sm_v(+ hml.',ee,1)).'   ;
    Mep=      - 0*hel            ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix (k)
    % first line of K
    Kuu=(tm2sm_v(dot(k,hsF,k),ee,2));
    Kup=(1i*tm2sm_v(dot(-hmH,k)+dot(k,hsG),ee,1));
    Kue= O_Q_2.' ;
    
    % second line of K
    Kpu=full((0*tm2sm_v(heE+dot(k,hvF,k),ee,1))).';
    Kpp=       dot(k,hvH,k)   ;
    Kpe= O_Q_1.' ;
    
    % third line of K
    Keu= O_Q_2 ; Kep= O_Q_1 ; Kee= LAM_Q_Q ;                                                                 

    % density matrix 
    M=[Muu Mup Mue;
       Mpu Mpp Mpe;
       Meu Mep Mee];
    % elasticity matrix 
    K=[Kuu Kup Kue;
       Kpu Kpp Kpe;
       Keu Kep Kee];
    
  % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------
    


    % are the eigenvectors orthonormal?
    %--------------------------------------------------------------------------
    isorthogonal = [isorthogonal ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,2))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,3))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,4))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,5))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,3))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,4))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,5))) + ...
                     abs(real(eigenvectors(:,3))'*real(eigenvectors(:,4))) + ...
                     abs(real(eigenvectors(:,3))'*real(eigenvectors(:,5))) + ...
                     abs(real(eigenvectors(:,4))'*real(eigenvectors(:,5))) ];
    %--------------------------------------------------------------------------




    % split displacement and pressure/omega2
    %--------------------------------------------------------------------------
    us               = eigenvectors([1 2],:);
    p_over_omega2    = eigenvectors(3,:) ./omega2.';
    eta              = eigenvectors(4:end,:);
    %--------------------------------------------------------------------------            

    % compute fluid-disp/aggregate using constitutive relatio
    %--------------------------------------------------------------------------
    uf =1/phi* ...
        [tm2sm_v(hvA,ee,2) 1i*tm2sm_v(dot(hvH,k),ee,1) tm2sm_v(hvL.',ee,1)]*...
                                                            [us           ;
                                                             p_over_omega2;
                                                             eta];
    %--------------------------------------------------------------------------                                            

   
   % eigen frequencies
   %--------------------------------------------------------------------------
   lambda(:,i) = omega2;
   %--------------------------------------------------------------------------


   % total displacement 
   ut = (1-phi)*us + phi*uf;


   % wave born
   %--------------------------------------------------------------------------
   waveborn(:,i)   = sqrt(sum(uf.*conj(uf),1))./ sqrt(sum(us.*conj(us),1));
   %--------------------------------------------------------------------------


   % total displacement polarization
   %--------------------------------------------------------------------------
   % in k direction, not k is real so conj is useless
   ut_dot_ek = sum( ut.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);

   % wave polarization 
   wavepol_t(:,i)   = sqrt(sum(ut_dot_ek.*conj(ut_dot_ek),1))./ sqrt(sum(ut.*conj(ut),1));
   %--------------------------------------------------------------------------

   % solid displacement polarization
   %--------------------------------------------------------------------------
   % in k direction, not k is real so conj is useless
   us_dot_ek = sum( us.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);

   % wave polarization
   wavepol_s(:,i)   = sqrt(sum(us_dot_ek.*conj(us_dot_ek),1))./ sqrt(sum(us.*conj(us),1));
   %--------------------------------------------------------------------------

   % fluid displacement polarization
   %--------------------------------------------------------------------------
   % in k direction, not k is real so conj is useless
   uf_dot_ek = sum( uf.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);

   % wave polarization
   wavepol_f(:,i)   = sqrt(sum(uf_dot_ek.*conj(uf_dot_ek),1))./ sqrt(sum(uf.*conj(uf),1));
   %-------------------------------------------------------------------------- 
    
    i=i+1;
    end
%--------------------------------------------------------------------------



% % plot log10 imag part of lambda
% %--------------------------------------------------------------------------
% figure
% hold on
% grid on
% log10_im_lambda = log10(abs(imag(lambda)));
% for i=1:number_of_branches
% % scatter(k_plot,real(frequencies(i,:)),'b.')
% plot(range_of_dimless_k,log10_im_lambda(i,:))
% end
% %--------------------------------------------------------------------------
% 
% 
% 
% % plot orthogonality ratio
% %--------------------------------------------------------------------------
% % c_gx =  2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-
% frequencies = sqrt(lambda)/2/pi;
% figure
% hold on;
% grid on;
% for i=1:number_of_branches
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(isorthogonal),'filled' );
% colorbar;
% end
% % ylim([0 9000])
% % title("slope \Gamma-X is c = " + num2str(c_gx)+"[m/s]") ; 
% grid on;
% %--------------------------------------------------------------------------



% % plot dispersion and wave polarization
% %--------------------------------------------------------------------------
% frequencies = sqrt(lambda)/2/pi;
% c_gx =  2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% 
% hold on;
% for i=1:number_of_branches
% scatter(range_of_dimless_k,real(frequencies(i,:)),3, wavepol_s(i,:),'filled') ;
% % scatter(k_plot,real(frequencies(i,:)),5, wavepol_t(i,:),'filled') ;
% colormap("cool");
% colorbar('Ticks',[0.01,0.99],...
%          'TickLabels',{'u_t \perp k', 'u_t || k'})
% end
% % ylim([0 9000])
% % title("slope \Gamma-X is c = " + num2str(c_gx)+"[m/s]") ; 
% grid on;
% %--------------------------------------------------------------------------

% % plot
% %--------------------------------------------------------------------------
% frequencies = sqrt(lambda)/2/pi;
% % c_x = 2*pi*(frequencies(:,2)-frequencies(:,1)) / (kstep*pi/ax); % speed of sound gamma-X
% % figure(24)
% % hold on
% % grid on
% for i=1:number_of_branches
% scatter(range_of_dimless_k,real(frequencies(i,:)),'r.')
% end
% % title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% ylim([0 9000])
% % hold off
% %--------------------------------------------------------------------------

% plot preferential propagation per branch
%--------------------------------------------------------------------------
% % figure(4)
% hold on;
% frequencies = sqrt(lambda)/2/pi;
% for i=1:number_of_branches
% % scatter(range_of_dimless_k,real(frequencies(i,:)), '.')
% % scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10( phi/(1-phi)*waveborn(i,:)),'filled' );
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(waveborn(i,:)),'filled' );
% end
% grid on;
% colorbar;
% title("wave amplitude log_{10} \phi u_f /((1-\phi)*u_s)") ; 
% hold off
%--------------------------------------------------------------------------


% % plot one color
% 
% %--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
% c_x = 2*pi*(frequencies(:,2)-frequencies(:,1)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(24)
hold on
grid on
for i=1:number_of_branches
plot(range_of_dimless_k,real(frequencies(i,:)),'r', "linewidth", 2)
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% xlim([0 1])
% hold off
% %-------------



%% DISPERSION CURVE ALL TERMS

% eigenvector orthogonality indicator
isorthogonal =[];

% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_branches = 3+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%-------------------------------------------------------figurefigure-------------------
kstep=0.01;
% range_of_dimless_k = -1:kstep:2;
% range_of_dimless_k = 0:kstep:1;
% range_of_dimless_k = 0:kstep:.05;
range_of_dimless_k = -1:kstep:1;warning("Gamma-X and Gamma-Y");
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_branches,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------

i=1;
% loop in k
%--------------------------------------------------------------------------
    for dimless_k = range_of_dimless_k


     % wavenumber
    if(dimless_k<0)
    % k =  - dimless_k   *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    k =                 0*pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end
    
    % Density matrix (k)
    % first line of M
    Muu=tm2sm_v(hmA+dot(k,hsB,k),ee,2)+1i*tm2sm_v(dot(-hmB,k)+dot(k,hsA),ee,2);
    Mup=tm2sm_v(hmC+dot(k,hsD,k),ee,1)+1i*tm2sm_v(dot(-hmD,k)+dot(k,hsC),ee,1);
    Mue=tm2sm_v(hmL.',ee,1)           +1i*tm2sm_v(dot(k,hsL.'),ee,1);
    
    % second line of M
    Mpu=(tm2sm_v(heA+dot(k,hvB,k),ee,1)+1i*tm2sm_v(dot(-heB,k)+dot(k,hvA),ee,1)).';
    Mpp=        heC+dot(k,hvD,k)      +1i*       (dot(-heD,k)+dot(k,hvC)     );
    Mpe=        heL.'                 +1i*        dot(k,hvL.')                ;
    
    % third line of M
    Meu=tm2sm_v(+ hml.',ee,1).'      +1i*tm2sm_v(dot( - hsl.',k),ee,1).'  ;
    Mep=        - hel                +1i*        dot( + hvl  ,k)          ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix (k)
    % first line of K
    Kuu=tm2sm_v(hmE+dot(k,hsF,k),ee,2)+1i*tm2sm_v(dot(-hmF,k)+dot(k,hsE),ee,2);
    Kup=tm2sm_v(hmG+dot(k,hsH,k),ee,1)+1i*tm2sm_v(dot(-hmH,k)+dot(k,hsG),ee,1);
    Kue= O_Q_2.' ;
    
    % second line of K
    Kpu=(tm2sm_v(heE+dot(k,hvF,k),ee,1)+1i*tm2sm_v(dot(-heF,k)+dot(k,hvE),ee,1)).';
    Kpp=        heG+dot(k,hvH,k)      +1i*       (dot(-heH,k)+dot(k,hvG)     );
    Kpe= O_Q_1.' ;
    
    % third line of K
    Keu= O_Q_2 ; Kep= O_Q_1 ; Kee= LAM_Q_Q ;                                                                    

    % density matrix 
    M=[ Muu    Mup   Mue;
        Mpu    Mpp   Mpe;
        Meu    Mep   Mee]; 
    % elasticity matrix 
    K=[ Kuu   Kup   Kue;
        Kpu   Kpp   Kpe;
        Keu   Kep   Kee];
    
    %     % compute eigevalues
    % %     [V,D]=eig(K,M)
    %     lambda(:,i)=sort(eig(K,M));       
    
    % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------
    
    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = omega2;
    %--------------------------------------------------------------------------


    % % are the eigenvectors orthonormal?
    % %--------------------------------------------------------------------------
    % isorthogonal = [isorthogonal ...
    %                  abs(eigenvectors(:,1)'*eigenvectors(:,2) + ...
    %                      eigenvectors(:,1)'*eigenvectors(:,3) + ...
    %                      eigenvectors(:,2)'*eigenvectors(:,3)) ];
    % %--------------------------------------------------------------------------    
    
    
    % % are the eigenvectors orthonormal?
    % %--------------------------------------------------------------------------
    % isorthogonal = [isorthogonal ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,2))) + ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,3))) + ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,4))) + ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,5))) + ...
    %                  abs(real(eigenvectors(:,2))'*real(eigenvectors(:,3))) + ...
    %                  abs(real(eigenvectors(:,2))'*real(eigenvectors(:,4))) + ...
    %                  abs(real(eigenvectors(:,2))'*real(eigenvectors(:,5))) + ...
    %                  abs(real(eigenvectors(:,3))'*real(eigenvectors(:,4))) + ...
    %                  abs(real(eigenvectors(:,3))'*real(eigenvectors(:,5))) + ...
    %                  abs(real(eigenvectors(:,4))'*real(eigenvectors(:,5))) ];
    % %--------------------------------------------------------------------------
    % are the eigenvectors orthonormal?
    %--------------------------------------------------------------------------
    isorthogonal = [isorthogonal ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,2))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,3))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,4))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,3))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,4))) + ...
                     abs(real(eigenvectors(:,3))'*real(eigenvectors(:,4))) ];

                     
                     % abs(real(eigenvectors(:,1))'*real(eigenvectors(:,5))) + ...
                     % abs(real(eigenvectors(:,2))'*real(eigenvectors(:,4))) + ...
                     % abs(real(eigenvectors(:,2))'*real(eigenvectors(:,5))) + ...
                     
                     % abs(real(eigenvectors(:,3))'*real(eigenvectors(:,5))) + ...
                     % abs(real(eigenvectors(:,4))'*real(eigenvectors(:,5))) ];
    %--------------------------------------------------------------------------

    % split displacement and pressure/omega2
    %--------------------------------------------------------------------------
    us               = eigenvectors([1 2],:);
    p_over_omega2    = eigenvectors(3,:) ./omega2.';
    eta              = eigenvectors(4:end,:);
    %--------------------------------------------------------------------------            

    % % compute fluid-disp/aggregate using constitutive relatio
    % %--------------------------------------------------------------------------
    % uf =1/phi* ...
    %     [tm2sm_v(hvA,ee,2) 1i*tm2sm_v(dot(hvH,k),ee,1) tm2sm_v(hvL.',ee,1)]*...
    %                                                         [us           ;
    %                                                          p_over_omega2;
    %                                                          eta];
    % %--------------------------------------------------------------------------                                            

  

   % % total displacement 
   % ut = (1-phi)*us + phi*uf;


   % % wave born
   % %--------------------------------------------------------------------------
   % waveborn(:,i)   = sqrt(sum(uf.*conj(uf),1))./ sqrt(sum(us.*conj(us),1));
   % %--------------------------------------------------------------------------


   % % total displacement polarization
   % %--------------------------------------------------------------------------
   % % in k direction, not k is real so conj is useless
   % ut_dot_ek = sum( ut.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);
   % 
   % % wave polarization 
   % wavepol_t(:,i)   = sqrt(sum(ut_dot_ek.*conj(ut_dot_ek),1))./ sqrt(sum(ut.*conj(ut),1));
   % %--------------------------------------------------------------------------
   

   % % solid displacement polarization
   % %--------------------------------------------------------------------------
   % % in k direction, not k is real so conj is useless
   % us_dot_ek = sum( us.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);
   % 
   % % wave polarization
   % wavepol_s(:,i)   = sqrt(sum(us_dot_ek.*conj(us_dot_ek),1))./ sqrt(sum(us.*conj(us),1));
   % %--------------------------------------------------------------------------

   % % fluid displacement polarization
   % %--------------------------------------------------------------------------
   % % in k direction, not k is real so conj is useless
   % uf_dot_ek = sum( uf.* conj(([dot(k,e1);dot(k,e2)]/norm(k).*ones(2,number_of_branches))) ,1);

   % % wave polarization
   % wavepol_f(:,i)   = sqrt(sum(uf_dot_ek.*conj(uf_dot_ek),1))./ sqrt(sum(uf.*conj(uf),1));
   % %-------------------------------------------------------------------------- 
    
    i=i+1;
    end
%--------------------------------------------------------------------------



%
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
%--------------------------------------------------------------------------


% % plot log10 imag part of lambda
% %--------------------------------------------------------------------------
% figure
% hold on
% grid on
% log10_im_lambda = log10(abs(imag(lambda)));
% for i=1:number_of_branches
% % scatter(k_plot,real(frequencies(i,:)),'b.')
% scatter(range_of_dimless_k,log10_im_lambda(i,:), '.')
% end
% %--------------------------------------------------------------------------
% 
% 
% 
% % plot orthogonality ratio
% %--------------------------------------------------------------------------
% % c_gx =  2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure
% hold on;
% grid on;
% for i=1:number_of_branches
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(isorthogonal),'filled' );
% colorbar;
% end
% % ylim([0 9000])
% % title("slope \Gamma-X is c = " + num2str(c_gx)+"[m/s]") ; 
% grid on;
% %--------------------------------------------------------------------------

%plot one color
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,2)-frequencies(:,1)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
hold on
grid on
for i=1:number_of_branches
% scatter(range_of_dimless_k,real(frequencies(i,:)),'r.')
plot(range_of_dimless_k,real(frequencies(i,:)),'r', "linewidth", 2)
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(isorthogonal),'filled' );
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 1600])
hold off
%--------------------------------------------------------------------------

% % plot preferential propagation per branch (CONSTITUTIVE VALUE OF uF is
% different
% %--------------------------------------------------------------------------
% % figure(4)
% hold on;
% frequencies = sqrt(lambda)/2/pi;
% for i=1:number_of_branches
% % scatter(range_of_dimless_k,real(frequencies(i,:)), '.')
% % scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10( phi/(1-phi)*waveborn(i,:)),'filled' );
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(  waveborn(i,:)),'filled' );
% end
% grid on;
% colorbar;
% title("wave amplitude log_{10} \phi u_f /((1-\phi)*u_s)") ; 
% hold off
% %--------------------------------------------------------------------------



