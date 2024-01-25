function [M,C,K,fnl,fext, outdof, MyAssembly] = buildModel(nDiscretization)
% we tune w to trigger 1:2 resonance

%% Finite Element Setup
% Geometry
startLIN = tic;
l = 1;  % width of domain [m]
b = 2;  % length of domain [m]
t = 1e-2; % thickness of plate [m]
w = 0.041; % curvature parameter (height of the midpoint relative to ends) [m]
% material properties
E       = 70e9;  % 70e9 % 200e9 % Young's modulus [Pa]
rho     = 2700; % 2700 % 7850 % density [kg/m^3]
nu      = 0.33;    % Poisson's ratio 
kappa   = 1e5; % material damping modulus 1e8

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor =  @()TriShellElement(t, myMaterial); % same element all across the domain

% Meshing the geometry
nl = nDiscretization;
nb = 2*nDiscretization; 
[nodes,elements,bnodes] = RectangularMesh(l,b,nl,nb,w);     % Rectangular Mesh definition

% creating Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

% Plot mesh
figure('Name','Mesh'); PlotMesh(nodes,elements,0);

%% Assemble linear stiffness, mass and damping
disp('Assembling M,C,K matrices')
% % parallelized assembly
% cluster = parcluster('local');
% cluster.NumWorkers = 4;
% parpool(cluster, 4)
% MyAssembly = Assembly(myMesh,true); 

MyAssembly = Assembly(MyMesh);
K = MyAssembly.stiffness_matrix();
M = MyAssembly.mass_matrix();
% C = MyAssembly.damping_matrix();
MyAssembly.DATA.K = K;
MyAssembly.DATA.M = M;
%% apply boundary conditions
disp('Applying boundary conditions')
MyMesh.set_essential_boundary_condition([bnodes{3}, bnodes{4}],1:3,0) % simply supported on opposite ends
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
% C = MyAssembly.constrain_matrix(C);



%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
[V0,omega2] = eigs(K,M,n_VMs,'SM');
omega = sqrt(diag(omega2));

V = MyAssembly.unconstrain_vector(V0);
mod = 1;
v1 = reshape(V(:,mod),6,[]);
figure;
PlotFieldonDeformedMesh(nodes,elements,v1(1:3,:).','factor',5);
title(['Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )
set(colorbar,'visible','off')

mod = 2;
v1 = reshape(V(:,mod),6,[]);
figure;
PlotFieldonDeformedMesh(nodes,elements,v1(1:3,:).','factor',5);
title(['Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )
set(colorbar,'visible','off')

%% Damping matrix
disp('Using Rayleigh damping')
W =   omega(1:2);
a = [W(1) 1/W(1);W(2) 1/W(2)]\[0.004;0.004]
C = a(2) * M + a(1) * K;
MyAssembly.DATA.C = a(2) * MyAssembly.DATA.M + a(1) * MyAssembly.DATA.K;
%% external force assembly
disp('Assembling external force vector')
outcoord = [l/2,b/4]; % output coordinate
outdir = 3; % transverse displacement
dist = vecnorm(MyMesh.nodes(:,1:2) - repmat(outcoord,[MyMesh.nNodes,1]),2,2);
[~,outnode] = min(dist);
outdof = (outnode-1)*MyAssembly.Mesh.nDOFPerNode+outdir;

outdofvec = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
outdofvec = MyAssembly.constrain_vector(outdofvec);
outdof = find(outdofvec);

fext = 100*outdofvec;

centcoord = [l/2,b/2];
dist = vecnorm(MyMesh.nodes(:,1:2) - repmat(centcoord,[MyMesh.nNodes,1]),2,2);
[~,outnode] = min(dist);
centdof = (outnode-1)*MyAssembly.Mesh.nDOFPerNode+outdir;
centdofvec = sparse(centdof,ones(size(centdof)),1,MyMesh.nDOFs,1);
centdofvec = MyAssembly.constrain_vector(centdofvec);
centdof = find(centdofvec);
outdof = [outdof;centdof];


% fext = MyAssembly.constrain_vector(MyAssembly.uniform_body_force());

computationTimeLIN = toc(startLIN);
%% Tensor Assembly
disp('Getting nonlinearity coefficients')
filename = ['tensors_' num2str(MyMesh.nElements) '.mat'];
try     
    load(filename,'fnl')
    disp('Loaded tensors from storage')
    load(filename, 'computationTimeTensors')
    disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])
catch
    fnl = cell(1,2);
    disp('Assembling Tensors')
    startTensors = tic;
    fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
    fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);       
    computationTimeTensors = toc(startTensors);
    disp('Saving Tensors')
    save(filename,'fnl','computationTimeTensors','-v7.3')
    disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])
end

% apply boundary conditions
for j = 1:length(fnl)
    fnl{j} = MyAssembly.constrain_tensor(fnl{j});
end

