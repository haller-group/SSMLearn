function [M,C,K,fnl,fext, outdof, PlotFieldonDefMesh] = build_model(nElements)
%% Finite Element Setup
% Geometry
startLIN = tic;

l = 1;
h = 1e-3;
b = 1; 
% Mesh parameters

% Material properties

E       = 70e9;  % 70e9 % 200e9 % Young's modulus
rho     = 2700; % 2700 % 7850 % density
nu      = 0.3;    % nu
kappa   = 2e7; % material damping modulus 1e8

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor = @()BeamElement(b, h, myMaterial); % same element all across the domain

% Meshing the geometry
dx = l/nElements;
x = (0:dx:l).';
nNodes = size(x,1);
nodes = [x, zeros(nNodes,1)];
elements = [1:nNodes-1;2:nNodes].';

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
C = MyAssembly.damping_matrix();

%% apply boundary conditions
disp('Applying boundary conditions')
MyMesh.set_essential_boundary_condition([1 nNodes],[1 2 3],0) % Clamped-Clamped
% MyMesh.set_essential_boundary_condition([1],[1 2 3],0) % Cantilevered beam
% MyMesh.set_essential_boundary_condition([1 nNodes],[1 2],0) % Pinned-Pinned
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
C = MyAssembly.constrain_matrix(C);



%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
[V0,omega2] = eigs(K,M,n_VMs,'SM');
omega = sqrt(diag(omega2));

PlotFieldonDefMesh = @(defor,factor) PlotFieldonDeformedMesh_wrap(defor,nodes,elements,MyAssembly,factor);

Info = struct('nodes',nodes,'elements',elements,'assembly',MyAssembly);
figure; mod = 1; PlotFieldonDefMesh(V0(:,mod),0.2);
title(['Mode ' num2str(mod) ', Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )

%% external force assembly
disp('Assembling external force vector')

outnode = ceil(MyMesh.nNodes/2);
outdof = outnode*3-1; % transverse direction

outdofvec = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
outdofvec = MyAssembly.constrain_vector(outdofvec);
outdof = find(outdofvec);

fext = outdofvec;


% weights = true(nElements,1); 
% fext = MyAssembly.constrain_vector(MyAssembly.uniform_body_force('weights',weights));


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
end

function PlotFieldonDeformedMesh_wrap(deformation,nodes,elements,MyAssembly,factor)
V = MyAssembly.unconstrain_vector(deformation);
v1 = reshape(V,3,[]);
PlotFieldonDeformedMesh(nodes,elements,v1(1:2,:).','factor',factor);
end