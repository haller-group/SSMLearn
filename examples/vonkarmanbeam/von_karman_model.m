function [M,C,K,fnl,fext,outdof] = von_karman_model(nElements, E, rho, nu, kappa, l, h, b)
%% Finite Element Setup
% Geometry
startLIN = tic;

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
% figure('Name','Mesh'); PlotMesh(nodes,elements,0);

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
MyMesh.set_essential_boundary_condition([1, nNodes],[1 2 3],0) % Clamped-clamped beam
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
C = MyAssembly.constrain_matrix(C);


%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 3; % first n_VMs modes with lowest frequency calculated 
[V0,omega2] = eigs(K,M,n_VMs,'SM');
omega = sqrt(diag(omega2));

V = MyAssembly.unconstrain_vector(V0);
mod = 1;
v1 = reshape(V(:,mod),3,[]);
figure;
PlotFieldonDeformedMesh(nodes,elements,v1(1:2,:).','factor',0.2);
title(['Mode ' num2str(mod) ', Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )

%% external force assembly
% disp('Assembling external force vector')

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
%     disp('Saving Tensors')
%     save(filename,'fnl','computationTimeTensors','-v7.3')
    disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])
end

% apply boundary conditions
for j = 1:length(fnl)
    fnl{j} = MyAssembly.constrain_tensor(fnl{j});
end
