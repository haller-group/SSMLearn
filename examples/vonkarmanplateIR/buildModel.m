function [M,C,K,fnl,fext,outdofs, MyAssembly] = buildModel(l,b,t,w,E,rho,nu,kappa,bc,nl,nb)

%% Finite Element Setup
% % Geometry
% l  = 1; % length of domain [m]
% b = 2; % breadth of domain [m]
% t = 1e-2; % thickness of plate [m]
% w = 1e-1; % curvature parameter (height of the midpoint relative to ends) [m]
% % material properties
% E       = 70e9;  % 70e9 % 200e9 % Young's modulus [Pa]
% rho     = 2700; % 2700 % 7850 % density [kg/m^3]
% nu      = 0.33;    % Poisson's ratio 
% kappa   = 1e5; % material damping modulus 1e8

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor =  @()TriShellElement(t, myMaterial); % same element all across the domain

% Meshing the geometry
% nl = 10;
% nb = 20; 
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
MyAssembly.DATA.K = K;
MyAssembly.DATA.M = M;
% C = MyAssembly.damping_matrix();

%% apply boundary conditions
disp('Applying boundary conditions')
switch bc
    case 'SFSF'
        bcnodes = [bnodes{1}, bnodes{2}];
    case 'FSFS'
        bcnodes = [bnodes{3}, bnodes{4}];
    case 'SSSS'
        bcnodes = cell2mat(bnodes);
        bcnodes = unique(bcnodes);
    otherwise
        error('Does not support such bc types');
end
MyMesh.set_essential_boundary_condition(bcnodes,1:3,0) % simply supported on opposite ends
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
% C = MyAssembly.constrain_matrix(C);



%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
[V0,omega2] = eigs(K,M,n_VMs,'SM');
omega = sqrt(diag(omega2));
disp('The circular natural frequencies of the first five modes:');
omega

V = MyAssembly.unconstrain_vector(V0);
% mod = 2;
% v1 = reshape(V(:,mod),6,[]);
% figure(1);
% PlotFieldonDeformedMesh(nodes,elements,v1(1:3,:).','factor',5);
% title(['Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )
% view(2);
% figure(2);
% mod = 2;
% v1 = reshape(V(:,mod),6,[]);
% PlotFieldonDeformedMesh(nodes,elements,v1(1:3,:).','factor',5);
% title(['Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )
% view(3);
% figure(3);
% mod = 3;
% v1 = reshape(V(:,mod),6,[]);
% PlotFieldonDeformedMesh(nodes,elements,v1(1:3,:).','factor',5);
% title(['Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )
% view(3);


for mod = 1:4
    figure;
    v1 = reshape(V(:,mod),6,[]);
    PlotFieldonDeformedMesh(nodes,elements,v1(1:3,:).','factor',5);
    title(['Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )
    set(colorbar,'visible','off')
    colormap turbo
    view(3);
end
%% Damping matrix
disp('Using Rayleigh damping')
% W =   omega(1:2);
% a = [W(1) 1/W(1);W(2) 1/W(2)]\[0.004;0.004];
a = [4e-6 1];
C = a(2) * M + a(1) * K;
MyAssembly.DATA.C = a(2) * MyAssembly.DATA.M + a(1) * MyAssembly.DATA.K;


%% Tensor Assembly
disp('Getting nonlinearity coefficients')
fileName = ['tensors_',num2str(nl),'_',num2str(nb),'.mat'];
try 
    load(fileName,'fnl')    
    disp('Loaded tensors from storage')
catch
    fnl = cell(1,2);
    disp('Assembling Tensors')
    fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
    fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);
    disp('Saving Tensors')
    save(fileName,'fnl','-v7.3')
end

% apply boundary conditions
for j = 1:length(fnl)
    fnl{j} = MyAssembly.constrain_tensor(fnl{j});
end

%% external force assembly
disp('Assembling external force vector')
% outcoord = [l/2,b/4]; % output coordinate
outcoords = [0.2*l,0.3*l; 0.7*l 0.3*b; 0.25*l 0.25*b; 0.75*l 0.25*b];
outdir = 3; % transverse displacement
outdofs = zeros(4,1);
for idof = 1:4
    outcoord = outcoords(idof,:);
    dist = vecnorm(MyMesh.nodes(:,1:2) - repmat(outcoord,[MyMesh.nNodes,1]),2,2);
    [~,outnode] = min(dist);
    outdof = (outnode-1)*MyAssembly.Mesh.nDOFPerNode+outdir;

    outdofvec = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
    outdofvec = MyAssembly.constrain_vector(outdofvec);
    outdofs(idof) = find(outdofvec);
end
outdofs = unique(outdofs);

% fext = outdofvec;

fext = MyAssembly.constrain_vector(MyAssembly.uniform_body_force());