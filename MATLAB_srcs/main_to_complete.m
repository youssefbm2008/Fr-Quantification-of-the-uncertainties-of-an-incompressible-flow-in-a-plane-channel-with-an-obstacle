clear all; % clears the workspace (removes all variables, etc...)
close all; % close all figures
clc;       % clears the command window

%%% PROBLEM SETUP (DATA, FILES,...) %%%

% définition des paramètres physiques - sources d'incertitude %%%
%-----------------------%

% Mean & std of Reynolds number
Re_mean = 50; 
Re_std = 10; 

% Mean & std of bump height
bump_height_mean = 0.5; 
bump_height_std = 0.1; 
%-----------------------%

% choisir la méthode à utiliser 'projection' 'montecarlo' ou 'regression'
%method = 'projection'; 

% l'utilisateur définit la méthode à chaque exécution
prompt = 'Saisir la méthode à utiliser: projection, montecarlo ou regression: ';
method=input(prompt,'s');





% Numerical parameters %%%
%------------------------%
Nx = 60; Ny = 30; % points pour la discrétisation spatiale

%Montecarlo
nnc = 5; 

%gPC
N = 2; % number of random parameters/dimensions
P = 3; % Polynomial chaos order
M = PCnumbterms(P,N) ;% use PCnumbterms
Nq_1D = P+1 ; % Number of quadrature grid points per random dimension
n1 = 13; n2 = 13; Ne = n1.*n2; % Number of (response surface) sampling points for VISUALIZATION purposes
Nq_ND = Nq_1D^N ;
%------------------------%

% On définitr les 3 méthodes de résolution

switch method
    
%------------------------    
% APPROCHE PAR PROJECTION
%------------------------

    case 'projection'




% Code and data structure parameters %%%
%--------------------------------------%
case_number = 0; % should be defined to differentiate among different cases
Project_path = '/auto/b/benmahmoud/TP1_AMS_V05/TP_Matlab'; % path to the location of your project

run_number = 0; % run number starts at 0
cd(Project_path);
addpath([Project_path,'/MATLAB_srcs']);  addpath([Project_path,'/BUMP']);
input_file = 'bump.inp'; % input file name
%--------------------------------------%

% Memory allocation %%%
%---------------------%
x1 = linspace(-1,1,n1); x2 = linspace(-1,1,n2); [X1,X2] = ndgrid(x1,x2);
x1t = linspace(Re_mean-Re_std,Re_mean+Re_std,n1);
x2t = linspace(bump_height_mean-bump_height_std,bump_height_mean+bump_height_std,n2);

gPC = struct('coeff',zeros(2*Nx-1,2*Ny-1,M),'data',zeros(2*Nx-1,2*Ny-1,Nq_1D*Nq_1D),...
    'Var',zeros(M,1),'data_viz',zeros(n1,n2),'poly',zeros(Ne,M));

U_gPC = struct('coeff',zeros(2*Nx-1,2*Ny-1,M),'mean', zeros(2*Nx-1,2*Ny-1),'var',zeros(2*Nx-1,2*Ny-1),...
    'data',zeros(2*Nx-1,2*Ny-1,Nq_1D^2),'cov',zeros(2*Nx-1,2*Ny-1));
V_gPC = struct('coeff',zeros(2*Nx-1,2*Ny-1,M),'mean', zeros(2*Nx-1,2*Ny-1),'var',zeros(2*Nx-1,2*Ny-1),...
    'data',zeros(2*Nx-1,2*Ny-1,Nq_1D^2),'cov',zeros(2*Nx-1,2*Ny-1));


%---------------------%

% Construction of the 1D Gauss-Legendre quadrature grid %%%
%---------------------------------------------------------%
[x1_qd,w1_qd] = GaussLegPtsPds(Nq_1D);% use GaussLegPtsPds
[x2_qd,w2_qd] = GaussLegPtsPds(Nq_1D); % use GaussLegPtsPds
x1t_qd = Mapping(x1_qd,Re_mean-Re_std,Re_mean+Re_std);% use Mapping
x2t_qd = Mapping(x2_qd,bump_height_mean-bump_height_std,bump_height_mean+bump_height_std);% use Mapping
%---------------------------------------------------------%

% Tensorization to get the N-dimensional quadrature grids %%%
%-----------------------------------------------------------%
[X1_qd,X2_qd] = ndgrid(x1_qd,x2_qd);% use ndgrid
[X1t_qd,X2t_qd] = ndgrid(x1t_qd,x2t_qd); % idem
%-----------------------------------------------------------%

% Creation of all input files for deterministic runs %%%
%------------------------------------------------------%
% use WriteFiles
WriteFiles(Nq_1D,Project_path,case_number,input_file,Nx,Ny,x1t_qd,x2t_qd);
display(' Directories and files are ready - Hit return or wait');
display('--'); pause(2)
%------------------------------------------------------%


%% RUNNING SIMULATIONS %%%

% Running of all deterministic simulations %%%
%--------------------------------------------%
% use Bumprun
Nr = Nq_1D^2 ; 
Bumprun(Nr,input_file,case_number,Project_path); %commenter, et lire
%directement quand on travaille sur les memes calculs. /!\ décommenter pour
%régression
%--------------------------------------------%


%% POLYNOMIAL CHAOS POST-PROCESSING %%%

cd(Project_path); pwd % moves back to right directory

% Reading the output data files %%%
%---------------------------------%
data = Dataread(Nq_1D,Project_path,case_number,Nx,Ny) ; % use ReadData
U_gPC.data = squeeze(data(:,:,:,1)); V_gPC.data = squeeze(data(:,:,:,2));
clear data;
%---------------------------------%


% Performing the projection to get the gPC modes %%%
%---------------------------------%
% use Legendre_poly; PCmultindex;
%---------------------------------%

mI = zeros(M,N); 
for k=0:M-1
    mI(k+1, :) = PCmultindex(N,P,k);
end 

for id_k =0:M-1
    id_run = 1; 
    for i1=1:Nq_1D
        for i2=1:Nq_1D
            U_gPC.coeff(:,:,id_k+1) =  U_gPC.coeff(:,:, id_k + 1) + U_gPC.data(:,:, id_run).* Legendre_poly([x1_qd(i1) x2_qd(i2)] ,N, id_k, mI)*w1_qd(i1)*w2_qd(i2);
            V_gPC.coeff(:,:,id_k+1) =  V_gPC.coeff(:,:, id_k + 1) + V_gPC.data(:,:, id_run).* Legendre_poly([x1_qd(i1) x2_qd(i2)] ,N, id_k, mI)*w1_qd(i1)*w2_qd(i2);
            id_run = id_run + 1; 
        end
    end

end 


% Computing the mean solution %%%
%---------------------------------%
% use definition from the class notes

U_gPC.mean(:,:) = U_gPC.coeff(:,:,1) ;
V_gPC.mean(:,:) = V_gPC.coeff(:,:,1) ;
%---------------------------------%


% Computing the std solution %%%
%---------------------------------%
% use definition from the class notes
for m=2:M
    U_gPC.var(:,:) = U_gPC.var(:,:) + U_gPC.coeff(:,:,m).^2 ; 
    V_gPC.var(:,:) = V_gPC.var(:,:) + V_gPC.coeff(:,:,m).^2 ; 
end 
%---------------------------------%

%computing the variation coefficient
U_gPC.cov = ((U_gPC.var).^(1/2))./(U_gPC.mean) ;
V_gPC.cov = ((V_gPC.var).^(1/2))./(V_gPC.mean) ;



% Getting the mesh coordinates to be able to plot the stochastic fields %%%
%---------------------------------%
meshfile = [Project_path,'/data','/CASE',num2str(case_number),'/RUNS/run',num2str(ceil((Nq_1D.^2)./2)-mod(Nq_1D+1,2)),'/xy_001.txt'];
[XmeshRep,YmeshRep] = GetMeshCoords(meshfile,Nx,Ny);

%représentation des champs moyens
figure;
subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), U_gPC.mean); colorbar; axis equal;  shading interp; %flat;
title('mean U','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;

subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), V_gPC.mean); colorbar; axis equal;  shading interp; %flat;
title('mean V','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;

%représentation des champs de vitesse de std de U
figure;
subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), sqrt(U_gPC.var)); colorbar; axis equal;  shading interp; %flat;
title('std U','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;


%représentation des champs de vitesse de std de V

subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), sqrt(V_gPC.var)); colorbar; axis equal;  shading interp; %flat;
title('std V','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;




%représentation des champs de coefficients de variation
figure;
subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), U_gPC.cov); colorbar; axis equal;  shading interp; %flat;
title('coV U','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;

subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), V_gPC.cov); colorbar; axis equal;  shading interp; %flat;
title('coV V','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;

% distributions spatiales des coeff gPC
%---------------------------------%
gPC_num = [3 5 7 10 13 15];
figure;
for gpc_num=1:6
    subplot(3,2,gpc_num); pcolor( XmeshRep(:,:),YmeshRep(:,:), U_gPC.coeff(:,:,gpc_num) ); colorbar; axis equal;  shading interp; %flat;
    title(['gPC coeff k=',num2str(gPC_num(gpc_num))],'fontname','Times','fontsize',14);
    xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
    axis tight; box on;
end 

%------------------------    
% MONTE CARLO
%------------------------



case 'montecarlo'


% sampling the parameters with respect to uniform density
re_test = (Re_mean - sqrt(3)*Re_std) + sqrt(12)*Re_std*rand(nnc,1);
bump_test = (bump_height_mean - sqrt(3)*bump_height_std) + sqrt(12)*bump_height_std.*rand(nnc,1); 

[Re_test,Bump_test] = ndgrid( re_test, bump_test); % use ndgrid
%---------------------------------------------------------%

case_number = 1; % should be defined to differentiate among different cases
Project_path = '/auto/b/benmahmoud/TP1_AMS_V05/TP_Matlab'; % path to the location of your project

run_number = 0; % run number starts at 0
cd(Project_path);
addpath([Project_path,'/MATLAB_srcs']);  addpath([Project_path,'/BUMP']);
input_file = 'bump.inp'; % input file name
%--------------------------------------%

% Memory allocation %%%
%---------------------%

U = struct('data',zeros(2*Nx-1,2*Ny-1,nnc^2), 'mean',zeros(2*Nx-1,2*Ny-1), 'var',zeros(2*Nx-1,2*Ny-1),'cov',zeros(2*Nx-1,2*Ny-1));
V = struct('data',zeros(2*Nx-1,2*Ny-1,nnc^2), 'mean',zeros(2*Nx-1,2*Ny-1), 'var',zeros(2*Nx-1,2*Ny-1),'cov',zeros(2*Nx-1,2*Ny-1));

%---------------------%

% Creation of all input files for deterministic runs %%%
%------------------------------------------------------%
% use WriteFiles
WriteFiles(nnc, Project_path, case_number, input_file, Nx, Ny, re_test, bump_test  );
display(' Directories and files are ready - Hit return or wait');
display('--'); pause(2)
%------------------------------------------------------%


%% RUNNING SIMULATIONS %%%

% Running of all deterministic simulations %%%
%--------------------------------------------%
% use Bumprun
Bumprun(nnc^2, input_file, case_number, Project_path);
%--------------------------------------------%


%% MC POST-PROCESSING %%%

cd(Project_path); pwd % moves back to right directory

% Reading the output data files %%%
%---------------------------------%
data = Dataread(nnc, Project_path, case_number, Nx, Ny) ; % use ReadData
U.data = squeeze(data(:,:,:,1)); V.data = squeeze(data(:,:,:,2));
clear data;
%---------------------------------%


% Computing the mean solution %%%
%---------------------------------%
% use definition from the class notes
%---------------------------------%
for run=1:nnc^2
   U.mean(:,:)=U.mean(:,:)+U.data(:,:,run) ;
   V.mean(:,:)=V.mean(:,:)+V.data(:,:,run) ; 
end

U.mean = 1/nnc^2*U.mean ; 
V.mean = 1/nnc^2*V.mean ; 

% Computing the std solution %%%
%---------------------------------%
% variance
for m=1:nnc^2
    U.var(:,:) = U.var(:,:) + (U.data(:,:,run)-U.mean(:,:)).^2 ; 
    V.var(:,:) = V.var(:,:) + (V.data(:,:,run)-V.mean(:,:)).^2 ; 
end 
U.var=(1/(nnc^2-1))*U.var ;
V.var=(1/(nnc^2-1))*V.var ;

%computing the variation coefficient
U.cov = ((U.var).^(1/2))./(U.mean) ;
V.cov = ((V.var).^(1/2))./(V.mean) ;

%---------------------------------%
% Generating response surfaces and solution pdfs %%%
%---------------------------------%
% use PCSampling; hist

%---------------------------------%


% Getting the mesh coordinates to be able to plot the stochastic fields %%%
%---------------------------------%
meshfile = [Project_path,'/data','/CASE',num2str(case_number),'/RUNS/run',num2str(ceil((nnc.^2)./2)-mod(nnc+1,2)),'/xy_001.txt'];
[XmeshRep,YmeshRep] = GetMeshCoords(meshfile,Nx,Ny);


%représentation des champs moyens 
figure;
subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), U.mean); colorbar; axis equal;  shading interp; %flat;
title('mean U','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;

subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), V.mean); colorbar; axis equal;  shading interp; %flat;
title('mean V','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;

%représentation des champs de std
figure;
subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), sqrt(U.var)); colorbar; axis equal;  shading interp; %flat;
title('std U','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;

subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), sqrt(V.var)); colorbar; axis equal;  shading interp; %flat;
title('std V','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;




%représentation des champs de coefficients de variation
figure;
subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), U.cov); colorbar; axis equal;  shading interp; %flat;
title('coV U','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;

subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), V.cov); colorbar; axis equal;  shading interp; %flat;
title('coV V','fontname','Times','fontsize',14);
xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
axis tight; box on;


%------------------------    
% APPROCHE PAR REGRESSION
%------------------------


case 'regression'
    

    % Code and data structure parameters %%%
    %--------------------------------------%
    case_number = 2; % should be defined to differentiate among different cases
    Project_path = '/auto/b/benmahmoud/TP1_AMS_V05/TP_Matlab'; % path to the location of your project

    run_number = 0; % run number starts at 0
    cd(Project_path);
    addpath([Project_path,'/MATLAB_srcs']);  addpath([Project_path,'/BUMP']);
    input_file = 'bump.inp'; % input file name
    %--------------------------------------%

    % Memory allocation %%%
    %---------------------%
    
    U_reg = struct('coeff',zeros(2*Nx-1,2*Ny-1,M),'mean', zeros(2*Nx-1,2*Ny-1),'var',zeros(2*Nx-1,2*Ny-1),...
        'data',zeros(2*Nx-1,2*Ny-1,Nq_1D^2),'cov',zeros(2*Nx-1,2*Ny-1));
    V_reg = struct('coeff',zeros(2*Nx-1,2*Ny-1,M),'mean', zeros(2*Nx-1,2*Ny-1),'var',zeros(2*Nx-1,2*Ny-1),...
        'data',zeros(2*Nx-1,2*Ny-1,Nq_1D^2),'cov',zeros(2*Nx-1,2*Ny-1));
    
    Phi_reg = zeros(Nq_ND, M); 

    %---------------------%

    % Construction des paramètres %
    %---------------------------------------------------------%
    re_test = (Re_mean - sqrt(3)*Re_std) + sqrt(12)*Re_std.*rand(Nq_ND,1);
    bump_test = (bump_height_mean - sqrt(3)*bump_height_std) + sqrt(12)*bump_height_std.*rand(Nq_ND,1); 
    
    %invMapping paramètres entre -1 et 1 
    re_un = InvMapping(re_test,Re_mean-Re_std,Re_mean+Re_std) ; 
    bump_un = InvMapping(bump_test,bump_height_mean-bump_height_std,bump_height_mean+bump_height_std) ;
        
    %-----------------------------------------------------------%

    % Creation of all input files for deterministic runs %%%
    %------------------------------------------------------%
    % use WriteFiles
    WriteFiles(Nq_1D,Project_path,case_number,input_file,Nx,Ny,re_test,bump_test);
    display(' Directories and files are ready - Hit return or wait');
    display('--'); pause(2)
    %------------------------------------------------------%


    %% RUNNING SIMULATIONS %%%

    % Running of all deterministic simulations %%%
    %--------------------------------------------%
    % use Bumprun
    Nr = Nq_ND ; 
    Bumprun(Nr,input_file,case_number,Project_path); 
    
    
    %--------------------------------------------%


    %% POLYNOMIAL CHAOS POST-PROCESSING %%%

    cd(Project_path); pwd % moves back to right directory

    % Reading the output data files %%%
    %---------------------------------%
    data = Dataread(Nq_1D,Project_path,case_number,Nx,Ny) ; % use ReadData
    U_reg.data = squeeze(data(:,:,:,1)); V_reg.data = squeeze(data(:,:,:,2));
    clear data;
    %---------------------------------%


    % Performing the projection to get the gPC modes %%%
    %---------------------------------%
    % use Legendre_poly; PCmultindex;
    %---------------------------------%
    mI = zeros(M,N); 
    for k=0:M-1
        mI(k+1, :) = PCmultindex(N,P,k);
    end 

    %construction de Phi
    for id_k =1:M
        id_run = 1; 
        for i1=1:Nq_1D
            for i2=1:Nq_1D
                Phi_reg(id_run,id_k) = Legendre_poly([re_un(i1) bump_un(i2)] ,N, id_k-1, mI);
                id_run = id_run + 1 ;
            end
        end
    end 

%construction des coeffs
    for x=1:2*Nx-1
        for y=1:2*Ny-1
            U = reshape(U_reg.data(x,y, :), [Nq_ND 1]); 
            V = reshape(V_reg.data(x,y, :), [Nq_ND 1]);
            U_reg.coeff(x,y, :) = (Phi_reg'*Phi_reg)\(Phi_reg'*U); 
            V_reg.coeff(x,y, :) = (Phi_reg'*Phi_reg)\(Phi_reg'*V); 
        end 
    end 


    % Computing the mean solution %%%
    %---------------------------------%
    % use definition from the class notes

    U_reg.mean(:,:) = U_reg.coeff(:,:,1) ;
    V_reg.mean(:,:) = V_reg.coeff(:,:,1) ;
    %---------------------------------%


    % Computing the std solution %%%
    %---------------------------------%
    % use definition from the class notes
    for m=2:M
        U_reg.var(:,:) = U_reg.var(:,:) + U_reg.coeff(:,:,m).^2 ; 
        V_reg.var(:,:) = V_reg.var(:,:) + V_reg.coeff(:,:,m).^2 ; 
    end 
    %---------------------------------%
    
    %computing the variation coefficient
    U_reg.cov = ((U_reg.var).^(1/2))./(U_reg.mean) ;
    V_reg.cov = ((V_reg.var).^(1/2))./(V_reg.mean) ; 

    


    % Getting the mesh coordinates to be able to plot the stochastic fields %%%
    %---------------------------------%
    meshfile = [Project_path,'/data','/CASE',num2str(case_number),'/RUNS/run',num2str(ceil((Nq_1D.^2)./2)-mod(Nq_1D+1,2)),'/xy_001.txt'];
    [XmeshRep,YmeshRep] = GetMeshCoords(meshfile,Nx,Ny);

    %représentation des champs moyens
    figure;
    subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), U_reg.mean); colorbar; axis equal;  shading interp; %flat;
    title('mean U','fontname','Times','fontsize',14);
    xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
    axis tight; box on;

    subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), V_reg.mean); colorbar; axis equal;  shading interp; %flat;
    title('mean V','fontname','Times','fontsize',14);
    xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
    axis tight; box on;
    
    %représentation des champs de std
    figure;
    subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), sqrt(U_reg.var)); colorbar; axis equal;  shading interp; %flat;
    title('std U','fontname','Times','fontsize',14);
    xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
    axis tight; box on;

    subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), sqrt(V_reg.var)); colorbar; axis equal;  shading interp; %flat;
    title('std V','fontname','Times','fontsize',14);
    xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
    axis tight; box on;
    
    
    
    %représentation des champs de coefficients de variation
    figure;
    subplot(2,1,1); pcolor(XmeshRep(:,:),YmeshRep(:,:), U_reg.cov); colorbar; axis equal;  shading interp; %flat;
    title('coV U','fontname','Times','fontsize',14);
    xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
    axis tight; box on;

    subplot(2,1,2); pcolor(XmeshRep(:,:),YmeshRep(:,:), V_reg.cov); colorbar; axis equal;  shading interp; %flat;
    title('coV V','fontname','Times','fontsize',14);
    xlabel('x','fontname','Times','fontsize',14); ylabel('y','fontname','Times','fontsize',14);
    axis tight; box on;

    

    % distributions spatiales coeff gPC
    %---------------------------------%
    Reg_num = [3 5 7 10 13 15];
    figure;
    for reg_num=1:6
        subplot(3,2,reg_num); pcolor( XmeshRep(:,:),YmeshRep(:,:), U_reg.coeff(:,:,reg_num) ); colorbar; axis equal;  shading interp; %flat;
        title(['reg coeff k=',num2str(Reg_num(reg_num))],'fontname','Times','fontsize',14);
        xlabel('x','fontname','Times','fontsize',13); ylabel('y','fontname','Times','fontsize',14);
        axis tight; box on;
    end 
  
end
