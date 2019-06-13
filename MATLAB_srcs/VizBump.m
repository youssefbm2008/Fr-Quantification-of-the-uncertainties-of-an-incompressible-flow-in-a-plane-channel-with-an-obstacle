function [XY, UV] = VizBump(varargin)
%
% VizBump   reads and displays the U- & V-velocity fields produced by the bump code
%            and outputs the X & Y coordinates and the U & V velocity as
%            matrices
% Synopsis:  [XY, UV] = VizBump(run_idnumber,case_number,PP);
%
% Inputs (optional):
%            run_idnumber = id number of the run to visualize
%            case_number = case under consideration
%            PP = Project path
% Output:    XY = coordinates matrix of size [2*nx-1,2*ny-1,2]
%            UV = velocity matrix of size [2*nx-1,2*ny-1,2]
%

%% Files names

if nargin>0
    run_idnumber = varargin{1}; case_number = varargin{2}; PP = varargin{3};
    filepath = strcat(PP,'/data/CASE',num2str(case_number),'/RUNS/run',num2str(run_idnumber),'/');
else
    filepath = './';
end

connectivity_file = strcat(filepath,'connect_001.txt');
coordinates_file = strcat(filepath,'xy_001.txt');
velocity_file = strcat(filepath,'uv_001.txt');
data_file = strcat(filepath,'data_quad.mat');

%% Files loading

if (exist(connectivity_file,'file'))
load(connectivity_file);
else error('Viz_Bump.m: sorry this file does not exist'); end

if (exist(coordinates_file,'file'))
load(coordinates_file);
else error('Viz_Bump.m: sorry this file does not exist'); end

if (exist(velocity_file,'file'))
load(velocity_file);
else error('Viz_Bump.m: sorry this file does not exist'); end

if (exist(data_file,'file'))
load(data_file); % this small data file is used to load the correct labels in the figures
else Re = 30; BH = 1; nx = 100; ny = 50; end

TRI6 = connect_001; TRI3 = TRI6(:,2:4); clear connect_001; clear TRI6;
X6 = xy_001(:,1); Y6 = xy_001(:,2); clear xy_001;
U6 = uv_001(:,1); V6 = uv_001(:,2); clear uv_001;

XY = zeros(2*nx-1,2*ny-1,2); XY(:,:,1) = reshape(X6,2*ny-1,2*nx-1)'; XY(:,:,2) = reshape(Y6,2*ny-1,2*nx-1)';
UV = zeros(2*nx-1,2*ny-1,2); UV(:,:,1) = reshape(U6,2*ny-1,2*nx-1)'; UV(:,:,2) = reshape(V6,2*ny-1,2*nx-1)';

ftname = 'Helvetica'; ftsize = 9;

figure;
subplot(2,1,1); trisurf(TRI3,X6,Y6,U6); view(0,89.99); colorbar; axis equal;  shading interp; %flat;
title('U-velocity','fontname','Helvetica','fontsize',ftsize);
xlabel('x','fontname',ftname,'fontsize',ftsize); ylabel('y','fontname',ftname,'fontsize',ftsize); 
axis tight; box on;

subplot(2,1,2); trisurf(TRI3,X6,Y6,V6); view(0,89.99); colorbar; axis equal;  shading interp; %flat;
title('V-velocity','fontname','Helvetica','fontsize',ftsize);
xlabel('x','fontname',ftname,'fontsize',ftsize); ylabel('y','fontname',ftname,'fontsize',ftsize);
axis tight; box on;

text(0,3.5,['Reynolds = ',num2str(Re),'; Bump height = ',num2str(BH),'; nx = ',num2str(nx),'; ny = ',num2str(ny),' '],'fontname',ftname,'fontsize',ftsize);