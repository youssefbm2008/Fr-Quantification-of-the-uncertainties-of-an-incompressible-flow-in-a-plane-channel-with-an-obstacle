function [data] = Dataread(Nq,PP,case_number,nx,ny)
%
% Dataread   reads and collects data output files generated from the bump code in the
% appropriate directories
%
% Synopsis:  [data] = Dataread(Nq,PP,case_number,nx,ny);
%
% Inputs:    Nq = total number of quadrature points in 1D
%            PP = Project path 
%            case_number = case under consideration
%            nx = number of grid points in the streamwise direction
%            ny = number of grid points in the crossflow direction
% Output:    data = concatenated UV velocity data from all runs
%

data = zeros(2*nx-1,2*ny-1,Nq*Nq,2); run_number = 0;

for i=1:Nq,
    for j=1:Nq,
        run_number = run_number + 1;
        dirname = [PP,'/data','/CASE',num2str(case_number),'/RUNS/run',num2str(run_number),'/'];
        cd(dirname);
        pwd
        load uv_001.txt;
        data(:,:,run_number,1) = reshape(uv_001(:,1),2*ny-1,2*nx-1)';
        data(:,:,run_number,2) = reshape(uv_001(:,2),2*ny-1,2*nx-1)';
        clear uv_001;
    end
end
cd(PP);