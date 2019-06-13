function WriteFiles(Nq,PP,case_number,iofile,nx,ny,x1t_qd,x2t_qd)
%
% WriteFiles writes and saves input files for the bump code in
% appropriately created directories
%
% Synopsis:  WriteFiles(Nq,PP,case_number,iofile,nx,ny,x1t_qd,x2t_qd);
%
% Inputs:    Nq = total number of quadrature points in 1D
%            PP = Project path 
%            case_number = case under consideration
%            iofile = name of the input file to the bump program
%            nx = number of grid points in the streamwise direction
%            ny = number of grid points in the crossflow direction
%            x1t_qd,x2t_qd = input quadrature values
% Output:    none
%
% Remark:   only writes two-dimensional data with a full tensor parametric structure

run_number = 0;

for i=1:Nq,
    for j=1:Nq,
        run_number = run_number + 1;
        dirname = [PP,'/data','/CASE',num2str(case_number),'/RUNS/run',num2str(run_number),'/'];
        if (~exist(dirname,'dir')), mkdir(dirname); end
        Re = x1t_qd(i); BH = x2t_qd(j);
        WriteInputParameters(iofile,nx,ny,Re,BH);
        movefile(iofile,dirname);
        save data_quad Re BH nx ny;
        movefile('data_quad.mat',dirname);
    end
end