function WriteInputParameters(sfile,nx,ny,Re,bs)
% 
% WriteInputParameters     generates an ascii input file with grid,
%            Reynolds number and bump height information
%
% Synopsis:  WriteInputParameters(sfile,nx,ny,Re,bs);
%
% Inputs:    sfile = name of the file generated
%            nx = number of grid points in the streamwise direction
%            ny = number of grid points in the crossflow direction
%            Re = value of the Reynolds number
%            bs = bump size in the vertical direction
% Output:    none
%
	fprintf('\n writing input parameters [file=%s]....',sfile)
	fid = fopen(sfile,'w');
	fprintf(fid,'%s \n','#  -------------------------------------------------');
	fprintf(fid,'%s \n','#      INPUT PARAMETER FILE FOR THE BUMPCODE');
	fprintf(fid,'%s \n','#  -------------------------------------------------');
	fprintf(fid,'%s \n','#Inputs');
	fprintf(fid,'\t %s \t %u \n','Nx',nx);
	fprintf(fid,'\t %s \t %u \n','Ny',ny);
	fprintf(fid,'\t %s \t %u \n','Reynolds',Re);
	fprintf(fid,'\t %s \t %u \n','Bumpsize',bs);
	fprintf(fid,'%s \n','#END_Inputs');
	fprintf(fid,'%s \n','#  -------------------------------------------------');
	fprintf(fid,'%s \n','#END_OF_FILE');
	fprintf(fid,'%s \n','#  -------------------------------------------------');    
	fclose(fid);
	fprintf(' [done]\n')
