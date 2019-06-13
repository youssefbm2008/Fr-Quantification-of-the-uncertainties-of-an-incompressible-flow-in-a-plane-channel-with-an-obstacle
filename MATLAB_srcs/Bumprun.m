function Bumprun(Nr,iofile,case_number,PP)
% 
% Bumprun    runs Nr realizations of the bump program in the appropriate
% directories
%
% Synopsis:  Bumprun(Nr,iofile,case_number,PP);
%
% Inputs:    Nr = total number of deterministic runs
%            iofile = name of the input file to the bump program
%            case_number = case under consideration
%            PP = Project path
% Output:    none
%
for k=1:Nr
    cd([PP,'/data','/CASE',num2str(case_number),'/RUNS/run',num2str(k),'/']);
    display('Current directory:'); pwd
   
    if (exist(iofile,'file')),
    execute_command = [PP,'/BUMP/BumpCode/bump ',iofile];
    system(execute_command);
    
    display(' ');
    display([' Run',num2str(k),' completed - Hit return or wait']);
    display('--');
    pause(3);
    else
        error('Bumprun.m: input file is missing!');
    end
end
