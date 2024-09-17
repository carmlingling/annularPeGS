function out = diskSolve(p, f, verbose)

%% Main function for diskSolve module

%% Set default parameters if no user input is given
%% Least squares fit options
%Algorithm to use. Other options: 'trust-region-reflective', 'interior-point'
if ~isfield(f,'algorithm')
    f.algorithm = 'levenberg-marquardt';
end

% Function evalution limits
if ~isfield(f,'maxIterations')
    f.maxIterations = 200;
end
if ~isfield(f,'maxFunctionEvaluations')
    f.maxFunctionEvaluations = 400;
end
if ~isfield(f,'functionTolerance')
    f.functionTolerance = 0.01;
end

%% Scaling of image for fit
if ~isfield(f,'scaling')
    f.scaling = 0.5;
end

%% Masking radius
% How much of particle edge to remove before fit
if ~isfield(f,'maskradius')
    f.maskradius = 0.96;
end

%% Which version to use?
% Run original unvectorised version of disk solver
if ~isfield(f,'original')
    f.original = 1;
end
%Run vectorised version of disk solver (coming soon)
if ~isfield(f,'vectorise')
    f.vectorise = 0;
end

%% Create parameter list
parameterlist = [fieldnames(f) struct2cell(f)];

%% Add verbose input
% show fit results if verbose 
if verbose
    display = 'final-detailed';
else
    display = 'none';
end
fitoptions = optimoptions('lsqnonlin','Algorithm',f.algorithm,'MaxIter',f.maxIterations,'MaxFunEvals',f.maxFunctionEvaluations,'TolFun',f.functionTolerance,'Display',display);
f.fitoptions = fitoptions;

%% Load particle data structure
% directory: particle data location
particledirectory = dir( fullfile(p.topDir,'particle/','*_preprocessing.mat'))  ;
noparticlestructs = length(particledirectory);
for frame = 1 : noparticlestructs
    data = load( fullfile(p.topDir , particledirectory(frame).name ));
    particle = data.particle;

    if f.original == 1
        particle = solver_original(particle, f);
    elseif f.original == 0 && f.vectorise == 1
        particle = solver_vectorised(particle, f);
    end

    %% Save output
    savename = strrep(particledirectory(frame).name , '_preprocessing','_solved');
    save(fullfile(topDirectory.outDir , savename), 'particle')

end

%% Save user input parameter list
cd(topDirectory.outDir)
writecell(parameterlist,'diskSolve_params.txt','Delimiter','tab')

%% Output variable if everything is finished
out = 1;


