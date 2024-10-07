%Updated to first release version of PeGS2 by Carmen Lee 29/9/24
% Written for PeGS 2.0 by Ben MacMillan 8/24 
%Adapted from Carmen Lee's adaptation of Jonathan Kollmer's PeGS 1.0
%
% # diskSolve
% This module uses a nonlinear least-squares solver to match the photoelastic fringe pattern of each *experimental* particle to a *theoretical* fringe pattern. The best fit theoretical fringe pattern provides the magnitude and direction of each force acting on a particle.
% 
% **I**The input of this module is the standardised MATLAB **particle** structure array. Each particle in the image has the following fields:
% 
% - `id` : assigned ID number of the particle
% - `x`: x coordinate of the particle centre, in *pixels*
% - `y` : y coordinate of the particle centre, in *pixels*
% - `r` : radius of the particle, in *pixels*
% - `rm` : radius of the particle, in *metres*
% - `color` : assigned particle color (for plotting)
% - `fsigma` : photoelastic stress coefficient of the particle
% - `z` : number of contacts on the particle
% - `f` : average force on the particle calculated using the G<sup>2</sup> method
% - `g2` : G<sup>2</sup> value of the particle image
% - `forces` : array of contact force magnitudes 
% - `betas` : array of contact force azimuthal angles
% - `alphas` : array of contact force 'contact angles'. 0 is a purely normal force, $\pm \pi/2$ is a purely tangential force
% - `neighbours` : array containing particle IDs of contacting particles
% - `contactG2s` : array of G<sup>2</sup> values, corresponding to area around each contact force
% - `forceImage` : cropped image of particle from original experimental image
% 
% diskSolve requires all fields of **particle** to be populated *other than* `forces` and `alphas`. diskSolve also adds the following fields:
% 
% - `fitError` : error of the least-squares fit
% - `synthImg` : fitted theoretical fringe pattern
% 
% In addition to this, the following user inputs are required:
% 
% 
% - `algorithm` : algorithm used for the least-squares fit. Default set to `levenberg-marquardt`
% - `maxIterations` : maximum number of iterations for solver. Default set to 200
% - `functionTolerance` : error tolerance for solver. Default set to 0.01
% - `scaling` : scale factor to change size of experimental image for solver input. Default set to 0.5
% - `maskradius` : percentage of particle radius to use for fit. i.e. 1 - `maskradius` is fraction of radius that is discarded
% 
%  The other inputs can be tuned to improve accuracy and/or run time.
% 
% **O**
% The output of this module is an updated **particle** structure with the force values and a synthetic image of each particle. If verbose is selected, it will also save an image of the entire packing composed of the synthetic particles (good for troubleshooting and quality control)
% 
% ---
% 
% ## Usage
% To run this module, input the `directory` variable on line 6 of `diskSolve.m` and run the script. This will call the following functions:
% 
% - `solver_original.m` : calculates initial conditions for least-squares fit and performs fit
% - `forceBalance.m` : applies force balance to the particle. **UNMODIFIED FROM PEGSV1 VERSION**
% - `fringe_pattern_original.m` : creates fringe pattern for given force magnitudes and angles
% - `stress_engine_original.m` : principal stress calculation to provide pixel intensities for fringe pattern 
%  
% --- 
% ### PLEASE READ
% This version functions the same as PeGS Version 1.0. It includes C Lee's parallelisation modification and B McMillan's fix to the stress tensor calculation. The following features will be added soon:
% 
% - Vectorised version of the fringe pattern calculation to improve efficiency
% - Line contact solution for forces. Currently all forces are assumed to be point contacts
% - The ability to switch between original version, vectorised version, line contact and point contact solutions
% - The ability to run diskSolve on multiple video frames
% - The ability to turn on/off force balance requirement
% 
% The ability to watch the fringe pattern converge on screen has been removed due to the parallelisation requirements. Comparing the experimental fringes to theoretical ones should now be done in post-processing.
% 
% **I have not validated or modified the force balance calculations** 

function out = diskSolve(fileParams, dsParams, verbose)

%% Main function for diskSolve module


if ~exist(fullfile(fileParams.topDir, fileParams.solvedDir) , 'dir')
    mkdir(fullfile(fileParams.topDir, fileParams.solvedDir))
end


if verbose
    disp('starting disksolve() to find all particle centroids and save results in particleDir')
end

%goes to function to set up default disksolve parameters at the bottom of
%this file
dsParams = setupParams(dsParams, verbose);


%% Load particle data structure
% directory: particle data location
particledirectory = dir( fullfile(fileParams.topDir, fileParams.contactDir,'*_contacts.mat'))  ;
nFrames = length(particledirectory);


imgp = dir(fullfile(fileParams.topDir, fileParams.warpedImgDir, '*warped.tif'));
img = imread(fullfile(imgp(1).folder, imgp(1).name)); % we need this for setting up the synthetic image to show the forces at the end of this


if verbose
    disp(['starting disksolve() to fit ', num2str(nFrames), 'images'])
end

%% begin disksolving
for frame = 1 : nFrames
    data = load( fullfile(particledirectory(frame).folder , particledirectory(frame).name ));
    particle = data.particle;

    if dsParams.original == 1
        particle = solver_original(particle, dsParams);
    elseif dsParams.original == 0 && dsParams.vectorise == 1
        particle = solver_vectorised(particle, dsParams);
    end

    %% Save output
    savename = strrep(particledirectory(frame).name , '_contacts','_solved');
    save(fullfile(fileParams.topDir, fileParams.solvedDir , savename), 'particle')

    h3 = figure(1);
    hAx1 = subplot(1,1,1,'Parent', h3);
    NN = length(particle);
    bigSynthImg = zeros(size(img,1),size(img,2)); %make an empty image with the same size as the camera image
    for n=1:NN %for all particles
        %display(['fitting force(s) to particle ',num2str(n)]); %status indicator
        if (particle(n).z > 0 )
            %Add the syntetic peImage for the particle to the
            %synthetic image of our whole packing 
            x = floor(particle(n).x); %interger rounded x coordinate of the current particle
            y = floor(particle(n).y); %interger rounded y coordinate of the current particle
            sImg = particle(round(n)).synthImg; %synthetic force image for the current particle
            sx = size(sImg,1)/2; %width of the synthetic force image of the current particle
            sy = size(sImg,2)/2; %heights of the synthetic force image of the current partice
            bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx)) = bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx))+sImg; %Add the syntetic Force Image of the current particle to the appropriate location
            
        end
    
    end
   if verbose
    imshow(bigSynthImg, 'Parent', hAx1);
    
   end
   hold off;
    drawnow;
    savename = strrep(particledirectory(frame).name , '_contacts.mat','_Synth.jpg');
    
    imwrite(bigSynthImg, fullfile(fileParams.topDir, fileParams.synthImgDir, savename), "jpg"); %save synthetic image

end

%% Create parameter list and module wrap up
dsParams = rmfield(dsParams, 'fitoptions');
fields = fieldnames(dsParams);
for i = 1:length(fields)
    fileParams.(fields{i}) = dsParams.(fields{i});
end


fileParams.time = datetime("now");
fields = fieldnames(fileParams);
C=struct2cell(fileParams);
dsParams = [fields C];

writecell(dsParams,fullfile(fileParams.topDir, fileParams.solvedDir,'diskSolve_params.txt'),'Delimiter','tab');


if verbose 
        disp('done with diskSolve()');
end

out = true;


end



%% Set default parameters if no user input is given
function params = setupParams(params, verbose)

%%Least squares fit options
%Algorithm to use. Other options: 'trust-region-reflective', 'interior-point'
if ~isfield(params,'algorithm')
    params.algorithm = 'levenberg-marquardt';
end

% Function evalution limits
if ~isfield(params,'maxIterations')
    params.maxIterations = 200;
end
if ~isfield(params,'maxFunctionEvaluations')
    params.maxFunctionEvaluations = 400;
end
if ~isfield(params,'functionTolerance')
    params.functionTolerance = 0.01;
end

%% Scaling of image for fit
if ~isfield(params,'scaling')
    params.scaling = 0.5;
end

%% Masking radius
% How much of particle edge to remove before fit
if ~isfield(params,'maskradius')
    params.maskradius = 0.96;
end

%% Which version to use?
% Run original unvectorised version of disk solver
if ~isfield(params,'original')
    params.original = 1;
end
%Run vectorised version of disk solver (coming soon)
if ~isfield(params,'vectorise')
    params.vectorise = 0;
end



%% Add verbose input
% show fit results if verbose 
if verbose
    display = 'final-detailed';
else
    display = 'none';
end
fitoptions = optimoptions('lsqnonlin','Algorithm',params.algorithm,'MaxIter',params.maxIterations,'MaxFunEvals',params.maxFunctionEvaluations,'TolFun',params.functionTolerance,'Display',display);
params.fitoptions = fitoptions;
end