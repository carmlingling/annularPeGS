function out = diskSolve(p, f, verbose)

%% Main function for diskSolve module
%% file housekeeping
if not(isfolder(append(p.topDir,p.solvedDir))) %make a new folder with particle centers
        mkdir(append(p.topDir,p.solvedDir));
end
if not(isfolder(append(p.topDir,p.synthImgDir))) %make a new folder with particle centers
        mkdir(append(p.topDir,p.synthImgDir));
end
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
particledirectory = dir( fullfile(p.topDir,p.contactDir,'*_contacts.mat'))  ;
noparticlestructs = length(particledirectory);


 imgp = dir(fullfile(p.topDir, p.warpedImgDir, '*.tif'));
 img = imread(fullfile(imgp(1).folder, imgp(1).name));
 bigSynthImg = zeros(size(img,1),size(img,2)); %make an empty image with the same size as the camera image
    

for frame = 1 : noparticlestructs
    data = load( fullfile(particledirectory(frame).folder , particledirectory(frame).name ));
    particle = data.particle;

    if f.original == 1
        particle = solver_original(particle, f);
    elseif f.original == 0 && f.vectorise == 1
        particle = solver_vectorised(particle, f);
    end

    %% Save output
    savename = strrep(particledirectory(frame).name , '_contacts','_solved');
    save(fullfile(p.topDir, p.solvedDir , savename), 'particle')
    
    h3 = figure(1);
    hAx1 = subplot(1,1,1,'Parent', h3);
    NN = length(particle);
    
    for n=1:NN %for all particles
        %display(['fitting force(s) to particle ',num2str(n)]); %status indicator
        if (particle(n).z > 0 )
            %Add the syntetic peImage for the particle to the
            %synthetic image of our whole packing 
            x = floor(particle(n).x); %interger rounded x coordinate of the current particle
            y = floor(particle(n).y); %interger rounded y coordinate of the current particle
            r = particle(n).r; %radius (in pixels) of the current particle
            sImg = particle(round(n)).synthImg; %synthetic force image for the current particle
            sx = size(sImg,1)/2; %width of the synthetic force image of the current particle
            sy = size(sImg,2)/2; %heights of the synthetic force image of the current partice
            bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx)) = bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx))+sImg; %Add the syntetic Force Image of the current particle to the appropriate location
            
        end
    
    end
   if verbose
    imshow(bigSynthImg, 'Parent', hAx1);
%     hold (hAx1, 'on');
%     for n=1:NN
%         viscircles([particle(n).x, particle(n).y], particle(n).r)
%         text(particle(n).x, particle(n).y, num2str(particle(n).fitError))
%     end     
   end
    drawnow;
    savename = strrep(particledirectory(frame).name , '_contacts.mat','_Synth.jpg');
    
    imwrite(bigSynthImg, fullfile(p.topDir, p.synthImgDir, savename), "jpg");

end

%% Save user input parameter list


p.time = datetime("now");
fields = fieldnames(p);
C=struct2cell(p);
params = [fields C];
%writecell(params,[p.topDir, 'particles/particleDetect_params.txt'],'Delimiter','tab')
%f.time = datetime; %add time to writeout

writecell(params,[p.topDir,'solved/','diskSolve_params.txt'],'Delimiter',',')


%% Output variable if everything is finished
out = 1;


