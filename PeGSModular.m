% PeGSv2 uses a modular format
%      all modules to be called using the following 3 inputs and 1 output
%      all other returned values are outputs to disk
%      document the output formats used in each module

% standard file param structure is: (with example values)
%      FileParams.inDir = '/home/username1/DATA/';
%      FileParams.outDir = '/home/username2/analysis/';
%      FileParams.imgReg = '240726a.1*.jpg';

% standard module parameters structure looks like this: 
%      ModNameParams.diameter = 100; %pixels
%      ModNameParmas.filter = 7; 
% for all parameters the user will set to run the module
% use if ~exist() statements to provide default values inside function

% all modules should be functions called as follows:
%      ModuleName(MyFileParams, MyModParams, verbose)
%      and return(True) or return(False) depending on their success
%      setting verbose = True/False will turn on/off outputs

% by default, you can run the function from the directory where your
% images are stored, and it will write the output to that same
% directory in a subdirectory 'output' using all .jpg files in that directory

function PeGSModular(fileParams, moduleParams,verbose)
   

if exist('fileParams', 'var' ) == 0
    fileParams=struct;
end
if exist('verbose', 'var' ) == 0
   verbose = true;
end
fileParams = paramsSetUp(fileParams, verbose);


if ~exist('moduleParams', 'var')
    pdParams = struct;
    ptParams = struct;
    cdParams = struct;
    dsParams = struct;
    nwParams = struct;
    amParams = struct;
else
    if isfield(moduleParams, 'pdParams') ==1
        pdParams = moduleParams.pdParams;
    else
        pdParams = struct;
    end

    if isfield(moduleParams, 'ptParams') ==1
        ptParams = moduleParams.ptParams;
    else
        ptParams = struct;
    end

    if isfield(moduleParams, 'cdParams') ==1
        cdParams = moduleParams.cdParams;
    else
        cdParams = struct;
    end

    if isfield(moduleParams, 'dsParams') ==1
        dsParams = moduleParams.dsParams;
    else
        dsParams = struct;
    end
    
    if isfield(moduleParams, 'nwParams') ==1
        nwParams = moduleParams.nwParams;
    else
        nwParams = struct;
    end

    if isfield(moduleParams, 'amParams') ==1
        amParams = moduleParams.amParams;
    else
        amParams = struct;
    end

end

% these are basic steps to run PeGS on the sample images

%% module to detect contacts between particles. Set parameters in cdParams structure
%preprocess(fileParams, pdParams, verbose);
particleDetect(fileParams, pdParams, verbose);


%% particleTrack is optional. Comment out if you do not want the particles tracked from frame to frame
% module to track particles and assign them from frame to frame with the
% same id, set parameters with ptParams


%particleTrack(fileParams, ptParams, verbose);


%% module to detect contacts between particles. Set parameters in cdParams structure


%contactDetect(fileParams, cdParams, verbose);

%% module to solve the forces on the particles. Set parameters in dsParams structure


%diskSolve(fileParams, dsParams, verbose);

%% newtonize

%newtonize(fileParams, nwParams, verbose);
%% module create an adjacency matrix for all images in the data file. Set parameters in amParams structure

%
%adjacencyMatrix(fileParams, amParams, verbose);

return

end

%%%%%%%%%%%%%%%%%%%%%


function [fileParams] = paramsSetUp(fileParams, verbose)

if isfield(fileParams,'topDir') == 0
    fileParams.topDir = './testdata/'; % where the images are stored
end

if isfield(fileParams,'imgDir') == 0
    fileParams.imgDir = 'images'; % where the output is saved
    imgDirPath = fullfile(fileParams.topDir,fileParams.imgDir);
    if ~exist(imgDirPath , 'dir')
        mkdir(imgDirPath)
        error(['error: put your images in ', imgDirPath, ' or change your path'])
    end
end

if isfield(fileParams,'particleDir') == 0
    fileParams.particleDir = 'particles';
end

if isfield(fileParams,'contactDir') == 0
   fileParams.contactDir = 'contacts';
end

if isfield(fileParams,'solvedDir') == 0
   fileParams.solvedDir = 'solved';
end

if isfield(fileParams,'adjacencyDir') == 0
   fileParams.adjacencyDir = 'adjacency';
end

if isfield(fileParams,'imgReg') == 0
    fileParams.imgReg = '*.jpg'; %image format and regex
    % several Step*.jpg files are present on GitHub as sample data
end



if verbose
    disp(fileParams)
end



end