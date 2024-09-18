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

function PeGSModular(FileParams, particleDetectParams,particleTrackParams, cdParams, dsParams, verbose)
if ~exist('FileParams')
    FileParams=struct;
end


if isfield(FileParams,'imgReg') == 0
    FileParams.imgReg = '*.jpg'; %image format and regex
    % several Step*.jpg files are present on GitHub as sample data
end

if exist('verbose') == 0
    verbose = true;
end

if verbose
    disp(FileParams)
end

% these are basic steps to run PeGS on the sample images

%% set parameters for particleDetect() before running

%preprocess(FileParams, particleDetectParams, verbose)

%particleDetect(FileParams, particleDetectParams, verbose);

%particleTrack(FileParams, particleTrackParams, verbose);
%% set parameters for contactDetection() before running
%cdParams = struct;
contactDetection(FileParams, cdParams, verbose);
% if verbose 
%         disp('done with contactDetection()');
% end
% 
% userOptions = struct;
% diskSolve(FileParams, userOptions, verbose);
% if verbose 
%         disp('done with diskSolve()');
% end
%%
%newtonize(topDirectory, imageNames, boundaryType, verbose)
%'newtonized and edges handled'
%%
%adjacencyMatrix(topDirectory, imageNames, boundaryType, frameidind,verbose)
%'Adjacency matrix built'

%%
%particleTrack(topDirectory, imageNames, boundaryType, frameidind, verbose, false);
%'trajectories connected'
return
