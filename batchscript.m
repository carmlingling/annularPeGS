%batch scripting


FileParams.topDir = ['/eno/cllee3/matlab/annularPeGS', '/']; % where the images are stored
FileParams.imgReg = '*.jpg'; %image format and regex
% several Step*.jpg files are present on GitHub as sample data
FileParams.frameIdInd = 16;

 
particleDetectParams.cen = [112+5304/2, 112+5304/2]; %measure the center of annulus in images taken by the camera
particleDetectParams.rad = [2810/2, 5304/2];
    % cen = [71+5313/2, 110+5313/2];
    %rad = [2783/2, 5313/2]; %measured in imageJ, pixels in untransformed image should be same as preprocess
particleDetectParams.R_overlap = 15; %how much do the particle have to overlap by to remove double detected particles
particleDetectParams.dtol = 25; %distance from edge
particleDetectParams.radiusRange = [40,60];
particleDetectParams.boundaryType = 'annulus';
particleDetectParams.sigma = 50; % chosen by visual inspection

particleTrackParams.skipvalue = 200

cdParams.radiusRange = particleDetectParams.radiusRange;
cdParams.metersperpixel = 0.015/939;
cdParams.fsigma = 141;%photoelastic stress coefficient
cdParams.g2cal = 145;%Calibration Value for the g^2 method, can be computed by joG2cal.m
%how far away can the outlines of 2 particles be to still be considered Neighbors
cdParams.dtol = 30;
cdParams.contactG2Threshold = 0.5;%sum of g2 in a contact area larger than this determines a valid contact
cdParams.CR = 15;
cdParams.boundaryType = particleDetectParams.boundaryType;
cdParams.rednormal = 8;

cdParams.minpeakheight = 0.10;
cdParams.minpeakprominence = 0.06;
cdParams.minpeakprom_main = 0.06;

cdParams.imadjust_limits = [0, 0.6];
cdParams.fineimadjust_limits = [0/255, 30/255];%[13/255, 39/255]

padding = 1;
cdParams.sigma = particleDetectParams.sigma; %for blurring large scale features
cdParams.polarizerstrip = [[2731,2719,3643,3666];[212,212,6099,6100]];
cdParams.calibrate =false; 
cdParams.figverbose = true;

dsParams.algorithm = 'levenberg-marquardt';
% Function evalution limits

dsParams.maxIterations = 200;
dsParams.maxFunctionEvaluations = 400;
dsParams.functionTolerance = 0.01;

%% Scaling of image for fit
dsParams.scaling = 0.5;

%% Masking radius
% How much of particle edge to remove before fit, 1 is the entire particle
dsParams.maskradius = 0.96;

verbose = true

PeGSModular(FileParams, particleDetectParams,cdParams, dsParams, verbose)
