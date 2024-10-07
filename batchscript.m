%batch scripting

%%file parameters
fileParams.topDir = './testdata/'; % where the images are stored
fileParams.imgReg = '*.jpg'; %image format and regex
fileParams.imgDir = 'images/';
fileParams.warpedImgDir = 'warpedimg/';
fileParams.particleDir = 'particles/';
fileParams.contactDir = 'contacts/';
fileParams.solvedDir = 'solved/';
fileParams.synthImgDir = 'synthimg/';
fileParams.adjacencyDir = 'adjacency/';
fileParams.frameIdInd = 15;

%% verbose
verbose = true;
 
%%particledetect parameters
pdParams.cen = [120+5304/2, 113+5304/2]; %measure the center of annulus in images taken by the camera
pdParams.rad = [2840/2, 5304/2];
    % cen = [71+5313/2, 110+5313/2];
    %rad = [2783/2, 5313/2]; %measured in imageJ, pixels in untransformed image should be same as preprocess
pdParams.R_overlap = 15; %how much do the particle have to overlap by to remove double detected particles
pdParams.dtol = 25; %distance from edge
pdParams.radiusRange = [40,60];
pdParams.boundaryType = 'annulus';
pdParams.sigma = 50; % chosen by visual inspection

% particleTrack
ptParams.skipvalue = 200;

%contactdetect
cdParams.radiusRange = pdParams.radiusRange;
cdParams.metersperpixel = 0.015/939;
cdParams.fsigma = 141;%photoelastic stress coefficient
cdParams.g2cal = 145;%Calibration Value for the g^2 method, can be computed by joG2cal.m
%how far away can the outlines of 2 particles be to still be considered Neighbors
cdParams.dtol = 30;
cdParams.contactG2Threshold = 0.5;%sum of g2 in a contact area larger than this determines a valid contact
cdParams.CR = 15;
cdParams.boundaryType = pdParams.boundaryType;
cdParams.rednormal = 8;
cdParams.minpeakheight = 0.10;
cdParams.minpeakprominence = 0.06;
cdParams.minpeakprom_main = 0.06;
cdParams.padding = 1;
cdParams.roach = false;
cdParams.imadjust_limits = [0, 0.6];
cdParams.fineimadjust_limits = [0/255, 30/255];%[13/255, 39/255]
cdParams.sigma = pdParams.sigma; %for blurring large scale features
cdParams.polarizerstrip = [[2900,2930,3560,3593];[220,220,6099,6100]];
cdParams.calibrate =false; 
cdParams.figverbose = true;

%% disksolve 
dsParams.algorithm = 'levenberg-marquardt';
dsParams.maxIterations = 200;
dsParams.maxFunctionEvaluations = 400;
dsParams.functionTolerance = 0.01;
dsParams.scaling = 0.5;
dsParams.maskradius = 0.96;

%% newtonizer params
nwParams.boundaryType = pdParams.boundaryType;

%% adjacency matrix parameters
amParams.go = false; %compile adjaceceny lists for each image, set to false if you have already done this and only want to have the big list
amParams.fmin = 0.000001;%minimum force (in Newtons) to consider a contact a valid contact
amParams.fmax = 1000; %maximum force (in Newtons) to consider a contact a valid contact
amParams.emax = 2800; %maximum fit error/residual to consider a contact a valid contact
amParams.skipvalue = 2000; %I chose this as a result of my system size, could and should be altered based on your specific system and variability in finding particles


%%
moduleParams.pdParams = pdParams;
moduleParams.ptParams = ptParams;
moduleParams.cdParams = cdParams;
moduleParams.dsParams = dsParams;
moduleParams.nwParams = nwParams;
moduleParams.amParams = amParams;
PeGSModular(fileParams, moduleParams, verbose)
