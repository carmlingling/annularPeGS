%Updated to first release version of PeGS2 by Carmen Lee 29/9/2024
%function to track particles and to assign them to matching IDs, adapted by
% %Carmen Lee from Abrar Nasser's code on Particle_IDtracing
% 
% # ParticleID_Tracing
% Reorganising the random IDs post particle detection. (optional)
% 
% **I**
% Takes in the particle positions from particleDetect and organizes them based on their position (sets the same id)
% 
% **O**
% Exports a text file organized as
% 
% frame number | particle id | x | y | r | edge
% 
%  Allots unique ID to each particle in the ensemble. The alloted ID remains persistant and unique throughout the 
%  strain steps. The fields in kth ID in two different strain steps will correspond to the same particle but in two 
%  different strain steps.
%  
% Useful for Quasi-static problems where focus is on evolution of parameters linked to each particle respectively. 
% 
% Also exports a particleTrack_params.txt storing the tracking parameters into the particles directory (which is really just the padding value, called p.skipValue)


function out = particleTrack(fileParams, ptParams, verbose)

if ~isfield(ptParams, 'skipValue')
    ptParams.skipValue = 20; %for initializing the data structure, will depend on the amount of particles you are looking for and the flux in and out of frame
end

%% File Management

particledirectory = fullfile(fileParams.topDir, fileParams.particleDir);
datafiles = dir(fullfile(particledirectory,'*centers.txt')); %output from particleDetect
nFrames = size(datafiles,1);    %how many files

%% Initializing data array

%%loading data from first in series to initialize the datastructure
par_ref = load(fullfile(datafiles(1).folder, datafiles(1).name)); %load the first set

%if you want to pull id from the filename (i.e if named sequentially and
%not starting from 1)
if isfield(fileParams, 'frameIdInd')
    frameId = str2double(datafiles(1).name(fileParams.frameIdInd:fileParams.frameIdInd+3));
else
    frameId = 1;
end

skipamount = length(par_ref)+ptParams.skipValue; %I chose this as a result of my system size, could and should be altered based on your specific system and variability in finding particles


%centers has the format [frame, particleID, x position, y position, r, edge] 
centers = nan(nFrames*skipamount, 6); %this is the big array where we store the tracked ids for all of the frames

%initialize array and particle ids with first dataset
ids = (1:length(par_ref))';
par_ref = cat(2,ids, par_ref); 
centers(1:length(par_ref),:)=cat(2,ones(length(par_ref),1)*frameId, par_ref); 




%% track from the previous frame and populate data array
for i = 1:numel(datafiles)-1

    par_curr = load(fullfile(datafiles(i+1).folder,datafiles(i+1).name));
    [tracked] = particletracingfn(par_ref,par_curr);
    
  
    %find the particles that were missed
    biggest = max(max(centers(:,2)));
    idx=find(tracked==0);
    
    for m = 1:length(idx)
        tracked(idx(m)) = biggest+m;
    end
    
    par_curr = cat(2, tracked, par_curr);
    
    %draw motion per frame
%     if verbose
%         figure(1);
%         viscircles(par_ref(:,2:3), par_ref(:,4));
%         hold on;
%         viscircles(par_curr(:,2:3), par_curr(:,4), 'Color', 'b');
%         for z = 1:size(par_curr, 1)
%             text(par_curr(z, 2), par_curr(z, 3), num2str(tracked(z)), 'Color', 'white');
%         end
%         axis('equal')
%         drawnow;
%         hold off
%     end

    
    
    %if it is wanted to go by some regex in the filename
    if isfield(fileParams, 'frameIdInd')
        frameId = str2double(datafiles(i+1).name(fileParams.frameIdInd:fileParams.frameIdInd+3));
    else
        frameId = i+1;
    end
    
    if size(par_curr, 1) > skipamount
        error('skipvalue is not large enough, please increase')
    end %double check that the skip amount is large enough and not going to overwrite other data
    
    %reorganize
    par_curr = sortrows(par_curr,1);
    par_ref = par_curr;
    update = cat(2,ones(length(par_ref),1)*(frameId), par_ref);
    centers(skipamount*i:skipamount*i +length(par_ref)-1,:)=update; %add to big matrix
end

centers(any(isnan(centers),2),:)=[]; %remove spacers
writematrix(centers, fullfile(fileParams.topDir,'particle_positions.txt'))



%% final documentation

if verbose %visualize particle tracks
    figure(2);
    refimages = dir(fullfile(fileParams.topDir, fileParams.warpedImgDir, '*tif'));
    ref =imread(fullfile(refimages(1).folder, refimages(1).name));
    imshow(ref)
    hold on;
    
    scatter(centers(:,3), centers(:,4),30, centers(:,1), 'filled');
    colormap(parula(size(centers,1)))
    axis 'equal'
    savefig(fullfile(fileParams.topDir, fileParams.particleDir, 'trajectories')) 
end


%%save parameters 
fields = fieldnames(ptParams);
for i = 1:length(fields)
    fileParams.(fields{i}) = ptParams.(fields{i});
end


fileParams.time = datetime("now");
fields = fieldnames(fileParams);
C=struct2cell(fileParams);
ptParams = [fields C];
writecell(ptParams,fullfile(fileParams.topDir, fileParams.particleDir,'particleTrack_params.txt'),'Delimiter','tab')


if verbose 
        disp(['Particles tracked for ' num2str(nFrames),' images','done with contactTrack()']);
end

out = true;
end




function [tracked] = particletracingfn(centers1,centers2)
% input : refernce centers --> centers1
%          current centers  --> centers2
%          reference radius --> r1
%          total number of particles --> num_par
% output: tracked --> 1D matrix of length equal to num_par
%                     Its kth index stores a value p
%                     the fields in the pth ID of current struct
%                     is mapped to the fields in the kth ID in reference struct
%                     -----> kth ID particle in reference struct is 
%                            pth ID particle in current struct



num_par = size(centers1,1);
r1 = centers1(:,4);
tracked = zeros(size(centers2,1),1);



for i = 1:num_par %could probably vectorize this
    shift = r1(i)/1;
    x1 = centers1(i,2)-shift;
    x2 = centers1(i,2)+shift;
    y1 = centers1(i,3)-shift;
    y2 =  centers1(i,3)+shift;
    xv = [x1 x1 x2 x2];
    yv = [y1 y2 y2 y1];
    detect = inpolygon(centers2(:,1),centers2(:,2),xv,yv); %checks for centers within polygon with axes xv,yv 
    
    if sum(detect) == 1
        tracked(detect == true) = centers1(i,1);
        
    end
end

end



