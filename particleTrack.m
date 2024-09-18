%function to track particles and to assign them to matching IDs, adapated
%from Abrar Nasser's code on Particle_IDtracing
function particleTrack(p, f, verbose)


directory = p.topDir;

particledirectory = [directory, p.particleDir;];
datafiles = dir([directory,p.particleDir,'*centers.txt']);
skipvalue = f.skipvalue; %for initialising the data structure

nFrames = length(datafiles);    %how many files
par_ref = load([particledirectory, datafiles(1).name]); %load the first set
frameId = str2double(datafiles(1).name(p.frameIdInd:p.frameIdInd+3));
skipamount = length(par_ref)+skipvalue; %I chose this as a result of my system size, could and should be altered based on your specific system and variability in finding particles

%centers has the format in 
%frame, particleID, x position, y position, r, edge 

centers = nan(nFrames*skipamount, 6); %this is the big array where we store the tracked ids for all of the frames

%initialize array and particle ids with first dataset
ids = (1:length(par_ref))';
par_ref = cat(2,ids, par_ref); 
centers(1:length(par_ref),:)=cat(2,ones(length(par_ref),1)*frameId, par_ref); 

%loop over datasets
for i = 1:numel(datafiles)-1

    par_curr = load([particledirectory,datafiles(i+1).name]);
    [tracked] = particletracingfn(par_ref,par_curr);
    
  
    %find the particles that were missed
    biggest = max(max(centers(:,2)));
    idx=find(tracked==0);
    
    for m = 1:length(idx)
        tracked(idx(m)) = biggest+m;
    end
    
    par_curr = cat(2, tracked, par_curr);
    
    if verbose
        figure(1);
        viscircles(par_ref(:,2:3), par_ref(:,4));
        
        viscircles(par_curr(:,2:3), par_curr(:,4), 'Color', 'b');
        for z = 1:size(par_curr, 1)
            text(par_curr(z, 2), par_curr(z, 3), num2str(tracked(z)), 'Color', 'white');
        end
        axis('equal')
        drawnow;
    end
    par_curr = sortrows(par_curr,1);
    par_ref = par_curr;
    frameId = str2double(datafiles(i+1).name(p.frameIdInd:p.frameIdInd+3));
    update = cat(2,ones(length(par_ref),1)*(frameId), par_ref);
    centers(skipamount*i:skipamount*i +length(par_ref)-1,:)=update;
end

centers(any(isnan(centers),2),:)=[];
writematrix(centers, [directory,'particle_positions.txt'])





if verbose
figure(1);
    refimages = dir([directory, p.imageDir, p.imgReg]);
    ref =imread([directory, p.imageDir, refimages(1).name]);
    imshow(ref)
    N = unique(centers(:,2));
    
    
    for particle = 1:length(N)
        ind = find(centers(:,2) ==particle);
        if length(ind)>1
            xx = [centers(ind, 3),centers(ind, 3)];
            yy = [centers(ind, 4),centers(ind, 4)];
            ff = [centers(ind,1),centers(ind,1)];
            zz=zeros(size(xx));
            surf(xx,yy,zz,ff,'EdgeColor','interp');
            %plot(centers(ind,3), centers(ind,4),'Color',cm(frame,:));
            hold on;
        end
    end
    axis('equal')
end


disp(['Particles tracked for ' num2str(nFrames),' images'])
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
%         store_esc --> sanity check
%                       should be zero if all particles have been detected
% SET STRINGENT AS 'TRUE' IF STORE_ESC HAS NON-ZERO VALUE

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
        tracked(find(detect == true)) = centers1(i,1);
        
    end
end

end


