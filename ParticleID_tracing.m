
%function particleTrack(p, f, verbose)


%directory = p.topDir;
directory = './';
particledirectory = [directory, 'particles/'];
datafiles = dir([directory,'particles/','*_centers.txt']);
stringent = false; %% set value as per info in the function at the end
skipvalue = 200;
verbose = false

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%No user input required beyond this line%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nFrames = length(datafiles);    %how many files

par_ref = load([particledirectory, datafiles(1).name]); %load the first set

skipamount = length(par_ref)+skipvalue; %I chose this as a result of my system size, could and should be altered based on your specific system and variability in finding particles
centers = nan(nFrames*skipamount, 6); %this is the big array where we store the tracked ids for all of the frames

store_save = [];

%initialize array and particle ids with first dataset
ids = (1:length(par_ref))';
par_ref = cat(2,ids, par_ref);

centers(1:length(par_ref),:)=cat(2,ones(length(par_ref),1), par_ref); 

%loop over particles
for i = 1:numel(datafiles)-1

    par_curr = load([particledirectory,datafiles(i+1).name]);
    
    
    num_par = length(par_ref); % number of particles from the previous frame
    
    [tracked, store_esc] = particletracingfn(par_ref,par_curr,stringent);
    store_save = [store_save; store_esc];
    
    n = length(tracked);
    [~,ind,~]=unique(tracked);
    out = unique(tracked(setdiff((1:n),ind)));
    biggest = max(max(centers(:,2)));
    idx=find(tracked==0);
    
    for m = 1:length(idx)
        tracked(idx(m)) = biggest+m;
    end
    
    par_curr = cat(2, tracked, par_curr);
    %particle = par_new;
    %save([directory,datafiles(i+1).name],'particle');
    if verbose
    figure(1);
    ref =imread([directory, 'warpedimg/', datafiles(i).name(1:end-12), '.tif']);
    imshow(ref);
    viscircles(par_ref(:,2:3), par_ref(:,4));
    for z = 1:num_par
        text(par_ref(z, 2), par_ref(z, 3), num2str(par_ref(z,1)), 'Color', 'white')
    end
    

    figure(2);
    curr =imread([directory, 'warpedimg/', datafiles(i+1).name(1:end-12), '.tif']);
    imshow(curr);
    viscircles(par_curr(:,2:3), par_curr(:,4));
    for z = 1:size(par_curr, 1);
        text(par_curr(z, 2), par_curr(z, 3), num2str(tracked(z)), 'Color', 'white');
    end
    drawnow;
    end
    par_curr = sortrows(par_curr,1);
    par_ref = par_curr;
    centers(skipamount*i:skipamount*i +length(par_ref)-1,:)=cat(2,ones(length(par_ref),1)*(i+1), par_ref);
end

centers(any(isnan(centers),2),:)=[];


figure(1);
    ref =imread([directory, 'warpedimg/', datafiles(i).name(1:end-12), '.tif']);
    
    imshow(ref)
    N = unique(centers(:,2));
    f = unique(centers(:,1));
    cm = colormap(parula(size(f,1))); 
    %for particle =1:100
    for particle = 1:length(N)
        ind = find(centers(:,2) ==particle);
        if length(ind)>1;
            xx = [centers(ind, 3),centers(ind, 3)];
            yy = [centers(ind, 4),centers(ind, 4)];
            ff = [centers(ind,1),centers(ind,1)];
            zz=zeros(size(xx));
            hs=surf(xx,yy,zz,ff,'EdgeColor','interp');
            %plot(centers(ind,3), centers(ind,4),'Color',cm(frame,:));
            hold on;
        end
    end

function [tracked, store_esc] = particletracingfn(centers1,centers2,stringent)
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
detect_store = false(size(centers2,1),1);

store_esc = 0;
for i = 1:num_par
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
        detect_store = detect_store | detect;
    end

    if sum(detect) ~=1
        store_esc = store_esc + 1;
    end


end
% if stringent
% store_esc = sum(tracked == 0);
% for i = 1:.25:5
%     num_esc = sum(tracked == 0);
%     if num_esc == 0
%         break
%     end
%     centers2_esc = centers2(detect_store == 0,:);
%     esc_ind = find(tracked == 0,num_esc,"first");
%     for j = 1:num_esc
%         shift = i * r1(esc_ind(j));
%         x1 = centers1(esc_ind(j),1)-shift;
%         x2 = centers1(esc_ind(j),1)+shift;
%         y1 = centers1(esc_ind(j),2)-shift;
%         y2 =  centers1(esc_ind(j),2)+shift;
%         xv = [x1 x1 x2 x2];
%         yv = [y1 y2 y2 y1];
%         detect = inpolygon(centers2_esc(:,1),centers2_esc(:,2),xv,yv);
% 
%         if sum(detect) == 1
%             store_esc = store_esc -1 ;
%             detect_store(centers2(:,1) == centers2_esc(detect,1) & centers2(:,2) == centers2_esc(detect,2)) = true;
%             tracked(esc_ind(j)) = find(centers2 == centers2_esc(detect));
%         end
% 
% 
% 
%     end
% end
%end
end


