%Written for PeGS 2.0 by Kerstin Nordstrom 8/24
%Adapted from Carmen Lee's adaptation of Jonathan Kollmer's PeGS 1.0
%CL's version copied below in entirety

%params: user must set radius range and boundary type as inputs. Other parameters to adjust
%may be defined in wrapper (p.sensitivity, p.dtol, p.edgethresh) otherwise will assume default
%values

%centers are saved as text files
%parameters are saved as text file
%verbose will display detected circles and save a sample image

% function out = particleDetect(f, p, verbose)
% %f = directories and image name pattern
% %p = parameters
% f.outDir = [f.topDir, '/particles/']; % where the output is saved
% if ~exist(f.outDir , 'dir')
%     mkdir(f.outDir);
% end



%sensitivity for edge finding is inside loop


% %%CL's ParticleDetect function%%
% 
% %function particle_detect(directory)
% % A script to find particle locations
function out = particleDetect(fileParams, pdParams, verbose)
%fileParams = directories and image name pattern
%pdParams = parameters for particle detect

%% FILE MANAGEMENT


    if ~exist(fullfile(fileParams.topDir, fileParams.particleDir) , 'dir')
        mkdir(fullfile(fileParams.topDir, fileParams.particleDir))
    end

    
    if verbose
        disp('starting particleDetect() to find all particle centroids and save results in particleDir')
    end

%% Set up for parameters is found in the function below
    pdParams = paramsSetUp(pdParams); %see below for defaults for params

if pdParams.boundaryType == "annulus"
    
    images=dir(fullfile(fileParams.topDir, fileParams.warpedImgDir,'*.tif'));
    if verbose
        disp([num2str(length(images)), ' images starting']);
    end
    nFrames = length(images);
    


    for frame = 1:nFrames
        
        im = imread(fullfile(images(frame).folder, images(frame).name));
        red = im(:,:,1);
        green = im(:,:,2);
        red = imsubtract(red, green*0.2); %this works for the annulus images, removes excess green
        red = imadjust(red, [0.20,0.60]); %this works for annulus, might need to tweak, brightens image

        sigma = pdParams.sigma; % chosen by visual inspection
        G = fspecial('gaussian', 3*sigma+1, sigma);
        yb = imfilter(red, G, 'replicate'); %removes large scale image features like bright spots
        red = bsxfun(@minus, red,yb);



        [centers,radii,metrics]=imfindcircles(red,pdParams.radiusRange,'objectpolarity','dark','sensitivity',0.945,'method','twostage','EdgeThreshold',0.02);%values found by tweaking
        % if you want to check out the images
        if verbose == true
            figure(1);
            imshow(red);
            viscircles(centers, radii);
            hold on;
            
        end
        %%
        xt = centers(:,1);
        yt = centers(:,2);
        rt = radii;


        %binarize radius
        rt(rt<49) = 44;
        rt(rt>49) = 55;
        %%      %beginning cleaning section

        %convert back to real space
        [midx,midy] = size(red);
        [theta,r] = cart2pol(xt-midx/2,yt-midy/2);

        d = -6.5*r.^2/(200*(925+6.5)); %6.5 is the thickness of the particles in mm, 925 is distance between particles and camera lens in mm
        
        s1 = d+r;
        [ut,vt] = pol2cart(theta,s1);
        ut = ut + midx/2;
        vt = vt + midy/2;

        ifcn = @(c) [ut(:) vt(:)];
        tform = geometricTransform2d(ifcn);
        [uv] = transformPointsInverse(tform, [0,0]); %particle original coordinates
        u = uv(:,1)-400;
        v = uv(:,2)-400;
        

        radialPos = sqrt((u-pdParams.cen(1)).^2+(v-pdParams.cen(2)).^2);
        closeind = find(radialPos <= pdParams.rad(1)+15 );
        closeind = sortrows(closeind, 'descend');


        xt(closeind) = [];
        yt(closeind) = [];
        
        rt(closeind) = [];
        metrics(closeind)= [];
        u(closeind) = [];
        v(closeind) = [];

        if verbose
            viscircles([xt, yt], rt,'EdgeColor', 'b');
        end
        %%
        %now we look for particles with a dramatic overlap
        dmat = pdist2([u,v],[u,v]); %Creates a distance matrix for particle center locations
        rmat = rt + rt'; %Makes a combination of radii for each particle

        friendmat = dmat < (rmat - 25) & dmat~=0; %Logical "friend" matrix
        [f1, f2] = find(friendmat == 1);


        badind = zeros(length(f1),1);

        M = length(f1);
        %this picks out the worse circle
        for n=1:M
            if metrics(f1(n)) > metrics(f2(n))
                badind(n) = f2(n);
            else
                badind(n) = f1(n);
            end
        end
        badind = badind(badind~=0);
        badind = unique(badind);
        badind = sortrows(badind, 'descend');


        xt(badind) = [];
        yt(badind) = [];
        
        rt(badind) = [];
        
        u(badind) = [];
        v(badind) = [];

        if verbose
            viscircles([xt, yt], rt,'EdgeColor','g'); %draw particle outline
            
        end
        %%
        dmat = pdist2([u,v],[u,v]);
        rmat = rt + rt';
        friendmat = dmat < (rmat -8 ) & dmat~=0; %Logical "friend" matrix

        % %friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
        [~, f2] = find(friendmat == 1);
        % [f3, f4] = find(friendmat2 == 1);
        badind = unique(f2);
        M = length(badind);
        toobig = zeros(M, 1);
        %
        badin2 = unique(f2);

        for n=1:M
            sum(f2 == badin2(n));
            if  sum(f2 == badin2(n))>2
                if rt(badind(n)) > 49 %cut off for middle of particle size
                    rt(badind(n)) = 44;
                    toobig(n) = badind(n);

                end
            end
        end
        




        %%

        radialPos = sqrt((u-pdParams.cen(1)).^2+(v-pdParams.cen(2)).^2);
        
        owi= radialPos <= pdParams.rad(2)+2.5*pdParams.dtol &radialPos >=pdParams.rad(2)-2.5*pdParams.dtol;
        iwi = find(radialPos <= pdParams.rad(1)+3.5*pdParams.dtol &radialPos >=pdParams.rad(1)-2.5*pdParams.dtol);
        edges = zeros(length(u), 1);
        edges(owi) = 1;
        edges(iwi) = -1;
        % if you want to check out the images
        if verbose == true
            viscircles([xt(owi), yt(owi)], rt(owi), 'Color', 'y');
            viscircles([xt(iwi), yt(iwi)], rt(iwi), 'Color', 'y');
            hold off
            drawnow;
        end
        particle = [xt, yt, rt, edges];
        if verbose
            disp(['image ', num2str(frame),' found ' num2str(length(particle)), ' particles'])
        end
        writematrix(particle,[fileParams.topDir,fileParams.particleDir, images(frame).name(1:end-4),'_centers.txt'])
    end


end




%% saving parameters in and finishing module

fields = fieldnames(pdParams);
for i = 1:length(fields)
    fileParams.(fields{i}) = pdParams.(fields{i});
end

fileParams.lastimagename=images(frame).name;
fileParams.time = datetime("now");
fields = fieldnames(fileParams);
C=struct2cell(fileParams);
pdParams = [fields C];
writecell(pdParams,fullfile(fileParams.topDir, fileParams.particleDir,'particleDetect_params.txt'),'Delimiter','tab')


if verbose 
        disp('done with particleDetect()');
end

out = true;

end

%% defaults for particleDetectModule
function p = paramsSetUp(p)
%classify boundary type
if isfield(p,'boundaryType') == 0
    p.boundaryType = "rectangle";
end

%set radius range
if isfield(p,'radiusRange') == 0
    p.radiusRange = [45 80];
end


%classify edge particles with tolerance
if isfield(p,'dtol') == 0
    p.dtol = 10;
end
end