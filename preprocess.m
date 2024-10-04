%function particle_detect(directory)
% A script to find particle locations
function preprocess(p, f, verbose)
% directory = './';
% imname = 'test.jpg';
% boundaryType = "annulus";
% verbose = false;
directory = [p.topDir, 'images/'];
images=dir([directory,p.imgReg]);
nFrames = length(images);
%nFrames = 1

%cen = [2710, 2768]; %measure the center of annulus in images taken by the camera
%rad = [2830/2, 5330/2];
%cen = [58+5313/2, 110+5313/2]; %before 240906
%rad = [2783/2, 5313/2]; %before 240906
%cen = [112+5304/2, 112+5304/2]; %measure this in imagej
%rad = [2790/2, 5304/2];

%line 58 contains experiment specific values)
xLimitsOut = [1,6304];
yLimitsOut = [1,6304]; %can be calculated just before imwarp is called but it is very slow. Recommend doing once then hardcoding in for speed)
%% for doing camera transforms for both camera lens distortion and particle distance correction on annulus


if f.boundaryType == "annulus"
    if not(isfolder(append(p.topDir,'warpedimg'))) %make a new folder with warped images
    mkdir(append(p.topDir,'warpedimg'))
    end
    im1=imread([directory,images(1).name]);
    
        %trim the image to get rid of edge and center points 
        ncols = (1:size(im1,2));
        nrows = (1:size(im1,1))';
        
        CC=sqrt((nrows-f.cen(2)).^2+(ncols-f.cen(1)).^2);
          
        Q=(CC>=f.rad(1)) & (CC<=f.rad(2));
        mask = Q;
        mask(:,:,2) = Q;
        %mask(:,:,3) = 0;
    for frame = 1:nFrames
    %frame = frame+130
   
    
    if ~isfile(append(p.topDir,'warpedimg/',images(frame).name(1:end-4),'warped.tif')) %check to see if images are already warped
        disp(frame)
        im1=imread([directory,images(frame).name]);
        im1 = im1(:,:,1:2);
        im1 = immultiply(im1, mask);
        
        %'trimmed'
        if verbose
            figure(1)
            im1(:,:,3) = 0 ;
            imshow(im1)
            drawnow;
            hold off
        end
        im1 = padarray(im1,[400 400], 'both'); % have to pad image with some amount of extra pixels to not lose any edge material when warping
        [nrows, ncols,~] = size(im1);

        [xi,yi] = meshgrid(1:ncols,1:nrows); %transform the points to polar coordinates
        xt = xi - ncols/2;
        yt = yi - nrows/2;
        [theta,r] = cart2pol(xt,yt);

        d = -6.5*r.^2/(200*(925+6.5)); %6.5 is the thickness of the particles in mm, 925 is distance between particles and camera lens in mm
        %rmax = max(r(:));
        s1 = d+r;
        [ut,vt] = pol2cart(theta,s1);
        ui = ut + ncols/2;
        vi = vt + nrows/2;
        %fill = [187, 22, 54];
        ifcn = @(c) [ui(:) vi(:)];
        tform = geometricTransform2d(ifcn);
        %'done transform' %finishes creating the geometric transform
    
        %[h,w,~] = size(im1);
        %[xLimitsOut,yLimitsOut] = outputLimits(tform,[1,w],[1,h]);
        %%uncomment to calculate for the first time
        
        Rout = imref2d(size(im1));
        Rout.XWorldLimits = xLimitsOut;      
        Rout.YWorldLimits = yLimitsOut; 
        %'limitsmade'
        im =imwarp(im1,tform, "OutputView",Rout);
        im(:,:,3) = 0;
    
        imwrite(im, append(p.topDir,'warpedimg/',images(frame).name(1:end-4),'warped.tif'))

        %'done image warp'
    
%     else
%     im = imread(append(directory,'warpedimg/',images(frame).name(1:end-4),'warped.tif'));

    end
    end
end
disp('files done preprocessing')
