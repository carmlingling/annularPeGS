% update Lori S McCabe 08/2024 for PEGS 2.0 wrapper compatibility
% original version Carmen Lee modular and converted as basis for PEGS 2.0

function out = contactDetection(p, f, verbose)
%p passed in from wrapper
%requires input/output directories, image names

%contactDectectParams (f) are parameters specific for contact detection
%this includes all thresholding limits (listed below)


% directory = './';
% fileNames = 'Step09-0001.jpg';
% radiusRange = [40, 90];
% verbose = true;



%% thresholding limits (setting defaults if none listed in wrapper)
if ~isfield(f,'radiusRange') %range of radii for particles (in pixels)
    f.radiusRange = [40 90];
end

if ~isfield(f,'metersperpixel') 
    f.metersperpixel = .007/160;
end

if ~isfield(f,'fsigma') %photoelastic stress coefficient
    f.fsigma = 140;
end

if ~isfield(f,'g2cal') %calibration Value for the g^2 method, can be computed by joG2cal.m (PEGS 1.0 version)
    f.g2cal = 100;
end

if ~isfield(f,'dtol') %how far away can the outlines of 2 particles be to still be considered Neighbors
    f.dtol = 10;
end

if ~isfield(f,'contactG2Threshold') %sum of g2 in a contact area larger than this determines a valid contact
    f.contactG2Threshold = 0.5; 
end

if ~isfield(f,'CR') %contact radius over which contact gradient is calculated
    f.CR = 10;
end

if ~isfield(f,'imadjust_limits') %adjust contrast in green channel
    f.imadjust_limits = [0,1];
end

if ~isfield(f,'rednormal') %fractional amount to subtract the red channel from the green channel (Rimg/rednormal)
    f.rednormal = 2;
end


if ~isfield(f,'figverbose') %show all the figures as you go
    f.figverbose = false;
end


%% directory buisness and importing files
%% file housekeeping
if not(isfolder(append(p.topDir,'particles'))) %make a new folder with particle centers
        mkdir(append(p.topDir,'particles'));
end





disp([p.topDir, 'warpedimg/', p.imgReg]) %print the filename
files = dir([ p.topDir, 'warpedimg/', p.imgReg]); %images

centersfile = dir([p.topDir, 'particles/','*_centers.txt']);
files

%% setting up mask
mask = abs(-f.CR:f.CR);
mask = mask.^2 + mask.^2';
maskCR = double(sqrt(mask) <= f.CR-1);



%% image manipulation
for imgnumb = 1:size(files,1)

    clear particle %reinitialize the particle structure

    %read in image
    Img = imread([p.topDir,'warpedimg/',files(imgnumb).name]);
    Rimg = Img(:,:,1);
    Gimg = Img(:,:,2); %force image


    %adjust green image contrast by subtracting red channel and adjusting
    %contrast levels as set by imadjust_limits. This will need to be
    %tweaked for different lighting and transmission or reflection
    %photoelasticimetry
    
    if f.boundaryType == "annulus"
        
        bckgnd = poly2mask(f.polarizerstrip(1,:),f.polarizerstrip(2,:), length(Gimg), length(Gimg));
        Gimg = inpaintCoherent(Gimg,bckgnd,'SmoothingFactor',5,'Radius',15);
        
        if f.calibrate == true
           figure;
           imshow(Gimg);
           hold on
           plot([2731,2719,3643,3666, 2731],[212,212,6099,6100, 212],'b','LineWidth',1)
           drawnow;
        end
        
        Gimg=imsubtract(Gimg,Rimg./f.rednormal);
   
        G = fspecial('gaussian', 3*f.sigma+1, f.sigma);
        yb = imfilter(imcomplement(Rimg), G, 'replicate');
        Gimg = bsxfun(@minus, Gimg,yb*.09);

        Gimg= im2double(Gimg);
        Gimg = Gimg.*(Gimg > 0);
    
    
        Gimgd = imadjust(Gimg,f.imadjust_limits); %regular contrast
    
        Gimgfine = imadjust(Gimg, f.fineimadjust_limits); %super boosted contrast
        frame = str2double(files(imgnumb).name(p.frameIdInd:p.frameIdInd+3));
        disp(files(imgnumb).name(p.frameIdInd:p.frameIdInd+3))
        disp(frame)

    end
   
%% initialize data structure
   
    if f.figverbose
        figure(1); %makes a lot of figures, override onto 1
        
        imshow(Gimg)
        title('Gimg')

    end
    
    %frame = str2double(files(imgnumb).name(frameidind:frameidind+3))


    %% initialize data structure
    pData = readmatrix([p.topDir, 'particles/',centersfile(imgnumb).name]); %,"NumHeaderLines", 1); %Read Position data from centers file
    
   
    if ~isempty(pData)

        N = size(pData,1);

        particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[],  'edge', 0);
        for n = 1:N %Bookkeeping from centers-tracked
            particle(n).id= n;
            particle(n).x = pData(n,1);
            particle(n).y = pData(n,2);
            particle(n).r = round(pData(n,3));
            particle(n).edge = pData(n, 4);
            particle(n).rm = particle(n).r*f.metersperpixel;
            particle(n).fsigma = f.fsigma;
        end

        if f.figverbose
            imshow(Gimg);
            viscircles([pData(:,1),pData(:,2)],pData(:,3));
            hold on;
        end

        for n=1:N %loop over particles

            %create a circular mask

            r = particle(n).r;
            if round(particle(n).y+r)<size(Gimg, 1)&&round(particle(n).x+r)<size(Gimg,2)&&round(particle(n).y-r)>1&&round(particle(n).x-r)>1 %double check to make sure the bounds are within the image

                mask = abs(-r:r);
                mask = mask.^2 + mask.^2';
                mask1 = (sqrt(mask) <= r);

                %This crops out a particle
                cropXstart = round(particle(n).x-r);
                cropXstop = round(particle(n).x-r)+ size(mask1,1)-1;
                cropYstart = round(particle(n).y-r);
                cropYstop = round(particle(n).y-r)+ size(mask1,2)-1;


                particleImg= Gimg(cropYstart:cropYstop, cropXstart:cropXstop).*mask1;
                particle(n).forceImage=particleImg; %save this so we can fit to this image later in diskSolve


                %create a circular mask with a radius that is one pixel smaller
                %for cropping out the relevant gradient

                mask2 = double(sqrt(mask) <= r-1);

                %Compute G^2 for each particle
                [gx,gy] = gradient(particleImg);
                g2 = (gx.^2 + gy.^2).*mask2;
                particle(n).g2 = sum(sum(g2));
                particle(n).f = particle(n).g2/f.g2cal; %saving some particle scale features
            else
                error('badimage!!')

            end
        end
        %% look at neighbours

        xmat = pData(:,1);
        ymat = pData(:,2);
        rmat = pData(:,3);

        rmats = rmat; %Saves our radius matrix for later

        dmat = pdist2([xmat,ymat],[xmat,ymat]); %Creates a distance matrix for particle center locations
        rmat = rmat + rmat'; %Makes a combination of radii for each particle

        friendmat = dmat < (rmat + f.dtol) & dmat~=0; %Logical "friend" matrix

        friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
        [f1, f2] = find(friendmat == 1); %Creates an index of particles that are considered touching
        %%
        xpairs = [xmat(f1),xmat(f2)]; %set up an array of pairs of x, y and r for easy manipulation later
        ypairs = [ymat(f1),ymat(f2)];
        rpairs = [rmats(f1),rmats(f2)];

        %% loop over friends

        for l = 1:length(f1)

            x = xpairs(l,:);
            y = ypairs(l,:);
            r = rpairs(l,:);

                        
            if f.figverbose
            plot(x, y, 'LineWidth', 2)
            title('neighbour candidates')
            end
            
            [contactG2p, contactIp] = contactspot(x,y,r, f.CR, Gimg, maskCR);

            if(contactG2p(1) > f.contactG2Threshold && contactG2p(2) > f.contactG2Threshold)

                %this is a valid contact, remember it
                particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
                particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p(1); %remember the g2 value of the current contact area
                particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp(1); %changes to color
                particle(f1(l)).color(particle(f1(l)).z)='r'; %changes to color
                particle(f1(l)).neighbours(particle(f1(l)).z) = particle(f2(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
                particle(f1(l)).betas(particle(f1(l)).z) = atan2(y(2)-y(1),x(2)-x(1)); %the contact angle to particle m is now noted in the particle l datastructure
                particle(f2(l)).z= particle(f2(l)).z+1; %increase coordination number
                particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p(2); %remember the g2 value of the current contact area
                particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp(2);
                particle(f2(l)).color(particle(f2(l)).z)='r'; %changes to color
                particle(f2(l)).neighbours(particle(f2(l)).z) = particle(f1(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
                particle(f2(l)).betas(particle(f2(l)).z) = atan2(y(1)-y(2),x(1)-x(2));


            end




        end
        %%

        %Check if any of the walls is a neighbour as well

        circs = [[particle.y]', [particle.x]', [particle.r]', [particle.edge]']; %Makes a circs matrix from old matrices

        
        for ep = 1:length(particle)
            if circs(ep,4) == 1
                contacts = 0;
            elseif circs(ep,4) == -1
                contacts = pi;
            elseif circs(ep,4) == 2
                contacts = pi/2;
            elseif circs(ep,4) == -2
                contacts = -pi/2;
            end
            if particle(ep).edge ~=0
                x = particle(ep).x;
                y = particle(ep).y;
                r = particle(ep).r;
                for c =1:length(contacts) %technically doesn't need to be a loop but we will keep it for alternate scenarios
                    [contactG2p, contactIp]= contactspotwall(x, y, r, f.CR, contacts(c),Gimg, maskCR);
                    if(contactG2p > f.contactG2Threshold)
                        particle(ep).z= particle(ep).z +1; %increase coordination number
                        particle(ep).contactG2s(particle(ep).z)=contactG2p;
                        particle(ep).contactIs(particle(ep).z)=contactIp;
                        particle(ep).neighbours(particle(ep).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
                        particle(ep).betas(particle(ep).z) = contacts(c); %the contact angle to the wall is now noted in the particle l datastructure
                        particle(ep).color(particle(ep).z)='g';
                        %     else
                    end
                end
            end
        end


    end




if f.figverbose 
    h3 = figure(20);
    hAx1 = subplot(1,1,1,'Parent', h3);
    imshow(Gimg, 'Parent', hAx1);
    hold (hAx1, 'on');
    for n = 1:length(particle)
        particle(n).id;
        %viscircles([particle(n).x, particle(n).y], particle(n).r, 'EdgeColor', particle(n).color);
        z = particle(n).z; %get particle coordination number
        if (z>0) %if the particle does have contacts
            for m = 1:z %for each contact
                %draw contact lines
                lineX(1)=particle(n).x;
                lineY(1)=particle(n).y;
                lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
                lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));
                viscircles([lineX(1), lineY(1)], particle(n).r, 'color', 'blue');
                viscircles([lineX(1) + (particle(n).r-f.CR) * cos(particle(n).betas(m)) lineY(1) + (particle(n).r-f.CR) * sin(particle(n).betas(m))], f.CR, 'color', 'white');

                plot(hAx1, lineX, lineY,particle(n).color(m),'LineWidth',2);
            end
        end
        text(hAx1, particle(n).x, particle(n).y, num2str(particle(n).id), 'Color', 'y')
    end
    drawnow;

    saveas(h3,[p.topDir,'Contacts_',files(imgnumb).name(1:end-4)],'jpg') %save fig(20) 
end %PEGS wrapper verbose


disp([num2str(sum([particle.z])), ' contacts detected'])

%save updated particle contact info
save([p.topDir,'particles/', files(imgnumb).name(1:end-4),'_preprocessing.mat'],'particle')

end %loop imgnumb

%save parameters 
fields = fieldnames(f);
for i = 1:length(fields)
    p.(fields{i}) = f.(fields{i})
end
%p = rmfield(p,'radiusRange');
p.lastimagename=files(imgnumb).name;
p.time = datetime("now");
fields = fieldnames(p);
C=struct2cell(p);
params = [fields C];
%writecell(params,[p.topDir, 'particles/particleDetect_params.txt'],'Delimiter','tab')
%f.time = datetime; %add time to writeout

writecell(params,[p.topDir,'particles/','contactDetect_params.txt'],'Delimiter',',')

% writetable(f,[p.outDir, files(imgnumb).name(1:end-4),'_parameters.txt'])


%% wrap up contact Detect function
out = true; %if function ran completely, return true to PEGS wrapper


%% other called functions used in the contact detect algorithm 
function contactG2 = gradientcalculator(imgchunk)
[gx,gy] = gradient(imgchunk);
g2 = (gx.^2 + gy.^2);
contactG2 = sum(sum(g2));

end

function [contactG2p, contactIp]=contactspot(x, y, r, CR, Gimgd, maskCR)
contactangle = [atan2(y(2)-y(1),x(2)-x(1)), atan2(y(1)-y(2), x(1)-x(2))];
contactXp = round(x + (r -  1 - CR).* cos(contactangle));
contactYp = round(y + (r - 1 - CR).* sin(contactangle));

contactImg = im2double(imcrop(Gimgd,[contactXp(1)-CR contactYp(1)-CR CR*2 CR*2]));
contactImg = contactImg.*maskCR;

contactG2p = [gradientcalculator(contactImg)];
contactIp = [sum(sum(contactImg))];

contactImg = im2double(imcrop(Gimgd,[contactXp(2)-CR contactYp(2)-CR CR*2 CR*2]));
contactImg = contactImg.*maskCR;
contactG2p(2,:)= gradientcalculator(contactImg);
contactIp(2,:) = sum(sum(contactImg));
%contactG2p = [G1 G2]
end

function [contactG2p, contactIp]=contactspotwall(x, y, r, CR, angle,Gimgd, maskCR)

contactX = round(x + (r -  1 - CR).* cos(angle));
contactY = round(y + (r -1- CR).* sin(angle));


contactImg = im2double(imcrop(Gimgd,[contactX-CR contactY-CR CR*2 CR*2]));
contactImg = contactImg.*maskCR;

contactG2p = gradientcalculator(contactImg);

contactIp = sum(sum(contactImg));

end

%% end of other functions block


end %function contact dectect


%% CL original code below 
% 
% 
% %function contactDetection(directory, fileNames,verbose)
% 

% 
% 
% 
% %% image manipulation
% for imgnumb = 1:length(files)
% 
%     clear particle %reinitialize the particle structure
% 
%     %read in image
%     Img = imread([directory,files(imgnumb).name]);
%     Rimg = Img(:,:,1);
%     Gimg = Img(:,:,2); %force image
% 
% 
%     %adjust green image contrast by subtracting red channel and adjusting
%     %contrast levels as set by imadjust_limits. This will need to be
%     %tweaked for different lighting and transmission or reflection
%     %photoelasticimetry
% 
%     Gimg = Gimg-Rimg./rednormal;
%     Gimg= im2double(Gimg);
% 
%     Gimg = Gimg.*(Gimg > 0);
%     Gimg = imadjust(Gimg,imadjust_limits);
% 
%     if verbose
%         figure;
% 
%         imshow(Gimg)
%         title('Gimg')
% 
%     end
% 
%     %frame = str2double(files(imgnumb).name(frameidind:frameidind+3))
% 
% 
%     %% initialize data structure
%     pData = readmatrix([datadirectory,centersfile.name]); %,"NumHeaderLines", 1); %Read Position data from centers file
% 
% 
%     if ~isempty(pData)
% 
%         N = size(pData,1);
% 
%         particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[],  'edge', 0);
%         for n = 1:N %Bookkeeping from centers-tracked
%             particle(n).id= n;
%             particle(n).x = pData(n,1);
%             particle(n).y = pData(n,2);
%             particle(n).r = round(pData(n,3));
%             particle(n).edge = pData(n, 4);
%             particle(n).rm = particle(n).r*metersperpixel;
%             particle(n).fsigma = fsigma;
%         end
% 
%         if verbose
%             imshow(Gimg);
%             viscircles([pData(:,1),pData(:,2)],pData(:,3));
%             hold on;
%         end
% 
%         for n=1:N %loop over particles
% 
%             %create a circular mask
% 
%             r = particle(n).r;
%             if round(particle(n).y+r)<size(Gimg, 1)&&round(particle(n).x+r)<size(Gimg,2)&&round(particle(n).y-r)>1&&round(particle(n).x-r)>1 %double check to make sure the bounds are within the image
% 
%                 mask = abs(-r:r);
%                 mask = mask.^2 + mask.^2';
%                 mask1 = (sqrt(mask) <= r);
% 
%                 %This crops out a particle
%                 cropXstart = round(particle(n).x-r);
%                 cropXstop = round(particle(n).x-r)+ size(mask1,1)-1;
%                 cropYstart = round(particle(n).y-r);
%                 cropYstop = round(particle(n).y-r)+ size(mask1,2)-1;
% 
% 
%                 particleImg= Gimg(cropYstart:cropYstop, cropXstart:cropXstop).*mask1;
%                 particle(n).forceImage=particleImg; %save this so we can fit to this image later in diskSolve
% 
% 
%                 %create a circular mask with a radius that is one pixel smaller
%                 %for cropping out the relevant gradient
% 
%                 mask2 = double(sqrt(mask) <= r-1);
% 
%                 %Compute G^2 for each particle
%                 [gx,gy] = gradient(particleImg);
%                 g2 = (gx.^2 + gy.^2).*mask2;
%                 particle(n).g2 = sum(sum(g2));
%                 particle(n).f = particle(n).g2/g2cal; %saving some particle scale features
%             else
%                 error('badimage!!')
% 
%             end
%         end
%         %% look at neighbours
% 
%         xmat = pData(:,1);
%         ymat = pData(:,2);
%         rmat = pData(:,3);
% 
%         rmats = rmat; %Saves our radius matrix for later
% 
%         dmat = pdist2([xmat,ymat],[xmat,ymat]); %Creates a distance matrix for particle center locations
%         rmat = rmat + rmat'; %Makes a combination of radii for each particle
% 
%         friendmat = dmat < (rmat + dtol) & dmat~=0; %Logical "friend" matrix
% 
%         friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
%         [f1, f2] = find(friendmat == 1); %Creates an index of particles that are considered touching
%         %%
%         xpairs = [xmat(f1),xmat(f2)]; %set up an array of pairs of x, y and r for easy manipulation later
%         ypairs = [ymat(f1),ymat(f2)];
%         rpairs = [rmats(f1),rmats(f2)];
% 
%         %% loop over friends
% 
%         for l = 1:length(f1)
% 
%             x = xpairs(l,:);
%             y = ypairs(l,:);
%             r = rpairs(l,:);
% 
% 
%             if verbose
%             plot(x, y, 'LineWidth', 2)
%             title('neighbour candidates')
%             end
% 
%             [contactG2p, contactIp] = contactspot(x,y,r, CR, Gimg, maskCR);
% 
%             if(contactG2p(1) > contactG2Threshold && contactG2p(2) > contactG2Threshold)
% 
%                 %this is a valid contact, remember it
%                 particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
%                 particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p(1); %remember the g2 value of the current contact area
%                 particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp(1); %changes to color
%                 particle(f1(l)).color(particle(f1(l)).z)='r'; %changes to color
%                 particle(f1(l)).neighbours(particle(f1(l)).z) = particle(f2(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
%                 particle(f1(l)).betas(particle(f1(l)).z) = atan2(y(2)-y(1),x(2)-x(1)); %the contact angle to particle m is now noted in the particle l datastructure
%                 particle(f2(l)).z= particle(f2(l)).z+1; %increase coordination number
%                 particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p(2); %remember the g2 value of the current contact area
%                 particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp(2);
%                 particle(f2(l)).color(particle(f2(l)).z)='r'; %changes to color
%                 particle(f2(l)).neighbours(particle(f2(l)).z) = particle(f1(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
%                 particle(f2(l)).betas(particle(f2(l)).z) = atan2(y(1)-y(2),x(1)-x(2));
% 
% 
%             end
% 
% 
% 
% 
%         end
%         %%
% 
%         %Check if any of the walls is a neighbour as well
% 
%         circs = [[particle.y]', [particle.x]', [particle.r]', [particle.edge]']; %Makes a circs matrix from old matrices
% 
% 
%         for p = 1:length(particle)
%             if circs(p,4) == 1
%                 contacts = 0;
%             elseif circs(p,4) == -1
%                 contacts = pi;
%             elseif circs(p,4) == 2
%                 contacts = pi/2;
%             elseif circs(p,4) == -2
%                 contacts = -pi/2;
%             end
%             if particle(p).edge ~=0
%                 x = particle(p).x
%                 y = particle(p).y
%                 r = particle(p).r
%                 for c =1:length(contacts) %technically doesn't need to be a loop but we will keep it for alternate scenarios
%                     [contactG2p, contactIp]= contactspotwall(x, y, r, CR, contacts(c),Gimg, maskCR)
%                     if(contactG2p > contactG2Threshold)
%                         particle(p).z= particle(p).z +1; %increase coordination number
%                         particle(p).contactG2s(particle(p).z)=contactG2p;
%                         particle(p).contactIs(particle(p).z)=contactIp;
%                         particle(p).neighbours(particle(p).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
%                         particle(p).betas(particle(p).z) = contacts(c); %the contact angle to the wall is now noted in the particle l datastructure
%                         particle(p).color(particle(p).z)='g';
%                         %     else
%                     end
%                 end
%             end
%         end
% 
% 
%     end
% 
% 
% 
% 
% if verbose
%     h3 = figure(20);
%     hAx1 = subplot(1,1,1,'Parent', h3);
%     imshow(Gimg, 'Parent', hAx1);
%     hold (hAx1, 'on');
%     for n = 1:length(particle)
%         particle(n).id;
%         %viscircles([particle(n).x, particle(n).y], particle(n).r, 'EdgeColor', particle(n).color);
%         z = particle(n).z; %get particle coordination number
%         if (z>0) %if the particle does have contacts
%             for m = 1:z %for each contact
%                 %draw contact lines
%                 lineX(1)=particle(n).x;
%                 lineY(1)=particle(n).y;
%                 lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
%                 lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));
%                 viscircles([lineX(1), lineY(1)], particle(n).r, 'color', 'blue')
%                 viscircles([lineX(1) + (particle(n).r-CR) * cos(particle(n).betas(m)) lineY(1) + (particle(n).r-CR) * sin(particle(n).betas(m))], CR, 'color', 'white')
% 
%                 plot(hAx1, lineX, lineY,particle(n).color(m),'LineWidth',2);
%             end
%         end
%         text(hAx1, particle(n).x, particle(n).y, num2str(particle(n).id), 'Color', 'y')
%     end
%     drawnow;
% end
% disp([num2str(sum([particle.z])), 'detected'])
% save([directory, 'particles/', files(imgnumb).name(1:end-4),'_preprocessing.mat'],'particle')
% end
% 
% 
% 
% function contactG2 = gradientcalculator(imgchunk)
% [gx,gy] = gradient(imgchunk);
% g2 = (gx.^2 + gy.^2);
% contactG2 = sum(sum(g2));
% 
% end
% 
% 
% function [contactG2p, contactIp]=contactspot(x, y, r, CR, Gimgd, maskCR)
% contactangle = [atan2(y(2)-y(1),x(2)-x(1)), atan2(y(1)-y(2), x(1)-x(2))];
% contactXp = round(x + (r -  1 - CR).* cos(contactangle));
% contactYp = round(y + (r - 1 - CR).* sin(contactangle));
% 
% contactImg = im2double(imcrop(Gimgd,[contactXp(1)-CR contactYp(1)-CR CR*2 CR*2]));
% contactImg = contactImg.*maskCR;
% 
% contactG2p = [gradientcalculator(contactImg)];
% contactIp = [sum(sum(contactImg))];
% 
% contactImg = im2double(imcrop(Gimgd,[contactXp(2)-CR contactYp(2)-CR CR*2 CR*2]));
% contactImg = contactImg.*maskCR;
% contactG2p(2,:)= gradientcalculator(contactImg);
% contactIp(2,:) = sum(sum(contactImg));
% %contactG2p = [G1 G2]
% end
% 
% function [contactG2p, contactIp]=contactspotwall(x, y, r, CR, angle,Gimgd, maskCR)
% 
% contactX = round(x + (r -  1 - CR).* cos(angle));
% contactY = round(y + (r -1- CR).* sin(angle));
% 
% 
% contactImg = im2double(imcrop(Gimgd,[contactX-CR contactY-CR CR*2 CR*2]));
% contactImg = contactImg.*maskCR;
% 
% contactG2p = gradientcalculator(contactImg);
% 
% contactIp = sum(sum(contactImg));
% 
% end
