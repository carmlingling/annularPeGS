% update Lori S McCabe 08/2024 for PEGS 2.0 wrapper compatibility
% original version Carmen Lee modular and converted as basis for PEGS 2.0

function out = contactDetect(fileParams, cdParams, verbose)
%p passed in from wrapper
%requires input/output directories, image names

%contactDectectParams (f) are parameters specific for contact detection
%this includes all thresholding limits (listed below)


% directory = './';
% fileNames = 'Step09-0001.jpg';
% radiusRange = [40, 90];
% verbose = true;



%% thresholding limits (setting defaults if none listed in wrapper)
cdParams = setupParams(cdParams);

%% directory buisness and importing files
warning('off','signal:findpeaks:largeMinPeakHeight')
global particleNumber1 particleNumber2
particleNumber1 = 222;
particleNumber2 = 1060;
cdParams.roach = false;
%% file housekeeping
if not(isfolder(append(fileParams.topDir,fileParams.contactDir))) %make a new folder with particle centers
        mkdir(append(fileParams.topDir,fileParams.contactDir));
end



%disp([p.topDir, 'warpedimg/', p.imgReg]) %print the filename
files = dir([ fileParams.topDir,fileParams.warpedImgDir, '*.tif']); %images

centersfile = dir([fileParams.topDir, 'particle_positions.txt']);
centersdata = readmatrix(fullfile(centersfile.folder, centersfile.name));

%% setting up mask
mask = abs(-cdParams.CR:cdParams.CR);
mask = mask.^2 + mask.^2';
maskCR = double(sqrt(mask) <= cdParams.CR-1);



%% image manipulation
for imgnumb = 1:size(files,1)

    clear particle %reinitialize the particle structure

    %read in image
    Img = imread([fileParams.topDir,fileParams.warpedImgDir,files(imgnumb).name]);
    Rimg = Img(:,:,1);
    Gimg = Img(:,:,2); %force image


    %adjust green image contrast by subtracting red channel and adjusting
    %contrast levels as set by imadjust_limits. This will need to be
    %tweaked for different lighting and transmission or reflection
    %photoelasticimetry
    
    if cdParams.boundaryType == "annulus"
        
        bckgnd = poly2mask(cdParams.polarizerstrip(1,:),cdParams.polarizerstrip(2,:), length(Gimg), length(Gimg));
        Gimg = inpaintCoherent(Gimg,bckgnd,'SmoothingFactor',5,'Radius',15);
        
        if cdParams.calibrate == true
           figure;
           imshow(Gimg);
           axis on
           hold on
           plot([cdParams.polarizerstrip(1,:) cdParams.polarizerstrip(1:1)],[cdParams.polarizerstrip(2,:) cdParams.polarizerstrip(2,1)],'b','LineWidth',1)
           drawnow;
        end
        
        Gimg=imsubtract(Gimg,Rimg./cdParams.rednormal);
   
        G = fspecial('gaussian', 3*cdParams.sigma+1, cdParams.sigma);
        yb = imfilter(imcomplement(Rimg), G, 'replicate');
        Gimg = bsxfun(@minus, Gimg,yb*.09);

        Gimg= im2double(Gimg);
        Gimg = Gimg.*(Gimg > 0);
    
    
        Gimgd = imadjust(Gimg,cdParams.imadjust_limits); %regular contrast
    
        Gimgfine = imadjust(Gimg, cdParams.fineimadjust_limits); %super boosted contrast
        %frame = str2double(files(imgnumb).name(p.frameIdInd:p.frameIdInd+3));
        

    end
   
%% initialize data structure
   
    if cdParams.figverbose
        figure(1); %makes a lot of figures, override onto 1
        
        imshow(Gimg)
        title('Gimg')

    end
    
    frame = str2double(files(imgnumb).name(fileParams.frameIdInd:fileParams.frameIdInd+3));
    

    %% initialize data structure
    %pData = readmatrix([p.topDir, 'adjacencyMatrix.m',centersfile(imgnumb).name]); %,"NumHeaderLines", 1); %Read Position data from centers file
    %ind = centersdata(:,1)==frame
    pData = centersdata(centersdata(:,1)==frame,:);
   
    if ~isempty(pData)

        N = size(pData,1);

        particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[],  'edge', 0);
        for n = 1:N %Bookkeeping from centers-tracked
            particle(n).id= pData(n,2);
            particle(n).x = pData(n,3);
            particle(n).y = pData(n,4);
            particle(n).r = round(pData(n,5));
            particle(n).edge = pData(n, 6);
            particle(n).rm = particle(n).r*cdParams.metersperpixel;
            particle(n).fsigma = cdParams.fsigma;
        end

        if cdParams.figverbose
            imshow(Gimgfine);
            %viscircles([pData(:,3),pData(:,4)],pData(:,5));
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
                particle(n).f = particle(n).g2/cdParams.g2cal; %saving some particle scale features
            else
                error('badimage!!')

            end
        end
        %% look at neighbours

        xmat = pData(:,3);
        ymat = pData(:,4);
        rmat = pData(:,5);

        rmats = rmat; %Saves our radius matrix for later

        dmat = pdist2([xmat,ymat],[xmat,ymat]); %Creates a distance matrix for particle center locations
        rmat = rmat + rmat'; %Makes a combination of radii for each particle

        friendmat = dmat < (rmat + cdParams.dtol) & dmat~=0; %Logical "friend" matrix

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

                        
            if cdParams.figverbose
            plot(x, y, 'LineWidth', 2)
            title('neighbour candidates')
            hold off
            end
            
            [contactG2p, contactIp] = contactspot(x,y,r, cdParams.CR, Gimg, maskCR);

            if(contactG2p(1) > cdParams.contactG2Threshold && contactG2p(2) > cdParams.contactG2Threshold)

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

            else %we try the more refined method of contact detection
       
            %find peaks in intensity for each particle, record the value of
            %the peak and the angular location relative to the x axis of
            %each particle
            [~, locs] = peakfinder(x(1), y(1), r(1), f1(l), Gimgfine, cdParams.minpeakheight, cdParams.minpeakprominence, cdParams.minpeakprom_main, cdParams.padding, cdParams.roach);
            [~, locs2] = peakfinder(x(2), y(2), r(2), f2(l), Gimgfine, cdParams.minpeakheight, cdParams.minpeakprominence, cdParams.minpeakprom_main, cdParams.padding, cdParams.roach);
    
            %compare the locations to whatever the nominal angle is (center
            %to center)
            
            nominalAngle = atan2(y(2)-y(1), x(2)-x(1));
            nominalAngleEXP = [nominalAngle-2*pi, nominalAngle, nominalAngle+2*pi];

            [angle] = ismembertol(nominalAngleEXP,locs',  pi/6,'DataScale', 1);

            nominalAngle2 = atan2(y(1)-y(2), x(1)-x(2));
            nominalAngleEXP = [nominalAngle2-2*pi, nominalAngle2, nominalAngle2+2*pi];

            [angle2] = ismembertol(nominalAngleEXP,locs2',  pi/6,'DataScale', 1);
        
%             if f1(l) == particleNumber1 && f2(l) == particleNumber2
%                 error()
%             end
            if ~isempty(angle(angle==1)) && ~isempty(angle2(angle2==1))

                particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
                particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p(1); %remember the g2 value of the current contact area
                particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp(1);
                particle(f1(l)).color(particle(f1(l)).z)='y';
                particle(f1(l)).neighbours(particle(f1(l)).z) = particle(f2(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
                particle(f1(l)).betas(particle(f1(l)).z) = nominalAngle; %the contact angle to particle m is now noted in the particle l datastructure
                particle(f2(l)).z= particle(f2(l)).z +1; %increase coordination number
                particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p(2); %remember the g2 value of the current contact area
                particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp(2);
                particle(f2(l)).color(particle(f2(l)).z)='y';
                particle(f2(l)).neighbours(particle(f2(l)).z) = particle(f1(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
                particle(f2(l)).betas(particle(f2(l)).z) = nominalAngle2;
            
            
            end
        
        
        
        
%             if isempty(particle(f1(l)).contactPos)
%                 particle(f1(l)).contactPos = locs;
%                 particle(f1(l)).contactInt = pks;
%             end
%             if isempty(particle(f2(l)).contactPos)
%                 particle(f2(l)).contactPos = locs2;
%                 particle(f2(l)).contactInt = pks2;
%             end
        
            end
    

    
    
        
            




        end
        %%

%         %Check if any of the walls is a neighbour as well
% 
%         circs = [[particle.y]', [particle.x]', [particle.r]', [particle.edge]']; %Makes a circs matrix from old matrices
% 
%         
%         for ep = 1:length(particle)
%             if circs(ep,4) == 1
%                 contacts = atan2(circs(ep, 1), circ(ep,2));
%             elseif circs(ep,4) == -1
%                 contacts = pi;
%             end
%             if particle(ep).edge ~=0
%                 x = particle(ep).x;
%                 y = particle(ep).y;
%                 r = particle(ep).r;
%                 for c =1:length(contacts) %technically doesn't need to be a loop but we will keep it for alternate scenarios
%                     [contactG2p, contactIp]= contactspotwall(x, y, r, cdParams.CR, contacts(c),Gimg, maskCR);
%                     if(contactG2p > cdParams.contactG2Threshold)
%                         particle(ep).z= particle(ep).z +1; %increase coordination number
%                         particle(ep).contactG2s(particle(ep).z)=contactG2p;
%                         particle(ep).contactIs(particle(ep).z)=contactIp;
%                         particle(ep).neighbours(particle(ep).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
%                         particle(ep).betas(particle(ep).z) = contacts(c); %the contact angle to the wall is now noted in the particle l datastructure
%                         particle(ep).color(particle(ep).z)='g';
%                         %     else
%                     end
%                 end
%              %end
%             end


        %end
    end


if cdParams.figverbose 
    h3 = figure(20);
    hAx1 = subplot(1,1,1,'Parent', h3);
    imshow(Gimg, 'Parent', hAx1);
    hold (hAx1, 'on');
    for n = 1:length(particle)
        particle(n).id;
        
        z = particle(n).z; %get particle coordination number
        if (z>0) %if the particle does have contacts
            for m = 1:z %for each contact
                %draw contact lines
                lineX(1)=particle(n).x;
                lineY(1)=particle(n).y;
                lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
                lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));
                viscircles([lineX(1), lineY(1)], particle(n).r, 'color', 'blue');
                viscircles([lineX(1) + (particle(n).r-cdParams.CR) * cos(particle(n).betas(m)) lineY(1) + (particle(n).r-cdParams.CR) * sin(particle(n).betas(m))], cdParams.CR, 'color', 'white');
                plot(hAx1, lineX, lineY,particle(n).color(m),'LineWidth',2);
            end
        end
        %text(hAx1, particle(n).x, particle(n).y, num2str(particle(n).id), 'Color', 'y')
    end
    drawnow;
    hold off;
    if verbose
    saveas(h3,fullfile(fileParams.topDir, fileParams.contactDir,['Contacts_' files(imgnumb).name(1:end-4)]),'jpg') %save fig(20) 
    end
end %PEGS wrapper verbose

if verbose
disp([num2str(sum([particle.z])), ' contacts detected'])
end

%save updated particle contact info
savename = strrep(files(imgnumb).name, 'warped.tif', '_contacts.mat');
save(fullfile(fileParams.topDir, fileParams.contactDir, savename),'particle')

end %loop imgnumb

%% save parameters 

fields = fieldnames(cdParams);
for i = 1:length(fields)
    fileParams.(fields{i}) = cdParams.(fields{i});
end


fileParams.time = datetime("now");
fields = fieldnames(fileParams);
C=struct2cell(fileParams);
cdParams = [fields C];

writecell(cdParams,fullfile(fileParams.topDir, fileParams.contactDir,'contactDetect_params.txt'),'Delimiter','tab')
if verbose 
        disp('done with contactDetect()');
end

out = true; %if function ran completely, return true to PEGS wrapper
end %function contact detect

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




function [profile] = contactfind(croppedImg, r)
    
    [a,b] = size(croppedImg);
    [X, Y] = meshgrid( (1:b)-r, (1:a)-r);
    R = sqrt(X.^2 + Y.^2);
    m2 = (3*r/4<R&R<r-1);

    %maskedImg = double(croppedImg).*m2; %if you want to see what the
    %masked image looks like
    values = double(croppedImg(m2));


    theta = atan2(Y, X);
    angles = theta(m2);
    combo = [angles, values];
    profile = sortrows(combo);

    profile = smoothdata(profile,1,'sgolay', length(profile)/25);

    profX = [profile(1:uint16(length(profile)/4),1)+(2*pi), profile(1:uint16(length(profile)/4), 2)];
    profile = [profile; profX];
    profile = sortrows(profile);
end

function [finalpks, finallocs] = peakfinder(x, y, r, ind, Gimgfine, minpeakheight, minpeakprominence,minpeakprom_main, padding, roach )
        global particleNumber1 particleNumber2
        croppedImg = (Gimgfine(round(y-r-padding):round(y+r+padding),round(x-r-padding):round(x+r+padding)));
        
        [profile] = contactfind(croppedImg, r-1);
        
        if any(profile(:,1)>minpeakheight)
            [pkints, locints] = findpeaks(profile(:,2),profile(:,1), "MinPeakHeight", minpeakheight,'MinPeakProminence', minpeakprominence, 'MaxPeakWidth', pi);
        else
            pkints = [];
            locints = [];
        end
        if any(profile(:,1)>minpeakprom_main)
            [pks,locs] =findpeaks(profile(:,2), profile(:,1),'MinPeakProminence',minpeakprom_main,'MaxPeakWidth', pi);
        else
            pks = [];
            locs = [];
        end
        pks = [pks;pkints];
        locs = [locs;locints];
        inds = find(locs>(pi));
        locs(inds) = locs(inds)-2*pi;

        [a, in] = uniquetol(locs, pi/12);
        [b] = uniquetol(locs, pi/12,'highest');
        finallocs = mean([a, b], 2);
        finalpks = pks(in);
        if roach == true
        if ind == particleNumber1 || ind == particleNumber2
            figure1 = figure;
            subplot(2,1,1, 'Parent', figure1)

       
        
            imshow(croppedImg);
            title(num2str(ind));
        
            hold on;
            contactLoc2 = [r+r.*cos(locs), r+r.*sin(locs)];
            contactLoc = [r+r.*cos(finallocs), r+r.*sin(finallocs)];
            %plot(r+r.*cos(nominalAngle), r+r.*sin(nominalAngle-pi),'yo')
            if ~isempty(contactLoc2)
            plot(contactLoc2(:,1), contactLoc2(:,2), 'ro');
            plot(contactLoc(:,1), contactLoc(:,2), 'bo');
            plot(r, r, 'gx')
            end
            subplot(2,1,2, 'Parent', figure1);
            plot(profile(:,1), profile(:,2));
            title(num2str(ind))
            hold on;
            plot(locs, pks, 'o');
            plot(finallocs, pks(1:length(finallocs)), 'o');
        
        
        end
        end
%         subplot(2,3,6, 'Parent', figure1)
%         plot(c);
%         hold on;
%         plot(loc, peaks, 'o');
        
end

function params = setupParams(params)
%% thresholding limits (setting defaults if none listed in wrapper)

if ~isfield(params,'metersperpixel') 
    params.metersperpixel = .007/160;
end

if ~isfield(params,'fsigma') %photoelastic stress coefficient
    params.fsigma = 140;
end

if ~isfield(params,'g2cal') %calibration Value for the g^2 method, can be computed by joG2cal.m (PEGS 1.0 version)
    params.g2cal = 100;
end

if ~isfield(params,'dtol') %how far away can the outlines of 2 particles be to still be considered Neighbors
    params.dtol = 10;
end

if ~isfield(params,'contactG2Threshold') %sum of g2 in a contact area larger than this determines a valid contact
    params.contactG2Threshold = 0.5; 
end

if ~isfield(params,'CR') %contact radius over which contact gradient is calculated
    params.CR = 10;
end

if ~isfield(params,'imadjust_limits') %adjust contrast in green channel
    params.imadjust_limits = [0,.65];
end

if ~isfield(params,'rednormal') %fractional amount to subtract the red channel from the green channel (Rimg/rednormal)
    params.rednormal = 2;
end


if ~isfield(params,'figverbose') %show all the figures as you go
    params.figverbose = true;
end
end
