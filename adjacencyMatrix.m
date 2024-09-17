function adjacencyMatrix(p, f, verbose)
directory = p.topDir
if not(isfolder(append(directory,'adjacency'))) %make a new folder with warped images
    mkdir(append(directory,'adjacency'));
end


files = dir([directory,'solved/',p.imgReg(1:end-4),'solved.mat']) %which files are we processing ?
nFrames = length(files) %how many files are we processing ?%snFrames = 92
if nFrames ==0
disp(['wrong spot:',directory,'solved/',p.imgReg(1:end-4),'solved.mat'])
return
end

%PARAMETERS NEEDED TO RUN THIS SCRIPT ARE SET HERE
go = true;
fmin = 0.000001; %minimum force (in Newton) to consider a contact a valid contact
fmax = 1000; %maximum force (in Newton) to consider a contact a valid contact
emax = 2800; %maximum fit error/residual to consider a contact a valid contact
fs=16; %plot font size
%verbose = False; %make lots of plots as we go



%%
%Global Metrics will be stored in these structures
%allContacts = struct('fAbs',0,'fNorm',0,'fTan',0); %data structure to store information about contacts
%aID = 1; %Global contact counter over all contacts in all cycles.

if go==true
for cycle = 1:nFrames %loop over these cycles 
    
    clearvars particle;
    clearvars contact;
    
    %input filnames
    peOutfilename = files(cycle).name %input filename 
    camImageFileName = [directory, 'warpedimg/',peOutfilename(1:end-19),'.tif'];  %adjusted force image filename
    
    % NO PARAMETERS SHOULD BE SET BY HAND BELOW THIS LINE

    %check if the data we want to read exists
    %if it does, load it, else abort
    if ~(exist([directory, 'solved/',peOutfilename], 'file') == 2) %if the file we try to open does not exist
        disp(['File not Found:', peOutfilename]); %complain about it
        return %and end the execution of this script
    else
        pres = load([directoryini, 'solved/', peOutfilename]); %read peDiscSolve ouput
        particle = pres.pres;
        NN = length(particle);
        IDN = max([particle.id]);
    end


    %DATA EVALUATION AND ANALYSIS STARTS HERE

    %particle(1:size(data,1)) = struct('id',0,'x',0,'y',0,'z',0,'fx',0,'fy',0); %data structure to store particle information
    contact = struct('id1',0,'id2',0,'x',0,'y',0,'fAbs',0,'fNorm',0,'fTan',0,'alpha',0,'beta',0,'contactX',0,'contactY',0,'error',0); %data structure to store information about contacts
    cID = 1; %contact counter
    %A = zeros(NN); %empty binary adjacency matrix
    W = NaN(IDN); %empty force weighted adjacency matrix
    N = NaN(IDN); %empty normal force weighted adjacency matrix
    T = NaN(IDN); %empty tangential force weighted adjacency matrix
    %P = NaN(NN);
    for n = 1:NN %for each particle
        err = particle(n).fitError; %get fit error 
        
        r = particle(n).r; %get particle radius in pixel
        
        if ~isempty(particle(n).neighbours) % particle is in contact
            contacts = particle(n).neighbours; %get IDs of all contacting particles
            betas = particle(n).betas+pi; %get the beta angle (position of contact point) associated with each contact
            forces = particle(n).forces; %get the force associated with each contact
            alphas = particle(n).alphas; %get the alpha angle (direction of force) associated with each contact
            
            for m=1:length(forces) %for each contact
		
                %if(forces(m) > fmin && err < emax && forces(m) < fmax) %is this a valid contact ?
                if (forces(m) < fmax && forces(m)> 0 && particle(n).color(m) ~= 'y')% && err <emax)
                    %put information about the first particle involved in this
                    %contact in the corresponding particle structure vector

                    %ideally the accumulated fx and fy should be zero, that is
                    %the particle is in force balance
                    %particle(n).fx = particle(n).fx + forces(m) * cos(betas(m)-pi); %x component of total force vector %CHECK AGAIN IF THIS IS GEOMETRICALLY CORRECT
                    %particle(n).fy = particle(n).fy + forces(m) * sin(betas(m)-pi); %y component of total force vector %CHECK AGAIN IF THIS IS GEOMETRICALLY CORRECT
                    %particle(n).z = particle(n).z+1; %increment the real contact number for the current particle

                    %put all the information about this contact
                    %into the contact struct vector
%                     contact(cID).id1 = particle(n).id; %first particle involved in this contact
                    targetid = contacts(m);
                    ids = [particle.id];
                    tind1m = ids == targetid;
                    tind1 = find(tind1m);
                    contact(cID).id2 = targetid; %second particle involved in this contact 
                    contact(cID).x = particle(n).x;
                    contact(cID).y = particle(n).y;
                    contact(cID).fAbs = forces(m); %absolute force
                    contact(cID).fNorm = forces(m)*cos(alphas(m)); %normal force (see Eq. 4.16)
                    contact(cID).fTan = forces(m)*sin(alphas(m)); %tangential force
                    contact(cID).alpha = alphas(m); %the alpha angle (direction of force) associated with this contact
                    contact(cID).beta = betas(m); %the beta angle (position of contact point) associated with this contact  
                    contact(cID).contactX =  r * cos(betas(m)-pi); %x component of vector to contact point
                    contact(cID).contactY =  r * sin(betas(m)-pi); %y component of vector to contact point
                    contact(cID).error = err; %fit error for this particle (the first particle in the contact)
% 
                    cID = cID + 1; %increment contact counter
%                     
%                     allContacts(aID).fAbs = forces(m); %absolute force
%                     allContacts(aID).fNorm = forces(m)*cos(alphas(m)); %normal force (see Eq. 4.16)
%                     allContacts(aID).fTan = forces(m)*sin(alphas(m)); %tangential force
%                     
%                     aID = aID+1;

                    %build some adjacency matrices
                    %if (contacts(m)>0) %correct for negative contact IDs in peDiscsolve, i.e.non-wall contacts only
                         %A(n,contacts(m)) = 1; %mark contact in the binary adjacency matrix
                         W(particle(n).id,particle(tind1).id) = real(forces(m)); %write the corrsponding force as a weight into an adjacency matrix
                         N(particle(n).id,particle(tind1).id) = real(forces(m))*cos(alphas(m)); %write the corrsponding normal force as a weight into an adjacency matrix
                         T(particle(n).id,particle(tind1).id) = real(forces(m))*sin(alphas(m)); %write the corrsponding tangential force as a weight into an adjacency matrix
                         %P(n, tind1) = [x, y, rm]
                    %end
                end
            end
        end
       
    end
    list = [];
    frameid = str2num(files(cycle).name(frameidind:frameidind+3));
    d = ~isnan(T);
    [row , col] = find(d==1);
    ind=sub2ind(size(T),row,col);
    
%     list_T = [row , col , T(row,col) , N(row,col) , Theta(row,col)];
    list = [ones(length(row),1).*frameid,row , col , T(ind) , N(ind)];
    writematrix(list, [directoryini, 'adjacency/',files(cycle).name(1:end-4),'-Adjacency.txt'] );
    %%
% length(nonzeros(W))
% edges = 10.^(-5:0.1:2);
% [N,edges] = histcounts(W,edges,'Normalization','countdensity');
% figure;
% g = histogram('BinEdges',edges,'BinCounts',N);
% set(gca, "Xscale", "log")
%     figure;
%     plot(nonzeros(W), '.')
%      drawnow;


        if verbose
        %figure(1)
            %read and display the original image used as input to peDisc
            %img = imcrop(imread(camImageFileName),[xoffset, yoffset, xsize, ysize]); %force image
            img = imread(camImageFileName); %force image
            imshow(img); hold on;
            colormap(gray);
            %plot the centers of all particles associated with a contact
            plot([contact.x],[contact.y],'or')
            %plot arrows from the centers of all particles associated with a contact to
            %the contact point
            f = [contact.fAbs];
            norm = max(max(f));
            shift = min(min(f));
            linewidths = 10*(f-shift+0.001)/norm;
            for m= 1:length(f)
                quiver([contact(m).x],[contact(m).y],[contact(m).contactX],[contact(m).contactY],0,'LineWidth',linewidths(m), Color='b')
            end
                
            %set font sizes and labels
            
            set(gca,'FontSize',fs);
            title('camera image','FontSize',fs);
            drawnow;
        end
    
 
end
end
%%
AdjFiles = dir([directoryini, 'adjacency/', fileNames(1:end-4), '-Adjacency.txt'])
nFrames = length(AdjFiles);   
posData = load([AdjFiles(nFrames-1).folder,'/', AdjFiles(nFrames -1).name]);
    
skipamount = length(posData)+2000; %I chose this as a result of my system size, could and should be altered based on your specific system and variability in finding particles
Adj_list = nan(nFrames*skipamount, 5);

for frame = 1:nFrames
        frame
        %posData = dlmread([directory, datafiles(n).name]);
    
        posData = load([AdjFiles(frame).folder, '/', AdjFiles(frame).name]);
        frameid = frame
        if length(posData) > skipamount
            error(['up the skipamount by', num2str(length(posData)-skipamount)])
            break
        end
        Adj_list((frame-1)*skipamount+1:(frame-1)*skipamount +length(posData),:) = posData;
end        
Adj_list(any(isnan(Adj_list),2),:)=[];


%fram number, particle 1, particle id 2, tangential force, normal force
dlmwrite([directoryini,'Adjacency_list.txt'],Adj_list);
disp('Adjacency matrix built')


end

