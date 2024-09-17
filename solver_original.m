function particle = solver_original(particle, userOptions )

maskradius = userOptions.maskradius;
scaling = userOptions.scaling;
fitoptions = userOptions.fitoptions;

% Number of particles in this frame
N = length(particle); 
    parfor n = 1:N
        display(['fitting force(s) to particle ',num2str(n)]); %status indicator
        if (particle(n).z > 0 )
            % Extract necessary information from particle structure
            fsigma = particle(n).fsigma;
            z = particle(n).z;
            forces = zeros(z,1);
            cg2s = sum(particle(n).contactG2s);
            beta = particle(n).betas;
            rm = particle(n).rm;
            template = particle(n).forceImage;
            template = imresize(template,scaling);
            px = size(template,1); 
            
            % Initial force and alpha values
            for i=1:z
                forces(i) = 2*particle(n).f*particle(n).contactG2s(i)/cg2s;
            end
            alphas = zeros(z,1);
            
            % Apply force balance to the initial guesses - unchanged from
            % PeGS1.0
            [alphas,forces] = forceBalance(forces,alphas,beta);

            % Create a circular mask
            cx=px/2;
            cy=px/2;
            ix=px;
            iy=px;
            r=maskradius*px;
            [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
            c_mask=((x.^2+y.^2)<=r^2);   

            %% Least squares fitting
            % Set up initial values
            p0 = zeros(2*z, 1);
            p0(1:z) = forces;
            p0(z+1:2*z) = alphas;

            % Fitting functions
            func = @(par) fringe_pattern_original(z, par(1:z),par(z+1:z+z), beta(1:z), fsigma, rm, px); 
            err = @(par) real(sum(sum( ( c_mask.*(template-func(par)).^2) ))); 
            p = lsqnonlin(err,p0,[],[],fitoptions);

            % Extract fitting outputs
            forces = p(1:z);
            alphas = p(z+1:z+z);
            fitError = err(p);
            
            % Generate an image with the fitted parameters
            imgFit = fringe_pattern_original(z, forces, alphas, beta, fsigma, rm, px*(1/scaling));
            
            % Redo force balance
            [alphas,forces] = forceBalance(forces,alphas,beta);
            
            % Store the new information in particle 
            particle(n).fitError = fitError;
            particle(n).forces = forces;
            particle(n).alphas = alphas;
            particle(n).synthImg = imgFit;
        end
    end
end