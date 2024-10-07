function img = joForceImgFixed (z, f, alpha, beta, fsigma, rm, px, verbose, fixedF, fixedAlpha, fixedbeta, fb)
     if nargin < 12
        fb = true;
     elseif nargin ==12
         fb = false;

     end
    
     f = [f; fixedF];
     alpha = [alpha; fixedAlpha];
     beta = [beta, fixedbeta];
    %make sure the forces are balanced
    if fb == true
        [alpha,f] = forceBalance(f,alpha,beta);
    end
    %create an empty placeholder image for our result/return value
    img = zeros(px);

    %Create a scale that maps the diameter of the particle onto the image size
    xx = linspace(-rm, rm, px); 
    for x=1:px %loop throgh image width
        xRow=zeros(px,1); %placeholder for the current row beeing prcessed, parallel for makes it necessary to split up the result data and later combine it
        for y=1:px  %loop throgh image height
            if ((xx(x)^2+xx(y)^2)<=rm^2) %check if we are actually inside
                xRow(y) = StressEngine(xx(x), xx(y), z, f, alpha, beta, fsigma, rm);  %call the StressEngine to compute the photoelastic response at each pixel
            end
        end
        img(x,:)=xRow; %consolidate processed row into the output image, necessary data shuffling to use parallel for
    end

    if(verbose)
        %plot the synthetic image each time it is generated (i.e. you will
        %see the fit converge on screen)
        subplot(1,2,2)
        imshow(img);
        drawnow
    end
end
