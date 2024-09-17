function img = fringe_pattern_original(z, f, alpha, beta, fsigma, rm, px)

    %make sure the forces are balanced
    [alpha,f] = forceBalance(f,alpha,beta);

    %create an empty placeholder image for our result/return value
    img = zeros(px);

    %Create a scale that maps the diameter of the particle onto the image size
    xx = linspace(-rm, rm, px); 
    for x=1:px 
        xRow=zeros(px,1); 
        for y=1:px  
            if ((xx(x)^2+xx(y)^2)<=rm^2) 
                xRow(y) = stress_engine_original(xx(x), xx(y), z, f, alpha, beta, fsigma, rm);
            end
        end
        img(x,:)=xRow; 
    end
end