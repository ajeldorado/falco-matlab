% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate a traditional grayscale apodizer
%
% Written by G. Ruane on 2018-05-19. 
%
% INPUTS: 
%  PUPIL: 2-D square array with the pupil
%  apDiaSamps: Aperture diameter in samples
%  effIWA: Effective inner working angle 
%  effOWA: Effective outer working angle 
%  useGPU: (logical) Run FFTs on GPUs
%
% OUTPUTS:
%  APOD:     2-D square array for the apodizer mask. Cropped down to the smallest even-sized array with no extra zero padding. 

function APOD = falco_gen_tradApodizer(PUPIL,apDiaSamps,effIWA,effOWA,useGPU)

    lamOverD = 4; 
    Narray = 2^nextpow2(lamOverD*apDiaSamps);
    
	[rows0,cols0] = size(PUPIL);
    
    PUPIL = padOrCropEven(PUPIL,Narray);
    
    if(useGPU)
        PUPIL = gpuArray(PUPIL);
    end
    
    [X,Y]=meshgrid(-Narray/2:Narray/2-1);
    [~,RHO]=cart2pol(X,Y);
    
%     Q = 1 - exp(-(RHO/(effIWA*lamOverD)).^1000);
    Q = exp(-(RHO/(effOWA*lamOverD)).^1000) - exp(-(RHO/(effIWA*lamOverD)).^1000);
    Q = fftshift(Q);
    

    prevMetric = 1; 
    minThpt = 0.3;
    maxIts = 500;
    
    A = PUPIL;
    normI = max(max(abs(fft2(PUPIL)).^2));
    
    
    % Initial conditions 
    B = fft2(A);
    peakPSFthpt = max(abs(B(:)).^2/normI);
    currMetric = mean(abs(B(logical(Q))).^2/normI)/peakPSFthpt.^2;
    
    
	disp('Optimizing traditional grayscale apodizer...');
    it = 0;
    while(it < maxIts && peakPSFthpt>minThpt && currMetric<prevMetric )

        B = B.*(1-Q);
        A = abs(ifft2(B)).*round(PUPIL);
        
        A(A>1) = 1; 
%         A = A/max(A(:));
       
        B = fft2(A);
        
        peakPSFthpt = max(abs(B(:)).^2/normI);
        
        prevMetric = currMetric;
        currMetric = mean(abs(B(logical(Q))).^2/normI)/peakPSFthpt.^2;
        
        it = it + 1;
    end
        
	APOD = padOrCropEven(A,rows0);
        
    if(useGPU)
        APOD = gather(APOD);
    end
    
    disp(['Traditional grayscale apodizer returned with T=',num2str(peakPSFthpt*100,2),'% and M=',num2str(currMetric)]);


end %--END OF FUNCTION
