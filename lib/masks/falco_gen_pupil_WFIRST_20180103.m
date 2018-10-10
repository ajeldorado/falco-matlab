% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to generate a WFIRST CGI input pupil, ID # 20180103. (ID # is the
% date it was received.)
%  -Has non-circular, offset secondary mirror.
%  -Has 3.22% OD) strut scraper width.

function mask = falco_gen_pupil_WFIRST_20180103(Nbeam,centering,varargin)

% % % %--DEBUGGING ONLY: HARD-CODED INPUTS
% centering = 'pixel'; %--Output centering
% Nbeam = 250; % 250 %--Number of points across the resized pupil


% Set default values of input parameters
flagRot180deg = false;
%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'rot180'}
            flagRot180deg = true; % For even arrays, beam center is in between pixels.
        otherwise
            error('falco_gen_pupil_WFIRST_20180103: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

% centering_orig = 'interpixel'; % centering of the output array: pixel centering or interpixel centering

%--Load the WFIRST Pupil
% if(exist('CGI Entrance Pupil 76mm 4096 180103.bmp','file')~=2)
%     unzip('CGI Entrance Pupil 76mm 4096 180103.bmp.zip')
% end
% pupil0 = (imread('CGI Entrance Pupil 76mm 4096 180103.bmp'));
pupil0 = imread('pupil_WFIRST_CGI_20180103.png'); pupil0 = rot90(pupil0,2);


pupil1 = sum(pupil0,3);
pupil1 = pupil1/max(pupil1(:));
if(flagRot180deg)
    pupil1 = rot90(pupil1,2);
end

% figure(1); imagesc(pupil1); axis xy equal tight; colorbar;
% figure(11); imagesc(pupil1-fliplr(pupil1)); axis xy equal tight; colorbar;

%%--Resize
Npup = length(pupil1);
xs0 = ( -(Npup-1)/2:(Npup-1)/2 )/Npup; %--original coordinates, normalized to the pupil diameter. True for the 20180103 design, which is interpixel centered.
Xs0 = meshgrid(xs0);


switch centering
    case{'interpixel','even'}
%         mask = imresize(pupil1,[1 1]*(Nbeam),'bilinear');
        xs1 = ( -(Nbeam-1)/2:(Nbeam-1)/2 )/Nbeam;
        Xs1 = meshgrid(xs1);
        mask = interp2(Xs0,Xs0.',pupil1,Xs1,Xs1.','spline',0); %--interp2 does not get the gray edges as well, but imresize requires the Image Processing Toolbox
    case{'pixel','odd'}
%         mask = imresize(pupil1,[1 1]*(Nbeam+1),'bilinear'); %--Pupil is slightly larger than it should be because it is resized to (Nbeam+1) pixels across instead of being Nbeam pixels across with a half-pixel centering offset.
        xs1 = ( -(Nbeam)/2:(Nbeam)/2 )/Nbeam;
        Xs1 = meshgrid(xs1);
        temp = interp2(Xs0,Xs0.',pupil1,Xs1,Xs1.','linear',0); %--interp2 does not get the gray edges as well, but imresize requires the Image Processing Toolbox
        mask = zeros(Nbeam+2,Nbeam+2); %--Initialize
        mask(2:end,2:end) = temp; %--Fill in values
        %mask = padarray(mask,[1 1],0,'pre'); % requires the Image Processing Toolbox
end
% figure(12); imagesc(mask); axis xy equal tight; colormap gray; colorbar;

end %--END OF FUNCTION
