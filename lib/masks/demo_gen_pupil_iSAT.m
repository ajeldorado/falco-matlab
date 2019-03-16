% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear; 
addpath('segMirrorFunctions');

% With three rings of rafts and zero gap width
% the circumscribed diameter is 18.8 segment widths 
% the inscribed diameter is 15.3 segment widths 
% The 20-meter design is 20 meters from flat-to-flat along the vertical 
% with (17 segments) and raft diameter (flat-to-flat) of 3.5 meters. 

% With small gaps, the circumscribed diameter is ~19 segment widths 

input.Nbeam = 1000;
input.Npad = 2^10;

raft_apDia_m = 3.5;% raft diameter in meters (20m design has 3.5m rafts) 

input.numRaftRings = 3;
input.raftDia = input.Nbeam/19*3;%samples
input.wGap = 6e-3/raft_apDia_m*input.raftDia;
input.raftGap = 30e-3/raft_apDia_m*input.raftDia;

% input.telDia = input.Nbeam;

% input.raftOffsets = zeros(hexSegMirror_numSegments( input.numRaftRings ),2);
% input.raftOffsets(2,:) = input.raftDia/raft_apDia_m*[-10e-3,10e-3];


RAFT = falco_gen_pupil_iSAT( input );

%%

[rows,cols]=size(RAFT);
[X,Y] = meshgrid(-rows/2:rows/2-1);
[~,RHO] = cart2pol(X,Y);

%%

figure(1);
imagesc(RAFT);
% imagesc(RAFT+(RHO<raft_apDia/2));
axis image;
set(gca,'ydir','normal');
