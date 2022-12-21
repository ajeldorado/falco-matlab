% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function handles = falco_plot_progress_testbed(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf)
    
    
    switch upper(mp.testbed)
        case 'HCST'
            handles = falco_plot_progress_hcst(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf);
        case 'GPCT'
            handles = falco_plot_progress_gpct(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf); 
        case 'DST'
            handles = falco_plot_progress_dst(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf); 
        case 'SIM'
            handles = falco_plot_progress(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf); 
        case 'OMC'
            handles = falco_plot_progress_omc(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf); 
        case 'IACT'
            handles = falco_plot_progress_iact(handles,mp,Itr,Inorm,Im_tb,DM1surf,DM2surf); 
        %Can put other testbeds or fancy models here
        %case 'WHATEVER'
            %normI = falco_plot_progress_whatever(mp,si); 
        otherwise
            error('Case not recognized.  Check falco_plot_progress_testbed.m to make sure this case is specified.');
    end


end %--END OF FUNCTION
