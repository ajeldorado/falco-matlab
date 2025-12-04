% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%     Calculate the complex-valued transmission of a 2-D mask.
% 
%     Calculates the thin-film complex transmission for the provided 2-D maps
%     of metal and dielectric thicknesses for a single wavelength.
% 
%     This function is a wrapper around falco_thin_film_material_def().
% 
%     Parameters
%     ----------
%     substrate : str
%         Name of the substrate material.
%     metal : str
%         Name of the metal used in the mask.
%     dielectric : str
%         Name of the dielectric used in the mask.
%     lam : float
%         Wavelength in meters.
%     aoi : flaot
%         Angle of incidence in degrees.
%     t_Ti : float
%         Titanium layer thickness in meters. Titanium is used in a uniform
%         thickness only between the substrate and the main metal to help
%         adhesion.
%     t_metal_map : array_like
%         2-D array of metal thicknesses in meters. This metal goes between the
%         titanium and dielectric layers.
%     t_diel_map : array_like
%         2-D array of dielectric thicknesses in meters.
%     d0 : float
%         Reference height for all phase offsets. Must be larger than the stack
%         of materials, not including the substrate. Units of meters.
%     pol : {0, 1, 2}
%         Polarization state to compute values for.
%         0 for TE(s) polarization,
%         1 for TM(p) polarization,
%         2 for mean of s and p polarizations
%     flagOPD : bool, optional
%         Flag to use the OPD convention. The default is False.
% 
%     Returns
%     -------
%     out_map : numpy ndarray
%         2-D complex transmission map for the provided layer thicknesses.

function out_map = falco_calc_complex_occulter(substrate, metal, dielectric, lam, aoi, t_Ti, ...
                                               t_metal_map, t_diel_map, d0, pol, varargin)

    % Optional Keyword Inputs
    flagOPD = false; %--Default value for OPD phase sign convention is false.
    icav = 0;             % index in cell array varargin
    while icav < size(varargin, 2)
        icav = icav + 1;
        switch lower(varargin{icav})
            case {'opd'}
                flagOPD  = true;       % Use the OPD phase sign convention.
            otherwise
                error('falco_thin_film_material_def: Unknown keyword: %s\n', varargin{icav});
        end
    end
    
    % Input checks
    Check.real_nonnegative_scalar(t_Ti)
    Check.two_dim_array(t_metal_map)
    Check.two_dim_array(t_diel_map)
    
    t_Ti_map = zeros(size(t_metal_map));
    t_Ti_map(t_metal_map > 10*eps) = t_Ti;
    
    % Put each vector as a column in a matrix and then keep only unique
    % rows to save computation time.
    t_mat = [t_diel_map(:), t_metal_map(:), t_Ti_map(:)];
    t_unique_mat = unique(t_mat, 'rows');

    t_diel_vec_short = t_unique_mat(:, 1);
    t_metal_vec_short = t_unique_mat(:, 2);
    t_Ti_vec_short = t_unique_mat(:, 3);

    Nshort = size(t_unique_mat, 1); 
    % tCoefShort = np.zeros(Nshort)
    % rCoefShort = np.zeros(Nshort)

    out_map = zeros(size(t_metal_map));
    for ii = 1:Nshort

        t_diel = t_diel_vec_short(ii);
        t_metal = t_metal_vec_short(ii);
        t_Ti_here = t_Ti_vec_short(ii);

        if flagOPD
            [tCoef, ~] = falco_thin_film_material_def(substrate, metal, dielectric, lam, aoi, ...
                t_Ti_here, t_metal, t_diel, d0, pol);
        else
            [tCoef, ~] = falco_thin_film_material_def(substrate, metal, dielectric, lam, aoi, ...
                t_Ti_here, t_metal, t_diel, d0, pol, 'opd');
        end

        thisRegion = (t_Ti_map == t_Ti_here) & (t_diel_map == t_diel) & (t_metal_map == t_metal);
        out_map(thisRegion) = tCoef;
        
    end

end
