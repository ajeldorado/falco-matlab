% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Constrain DM commands with tied actuators and neighbor restrictions.
%

classdef ConstrainDM
   
methods (Static)

    
    function smoothed = constrain_dm(volts, flatmap, tie, vmax, vlat, vdiag, vquant, maxiter)
        % Given a DM setting, return one consistent with physical and user
        % mandated constraints.
        % 
        % Constrains each individual voltage to be in 0 <= v <= `vmax`.  Use of
        % defaults (min voltage=0V, `vmax`=100V) enforces clip requirements in DNG
        % 884740 and 884741.
        % 
        % Constrains each pair of laterally-adjacent actuators to be <= `vlat` after
        % subtraction of the DM flat map in `flatmap`. Use of default (`vlat`=50V)
        % enforces neighbor-rule requirements in DNG 884742 and 884743.
        % 
        % Constrains each pair of diagonally-adjacent actuators to be <= `vdiag`
        % after subtraction of the DM flat map in `flatmap`. Use of default
        % (`vdiag`=75V) enforces neighbor-rule requirements in DNG 1073291 and
        % 1073292.
        % 
        % Constrains all tied actuators (groups in the `tie` matrix with value > 0)
        % to have the same voltage.  Constrains all dead actuators (groups in the
        % `tie` matrix with value = -1) to be 0V.
        % 
        % Note that in practical use, the flatmap and tie matrices have certain
        % expectations on their format:
        %  - flatmap should be >= 0, <= vmax, have all ties at the same voltage,
        %   and have all dead actuators at 0V voltage.
        %  - tie should be -1, 0, or the integer range 1-N for some N (with no gaps)
        % These expectations will be enforced.
        % 
        % NOTE: this function assumes that the minimum commandable voltage is equal
        % to the dead actuator voltage == 0V, and this also corresponds to 0x0000 in
        % the DAC (so exactly the bottom value).  All of these things are true for
        % CGI in its current implementation, but if this were to be revisited in
        % another use case where these assumptions are not valid, then the logic here
        % needs to be revisited.  (In particular, dead actuator voltage < minimum
        % commandable voltage makes strange edge effects, as does the case where the
        % voltage corresponding to the DAC value 0x0000 at the lower edge is not 0V
        % but a dead actuator is.)
        % 
        % Arguments:
        %  volts: a 2D array of floating-point voltages.  This is the set of voltages
        %   which we are fixing before sending to the DM.
        %  flatmap: a 2D array of floating-point voltages, of the same size as
        %   `volts`.  This array represents a physically-flat DM surface.
        %  tie: a 2D array of integers, of the same size as `volts`, which can take
        %   on the values 0, -1, or consecutive integers 1 -> N.
        % 
        % Keyword Arguments:
        %  vmax: maximum commandable voltage, in volts.  Floating-point scalar, must
        %   be > 0. Defaults to 100V, which is the max voltage
        %   for CGI.
        %  vlat: maximum allowable voltage differential between laterally-adjacent
        %   actuators, in volts.  Floating-point scalar > 0.  Defaults to 50V, which
        %   is the CGI requirement.
        %  vdiag: maximum allowable voltage differential between diagonally-adjacent
        %   actuators, in volts.  Floating-point scalar > 0.  Defaults to 75V, which
        %   is the CGI requirement.
        %  vquant: smallest voltage step (1 LSB) which the DME electronics can
        %   produce.  Used to keep the constraints from being broken after EU->DN
        %   conversion. Floating-point scalar >= 0; using 0 is equivalent to not
        %   accounting for DAC discretization effects at all.  Defaults to 110/2^16,
        %   the CGI LSB.
        %  maxiter: number of times to iterate between smoothing and tying before
        %   giving up.  Smoothing and tying are both convergent and we do not
        %   expect to need this, but it seemed a reasonable safety measure
        %   against a having the 'while True' loop repeat indefinitely for some
        %   unforeseen corner case.  Defaults to 1000.
        % 
        % Returns:
        %  a constrained 2D array of voltages

        % Check argumemts
        Check.two_dim_array(volts)
        Check.two_dim_array(flatmap)
        Check.two_dim_array(tie)
        if ~isequal(size(volts), size(flatmap))
            error('volts and flatmap must have the same shape')
        end
        if ~isequal(size(volts), size(tie))
            error('volts and tie must have the same shape')
        end

        % Check keyword arguments
        Check.real_positive_scalar(vmax)
        Check.real_positive_scalar(vlat)
        Check.real_positive_scalar(vdiag)
        Check.real_nonnegative_scalar(vquant)
        Check.positive_scalar_integer(maxiter)

        VMIN = 0;
        if vmax <= VMIN
            error('VMIN must be < vmax')
        end

        % enforce tie and flat formatting
        if ~ConstrainDM.checktie(tie)
            error('tie must have values 0, -1, or consecutive integers 1 to N')
        end
        if ~ConstrainDM.checkflat(flatmap, VMIN, vmax, tie)
            error(['flatmap must be <= vmax, >= VMIN, have all tied ' 
                    'actuators tied already, and have all dead ' 
                     'actuators = 0V'])
        end


        % Run initial smoothing
        dmflat = flatmap;
        smoothed = ConstrainDM.dmsmooth(volts, vmax, vquant, vlat, vdiag, dmflat);

        % Dummy array to initialize while loop
        tied = []; % None

        % Loop tie and smooth unitl they converge
        ii = 0;
        while ~isequal(tied, smoothed) && (ii < maxiter)
            tied = ConstrainDM.tie_with_matrix(smoothed, tie);
            smoothed = ConstrainDM.dmsmooth(tied, vmax, vquant, vlat, vdiag, dmflat);
            ii = ii + 1;
        end
        if ii >= maxiter
            error('maxiter exceeded for constrain_dm')
        end
    end
    
    
    function dmout = dmsmooth(dmin, vmax, vquant, vneighbor, vcorner, varargin)
        %     Modify a DM setting so that it obeys neighbor rules and voltages limits
        % 
        %     This implementation is conservative: it will produce settings that always
        %     obey neighbor rules after the floating-point DM specification has been
        %     quantized onto the n-bit DAC.  All neighbor-rule and voltage bounds will be
        %     at least 2 LSBs away from their threshold.
        % 
        %     This method will produce a result that is safe, but not necessarily
        %     optimal for any particular definition of optimality.
        % 
        %     This method is idempotent: running the output through the function a second
        %     time with identical settings will give an identical output to the original
        %     function.
        % 
        %     NOTE: this function assumes that the minimum commandable voltage is equal
        %     to the dead actuator voltage == 0V, and this also corresponds to 0x0000 in
        %     the DAC (so exactly the bottom value).  All of these things are true for
        %     CGI in its current implementation, but if this were to be revisited in
        %     another use case where these assumptions are not valid, then the logic here
        %     needs to be revisited.  (In particular, dead actuator voltage < minimum
        %     commandable voltage makes strange edge effects, as does the case where the
        %     voltage corresponding to the DAC value 0x0000 at the lower edge is not 0V
        %     but a dead actuator is.)
        % 
        %     Arguments:
        %      dmin: 2D square array of actuator settings to be fixed, in floating-point
        %       volts.
        %      vmax: Maximum voltage permitted, in floating-point volts.  The output
        %       will be no larger than the largest number <= vmax which is an integer
        %       multiple of vquant (unless vquant == 0, in which case it will be exact)
        %      vquant: quantization step in the DAC (LSB).  Used to determine how close
        %       to the vmin/vmax/vneighbor/vcorner thresholds a setting can be before it
        %       triggers a correction.  Must be >= 0; if vquant=0, thresholds will be
        %       exact to floating-point error
        %      vneighbor: Permitted voltage difference between lateral neighbors, in
        %       floating-point volts.  Must be >= 0.
        %      vcorner: Permitted voltage difference between diagonal neighbors, in
        %       floating-point volts.  Must be >= 0.
        % 
        %     Keyword Arguments:
        %      dmflat: 2D array of the same size as dmin, or None. Defines the phase-flat
        %       voltages for neighbor rule checking.  If None, assumes 0V, consistent
        %       with unpowered polish (CGI baseline).  Neighbor rule must only be
        %       maintained with respect to a phase-flat array.  dmflat must be <= vmax,
        %       >= vmin at all points.
        % 
        %     Returns:
        %      a 2D array of the same size as dmin

        % Check inputs
        Check.two_dim_square_array(dmin)
        Check.real_scalar(vmax)
        Check.real_nonnegative_scalar(vquant)
        Check.real_nonnegative_scalar(vneighbor)
        Check.real_nonnegative_scalar(vcorner)

        % Optional arguments
        dmflat = zeros(size(dmin)); % initial DM command map
        if length(varargin) == 1
            dmflat = varargin{1};
            
        elseif length(varargin) > 1
            error('Too many inputs.')
        end
        if isempty(dmflat)
            dmflat = zeros(size(dmin));
        end
        Check.two_dim_square_array(dmflat);
        

        VMIN = 0;
        N = size(dmin, 1);

        % Value Checks
        if vmax <= VMIN
            error('VMIN must be less than vmax')
        end
        if ~isequal(size(dmflat), size(dmin))
            error('dmflat should be same size as dmin')
        end
        if any(dmflat < VMIN)
            error('dmflat should be >= VMIN')
        end
        if any(dmflat > vmax)
            error('dmflat should be <= vmax')
        end

        dmout = double(dmin);
        margin = 2 * vquant;
        if vquant == 0
            vmax_q = vmax;
        else
            vmax_q = vquant * floor(vmax/vquant);
        end

        % Fix DM upper and lower bounds
        dmout(dmout < VMIN) = VMIN;
        dmout(dmout > vmax_q) = vmax_q;

        % If nside is even, exclude the last member in the second correction step.
        if mod(size(dmout(1)), 2) == 0
            endVal = N-1;
        else
            endVal = N;
        end

        % Fix internal neighbor rule violations
        isDone = false;
        dmout = dmout - dmflat; % subtract only while running
        while ~isDone
            % Loop has four essentially similar sets of operations
            % Detailed comments only for first, but rest follow same pattern

            % Subtract adjacent columns (horizontal)
            % Fix the array in two steps.
            % First, subtract every other column so that difference is calculated
            % in pairs. For example, for [1, 2, 3, 4, 5, 6], take 2-1, 4-3, 6-5.
            c1x = 1:2:(N-1); % slice(0, -1, 2)
            c2x = 2:2:N; % slice(1, None, 2)
            dmx = dmout;
            [dmx, maskx] = ConstrainDM.fix_nr(dmx, vneighbor, margin, 1:N, c1x, 1:N, c2x);

            % Next, follow the same procedure but excluding the first column, and
            % if nside is even exclude the last column also.
            % For example, for [1, 2, 3, 4, 5, 6] take 3-2, 5-4.
            dmx_ = dmout;
            dmTemp = dmx_(:, 2:endVal);
            [nRows, nCols] = size(dmTemp);
            c1x = 1:2:(nCols-1); % slice(0, -1, 2)
            c2x = 2:2:nCols; % slice(1, None, 2)
            [temp, maskx_] = ConstrainDM.fix_nr(dmTemp, vneighbor, margin, 1:nRows, c1x, 1:nRows, c2x);
            dmx_(:, 2:endVal) = temp;
            % Combine the two adjustements together by applying them one at a time.
            dmout(maskx) = dmx(maskx);
            % Only apply second step adjustments to elements that were not adjusted
            % in the first step. This will avoid adjusting the same elements twice.
            dmout(~maskx) = dmx_(~maskx);


            % Subtract adjacent rows (vertical)
            r1y = 1:2:(N-1); % slice(0, -1, 2)
            r2y = 2:2:N; % slice(1, None, 2)
            dmy = dmout;
            [dmy, masky] = ConstrainDM.fix_nr(dmy, vneighbor, margin, r1y, 1:N, r2y, 1:N);

            
            dmy_ = dmout;
            dmTemp = dmy_(2:endVal, :);
            [nRows, nCols] = size(dmTemp);
            r1y = 1:2:(nRows-1); % slice(0, -1, 2)
            r2y = 2:2:nRows; % slice(1, None, 2)
            
            [temp, masky_] = ConstrainDM.fix_nr(dmTemp, vneighbor, margin, r1y, 1:nCols, r2y, 1:nCols);
            dmy_(2:endVal, :) = temp;
            dmout(masky) = dmy(masky);
            dmout(~masky) = dmy_(~masky);


            % Subtract adjacent rows (diagonal right)
            r1xy = 1:2:(N-1); % slice(0, -1, 2)
            r2xy = 2:2:N; % slice(1, None, 2)
            c1xy = 1:(N-1); % slice(0, -1)
            c2xy = 2:N; % slice(1, None)
            dmxy = dmout;
            [dmxy, maskxy] = ConstrainDM.fix_nr(dmxy, vcorner, margin, r1xy, c1xy, r2xy, c2xy);

            dmxy_ = dmout;
            dmTemp = dmxy_(2:endVal, :);
            [nRows, nCols] = size(dmTemp);
            r1xy = 1:2:(nRows-1); % slice(0, -1, 2)
            r2xy = 2:2:nRows; % slice(1, None, 2)
            c1xy = 1:(nCols-1); % slice(0, -1)
            c2xy = 2:nCols; % slice(1, None)
            [temp, maskxy_] = ConstrainDM.fix_nr(dmTemp, vcorner, margin, r1xy, c1xy, r2xy, c2xy);
            dmxy_(2:endVal, :) = temp;
            dmout(maskxy) = dmxy(maskxy);
            dmout(~maskxy) = dmxy_(~maskxy);


            % Subtract adjacent rows (diagonal left)
            r1yx = 1:2:(N-1); % slice(0, -1, 2)
            r2yx = 2:2:N; % slice(1, None, 2)
            c1yx = 2:N; % slice(1, None)
            c2yx = 1:(N-1); % slice(0, -1)
            dmyx = dmout;
            [dmyx, maskyx] = ConstrainDM.fix_nr(dmyx, vcorner, margin, r1yx, c1yx, r2yx, c2yx);

            dmyx_ = dmout;
            dmTemp = dmyx_(2:endVal, :);
            r1yx = 1:2:(nRows-1); % slice(0, -1, 2)
            r2yx = 2:2:nRows; % slice(1, None, 2)
            c1yx = 2:nCols; % slice(1, None)
            c2yx = 1:(nCols-1); % slice(0, -1)
            [temp, maskyx_] = ConstrainDM.fix_nr(dmTemp, vcorner, margin, r1yx, c1yx, r2yx, c2yx);
            dmyx_(2:endVal, :) = temp;
            dmout(maskyx) = dmyx(maskyx);
            dmout(~maskyx) = dmyx_(~maskyx);


            % If any of them had violations, go around again
            isDone = ...
                (sum(maskx(:)) == 0) && (sum(maskx_(:)) == 0) && ...
                (sum(masky(:)) == 0) && (sum(masky_(:)) == 0) && ...
                (sum(maskxy(:)) == 0) && (sum(maskxy_(:)) == 0) && ...
                (sum(maskyx(:)) == 0) && (sum(maskyx_(:)) == 0);

        end

        dmout = dmout + dmflat; % reinsert

        % recheck clip, clip is applied without flat
        dmout(dmout < VMIN) = VMIN;
        dmout(dmout > vmax_q) = vmax_q;

    end


    function out = sign_array(x)
        % Fast scalar sign function for arrays
        out = zeros(size(x));
        out(x > 0) = 1;
        out(x < 0) = -1;
    end


    function [dmout, fix_mask] = fix_nr(dmin, vneighbor, margin, r1, c1, r2, c2)
        % Take a pair of columns and/or rows,
        % check for neighbor rule violations between them,
        % and fix them.

        dmout = dmin;

        % Subtract adjacent rows/cols
        diff = dmout(r2, c2) - dmout(r1, c1);

        % Find if any neighbor rule violations exist. Use mask to exclude good
        % neighbors from correction.
        diff_mask = abs(diff) > vneighbor - margin;
        diff = diff .* diff_mask;

        % Split excess in half, with a little extra margin to keep
        % numerical errors from making this reappear.  2x margin to
        % overshoot the correction slightly so this doesn't pop up
        % again.
        delta = (diff - ConstrainDM.sign_array(diff) .* (vneighbor - 2*margin))/2;

        % Half on each side of the violation
        dmout(r1, c1) = dmout(r1, c1) + delta;
        dmout(r2, c2) = dmout(r2, c2) - delta;

        % Return a mask for all the elements that were fixed
        fix_mask = dmout ~= dmin;

    end    
    
    
    function dmtied = tie_with_matrix(volts, tie)
        % Tie specified actuators to single value. The value is the mean value of
        % all actuators.  This uses a matrix with to indicate which actuators are
        % tied together, with a specific format:
        %  - 0 indicates no tie
        %  - -1 indicates dead (0V)
        %  - any other integer 1->N indicates a tie group; all actuators in that
        %   group will be assigned the mean voltage across that set of actuators.
        % 
        % Arguments:
        %  volts: a 2D array of floating-point voltages.  This is the set of voltages
        %   which we are fixing before sending to the DM.
        %  tie: a 2D array of integers, of the same size as `volts`, which can take
        %   on the values 0, -1, or consecutive integers 1 -> N.
        % 
        % Returns:
        %  Output DM array with tied values.

        % dimensionality checks
        Check.two_dim_array(volts)
        Check.two_dim_array(tie)
        if ~isequal(size(volts), size(tie))
            error('ValueError:SizeMismatch', 'volts and tie must have the same shape')
        end

        % enforce tie formatting
        if ~ConstrainDM.checktie(tie)
            error('tie must have values 0, -1, or consecutive integers 1 -> N')
        end

        % Loop through each DM mask in data cube
        dmtied = volts;
        tienumset = sort(unique(tie(:).')); % set(tie.ravel())
        for tienum = tienumset
            if tienum == 0
                % not tied
            elseif tienum == -1
                % dead actuators
                dmtied(tie == tienum) = 0;
            else
                % Grab mean value of all actuators that are tied together
                mean_val = mean(dmtied(tie == tienum));
                % Assign indices in DM mask to mean value
                dmtied(tie == tienum) = mean_val;
            end

        end
    end
    
    
    function isCompliant = checktie(tie)
        % Check whether a tie matrix satisfies the format constraints
        %  - tie should be -1, 0, or the integer range 1-N for some N (with no gaps)
        % 
        % Arguments:
        %  tie: 2D tie matrix to check
        % 
        % Returns:
        %  True if valid, False if not
        Check.two_dim_array(tie)

        tienumset = sort(unique(tie(:).')); % set(tie.ravel())
        tmp = tienumset(tienumset~=0 & tienumset~=-1); % Remove -1 and 0.
        vs = 1:length(tmp);
        if isempty(tmp) && isempty(vs)
            isCompliant = true;
        elseif isequal(tmp, vs)
            isCompliant = true;
        else
            isCompliant = false;
        end
    end


    function isCompliant = checkflat(flatmap, vmin, vmax, tie)
        % Check whether a flatmap matrix satisfies the format constraints
        %  - flatmap should be >= vmin, <= vmax, have all ties at the same voltage,
        %   and have all dead actuators at 0V voltage.
        % 
        % Arguments:
        %  flatmap: 2D flat matrix to check
        %  vmin: min voltage, used in check
        %  vmax: max voltage, used in check
        %  tie: 2D tie matrix to use in check, same size as flatmap
        % 
        % Returns:
        %  True if valid, False if not

        Check.two_dim_array(flatmap)
        Check.real_scalar(vmin)
        Check.real_scalar(vmax)
        Check.two_dim_array(tie)
        if ~isequal(size(tie), size(flatmap))
            error('flatmap and tie must be the same shape')
        end
        if ~ConstrainDM.checktie(tie)
            error('tie matrix contents do not match tie spec')
        end

        isCompliant = true;
        if any(any(flatmap < vmin))
            isCompliant = false;
        end
        if any(any(flatmap > vmax))
            isCompliant = false;
        end
        if any(any(flatmap(tie == -1) ~= 0))
            isCompliant = false;
        end

        tienumset = sort(unique(tie(:).')); % set(tie.ravel())
        tmp = tienumset(tienumset~=0 & tienumset~=-1); % Remove -1 and 0.
        for t = tmp
            if length(unique(flatmap(tie == t))) ~= 1
                % non-identical values implies more than one element
                isCompliant = false;
            end
        end        
    end

    
    function isCompliant = check_valid_neighbors(array, plus_limit, diag_limit, high_limit, low_limit, dmflat)
        % Voltage constraints (caps/neighbor rules) only; does not check ties.  This
        % function is also used for dmsmooth checkout, so not worth trying to force
        % them in.  check_tie_dead() will cover tied and dead actuators.
        % 
        % Arguments:
        %  array: 2D voltage array of interest
        %  plus_limit: the limit on the absolute value difference of elements
        %   along the plus dimension.
        %  diag_limit: the limit on the absolute value difference of elements
        %   along the diagonal dimension.
        % 
        % Keyword Arguments:
        %  high_limit: the maximum value of the arrays, or None.  If None, assumes
        %   there are no upper voltage caps
        %  low_limit: the minimum value of the arrays, or None.  If None, assumes
        %   there are no lower voltage caps
        %  dmflat: 2D voltage array of the same size as 'array' representing the set
        %   of voltages which make the DM surface flat, or None.  These may not
        %   necessarily be uniform.  The requirement is that neighbor rules are
        %   obeyed after the subtraction of a flatmap.  If None, the flatmap is
        %   treated as a uniform array of 0V.
        % 
        % Returns:
        %  True if all elements pass the test, False otherwise

        isCompliant = true;

        % Check input formats
        Check.two_dim_array(array)
        Check.real_positive_scalar(plus_limit)
        Check.real_positive_scalar(diag_limit)

        if ~isempty(high_limit)
            Check.real_scalar(high_limit)
        end
        if ~isempty(low_limit)
            Check.real_scalar(low_limit)
        end

        if ~isempty(dmflat)
            Check.two_dim_array(dmflat)
        else
            dmflat = zeros(size(array));
        end

        % Raise exceptions if out of bounds high/low
        if ~isempty(high_limit) && ~isempty(low_limit)
            lims = ConstrainDM.check_high_low_limit(array, high_limit, low_limit);
        elseif ~isempty(high_limit)
            lims = ConstrainDM.check_high_low_limit(array, high_limit, -Inf);
        elseif ~isempty(low_limit)
            lims = ConstrainDM.check_high_low_limit(array, Inf, low_limit);
        else
            % No voltage bounds supplied, so they can't fail
            lims = true;
        end

        if ~lims
            isCompliant = false;
        end

        fromflat = array - dmflat;

        % Raise exceptions if we break neighbor rules
        for irow = 1:size(fromflat, 1)
            for icol = 1:size(fromflat, 2)
                % index is the (i, j) coordinate,
                % testval is value at that coordinate
                testval = fromflat(irow, icol);

                pn_vals = ConstrainDM.neighbor_values_plus(fromflat, irow, icol);
                for q = pn_vals
                    if any(abs(q - testval) > plus_limit)
                        isCompliant =  false;
                    end
                end

                dn_vals = ConstrainDM.neighbor_values_diag(fromflat, irow, icol);
                for q = dn_vals
                    if any(abs(q - testval) > diag_limit)
                        isCompliant = false;
                    end
                end
            end
        end

    end
    
    function isCompliant = check_tie_dead(volts, tie)
        % Checks that all tied actuator groups are at the same voltage, and that all
        % the dead actuators are at zero volts.
        %
        % Arguments:
        %  volts: 2D array of floating-point voltages to check
        %  tie: 2D array of integers of the same size as volts.  tie should be -1, 0,
        %   or the integer range 1-N for some N (with no gaps).  -1 indicates a dead
        %   actuator; integers > 0 denote tie groups.
        % 
        % Returns:
        %  True if all tied groups are at the same voltage, and all dead actuators
        %   are 0V.  False otherwise.
        Check.two_dim_array(volts)
        Check.two_dim_array(tie)
        if ~isequal(size(volts), size(tie))
            error('ValueError:SizeMismatch','volts and tie must be the same shape')
        end
        
        isCompliant = true;
        
        % Get list of groups
        tienumset = sort(unique(tie(:).')); % row vectorize
        for t = tienumset
            %fprintf('t = %s\n', string(t))
            if t == -1 % dead
                if ~all(volts(tie == t) == 0)
                    isCompliant = false;
                end
                
            elseif t == 0
                % not tied or dead, skip
            else
                tmp = sort(unique(volts(tie == t)));
                %fprintf(' Length of tmp = %d\n', length(tmp))
                if length(tmp) ~= 1
                    isCompliant = false;
                end
            end
        end

    end
    
    function isCompliant = check_high_low_limit(array, high_limit, low_limit)
        % Checks if an array is between two values, inclusive
        % 
        % Inputs:
        %  array: 2D voltage array
        %  high_limit: maximum voltage to compare against
        %  low_limit: minimum voltage to compare against.  Must be lower than
        %   high_limit
        % 
        % Outputs:
        %  Returns True if array is in bounds, otherwise returns False

        % Input Checks
        Check.two_dim_array(array)
        Check.real_scalar(high_limit)
        Check.real_scalar(low_limit)
        if low_limit > high_limit
            error('low_limit must be lower than high_limit')
        end

        if (array >= low_limit) & (array <= high_limit)
            isCompliant = true;
        else
            isCompliant = false;
        end

    end

    
    function isCompliant = check_index_validity(index, nRows, nCols)
        % Check if an index is within the rows and columns.
        %
        % Arguments:
        %  index: a 2-element iterable of integers with the row and column index,
        %   respectively, to be tested.  Maybe be negative, but negative index
        %   elements will always produce a False return.
        %  nRows: total number of rows in array (integer > 0)
        %  nCols: total number of cols in array (integer > 0)
        % 
        % Returns:
        %  returns True if this index would fall in a rows x cols array, False
        %   otherwise

        % Input checks
        if numel(index) ~= 2
            error('index must have 2 elements')
        else
            x = index(1);
            y = index(2);
            Check.scalar_integer(x)
            Check.scalar_integer(y)
        end
        Check.positive_scalar_integer(nRows)
        Check.positive_scalar_integer(nCols)

        isCompliant = false;
        if x >= 1
            if x <= nRows
                if y >= 1
                    if y <= nCols
                        isCompliant = true;
                    end
                end
            end
        end
        
    end

    
    function neighbor_vals = neighbor_values_plus(array, i, j)
        % For a given index in an array, returns the VALUES of the
        % neighboring points according to the 'PLUS' geometry
        % [ ] [x] [ ]
        % [x] [ ] [x]
        % [ ] [x] [ ]
        % array:  a 2d array representing the DM
        % i    :  the x coordinate of the index
        % j    :  the y coordinate of the index
        %
        % Returns:
        %    a list of values of neighboring points
        Check.two_dim_array(array)
        Check.nonnegative_scalar_integer(i)
        Check.nonnegative_scalar_integer(j)

        [nRows, nCols] = size(array);
        if  ~ConstrainDM.check_index_validity([i, j], nRows, nCols)
            error([str(i) + ", " + str(j) " not a valid index for array size!"])
        end
        neighbors = [i-1, j; ...
                     i+1, j; ...
                     i, j-1; ...
                     i, j+1];

        % Remove neighbor locations outside of the array
        for inb = size(neighbors, 1):-1:1
            valid_neighbors(inb) = ConstrainDM.check_index_validity(neighbors(inb, :), nRows, nCols);
        end
        neighbors = neighbors(valid_neighbors, :);
        
        % Retrieve neighbor values
        for inb = size(neighbors, 1):-1:1
            neighbor_vals(inb) = array(neighbors(inb, 1), neighbors(inb, 2));
        end
                
    end
    
    function neighbor_vals = neighbor_values_diag(array, i, j)
        % For a given index in an array, returns the VALUES of the
        % neighboring points according to the 'X' geometry:
        % [x] [ ] [x]
        % [ ] [ ] [ ]
        % [x] [ ] [x]
        % array:  a 2d array representing the DM
        % i    :  the x coordinate of the index
        % j    :  the y coordinate of the index
        %
        % Returns:
        %     a list of values of coordinates of neighboring points
        Check.two_dim_array(array)
        Check.nonnegative_scalar_integer(i)
        Check.nonnegative_scalar_integer(j)

        [nRows, nCols] = size(array);
        if  ~ConstrainDM.check_index_validity([i, j], nRows, nCols)
            error([str(i) + ", " + str(j) " not a valid index for array size!"])
        end
        neighbors = [i-1, j+1;...
                     i+1, j+1;...
                     i+1, j-1;...
                     i-1, j-1];

        % Remove neighbor locations outside of the array
        for inb = size(neighbors, 1):-1:1
            valid_neighbors(inb) = ConstrainDM.check_index_validity(neighbors(inb, :), nRows, nCols);
        end
        neighbors = neighbors(valid_neighbors, :);
        
        % Retrieve neighbor values
        for inb = size(neighbors, 1):-1:1
            neighbor_vals(inb) = array(neighbors(inb, 1), neighbors(inb, 2));
        end

    end
    
     
end

end
