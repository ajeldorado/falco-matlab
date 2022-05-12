% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute which actuators are frozen or co-moving.
%

classdef ActLimits
   
methods (Static)
    
    function  [freeze, link] = maplimits(dmv, dmObject, varargin)
    %     """
    %     Create lists of actuators to freeze and link for a single DM for HOWFSC.
    % 
    %     This includes an input matrix (tiemap) which indicates which actuators have
    %     fixed electrical constraints (dead or driven together)
    % 
    %     This map incorporates 5 effects:
    %      - dead actuators are frozen
    %      - actuators at or below low-voltage constraints are frozen
    %      - actuators at or above high-voltage constraints are frozen
    %      - electrically-connected actuators groups are linked together to move
    %        together
    %      - actuators at the boundaries of neighbor rule violations (both lateral
    %       and diagonal) are linked togther to move together
    %     Since HOWFSC computes changes in DM settings, constraining them to move
    %     together is equivalent to forcing their delta-DM settings to be the same.
    %     It doesn't force the underlying absolute voltages to be the same.
    % 
    %     Does not fix neighbor rule violations or cap violations! It only identifies
    %      ones that do violate; the handling is elsewhere.
    % 
    %     Regions can be marked as candidates for either freezing or linking, or
    %     neither. They will never be a member of both, and if an actuator would be
    %     in both, it is frozen.  Indices will be stored in ndarrays.
    % 
    %     For each of x, y, diag, other diag:
    %      find all the violations, pick the first one, and check it for adjacent
    %      violations until you run out nearby acts.  Go to the next one.  Repeat
    %      until x is done, then y, then diag, then other diag.
    % 
    %     Based on original HCIT implementation by Brian Kern.
    % 
    %     Arguments:
    %      dmv: array with current absolute voltages
    %      dmObject: DM object containing voltage information for a specific DM
    % 
    %     Keyword Arguments:
    %      tiemap: a 2D array of integers, of the same size as dmv, which can take
    %       on the values 0, -1, or consecutive integers 1 -> N, or None.  If not
    %       None, these values will encode dead and electrically-connected actuators
    %       (dead = -1, connected = 1, 2, etc., neither = 0).  If None, this function
    %       will behave as though there are no electrical constraints.
    % 
    %       The tiemap encodes actuators that not only move together, but also must
    %       have the same absolute voltage.  Generally it represents physical
    %       hardware constraints.
    % 
    %     Returns:
    %      dictionary with two keys: ['freeze', 'link'], which are lists of
    %       indices of actuators which are low/high (and need to be frozen) and
    %       groups that are outside the neighbor-rule limit (and need to be link)
    % 
    %     """
        ngval = 0; % value for "no group"
        deadval = -1; % special code for dead actuators
        invalid = -2; % value not otherwise permitted, must be < ngval

        % Check inputs
        Check.two_dim_array(dmv)
        
        % Optional inputs
        if length(varargin) == 1
            tiemap = varargin{1};
        elseif length(varargin) > 1
            error('Too many inputs.')
        else
            tiemap = ngval*ones(size(dmv));
        end

        Check.two_dim_array(tiemap)
        if ~ConstrainDM.checktie(tiemap)  % enforce tie formatting
            error('tie must have values 0, -1, or consecutive integers 1 -> N')
        end
        if ~all(size(dmv) == size(tiemap))
            error('dmv shape must match tiemap')
        end

        ddm = dmv - dmObject.facesheetFlatmap;

        margin = 2*dmObject.marginNbrRule; % use 2LSB to match dmsmooth implementation

        % Call it a match if it's within margin to avoid rounding issues
        vlo = dmObject.Vmin + margin;
        vhi = dmObject.Vmax - margin;
        vnb = dmObject.dVnbrLat - margin;
        vco = dmObject.dVnbrDiag - margin;

        % freeze if *at* threshold as well as beyond, so we can fix values to
        % a threshold elsewhere and they stick.
        mlohi = (dmv <= vlo) | (dmv >= vhi);
    %     mlohi = np.logical_or(dmv <= vlo, dmv >= vhi)

        % freeze dead actuators
        mlohi = mlohi | (tiemap == deadval);
    %     mlohi = np.logical_or(mlohi, tiemap == deadval)

        delx = ddm(:, 2:end) - ddm(:, 1:end-1); % x neighbors
        dely = ddm(2:end, :) - ddm(1:end-1, :); % y neighbors
        delxy = ddm(2:end, 2:end) - ddm(1:end-1, 1:end-1); % +x+y neighbors
        delyx = ddm(2:end, 1:end-1) - ddm(1:end-1, 2:end); % -x+y neighbors
    %     delx = ddm[:, 1:] - ddm[:, :-1] % x neighbors
    %     dely = ddm[1:, :] - ddm[:-1, :] % y neighbors
    %     delxy = ddm[1:, 1:] - ddm[:-1, :-1] % +x+y neighbors
    %     delyx = ddm[1:, :-1] - ddm[:-1, 1:] % -x+y neighbors

        mx = abs(delx) > vnb;
        my = abs(dely) > vnb;
        mxy = abs(delxy) > vco;
        myx = abs(delyx) > vco;
    %     mx, my = [np.abs(d) > vnb for d in [delx, dely]]
    %     mxy, myx = [np.abs(d) > vco for d in [delxy, delyx]]
        ms = {mx, my, mxy, myx}; % pack for passing to _growgrp compactly
        
%         figure(11); imagesc(mx); axis xy equal tight; colorbar;
%         figure(12); imagesc(my); axis xy equal tight; colorbar;
%         figure(13); imagesc(mxy); axis xy equal tight; colorbar;
%         figure(14); imagesc(myx); axis xy equal tight; colorbar;
%         drawnow;
        

        % igrp is a map of group members. 0 is no group, 1+ by int are individual
        % groups, -1 is the special group of dead actuators.  -2 is used by
        % igrploop only to get through first call to while() since it will always
        % fail
        igrp = ngval * ones(size(tiemap));
%         figure(11); imagesc(igrp); axis xy equal tight; colorbar; drawnow;

        ig = ngval;
        igrploop = invalid * ones(size(igrp));
        while any(any(igrp ~= igrploop))
            igrploop = igrp;
            mnog = (igrp == ngval);
            mnew = mx & (mnog(:, 2:end) & mnog(:, 1:end-1));
%             figure(20); imagesc(mnew); axis xy equal tight; colorbar; title('mnew'); drawnow; pause(0.2); 

    %         mnew = np.logical_and(mx, np.logical_and(mnog[:, 1:], mnog[:, :-1]))
            % new groups that contain x neighbors, as (row, col) tuples
            % from the difference matrix mx
            wnew = find(mnew);
    %         wnew = np.transpose(np.nonzero(mnew))
            if numel(wnew) > 0
                firstw = wnew(1);
                ig = ig + 1;
                [row, col] = ind2sub(size(mnew), firstw);
                igrp(row, col) = ig;
                igrp(row, col+1) = ig;
%                 figure(31); imagesc(igrp); axis xy equal tight; colorbar; title('Before'); drawnow;

    %             igrp[firstw[0], firstw[1]] = ig
    %             igrp[firstw[0], firstw[1]+1] = ig
                igrp = ActLimits.grow_group(ig, igrp, ms);
%                 figure(32); imagesc(igrp); axis xy equal tight; colorbar; title('After'); drawnow;
%                 pause(1);

            end
        end
%         figure(12); imagesc(igrp); axis xy equal tight; colorbar; drawnow;


        igrploop = invalid * ones(size(igrp));
        while any(any(igrp ~= igrploop))
            igrploop = igrp;
            mnog = igrp == ngval;
            mnew = (my & (mnog(2:end, :) & mnog(1:end-1, :)));
    %         mnew = np.logical_and(my, np.logical_and(mnog[1:, :], mnog[:-1, :]))
            % new groups that contain y neighbors, as (row, col) tuples
            % from the difference matrix my
            wnew = find(mnew);
    %         wnew = np.transpose(np.nonzero(mnew))
            if numel(wnew) > 0
                firstw = wnew(1);
                ig = ig + 1;
                [row, col] = ind2sub(size(mnew), firstw);
                igrp(row, col) = ig;
                igrp(row+1, col) = ig;
    %             igrp[firstw[0], firstw[1]] = ig
    %             igrp[firstw[0]+1, firstw[1]] = ig
                igrp = ActLimits.grow_group(ig, igrp, ms);
            end
        end

        igrploop = invalid * ones(size(igrp));
        while any(any(igrp ~= igrploop))
            igrploop = igrp;
            mnog = igrp == ngval;
            mnew = (mxy & (mnog(2:end, 2:end) & mnog(1:end-1, 1:end-1)));
    %         mnew = np.logical_and(mxy, np.logical_and(mnog[1:, 1:], mnog[:-1, :-1]))
            % new groups that contain +x+y neighbors, as (row, col) tuples
            % from the difference matrix mxy
            wnew = find(mnew);
    %         wnew = np.transpose(np.nonzero(mnew))
            if numel(wnew) > 0
                firstw = wnew(1);
                ig = ig + 1;
                [row, col] = ind2sub(size(mnew), firstw);
                igrp(row, col) = ig;
                igrp(row+1, col+1) = ig;
    %             igrp[firstw[0], firstw[1]] = ig
    %             igrp[firstw[0]+1, firstw[1]+1] = ig
                igrp = ActLimits.grow_group(ig, igrp, ms);
            end
        end

        igrploop = invalid * ones(size(igrp));
        while any(any(igrp ~= igrploop))
            igrploop = igrp;
            mnog = igrp == ngval;
            mnew = (myx & (mnog(2:end, 1:end-1) & mnog(1:end-1, 2:end)));
    %         mnew = np.logical_and(myx, np.logical_and(mnog[1:, :-1], mnog[:-1, 1:]))
            % new groups that contain -x+y neighbors, as (row, col) tuples
            % from the difference matrix myx
            wnew = find(mnew);
    %         wnew = np.transpose(np.nonzero(mnew))
            if numel(wnew) > 0
                firstw = wnew(1);
                ig = ig + 1;
                [row, col] = ind2sub(size(mnew), firstw);
                igrp(row, col+1) = ig;
                igrp(row+1, col) = ig;
    %             igrp[firstw[0], firstw[1]+1] = ig
    %             igrp[firstw[0]+1, firstw[1]] = ig
                igrp = ActLimits.grow_group(ig, igrp, ms);
            end
        end

        % Combine tiemap with neighbor-rule-only map at the very end
        filt = tiemap;
        filt(filt == deadval) = 0; % don't try to treat dead as linked together here
        fgrps = setdiff(filt, ngval);
    %     fgrps = set(filt.ravel()) - {ngval}

        if ~isempty(fgrps)
            for f = fgrps(:).'
                inds = (filt == f);
                igrp_num = setdiff(igrp(inds), ngval);

                % Pick a new number not used in igrp, and set 1) all igrp indices that
                % match the tiemap group (inds) to that number and 2) all igrp members
                % that share a number with a element in igrp[inds], to spread out and
                % catch nonlocal linkages
                nextval = max(igrp(:)) + 1;
                for n = igrp_num(:).'
                    igrp(igrp == n) = nextval;
                end
                igrp(inds) = nextval;
            end
        end

        % Separate low and high for freezing; rest to be linked
        igrplohi = unique(igrp(mlohi));
        for ii = igrplohi(:).' % if a group contains any low or high limits,
                          % freeze whole group
            if ii > ngval
                mlohi(igrp == ii) = true;
                igrp(igrp == ii) = ngval; % If it's frozen, don't link it together
            end
        end

        freeze = find(mlohi ~= 0);
    %     freeze = np.flatnonzero(mlohi)
        link = {}; % cell array of vectors
        igrp_num = setdiff(igrp, ngval);
        counter = 0;
        if ~isempty(igrp_num)
            for ii = igrp_num(:).'
                matches = (igrp == ii);
                igrplist = find(matches ~= 0);
    %             igrplist = np.flatnonzero(igrp == ii)
                if numel(igrplist) ~= 0 % strip empty arrays left after freeze
                    counter = counter + 1;
                    link{counter} = igrplist;
    %                 link.append(igrplist)
                end
            end
        end
    %     return {freeze, 'link':link}
    end

    function igrp = grow_group(ig, igrp, ms)

    %     Given a seed group of neighbor-rule-violating actuators, expand that group
    %     until it includes all NR violations linked to that seed group.
    % 
    %     This function modifies the igrp array in place! As such it returns nothing.
    % 
    %     Arguments:
    %      ig: group index, will be an integer >= 0
    %      igrp: 2D array of integers of the same size as the DM.  This matrix
    %       assigns a group number to each actuator (integer >= 0) or -1 if that
    %       actuator is not part of a group.  Actuators in a group all move together.
    %      ms: a cell array with four arrays of neighbor-rule violators -> mx (horizontal)
    %       my (vertical), mxy (diagonal), myx (other diagonal)

        Check.nonnegative_scalar_integer(ig)
        Check.two_dim_array(igrp)
        if ~iscell(ms)
            error('ms must be a cell array with length 4.')
        elseif length(ms) ~= 4
            error('ms must be a cell array with length 4.')
        end

        mx = ms{1};
        my = ms{2};
        mxy = ms{3};
        myx = ms{4};

        % Loop through and check all left/right/up/down/diagonal neighbors to
        % existing group to see if they are neighbor-rule violators with group
        % members.  If they are, expand the group around them and recheck to see
        % if the expansion picked up new violators.  If not, call it done.
        igrplast = (zeros(size(igrp)));
%         igrplast = logical(zeros(size(igrp)));
%         figure(61); imagesc(igrp); axis xy equal tight; colorbar; title('igrp'); drawnow;
%         figure(62); imagesc(igrplast); axis xy equal tight; colorbar; title('igrplast'); drawnow;
        while any(any(igrplast ~= igrp))

            igrplast = igrp;

            % horizontal
            ig0 = (igrp == ig);
            mnew = (mx & (ig0(:, 2:end) | ig0(:, 1:end-1)));
            wnew = find(mnew);
            for w = wnew(:).'
                [row, col] = ind2sub(size(mnew), w);
                igrp(row, col) = ig;
                igrp(row, col+1) = ig;
            end
%             figure(41); imagesc(igrp); axis xy equal tight; colorbar; title('post-horizontal'); drawnow;

    %         mnew = np.logical_and(mx, np.logical_or(ig0[:, 1:], ig0[:, :-1]))
    %         wnew = np.transpose(np.nonzero(mnew))
    %         for w = wnew
    %             igrp[w[0], w[1]] = ig
    %             igrp[w[0], w[1]+1] = ig
    %         end

            % vertical
            ig0 = (igrp == ig);
            mnew = (my & (ig0(2:end, :) | ig0(1:end-1, :)));
            wnew = find(mnew);
            for w = wnew(:).'
                [row, col] = ind2sub(size(mnew), w);
                igrp(row, col) = ig;
                igrp(row+1, col) = ig;
            end
%             figure(42); imagesc(igrp); axis xy equal tight; colorbar; title('post-vertical'); drawnow;
    %         wnew = np.transpose(np.nonzero(mnew))
    %         for w in wnew:
    %             igrp[w[0], w[1]] = ig
    %             igrp[w[0]+1, w[1]] = ig
    %         end

            % one diagonal
            ig0 = (igrp == ig);
            mnew = (mxy & (ig0(2:end, 2:end) | ig0(1:end-1, 1:end-1)));
            wnew = find(mnew);
%             figure(100); imagesc(mnew); axis xy equal tight; colorbar; title('mnew internal'); drawnow; %pause(5);
            for w = wnew(:).'
                [row, col] = ind2sub(size(mnew), w);
                igrp(row, col) = ig;
                igrp(row+1, col+1) = ig;
            end       
%             figure(43); imagesc(igrp); axis xy equal tight; colorbar; title('post-diagonal'); drawnow;
    %         wnew = np.transpose(np.nonzero(mnew))
    %         for w in wnew:
    %             igrp[w[0], w[1]] = ig
    %             igrp[w[0]+1, w[1]+1] = ig
    %         end

            % other diagonal
            ig0 = (igrp == ig);
            mnew = (myx & (ig0(2:end, 1:end-1) | ig0(1:end-1, 2:end)));
            wnew = find(mnew);
            for w = wnew(:).'
                [row, col] = ind2sub(size(mnew), w);
                igrp(row, col+1) = ig;
                igrp(row+1, col) = ig;
            end   
%             figure(44); imagesc(igrp); axis xy equal tight; colorbar; title('post-other-diagonal'); drawnow;
%             pause(3);
    %         wnew = np.transpose(np.nonzero(mnew))
    %         for w = wnew
    %             igrp[w[0], w[1]+1] = ig
    %             igrp[w[0]+1, w[1]] = ig
    %             pass
    %         end

        end

    end

end

end

