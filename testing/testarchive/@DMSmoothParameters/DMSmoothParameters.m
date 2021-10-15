% DMSmoothParameters Function

function self = DMSmoothParameters()

    VMIN = 0;
    self.vmin = VMIN;
    self.vmax = 100.;
    self.vquant = 100./(2^16 - 1.);
    self.vneighbor = 30.;
    self.vcorner = 30.;
    self.nact = 48;

    self.vmid = (self.vmin + self.vmax)/2.;

    % Uniform
    self.dm0 = self.vmid * ones(self.nact, self.nact);
%     figure; imagesc(self.dm0); axis xy equal tight; colorbar; title('Uniform');

    % Checkerboard
    self.dmcheck = self.vmin * ones(self.nact, self.nact);
    self.dmcheck(1:2:end, 1:2:end) = self.vmax;
    self.dmcheck(2:2:end, 2:2:end) = self.vmax;
%     figure; imagesc(self.dmcheck); axis xy equal tight; colorbar; title('checkerboard');

    % Up-down checkerboard
    self.dmud = self.vmid*ones(self.nact, self.nact);
    self.dmud(1:2:end, 1:2:end) = self.vmin;
    self.dmud(2:2:end, 2:2:end) = self.vmax;
%     figure; imagesc(self.dmud); axis xy equal tight; colorbar; title('up-down checkerboard');

    % Out-of-range checkerboard
    self.dmoor = (self.vmin - 1.) * ones(self.nact, self.nact);
    self.dmoor(1:2:end, 1:2:end) = self.vmax + 1.;
    self.dmoor(2:2:end, 2:2:end) = self.vmax + 1.;
%     figure; imagesc(self.dmoor); axis xy equal tight; colorbar; title('OOR checkerboard');

    % Horizontal stripes
    self.dmhstripes = self.vmin*ones(self.nact, self.nact);
    self.dmhstripes(1:2:end) = self.vmax;
%     figure; imagesc(self.dmhstripes); axis xy equal tight; colorbar; title('horizontal stripes');

    % Vertical stripes
    self.dmvstripes = self.vmin*ones(self.nact, self.nact);
    self.dmvstripes(:, 1:2:end) = self.vmax;
%     figure; imagesc(self.dmvstripes); axis xy equal tight; colorbar; title('vertical stripes');

    % Uniform and non-uniform flatmaps (which don't violate NR themselves)
    self.uniform_flat = self.vmid*ones(self.nact, self.nact);
    self.nonuniform_flat = self.vmid*ones(self.nact, self.nact) + min([self.vneighbor, self.vcorner])/10.*eye(self.nact);
%     figure; imagesc(self.uniform_flat); axis xy equal tight; colorbar; title('uniform flatmap');
%     figure; imagesc(self.nonuniform_flat); axis xy equal tight; colorbar; title('nonuniform flatmap');

end
