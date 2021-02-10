%   Copyright 2019 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology


function prop_free_threads
% Shut down processes created by prop_run_multi
% No inputs or outputs.
%   History:
%   Created 31 Dec 2019 by John Krist
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

  delete( prop_pool );
  n_prop_pool = 0;

end

