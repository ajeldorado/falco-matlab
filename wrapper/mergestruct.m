function into = mergestruct(into, from)
   %MERGESTRUCT merge all the fields of scalar structure from into scalar structure into
   validateattributes(from, {'struct'}, {'scalar'});
   validateattributes(into, {'struct'}, {'scalar'});
   fns = fieldnames(from);
   for fn = fns.'
      if isstruct(from.(fn{1})) && isfield(into, fn{1})
           %nested structure where the field already exist, merge again
           into.(fn{1}) = mergestruct(into.(fn{1}), from.(fn{1}));
      else
          %non structure field, or nested structure field that does not already exist, simply copy
          into.(fn{1}) = from.(fn{1});
      end
   end
end
