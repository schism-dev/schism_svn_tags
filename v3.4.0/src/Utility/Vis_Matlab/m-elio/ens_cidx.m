function idxOut=er_cidx(idxIn)
% idxOut=er_cidx(idxIn)
% helper function to compute idx from 
% string form '[1:456]' or array form [1 ... 456] to strictly array form [1 ... 456]
% idxIn  -- string or num array defineinig an index
% idxOut -- num array
%
% Sergey Frolov Nov, 2004

if ischar(idxIn)
  idxOut = eval(idxIn);
else
  idxOut = idxIn;
end
