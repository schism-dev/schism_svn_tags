function [idxD, flagSet]=setFind(idxA, idxB)
%[idxD, flagSet]=setFind(idxA, idxB)
% given the set of numbers idxB find their indexes in the ordered set idxA
% returns idxD set of indexes and 
% flagSet - set of flags for each memebr of idxA
%
% example idxA = [1 2 3 4 5 5], idxB=[2 5 6] -> flagSet=[0 1 0 0 1 1], idxD=[2 5 6]
%
% Sergey Frolov, August 2004

%     idxnA               = salt.h.idx.idxNodes;
%     idxnB               = state.dry.idxE;
    flagSet               = uint8(zeros(size(idxA)));
    for i=1:length(idxB)
        flagSet(find(idxA == idxB(i))) = 1;
    end
    idxTmp              = uint32([1:1:length(flagSet)]');
    idxD                = idxTmp(find(flagSet));  
