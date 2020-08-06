function [cnap, errval] = CNAremoveGenericReactionData(cnap,idx,varargin)
%
% CellNetAnalyzer API function 'CNAremoveGenericReactionData'
% ---------------------------------------------
% --> removes a Generic-Reaction-Data field
%
% Usage: err = CNAremoveGenericReactionData(cnap,idx,key);
%
% cnap: CNA project with updated reaction notes / fields
% errval : 0 if a field was deleted, 1 if nothing was deleted
%

if idx > cnap.numr
   errval = 1;
   return;
end

[aos, errval] = deserializeReactionNotes(cnap);
if errval ~= 0
   return;
end

for key = varargin
    aos(idx) = cellfun(@(x) rmfield(x,key),aos(idx),'UniformOutput',false);
end
[cnap,errval] = serializeReactionNotes(cnap,aos);
