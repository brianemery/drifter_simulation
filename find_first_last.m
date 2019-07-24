function [iStart,iEnd,row,cStart,cEnd] = find_first_last(x)
% FIND FIRST LAST - for each row find start and end matrix addresses of not-NaN elements 
% [iStart,iEnd,row,cStart,cEnd] = find_first_last(x)
% 
% OUTPUT
% iStart    - matrix address start of non-nan
% iEnd      - matrix address of the end of non-nan
% row       - row index for cStart, cEnd
% cStart    - column index of start of nan -nan
% cEnd      - column index of end of non-nan
%
% Finds the start and end Matrix addresses of non-NaN data in a matrix. The outputs
% will have the same number of rows as the input.
% 
% Typically used for matricies of drifter/trajectory data. 
%
% To convert the matrix addresses to row,column indicies:
% [rStart,cStart] = ind2sub(size(x),iStart);
%
% NOTE 
% these may need to be recomputed after using subsref_struct, which will
% not apply indexing changes to the vectors. (possible limitation of this
% kind of data structure - maybe cell arrays of vectors would be better?)
%
% EXAMPLE
% [DRFT.iStart,DRFT.iEnd,DRFT.row,DRFT.cStart,DRFT.cEnd] = ...
%                                           find_first_last(DRFT.Lon);

% Brian Emery 7 Jan 2010

if strcmp(x,'--t')
    test_case, return
end

[iStart,iEnd,cStart,cEnd] = deal(NaN(size(x,1),1));

for i = 1:size(x,1)
    n = find(~isnan(x(i,:)));
    if ~isempty(n)
        
        iStart(i) = sub2ind(size(x),i,min(n));
        iEnd(i) = sub2ind(size(x),i,max(n));
        
        cStart(i) = min(n);
        cEnd(i) = max(n);
        
    end
end

row = (1:size(x,1))';

end

%% -----------------------------------------------------------------------
function test_case

% create a matrix with a bunch of NaN's in it
x=NaN(4,7);
x(1,3:end) = 1;
x(2,2:5) = 1;
x(3,5)=1;
x(4,:)=1;

% create a matrix of addresses
y = reshape(1:28,4,7);
 
[iStart,iEnd,row,cStart,cEnd] = find_first_last(x);

keyboard


iStart
iEnd

disp('The Matrix'), x
disp('Index '), y

disp('row, col indecies of iStart:')
[rStart,cStart] = ind2sub(size(x),iStart)

end
