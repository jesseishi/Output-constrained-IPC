function FAST_indexing = GetOutDataIndices
% I haven't found a perfect way for indexing yet, but I like this the most
% so far.
% If OutList doesn't exist, do a 1 sec simulation (which will fail) but
% will generate OutList, which you can then load here.
load('OutList.mat', 'OutList')
FAST_indexing.n_vars = length(OutList);
FAST_indexing.Azimuth = find(matches(OutList, 'Azimuth'));
FAST_indexing.GenSpeed = find(matches(OutList, 'GenSpeed'));
temp = find(matches(OutList, 'RootMyc1'));
FAST_indexing.RootMyci = temp:temp+2;
end