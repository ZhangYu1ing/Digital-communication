%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the minimum distance of 2 vectors
%
% Parameters : vect1.
%              vect2.
% 
% Return : dmin : the minimum distance.
%
% Example :  vect1 = [0 1 0]
%            vect2 = [1 1 1]
%            dmin = distanceMin(vect1, vect2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dmin = distanceMin(vect1, vect2)
    vect_sum = [];
    for ii=1:length(vect1)
        vect_sum = [vect_sum mod(vect1(ii)+vect2(ii),2)];
    end
    dmin = numel(find(vect_sum));
end