%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the euclidean distance of 2 vectors
%
% Parameters : vect1.
%              vect2.
% 
% Return : d : the euclidean distance.
%
% Example :  vect1 = [0 1 0]
%            vect2 = [1 1 1]
%            d = distanceEucl(vect1, vect2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = distanceEucl(vect1, vect2)
    d = sum((vect1-vect2).^2).^(1/2);
end