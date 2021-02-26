%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function converts a decimal number to a binary sequence
%
% Parameters : dec : the decimal number to convert.
%              nb : the number of wanted bits for the binary sequence.
% 
% Return : bits : the binary sequence.
%
% Example :  bits = optiDe2Bi(10, 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bits = optiDe2Bi(dec, nb)
    bits = zeros(1,nb);
    for i=1:nb
        if dec - 2^(nb-i) >= 0
            bits(i) = 1;
            dec = dec-2^(nb-i);
        end
    end
end