%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the demapping of a given binary vector
%
% Parameters : vect.
% 
% Return : demap : the vector containing the real part and imaginary part
% of the demapped binary vector.
%
% Example :  vect = [0 1 0]
%            demap = demapping(vect)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function demap = demapping(vect)
    if length(vect) == 1    % BPSK
        temp = [1-2*vect];
    elseif length(vect) == 2    % QPSK
        temp = [1-2*vect(1) 1-2*vect(2)];
    else     % AMPM
        if vect == [0 0 0]
            temp = [1 -1];
        elseif vect == [0 0 1]
            temp = [-3 3];
        elseif vect == [0 1 0]
            temp = [1 3];
        elseif vect == [0 1 1]
            temp = [-3 -1];
        elseif vect == [1 0 0]
            temp = [3 -3];
        elseif vect == [1 0 1]
            temp = [-1 1];
        elseif vect == [1 1 0]
            temp = [3 1];
        elseif vect == [1 1 1]
            temp = [-1 -3];
        end
    end
    demap = temp;
end