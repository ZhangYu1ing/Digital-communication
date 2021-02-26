%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the trellis of a convolutional code
% given its generator polynomials
%
% Parameters : k : number of input bits.
%              n : number of output bits.
%              poly : vector of the generator polynomials in octal.
% 
% Return : structure trellis :  numInputSymbols : number of symbols at the
%                               encoder input
%
%                               numOutputSymbols : number of symbols at the
%                               encoder output.
%                               
%                               numStates : number of possible states
%                               in the encoder.
%
%                               nextStates : matrix of all the possible
%                               states.
%                               nextStates(initial state, input symbol) =
%                               final state.
%
%                               outputs : matrix containing the codewords
%                               (in octal).
%                               outputs(initial state, input symbol) =
%                               codeword.
%
% Examples : trellis = generateTrellis(1,2,[5 7])
%           -> Rate 1/2 convolutional code with G = (1+D^2,1+D+D^2)
%           trellis = generateTrellis(1,2,[27 26])
%           -> Rate 1/2 convolutional code with G = (1+D^2+D^3+D^4,1+D^2+D^3)
%           trellis = generateTrellis(1,2,[23 33])
%           -> Rate 1/2 convolutional code with G = (1+D^3+D^4,1+D+D^3+D^4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trellis = generateTrellis(k, n, poly)
    
    CL = [];
    poly_binary = [];   % a line of poly_binary = an element of poly in binary form
    for i=1:length(poly)
       poly_binary = [poly_binary ; fliplr(de2bi(oct2dec(poly(i))))];
       CL = [CL length(poly_binary)];
    end
    constraintLength = max(CL);     % The constraint length of the encoder
    
    trellis.numInputSymbols = 2^k;
    trellis.numOutputSymbols = 2^n;
    trellis.numStates  = 2^(constraintLength-1);

    trellis.nextStates = zeros(trellis.numStates,trellis.numInputSymbols);
    trellis.outputs = zeros(trellis.numStates,trellis.numInputSymbols);
    for ii=1:trellis.numStates
        currentState = optiDe2Bi(ii-1,constraintLength-1);  % all the different states in binary form
        for jj=1:trellis.numInputSymbols
            trellis.nextStates(ii,jj) = bi2de([jj-1 currentState(1:end-1)],'left-msb'); % New state = [input_bit currentState(1:end-1)]
            for kk=1:length(poly)   % Here, number of generator poly = number of output bits
                output(kk) = mod(sum(poly_binary(kk,1:constraintLength).*[jj-1 optiDe2Bi(ii-1,constraintLength-1)]),2); % Output = [generator poly] .* [input_bit state]
                trellis.outputs(ii,jj) = dec2oct(bi2de(output,'left-msb'));     % Result in octal
            end
        end
    end
    
end