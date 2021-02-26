%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function procedes the decoding of a message using the Viterbi
% algorithm (soft)
%
% Parameters : message : the encoded message to be decoded.
%              trellis : the trellis of the convolutional code
% 
% Return : decoded : the decoded message.
%
% Example : message = [0 1 0 1 1 1 0 0 1 0]
%            trellis = generateTrellis(1,2,[5 7])
%            -> Rate 1/2 convolutional code with G = (1+D^2,1+D+D^2)
%            decoded = viterbiDecodingSoft(message, trellis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function decoded = viterbiDecodingSoft(message, trellis)

    k = log2(trellis.numInputSymbols);  % Number of input bits
    n = log2(trellis.numOutputSymbols); % Number of output bits
    paths = -1*ones(trellis.numStates,length(message)/2,3);	% Contains the total metrics of each state at each transition, the transmitted symbol and its position
    
    % Initialization (1st transition of the trellis)
    initial_state = 0;  % Start from the 0-state
    cpt = 1;    % Transition n°1 of the trellis
    for ii=1:trellis.numInputSymbols
       paths(trellis.nextStates(initial_state+1,ii)+1,cpt,1) = distanceEucl(message(1:2),demapping(optiDe2Bi(trellis.outputs(initial_state+1,ii),n)));  % Branch metric stored in the first dim of paths
       paths(trellis.nextStates(initial_state+1,ii)+1,cpt,2) = ii-1;    % Transmitted symbol stored in the second dim of paths
    end
    
    % Generalization
    for jj=3:2:length(message)
        cpt = cpt+1;    % Indicates the current transition
        for kk=1:trellis.numStates
           if (paths(kk,cpt-1,1) ~= -1)
               previousMetric = paths(kk,cpt-1,1);  % Total metric of the last transition
               for ii=1:trellis.numInputSymbols
                   branchMetric = distanceEucl(message(jj:jj+1),demapping(optiDe2Bi(trellis.outputs(kk,ii),n)));  % Current branch metric
                   if (paths(trellis.nextStates(kk,ii)+1,cpt,1) == -1)  % If we haven't stored any metric yet
                       paths(trellis.nextStates(kk,ii)+1,cpt,1) = previousMetric + branchMetric;    % We store the total metric
                       paths(trellis.nextStates(kk,ii)+1,cpt,2) = ii-1; % We store the transmitted symbol associated to this metric
                       paths(trellis.nextStates(kk,ii)+1,cpt,3) = kk;   % We store the position of the transmitted symbol
                   else
                       if (previousMetric + branchMetric < paths(trellis.nextStates(kk,ii)+1,cpt,1))    % We replace if we find a lower total metric
                           paths(trellis.nextStates(kk,ii)+1,cpt,1) = previousMetric + branchMetric;
                           paths(trellis.nextStates(kk,ii)+1,cpt,2) = ii-1;
                           paths(trellis.nextStates(kk,ii)+1,cpt,3) = kk;
                       end
                   end
               end
           end
        end
    end
    
    for i=1:trellis.numStates
       paths(i,:,:) = fliplr(paths(i,:,:)); % We flip the matrix for easier processing (traceback from left to right)
    end
    
    decoded = [optiDe2Bi(paths(find(min(paths(:,1,1)) == paths(:,1,1)),1,2),k)];    % Find the transmitted symbol associated to the minimum total metric obtained, i.e for the last transition (located in the first column of paths)
    traceback_pos = paths(find(min(paths(:,1,1)) == paths(:,1,1)),1,3); % Position of the next symbol
    for j=2:length(message)/2
        decoded = [decoded fliplr(optiDe2Bi(paths(traceback_pos,j,2),k))];  % We append the transmitted symbols for each transition of the good path
        traceback_pos = paths(traceback_pos,j,3);   % Update of the position of the next symbol
    end
    
    decoded = fliplr(decoded);  % We flip again in order to restore the good order
    
end