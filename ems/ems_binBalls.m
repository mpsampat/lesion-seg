function list = ems_binBalls(nrOfBins,nrOfBalls)
% FORMAT list = ems_binBalls(nrOfBins,nrOfBalls)
%
% Returns all possible combinations for dividing nrOfBalls identical
% balls over nrOfBins bins. 
%
% ------------------------------------------------------------------------
% ems_binBalls.m      Koen Van Leemput - August 17, 2001


if (nrOfBins<1)
  error('nrOfBins must be greater than zero');
end


if (nrOfBins==1)
  list = nrOfBalls;
else
  list = []; 
  for i=0:nrOfBalls
    tmplist = ems_binBalls(nrOfBins-1, nrOfBalls-i);
    for j=1:size(tmplist,1)
      list = [list; [i tmplist(j,:)]];
    end
  end
end




