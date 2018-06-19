function [tfact] = defac(tair)

  if(tair .lt. -10.0) then
     tair = -10.0;
  endif

  tfact = (47.91 / (1.0 + exp(106.06/(tair+18.27))));
