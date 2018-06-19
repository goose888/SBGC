function [topsoil, botsoil, nodesoil] = getsoilinfo(csoil, depsoil)
numlayer = length(csoil);
for i = 1:numlayer
   if(i == 1)
     topsoil(i) = 0.0;
     botsoil(i) = depsoil(i);
   else
     topsoil(i) = botsoil(i-1);
     botsoil(i) = botsoil(i-1) + depsoil(i);
   end
end
nodesoil = topsoil + depsoil/2;


