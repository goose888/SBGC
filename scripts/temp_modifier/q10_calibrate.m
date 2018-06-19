clear;clc

%% % Soil profile Samples
%% % S04AK-090-001 (Orthel)
%% csoil_s04 = [153.2576995, 34.743, 22.143, 18.544, 42.6, 26.656];  % kgC/m3
%% depsoil_s04 = [0.06, 0.07, 0.12, 0.26, 0.25, 0.44];    % m
%% [topsoil_s04, botsoil_s04, nodesoil_s04] = getsoilinfo(csoil_s04, depsoil_s04);
%% 
%% % S07AK001 (Turbel)
%% csoil_s07 = [104.857, 94.851, 49.28, 15.105, 8.968, 7.315653163];    % kgC/m3
%% depsoil_s07 = [0.14, 0.19, 0.19, 0.16, 0.44, 0.40];
%% [topsoil_s07, botsoil_s07, nodesoil_s07] = getsoilinfo(csoil_s07, depsoil_s07);
%% 
%% % S09AK290008 (Histel)
%% csoil_s09 = [172.8862411, 170.7165193, 172.3230764, 172.9122776, 128.089266, 34.216, 9.6];  % kgC/m3
%% depsoil_s09 = [0.18, 0.12, 0.27, 0.37, 0.32, 0.30, 0.44];
%% [topsoil_s09, botsoil_s09, nodesoil_s09] = getsoilinfo(csoil_s09, depsoil_s09);

% Call Q10 function here
q10const = 1.5;
tsoil = 17+25;
q10len = length(q10const);
tlen = length(tsoil);
tfact = zeros(q10len, tlen);
for i = 1:q10len
  for j = 1:tlen
    tfact(i,j) = q10(q10const(i), tsoil(j));
  end
end
tfact

