clear;clc

fname = "S04AK-090-001.csv";
socprof = csvread(fname, 1, 1);
filled_prof = zeros(1,100);
filled_prof = socprof(1:100);
for i=1:100
  if(filled_prof(i)<=1e-5)
    filled_prof(i) = filled_prof(i-1);
  end
end

% layer thickness in cm (ceiled to integer)
tk = [2, 3, 5, 8, 13, 21, 34, 56];

% transfer kgC/m3 to kgC/m2
c_profile = zeros(8,1);
st = 1;
for i = 1:7
   c_profile(i) = sum(filled_prof(st:st+tk(i)-1)) * 0.01;
   st = st+tk;
end
c_profile(8) =  filled_prof(100) * tk(8) * 0.01;
c_profile
