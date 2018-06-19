clear;clc


 prompt = 'What is the filename? ';
 fname = input(prompt, 's');
 prompt = 'How many lines woud you like to merge into one line? (Must make sure to be a dividable number to the total line number of orginal data):';
 perline = input(prompt);
 org_data = csvread(fname);
 org_data(1:2,:)
 [x y] = size(org_data);
 new_data = zeros(round(x/perline), y*perline);

 for i=1:ceil(x/perline)
   for j=1:perline
     new_data(i, (j-1)*y+1:j*y) = org_data(2*(i-1)+j,1:y);
   end
 end

 prompt = 'Finished transformation, please specify new file to be saved.'
 newfname = input(prompt, 's');
 csvwrite(newfname,new_data);

 %% Done
