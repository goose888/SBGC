 clear;clc

 transfer_id_2_latlon = false
 if(transfer_id_2_latlon)
   prompt = 'What is the lonid? ';
   lonid = input(prompt);
   prompt = 'What is the latid? ';
   latid = input(prompt);
   if(lonid > 360)
     lon = (lonid-360-0.5)/2 - 180;
   else
     lon = (lonid-0.5)/2;
   end
   lat = (latid-0.5)/2 - 90;
   % Display
   lon
   lat
 else
   prompt = 'What is the lon (transfer to central lon first)? ';
   lon = input(prompt);
   prompt = 'What is the lat (transfer to central lat first)? ';
   lat = input(prompt);
   if(lon < 0)
     lonid = 720 + lon*2 + 0.5;
   else
     lonid = lon*2 + 0.5;
   end
   latid = lat*2 + 0.5 + 180;
   % Display
   lonid
   latid
 end
