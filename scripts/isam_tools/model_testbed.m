 clear;clc

 % Put flags for test sub model component
 Get_node_depth = false;
 Get_layer_thickness = false;
 Get_interface_depth = false;
 Temp_mod = true;
 Moist_mod = false;
 Depth_mod = false;
 Turnover = false;
 CEC_Rh = false;
 Therm_conduct = false;
 Percolation = false;
 Cryoturbation = true;

 % Put your code here to test
 if(Get_node_depth)
   prompt = 'Please specify the total number of layers? ';
   nlevgrnd = input(prompt);
   z = zeros(nlevgrnd,1);
   for j = 1:nlevgrnd
     z(j) = 0.025*(exp(0.5*(j-0.5))-1.);  % node depths
   end
   z
 end

 if(Get_node_depth && Get_layer_thickness)
   dz = zeros(nlevgrnd,1);
   dz(1) = 0.5*(z(1)+z(2));
   for j = 2:nlevgrnd-1
     dz(j)= 0.5*(z(j+1)-z(j-1));          % thickness b/n two interfaces
   end
   dz(nlevgrnd)= z(nlevgrnd)-z(nlevgrnd-1);
   dz
 end

 if(Get_node_depth && Get_layer_thickness && Get_interface_depth)
   zsoih = zeros(nlevgrnd,1);
   for j = 1:nlevgrnd-1
      zsoih(j)= 0.5*(z(j)+z(j+1));        % interface depths
   end
   zsoih(nlevgrnd) = z(nlevgrnd) + 0.5*dz(nlevgrnd);
   zsoih
 end

 if(Temp_mod)
   tair = -10:0.2:40;
   tsoil = -10:0.2:40;
   % Q10 method
   tref = 25;
   q10const = 1.7;   % Machea 2010
   tfact = q10const.^((tsoil-tref)/10);
%   tfact2 = 2.4.^((tsoil-tref)/10);
%   tfact3 = 2.5.^((tsoil-tref)/10);
   % Century Old
   t1 = (45-tsoil)/10;
   t2 = exp(0.076*(1-t1.^2.63));
   ft_old = (t1.^0.2).*t2;
   % Century method, defac
   tmod=0.73;
   ft=tmod*0.09*exp(0.095*tair);
   % Rothc method
   if(tair < -10.0)
      tair = -10.0;
   end
   t_rothc = 47.91 ./ (1.0 + exp(106.06./(tair+18.27)));
   t_rothc2 = 1.3*2.0.^((tsoil-10)/10)-1.3*2^-2;
 %  t_rothc = 1-exp(-47.91 ./ (1.0 + exp(106.06./(tair+18.27))));
 %  t_rothc2 = 1-exp(-2.0.^((tsoil-10)/10)+2^-2);
   plot(tair, tfact, 'r-');
   hold on
 %  plot(tair, tfact2, 'r--');
 %  plot(tair, tfact3, 'ro');
   plot(tair, ft ,'g-');
   plot(tair, t_rothc, 'b-');
   plot(tair, t_rothc2, 'b--');
   title('Comparison of different temperature modifiers');
   legend('Q10-based', 'Parton 1987', 'Jenkinson 1990', 'Q10-based derived from Gu 2004');
   xlabel('Temperature, degreeC', 'FontSize', 14);
   ylabel('Dimensionless Modifier', 'FontSize', 14);
 %  ylim([0,2]);
   hold off
 end

%% if(Moist_mod)
%%   moist = -10:0.2:40;
%%   % PET based modifier
%%   smod=0.5
%%
%%   % Treatment for negative PET
%%   if(pet<0.0) pet = 0.0
%%
%%   if(pet.eq.0.0)then
%%      defac_func_w=0.1
%%   else
%%   % wrat, which is the ratio of rainfall plus stored water to the pet
%%      wrat=(precip+soilwater)/pet
%%      fw= 1.0/( 1.0+30.*smod*exp(-8.5*wrat) )
%%      defac_func_w=fw
%%   end if
%%
%%   % Q10 method
%%   tref = 25;
%%   q10const = 2.0;   % Machea 2010
%%   tfact = q10const.^((tsoil-tref)/10);
%%%   tfact2 = 2.4.^((tsoil-tref)/10);
%%%   tfact3 = 2.5.^((tsoil-tref)/10);
%%   % Century Old
%%   t1 = (45-tsoil)/10;
%%   t2 = exp(0.076*(1-t1.^2.63));
%%   ft_old = (t1.^0.2).*t2;
%%   % Century method, defac
%%   tmod=0.73;
%%   ft=tmod*0.09*exp(0.095*tair);
%%   % Rothc method
%%   if(tair < -10.0)
%%      tair = -10.0;
%%   end
%%   t_rothc = 47.91 ./ (1.0 + exp(106.06./(tair+18.27)));
%%   t_rothc2 = 1.3*2.0.^((tsoil-10)/10)-1.3*2^-2;
%% %  t_rothc = 1-exp(-47.91 ./ (1.0 + exp(106.06./(tair+18.27))));
%% %  t_rothc2 = 1-exp(-2.0.^((tsoil-10)/10)+2^-2);
%%   plot(tair, tfact, 'r-');
%%   hold on
%% %  plot(tair, tfact2, 'r--');
%% %  plot(tair, tfact3, 'ro');
%%   plot(tair, ft ,'g-');
%%   plot(tair, t_rothc, 'b-');
%%   plot(tair, t_rothc2, 'b--');
%%   title('Comparison of different temperature modifiers');
%%   legend('Q10-based', 'Parton 1987', 'Jenkinson 1990', 'Q10-based derived from Gu 2004');
%%   xlabel('Temperature, degreeC', 'FontSize', 14);
%%   ylabel('Dimensionless Modifier', 'FontSize', 14);
%% %  ylim([0,2]);
%%   hold off
%% end

 if(Depth_mod)
   % Jenkinson 2008
   s = -8;
   f = z(1);
   F = (z - z(1));
   numerator = -1./(1+exp(-s.*(F-f)));
   denominator = -1/(1+exp(-s*(-f)));
   dfact = numerator/denominator;
   % Koven 2013
   zt = 1.0;
   dfact_koven = exp(-z./zt);

   plot(dfact, 'r-');
   hold on
   plot(dfact_koven, 'b-');
%   title('Comparison of different temperature modifiers');
%   legend('Q10-based', 'Parton 1987', 'Jenkinson 1990', 'Q10-based derived from Gu 2004');
%   xlabel('Temperature, degreeC', 'FontSize', 14);
%   ylabel('Dimensionless Modifier', 'FontSize', 14);
 %  ylim([0,2]);
   hold off
 end

 if(Turnover)
   % Bassil new *calibrated*
   % ISAM belowground (RothC), why these values needed to be modified into biome specific values?
   k_hum2=[0.01,0.007,0.022,0.012,0.01,0.015,0.009,0.011,0.012,0.0,0.0,0.01,0.009,0.01,0.007,0.02,0.01,0.01,0.0,0.012,0.009,0.01,0.01,0.012];
   k_dpm2=[5,5,3,10,20,10,5,5,20,0,0,5,5,5,5,3,10,20,0,20,5,5,5,20];
   k_rpm2=[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.0,0.0,0.3,0.3,0.3,0.3,0.3,0.3,0.2,0.0,0.2,0.3,0.3,0.3,0.2];
   k_bio2=0.66;
   t_hum2= 1./k_hum2;
   t_dpm2= 1./k_dpm2;
   t_rpm2= 1./k_rpm2;
   t_bio2= 1/k_bio2;

   % ISAM aboveground (CENTURY)
   k_str2 = 4.5;
   k_met2 = 18.5;
   k_bio_ag2 = 7.3;
   k_hum_ag2 = 2.0;
   t_str2 = 1/k_str2;
   t_met2 = 1/k_met2;
   t_bio_ag2 = 1/k_bio_ag2;
   t_hum_ag2 = 1/k_hum_ag2;

   % Xiaojuan's model
   % ISAM belowground (RothC)
   k_hum1=0.02;
   k_dpm1=10;
   k_rpm1=0.3;
   k_bio1=0.66;
   t_hum1= 1/k_hum1;
   t_dpm1= 1/k_dpm1;
   t_rpm1= 1/k_rpm1;
   t_bio1= 1/k_bio1;

   % ISAM aboveground (CENTURY)
   k_str1 = 4.5;
   k_met1 = 18.5;
   k_bio_ag1 = 7.3;
   k_hum_ag1 = 2.0;
   t_str1 = 1/k_str1;
   t_met1 = 1/k_met1;
   t_bio_ag1 = 1/k_bio_ag1;
   t_hum_ag1 = 1/k_hum_ag1;

 end

 if(CEC_Rh)
   clay = 0.01:0.01:1.0;
   cec = 46.0*clay;
   f_biohum =  1.0 ./ ( 2.67 * ( 1.21 + 2.24*exp(-0.085*cec) ) );
   f_co2 = 1-f_biohum;
   plot(clay,f_co2);
 end

 if(Therm_conduct)
   sand_vol = 10;
   porsl_min = 0.489 - 0.00126*sand_vol;
   bd_min = (1.0 - porsl_min)*2700;
   dkdry_min = (.135*bd_min + 64.7) / (2700 - 0.947*bd_min)
 end

 if(Percolation)
   sand_vol = 1:1:40;
   om_frac = 0.01:0.01:0.4;
   hksati  = zeros(length(om_frac),1);
   pc = 0.5;         % perlocation threshold
   pcbeta = 0.139;   % perlocation exponent
   om_hksat = 0.1;   % saturated hydraulic conductivity of organic soil [mm/s]

   xksat = 0.0070556 *( 10.^(-0.884+0.0153*sand_vol) ); % mm/s
 
   for i = 1:length(sand_vol)
     for j = 1:length(om_frac)
       if (om_frac(j) > pc)
          perc_norm=(1.0 - pc)^(-pcbeta);
          perc_frac=perc_norm*(om_frac(j) - pc)^pcbeta;
       else
          perc_frac=0.;
       end

       uncon_frac = (1.0-om_frac(j))+(1.0-perc_frac)*om_frac(j);
       % uncon_hksat is series addition of mineral/organic conductivites
       if (om_frac(j) < 1.0)
          uncon_hksat = uncon_frac/((1.0-om_frac(j))/xksat(i) + ((1.0-perc_frac)*om_frac(j))/om_hksat);
       else
          uncon_hksat = 0;
       end
       hksati(j)  = uncon_frac*uncon_hksat + (perc_frac*om_frac(j))*om_hksat;
     end

     plot(om_frac,hksati);
     if(i==1)
       hold on
     else
       if(i==length(sand_vol))
         hold off
       end
     end
   end
 end

 if(Cryoturbation)
   % Jenkinson 2008
   cryo = 1.;
   ald_perma = 0.5;
   mod = zeros(6,1);
   for i = 1:6
      mod(i) = cryo * (1 - (ald_perma-z(i))/ald_perma );
   end
   for i = 1:6
      mod_exp(i) = cryo * ( exp( log(2)*z(i)/ald_perma ) - 1);
   end

   plot(z(1:6), mod, 'r-');
   hold on
   plot(z(1:6), mod_exp, 'b-');
%   title('Comparison of different temperature modifiers');
%   legend('Q10-based', 'Parton 1987', 'Jenkinson 1990', 'Q10-based derived from Gu 2004');
%   xlabel('Temperature, degreeC', 'FontSize', 14);
%   ylabel('Dimensionless Modifier', 'FontSize', 14);
 %  ylim([0,2]);
   hold off
 end


