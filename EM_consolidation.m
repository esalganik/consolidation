% Figure 1: Spatial variability of ridge consolidation
% Alli's Ridge data
close all; clear; clc;
% Data from ice mass balance buoys (IMBs)
load('T61.mat',"t_T61","T_T61"); % data, 2020T61, doi:10.1594/PANGAEA.926580
load('PS_meteo.mat',"Ta_PS","t_PS","t_ASFS","Ta_ASFS_clean"); % data, meteo, doi:10.18739/A2FF3M18K
load('DTC26.mat',"t_DTC26","T_DTC_all_new","z_DTC_new"); % data, DTC26, doi:10.1594/PANGAEA.951780
z = 0:-0.02:-4.78; z = z + 0.8; z_T61 = z; % depth SIMBA T61
c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors
t_hc =  [  1   9  50 100 140 150 175 200 225 250 275 290 300 310 315 325 330 340 350 375 400 425 450 475 525 550 689];
hc_T = -[104 106 118 130 132 134 142 144 149 156 166 167 172 181 187 205 217 222 227 251 273 306 329 352 376 384 392]/100;
hc_T_int = interp1(datenum(t_T61(t_hc)),hc_T,datenum(t_T61),'linear'); % consolidated layer depth, T61
t_hs =  [ 1 50 100 200 250 275 325 350 385 485 600 620 640 650 659 675 689];
hs_T =  [70 70  66  70  68  68  74  74  74  70  66  36  36  32  30  14   8]/100;
hs_T_int = interp1(datenum(t_T61(t_hs)),hs_T,datenum(t_T61),'linear'); % snow-ice interface depth, T61
t_hc_DTC =  [  1  15  30  60 160 200 240 300 330 375 440 500 550 600 630 656];
hc_DTC =   -[210 210 210 210 214 222 224 238 238 242 238 256 242 238 224 214]/100;
hc_DTC_int = interp1(datenum(t_DTC26(t_hc_DTC)),hc_DTC,datenum(t_DTC26),'linear','extrap'); % consolidated layer depth, DTC
t_hs_DTC_dT = [836 840 842 844 847 859 865 884 909 936 940 950 953 965 969 973 977 980 984 989 992 994 997 998] + 737000;
hs_DTC_dT =   [052 054  62  72  88 090 086 084 086 088 080 078 080 078 056 054 044 030 030 022 016 014 014 010]/100;
hs_DTC_dT_int = interp1(t_hs_DTC_dT,hs_DTC_dT,datenum(t_DTC26),'linear','extrap'); % air-snow interface depth, DTC
t_hsi_DTC_dT = [836 898 939 951 959 962 965 972 973 977 983 989 998] + 737000;
hsi_DTC_dT =   [ 28  28  28  30  36  52  52  36  26  26  16   8   0]/100;
hsi_DTC_dT_int = interp1(t_hsi_DTC_dT,hsi_DTC_dT,datenum(t_DTC26),'linear','extrap'); % snow-ice interface depth, DTC
t_hk_DTC_dT = [836 889 907 925 946 958 965 971 983 994 998] + 737000;
hk_DTC_dT =  -[594 594 594 594 594 590 584 580 574 566 562]/100;
hk_DTC_dT_int = interp1(t_hk_DTC_dT,hk_DTC_dT,datenum(t_DTC26),'pchip','extrap'); % keel-water interface depth, DTC
for i = 1:length(t_DTC26) % time
    T_si_DTC(i) = T_DTC_all_new{60}(-1 + round((z_DTC_new(1) - hsi_DTC_dT_int(i))*50),i); % T at the snow-ice interface
end
load('T66.mat'); % Thermistor 2019T66 data import (FYI)
cut = 1097; t_T66 = t_T66(1:cut); T_T66 = T_T66(1:cut,:);
td_0_T66 = (datenum(t_T66)- datenum(t_T66(1))); % time of experiment [d]
bot_t_dir_T66 = [ 0  4  5 24 45 64  99 154 195 212 244 261 272]; % Direct time of measurements [d]
bot_dir_T66 =  -[37 44 46 58 70 86 108 144 158 158 150 142 130]/100 + 0.04; % Thickness for interpolation (m)
bot_int_T66 = interp1(bot_t_dir_T66,bot_dir_T66,td_0_T66,'pchip'); % Experimental bottom interpolation
top_t_dir_T66 = [0 212 215 223 225 233 241 246 254 262 266 269 272]; % Direct time of measurements [d]
top_dir_T66 =  -[0   0   6  12   0   2   0  12  22  40  42  52  54]/100; % Thickness for interpolation [m]
top_int_T66 = interp1(top_t_dir_T66,top_dir_T66,td_0_T66,'pchip'); % Experimental bottom interpolation

% new processing
sn_DTC = hs_DTC_dT_int-hsi_DTC_dT_int; % snow thickness, from DTC
t_MP = datetime(['17-Jan-2020';'31-Jan-2020';'28-Feb-2020';'10-Apr-2020';'28-Jun-2020';'13-Jul-2020';'28-Jul-2020']); % time
h_sn_MP_DTC = [0.02 0.02 0.56 0.61 0.24 0.15 0.04]; % snow thickness, DTC site from Magnaprobe, doi:10.1594/PANGAEA.937781
h_sn = interp1([datenum(t_T66(270)); datenum(t_MP(1:2)); datenum(t_DTC26(4:end)); datenum(t_T66(end))],[0; h_sn_MP_DTC(1:2)'; sn_DTC(4:end); 0],datenum(t_T66),'linear');
h_sn_MP_T61 = [0.09 0.19 0.45 0.98 0.44 0.07 0.07]; % snow thickness, T61 site from Magnaprobe
h_sn_T61 = interp1([datenum(t_T66(270)); datenum(t_MP); datenum(t_T66(end))],[0; h_sn_MP_T61'; 0],datenum(t_T66),'linear');

% Data from ridge coring: salinity and density
% Density for Alli's ridge, T61 temperatures
load('Coring_AR.mat'); % data, ridge coring, doi:10.1594/PANGAEA.953865
z = z_T61; % depth SIMBA T61
S_rho = interp1(z_S,S,z,'linear','extrap');
rho_int = interp1(z_S,rho_lab,z,'linear','extrap');
T_lab = ones(size(T_T61)) * -15; % Lab temperature
F1_pr_rho = -4.732-22.45*T_lab - 0.6397*T_lab.^2 - 0.01074*T_lab.^3;    
F2_pr_rho = 8.903*10^-2 - 1.763*10^-2*T_lab - 5.33*10^-4*T_lab.^2 - 8.801*10^-6*T_lab.^3;
vb_pr_rho = rho_int .* S_rho ./ F1_pr_rho; % brine volume for T_lab
rhoi_pr = (917-1.403*10^-1*T_lab); % pure ice density, Pounder (1965)
vg_pr = max((1-rho_int.*(F1_pr_rho-rhoi_pr.*S_rho/1000.*F2_pr_rho)./(rhoi_pr.*F1_pr_rho))*1000,0,'includenan'); % gas volume for T_lab    
F3_pr = rhoi_pr.*S_rho/1000./(F1_pr_rho-rhoi_pr.*S_rho/1000.*F2_pr_rho);
rhoi_rho = (917-1.403*10^-1*T_T61); % pure ice density, Pounder (1965) for T_insitu
F1_rho = -4.732-22.45*T_T61 - 0.6397*T_T61.^2 - 0.01074*T_T61.^3; % Cox and Weeks (1983)  
F1_rho(T_T61>-2) = -4.1221*10^-2 + -1.8407*10^1*T_T61(T_T61>-2).^1 + 5.8402*10^-1*T_T61(T_T61>-2).^2 + 2.1454*10^-1*T_T61(T_T61>-2).^3; % F1 from Lepparanta and Manninen (1988)
F2_rho = 8.903*10^-2 - 1.763*10^-2*T_T61 - 5.33*10^-4*T_T61.^2 - 8.801*10^-6*T_T61.^3; % Cox and Weeks (1983)
F2_rho(T_T61>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_T61(T_T61>-2).^1 + 1.2291*10^-4*T_T61(T_T61>-2).^2 + 1.3603*10^-4*T_T61(T_T61>-2).^3; % F2 from Lepparanta and Manninen (1988)
F3_rho = rhoi_rho.*S_rho/1000./(F1_rho-rhoi_rho.*S_rho/1000.*F2_rho);
vb_rho = vb_pr_rho .* F1_pr_rho ./ F1_rho; % Brine volume for T_insitu, CW + LM
vg = max(0,(1-(1-vg_pr/1000).*(rhoi_rho./rhoi_pr).*(F3_pr.*F1_pr_rho./F3_rho./F1_rho))*1000,'includenan'); % Gas volume for T_insitu, CW + LM
rho_si = (-vg/1000+1).*rhoi_rho.*F1_rho./(F1_rho-rhoi_rho.*S_rho/1000.*F2_rho); % Ice density for T_insitu

% Density for Alli's ridge, DTC temperatures
T_DTC = T_DTC_all_new{60}(1:368,:);
T_DTC_clean = T_DTC_all_new{60}(1:368,4:end);
indices = find(T_DTC_clean>10); T_DTC_clean(indices) = NaN;
filter = ~(any(isnan(T_DTC_clean(:,:)),1)); filter(422) = 0; filter(1:3) = 0;
t_DTC_clean = t_DTC26(filter); T_DTC_clean = T_DTC_clean(:,filter); T_DTC_clean = T_DTC_clean';
z = z_DTC_new; % depth DTC 26 (368 elements)
S_rho = interp1(z_S,S,z,'linear','extrap');
rho_int = interp1(z_S,rho_lab,z,'linear','extrap');
T_lab = ones(size(T_DTC_clean)) * -15; % Lab temperature
F1_pr_rho = -4.732-22.45*T_lab - 0.6397*T_lab.^2 - 0.01074*T_lab.^3;    
F2_pr_rho = 8.903*10^-2 - 1.763*10^-2*T_lab - 5.33*10^-4*T_lab.^2 - 8.801*10^-6*T_lab.^3;
vb_pr_rho = rho_int .* S_rho ./ F1_pr_rho; % brine volume for T_lab
rhoi_pr = (917-1.403*10^-1*T_lab); % pure ice density, Pounder (1965)
vg_pr = max((1-rho_int.*(F1_pr_rho-rhoi_pr.*S_rho/1000.*F2_pr_rho)./(rhoi_pr.*F1_pr_rho))*1000,0,'includenan'); % gas volume for T_lab    
F3_pr = rhoi_pr.*S_rho/1000./(F1_pr_rho-rhoi_pr.*S_rho/1000.*F2_pr_rho);
rhoi_rho = (917-1.403*10^-1*T_DTC_clean); % pure ice density, Pounder (1965) for T_insitu
F1_rho = -4.732-22.45*T_DTC_clean - 0.6397*T_DTC_clean.^2 - 0.01074*T_DTC_clean.^3; % Cox and Weeks (1983)  
F1_rho(T_DTC_clean>-2) = -4.1221*10^-2 + -1.8407*10^1*T_DTC_clean(T_DTC_clean>-2).^1 + 5.8402*10^-1*T_DTC_clean(T_DTC_clean>-2).^2 + 2.1454*10^-1*T_DTC_clean(T_DTC_clean>-2).^3; % F1 from Lepparanta and Manninen (1988)
F2_rho = 8.903*10^-2 - 1.763*10^-2*T_DTC_clean - 5.33*10^-4*T_DTC_clean.^2 - 8.801*10^-6*T_DTC_clean.^3; % Cox and Weeks (1983)
F2_rho(T_DTC_clean>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_DTC_clean(T_DTC_clean>-2).^1 + 1.2291*10^-4*T_DTC_clean(T_DTC_clean>-2).^2 + 1.3603*10^-4*T_DTC_clean(T_DTC_clean>-2).^3; % F2 from Lepparanta and Manninen (1988)
F3_rho = rhoi_rho.*S_rho/1000./(F1_rho-rhoi_rho.*S_rho/1000.*F2_rho);
vb_rho = vb_pr_rho .* F1_pr_rho ./ F1_rho; % Brine volume for T_insitu, CW + LM
vb_rho(vb_rho < 0) = 0; vb_rho(vb_rho > 400) = 400;
vg = max(0,(1-(1-vg_pr/1000).*(rhoi_rho./rhoi_pr).*(F3_pr.*F1_pr_rho./F3_rho./F1_rho))*1000,'includenan'); % Gas volume for T_insitu, CW + LM
rho_si_DTC = (-vg/1000+1).*rhoi_rho.*F1_rho./(F1_rho-rhoi_rho.*S_rho/1000.*F2_rho);
Sb = -1.2 - 21.8*T_DTC_clean - 0.919*T_DTC_clean.^2 - 0.0178*T_DTC_clean.^3; % brine salinity for seawater (Assur, 1958)
k_i = 2.21 - 10^-2*T_DTC_clean + 3.44*10^-5*T_DTC_clean.^2; % pure ice thermal conductivity (Yen at al., 1991)
k_b = 0.52325*(1-Sb./1000)+0.01256*T_DTC_clean + 5.8604*10^-5*T_DTC_clean.^2;  % seawater thermal conductivity (Schwerdtfeger, 1963) 
k_si = (1-vb_rho/1000-vg/1000).*k_i + vb_rho/1000.*k_b; % sea ice thermal conductivity (Batchelor, 1974)
T_DTC_cl = T_DTC_clean; % temperature of only frozen sesnors
for i = 1:size(T_DTC_clean,1) % time
    for j = 1:size(T_DTC_clean,2) % depth
        if  j > (z_DTC_new(1) - hc_DTC_int(i))*100/2
            T_DTC_cl(i,j) = NaN;
        elseif  j < (z_DTC_new(1) - z_DTC_new(65))*100/2
            T_DTC_cl(i,j) = NaN;
        end
    end
end
rho_DTC_cl = rho_si_DTC; rho_DTC_cl(isnan(T_DTC_cl)) = NaN; rho_bulk = nanmean(rho_DTC_cl,2); % bulk density of CL
vb_DTC_cl = vb_rho; vb_DTC_cl(isnan(T_DTC_cl)) = NaN; vb_bulk = nanmean(vb_DTC_cl,2)/1000; % bulk brine volume of CL
k_si_DTC_cl = k_si; k_si_DTC_cl(isnan(T_DTC_cl)) = NaN; k_si_bulk = nanmean(k_si_DTC_cl,2); % bulk thermal conductivity of CL

% Atmospheric data
Ta_PS_cor = Ta_PS(14902:454102); t_PS_cor = t_PS(14902:454102); t_PS_cor(Ta_PS_cor == 9) = []; Ta_PS_cor(Ta_PS_cor == 9) = [];
Ta_PS_6h = accumarray(ceil((1:numel(Ta_PS_cor))/360)',Ta_PS_cor(:),[],@mean); t_PS_a_6h = t_PS_cor(1:360:end); % 6h mean Ta from Polarstern
t_ASFS_hour = t_ASFS(1:360:end);
Ta_ASFS_hour = accumarray(ceil((1:numel(Ta_ASFS_clean))/360)',Ta_ASFS_clean(:),[],@mean);
Ta_T66 = interp1(datenum(t_PS_a_6h),Ta_PS_6h,datenum(t_T66),'linear');
Ta_T66(800:923) = interp1(datenum(t_ASFS_hour),Ta_ASFS_hour,datenum(t_T66(800:923)),'linear');

% Ocean data
project = 'MOSAiC_dailyaverage_skalars.nc'; % MSS ocean data, doi:10.18739/A21J9790B
SST = ncread(project,'SST'); time = ncread(project,'day'); t_0 = (datetime('0000-01-01 00:00:00')); t_MSS = t_0 + days(time-1);
Tw_mss = interp1(datetime(t_MSS),SST,datetime(t_T66),'linear');
clearvars -except c t_T66 Ta_T66 t_T61 t_DTC_clean T_DTC_clean hs_DTC_dT_int hsi_DTC_dT_int t_DTC26 ...
    T_si_DTC hc_T_int hc_DTC_int Tw_mss h_sn h_sn_T61 t_MP h_sn_MP_DTC h_sn_MP_T61 rho_bulk vb_bulk k_si_bulk vg_bulk

% Consolidation modelling, model description doi:10.1016/j.coldregions.2020.103194
% DTC site, fresh ice model
por = 0.29; % ridge macroporosity (from drilling)
Ta = Ta_T66(270:1097); Tf = Tw_mss(270:1097); % air and water temperatures
Hia = 21; % heat transfer coefficient (function of wind speed)
hs = h_sn(270:1097); % snow thickness
ks = 0.26; % snow thermal conductivity, Macfarlane et al., doi:10.5194/tc-17-5417-2023
ki = 2.2; rhoi = 917; Li = 333400; % fresh ice thermodynamic parameters
dc_0 = 0.20; % initial ice thickness = initial ice freeboard from drilling
td = datenum(t_T66(270:1097))-datenum(t_T66(270)); t = (td)*24*3600; t_sim = t_T66(270:1097); % time
% fresh ice model
dt = diff(t); n = length(t); [dc,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc(1) = dc_0;
for i = 1:n-1
     R1(i) = 1./Hia; R2(i) = hs(i)/ks; R3(i) = dc(i)/ki;
     Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
     ddc(i) = -ki/rhoi/(Li*por)*(Tsi(i) - Tf(i))/dc(i)*dt(i);
     dc(i+1) = dc(i) + ddc(i);
end
dc = dc - dc_0; % CL thickness = ice thickness - initial freeboard
% DTC site, saline ice model
vb_bulk_int = interp1([datenum(t_T66(270)); datenum(t_DTC_clean); datenum(t_T66(end));],[vb_bulk(1); vb_bulk; vb_bulk(end)],datenum(t_T66(270:1097)),'linear','extrap');
rho_bulk_int = interp1([datenum(t_T66(270)); datenum(t_DTC_clean); datenum(t_T66(end));],[rho_bulk(1); rho_bulk; rho_bulk(end)],datenum(t_T66(270:1097)),'linear','extrap');
k_si_bulk_int = interp1([datenum(t_T66(270)); datenum(t_DTC_clean); datenum(t_T66(end));],[k_si_bulk(1); k_si_bulk; k_si_bulk(end)],datenum(t_T66(270:1097)),'linear','extrap');
Lsi = 333400*(1-vb_bulk_int); rhoi = rho_bulk_int; ki = k_si_bulk_int; % saline ice parameters
dt = diff(t); n = length(t); [dc_si,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc_si(1) = dc_0;
for i = 1:n-1
     R1(i) = 1./Hia; R2(i) = hs(i)/ks; R3(i) = dc_si(i)./ki(i);
     Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
     ddc(i) = -ki(i)./rhoi(i)./(Li*por).*(Tsi(i) - Tf(i))/dc_si(i)*dt(i);
     dc_si(i+1) = dc_si(i) + ddc(i);
end
dc_si = dc_si - dc_0; % CL thickness = ice thickness - initial freeboard
dc_si = dc_si*Li./Lsi'; % correction for CL microporosity

% Fresh ice model, T61 site
ki = 2.2; rhoi = 917; Li = 333400; % fresh ice thermodynamic parameters
hs = h_sn_T61(270:1097); % snow thickness (from Magnaprobe)
dc_0 = 0.72; % initial ice thickness = initial ice freeboard from drilling
td = datenum(t_T66(270:1097))-datenum(t_T66(270)); t = (td)*24*3600; t_sim = t_T66(270:1097); % time
dt = diff(t); n = length(t); [dc_2,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc_2(1) = dc_0;
for i = 1:n-1
     R1(i) = 1./Hia; R2(i) = hs(i)/ks; R3(i) = dc_2(i)/ki;
     Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
     ddc(i) = -ki/rhoi/(Li*por)*(Tsi(i) - Tf(i))/dc_2(i)*dt(i);
     dc_2(i+1) = dc_2(i) + ddc(i);
end
dc_2 = dc_2 - dc_0;
% saline ice model, T61 site
Lsi = 333400*(1-vb_bulk_int); rhoi = rho_bulk_int; ki = k_si_bulk_int;
dc_0 = 0.72; td = datenum(t_T66(270:1097))-datenum(t_T66(270)); t = (td)*24*3600; t_sim = t_T66(270:1097);
dt = diff(t); n = length(t); [dc_2_si,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc_2_si(1) = dc_0;
for i = 1:n-1
     R1(i) = 1./Hia; R2(i) = hs(i)/ks; R3(i) = dc_2_si(i)./ki(i);
     Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
     ddc(i) = -ki(i)./rhoi(i)/(Li*por)*(Tsi(i) - Tf(i))/dc_2_si(i)*dt(i);
     dc_2_si(i+1) = dc_2_si(i) + ddc(i);
end
dc_2_si = dc_2_si - dc_0;
dc_2_si = dc_2_si*333400./Lsi';
clearvars ddc dc_0 Li ks rhoi por R1 R2 R3 Hia dt n i Ta Tsi ki Lsi t

% Ridge drilling data, doi:10.1594/PANGAEA.950192
x1 = [ -2.5     0   2.5     5   7.5    10  12.5    15  17.5];
k1 = [-4.25 -4.89 -5.53 -6.16 -6.79 -6.37 -5.95 -5.55 -5.15];
c1 = [-2.75 -2.29 -1.82 -1.34 -0.89 -1.12 -1.35 -1.70 -2.05];
s1 = [ 1.25  1.12  0.99  0.85  0.72  0.61  0.50  0.35  0.20];
sn1 =[ 1.28  1.16  1.04  0.91  0.79  0.90  1.00  0.61  0.21];
x2 = [ -2.5     0   2.5     5   7.5    10  12.5    15  17.5]; x0 = [17.5    20  22.5];
c2 = [-4.70 -3.00 -2.90 -2.55 -4.76 -5.40 -5.10 -4.05 -2.20]; k0 = [-5.15 -0.97 -0.86];
k2 = [-4.70 -3.00 -3.38 -5.80 -5.60 -5.40 -5.10 -4.90 -4.53]; c0 = [-2.05 -0.97 -0.86];
s2 = [ 0.80  1.75  1.06  0.80  0.65  0.53  0.35  0.15  0.00]; s0 = [ 0.20  0.13  0.09];
sn2= [ 1.17  1.83  1.11  0.82  0.68  0.57  0.37  0.20  0.00]; sn0 =[ 0.21  0.14  0.16];

% Figure 1
figure
tile = tiledlayout(2+2,3); tile.TileSpacing = 'compact'; tile.Padding = 'none';
nexttile([2 1]) % ridge morphology from drilling
msz = 2.5; L1 = ([1 5 7 9]);
p = plot(x1,-s1,'b-.','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = msz; hold on
p = plot(x1,-k1,'b--','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = msz;
p = plot(x1(L1),-c1(L1),'bo-','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = msz;
p = plot(x2,-c2,'ro-','color',c{2},'MarkerFaceColor','w'); p.MarkerSize = msz; 
p = plot(x1(L1),-k1(L1),'bo--','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = msz;
p = plot(x2,-k2,'ro--','color',c{2},'MarkerFaceColor','w'); p.MarkerSize = msz; 
p = plot(x1(L1),-s1(L1),'bo-.','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = msz; 
p = plot(x2,-s2,'ro-.','color',c{2},'MarkerFaceColor','w'); p.MarkerSize = msz;
p = plot(x0,-s0,'bo-.','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = msz;
p = plot(x0,-k0,'bo--','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = msz;
p = plot(x0,-c0,'bo-','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = msz;
p = fill([[x1 x0] fliplr([x1 x0])],[[-c1 -c0]*0 fliplr([-c1 -c0])],'b','facecolor',c{5},'edgecolor','none'); set(p,'facealpha',.1); % shaded initial cons.
p = fill([x1 fliplr(x2)],[-c1 fliplr(-c2)],'b','facecolor',c{1},'edgecolor','none'); set(p,'facealpha',.1); % shaded cons.
p = fill([x1 fliplr(x2)],[-k1 fliplr(-k2)],'r','facecolor',c{2},'edgecolor','none'); set(p,'facealpha',.1); % shaded melt
p = text(7.5,-1.7+0.5,'T61'); set(p,'Color',c{3},'HorizontalAlignment','center','FontSize',8);
p = text(17.5,-1.7+0.9,'DTC26'); set(p,'Color',c{4},'HorizontalAlignment','center','FontSize',8);
p = text(-2.0,0.4,'Consolidation'); set(p,'Color',c{5},'HorizontalAlignment','left','FontSize',8);
p = text(-2.0,1.0,'Jan-Feb'); set(p,'Color',c{5},'HorizontalAlignment','left','FontSize',8);
p = text(11.0,2.7,'Consolidation'); set(p,'Color',c{1},'HorizontalAlignment','center','FontSize',8);
p = text(11.0,3.3,'Feb-Jul'); set(p,'Color',c{1},'HorizontalAlignment','center','FontSize',8);
p = text(8.8,5.95,'Melt'); set(p,'Color',c{2},'HorizontalAlignment','center','FontSize',8);
leg = legend('Sail, Jan-Feb','Sail, July','CL, Jan-Feb','CL, July','Keel, Jan-Feb','Keel, July','box','off','NumColumns',2);set(leg,'FontSize',6,'Location','south');
xlim([-5 25]); ylim([-2 10]); leg.ItemTokenSize = [30*0.72,18*0.72];
hXLabel = xlabel('x (m)'); hYLabel = ylabel('Ice draft (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');
clearvars x1 k1 c1 s1 sn1 x2 k2 c2 s2 sn2 x0 k0 c0 s0 sn0 L1 msz

nexttile % Temperatures
plot(t_T66(270:1097),movmean(Tw_mss(270:1097),20),'color',c{5}); hold on  % water temperature
plot(t_DTC26(4:end),movmean(T_si_DTC(4:end),20),'color',c{2}); % snow-ice temperature
plot(t_T66(270:1097),movmean(Ta_T66(270:1097),20),'color',c{1}); hold on % air temperature
leg = legend('Water','Snow-ice, DTC26','Air','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.2,18*0.2];
hYLabel = ylabel('Temperature (°C)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([-40 10]); yticks(-40:10:10);
ax = gca; ax.XTick = datetime(['01-Jan-2020';'01-Feb-2020';'01-Mar-2020';'01-Apr-2020';'01-May-2020';'01-Jun-2020';'01-Jul-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile([2 1]) % Consolidation: model and observations
plot(t_T61,-hc_T_int,':','linewidth',3,'color',c{3}); hold on % observations, T61
plot(t_sim,dc_2_si,'-','linewidth',1,'color',c{3}); % saline ice model, T61
plot(t_DTC26(4:end),-hc_DTC_int(4:end),':','linewidth',3,'color',c{4});  % observations, DTC
plot(t_sim,dc_si,'-','linewidth',1,'color',c{4}); % saline ice model, DTC
leg = legend('CL, obs., T61','CL, model, T61','CL, obs., DTC26','CL, model, DTC26','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('Ice draft (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');
ax = gca; ax.XTick = datetime(['01-Jan-2020';'01-Feb-2020';'01-Mar-2020';'01-Apr-2020';'01-May-2020';'01-Jun-2020';'01-Jul-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile % Snow thickness
p = plot(t_MP,h_sn_MP_T61,'o--','color',c{3}); p.MarkerSize = 2.3; set(p,'markerfacecolor',get(p,'color')); hold on % snow depth from MP, T61 site
p = plot(t_MP,h_sn_MP_DTC,'o','color',c{4}); p.MarkerSize = 2.3; set(p,'markerfacecolor',get(p,'color')); % snow depth from MP, DTC site
plot(t_T66,h_sn,'color',c{4}); % snow depth from DTC
leg = legend('T61, Magnaprobe','DTC26, Magnaprobe','DTC26, temperature','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.6,18*0.6];
hYLabel = ylabel('Snow thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([0 2]); 
ax = gca; ax.XTick = datetime(['01-Jan-2020';'01-Feb-2020';'01-Mar-2020';'01-Apr-2020';'01-May-2020';'01-Jun-2020';'01-Jul-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time
clearvars leg hXLabel hYLabel tile p ax

nexttile([2 1]) % ridge morphology from drilling, Fort Ridge
x = [+12     +17   +19   +22   +27   +30];
k = [-0.74 -5.73 -6.78 -5.60 -3.91 -3.86];
cl = [-0.74 -1.33 -2.23 -1.00 -1.01 -1.06];
% cl_T = [NaN -1.6 NaN -4.0 -2.5 NaN];
cl_T = [NaN -1.33-0.62 NaN -1.0-3.03 -1.01-1.07 NaN]; % consolidation from IMBs
s = [0.09   0.32  0.32  0.50  0.19  0.19];
sn =[0.09   0.82  0.47  0.60  1.00  0.78];
p = plot(x,-s,'o-.','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = 2.5; hold on
p = plot(x,-cl,'o-','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = 2.5;
p = plot(x(~isnan(cl_T)),-cl_T(~isnan(cl_T)),'o-','color',c{2},'MarkerFaceColor','w'); p.MarkerSize = 2.5;
p = plot(x,-k,'o--','color',c{1},'MarkerFaceColor','w'); p.MarkerSize = 2.5;
p = fill([x fliplr(x)],[-cl*0 fliplr(-cl)],'b','facecolor',c{5},'edgecolor','none'); set(p,'facealpha',.1); % shaded initial cons.
p = fill([x(2:5) fliplr(x(~isnan(cl_T)))],[-cl(2:5) fliplr(-cl_T(~isnan(cl_T)))],'b','facecolor',c{1},'edgecolor','none'); set(p,'facealpha',.1); % shaded cons.
p = text(12.5,0.3,'Consolidation, Oct-Jan'); set(p,'Color',c{5},'HorizontalAlignment','left','FontSize',8);
p = text(21.2,1.5,'Consolidation'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',8);
p = text(21.2,2.1,'Jan-Apr'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',8);
p = text(17,-1.4,'DTC'); set(p,'Color',c{3},'HorizontalAlignment','center','FontSize',8);
p = text(17,-0.8,'25'); set(p,'Color',c{3},'HorizontalAlignment','center','FontSize',8);
p = text(22,-1.4,'T60'); set(p,'Color',c{4},'HorizontalAlignment','center','FontSize',8);
p = text(27,-1.4,'DTC'); set(p,'Color',c{6},'HorizontalAlignment','center','FontSize',8);
p = text(27,-0.8,'24'); set(p,'Color',c{6},'HorizontalAlignment','center','FontSize',8);
leg = legend('Sail, Jan','CL, Jan','CL, Mar-Apr','Keel, Jan','box','off','NumColumns',2);set(leg,'FontSize',6,'Location','south');
ylim([-2 10]); leg.ItemTokenSize = [30*0.72,18*0.72]; % xlim([-5 25]); 
hXLabel = xlabel('x (m)'); hYLabel = ylabel('Ice draft (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');

nexttile % temperatures, Fort Ridge
plot(t_T66(1:1097),movmean(Tw_mss(1:1097),20),'color',c{5}); hold on  % water temperature
plot(t_T66(1:1097),movmean(Ta_T66(1:1097),20),'color',c{1}); hold on % air temperature
leg = legend('Water','Air','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.2,18*0.2];
hYLabel = ylabel('Temperature (°C)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([-40 10]); yticks(-40:10:10);
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Nov-2019';'01-Dec-2019';'01-Jan-2020';'01-Feb-2020';'01-Mar-2020';'01-Apr-2020';'01-May-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile([2 1]); % consolidation
load('DTC_data.mat',"t_DTC","fb_DTC","hs_DTC","hi_DTC");
load('T60.mat',"t_T60");
t_hc =  [  1   9  50 100 150 200 250 275 285 290 310 315 325 350 375 395 405 410 420];
hc_T = -[147 147 177 200 214 233 256 280 319 323 354 374 412 429 448 450 450 450 450]/100;
hc_T_int = interp1(datenum(t_T60(t_hc)),hc_T,datenum(t_T60),'linear');
i = 12; plot(datetime(t_DTC{i},'ConvertFrom','datenum'),-hi_DTC{i},':','linewidth',3,'Color',c{3}); hold on % DTC 25
plot(t_T60,-hc_T_int-0.5,':','linewidth',3,'Color',c{4}); % T60
i = 11; plot(datetime(t_DTC{i},'ConvertFrom','datenum'),-hi_DTC{i},':','linewidth',3,'Color',c{6}); % DTC 24
leg = legend('CL, obs., DTC25','CL, obs., T60','CL, obs., DTC24','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('Ice draft (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Nov-2019';'01-Dec-2019';'01-Jan-2020';'01-Feb-2020';'01-Mar-2020';'01-Apr-2020';'01-May-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time
ylim([0 4.5]); 

nexttile % Snow thickness
i = 12; plot(datetime(t_DTC{i},'ConvertFrom','datenum'),hs_DTC{i}-fb_DTC{i},'Color',c{3}); hold on % snow depth from DTC 25
i = 11; plot(datetime(t_DTC{i},'ConvertFrom','datenum'),hs_DTC{i}-fb_DTC{i},'Color',c{6}); % snow depth from DTC 24
leg = legend('DTC25','DTC24','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.6,18*0.6];
hYLabel = ylabel('Snow thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([0 2]); 
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Nov-2019';'01-Dec-2019';'01-Jan-2020';'01-Feb-2020';'01-Mar-2020';'01-Apr-2020';'01-May-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time
clearvars leg hXLabel hYLabel tile p ax

% annotation('textbox',[0.005 .51 0.01 .51],'String','(a)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
% annotation('textbox',[0.220 .51 0.22 .51],'String','(b)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
% annotation('textbox',[0.440 .51 0.45 .51],'String','(c)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
% annotation('textbox',[0.220 .28 0.22 .28],'String','(d)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.005 .50 0.01 .51],'String','(a)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.215 .50 0.22 .51],'String','(b)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.440 .50 0.45 .51],'String','(c)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.005 .25 0.01 .26],'String','(d)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.215 .25 0.22 .26],'String','(e)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.440 .25 0.45 .26],'String','(f)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');

%% Figure 2: Comparison of consolidated layer thickness evolution from EMI, model, and IMBs
% Ridge EMI data from Polona, Alli's Ridge, transect A1
clear; clc;
load('PS_meteo.mat',"Tw_mss","Ta_T66"); % data, meteo, doi:10.18739/A2FF3M18K
load('T66.mat',"t_T66"); t_T66 = t_T66(1:1097);
t_sim = t_T66(270:1097); Tw_sim = Tw_mss(270:1097); Ta_sim = Ta_T66(270:1097);
t_sim = t_sim(53:end); Tw_sim = Tw_sim(53:end); Ta_sim = Ta_sim(53:end); % only after first EM survey

% import of EMI data
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA1_cc_consolidatedlayer.csv";
T = readtable(project); A = table2array(T); x = A(1,2:end); cl = A(2:size(A,1),2:end); time = num2str(A(2:size(A,1),1)); t = datetime(time,'InputFormat','yyyyMMdd'); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA1_cc_freeboard.csv";
T = readtable(project); A = table2array(T); fb = A(2,2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA1_cc_snow.csv";
T = readtable(project); A = table2array(T); sn = A(2:size(A,1),2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA1_cc_totalthickness.csv";
T = readtable(project); A = table2array(T); k = A(2:size(A,1),2:end); clearvars project T A

for i = 1:length(x)
    hs_int(:,i) = interp1([datenum(t_sim(1)); datenum(t); datenum(t_sim(end))],[0; sn(:,i); 0],datenum(t_sim),'linear');
end

% Alli's Ridge EMI transect, fresh ice model
por = 0.29; % ridge macroporosity (from drilling)
Ta = Ta_sim; Tf = Tw_sim; % air and water temperatures
Hia = 21; % heat transfer coefficient (function of wind speed)
hs = hs_int; % snow thickness
ks = 0.26; % snow thermal conductivity, Macfarlane et al., doi:10.5194/tc-17-5417-2023
ki = 2.2; rhoi = 917; Li = 333400; % fresh ice thermodynamic parameters
% dc_0 = max(fb,0.01); % initial ice thickness = initial ice freeboard from EMI
dc_0 = cl(1,:); % initial ice thickness = initial CL thickness from EMI
td = datenum(t_sim)-datenum(t_sim(1)); t_sec = (td)*24*3600; % time
% fresh ice model
dt = diff(t_sec); n = length(t_sec); [dc,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc_em = hs_int*0;
for j = 1:length(x)
    dc(1) = dc_0(j);
    for i = 1:n-1
        R1(i) = 1./Hia; R2(i) = hs(i,j)/ks; R3(i) = dc(i)/ki;
        Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
        ddc(i) = -ki/rhoi/(Li*por)*(Tsi(i) - Tf(i))/dc(i)*dt(i);
        dc(i+1) = dc(i) + ddc(i);
    end
    % dc = dc - dc_0(j); % CL thickness = ice thickness - initial freeboard
    dc = dc - fb(j); % CL thickness = ice thickness - initial freeboard
    dc_em(:,j) = dc(:);
end

c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors
for i = 1:length(t); [~,t_em(i)] = min(abs(datenum(t_sim)-datenum(t(i)))); end

figure
tile = tiledlayout(1,5); tile.TileSpacing = 'compact'; tile.Padding = 'none';
nexttile
for i = 1:length(t)
    di_em_avg(i) = mean(cl(i,:)-fb(1,:)); di_em_std(i) = std(cl(i,:));
end
for i = 1:length(t_sim)
    di_mdl_avg(i) = mean(dc_em(i,:)); di_mdl_std(i) = std(dc_em(i,:));
end
p = fill([t_sim; flipud(t_sim)],[di_mdl_avg+di_mdl_std fliplr(di_mdl_avg-di_mdl_std)],1,'FaceColor',c{1},'edgecolor','none'); set(p,'facealpha',0.1); hold on
errorbar(t,di_em_avg,di_em_std,'color',c{4}); hold on
p = plot(t,di_em_avg,'o','color',c{4}); p.MarkerSize = 3; set(p,'markerfacecolor',get(p,'color'));
plot(t_sim,di_mdl_avg,'color',c{1});

% IMB data
load('T61.mat',"t_T61"); load('DTC26.mat',"t_DTC26");
t_hc =  [  1   9  50 100 140 150 175 200 225 250 275 290 300 310 315 325 330 340 350 375 400 425 450 475 525 550 689];
hc_T = -[104 106 118 130 132 134 142 144 149 156 166 167 172 181 187 205 217 222 227 251 273 306 329 352 376 384 392]/100;
hc_T_int = interp1(datenum(t_T61(t_hc)),hc_T,datenum(t_T61),'linear'); % consolidated layer depth, T61
t_hc_DTC =  [  1  15  30  60 160 200 240 300 330 375 440 500 550 600 630 656];
hc_DTC =   -[210 210 210 210 214 222 224 238 238 242 238 256 242 238 224 214]/100;
hc_DTC_int = interp1(datenum(t_DTC26(t_hc_DTC)),hc_DTC,datenum(t_DTC26),'linear','extrap'); % consolidated layer depth, DTC
plot(t_T61,-hc_T_int,'--','linewidth',0.5,'color',c{3}); % observations, T61
plot(t_DTC26(4:end),-hc_DTC_int(4:end),'-','linewidth',0.5,'color',c{3}); % observations, DTC

hYLabel = ylabel('Thickness of consolidated layer (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([1 5]);
title('A1','FontSize',8,'FontWeight','normal');
t_start = datetime('01-Jan-2020'); t_end = datetime('01-Aug-2020'); xlim([t_start t_end]); datetick('x','mmm','keeplimits'); xtickangle(0); set(gca,'YDir','reverse');
% leg = legend('','EMI','','Model','IMB','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','southeast'); leg.ItemTokenSize = [30*0.66,18*0.66];

% figure
% tile = tiledlayout(2,3); tile.TileSpacing = 'compact'; tile.Padding = 'none';
% nexttile
% for i = 2:length(t)
%     p = plot(cl(i,:)-cl(1,:),dc_em(t_em(i),:)-dc_em(t_em(1),:),'o','color',c{i}); p.MarkerSize = 3.0; set(p,'markerfacecolor',get(p,'color')); hold on
% end
% plot([0 max(cl(end,:)-cl(1,:))],[0 max(cl(end,:)-cl(1,:))],'k--');
% leg = legend(datestr(t(2:end),'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% hXLabel = xlabel('EMI consolidation (m)'); hYLabel = ylabel('Model consolidation (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
% title('A1','FontSize',8,'FontWeight','normal');

% figure
% for i = 1:length(t); [~,t_em(i)] = min(abs(datenum(t_sim)-datenum(t(i)))); end
% for i = 1:length(t)
%     plot(x,dc_em(t_em(i),:),'color',c{i}); hold on
% end
% for i = 1:length(t)
%     plot(x,cl(i,:)-fb,':','linewidth',3.5,'color',c{i});
% end
% plot(x,mean(k(1,:),1)-fb,'-.','linewidth',1,'color','k');
% plot(x,-mean(sn,1)-fb,'--','linewidth',1,'color','k');
% leg = legend(datestr(t,'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% title('Allis Ridge Central Transect (A1)','FontSize',8,'FontWeight','normal');
% hXLabel = xlabel('x (m)'); hYLabel = ylabel('Ice draft (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');

% Ridge EMI data from Polona, Alli's Ridge, transect A2
clear; clc; 
load('PS_meteo.mat',"Tw_mss","Ta_T66"); % data, meteo, doi:10.18739/A2FF3M18K
load('T66.mat',"t_T66"); t_T66 = t_T66(1:1097);
t_sim = t_T66(270:1097); Tw_sim = Tw_mss(270:1097); Ta_sim = Ta_T66(270:1097);
t_sim = t_sim(155:end); Tw_sim = Tw_sim(155:end); Ta_sim = Ta_sim(155:end); % only after first EM survey

% import of EMI data
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA2_cc_consolidatedlayer.csv";
T = readtable(project); A = table2array(T); x = A(1,2:end); cl = A(2:size(A,1),2:end); time = num2str(A(2:size(A,1),1)); t = datetime(time,'InputFormat','yyyyMMdd'); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA2_cc_freeboard.csv";
T = readtable(project); A = table2array(T); fb = A(2,2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA2_cc_snow.csv";
T = readtable(project); A = table2array(T); sn = A(2:size(A,1),2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA2_cc_totalthickness.csv";
T = readtable(project); A = table2array(T); k = A(2:size(A,1),2:end); clearvars project T A

for i = 1:length(x)
    % hs_int(:,i) = interp1([datenum(t_sim(1)); datenum(t); datenum(t_sim(end))],[0; sn(:,i); 0],datenum(t_sim),'linear');
    hs_int(:,i) = interp1([datenum(t); datenum(t_sim(end))],[sn(:,i); 0],datenum(t_sim),'linear');
end

% Alli's Ridge EMI transect, fresh ice model
por = 0.29; % ridge macroporosity (from drilling)
Ta = Ta_sim; Tf = Tw_sim; % air and water temperatures
Hia = 21; % heat transfer coefficient (function of wind speed)
hs = hs_int; % snow thickness
ks = 0.26; % snow thermal conductivity, Macfarlane et al., doi:10.5194/tc-17-5417-2023
ki = 2.2; rhoi = 917; Li = 333400; % fresh ice thermodynamic parameters
% dc_0 = max(fb,0.01); % initial ice thickness = initial ice freeboard from EMI
dc_0 = cl(1,:); % initial ice thickness = initial CL thickness from EMI
td = datenum(t_sim)-datenum(t_sim(1)); t_sec = (td)*24*3600; % time
% fresh ice model
dt = diff(t_sec); n = length(t_sec); [dc,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc_em = hs_int*0;
for j = 1:length(x)
    dc(1) = dc_0(j);
    for i = 1:n-1
        R1(i) = 1./Hia; R2(i) = hs(i,j)/ks; R3(i) = dc(i)/ki;
        Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
        ddc(i) = -ki/rhoi/(Li*por)*(Tsi(i) - Tf(i))/dc(i)*dt(i);
        dc(i+1) = dc(i) + ddc(i);
    end
    % dc = dc - dc_0(j); % CL thickness = ice thickness - initial freeboard
    dc = dc - fb(j); % CL thickness = ice thickness - initial freeboard
    dc_em(:,j) = dc(:);
end
for i = 1:length(t); [~,t_em(i)] = min(abs(datenum(t_sim)-datenum(t(i)))); end
c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors

% figure
% for i = 1:length(t)
%     plot(x,dc_em(t_em(i),:),'color',c{i}); hold on
% end
% for i = 1:length(t)
%     plot(x,cl(i,:)-fb,':','linewidth',3.5,'color',c{i});
% end
% plot(x,mean(k(1,:),1)-fb,'-.','linewidth',1,'color','k');
% plot(x,-mean(sn,1)-fb,'--','linewidth',1,'color','k');
% leg = legend(datestr(t,'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% title('Allis Ridge Central Transect (A2)','FontSize',8,'FontWeight','normal');
% hXLabel = xlabel('x (m)'); hYLabel = ylabel('Ice draft (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');

% nexttile
% for i = 2:length(t)
%     p = plot(cl(i,:)-cl(1,:),dc_em(t_em(i),:)-dc_em(t_em(1),:),'o','color',c{i}); p.MarkerSize = 3.0; set(p,'markerfacecolor',get(p,'color')); hold on
% end
% plot([0 max(cl(end,:)-cl(1,:))],[0 max(cl(end,:)-cl(1,:))],'k--');
% leg = legend(datestr(t(2:end),'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% hXLabel = xlabel('EMI consolidation (m)'); hYLabel = ylabel('Model consolidation (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
% title('A2','FontSize',8,'FontWeight','normal');

nexttile
for i = 1:length(t)
    di_em_avg(i) = mean(cl(i,:)-fb(1,:)); di_em_std(i) = std(cl(i,:));
end
for i = 1:length(t_sim)
    di_mdl_avg(i) = mean(dc_em(i,:)); di_mdl_std(i) = std(dc_em(i,:));
end
p = fill([t_sim; flipud(t_sim)],[di_mdl_avg+di_mdl_std fliplr(di_mdl_avg-di_mdl_std)],1,'FaceColor',c{1},'edgecolor','none'); set(p,'facealpha',0.1); hold on
errorbar(t,di_em_avg,di_em_std,'color',c{4}); hold on
p = plot(t,di_em_avg,'o','color',c{4}); p.MarkerSize = 3; set(p,'markerfacecolor',get(p,'color'));
plot(t_sim,di_mdl_avg,'color',c{1});
set(gca,'FontSize',8,'FontWeight','normal'); ylim([1 5]);
title('A2','FontSize',8,'FontWeight','normal');
t_start = datetime('01-Jan-2020'); t_end = datetime('01-Aug-2020'); xlim([t_start t_end]); datetick('x','mmm','keeplimits'); xtickangle(0); set(gca,'YDir','reverse');
% leg = legend('','EMI','','Model','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','northeast'); leg.ItemTokenSize = [30*0.66,18*0.66];

% Ridge EMI data from Polona, Alli's Ridge, transect A3
clear; clc;
load('PS_meteo.mat',"Tw_mss","Ta_T66"); % data, meteo, doi:10.18739/A2FF3M18K
load('T66.mat',"t_T66"); t_T66 = t_T66(1:1097);
t_sim = t_T66(270:1097); Tw_sim = Tw_mss(270:1097); Ta_sim = Ta_T66(270:1097);
t_sim = t_sim(155:end); Tw_sim = Tw_sim(155:end); Ta_sim = Ta_sim(155:end); % only after first EM survey

% import of EMI data
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA3_cc_consolidatedlayer.csv";
T = readtable(project); A = table2array(T); x = A(1,2:end); cl = A(2:size(A,1),2:end); time = num2str(A(2:size(A,1),1)); t = datetime(time,'InputFormat','yyyyMMdd'); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA3_cc_freeboard.csv";
T = readtable(project); A = table2array(T); fb = A(2,2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA3_cc_snow.csv";
T = readtable(project); A = table2array(T); sn = A(2:size(A,1),2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeA3_cc_totalthickness.csv";
T = readtable(project); A = table2array(T); k = A(2:size(A,1),2:end); clearvars project T A

for i = 1:length(x)
    % hs_int(:,i) = interp1([datenum(t_sim(1)); datenum(t); datenum(t_sim(end))],[0; sn(:,i); 0],datenum(t_sim),'linear');
    hs_int(:,i) = interp1([datenum(t); datenum(t_sim(end))],[sn(:,i); 0],datenum(t_sim),'linear');
end

% Alli's Ridge EMI transect, fresh ice model
por = 0.29; % ridge macroporosity (from drilling)
Ta = Ta_sim; Tf = Tw_sim; % air and water temperatures
Hia = 21; % heat transfer coefficient (function of wind speed)
hs = hs_int; % snow thickness
ks = 0.26; % snow thermal conductivity, Macfarlane et al., doi:10.5194/tc-17-5417-2023
ki = 2.2; rhoi = 917; Li = 333400; % fresh ice thermodynamic parameters
% dc_0 = max(fb,0.01); % initial ice thickness = initial ice freeboard from EMI
dc_0 = cl(1,:); % initial ice thickness = initial CL thickness from EMI
td = datenum(t_sim)-datenum(t_sim(1)); t_sec = (td)*24*3600; % time
% fresh ice model
dt = diff(t_sec); n = length(t_sec); [dc,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc_em = hs_int*0;
for j = 1:length(x)
    dc(1) = dc_0(j);
    for i = 1:n-1
        R1(i) = 1./Hia; R2(i) = hs(i,j)/ks; R3(i) = dc(i)/ki;
        Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
        ddc(i) = -ki/rhoi/(Li*por)*(Tsi(i) - Tf(i))/dc(i)*dt(i);
        dc(i+1) = dc(i) + ddc(i);
    end
    % dc = dc - dc_0(j); % CL thickness = ice thickness - initial freeboard
    dc = dc - fb(j); % CL thickness = ice thickness - initial freeboard
    dc_em(:,j) = dc(:);
end
c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors
for i = 1:length(t); [~,t_em(i)] = min(abs(datenum(t_sim)-datenum(t(i)))); end

nexttile
for i = 1:length(t)
    di_em_avg(i) = mean(cl(i,:)-fb(1,:)); di_em_std(i) = std(cl(i,:));
end
for i = 1:length(t_sim)
    di_mdl_avg(i) = mean(dc_em(i,:)); di_mdl_std(i) = std(dc_em(i,:));
end
p = fill([t_sim; flipud(t_sim)],[di_mdl_avg+di_mdl_std fliplr(di_mdl_avg-di_mdl_std)],1,'FaceColor',c{1},'edgecolor','none'); set(p,'facealpha',0.1); hold on
errorbar(t,di_em_avg,di_em_std,'color',c{4}); hold on
p = plot(t,di_em_avg,'o','color',c{4}); p.MarkerSize = 3; set(p,'markerfacecolor',get(p,'color'));
plot(t_sim,di_mdl_avg,'color',c{1});
plot(t_sim(1),di_mdl_avg(1)*0,'--','color',c{3});
set(gca,'FontSize',8,'FontWeight','normal'); ylim([1 5]);
title('A3','FontSize',8,'FontWeight','normal');
t_start = datetime('01-Jan-2020'); t_end = datetime('01-Aug-2020'); xlim([t_start t_end]); datetick('x','mmm','keeplimits'); xtickangle(0); set(gca,'YDir','reverse');
% leg = legend('','EMI','','Model','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','northeast'); leg.ItemTokenSize = [30*0.66,18*0.66];
leg = legend('','EMI','','Model','IMB','box','off','NumColumns',3); set(leg,'FontSize',7,'Location','southoutside'); leg.ItemTokenSize = [30*0.66,18*0.66];

% nexttile
% for i = 2:length(t)
%     p = plot(cl(i,:)-cl(1,:),dc_em(t_em(i),:)-dc_em(t_em(1),:),'o','color',c{i}); p.MarkerSize = 3.0; set(p,'markerfacecolor',get(p,'color')); hold on
% end
% plot([0 max(cl(end,:)-cl(1,:))],[0 max(cl(end,:)-cl(1,:))],'k--');
% leg = legend(datestr(t(2:end),'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% hXLabel = xlabel('EMI consolidation (m)'); hYLabel = ylabel('Model consolidation (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
% title('A3','FontSize',8,'FontWeight','normal');

% figure
% c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors
% for i = 1:length(t); [~,t_em(i)] = min(abs(datenum(t_sim)-datenum(t(i)))); end
% for i = 1:length(t)
%     plot(x,dc_em(t_em(i),:),'color',c{i}); hold on
% end
% for i = 1:length(t)
%     plot(x,cl(i,:)-fb,':','linewidth',3.5,'color',c{i});
% end
% plot(x,mean(k(1,:),1)-fb,'-.','linewidth',1,'color','k');
% plot(x,-mean(sn,1)-fb,'--','linewidth',1,'color','k');
% leg = legend(datestr(t,'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% title('Allis Ridge Central Transect (A2)','FontSize',8,'FontWeight','normal');
% hXLabel = xlabel('x (m)'); hYLabel = ylabel('Ice draft (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');

% Ridge EMI data from Polona, Fort Ridge, transect 1
clear; clc;
load('PS_meteo.mat',"Tw_mss","Ta_T66"); % data, meteo, doi:10.18739/A2FF3M18K
load('T66.mat',"t_T66"); t_T66 = t_T66(1:1097);
t_sim = t_T66(286:1097); Tw_sim = Tw_mss(286:1097); Ta_sim = Ta_T66(286:1097); % only after first EM survey

% import of EMI data
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeFR1_cc_consolidatedlayer.csv";
T = readtable(project); A = table2array(T); x = A(1,2:end); cl = A(2:size(A,1),2:end); time = num2str(A(2:size(A,1),1)); t = datetime(time,'InputFormat','yyyyMMdd'); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeFR1_cc_freeboard.csv";
T = readtable(project); A = table2array(T); fb = A(2,2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeFR1_cc_snow.csv";
T = readtable(project); A = table2array(T); sn = A(2:size(A,1),2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeFR1_cc_totalthickness.csv";
T = readtable(project); A = table2array(T); k = A(2:size(A,1),2:end); clearvars project T A

for i = 1:length(x)
    hs_int(:,i) = interp1([datenum(t); datenum(t_sim(end))],[sn(:,i); 0],datenum(t_sim),'linear');
end

% Fort Ridge EMI transect, fresh ice model
por = 0.295; % ridge macroporosity (from drilling)
Ta = Ta_sim; Tf = Tw_sim; % air and water temperatures
Hia = 21; % heat transfer coefficient (function of wind speed)
hs = hs_int; % snow thickness
ks = 0.26; % snow thermal conductivity, Macfarlane et al., doi:10.5194/tc-17-5417-2023
ki = 2.2; rhoi = 917; Li = 333400; % fresh ice thermodynamic parameters
% dc_0 = max(fb,0.01); % initial ice thickness = initial ice freeboard from EMI
dc_0 = cl(1,:); % initial ice thickness = initial CL thickness from EMI
td = datenum(t_sim)-datenum(t_sim(1)); t_sec = (td)*24*3600; % time
% fresh ice model
dt = diff(t_sec); n = length(t_sec); [dc,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc_em = hs_int*0;
for j = 1:length(x)
    dc(1) = dc_0(j);
    for i = 1:n-1
        R1(i) = 1./Hia; R2(i) = hs(i,j)/ks; R3(i) = dc(i)/ki;
        Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
        ddc(i) = -ki/rhoi/(Li*por)*(Tsi(i) - Tf(i))/dc(i)*dt(i);
        dc(i+1) = dc(i) + ddc(i);
    end
    % dc = dc - dc_0(j); % CL thickness = ice thickness - initial freeboard
    dc = dc - fb(j); % CL thickness = ice thickness - initial freeboard
    dc_em(:,j) = dc(:);
end
for i = 1:length(t); [~,t_em(i)] = min(abs(datenum(t_sim)-datenum(t(i)))); end
c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors

nexttile
for i = 1:length(t)
    di_em_avg(i) = mean(cl(i,:)-fb(1,:)); di_em_std(i) = std(cl(i,:));
end
for i = 1:length(t_sim)
    di_mdl_avg(i) = mean(dc_em(i,:)); di_mdl_std(i) = std(dc_em(i,:));
end
p = fill([t_sim; flipud(t_sim)],[di_mdl_avg+di_mdl_std fliplr(di_mdl_avg-di_mdl_std)],1,'FaceColor',c{1},'edgecolor','none'); set(p,'facealpha',0.1); hold on
errorbar(t,di_em_avg,di_em_std,'color',c{4}); hold on
p = plot(t,di_em_avg,'o','color',c{4}); p.MarkerSize = 3; set(p,'markerfacecolor',get(p,'color'));
plot(t_sim,di_mdl_avg,'color',c{1});

load('DTC_data.mat',"t_DTC","fb_DTC","hs_DTC","hi_DTC"); load('T60.mat',"t_T60");
t_hc =  [  1   9  50 100 150 200 250 275 285 290 310 315 325 350 375 395 405 410 420];
hc_T = -[147 147 177 200 214 233 256 280 319 323 354 374 412 429 448 450 450 450 450]/100;
hc_T_int = interp1(datenum(t_T60(t_hc)),hc_T,datenum(t_T60),'linear');
i = 12; plot(datetime(t_DTC{i},'ConvertFrom','datenum'),-hi_DTC{i},'-','linewidth',0.5,'Color',c{3}); hold on % DTC 25
plot(t_T60,-hc_T_int-0.5,'-.','linewidth',0.5,'Color',c{3}); % T60
i = 11; plot(datetime(t_DTC{i},'ConvertFrom','datenum'),-hi_DTC{i},'--','linewidth',0.5,'Color',c{3}); % DTC 24

set(gca,'FontSize',8,'FontWeight','normal'); ylim([0.5 4.0]);
title('FR1','FontSize',8,'FontWeight','normal');
t_start = datetime('01-Jan-2020'); t_end = datetime('01-Apr-2020'); xlim([t_start t_end]); datetick('x','mmm','keeplimits'); xtickangle(0); set(gca,'YDir','reverse');
% leg = legend('','EMI','','Model','IMB','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','southwest'); leg.ItemTokenSize = [30*0.66,18*0.66];

% nexttile
% for i = 2:length(t)
%     p = plot(cl(i,:)-cl(1,:),dc_em(t_em(i),:)-dc_em(t_em(1),:),'o','color',c{i}); p.MarkerSize = 3.0; set(p,'markerfacecolor',get(p,'color')); hold on
% end
% plot([0 max(cl(end,:)-cl(1,:))],[0 max(cl(end,:)-cl(1,:))],'k--');
% leg = legend(datestr(t(2:end),'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','northwest'); leg.ItemTokenSize = [30*0.5,18*0.5];
% hXLabel = xlabel('EMI consolidation (m)'); hYLabel = ylabel('Model consolidation (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
% title('FR1','FontSize',8,'FontWeight','normal');

% figure
% for i = 1:length(t)
%     plot(x,dc_em(t_em(i),:),'color',c{i}); hold on
% end
% for i = 1:length(t)
%     plot(x,cl(i,:)-fb,':','linewidth',3.5,'color',c{i});
% end
% plot(x,mean(k(1:length(t),:),1)-fb,'-.','linewidth',1,'color','k');
% plot(x,-mean(sn,1)-fb,'--','linewidth',1,'color','k');
% leg = legend(datestr(t,'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% title('Allis Ridge Central Transect (A2)','FontSize',8,'FontWeight','normal');
% hXLabel = xlabel('x (m)'); hYLabel = ylabel('Ice draft (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');

% Ridge EMI data from Polona, Fort Ridge, transect 2
clear; clc;
load('PS_meteo.mat',"Tw_mss","Ta_T66"); % data, meteo, doi:10.18739/A2FF3M18K
load('T66.mat',"t_T66"); t_T66 = t_T66(1:1097);
t_sim = t_T66(294:1097); Tw_sim = Tw_mss(294:1097); Ta_sim = Ta_T66(294:1097); % only after first EM survey

% import of EMI data
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeFR2_cc_consolidatedlayer.csv";
T = readtable(project); A = table2array(T); x = A(1,2:end); cl = A(2:size(A,1),2:end); time = num2str(A(2:size(A,1),1)); t = datetime(time,'InputFormat','yyyyMMdd'); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeFR2_cc_freeboard.csv";
T = readtable(project); A = table2array(T); fb = A(2,2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeFR2_cc_snow.csv";
T = readtable(project); A = table2array(T); sn = A(2:size(A,1),2:end); clearvars project T A
project = "C:\Users\evgenii.salganik\Documents\MATLAB\datasets\ridge_model_input_EM\ridgeFR2_cc_totalthickness.csv";
T = readtable(project); A = table2array(T); k = A(2:size(A,1),2:end); clearvars project T A

for i = 1:length(x)
    hs_int(:,i) = interp1([datenum(t); datenum(t_sim(end))],[sn(:,i); 0],datenum(t_sim),'linear');
end

% Fort Ridge EMI transect, fresh ice model
por = 0.295; % ridge macroporosity (from drilling)
Ta = Ta_sim; Tf = Tw_sim; % air and water temperatures
Hia = 21; % heat transfer coefficient (function of wind speed)
hs = hs_int; % snow thickness
ks = 0.26; % snow thermal conductivity, Macfarlane et al., doi:10.5194/tc-17-5417-2023
ki = 2.2; rhoi = 917; Li = 333400; % fresh ice thermodynamic parameters
% dc_0 = max(fb,0.01); % initial ice thickness = initial ice freeboard from EMI
dc_0 = cl(1,:); % initial ice thickness = initial CL thickness from EMI
td = datenum(t_sim)-datenum(t_sim(1)); t_sec = (td)*24*3600; % time
% fresh ice model
dt = diff(t_sec); n = length(t_sec); [dc,R1,R2,R3,ddc,Tsi] = deal(zeros(1,n)); dc_em = hs_int*0;
for j = 1:length(x)
    dc(1) = dc_0(j);
    for i = 1:n-1
        R1(i) = 1./Hia; R2(i) = hs(i,j)/ks; R3(i) = dc(i)/ki;
        Tsi(i) = (Ta(i) - Tf(i))*R3(i)./(R1(i)+R2(i)+R3(i)) + Tf(i);
        ddc(i) = -ki/rhoi/(Li*por)*(Tsi(i) - Tf(i))/dc(i)*dt(i);
        dc(i+1) = dc(i) + ddc(i);
    end
    % dc = dc - dc_0(j); % CL thickness = ice thickness - initial freeboard
    dc = dc - fb(j); % CL thickness = ice thickness - initial freeboard
    dc_em(:,j) = dc(:);
end
for i = 1:length(t); [~,t_em(i)] = min(abs(datenum(t_sim)-datenum(t(i)))); end
c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors

nexttile
for i = 1:length(t)
    di_em_avg(i) = mean(cl(i,:)-fb(1,:)); di_em_std(i) = std(cl(i,:));
end
for i = 1:length(t_sim)
    di_mdl_avg(i) = mean(dc_em(i,:)); di_mdl_std(i) = std(dc_em(i,:));
end
p = fill([t_sim; flipud(t_sim)],[di_mdl_avg+di_mdl_std fliplr(di_mdl_avg-di_mdl_std)],1,'FaceColor',c{1},'edgecolor','none'); set(p,'facealpha',0.1); hold on
errorbar(t,di_em_avg,di_em_std,'color',c{4}); hold on
p = plot(t,di_em_avg,'o','color',c{4}); p.MarkerSize = 3; set(p,'markerfacecolor',get(p,'color'));
plot(t_sim,di_mdl_avg,'color',c{1});
set(gca,'FontSize',8,'FontWeight','normal'); ylim([0.5 4.0]);
title('FR2','FontSize',8,'FontWeight','normal');
t_start = datetime('01-Jan-2020'); t_end = datetime('01-Apr-2020'); xlim([t_start t_end]); datetick('x','mmm','keeplimits'); xtickangle(0); set(gca,'YDir','reverse');
% leg = legend('','EMI','','Model','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','southeast'); leg.ItemTokenSize = [30*0.66,18*0.66];

annotation('textbox',[0.005 .51 0.01 .51],'String','(a)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.14  .51 0.15 .51],'String','(b)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.27  .51 0.28 .51],'String','(c)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.40  .51 0.41 .51],'String','(d)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.535 .51 0.54 .51],'String','(e)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');

% nexttile
% for i = 2:length(t)
%     p = plot(cl(i,:)-cl(1,:),dc_em(t_em(i),:)-dc_em(t_em(1),:),'o','color',c{i}); p.MarkerSize = 3.0; set(p,'markerfacecolor',get(p,'color')); hold on
% end
% plot([0 max(cl(end,:)-cl(1,:))],[0 max(cl(end,:)-cl(1,:))],'k--');
% leg = legend(datestr(t(2:end),'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% hXLabel = xlabel('EMI consolidation (m)'); hYLabel = ylabel('Model consolidation (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
% title('FR2','FontSize',8,'FontWeight','normal');

% figure
% c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors
% for i = 1:length(t)
%     plot(x,dc_em(t_em(i),:),'color',c{i}); hold on
% end
% for i = 1:length(t)
%     plot(x,cl(i,:)-fb,':','linewidth',3.5,'color',c{i});
% end
% plot(x,mean(k(1:3,:),1)-fb,'-.','linewidth',1,'color','k');
% plot(x,-mean(sn,1)-fb,'--','linewidth',1,'color','k');
% leg = legend(datestr(t,'dd mmm'),'box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
% title('Allis Ridge Central Transect (A2)','FontSize',8,'FontWeight','normal');
% hXLabel = xlabel('x (m)'); hYLabel = ylabel('Ice draft (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');