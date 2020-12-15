function [avgACCW,CCE,ACE] = calcACE(AP,mtlmass)


%% Define variables

locations = ["Alaska";                 % Headings for scalingFactor
             "Washington";
             "Northern Oregon";
             "Oregon";
             "Northern California";
             "Southern California";
             "Hawaii"];

scalingFactor = [0.2430 0.3320 0.0750 0.2000 0.0240 0.0120;  % These are given in
                 0.1370 0.2770 0.0410 0.3380 0.0220 0.0450;  % the WE Prize doc
                 0.1550 0.3070 0.0560 0.3440 0.0370 0.0420;
                 0.1750 0.2680 0.0580 0.2950 0.0340 0.0540;
                 0.2070 0.2300 0.0120 0.4660 0.0160 0.0640;
                 0.1520 0.2700 0.0140 0.3910 0.0100 0.0950;
                 0.3280 0.2450 0.0010 0.1330 0.0000 0.0130];

CP = [35.5; 32.7; 39.3; 37.9; 31.5; 31.2; 16.8];

MMCLow  = [2250 424 6000 7200 4900 7500 4630];
MMCMed  = [3000 510 7900 9500 5900 8200 5510];
MMCHigh = [4500 557 12000 13500 8000 9500 6620];

%% Calculate ACCW
% prompt_AP = "Enter a vector with average "; 
% prompt_AP = prompt_AP + newline + "power for sea states 1-6: ";

% AP = input(prompt_AP);   % Input a vector for AP, average power  
                         % absorbed by the WEC for each sea state   
tposecheck = isrow(AP);
if tposecheck == 1
   AP = AP.';
end
                         
ACCW = (scalingFactor*AP)./CP;   % From WE Prize doc
avgACCW = mean(ACCW);

%% Calculate CCE 
% prompt = "Enter mass (metric tonne) of [A36 steel, ";  
% prompt = prompt + newline + "concrete, HDPE, coated fabric, 5083 aluminum, " ;
% prompt = prompt + newline + "fiberglass, filament wound fiberglass] as a vector: ";
% 
% disp(" ");
% mtlmass = input(prompt); % Input masses of each material in WEC

tposecheck2 = isrow(mtlmass);
if tposecheck == 0
   mtlmass = mtlmass.';
end

CCE_low  = dot(MMCLow,mtlmass);   % material mass in kg * cost per kg
CCE_med  = dot(MMCMed,mtlmass);
CCE_high = dot(MMCHigh,mtlmass);
CCE_dollars = [CCE_low CCE_med CCE_high];

CCE_low_mill  = CCE_low./1000000;
CCE_med_mill  = CCE_med./1000000;
CCE_high_mill = CCE_high./1000000;
CCE = [CCE_low_mill CCE_med_mill CCE_high_mill];
out = " ";
out = out + newline + newline + "      CCE_low     CCE_med     CCE_high  ($M)";
% disp(out)
CCE;

%% Calculate ACE

ACE_low  = avgACCW./CCE_high_mill;
ACE_med  = avgACCW./CCE_med_mill;
ACE_high = avgACCW./CCE_low_mill;
ACE = [ACE_low ACE_med ACE_high];

out2 = " ";
out2 = out2 + newline;
out2 = out2 + newline + "      ACE_low     ACE_med     ACE_high  (m/$M)";
% disp(out2)
ACE;

end







