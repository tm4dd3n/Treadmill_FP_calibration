%% This uses ONLY FORCES to find the calibration matrix
%% Creates 12X12 calibration matrix to be applied to ANALOG forces
%THIS CODE does not rearrange LOAD CELL DATA TO 1XYZ 2XYZ 3XYZ 4XYZ 
%Maintains original order of analog data

%%
clear
clc

data_folder = '/Users/meganpeach/Documents/MATLAB/Processed4';
num_trials = 12;
slope = 50.1779; %from calwand load cell calibration for volts to lbs (see bottom of code)
% fpo_lab = [.782 .305 .015];  %fp origin in lab cs from cortex BSI02
% fpo_lab = [.568 .304 .0112];  %fp origin in lab cs from cortex BSI03
fpo_lab = [.769 .308 .0044];  %fp origin in lab cs from cortex BSI04
r_fp_to_lab = [0 1 0; 1 0 0; 0 0 -1];  %rotation matrix TM to LAB CS (from Cortex), X is AP, Y is ML

%% Read in raw voltage data from treadmill load cells and calwand

%build the file path
%reads in .anc excel file for each #1-12 trial into "analog" structure
for i = 1:num_trials
    disp(['Processing trial number: ', num2str(i)]);
    file_name = sprintf('Trimmed_PILS%02d.anc.xlsx',i); %formats filename for each iteration
    file_path = fullfile(data_folder, file_name);

    analog(i).data = readmatrix(file_path); %read data
    analog(i).data = analog(i).data(4:4:end,:); %downsample
end
num_frames = length(analog(1).data);
    clear file_name file_path i

    %% Pull out load cell and calwand data
        % original .anc output: F1Y1-TM F1X1-TM F1Z1-TM	F1Y4-TM	F1X4-TM	F1Z4-TM	 F1Y2-TM F1X2-TM F1Z2-TM	F1Y3-TM	F1X3-TM	F1Z3-TM
        % Column 1 is time
        % Column 2-13 are F1Y F1X F1Z.....etc
        % Column 14 is load cell analog data

%conversion factor must be applied prior to rearranging to XYZ

conversion = [73.8367; 73.7637; 315.8560; 73.8127; 73.9503; 317.6015; 73.5251; 73.4808; 315.8161; 73.7191; 73.7300; 316.9170];

for i = 1:num_trials
    lcdata(i).analog = analog(i).data(:,2:13); %ADC counts = raw digital units
    lcdata(i).analog_scaled = lcdata(i).analog * .000152587890625;  %to scale to VOLTS from Cortex scaling factor, YXZ and 1423 order
    %DOES THE CONVERSION FACTOR ROTATE TO LAB CS?
    lcdata(i).newtonYXZ = lcdata(i).analog_scaled .* conversion';  %scales from VOLTS to NEWTONS using conversion (Steve), 1423 and YXZ order, CS?
    
    caldata(i).analog = analog(i).data(:,14); %ADC counts, FP CS
    caldata(i).volts = caldata(i).analog * .000152587890625;  % ADC to volts, FP CS
    caldata(i).lb = caldata(i).volts * slope;  % Volts to lbs, FP CS
    caldata(i).newton = caldata(i).lb * 4.4482216; %lbs to newtons, FP CS
end

%% Define individual load cell analog data
%CS?
% for i = 1:num_trials
%     lcdata(i).lc1_newton = [lcdata(i).newtonYXZ(:,1:3)]'; 
%     lcdata(i).lc2_newton = [lcdata(i).newtonYXZ(:,7:9)]';
%     lcdata(i).lc3_newton = [lcdata(i).newtonYXZ(:,10:12)]';
%     lcdata(i).lc4_newton = [lcdata(i).newtonYXZ(:,4:6)]';
% 
%     % lcdata(i).newtonXYZ = [lcdata(i).lc1_newton lcdata(i).lc2_newton lcdata(i).lc3_newton lcdata(i).lc4_newton]; %units = newton, 1234 and XYZ order, %?? CS
% end

%% Call in markers for calwand
%build the file path for calwand .trc 
% LAB CS
for i = 1:num_trials
    disp(['Processing trial number: ', num2str(i)]);
    file_name = sprintf('Trimmed_PILS%02d-FP_Cal_Wand.trc.xlsx',i); %formats filename for each iteration
    file_path = fullfile(data_folder, file_name);
       %disp(['Reading file: ', file_path]);

    caldata(i).trc = readmatrix(file_path); %read in trc data
    caldata(i).trc = caldata(i).trc(4:end,:);  
    caldata(i).trc_m = caldata(i).trc/1000; %unit conversion from mm to m, LAB CS
end

%check calwand tip positions on treadmill
for i = 1:num_trials
    trc = caldata(i).trc_m;
    x=trc(:,24);
    y=trc(:,25);
scatter(x,y,10,'filled');
hold on
end
clear x y trc i file_name file_path trc

%% Calculate vector from tip to top as a unit vector
%V_top_mid X Y Z = row 18,19,20
%V_tip X Y Z = row 24, 25, 26
for i = 1:num_trials
caldata(i).tip = caldata(i).trc_m(:,24:26); %calwand tip marker coordinates
caldata(i).tiptotopvector = caldata(i).trc_m(:,18:20) - caldata(i).trc_m(:,24:26); %top - tip

rowNorms = vecnorm(caldata(i).tiptotopvector, 2, 2);     % Normalize the rows first: N×1 vector
caldata(i).vecnormed = caldata(i).tiptotopvector ./ rowNorms;    % Normalize vector per row, LAB CS
end

% Multiply normed vector with force vector to get calwand XYZ vector
for i = 1:num_trials
    for j = 1:num_frames
    force = caldata(i).newton; %FP CS
    dir = (r_fp_to_lab * caldata(i).vecnormed(j,:)')'; %Rotates normed vector to force place CS
    caldata(i).vector(j,:) = force(j)*dir;   %newtons, FPCS
    caldata(i).vector_lab(j,:) = (caldata(i).vector(j,:) * r_fp_to_lab); % Rotates Calwand Vector to Lab CS
    end
end

clear i j force dir rowNorms

% Can check normalized tip to top vector with: norm(normalized vector) = 1

%% Find offset of calwand in LAB CS with relation to center of FP (in LAB CS) then find CALWAND MOMENTS

for i = 1:num_trials
    caldata(i).offset =  caldata(i).tip - fpo_lab; % LAB CS
end
clear i

% Find moments from calwand forces and offset from FP origin
% Moment (calwand) = COP(calwand) x force(calwand) = (COP(calwand) - FPC)) x force(calwand)  THIS IS A CROSS-PRODUCT!
for i = 1:num_trials
    caldata(i).moments_lab = (cross( caldata(i).offset, caldata(i).vector_lab));  %in Nm, LAB CS  
end
clear i 

%% Load cell positions in LAB CS 

% Call in markers (LAB CS) for treadmill
for i = 1:num_trials
    disp(['Processing trial number: ', num2str(i)]);
    file_name = sprintf('Trimmed_PILS%02d-Treadmill.trc.xlsx',i); %formats filename for each iteration
    file_path = fullfile(data_folder, file_name);
       %disp(['Reading file: ', file_path]);

    treadmill(i).trc = readmatrix(file_path); %read in trc data
    treadmill(i).trc = treadmill(i).trc(4:end,:);  
    treadmill(i).trc_m = treadmill(i).trc/1000; %unit conversion from mm to m, LAB CS

    lc1_lab = mean(treadmill(i).trc_m(:,3:5));
    lc2_lab = mean(treadmill(i).trc_m(:,6:8));
    lc3_lab = mean(treadmill(i).trc_m(:,9:11));
    lc4_lab = mean(treadmill(i).trc_m(:,12:14));
end

% 3x4 matrix for four load cell positions with rows = XYZ and columns = 1234
% LAB CS
lc_positions_lab = [lc1_lab; lc2_lab; lc3_lab; lc4_lab]';
clear lc1_lab lc2_lab lc3_lab lc4_lab
%% Decompose calibration wand forces in to treadmill load cell forces

for i = 1:num_trials
    wand_pos = (mean(caldata(i).tip))'; %LAB CS
    wand_force = caldata(i).vector_lab'; %LAB CS

    % find distance from wand to each load cell
    % vecnorm computes Euclidian norm giving a 1x4 matrix of distances from calwand to each load cell
    distances = vecnorm(lc_positions_lab - wand_pos, 2, 1);
    inv_dist = 1 ./distances; %weights the distances
    weights = inv_dist / sum(inv_dist); %normalize the weights to sum to 1

    % assign scaled forces to each load cell
    % for lc = 1:4 %four load cells
    %     idx = (lc-1)*3 + [1:3]; %channel index for LCx LCy LCz
        for j = 1:num_frames
        lc_temp1 = (weights(1) * wand_force(:,j))'; %load cell 1
        lc_temp2 = (weights(2) * wand_force(:,j))'; %load cell 2, etc...
        lc_temp3 = (weights(3) * wand_force(:,j))';
        lc_temp4 = (weights(4) * wand_force(:,j))';
        caldata(i).lcexpected(j,:) = [lc_temp1 lc_temp4 lc_temp2 lc_temp3]; %LAB CS  1423 order
        end
end

clear lc idx S_temp distances inv_dist weights wand_pos wand_force i j lc_temp4 lc_temp3 lc_temp2 lc_temp1

%% Plot raw load cell data (N) and raw deconstructed wand data (N)
% figure
% subplot(3,1,1)          % 3 rows, 1 column, 1st plot     
% plot(lcdata(6).newtonYXZ(:,3:3:end),"Color",'r')
% hold on
% plot(caldata(6).lcexpected(:,3:3:end),"Color",'g')

%% Calculate calibation(C) matrix 12x12
%expected = caldata.lcexpected -> decomposed calwand forces to each load cell
%raw = lcdata.newtonYXZ
s_raw = (vertcat(lcdata(:).newtonYXZ))'; %1423 LC order
s_expected = (vertcat(caldata(:).lcexpected))';

%Equation: expected = C * raw
% C = 12x12 global calibration matrix
C = s_expected * pinv(s_raw);

%% Calibrate raw analog data
for i = 1:num_trials
    for j = 1:num_frames
    raw_temp = lcdata(i).newtonYXZ';
    calibrated(i).raw(j,:) = C * raw_temp(:,j); %newtons, YXZ and 1423 order, CS???
    end
end

%% RESIDUAL ERROR FROM KNOWN TO MEASURED FORCE BEFORE CALIBRATION (12 columns, YXZ format)
for i=1:num_trials
    for j = 1:num_frames
    force_calwand = caldata(i).lcexpected(j,3:3:end);
    force_measured = lcdata(i).newtonYXZ(j,3:3:end);
    error(i).forcez_original = sqrt(mean((force_calwand - force_measured).^2));

    force_calwand = caldata(i).lcexpected(j,3:3:end);
    force_calibrated = calibrated(i).raw(j,3:3:end);
    error(i).forcez_cal = sqrt(mean((force_calwand - force_calibrated).^2));

    force_calwand = caldata(i).lcexpected(j,2:3:end);
    force_measured = lcdata(i).newtonYXZ(j,2:3:end);
    error(i).forcey_original = sqrt(mean((force_calwand - force_measured).^2));

    force_calwand = caldata(i).lcexpected(j,2:3:end);
    force_calibrated = calibrated(i).raw(j,2:3:end);
    error(i).forcey_cal = sqrt(mean((force_calwand - force_calibrated).^2));
  
    force_calwand = caldata(i).lcexpected(j,1:3:end);
    force_measured = lcdata(i).newtonYXZ(j,1:3:end);
    error(i).forcex_original = sqrt(mean((force_calwand - force_measured).^2));

    force_calwand = caldata(i).lcexpected(j,1:3:end);
    force_calibrated = calibrated(i).raw(j,1:3:end);
    error(i).forcex_cal = sqrt(mean((force_calwand - force_calibrated).^2));
    end
end
clear force_calwand force_calibrated i j 

%% Find composite force from load cells in XYZ BEFORE calibration
%Equations from Steve's Spreadsheet
%Fy1 in equation corresponds to light green column (EX. F1Y1)
for i = 1:num_trials
    for j = 1:num_frames
    lcdata(i).grfx(j) = lcdata(i).newtonYXZ(j,2) - lcdata(i).newtonYXZ(j,8) - lcdata(i).newtonYXZ(j,11) + lcdata(i).newtonYXZ(j,5);

    lcdata(i).grfy(j) = -lcdata(i).newtonYXZ(j,1) + lcdata(i).newtonYXZ(j,7) + lcdata(i).newtonYXZ(j,10) - lcdata(i).newtonYXZ(j,4);

    lcdata(i).grfz(j) = lcdata(i).newtonYXZ(j,3) + lcdata(i).newtonYXZ(j,9) + lcdata(i).newtonYXZ(j,12) + lcdata(i).newtonYXZ(j,6);
    
    %%%%% ADDED NEG(-) TO GRFZ
    lcdata(i).forces = [lcdata(i).grfx; lcdata(i).grfy; -lcdata(i).grfz]; %CS?
    end
end

%rotate forces to lab CS
for i = 1:num_trials
    for j = 1:num_frames
        lcdata(i).forces_lab(j,:) = r_fp_to_lab * lcdata(i).forces(:,j);  %rotate to lab CS
    end
end

%% Find composite Force from LC in XYZ AFTER calibration
%Equations from Steve's Spreadsheet
%Fy1 in equation corresponds to light green column (EX. F1Y1)

for i = 1:num_trials
    for j = 1:num_frames
    calibrated(i).grfx(j) = calibrated(i).raw(j,2) - calibrated(i).raw(j,8) - calibrated(i).raw(j,11) + calibrated(i).raw(j,5);

    calibrated(i).grfy(j) = -calibrated(i).raw(j,1) + calibrated(i).raw(j,7) + calibrated(i).raw(j,10) - calibrated(i).raw(j,4);

    calibrated(i).grfz(j) = calibrated(i).raw(j,3) + calibrated(i).raw(j,9) + calibrated(i).raw(j,12) + calibrated(i).raw(j,6);

    %%%%% ADDED NEG(-) TO GRFZ
    calibrated(i).forces = [calibrated(i).grfx; calibrated(i).grfy; -calibrated(i).grfz]; %CS??
    end
end

for i = 1:num_trials
    for j = 1:num_frames
        calibrated(i).forces_lab(j,:) = r_fp_to_lab * calibrated(i).forces(:,j);  %rotate to lab CS
    end
end


%% Compare corrected FORCES from treadmill load cells to Calwand
%red = original lc data
%black = calibrated lc data
%green = calwand forces

r = randi([1, 12], 1, 1); % to generate a random trial from 1:12

figure
subplot(3,1,1)          % 3 rows, 1 column, 1st plot     
plot(lcdata(r).forces_lab(:,1),"Color",'r')
hold on
plot(caldata(r).vector_lab(:,1),"Color",'g')
plot(calibrated(r).forces_lab(:,1),"Color", "k") %forces post-calibration

subplot(3,1,2)          % 3 rows, 1 column, 1st plot     
plot(lcdata(r).forces_lab(:,2),"Color",'r')
hold on
plot(caldata(r).vector_lab(:,2),"Color",'g')
plot(calibrated(r).forces_lab(:,2),"Color", "k") %forces post-calibration

subplot(3,1,3)          % 3 rows, 1 column, 1st plot     
plot(lcdata(r).forces_lab(:,3),"Color",'r')
hold on
plot(caldata(r).vector_lab(:,3),"Color",'g')
plot(calibrated(r).forces_lab(:,3),"Color", "k") %forces post-calibration

%% RESIDUAL ERROR FROM KNOWN TO MEASURED FORCE AFTER CALIBRATION

for i=1:num_trials
    for j = 1:num_frames
    force_calwand = caldata(i).vector_lab(:,3);
    force_calibrated = calibrated(i).forces_lab(:,3);
    error(i).forcez = sqrt(mean((force_calwand - force_calibrated(j)).^2));

    force_calwand = caldata(i).vector_lab(:,2);
    force_calibrated = calibrated(i).forces_lab(:,2);
    error(i).forcey = sqrt(mean((force_calwand - force_calibrated(j)).^2));

    force_calwand = caldata(i).vector_lab(:,1);
    force_calibrated = calibrated(i).forces_lab(:,1);
    error(i).forcex = sqrt(mean((force_calwand - force_calibrated(j)).^2));
    end
end
clear force_calwand force_calibrated i j 

%% Find moments FROM ORIGINAL/PRE-CALIBRATED LOAD CELL FORCE 
% Moment (loadcell) = r(TIP to tm origin) x COMPOSITE FORCE FROM LOAD CELLS
a = .4;
b = .5588;
az0 = -.035;

for i = 1:num_trials
    for j = 1:num_frames
    lcdata(i).momentx(j,:) = -b*(lcdata(i).newtonYXZ(j,3) + lcdata(i).newtonYXZ(j,9) - lcdata(i).newtonYXZ(j,12) - lcdata(i).newtonYXZ(j,6)) - ...
                                (az0 * (lcdata(i).newtonYXZ(j,2) + lcdata(i).newtonYXZ(j,8) + lcdata(i).newtonYXZ(j,11) + lcdata(i).newtonYXZ(j,5)));  %in Nm, LAB CS?
    lcdata(i).momenty(j,:) = -a*(-lcdata(i).newtonYXZ(j,3) + lcdata(i).newtonYXZ(j,9) + lcdata(i).newtonYXZ(j,12) - lcdata(i).newtonYXZ(j,6)) + ...
                                (az0 * (lcdata(i).newtonYXZ(j,1) + lcdata(i).newtonYXZ(j,7) + lcdata(i).newtonYXZ(j,10) + lcdata(i).newtonYXZ(j,4)));  %in Nm, LAB CS?
    lcdata(i).momentz(j,:) = b*(-lcdata(i).newtonYXZ(j,1) - lcdata(i).newtonYXZ(j,7) + lcdata(i).newtonYXZ(j,10) + lcdata(i).newtonYXZ(j,4)) + ...
                                a*(lcdata(i).newtonYXZ(j,2) - lcdata(i).newtonYXZ(j,8) - lcdata(i).newtonYXZ(j,11) + lcdata(i).newtonYXZ(j,5));
    end
end

for i = 1:num_trials
    for j = 1:num_frames
        lcdata(i).moments(:,j) = [lcdata(i).momentx(j) lcdata(i).momenty(j) lcdata(i).momentz(j)]; %combine into a single matrix
        lcdata(i).moments_lab(j,:) = r_fp_to_lab * lcdata(i).moments(:,j);  %rotate to lab CS
    end
end
clear i 

%% Find moments FROM CALIBRATED FORCE 
% Moment (loadcell) = r(TIP to tm origin) x COMPOSITE FORCE 

for i = 1:num_trials
   for j = 1:num_frames
   calibrated(i).moments_lab(j,:) = -(cross(caldata(i).offset(j,:), calibrated(i).forces(:,j)));  %in Nm, LAB CS  
   end
end
clear i j

%% Compare corrected MOMENTS from treadmill load cells to Calwand

r = randi([1, 12], 1, 1); % to generate a random trial from 1:12

figure
subplot(3,1,1)          % 3 rows, 1 column, 1st plot     
plot(calibrated(r).moments_lab(:,1), "Color","k")
hold on
plot(lcdata(r).moments_lab(:,1),"Color",'r')
plot(caldata(r).moments_lab(:,1),"Color",'g')

subplot(3,1,2)            
plot(calibrated(r).moments_lab(:,2), "Color","k")
hold on
plot(lcdata(r).moments_lab(:,2),"Color",'r')
plot(caldata(r).moments_lab(:,2),"Color",'g')

subplot(3,1,3)             
plot(calibrated(r).moments_lab(:,3), "Color","k")
hold on
plot(lcdata(r).moments_lab(:,3),"Color",'r')
plot(caldata(r).moments_lab(:,3),"Color",'g')

%% Read in .forces files from CORTEX for each PILS trial
% % To be compared to TM forces 
for i = 1:num_trials
    disp(['Processing trial number: ', num2str(i)]);
    file_name = sprintf('Trimmed_PILS%02d.forces.xlsx',i); %formats filename for each iteration
    file_path = fullfile(data_folder, file_name);
       %disp(['Reading file: ', file_path]);

    cortex(i).data = readmatrix(file_path); %read data
    cortex(i).data = cortex(i).data(1:4:end,:); %downsample 
    cortex(i).force = cortex(i).data(:,2:4);
    cortex(i).cop = cortex(i).data(:,5:7)/1000; %unit conversion mm to m
end

%check cop positions fron Cortex on treadmill
for i = 1:num_trials
    x=cortex(i).cop(:,1);
    y=cortex(i).cop(:,2);
scatter(x,y,10,'filled');
hold on
end

clear Mx My Mz i x y j r

%% Find measured COP with ORIGINAL FORCES
% COPx = -My/ Fz
% COPy = Mx/ Fz

% a = fpo_lab(1);
% b = fpo_lab(2);
% az0 = fpo_lab(3);

a = .4;
b = .5588;
az0 = -.035;

%Steve's equations:
for i = 1:num_trials
    for j = 1:num_frames
    lcdata(i).copx(j,:) = ((az0 * lcdata(i).forces(1,j) - lcdata(i).moments(2,j)) / lcdata(i).forces(3,j)) + a;
    lcdata(i).copy(j,:) = ((lcdata(i).moments(1,j) + az0 * lcdata(i).forces(2,j)) / lcdata(i).forces(3,j)) + b;
    end
end

copz = [zeros(num_frames,1)];
for i = 1:num_trials
    for j = 1:num_frames
        lcdata(i).cop = [lcdata(i).copx lcdata(i).copy copz]; %combine into a single matrix
        lcdata(i).cop_lab(j,:) = r_fp_to_lab * lcdata(i).cop(j,:)';  %rotate to lab CS
    end
end
clear copz


%% check cop positions from measured (lc) COP on treadmill PRE CALIBRATION
for i = 1:num_trials
    x=lcdata(i).cop_lab(:,1);
    y=lcdata(i).cop_lab(:,2);
scatter(x,y,10,'filled','MarkerFaceColor', 'r');
hold on
    x2=cortex(i).cop(:,1);
    y2=cortex(i).cop(:,2);
scatter(x2,y2,10,'filled','MarkerFaceCOlor','c');
hold on
    x3=caldata(i).tip(:,1);
    y3=caldata(i).tip(:,2);
scatter(x3,y3,10,'filled','MarkerFaceCOlor','g');
end

%% Find COP with CALIBRATED FORCES
% COPx = -My/ Fz
% COPy = Mx/ Fz

% a = .4;
% b = .5588;
% az0 = -.035;

a = fpo_lab(1);
b = fpo_lab(2);
az0 = fpo_lab(3);

for i = 1:num_trials
    for j = 1:num_frames
    calibrated(i).copx(j,:) = ((az0 * calibrated(i).forces_lab(j,1) - calibrated(i).moments_lab(j,2)) / calibrated(i).forces_lab(j,3)) + a;
    calibrated(i).copy(j,:) = ((calibrated(i).moments_lab(j,1) + az0 * calibrated(i).forces_lab(j,2)) / calibrated(i).forces_lab(j,3)) + b;
    end
end

copz = [zeros(num_frames,1)];
for i = 1:num_trials
    for j = 1:num_frames
        calibrated(i).cop_lab = [calibrated(i).copx calibrated(i).copy copz]; %combine into a single matrix
    end
end
clear copz

%% check COP positions + calibrated COP
figure
for i = 1:num_trials
    x=lcdata(i).cop_lab(:,1);
    y=lcdata(i).cop_lab(:,2);
scatter(x,y,10,'filled','MarkerFaceColor', 'r');
hold on
    x2=cortex(i).cop(:,1);
    y2=cortex(i).cop(:,2);
scatter(x2,y2,10,'filled','MarkerFaceCOlor','c');
hold on
    x3=caldata(i).tip(:,1);
    y3=caldata(i).tip(:,2);
scatter(x3,y3,10,'filled','MarkerFaceCOlor','g');
hold on
    x4=calibrated(i).cop_lab(:,1);
    y4=calibrated(i).cop_lab(:,2);
scatter(x4,y4,10,'filled','MarkerFaceCOlor','k');
end

%% RESIDUAL ERROR FROM KNOWN TO MEASURED COP after CALIBRATION

for i=1:num_trials
    for j = 1:num_frames
    copx_calwand = caldata(i).tip(:,1);
    copx_calibrated = calibrated(i).copx;
    error(i).copx = sqrt(mean((copx_calwand - copx_calibrated(j)).^2));

    copy_calwand = caldata(i).tip(:,2);
    copy_calibrated = calibrated(i).copy;
    error(i).copy = rmse(copy_calwand(j),copy_calibrated(j));
    end
end

clear copx_calwand copx_calibrated copy_calwand copy_calibrated i j 


%% MAXES (FOR COMPARISON)

% %original data
% for i = 1:num_trials
%     cortexforce = cortex(i).force';
%     maxforce(i).cortex = max(cortexforce(:,3));
% 
%     measuredforce = tmdata(i).force_lab(:,3);
%     maxforce(i).measured = max(measuredforce);
% 
%     calwandforce = caldata(i).vector_lab(:,3);
%     maxforce(i).calwand = max(calwandforce);
% end
% clear cortexforce measuredforce calwandforce
% % plot forces, should be a straight line
% figure
% for i = 1:num_trials
% scatter(maxforce(i).measured, maxforce(i).calwand, 'filled');
% xlabel('Measured Force');
% ylabel('Calwand (Known) Force');
% title('Force Comparison: Measured vs Calwand');
% grid on;
% hold on
% end
% 
% 
% %corrected data
% for i = 1:num_trials
%     correctedforce = corrected(i).force_lab';
%     maxforce(i).measured = max(correctedforce(:,3));
% 
%     calwandforce = caldata(i).vector_lab(:,3);
%     maxforce(i).calwand = max(calwandforce);
% end
% % plot forces, should be a straight line
% figure
% for i = 1:num_trials
% scatter(maxforce(i).measured, maxforce(i).calwand, 'filled');
% xlabel('Measured Force');
% ylabel('Calwand (Known) Force');
% title('Force Comparison: Measured vs Calwand');
% grid on;
% hold on
% end
% clear cortexforce measuredforce calwandforce