% Functions and script developed by Chaitanya Mandala, Christina Reid,
% Francisco Pallares Chamorro, Giulia Ciabatti, Laura De Leo and Roman Horskov
%
% This script will pull EOP & TLE data from their respected .txt file, run
% the SGP4 propagator with WGS84 ellipsoid.

%% ---------------------- REQUIRED INPUTS ---------------------------------

%% ---------------------- PERSONAL SETTINGS -------------------------------
clc; clear all; close all;
format long
% EOP PATH
addpath('C:\Users\chaitu\Documents\space and astro\4th Semester\GNC\Code\Observator_P\EOP');
% TLE PATH
addpath('C:\Users\chaitu\Documents\space and astro\4th Semester\GNC\Code\Observator_P\TLE');
% FUNCTIONS PATH
addpath('C:\Users\chaitu\Documents\space and astro\4th Semester\GNC\Code\Observator_P\Functions');
% FUNCTIONS TLE PATH
addpath('C:\Users\chaitu\Documents\space and astro\4th Semester\GNC\Code\Observator_P\Functions TLE');
% VALLADO PATH
addpath('C:\Users\chaitu\Documents\space and astro\4th Semester\GNC\Code\Observator_P\Functions Vallado');
% TLEs GROUP file txt of all the satellites to analize
file_name_in1 = 'visual.txt';

%% -------------------- OBSERVATION SETTINGS ------------------------------
R_earth = 6378.137; %km, Earth Radius
%INPUT:
OBS_DATES(1,:) = '2021 04 23 19 54 36';
OBS_DATES(2,:) = '2021 04 23 19 57 36';
SITE_COORD(1,:) = 41.765278; % LATITUDE: deg (positive towards North) [GEODETIC WSG84]
SITE_COORD(2,:) = 13.375; %LONGITUDE: deg (positive towards East) [GEODETIC WSG84]
SITE_COORD(3,:) = 555; %ALTITUDE: meters

%% ------------------ SATELLITE PERIOD FILTERING --------------------------

T = 200;        % Limit period for acceptance  
option = -1;    % Period less than limit
file_name_out1 = 'Cleaned_TLE_1.txt';
PeriodFilter(T,option,file_name_in1,file_name_out1)

%% --------------- ORBIT PROPAGATOR 1 MINUTE - ALL THE LIST ---------------

delta_t = 60; % 1 minute time step

%% ------ EOP DATA: Pull in TLE Data and Extract Epoch in UTC
% grabbing EOP data.
% EOP at the time of the observation
date_vec = str2num(OBS_DATES(1,:));
UTC_date = datetime(date_vec(1:3));

[xEOP,yEOP,dUTC_EOP, LOD_EOP,dPsi_EOP,dEpsilon_EOP,dX_EOP,dY_EOP,DAT_EOP]...
    = EOPfromcelestrak(UTC_date);

DATA_EOP(1) = table2array(xEOP);
DATA_EOP(2) = table2array(yEOP);
DATA_EOP(3) = table2array(dUTC_EOP);
DATA_EOP(4) = table2array(LOD_EOP);
DATA_EOP(5) = table2array(dPsi_EOP);
DATA_EOP(6) = table2array(dEpsilon_EOP);
DATA_EOP(7) = table2array(dX_EOP);
DATA_EOP(8) = table2array(dY_EOP);
DATA_EOP(9) = table2array(DAT_EOP);

%% --------------------- SATELLITES SETTINGS ------------------------------
% while cycle for grabbing satellites and names

file_name_in2 = file_name_out1;
[SAT_names,SAT_TLEs1,SAT_TLEs2] = ReadTLE_v3(file_name_in2);

fprintf('\n     Satellite Group - 1 minute propagation: \n')
N = size(SAT_names,1);

for i = 1:N % data from every satellite
    tic
    fprintf('Satellite %d/%d propagation in progress...', i, N);
    
    % grabbing TLE data.
    SAT_TLE = [SAT_TLEs1(i,:); SAT_TLEs2(i,:)];
    
    %% ----------------------- PROPAGATION --------------------------------
    
    [DATA_SAT_aux] = Orbit_Propagator(delta_t,DATA_EOP,OBS_DATES,SAT_TLE,SITE_COORD);
    
    %% ----------------------- DATA STRUCTURE -----------------------------
    % Collecting column vectors of [Date, Ra, Dec, Az, El, Visibility] in
    % raws for every epoch. Meanwhile those matrixes are collected in 3rd
    % dimension for every satellite
    DATA_SAT(:,:,i) = DATA_SAT_aux;
    toc
end

%% Clean DATA_SAT:
% at this point we take off the list the satellites are not visible for
% more tran 1 minute

% TLEs GROUP file txt of all the satellites to analize
file_name_out2 = 'Cleaned_TLE_2.txt';
TLE_Clean_list_v2(DATA_SAT, file_name_in2, file_name_out2);

%% ------------ ORBIT PROPAGATOR 2 - 1 SEC   ONLY SELECTED LIST -----------
% SELECT THIS PROPAGATOR IF Cleaned TLE file exists

delta_t = 1; % 1 second time step

%% --------------------- SATELLITES SETTINGS ------------------------------
% while cycle for grabbing satellites and names

file_name_in3 = file_name_out2;
[SAT_names,SAT_TLEs1,SAT_TLEs2] = ReadTLE_v3(file_name_in3);
M = size(SAT_names,1);

fprintf('\n     Satellite Cleaned Group - 1 second propagation: \n')

for j = 1:M % data from every satellite
    tic
    fprintf('Satellite %d/%d propagation in progress...',j, M);
    
    % grabbing TLE data.
    SAT_TLE = [SAT_TLEs1(j,:); SAT_TLEs2(j,:)];
    
    %% ----------------------- PROPAGATION --------------------------------
    
    [DATA_SAT_aux] = Orbit_Propagator(delta_t,DATA_EOP,OBS_DATES,SAT_TLE,SITE_COORD);
    
    %% ----------------------- DATA STRUCTURE -----------------------------
    % Collecting column vectors of [Date, Ra, Dec, Az, El, Visibility] in
    % raws for every epoch. Meanwhile those matrixes are collected in 3rd
    % dimension for every satellite
    DATA_SAT_2(:,:,j) = DATA_SAT_aux;
    toc
end

%% -------------------------- VISIBILITY PLOT -----------------------------
% Visibility lot YES / NO
plot_visibility(DATA_SAT_2, SAT_names)
% Plot of Ra, Dec, Az, El for all satellites
plot_az_el_during_visibility(DATA_SAT_2, SAT_names)

%% --------------------- FILTERING FOR ORCHESTRATE ------------------------

max_El = 70;            % [deg]
min_OBS_interval = 3;   % [min]
[DATA_SAT_OK, SAT_names_OK, DATE_OBS_start, DATE_OBS_end,...
    START_idx_OK, END_idx_OK] = Orchestrate_Filter(DATA_SAT_2, SAT_names, max_El,...
    min_OBS_interval);

N = size(SAT_names_OK,1); % number of satellites
Nt = size(DATA_SAT_OK,2); % number of time-steps

%% ----------------------SATELLITE ORDERING--------------------------------
%Order the satellites on the basis of the time start of visiblity 

J_startdate = nan(1,N); %trasform each date in julian date 
for i = 1:N
    J_startdate(i) = date2jd(DATE_OBS_start(:,i)');
end

% ++++ COMPUTE DATE TO JULIAN DATE FOR EACH TIME STEP OF OBSERVATION
% ++++ SESSION
J_times = nan(1,Nt);

for i = 1:Nt
    J_times(i) = date2jd(DATA_SAT_OK(1:6,i,1)');
end

[J_sort,Index] = sort(J_startdate,'ascend'); %order dates in ascending order

DATE_start_sorted = DATE_OBS_start(:,Index);
DATE_end_sorted = DATE_OBS_end(:,Index); 
SAT_names_sorted = SAT_names_OK(Index,:); 
DATA_SAT_sorted = DATA_SAT_OK(:,:,Index);
START_idx_s = START_idx_OK(1,Index);
END_idx_s = END_idx_OK(1,Index); 

%% --------------------------GENERATE ORCHESTRATE--------------------------

% ++++ FILE DEFINITION
orchestrate_name = 'orchestrate.txt';
satellites_name = 'orchestrate_sat.txt';
orchestrate_folder = 'Orchestrate';

if exist(fullfile(orchestrate_folder,orchestrate_name), 'file') == 2
  delete(fullfile(orchestrate_folder,orchestrate_name));
end
if exist(fullfile(orchestrate_folder,satellites_name),'file') == 2
  delete(fullfile(orchestrate_folder,satellites_name));
end

% ++++ FILE CREATION
UTC2LOC = 2; %hours
t_exposure = 2; %sec
t_between = 60; % Required time between shots [sec]
t_consider = 180; % Minimum available shooting time for passage consideration [sec]
jd_between = t_between/(24*3600);

start_exposure = nan(1,N);
end_exposure = nan(1,N);

% Generate orchestrate for 1st satellites:
start_exposure(1) = J_sort(1); %convert to julian date 
end_exposure(1) = date2jd(DATE_end_sorted(:,1)');
counter1 = START_idx_s(1);
counter2 = END_idx_s(1)-1;
delta_t = ... %Time step between instants in the Data Structure [sec]
    (end_exposure(1)-start_exposure(1))/(counter2-counter1)*24*3600;
idx_between = round(t_between/delta_t);
DATA_VIS_SAT = DATA_SAT_sorted(:,counter1:counter2,1); 
 k=0;
 i=1;
 for j = 1:idx_between:(length(DATA_VIS_SAT(1,:))-1)
     Add_Orchestrate(orchestrate_name,orchestrate_folder,DATA_VIS_SAT(:,j),t_exposure,UTC2LOC);
     Orchestrate_sat(satellites_name,orchestrate_folder,SAT_names_sorted(i,:));
     k = k+1;
     if k>=6 %check on maximum number of pictures on each satellites
         end_exposure(i) = date2jd(DATA_VIS_SAT(1:6,j))';
         break;
     end
 end

 %main loop for the other satellites:
 for i = 2:N
     
     %find the interval of visibility of each satellite
     counter1 = START_idx_s(i);
     counter2 = END_idx_s(i)-1;
     
     %convert the datetime start and end of the visibility in julian date
     start_exposure(i) =  date2jd(DATA_SAT_sorted(1:6,counter1,i))';
     end_exposure(i) = date2jd(DATA_SAT_sorted(1:6,counter2,i))';
     delta_t = ... %Time step between instants in the Data Structure [sec]
         (end_exposure(i)-start_exposure(i))/(counter2-counter1)*24*3600;
     idx_between = round(t_between/delta_t);
     
     % if the start of the exposure time is less than the end of visibility 
     %time of the previous satellite plus the required time between shots,
     %we postpone the visibility: 
     while start_exposure(i) < end_exposure(i-1) + jd_between
         start_exposure(i) = ...  %postpone start of 60 sec = 1 min 
             end_exposure(i-1) + jd_between;
         
         % LOCALIZATION OF THE INDEX FOR START OF THE OBS         
         counter1 = find(start_exposure(i) < J_times,1,'first');
     end
     
     %if there is not enough visibility time between the start of the
     %observation and its end, the observation is discarded
     if start_exposure(i) + t_consider/(24*3600) > end_exposure(i)
         end_exposure(i) = end_exposure(i-1);
     else   % else, add the passage for shooting
         DATA_VIS_SAT = DATA_SAT_sorted(:,counter1:counter2,i);
         k=0;
         for j = 1:idx_between:(length(DATA_VIS_SAT(1,:))-1)
             Add_Orchestrate(orchestrate_name,orchestrate_folder,DATA_VIS_SAT(:,j),t_exposure,UTC2LOC);
             Orchestrate_sat(satellites_name,orchestrate_folder,SAT_names_sorted(i,:));
             k=k+1;
             if k >= 6 %check on maximum number of pictures on each satellites
                 end_exposure(i) = date2jd(DATA_VIS_SAT(1:6,j))';
                 break;
             end
         end
     end
 end