%% Convert the Matlab ADV data to universal txt
%
% Note:
% - shift the time by 189 s; 
% - project to East/North reference: vector rotation of -25º (i.e. 25º CW)
%   corresponds to a reference rotation of +25º (25º CCW).
%
% Written by P. DERIAN 2017-01-16

%% which day
day = 13;    %change the day here
year = 2014; %do not edit
month = 3;   %do not edit
% Note: data is at 8 Hz
switch day
    case 13
        % 2014-03-13 data span indices 1:460800
        % over 08:00 (?) to 17:00
        imin = 1;
        imax = 460800;    
    otherwise
        error('This day (%04d-%03d-%02d) is not yet configured', year, month, day);
end


%% parameters
input = 'vel_13-16mars2014_XYZ'; 
output = sprintf('ADVdata_%04d%02d%02d_refEastNorth.txt', year, month, day);
t_format = 'yyyy-mm-dd HH:MM:SS.FFF000'; % the format for time display
t_shift = 189.; % [s] the shift to apply: t' = t + t_shift
a_rotation = -25.; % [deg] the angle of the reference to project onto.

%% load
load(input);

%% write output
cos_a = cosd(a_rotation);
sin_a = sind(a_rotation);
dt = t_shift/(3600.*24.); % dates are expressed in days
% header
f = fopen(output, 'w');
fprintf(f, '#ADV data for Grand Popo, %04d-%02d-%02d\n', year, month, day);
fprintf(f, '#generated from "%s"\n', input);
fprintf(f, '#applying (1) a time shift of %.2f s\n', t_shift);
fprintf(f, '#(2) changing the sign of second column (UADV)\n');
fprintf(f, '#and (3) a rotation of %.2f deg to compensate compas misalignment.\n', a_rotation);
fprintf(f, '#Reference frame: Easting-Northing\n');
fprintf(f, '#date, time, easting, northing, W\n');
% data
for i=imin:imax
    tmp_t = datestr(tADV(i)+dt, t_format);
    % velocities were recorded in [cm/s]
    tmp_u = A(i,2)/100.; % theoretical South-North axis, positive towards S
    tmp_v = A(i,3)/100.; % theoretical West-East axis, positive towards E
    tmp_w = A(i,4)/100.;
    % Here we do
    % (u, v) => (e, n)=(v, -u) => rotation(a)
    tmp_e = tmp_v; %easting
    tmp_n = -tmp_u; %northing 
    tmp_l = cos_a*tmp_e - sin_a*tmp_n; % corrected "easting"
    tmp_c = sin_a*tmp_e + cos_a*tmp_n; % corrected "northing"
    fprintf(f, '%s %.3f %.3f %.3f\n', tmp_t, tmp_l, tmp_c, tmp_w);
end
fclose(f);
