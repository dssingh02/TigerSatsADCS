clear all

addpath('IMU_calibration_package_light\')
load imu_calibration_params_04_10_2024.mat

a = arduino('COM8', 'Uno', 'Libraries', 'I2C');
imu = lsm9ds1(a);

%% measure once on demand

num_avg = 50;

data = NaN(num_avg,3);

fprintf('Press any key to measure\n');
while (true)
    input('Measure');

    for i = 1:num_avg
        data(i,:) = imu.readMagneticField();
    end
    data = calibrate_data(data, params);
    fprintf('B: (%.3f, %.3f, %.3f) %s\n', mean(data), 'microTesla');
end

%% measure continuously (option to plot)

fname = 'logs/04-18-2024_log_01.csv';
plt_flag = false;

figure(1)
clf;

fid = fopen(fname, 'a+');

while (true)
    [mag_val, t_val] = readMagneticField(imu);
    mag_val = calibrate_data(mag_val, params);

    fprintf(fid,'%s,%.8f,%.8f,%.8f\n',t_val, mag_val);
    fprintf('%s,%.8f,%.8f,%.8f\n',t_val, mag_val);
end

%%

fclose(fid);
clear;

