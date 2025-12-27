clc; clear; close all;
load('hodnoty.mat');   % MUSÍ vytvoriť "data"

cfg.prediction.imu = "IMU1";

cfg.update.gps  = true;
cfg.update.baro = false;
cfg.update.baro_id = 2;
cfg.update.mag  = true;
cfg.update.mag_id  = 2;

[log, x_final, P_final] = runEKF(data, cfg);

%% 3D POROVNANIE TRAJEKTÓRIE: EKF vs GPS

% --- GPS prepočet do lokálneho rámca (rovnaký ako v EKF) ---
start_lat = data.GPS.Pos(1,1);
start_lon = data.GPS.Pos(1,2);
start_alt = data.GPS.Pos(1,3);

gps_x = (data.GPS.Pos(:,1) - start_lat) * 111132;
gps_y = (data.GPS.Pos(:,2) - start_lon) * 111132 * cosd(start_lat);
gps_z = data.GPS.Pos(:,3) - start_alt;

% --- EKF trajektória ---
ekf_x = log.pos(:,1);
ekf_y = log.pos(:,2);
ekf_z = log.pos(:,3);

%% --- VYKRESLENIE ---
figure('Name','3D Trajektória: EKF vs GPS','Color','w');

% GPS (meranie)
plot3(gps_x, gps_y, gps_z, ...
      'r--', 'LineWidth', 1.5, 'DisplayName', 'GPS (meranie)');
hold on;

% EKF (odhad)
plot3(ekf_x, ekf_y, ekf_z, ...
      'b', 'LineWidth', 2.5, 'DisplayName', 'EKF (odhad)');

% Štart a koniec
plot3(ekf_x(1), ekf_y(1), ekf_z(1), ...
      'go', 'MarkerSize', 10, 'MarkerFaceColor','g', ...
      'DisplayName','Štart');

plot3(ekf_x(end), ekf_y(end), ekf_z(end), ...
      'rs', 'MarkerSize', 10, 'MarkerFaceColor','r', ...
      'DisplayName','Koniec');

% Kozmetika
grid on;
axis equal;
xlabel('X [m] (Sever–Juh)');
ylabel('Y [m] (Východ–Západ)');
zlabel('Z [m] (Výška)');
title('3D Rekonštrukcia trajektórie letu');
legend('Location','best');
view(45,30);

%% chyba od realnej drahy
err = [ekf_x - gps_x, ...
       ekf_y - gps_y, ...
       ekf_z - gps_z];

rmse_total = sqrt(mean(sum(err.^2, 2)));

rmse_x = sqrt(mean(err(:,1).^2));
rmse_y = sqrt(mean(err(:,2).^2));
rmse_z = sqrt(mean(err(:,3).^2));

fprintf('RMSE [X Y Z] = %.2f %.2f %.2f m\n', rmse_x, rmse_y, rmse_z);
fprintf('RMSE celkova = %.2f m\n', rmse_total);
