function [log, x, P] = runEKF(data, cfg)
% runEKF - 16-stavový EKF pre dron (IMU + GPS)
%
% VSTUPY:
%   data - štruktúra s dátami (IMU1, GPS, ...)
%   cfg  - konfigurácia senzorov
%
% VÝSTUPY:
%   log  - štruktúra s uloženými výsledkami
%   x,P  - finálny stav a kovariancia

%% --- NASTAVENIA ---
dt = mean(diff(data.Time));
N  = length(data.Time);

start_lat = data.GPS.Pos(1,1);
start_lon = data.GPS.Pos(1,2);
start_alt = data.GPS.Pos(1,3);

%% --- INICIALIZÁCIA STAVU ---
x = zeros(16,1);
x(7) = 1; % quaternion

P = eye(16);

% Šumy
imuName = cfg.prediction.imu;
var_acc = data.(imuName).Var_Acc;
var_gyr = data.(imuName).Var_Gyr;

Q = diag([ ...
    1e-5 1e-5 1e-5, ...
    var_acc, ...
    1e-6 1e-6 1e-6 1e-6, ...
    1e-7 1e-7 1e-7, ...
    1e-7 1e-7 1e-7]);

var_gps_pos = data.GPS.Var_Pos(1);

%% --- LOGOVANIE ---
log.pos   = zeros(N,3);
log.euler = zeros(N,3);

%% =================== HLAVNÁ SLUČKA ===================
for k = 1:N

    %% --- 1. VSTUPY Z IMU (PREDIKCIA) ---
    acc = data.(imuName).Acc(k,:)';
    gyr = data.(imuName).Gyr(k,:)';

    [x_pred, F] = predictState(x, acc, gyr, dt);

    P = F*P*F' + Q;

    x = x_pred;

    % --- 2. KOREKCIA (GPS) ---
    if cfg.update.gps
        R = diag([var_gps_pos, var_gps_pos, var_gps_pos*2, 0.2, 0.2, 0.2]);
        [Z, H] = gpsMeasurement(data, k, start_lat, start_lon, start_alt);
        [x, P] = ekfUpdate(x, P, Z, H, R);
    end

    % --- 3. KOREKCIA (BAROMETR) ---
    if cfg.update.baro
        baroName = "Baro" + cfg.update.baro_id;
        
        Z = data.(baroName).Alt(k);
        
        H = zeros(1,16);
        H(3) = 1;   % p_z
        
        R = data.(baroName).Var_Alt;

        [x, P] = ekfUpdate(x, P, Z, H, R);
    end

    % --- 4. KOREKCIA (MAGNETOMETER – YAW) ---
    if cfg.update.mag
        magName = "Mag" + cfg.update.mag_id;
        
        m = data.(magName).Field(k,:)';
        m = m / norm(m);
        
        yaw_meas = atan2(m(2), m(1));   % zjednodušený heading
        yaw_est  = quat2eul(x(7:10)');
        yaw_est  = yaw_est(1);
        
        z = yaw_meas - yaw_est;
        z = wrapToPi(z);
        
        H = zeros(1,16);
        H(7:10) = [0 0 0 1];   % aproximácia vplyvu yaw
        
        R = data.(magName).Var_Field(1);
        
        [x, P] = ekfUpdate(x, P, z, H, R);
    end

    %% --- LOG ---
    log.pos(k,:)   = x(1:3)';
    log.euler(k,:) = quat2eul(x(7:10)');

end
end
