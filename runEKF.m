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

Rmultiplier = 250;
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

var_lat = data.GPS.Var_Pos(1);
var_lon = data.GPS.Var_Pos(2);
var_alt = data.GPS.Var_Pos(3);
% Prevod lat/lon do metru
scale_x = 111132;                 % m / deg (lat)
scale_y = 111132 * cosd(start_lat);% 

%% --- LOGOVANIE ---
log.pos   = zeros(N,3);
log.euler = zeros(N,3);
log.nis.gps  = nan(N,1);
log.nis.baro = nan(N,1);
log.nis.mag  = nan(N,1);

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
        R = diag([var_lat*scale_x^2, var_lon*scale_y^2, var_alt, 0.2, 0.2, 0.2]*Rmultiplier);
        [Z, H] = gpsMeasurement(data, k, start_lat, start_lon, start_alt);
        [x, P, nis_gps] = ekfUpdate(x, P, Z, H, R);
        log.nis.gps(k) = nis_gps;
    end

    % --- 3. KOREKCIA (BAROMETR) ---
    if cfg.update.baro
        baroName = "Baro" + cfg.update.baro_id;
        
        Z = data.(baroName).Alt(k);
        
        H = zeros(1,16);
        H(3) = 1;   % p_z
        
        R = data.(baroName).Var_Alt*Rmultiplier;

        [x, P, nis_baro] = ekfUpdate(x, P, Z, H, R);
        log.nis.baro(k) = nis_baro;

    end

    % --- 4. KOREKCIA (MAGNETOMETER – YAW) ---
    if cfg.update.mag
        magName = "Mag" + cfg.update.mag_id;
        
        m = data.(magName).Field(k,:)';
        m = m / norm(m);
        
        yaw_meas = atan2(m(2), m(1));   % zjednodušený heading
        yaw_est  = yawFromQuat(x(7:10));
        nu = wrapToPi(yaw_meas - yaw_est);

        % Numerická derivace yaw wrt quaternion
        eps = 1e-6;
        H = zeros(1,16);
        
        q0 = x(7:10);
        for i = 1:4
            dq = zeros(4,1); dq(i) = eps;
            qp = q0 + dq; qp = qp / norm(qp);
            qm = q0 - dq; qm = qm / norm(qm);
        
            hp = yawFromQuat(qp);
            hm = yawFromQuat(qm);
        
            d = wrapToPi(hp - hm) / (2*eps);
            H(6+i) = d;  % indexy 7:10
        end
        
        R = data.(magName).Var_Field(1)*Rmultiplier;
        
        [x, P, nis_mag] = ekfUpdate(x, P, nu, H, R);
        log.nis.mag(k) = nis_mag;

    end

    %% --- LOG ---
    log.pos(k,:)   = x(1:3)';
    log.euler(k,:) = quat2eul(x(7:10)');
end
end
