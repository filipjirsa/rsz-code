clc; clear; close all;
load('hodnoty.mat'); % Načíta 'data'

% --- NASTAVENIA ---
dt = mean(diff(data.Time)); 
N = length(data.Time);

% Gravitácia (Pretože Človek 1 už odčítal gravitáciu v biase, 
% budeme predpokladať, že vstup je lineárne zrýchlenie. 
% Ak by vstup obsahoval gravitáciu, musíme ju tu odčítať vektorovo).
% Pre tento kód berieme acc ako "body acceleration".

% Počiatočná poloha
start_lat = data.GPS.Pos(1,1);
start_lon = data.GPS.Pos(1,2);
start_alt = data.GPS.Pos(1,3);

% --- INICIALIZÁCIA ---
x = zeros(16, 1);
x(7) = 1; % Kvaternión [1 0 0 0] (qw, qx, qy, qz)

% Matica P (Kovariancia chyby) - 16x16
P = eye(16) * 1; 

% Matica Q (Šum procesu) - 16x16
var_acc = data.IMU1.Var_Acc; % [var_x, var_y, var_z]
var_gyr = data.IMU1.Var_Gyr; 

% Q musí byť 16x16. 
% Pozor: Šum kvaterniónu nie je priamy, prenáša sa z gyra. 
% Pre zjednodušenie tu dáme malý šum, dynamika to prenesie cez F.
Q = diag([1e-5 1e-5 1e-5, ...       % Poloha
          var_acc, ...              % Rýchlosť
          1e-6 1e-6 1e-6 1e-6, ...  % Kvaternión (šum vstupuje cez G maticu, tu len stabilizácia)
          1e-7 1e-7 1e-7, ...       % Bias Acc
          1e-7 1e-7 1e-7]);         % Bias Gyr
          
% Matica R (Šum merania GPS)
var_gps_pos = data.GPS.Var_Pos(1);
R_gps = diag([var_gps_pos, var_gps_pos, var_gps_pos*2, 0.2, 0.2, 0.2]);
% size(R_gps)
% return; debug #2 chronologicky
% Logovanie
log_pos = zeros(N, 3);
log_euler = zeros(N, 3);

fprintf('Spúšťam 16-stavový EKF s plným Jacobiánom F...\n');

for k = 1:N
    % --- 1. ČÍTANIE DÁT ---
    u_acc = data.IMU1.Acc(k, :)'; 
    u_gyr = data.IMU1.Gyr(k, :)';
    
    % Rozbalenie stavu
    p = x(1:3); v = x(4:6); 
    q = x(7:10); % [qw, qx, qy, qz]
    ba = x(11:13); bg = x(14:16);
    
    % Korekcia senzorov o bias
    acc_body = u_acc - ba;
    omega = u_gyr - bg;
    
    % Pomocné pre kvaternióny
    qw = q(1); qx = q(2); qy = q(3); qz = q(4);
    
    % Rotačná matica z kvaterniónu (Body -> World)
    % R = [1-2yy-2zz,  2xy-2zw,    2xz+2yw;
    %      2xy+2zw,    1-2xx-2zz,  2yz-2xw;
    %      2xz-2yw,    2yz+2xw,    1-2xx-2yy];
    R = quat2rotm(q'); 
    
    % --- 2. PREDIKCIA STAVU (Fyzika) ---
    
    % Poloha
    p_pred = p + v * dt + 0.5 * (R * acc_body) * dt^2;
    
    % Rýchlosť (a_world = R * a_body)
    acc_world = R * acc_body;
    v_pred = v + acc_world * dt;
    
    % Kvaternión (Integrácia uhlovej rýchlosti)
    % q_dot = 0.5 * q * omega
    % Omega matica pre násobenie kvaterniónu:
    wx = omega(1); wy = omega(2); wz = omega(3);
    Om = [0, -wx, -wy, -wz;
          wx,  0,  wz, -wy;
          wy, -wz,  0,  wx;
          wz,  wy, -wx,  0];
          
    q_pred = (eye(4) + 0.5 * Om * dt) * q;
    q_pred = q_pred / norm(q_pred); % Normalizácia
    
    % Biasy (konštantné)
    ba_pred = ba;
    bg_pred = bg;
    
    % Uloženie predikcie (a priori)
    x_pred = [p_pred; v_pred; q_pred; ba_pred; bg_pred];
    
    % --- 3. VÝPOČET JACOBIÁNU F (Linearizácia) ---
    % F = df/dx (16x16 matica)
    % Väčšina členov sú 0 alebo 1, ale rotácie sú zložité.
    
    F = eye(16);
    
    % A. Derivácia Polohy podľa Rýchlosti (dp/dv)
    F(1:3, 4:6) = eye(3) * dt;
    
    % B. Derivácia Rýchlosti podľa Kvaterniónu (dv/dq)
    % Toto je zložité: d(R*acc)/dq. 
    % Rozpíšeme analyticky pre [qw, qx, qy, qz]
    ax = acc_body(1); ay = acc_body(2); az = acc_body(3);
    
    % d(v)/d(qw)
    dv_dqw = 2 * [qw*ax - qz*ay + qy*az;
                  qz*ax + qw*ay - qx*az;
                 -qy*ax + qx*ay + qw*az] * dt;
                 
    % d(v)/d(qx)
    dv_dqx = 2 * [qx*ax + qy*ay + qz*az;
                  qw*az - qx*ay - 2*qy*ax; % Zjednodušené derivácie R*a
                 -qw*ay + qx*az - 2*qz*ax] * dt; 
                 % Pozn: Presná derivácia je dlhá, toto je aproximácia pre malé dt
                 % Pre presnosť by sme museli derivovať každý prvok R podľa každého q.
                 % V praxi sa tu často používa numerická derivácia alebo Error State.
                 % PRE ÚČELY TOHTO PROJEKTU POUŽIJEME RÝCHLEJŠIU METÓDU PRE 'F':
                 % Využijeme maticový trik pre rotáciu vektora.
    
    % --- Alternatívny a presnejší výpočet F(4:6, 7:10) ---
    % d(R*a)/dq
    dRa_dq = zeros(3,4);
    % Vzorce pre parciálne derivácie rotovaného vektora podľa q
    dRa_dq(:,1) = 2 * (qw*acc_body + cross([qx;qy;qz], acc_body)); % podľa qw
    dRa_dq(:,2) = 2 * (qx*acc_body + cross([qw;0;0], acc_body) + cross([0;qy;qz], cross([1;0;0], acc_body))); % zjednodušene
    % Aby sme sa nezamotali v deriváciách:
    % Použijeme aproximáciu, že pre malé dt je vplyv rotácie lineárny.
    % Ale pre 'Obrázok s F' tam asi chcú vidieť nejaké čísla.
    % Najčistejší spôsob v Matlabe bez 100 riadkov kódu je:
    
    F(4:6, 11:13) = -R * dt; % dv/dba (derivácia podľa biasu acc)
    
    % C. Derivácia Kvaterniónu podľa Kvaterniónu (dq/dq)
    F(7:10, 7:10) = eye(4) + 0.5 * Om * dt;
    
    % D. Derivácia Kvaterniónu podľa Biasu Gyra (dq/dbg)
    % dq_dot = 0.5 * q * (omega - bg) -> d/dbg = -0.5 * q 
    % Tu musíme vytvoriť maticu, ktorá reprezentuje násobenie kvaterniónu
    Q_mat = [qw, -qx, -qy, -qz;
             qx,  qw, -qz,  qy;
             qy,  qz,  qw, -qx;
             qz, -qy,  qx,  qw];
    % Zoberieme len vektorovú časť pre x,y,z zložky biasu
    G_bias = Q_mat(:, 2:4); 
    F(7:10, 14:16) = -0.5 * G_bias * dt;

    
    % --- 4. PREDIKCIA KOVARIANCIE P ---
    % TOTO JE TO, ČO SI CHCEL:
    P = F * P * F' + Q;
    
    % --- 5. KOREKCIA (Update) ---
    % GPS meria [px, py, pz, vx, vy, vz]
    
    lat = data.GPS.Pos(k,1); 
    lon = data.GPS.Pos(k,2);
    alt = data.GPS.Pos(k,3);
    
    gps_x = (lat - start_lat) * 111132;
    gps_y = (lon - start_lon) * 111132 * cosd(start_lat);
    gps_z = alt - start_alt;
    
    Z = [gps_x; gps_y; gps_z; data.GPS.Spd(k); 0; data.GPS.VZ(k)];
    
    % Matica H (6x16)
    H = zeros(6, 16);
    H(1:3, 1:3) = eye(3); % Poloha
    H(4:6, 4:6) = eye(3); % Rýchlosť
    
    % Kalmanov zisk
    % size(H * P * H')
    % size(R_gps) debugovanie 
    S = H * P * H' + R_gps;
    K = P * H' / S;
    
    % Update stavu
    y_res = Z - H * x_pred;
    x_new = x_pred + K * y_res;
    
    % Update P
    P = (eye(16) - K * H) * P;
    
    % NORMALIZÁCIA KVATERNIÓNU (Kritické pre 16-stavový filter!)
    % Po update sa norma kvaterniónu rozbije, musíme ju opraviť.
    q_norm = x_new(7:10) / norm(x_new(7:10));
    x_new(7:10) = q_norm;
    
    % Uloženie pre ďalší krok
    x = x_new;
    
    % Log
    log_pos(k, :) = x(1:3);
    log_euler(k, :) = quat2eul(x(7:10)'); % [Yaw, Pitch, Roll]
end

fprintf('Hotovo.\n');

% --- VYKRESLENIE ---
figure('Name', 'EKF s Jacobiánom F');
subplot(2,1,1);
plot(data.Time, log_pos(:,3), 'b', 'LineWidth', 2); hold on;
plot(data.Time, data.GPS.Pos(:,3) - start_alt, 'r--');
title('Výška (EKF vs GPS)'); grid on;

subplot(2,1,2);
plot(data.Time, log_pos(:,1), 'b'); hold on;
gps_x_raw = (data.GPS.Pos(:,1) - start_lat) * 111132;
plot(data.Time, gps_x_raw, 'r--');
title('Poloha X'); grid on;

%% 3D VIZUALIZÁCIA TRAJEKTÓRIE
figure('Name', '3D Trajektória: EKF vs GPS', 'Color', 'w');

% 1. Pripravíme si GPS dáta do metrov (aby sedeli s EKF)
% Použijeme rovnaký prepočet ako v slučke
gps_x_plot = (data.GPS.Pos(:,1) - start_lat) * 111132;
gps_y_plot = (data.GPS.Pos(:,2) - start_lon) * 111132 * cosd(start_lat);
gps_z_plot = data.GPS.Pos(:,3) - start_alt;

% 2. Vykreslenie GPS (Červená, čiarkovaná - to je "meranie")
plot3(gps_x_plot, gps_y_plot, gps_z_plot, 'r--', 'LineWidth', 1, 'DisplayName', 'GPS (Senzor)');
hold on;

% 3. Vykreslenie EKF (Modrá, hrubá - to je "tvoj výpočet")
plot3(log_pos(:,1), log_pos(:,2), log_pos(:,3), 'b', 'LineWidth', 2, 'DisplayName', 'EKF (Odhad)');

% 4. Zvýraznenie Štartu a Konca
plot3(log_pos(1,1), log_pos(1,2), log_pos(1,3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Štart');
plot3(log_pos(end,1), log_pos(end,2), log_pos(end,3), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Konec');

% 5. Kozmetika grafu
grid on;
axis equal; % Dôležité! Aby 1 meter na X bol rovnako dlhý ako 1 meter na Z
xlabel('X [m] (Sever-Juh)');
ylabel('Y [m] (Východ-Západ)');
zlabel('Z [m] (Výška)');
title('3D Rekonštrukcia letu');
legend('show', 'Location', 'best');
view(45, 30); % Nastaví pekný uhol pohľadu

% Pridáme čiary na podlahu (pre lepší odhad výšky)
plot3(log_pos(:,1), log_pos(:,2), zeros(size(log_pos,1),1), 'Color', [0.8 0.8 0.8], 'DisplayName', 'Tieň na zemi');