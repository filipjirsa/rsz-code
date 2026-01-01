function out = detectFaultySensor(results, cfgs)
    % Nastavení přesnost
    alpha = 0.999;   
    W = 50;        
    threshold_ratio = 0.5; 
    
    dof_gps  = 6;   
    dof_baro = 1;   
    dof_mag  = 1;   
    nFilt = numel(results);
    
    sensorList = ["IMU1","IMU2","IMU3","Baro1","Baro2","Mag1","Mag2","GPS"];
    nSens = numel(sensorList);
    M = zeros(nFilt, nSens);
    
    % Sestavení matice M 
    for i = 1:nFilt
        imuName = string(cfgs(i).prediction.imu);
        M(i, sensorList == imuName) = 1;
        baroName = "Baro" + cfgs(i).update.baro_id;
        M(i, sensorList == baroName) = 1;
        magName = "Mag" + cfgs(i).update.mag_id;
        M(i, sensorList == magName) = 1;
        M(i, sensorList == "GPS") = 1; 
    end

    %Upravené rozhodování o vadě (vektor F) 
    F = false(nFilt,1);
    for i = 1:nFilt
        % Získání residuí (NIS)
        nis_gps  = results(i).log.nis.gps;
        nis_baro = results(i).log.nis.baro;
        nis_mag  = results(i).log.nis.mag;
        
        % Teoretické prahy pro jeden vzorek 
        Tg = chi2inv(alpha, dof_gps);
        Tb = chi2inv(alpha, dof_baro);
        Tm = chi2inv(alpha, dof_mag);
        
        % Kontrola, zda residua překračují práh soustavně 
        is_bad_gps  = movmean(nis_gps  > Tg, W, 'omitnan') > threshold_ratio;
        is_bad_baro = movmean(nis_baro > Tb, W, 'omitnan') > threshold_ratio;
        is_bad_mag  = movmean(nis_mag  > Tm, W, 'omitnan') > threshold_ratio;
        
        % Pokud alespoň v jednom čase došlo k trvalému překročení prahu
        F(i) = any(is_bad_gps) || any(is_bad_baro) || any(is_bad_mag);
    end

    % Výpočet pravděpodobnosti P
    Frow = double(F(:))';
    denom = ones(1, nFilt) * M;
    numer = Frow * M;
    P = numer ./ max(denom, 1e-12);

    [maxP, idx] = max(P);
    % Pokud je maxP rovno 0, žádný senzor není vadný
    if maxP > 0
        winner = sensorList(idx);
    else
        winner = "None";
    end

    out.F = F;
    out.P = P;
    out.sensorList = sensorList;
    out.winner = winner;
end