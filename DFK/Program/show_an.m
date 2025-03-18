GVD = zeros(3, TimeEnd);
Psi = zeros(1, sum(k_arr));
for i = 1:TimeEnd
    for j = j_min : j_max
        for k = 1 : k_arr(j - j_min + 1)
           dT  = (2^j) * t_max / k_arr(j - j_min + 1);  % step of wavelet grid
           t_k = k * dT;                % knot in grid
           Psi(:, sum(k_arr(1: j - j_min)) + k) = wavel_trf(j,t_k,TimeWave(i));   
        end
    end
    Psi3d = [Psi, zeros(size(Psi)), zeros(size(Psi));
             zeros(size(Psi)),  Psi, zeros(size(Psi));
             zeros(size(Psi)), zeros(size(Psi)), Psi;];
    GVD(:, i) = Psi3d * X_j(13:end, i);
end

% figure(1)
% hold
% plot(time(1:time_end*10), GVD(1, :))
% plot(true_anomaly_res(:, 1),true_anomaly_res(:,2));

% figure(2)
% plot(time(1:time_end*10), GVD(2, :))
% plot(true_anomaly_res(:, 1),true_anomaly_res(:,3), LineWidth=1.5);
figure(3)
hold
plot(time(1:time_end*10), GVD(3, :))
plot(true_anomaly_all(:, 1),true_anomaly_all(:,4), LineWidth=1.5);