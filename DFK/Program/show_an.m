cfg;

D               = importdata('../../INS/Output_data/output_INS_data.txt');  
time            = D.data(:, 1);
clear D

t_val           = time(1:TimeEnd);
t_min           = min(t_val);
t_max           = max(t_val);
TimeWave = t_val - t_min;

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
    GVD(:, i) = Psi3d * X_j(13:end, end);
end

figure(1)
hold
plot(X_j(33:4:52, :)')

% figure(1)
% hold
% % plot(time(1:time_end*10), 10^5*(GVD(1, :)-true_anomaly_all(1:10:100*time_end,2)'),  LineWidth=1.5)
% plot(time(1:time_end*10), 10^5*GVD(1, :),  LineWidth=1.5)
% plot(time(1:time_end*10), 10^5*true_anomaly_all(1:10:100*time_end,2)',  LineWidth=1.5)
% % legend('est', 'true')
% 
% figure(2)
% hold
% % plot(time(1:time_end*10), 10^5*(GVD(2, :)-true_anomaly_all(1:10:100*time_end,3)'),  LineWidth=1.5)
% plot(time(1:time_end*10), 10^5*GVD(2, :),  LineWidth=1.5)
% plot(time(1:time_end*10), 10^5*true_anomaly_all(1:10:100*time_end,3)',  LineWidth=1.5)
% figure(3)
% hold
% % plot(time(1:time_end*10), 10^5*(GVD(3, :)-true_anomaly_all(1:10:100*time_end,4)'),  LineWidth=1.5)
% plot(time(1:time_end*10), 10^5*GVD(3, :),  LineWidth=1.5)
% plot(time(1:time_end*10), 10^5*true_anomaly_all(1:10:100*time_end,4)',  LineWidth=1.5)
% % plot(time(1:time_end*10), 10^5*(GVD(3, :)true_anomaly_all(1:10:100*time_end,4)), LineWidth=1.5)
% % plot(true_anomaly_all(1:10:time_end*100, 1),10^5*true_anomaly_all(1:10:100*time_end,4), LineWidth=1.5);
% legend('est', 'true')