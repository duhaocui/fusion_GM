clear
clc
close all

%% target
% parameter initialization
T = 0.1;
nSteps = 200;
nModels = 2;
nStates = [2; 3];
x_tgt_init = [3; 2; 1];
R = 0.1;
% variable initialization
ID_model = zeros(nSteps, 1);
sigma_v = zeros(nModels, 1);
F = cell(nModels, 1);
Gamma = cell(nModels, 1);
x_tgt = cell(nSteps, 1);
x_tgt_log = zeros(nSteps, length(x_tgt_init));
% `ID_model` update
for i = 1 : nSteps
    if i <= 50
        ID_model(i) = 1;
    elseif i <= 70
        ID_model(i) = 2;
    elseif i <=120
        ID_model(i) = 1;
    elseif i <= 150
        ID_model(i) = 2;
    elseif i <= 200
        ID_model(i) = 1;
    end
end
% update of `sigma_v`, `F`, `Gamma`
for i = 1 : nModels
    switch i
        case 1
            % discrete white noise acceleration model
            sigma_v(i) = 0.1;
            F{i} = [1, T; 0, 1];
            Gamma{i} = [T^2 / 2; T];
        case 2
            % discrete Wiener process acceleration model
            sigma_v(i) = 1;
            F{i} = [1, T, T^2 / 2; 0, 1, T; 0, 0, 1];
            Gamma{i} = [T^2 / 2; T; 1];
    end
end
% target `X` update
for i = 1 : nSteps
    if i == 1
        x_tgt{i} = F{ID_model(i)} * x_tgt_init(1:nStates(ID_model(i))) + Gamma{ID_model(i)} * normrnd(0, sigma_v(ID_model(i)));
    else
        if size(F{ID_model(i)}, 1) <= length(x_tgt{i - 1})
            x_tgt{i} = F{ID_model(i)} * x_tgt{i - 1}(1:nStates(ID_model(i))) + Gamma{ID_model(i)} * normrnd(0, sigma_v(ID_model(i)));
        else
            x_pseudo = [x_tgt{i - 1}; 0];
            x_tgt{i} = F{ID_model(i)} * x_pseudo + Gamma{ID_model(i)} * normrnd(0, sigma_v(ID_model(i)));
        end
    end
end
% archive of `X_log` 
for i = 1 : nSteps
    if length(x_tgt{i}) == 3
        x_tgt_log(i, :) = x_tgt{i}';
    else
        x_tgt_log(i, :) = [x_tgt{i}; 0]';
    end
end
%% measurement
% parameter initialization
H = cell(nModels, 1);
for i = 1 : nModels
    switch i
        case 1
            H{i} = [1, 0];
        case 2
            H{i} = [1, 0, 0];
    end
end
% variable initialization
Z = zeros(nSteps, 1);
% observation update
for i = 1 : nSteps
    Z(i) = H{ID_model(i)} * x_tgt{i} + normrnd(0, sqrt(R));
end
%% IMM
% parameter initialization
p = [0.98, 0.02; 0.02, 0.98]; % mode transition probability matrix
mu_init = [0.9; 0.1];
x_IMM_init = [3; 2; 1];
P_IMM_init = diag([0.1, 0.1, 0.5]);
% variable initialization
c_bar = zeros(nModels, 1);
mu = zeros(nSteps, nModels); % mode probability
mu_cond = zeros(nModels, nModels); % conditional mode probability
x_prev = cell(nModels, 1);
P_prev = cell(nModels, 1);
x_local_pred = cell(nModels, 1);
P_local_pred = cell(nModels, 1);
x_local_upd = cell(nModels, 1);
P_local_upd = cell(nModels, 1);
Lambda = zeros(nModels, 1); % likelihood of the measurement for each mode
for i = 1 : nModels
    x_local_pred{i} = zeros(nSteps, nStates(i));
    P_local_pred{i} = zeros(nSteps, nStates(i), nStates(i));
    x_local_upd{i} = zeros(nSteps, nStates(i));
    P_local_upd{i} = zeros(nSteps, nStates(i), nStates(i));
end
x_mixed_init = cell(nModels, 1);
P_mixed_init = cell(nModels, 1);
% main loop
for i = 1 : nSteps
    % step 1: calculation of the mixing probabilities
    % `mu_prev` update
    if i == 1
        mu_prev = mu_init;
    else
        mu_prev = (mu(i - 1, :))';
    end
    % calculation of normalizing constant `c_bar(j)` for j = 1, ... , nModels 
    for j = 1 : nModels
        for k = 1 : nModels
            if k == 1
                c_bar(j) = p(k, j) * mu_prev(k);
            else
                c_bar(j) = c_bar(j) + p(k, j) * mu_prev(k);
            end
        end
    end
    % calculation of conditioanl mode probability `mu_cond(j, k)` for j, k = 1, ..., nModels
    for j = 1 : nModels
        for k = 1 : nModels
            mu_cond(j, k) = p(j, k) * mu_prev(j) / c_bar(k);
        end
    end
    % step 2: mixing
    % update of `x_prev{j}` and `P_prev{j}` for j = 1, ... , nModels 
    for j = 1 : nModels
        if i == 1
            x_prev{j} = x_IMM_init(1 : nStates(j));
            P_prev{j} = P_IMM_init(1 : nStates(j), 1 : nStates(j));
        else
            x_prev{j} = (x_local_upd{j}(i - 1, :))';
            P_prev{j} = shiftdim(P_local_upd{j}(i - 1, :, :));
        end
    end
    % calculation of `x_mixed_init{j}` for j = 1, ... , nModels  
    for j = 1 : nModels
        for k = 1 : nModels
            if k == 1
                x_mixed_init{j} = x_prev{k} * mu_cond(k, j);                
            else
                x_mixed_init{j} = vector_sum(x_mixed_init{j}, x_prev{k} * mu_cond(k, j), nStates(j));
            end
        end
    end
    % calculation of `P_mixed_init{j}` for j = 1, ... , nModels  
    for j = 1 : nModels
        for k = 1 : nModels
            if k == 1
                P_mixed_init{j} = mu_cond(k, j) * square_sum(P_prev{k}, (vector_sum(x_prev{k}, -x_mixed_init{j}, nStates(j))) * (vector_sum(x_prev{k}, -x_mixed_init{j}, nStates(j)))', nStates(j));
            else
                P_mixed_init{j} = square_sum(P_mixed_init{j}, mu_cond(k, j) * square_sum(P_prev{k}, (vector_sum(x_prev{k}, -x_mixed_init{j}, nStates(j))) * (vector_sum(x_prev{k}, -x_mixed_init{j}, nStates(j)))', nStates(j)), nStates(j));
            end
        end
    end
    % step 3: Mode-matched filtering
    z = Z(i);
    for j = 1 : nModels
        Q = Gamma{j} * (Gamma{j})' * (sigma_v(j))^2;
        % Kalman filter, prediction
        [x_pred_KF, P_pred_KF] = KF_pred(x_mixed_init{j}, P_mixed_init{j}, F{j}, Q);
        % Kalman filter, update
        [x_upd_KF, P_upd_KF] = KF_upd(z, x_pred_KF, P_pred_KF, H{j}, R);
        % archive of `x_local{j}` and `P_local{j}` for j = 1, ..., nModels
        x_local_pred{j}(i, :) = x_pred_KF';
        P_local_pred{j}(i, :, :) = P_pred_KF;
        x_local_upd{j}(i, :) = x_upd_KF';
        P_local_upd{j}(i, :, :) = P_upd_KF;
    end
    % calculation of likelihood function `Lambda(j)` for j = 1, ..., nModels
    for j = 1 : nModels
        Lambda(j) = normpdf(z, H{j} * (x_local_pred{j}(i, :))', sqrt(H{j} * shiftdim(P_local_pred{j}(i, :, :)) * (H{j})' + R));
    end
    % step 4: mode probability update
    for j = 1 : nModels
        if j == 1
            c = Lambda(j) * c_bar(j);
        else
            c = c + Lambda(j) * c_bar(j);
        end
    end
    % update of mode probability `mu`
    for j = 1 : nModels
        mu(i, j) = Lambda(j) * c_bar(j) / c;
    end
end

%% test
figure
plot(1:nSteps, mu(:, 1), 1:nSteps, mu(:, 2))
grid on
legend('mode 1', 'mode 2')
