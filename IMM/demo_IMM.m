clear
clc
close all

%% target
% parameter initialization
T = 0.1;
nSteps = 200;
nModels = 2;
% `ID_model` update
ID_model = zeros(nSteps, 1);
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
% `nStates` update
nStates = [2; 3];
% `sigma_v`, `F`, `Gamma` update
sigma_v = zeros(nModels, 1);
F = cell(nModels, 1);
Gamma = cell(nModels, 1);
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
% target update 
X = cell(nSteps, 1);
x_0 = [0; 1; 0];
for i = 1 : nSteps
    if i == 1
        X{i} = F{ID_model(i)} * x_0(1:nStates(ID_model(i))) + Gamma{ID_model(i)} * normrnd(0, sigma_v(ID_model(i)));
    else
        if size(F{ID_model(i)}, 1) <= length(X{i - 1})
            X{i} = F{ID_model(i)} * X{i - 1}(1:nStates(ID_model(i))) + Gamma{ID_model(i)} * normrnd(0, sigma_v(ID_model(i)));
        else
            x_pseudo = [X{i - 1}; 0];
            X{i} = F{ID_model(i)} * x_pseudo + Gamma{ID_model(i)} * normrnd(0, sigma_v(ID_model(i)));
        end
    end
end
% log
X_log = zeros(nSteps, length(x_0));
for i = 1 : nSteps
    if length(X{i}) == 3
        X_log(i, :) = X{i}';
    else
        X_log(i, :) = [X{i}; 0]';
    end
end
%% measurement
% parameter initialization
Z = zeros(nSteps, 1);
H = cell(nModels, 1);
for i = 1 : nModels
    switch i
        case 1
            H{i} = [1, 0];
        case 2
            H{i} = [1, 0, 0];
    end
end
% observation update
for i = 1 : nSteps
    Z(i) = H{ID_model(i)} * X{i};
end