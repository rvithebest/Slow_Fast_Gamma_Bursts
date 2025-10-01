function generate_WC_OU_model_gamma(indicator,folder_name,theta_E,theta_I)
% folder_name='Model_SG_OU_noise_SOM_2';
basefolder = pwd();
savedir = fullfile(basefolder, folder_name); mkdir(savedir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niterations = 50;
iters = 1:niterations;
thetasigma_multipliers = [[1,1]];
uniqmultids = 1:size(thetasigma_multipliers,1);
uniqE0=8;
uniqI0=8;
%% unwrap parameter combinations
t = -2.047:2e-5:2.048;
% t = -0.8475:2e-5:1.1965;
paramnames = {'E0', 'I0', 'multiplierIDs', 'iterIDs'};
params = { uniqE0, uniqI0, uniqmultids, iters};
[varVals, varSels, nsims] = unwrapParameters(params);
E0 = varVals(:,1);
I0 = varVals(:,2);
iterIDs = varVals(:,4);
multIDs = varVals(:,3);
thetamultipliers = {thetasigma_multipliers(multIDs,1),thetasigma_multipliers(multIDs,1)};
sigmamultipliers = {thetasigma_multipliers(multIDs,2),thetasigma_multipliers(multIDs,2)};
%% save directories
inputfilename = fullfile(savedir,'inputs.mat');
simulationfilename = fullfile(savedir,'simulationResults.mat');

%%%%%%%%%%%%%%
inputVal = {E0, I0};
baselineinputAmp = {0*E0, 0*I0};

ipdc = {0*E0, 0*I0}; % dc offset value for sinusoid

% Consider 2 intervals : [baseline, stimulus]
% intervals = [-1/0, 0, 1/0];
intervals= [-1/0,0,1.5];
% intervals= [-1/0,0,0.8];
% DC component
inputDC = {[],[]}; % {1}-E, {2}-I
% Cosine component - to be set to 0 amplitude
inputPk2Base = {[],[]}; 
inputFreq = {[],[]}; 
inputOnsetPhase = {[],[]}; 
for i = 1:numel(inputPk2Base) % {1}-E, {2}-I
    inputDC{i} = baselineinputAmp{i}*[1 0] + inputVal{i}*[0 1]; % specifying DC component of input for each input combination in each interval ([1,0]-baseline; [0,1]-stimulus)
    inputPk2Base{i} = [0,0];
    inputFreq{i} = [0,0]; % constant contrast stimuli; frequency of cosine signal input
    inputOnsetPhase{i} = [0,0]; % initial phase of cosine contrast mask at onset instant (start of interval)
end
%%%%% OU noise %%%%%%%%%%%%%%
outhetas = {2*pi*theta_E, 2*pi*theta_I};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[input_t, input] = generate_OUinput_rectifiedsine(t, intervals, nsims, inputDC, inputPk2Base, inputFreq, inputOnsetPhase, thetamultipliers, sigmamultipliers, outhetas); 

save(inputfilename, 'simulationfilename','thetamultipliers', 'sigmamultipliers', 'intervals','inputDC', 'inputPk2Base','inputFreq', 'inputOnsetPhase', 'input_t','input','-v7.3');
%% model instantiation
if indicator
   taus = (10/13)*([20,10]*1e-3); % seconds- fast gamma
else 
   taus = (10/13)*([20,20]*1e-3); % seconds- slow gamma
end
if indicator
   Weights = [16, -26 ; 20, -1]; % PV
else 
   Weights = [16, -26 ; 20, 0]; % SOM
end
theta = [7.5 18];
JSpop = ISN_KkSR_JS_OUip(nsims, inputfilename, taus,theta, Weights); % or , taus);

JSpop.input([0,0], input_t); % set timestamps of simulation

%% model simulation and lfp proxy
solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);
JSpop.ode(solver);

t = JSpop.EIpairs.t;
rE = JSpop.EIpairs.R(1:end/2, :);
rI = JSpop.EIpairs.R(end/2+1:end, :);
lfp = - rE - rI;
save(simulationfilename, 'JSpop','lfp','inputfilename','niterations','uniqE0','uniqI0',"thetasigma_multipliers",'-v7.3');
clear JSpop