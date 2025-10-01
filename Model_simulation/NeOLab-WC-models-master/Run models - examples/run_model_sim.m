% Run ISN(KkSR) model with 21x21 mean inputs for 2 iterations for 2
% (thetamultiplier, sigmamultiplier) values (1/0,0) and (1,1) (Poissonlike)

basefolder = pwd();
savedir = fullfile(basefolder, 'Model_SG_1'); mkdir(savedir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niterations = 50;
iters = 1:niterations;
% thetasigma_multipliers = [[1/0,0]];
thetasigma_multipliers = [[1,1]];
uniqmultids = 1:size(thetasigma_multipliers,1);
uniqE0=8;
uniqI0=12;

%% unwrap parameter combinations
t = -1:0.5e-4:1; % seconds - for noisy inputs
% t = -0.3:0.5e-2:1; % seconds - sufficient for non-noisy input and code testrun
% t = -2.047:0.2e-5:2.048;
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
baselineinputAmp = {0*E0,0*I0};

ipdc = {0*E0, 0*I0}; % dc offset value for sinusoid

% Consider 2 intervals : [baseline, stimulus]
intervals = [-1/0, 0, 1/0];
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

[input_t, input] = generate_OUinput_rectifiedsine(t, intervals, nsims, inputDC, inputPk2Base, inputFreq, inputOnsetPhase, thetamultipliers, sigmamultipliers); 

save(inputfilename, 'simulationfilename','thetamultipliers', 'sigmamultipliers', 'intervals','inputDC', 'inputPk2Base','inputFreq', 'inputOnsetPhase', 'input_t','input','-v7.3');
%% model instantiation

JSpop = ISN_KkSR_JS_OUip(nsims, inputfilename); % or , taus);

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