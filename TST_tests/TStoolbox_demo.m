%% let's make some continuous "data", simulating the firing rate of 
duration = 5;

Fs = 30000; % 30 kHz (sampling rate for open-ephys) 
EEG0_times = (1:(Fs*duration))'; % five seconds of simulated EEG, in units of the sampling rate. Must be a column vetor
EEG0_data = rand(size(EEG0_times)); % EEG is white noise 


EEG1_times = EEG0_times;
EEG1_data = sin(2 * pi * 8 * EEG1_times / Fs);


%% package the data in a tsd object
EEG0 = tsd(EEG0_times * (1./Fs), EEG0_data); % a tsd object contains the times and the timestamps. 
                                              % the Tsd object expect times
                                              % in units of seconds!
EEG1 = tsd(EEG1_times * (1./Fs), EEG1_data ); % again, note unit conversion for the times (to seconds)

%% plot some EEG 

figure(1)
clf
plot(Range(EEG1), Data(EEG1))
hold on 
% plot(Range(EEG0), Data(EEG0))

%% let's make some "position data"...

% let's have a mouse running on a circular track, very regularly with a
% period of 20 s 

Pos_times = (1:(Fs*duration))'; % five minutes of simulated position data, in units of the sampling rate of the acquisition system
                          % that assumes synchronization and rescaling of
                          % video timestamps! (must be column vector)
                          
X_data = 300 * (0.5 + 0.5 * cos(2 * pi * (1/2) * Pos_times / Fs)) + 5 * rand(size(Pos_times));                      
Y_data = 300 * (0.5 + 0.5 * sin(2 * pi * (1/2) * Pos_times / Fs)) + 5 * rand(size(Pos_times));                      


%% ...and let's package them in tsd objects 
X = tsd(Pos_times * (1. / Fs), X_data);
Y = tsd(Pos_times * (1. / Fs), Y_data);


%% plot animal position 

figure(2)
clf
plot(Data(X), Data(Y))
axis equal tight

%% let's make random "spike trains"

Spikes = cell(1,10);
for i = 1:10 % 10 cells 
    n_spikes = exprnd(5); % mean firing rate will be exponentially distributed with a population mean of 5 Hz
    sp = sort(rand(floor(n_spikes * duration), 1) * duration * Fs);
    Spikes{i} = ts(sp * (1. / Fs)); % spike times in units of seconds, covering a 5 min period
end

Spikes = tsdArray(Spikes);
