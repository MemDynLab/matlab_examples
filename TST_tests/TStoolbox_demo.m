%% let's make some continuous "data", simulating the firing rate of 
duration = 5;

Fs = 30000; % 30 kHz (sampling rate for open-ephys) 
EEG0_times = (1:(Fs*duration))'; % five seconds of simulated EEG, in units of the sampling rate. Must be a column vector
EEG0_data = rand(size(EEG0_times)); % EEG is white noise 


EEG1_times = EEG0_times;
EEG1_data = sin(2 * pi * 8 * EEG1_times / Fs);


%% package the data in a tsd object
EEG0 = tsd(EEG0_times * (1./Fs), EEG0_data); % a tsd object contains the times and the timestamps. 
                                              % the Tsd object expect times
                                              % in units of seconds!
                                              % ts and tsd expect COLUMN
                                              % VECTORS as inputs!
EEG1 = tsd(EEG1_times * (1./Fs), EEG1_data ); % again, note unit conversion for the times (to seconds)

%% plot some EEG 

figure(1)
clf
plot(Range(EEG1), Data(EEG1))
hold on 
% plot(Range(EEG0), Data(EEG0))

%% let's make some "position data"...

% let's have a mouse running on a circular track, very regularly with a
% period of 2 s 

Fs_position = 30;
Pos_times = (1:(Fs_position*duration))'; % five minutes of simulated position data, in units of the sampling rate of the acquisition system
                          % that assumes synchronization and rescaling of
                          % video timestamps! (must be column vector)
                          
X_data = 300 * (0.5 + 0.5 * cos(2 * pi * (1/2) * Pos_times / Fs_position)) + 5 * rand(size(Pos_times));                      
Y_data = 300 * (0.5 + 0.5 * sin(2 * pi * (1/2) * Pos_times / Fs_position)) + 5 * rand(size(Pos_times));                      
 
% here position is given in units of pixels, you may want to convert to
% centimeters at this stage 

%% ...and let's package them in tsd objects 
X = tsd(Pos_times * (1. / Fs_position), X_data);
Y = tsd(Pos_times * (1. / Fs_position), Y_data);


%% plot animal position 

figure(2)
clf
plot(Data(X), Data(Y))
axis equal tight

%% let's make random "spike trains"

Fs_spikes = 10000; % this is the convention that should be implemented by MClust PLEASE VERIFY!!!!
Spikes = cell(1,10);
for i = 1:10 % 10 cells 
    firing_rate = exprnd(5); % mean firing rate will be exponentially distributed with a population mean of 5 Hz
    sp = floor(sort(rand(floor(firing_rate * duration), 1) * duration * Fs_spikes)); % these are the timestamps of the spikes 
                                                                                     % in units of samples at Fs_spikes rate
    Spikes{i} = ts(sp * (1. / Fs_spikes)); % spike times in units of seconds, covering a 5 s period
                                    % we encapsulate them in a ts object
                                    % (time series without accompanying
                                    % data) 
end

Spikes = tsdArray(Spikes);



%% let's compute "place fields"

XS = cell(1, 10); % X and Y coordinates of the mouse at the time each spike is emitted
YS = cell(1, 10); 
for i = 1:10 
    XS{i} = Align(X, Spikes{i});
    YS{i} = Align(Y, Spikes{i});
end

cell_ix = 1;
plot(Data(X), Data(Y))
axis equal tight
hold on 
plot(Data(XS{cell_ix}), Data(YS{cell_ix}), 'r.')



%% the mouse speed

XV = timeDeriv(X); % speed on the X and Y component
YV = timeDeriv(Y);

Vel = tsd(Range(XV), sqrt((Data(XV)).^ 2 + (Data(YV).^ 2))); 


