% this is a demo of basic place cell-like analyses, using the TS toolbox.
% The demo is divided in SIMULATE sections, that generate the data and
% COMPUTE section, that exemplifies how you would process real data. 

%% SIMULATE let's make some continuous "data", simulating the firing rate of 
duration = 50;

Fs = 30000; % 30 kHz (sampling rate for open-ephys) 
EEG0_times = (1:(Fs*duration))'; % five seconds of simulated EEG, in units of the sampling rate. Must be a column vector
EEG0_data = rand(size(EEG0_times)); % EEG is white noise 


EEG1_times = EEG0_times;
EEG1_data = sin(2 * pi * 8 * EEG1_times / Fs);


%% INPUT/OUTPUT package the data in a tsd object
EEG0 = tsd(EEG0_times * (1./Fs), EEG0_data); % a tsd object contains the times and the timestamps. 
                                              % the Tsd object expect times
                                              % in units of seconds!
                                              % ts and tsd expect COLUMN
                                              % VECTORS as inputs!
EEG1 = tsd(EEG1_times * (1./Fs), EEG1_data ); % again, note unit conversion for the times (to seconds)

%% COMPUTE plot some EEG 

figure(1)
clf
plot(Range(EEG1), Data(EEG1))
hold on 
% plot(Range(EEG0), Data(EEG0))

%% SIMULATE let's make some "position data"...

% let's have a mouse running on a circular track, very regularly with a
% period of 2 s 
center_arena = [150, 150];

Fs_position = 30;
Pos_times = (1:(Fs_position*duration))'; % five minutes of simulated position data, in units of the sampling rate of the acquisition system
                          % that assumes synchronization and rescaling of
                          % video timestamps! (must be column vector)
                          
X_data = 300 * (0.5 + 0.5 * cos(2 * pi * (1/2) * Pos_times / Fs_position)) + 20 * rand(size(Pos_times));                      
Y_data = 300 * (0.5 + 0.5 * sin(2 * pi * (1/2) * Pos_times / Fs_position)) + 20 * rand(size(Pos_times));                      
 
% here position is given in units of pixels, you may want to convert to
% centimeters at this stage 

%% ...and let's package them in tsd objects 
X = tsd(Pos_times * (1. / Fs_position), X_data);
Y = tsd(Pos_times * (1. / Fs_position), Y_data);


%% COPPUTE plot animal position 

figure(2)
clf
plot(Data(X), Data(Y))
axis equal tight

%% SIMULATE let's make random "spike trains"

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



%% COMPUTE let's compute and draw "place fields"

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



%% COMPUTE the mouse speed

XV = timeDeriv(X); % speed on the X and Y component
YV = timeDeriv(Y);

Vel = tsd(Range(XV), sqrt((Data(XV)).^ 2 + (Data(YV).^ 2))); 

%% SIMULATE let's simulate place cells 
Spikes = cell(1,10);
phi_tsd = tsd(Range(X), atan2(Data(Y) - center_arena(1), Data(X) - center_arena(2))); % extract the angular coordinate from the X/Y position data

t = Range(phi_tsd);
phi = Data(phi_tsd);
for i = 1:10
    pf_centre = rand(1,1) * 2 * pi;
    pf_size = exprnd(0.5);
    kappa = 1 / pf_size;
    peak_firing_rate = exprnd(10) ;
    firing_prob = peak_firing_rate * (1 / Fs_position) * exp(kappa * cos(phi-pf_centre)) / exp(kappa); % probability of firing a von Mises function of the angular deviation from 
                                                                       % place
                                                                       % field
                                                                       % centre
    sp = rand(size(phi)) < firing_prob;
    Spikes{i} = ts(t(sp));
end

Spikes = tsdArray(Spikes);

%% COMPUTE let's compute and draw "place fields"

XS = cell(1, 10); % X and Y coordinates of the mouse at the time each spike is emitted
YS = cell(1, 10); 
for i = 1:10 
    XS{i} = Align(X, Spikes{i});
    YS{i} = Align(Y, Spikes{i});
end

cell_ix = 10;
figure(2)
clf
plot(Data(X), Data(Y))
axis equal tight
hold on 
for cell_ix = 1:10
    plot(Data(XS{cell_ix}), Data(YS{cell_ix}), '.', 'MarkerSize', 15)
end


%% SIMULATE for a more realistic simulation: 
% speed is now variable and random 
duration = 200;
% let's have a mouse running on a circular track, very regularly with a
% period of 2 s 
center_arena = [150, 150];

Fs_position = 30;
Pos_times = (1:(Fs_position*duration))' / Fs_position; % five minutes of simulated position data, in units of the sampling rate of the acquisition system
                          % that assumes synchronization and rescaling of              
                        % video timestamps! (must be column vector)
Vel_factor = 0.01;                        
v = Vel_factor * (rand(size(Pos_times))-0.49);
v = conv(v, hanning(5 * Fs_position), 'same');
v(v < 0) = 0;



phi = cumsum(v);

X_data = 300 * (0.5 + 0.5 * cos(phi)) + 20 * rand(size(Pos_times));                      
Y_data = 300 * (0.5 + 0.5 * sin(phi)) + 20 * rand(size(Pos_times));                      
 
figure(5), clf
plot(Pos_times, X_data)
hold on
plot(Pos_times, Y_data)

X = tsd(Pos_times, X_data);
Y = tsd(Pos_times, Y_data);

%% SIMULATE let's simulate place cells  
% we pretend we don't know how X and Y were simulated 
Spikes = cell(1,10);
phi_tsd = tsd(Range(X), atan2(Data(Y) - center_arena(1), Data(X) - center_arena(2))); % extract the angular coordinate from the X/Y position data

t = Range(phi_tsd);
phi = Data(phi_tsd);
XV = timeDeriv(smooth(X, 10)); % speed on the X and Y component
YV = timeDeriv(smooth(Y, 10));

Vel = tsd(Range(XV), sqrt((Data(XV)).^ 2 + (Data(YV).^ 2))); 

not_running = Data(Vel) < 12;
not_running(end+1) = true;
not_running_FR = 3;

for i = 1:10
    pf_centre = rand(1,1) * 2 * pi;
    pf_size = exprnd(0.1);
    kappa = 1 / pf_size;
    peak_firing_rate = exprnd(20) ;
    firing_prob = (~not_running) .* (peak_firing_rate * (1 / Fs_position) * exp(kappa * cos(phi-pf_centre)) / exp(kappa)); % probability of firing a von Mises function of the angular deviation from 
                                                                       % place
                                                                       % field
                                                                       % centre
    sp = find(rand(size(phi)) < firing_prob);
    
    sp_non_running = find(rand(size(phi)) < not_running * not_running_FR / Fs_position);
   
    sp = sort([sp ; sp_non_running]);
    Spikes{i} = ts(t(sp));
end

Spikes = tsdArray(Spikes);

%% COMPUTE let's compute and draw "place fields"

XS = cell(1, 10); % X and Y coordinates of the mouse at the time each spike is emitted
YS = cell(1, 10); 
for i = 1:10 
    XS{i} = Align(X, Spikes{i});
    YS{i} = Align(Y, Spikes{i});
end

figure(2)
subplot(1, 2, 1)
cla
plot(Data(X), Data(Y))
axis equal tight
hold on 
for cell_ix = 1:10
    plot(Data(XS{cell_ix}), Data(YS{cell_ix}), '.', 'MarkerSize', 15)
end


% in the second subplot we will filter the spikes fired while the mouse
% wasn't running 

XV = timeDeriv(smooth(X, 10)); % speed on the X and Y component
YV = timeDeriv(smooth(Y, 10));

Vel = tsd(Range(XV), sqrt((Data(XV)).^ 2 + (Data(YV).^ 2))); 

running = thresholdIntervals(Vel, 12, 'Direction', 'Above');

XS = cell(1, 10); % X and Y coordinates of the mouse at the time each spike is emitted
YS = cell(1, 10); 
for i = 1:10 
    spikes_running = Restrict(Spikes{i}, running);
    XS{i} = Align(X, spikes_running);
    YS{i} = Align(Y, spikes_running);
end

figure(2)
subplot(1, 2, 2)
cla
plot(Data(X), Data(Y))
axis equal tight
hold on 
for cell_ix = 1:10
    plot(Data(XS{cell_ix}), Data(YS{cell_ix}), '.', 'MarkerSize', 15)
end



%% SIMULATE let's get some "theta" so that we can compute phase precession 


Fs = 30000; % 30 kHz (sampling rate for open-ephys) 
EEG0_times = (1:(Fs*duration))' / Fs; % five seconds of simulated EEG, in units of the sampling rate. Must be a column vector


EEG1_times = EEG0_times;


EEG1_phase = 2 * pi * 8.1 * EEG1_times + 1.2 * rand(size(EEG1_times)); 
EEG1_data = cos(EEG1_phase) ;

EEG1 = tsd(EEG1_times, EEG1_data); 
EEG1_phase = tsd(EEG1_times, EEG1_phase);

%% SIMULATE let's make phase precessing cells 

Spikes = cell(1,10);
phi_tsd = tsd(Range(X), atan2(Data(Y) - center_arena(1), Data(X) - center_arena(2))); % extract the angular coordinate from the X/Y position data

t = Range(phi_tsd);
phi = Data(phi_tsd);

theta_phase = Align(EEG1_phase, Range(X));
theta_phase = Data(theta_phase);
kappa_theta = 3;
phase_precession_center = pi;
phase_precession_alpha = -3;
for i = 1:10
    pf_centre = rand(1,1) * 2 * pi;
    pf_size = exprnd(0.1);
    kappa = 1 / pf_size;
    peak_firing_rate = exprnd(30) ;
    firing_prob = peak_firing_rate * (1 / Fs_position) * ...
        exp(kappa * cos(phi-pf_centre) + ...
        kappa_theta * cos(theta_phase - phase_precession_center - phase_precession_alpha*(phi-pf_centre))) / ...
        exp(kappa + kappa_theta); % probability of firing a von Mises function of the angular deviation from 
                                                                       % place
                                                                       % field
                                                                       % centre
    sp = rand(size(phi)) < firing_prob;
    Spikes{i} = ts(t(sp));
end

Spikes = tsdArray(Spikes);

%% COMPUTE now let's compute phase precession

resample_rate = 30;
t = Range(EEG1);
d = Data(EEG1);
t_resample = t(1:resample_rate:end);
d_resample = resample(d, 1, resample_rate);

EEG1_res = tsd(t_resample, d_resample);


EEG1_theta = filter_7hz(EEG1_res);

[phaseTsd, ph] = firingPhaseHilbert(EEG1_theta, Spikes);


%%  COMPUTE phase precession plots
phi_tsd = tsd(Range(X), atan2(Data(Y) - center_arena(1), Data(X) - center_arena(2))); % extract the angular coordinate from the X/Y position data


for i = 1:10
    
    phiS = Align(phi_tsd, Spikes{i});
    figure(8), clf
    plot(Data(phiS), Data(ph{i}), '.');
    axis([-pi, pi, 0, 2*pi]);
    xlabel('position')
    ylabel('theta phase')
    
    
    title('type dbcont to continue');
    keyboard
end


