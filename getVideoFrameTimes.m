function [pos_times, valid_interp_interval] = getVideoFrameTimes(ttl_data, vecbrilho, Fs_ttl, Fs_position)
% finds the times of each video frame in times of the data acquisition system recording LED-related TTL pulses
%
% INPUTS 
% ttl_data: a vector of TTL data reporting the LED electrical state
% vecbrilho: light intensity in each video frame in the LED region 
% Fs_ttl: sampling rate of ttl_data
% Fs_position: sampling rate of vecbriho
% OUTPUTS

video_threshold = 15;
ttl_threshold = 4;

%% Get the moments when the LED on transTTL is on
led_transition_rising_frame_idx = find(diff(vecbrilho > video_threshold) == 1);
led_transition_rising_frame_idx = led_transition_rising_frame_idx(:);

valid_interp_interval = intervalSet(led_transition_rising_frame_idx(1)/Fs_position, ...
                                    led_transition_rising_frame_idx(end)/Fs_position); %#ok<COLND>
                               
ttl_transition_rising_sample_idx=find(diff(ttl_data > ttl_threshold) == 1);
ttl_transition_rising_sample_idx = ttl_transition_rising_sample_idx(:);


%% find which TTL pulse corresponds to the first LED blink 

ttl_time = (ttl_transition_rising_sample_idx - ttl_transition_rising_sample_idx(1)) / Fs_ttl;
video_time = (led_transition_rising_frame_idx - led_transition_rising_frame_idx(1)) / Fs_position;


[X, Y] = meshgrid(0:99, 1:100 );
idx_matrix = X+Y;

align_mat = ttl_time(idx_matrix) - ttl_time(idx_matrix(1, :))';
align_mat2 = align_mat - video_time(1:100);

align_deviation = sum(abs(align_mat2));
[~, first_useful_ttl] = min(align_deviation);

ttl_transition_rising_sample_idx = ttl_transition_rising_sample_idx(first_useful_ttl:end);
led_transition_rising_frame_idx = led_transition_rising_frame_idx(1:length(ttl_transition_rising_sample_idx));

%% perform linear interpolation 

video_frame_times_sec = (1:length(vecbrilho))' / Fs_position;
led_transition_rising_times_sec = video_frame_times_sec(led_transition_rising_frame_idx); %#ok<FNDSB>
ttl_transition_rising_times_sec = ttl_transition_rising_sample_idx / Fs_ttl;

pos_times = interp1(led_transition_rising_times_sec, ...
                    ttl_transition_rising_times_sec, ...
                    video_frame_times_sec, ...
                    'linear', 'extrap');
                
