function Y = Filter_7hz(IN_tsd, lowlimit, highlimit)
%function Y = Filter_7hz(IN_tsd)
% Filter for Theta EEG
% 
% INPUT:
%        IN_tsd     = The EEG data to be filtered (a tsd object).
%        lowlimit   = lower limit of the band pass filter 
%                     (default is 6 Hz)
%        highlimit  = upper limit of the band pass filter 
%                     (default is 10 Hz) Values above  Hz don't work
% 
% OUTPUT: 
%        Y        = tsd of filtered data that is subsampled.

% cowen Wed Mar 24 14:08:41 1999

sFreq = 200; % Hz. The desired sample frequency. It can't be too large for theta.


input_sFreq = 1./median(diff(Range(IN_tsd)))

if input_sFreq > sFreq % only subsample if the original sampling rate greater than desired one
    interval = floor(input_sFreq/sFreq);
else
    interval = 1;
end

crdata = Data(IN_tsd);
crtt = Range(IN_tsd);

if nargin == 1
  % Default band pass
  lowlimit  = 6;  % Hz
  highlimit = 10; % Hz
end
if highlimit > 50
  error('Keep the upper limit below 50. The filter cannot handle larger valuse')
end


% Filter parameters
F_Ny = sFreq/2;           % Hz
lowcut = lowlimit/F_Ny;   % Hz
highcut = highlimit/F_Ny; % Hz
N = 4;             % Order of the filter
passband = [lowcut highcut];  
ripple = .1;

% Make sure the Cr tsd is OK.
if length(crdata) ~= length(crtt)
  error('input is misaligned');
end



%[Bb,Ab] = butter(N, passband);
[Bc,Ac] = cheby1(N, ripple, passband);
%h = [abs(hh) abs(freqz(Bb,Ab,n)) abs(freqz(Bc,Ac,n))];
% Use filtfilt instead of filter because it corrects for phase
% distortion due to filtering. It runs through the filter twice.
F = filtfilt(Bc,Ac,crdata(1:interval:end));

% The interval is used for subsampling the data
Y = tsd(crtt(1:interval:end),F);
