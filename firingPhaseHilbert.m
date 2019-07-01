function [phaseTsd, ph] = firingPhaseHilbert(CR, S) 
% [ph, phaseTsd] = thetaPhaseHilbert(S, CRtheta) compute firing phase of spikes with the Hilbert Transform 
% INPUTS:
% CR: tsd containing EEG filtered in the frequency domain of interest (e.g.
%        theta, or gamma
% S (optional): tsd array of cells, the firing phase of each spike will be calculated 
% OUTPUT: 
% phaseTsd: a tsd of phas efor each of the points in the CR input (defined
% between  0 and 2*pi
% ph: a tsdArray of firing phases for the spike trains in S


% copyright (c) 2006 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
%version 0.1 under construction


zr = hilbert(Data(CR));

phzr = atan2(imag(zr), real(zr));
phzr(phzr < 0) = phzr(phzr < 0) + 2 * pi;
phaseTsd = tsd(Range(CR), phzr);


if nargin > 1
    ph = {};
    for i = 1:length(S)
        ph{i} = Align(phaseTsd, S{i}, 'align', 'closest', 'time', 'align');
    end
    ph = ph(:);
    ph = tsdArray(ph);
end
