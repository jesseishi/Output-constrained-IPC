function [optimal_azimuth_offset_rad, dMtdThetat, dMydThetay] = getOptimalAzimuthOffset(windSpeed_mps)
% getOptimalAzimuthOffset
% Read the optimal azimuth offset from a precompiled table.

T = readtable('optimal_azimuth_offset_fullDOF.csv');

k = find(T.wind_speed_mps == windSpeed_mps);
if ~isempty(k)
    optimal_azimuth_offset_deg = T.psi_o_deg(k);
    dMtdThetat = T.dMtdThetat(k);
    dMydThetay = T.dMydThetay(k);
else
    optimal_azimuth_offset_deg = interp1(T.wind_speed_mps, T.psi_o, windSpeed_mps);
    dMtdThetat = interp1(T.wind_speed_mps, T.dMtdThetat, windSpeed_mps);
    dMydThetay = interp1(T.wind_speed_mps, T.dMydThetay, windSpeed_mps);
end

% Convert to the right unit.
optimal_azimuth_offset_rad = deg2rad(optimal_azimuth_offset_deg);

end
