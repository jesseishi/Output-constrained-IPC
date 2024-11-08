function theta_out_rad = MBC_inv(theta_in_rad, psi_rad, psi_offset, B, n)
% Reverse MBC (multiblade coordinate) transformation using azimuth offset
% psi_0.
% Using the convention from muldersAnalysisOptimalIndividual2019.
% Theta_in is [0, tilt, yaw].

% Calculate the azimuth angle of each blade, assuming that they are evenly
% spaced.
psi_1_rad = psi_rad;
psi_2_rad = wrapTo2Pi(psi_rad + (2-1) * 2*pi/B);
psi_3_rad = wrapTo2Pi(psi_rad + (3-1) * 2*pi/B);

% Define the transformation matrix.
T_n_inv = [1 cos(n*(psi_1_rad+psi_offset)) sin(n*(psi_1_rad+psi_offset));
           1 cos(n*(psi_2_rad+psi_offset)) sin(n*(psi_2_rad+psi_offset));
           1 cos(n*(psi_3_rad+psi_offset)) sin(n*(psi_3_rad+psi_offset))];

% Calculate the transformation.
theta_out_rad = T_n_inv * theta_in_rad;
end
