function [M_out_Nm] = MBC(M_in_Nm, psi_rad, psi_offset_rad, B, n)
% MBC (multiblade coordinate) transformation
% Using the convention from muldersAnalysisOptimalIndividual2019.

% Calculate the azimuth angle of each blade, assuming that they are evenly
% spaced.
% TODO: If I make this code open-source it'd be nice to not depend on the
% Mapping toolbox just for the wrapTo2Pi function.
psi_1_rad = psi_rad + psi_offset_rad;
psi_2_rad = wrapTo2Pi(psi_rad + (2-1) * 2*pi/B + psi_offset_rad);
psi_3_rad = wrapTo2Pi(psi_rad + (3-1) * 2*pi/B + psi_offset_rad);

% Define the transformation matrix.
% TODO: is the top row 1/2 (like in his PhD) or 1 (like in the journal
% paper)?
T_n = 2/B * [1/2          1/2          1/2;
             cos(n*psi_1_rad) cos(n*psi_2_rad) cos(n*psi_3_rad);
             sin(n*psi_1_rad) sin(n*psi_2_rad) sin(n*psi_3_rad)];

% Calculate the transformation.
M_out_Nm = T_n * M_in_Nm;
end
