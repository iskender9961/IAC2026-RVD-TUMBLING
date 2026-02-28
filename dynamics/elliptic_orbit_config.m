function cfg = elliptic_orbit_config(ecc, a_km, mu)
%ELLIPTIC_ORBIT_CONFIG  Build orbital config struct for YA dynamics.
%
%   cfg = elliptic_orbit_config(ecc, a_km)
%   cfg = elliptic_orbit_config(ecc, a_km, mu)
%
%   Inputs:
%     ecc    - orbital eccentricity (0 for circular)
%     a_km   - semi-major axis [km]
%     mu     - gravitational parameter [m^3/s^2] (default: Earth)
%
%   Output:
%     cfg    - struct compatible with YA_A_matrix.m

if nargin < 3
    mu = 3.986004418e14;   % Earth [m^3/s^2]
end

cfg.mu   = mu;
cfg.ecc  = ecc;
cfg.ecc2 = ecc^2;
cfg.a    = a_km * 1e3;            % [m]
cfg.p    = cfg.a * (1 - cfg.ecc2);  % semi-latus rectum [m]
cfg.n    = sqrt(mu / cfg.a^3);    % mean motion [rad/s]
cfg.rp   = cfg.p / (1 + ecc);     % periapsis radius [m]
cfg.v0   = sqrt(2*mu*(1/cfg.rp - 1/(2*cfg.a)));  % speed at periapsis
cfg.h    = cfg.rp * cfg.v0;       % specific angular momentum
cfg.ya_k2 = cfg.h / cfg.p^2;     % YA scaling constant

% Precomputed quantities for nu->M conversion
cfg.X1a = [sqrt(1 - ecc); sqrt(1 + ecc)];
cfg.X2  = ecc * sqrt(1 - cfg.ecc2);
end
