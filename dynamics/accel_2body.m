function a = accel_2body(r, mu)
%ACCEL_2BODY  Two-body gravitational acceleration in ECI.
%   a = accel_2body(r, mu)
%   r  : 3x1 position vector [m]
%   mu : gravitational parameter [m^3/s^2]
%   a  : 3x1 acceleration [m/s^2]
    rnorm = norm(r);
    a = -mu / rnorm^3 * r;
end
