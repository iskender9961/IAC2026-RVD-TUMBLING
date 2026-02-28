function a = accel_J2(r, mu, Re, J2)
%ACCEL_J2  J2 zonal harmonic perturbation acceleration in ECI.
%   a = accel_J2(r, mu, Re, J2)
%   r   : 3x1 ECI position [m]
%   mu  : gravitational parameter [m^3/s^2]
%   Re  : Earth equatorial radius [m]
%   J2  : J2 coefficient [-]
%   a   : 3x1 perturbation acceleration [m/s^2]
    x = r(1); y = r(2); z = r(3);
    rnorm = norm(r);
    r2 = rnorm^2;
    r5 = rnorm^5;
    coeff = -1.5 * J2 * mu * Re^2 / r5;
    z2_r2 = z^2 / r2;
    a = coeff * [ x * (1 - 5*z2_r2);
                  y * (1 - 5*z2_r2);
                  z * (3 - 5*z2_r2) ];
end
