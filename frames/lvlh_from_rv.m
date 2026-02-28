function R_eci_lvlh = lvlh_from_rv(r_eci, v_eci)
%LVLH_FROM_RV  Rotation matrix from LVLH to ECI given target orbit state.
%   R_eci_lvlh = lvlh_from_rv(r_eci, v_eci)
%
%   LVLH convention:
%     x_L = r/|r|            (radial outward)
%     z_L = h/|h|            (orbit normal = r x v / |r x v|)
%     y_L = z_L x x_L        (approximately along-track for near-circular)
%
%   Returns R_eci_lvlh such that  v_eci = R_eci_lvlh * v_lvlh.
    r_hat = r_eci / norm(r_eci);
    h = cross(r_eci, v_eci);
    z_hat = h / norm(h);
    y_hat = cross(z_hat, r_hat);
    % columns are LVLH basis vectors expressed in ECI
    R_eci_lvlh = [r_hat, y_hat, z_hat];
end
