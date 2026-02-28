function a_eci = control_TB_to_ECI(u_tb, R_eci_tb)
%CONTROL_TB_TO_ECI  Convert control acceleration from target body to ECI.
%   a_eci = control_TB_to_ECI(u_tb, R_eci_tb)
%   u_tb     : 3x1 control acceleration in target-body frame [m/s^2]
%   R_eci_tb : 3x3 rotation matrix (ECI from TB)
%   a_eci    : 3x1 control acceleration in ECI [m/s^2]
    a_eci = R_eci_tb * u_tb;
end
