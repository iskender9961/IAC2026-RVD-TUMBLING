function R = rotm_from_quat(q)
%ROTM_FROM_QUAT  Rotation matrix from scalar-first quaternion.
%   R = rotm_from_quat(q)
%   q : 4x1 quaternion [qw; qx; qy; qz]
%   R : 3x3 rotation matrix (R rotates body-frame vectors to reference frame)
%       v_ref = R * v_body
    qw = q(1); qx = q(2); qy = q(3); qz = q(4);
    R = [1 - 2*(qy^2+qz^2),   2*(qx*qy - qw*qz),   2*(qx*qz + qw*qy);
         2*(qx*qy + qw*qz),   1 - 2*(qx^2+qz^2),   2*(qy*qz - qw*qx);
         2*(qx*qz - qw*qy),   2*(qy*qz + qw*qx),   1 - 2*(qx^2+qy^2)];
end
