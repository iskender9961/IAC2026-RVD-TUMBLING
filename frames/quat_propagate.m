function q_next = quat_propagate(q, omega_body, dt)
%QUAT_PROPAGATE  Propagate attitude quaternion by one step assuming constant body rate.
%   q_next = quat_propagate(q, omega_body, dt)
%   q          : 4x1 quaternion [qw; qx; qy; qz] (scalar-first)
%   omega_body : 3x1 body angular velocity [rad/s]
%   dt         : time step [s]
%   q_next     : 4x1 propagated quaternion (normalized)
%
%   Uses exact exponential integration for constant omega over dt.
    w = omega_body;
    wnorm = norm(w);
    if wnorm < 1e-14
        q_next = q;
        return;
    end
    half_angle = 0.5 * wnorm * dt;
    e = w / wnorm;     % rotation axis
    % quaternion for the rotation increment
    dq = [cos(half_angle); sin(half_angle) * e];
    % quaternion multiply:  q_next = q * dq  (Hamilton convention)
    q_next = quat_mult(q, dq);
    q_next = q_next / norm(q_next);
end

function qr = quat_mult(q, p)
%QUAT_MULT  Hamilton quaternion multiplication (scalar-first).
    qw = q(1); qv = q(2:4);
    pw = p(1); pv = p(2:4);
    qr = [qw*pw - dot(qv,pv);
          qw*pv + pw*qv + cross(qv,pv)];
end
