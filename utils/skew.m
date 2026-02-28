function S = skew(w)
%SKEW  Skew-symmetric matrix from 3-vector.
%   S = skew(w)  returns the 3x3 skew-symmetric matrix such that
%   S*v = cross(w,v) for any 3-vector v.
    S = [  0   -w(3)  w(2);
          w(3)   0   -w(1);
         -w(2)  w(1)   0  ];
end
