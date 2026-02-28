function y = clamp(x, lo, hi)
%CLAMP  Element-wise clamping of x to [lo, hi].
    y = min(max(x, lo), hi);
end
