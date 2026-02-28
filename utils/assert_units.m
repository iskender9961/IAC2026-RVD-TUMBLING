function assert_units(value, expected_order, name)
%ASSERT_UNITS  Sanity-check that a scalar's order of magnitude is plausible.
%   assert_units(value, expected_order, name)
%   Warns if log10(abs(value)) differs from expected_order by more than 2.
    if value == 0
        return;
    end
    actual = log10(abs(value));
    if abs(actual - expected_order) > 2
        warning('assert_units:mismatch', ...
            '%s = %.4e  (order %.1f, expected ~%.1f). Check units.', ...
            name, value, actual, expected_order);
    end
end
