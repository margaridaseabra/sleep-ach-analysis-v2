function value = value_to_double(inputValue, defaultValue)
%VALUE_TO_DOUBLE Convert a scalar table value to double safely.

if nargin < 2
    defaultValue = NaN;
end

if isnumeric(inputValue) || islogical(inputValue)
    value = double(inputValue);
else
    value = str2double(string(inputValue));
end

if isempty(value) || ismissing(value) || isnan(value)
    value = defaultValue;
end
end
