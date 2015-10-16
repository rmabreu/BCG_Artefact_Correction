% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

function [ px, py ] = parabolic(f, x)

px = (1 / 2) .* (f(x - 1) - f(x + 1)) ./ (f(x - 1) - 2 .* f(x) + f(x + 1)) + x;

py = f(x) - (1 / 4) .* (f(x - 1) - f(x + 1)) .* (px - x);

end