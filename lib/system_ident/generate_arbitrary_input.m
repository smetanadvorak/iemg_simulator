function [ profile ] = generate_identification_input( len, pause, levels, amplitude )

function [ profile ] = generate_arbitrary_input( len, pause, levels, amplitude )

if nargin < 
% Part 1: rectangles
T1 = floor(T/2);
part1 = zeros(T1,1);
part1(floor(T1/6):end-floor(T1/6)) = 1;

% Part 2: trapezoids
T2 = T - T1;
part2 = sin(linspace(0,pi,T2))';
part2 = part2 + 0.2*sin(9*linspace(0,pi,T2))';
part2 = part2/max(part2);
profile = [part1; part2];

% Part 3: slow variations around set level


end

