% plotRealSearchg.m
% Plots objective function g(x) = ceil(log_2(x)) for x in [1/16, 16] for Figure 3(a)

linethickness = 2.5;

figure
axis([-0.5, 16.5, -4.5, 4.5]);
line([1/16, 1/16], [-4, -4], 'LineWidth', linethickness, 'Color', 'k');
line([1/16, 1/8], [-3, -3], 'LineWidth', linethickness, 'Color', 'k');
line([1/8, 1/4], [-2, -2], 'LineWidth', linethickness, 'Color', 'k');
line([1/4, 1/2], [-1, -1], 'LineWidth', linethickness, 'Color', 'k');
line([1/2, 1], [0, 0], 'LineWidth', linethickness, 'Color', 'k');
line([1, 2], [1, 1], 'LineWidth', linethickness, 'Color', 'k');
line([2, 4], [2, 2], 'LineWidth', linethickness, 'Color', 'k');
line([4, 8], [3, 3], 'LineWidth', linethickness, 'Color', 'k');
line([8, 16], [4, 4], 'LineWidth', linethickness, 'Color', 'k');

xlabel('System (x)', 'FontSize', 14);
ylabel('True Performance g(x)', 'FontSize', 14);
title('Realistic Example', 'FontSize', 14);