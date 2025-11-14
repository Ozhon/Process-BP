function draw_lines(points, end_points_idx, col)
    starts = points(end_points_idx(:, 1), :);
    ends = points(end_points_idx(:, 2), :);
    line([starts(:,1), ends(:,1)]', [starts(:,2), ends(:,2)]', [starts(:,3), ends(:,3)]', 'Color', col)
    axis equal
end