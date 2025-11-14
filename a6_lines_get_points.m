% Read "skeleton_connections_s1.mat" from a2_find_all_connections and
% "skeleton_connections_s2.mat", "all_lines_1.mat" from a5_find_all_lines.m
clear;
points = load('skeleton_coordinates_s1.mat').points;
skeleton_connections = load('skeleton_connections_s2.mat').skeleton_connections;
all_lines = load('all_lines_1.mat').all_lines;
all_edges = sort(cat(1, skeleton_connections.end_points), 2);

% Cropped points to remove due to the boundary effect
fcrop = 0.05;
dcrop = (max(points) - min(points)) * fcrop;
bounding_box = [min(points) + dcrop; max(points) - dcrop];

for i = 1: length(all_lines)
    edges = sort(all_lines(i).edges, 2);
    line_points = [];
    len = 0;
    for ii = 1: size(edges, 1)
        eg = edges(ii, :);
        [Lia, Locb] = ismember(eg, all_edges, 'rows');
        eg_idx = find(Lia);
        if length(eg_idx) > 1
            error('Wrong!!!')
        end
        a_idx = Locb(eg_idx);
        line_points = [line_points; skeleton_connections(a_idx).all_points];
        len = len + skeleton_connections(a_idx).distance;
    end
    l_points = points(sort(unique(line_points)), :);
    keep_p = all(l_points >= bounding_box(1,:) & l_points <= bounding_box(2,:), 2);
    l_points = l_points(keep_p, :);
    all_lines(i).points = l_points;
    all_lines(i).length = len;
end

keep_a = arrayfun(@(s) length(s.points) > 3, all_lines);
all_lines = all_lines(keep_a);

% Linearly fitted the lines and save the final found lines as
% "all_lines_last.mat"
std_v = [1; 1; 1];
for i = 1: length(all_lines)
    l_points = all_lines(i).points;
    all_lines(i).mean = mean(l_points);
    [bestline, inliers] = fitLineRANSAC(l_points, 20);
    all_lines(i).p0 = bestline(1: 3);
    if dot(std_v, bestline(4: 6)) > 0
        all_lines(i).vector = bestline(4: 6);
    else
        all_lines(i).vector = -bestline(4: 6);
    end
    all_lines(i).inliers = inliers;
end
save('all_lines_last.mat', 'all_lines');