% Read "skeleton_connections.mat" and "skeleton_coordinates_s1.mat" from
% a2_find_all_connections.m
clear;
filename = 'skeleton_connections.mat';
points = load('skeleton_coordinates_s1.mat').points;
skeleton_connections = load(filename).skeleton_connections;

% Save flagged edges in the graph
skeleton_connections = find_long_lines(points, skeleton_connections);
save('skeleton_connections_s1.mat', 'skeleton_connections');


% Visualize the used, unused, and deleted edge connections as red, blue,
% and green color, respectively
a_idx = arrayfun(@(s) s.chose == 1, skeleton_connections);
end_points = cat(1, skeleton_connections(a_idx).end_points);
starts = points(end_points(:, 1), :);
ends = points(end_points(:, 2), :);
line([starts(:,1), ends(:,1)]', [starts(:,2), ends(:,2)]', [starts(:,3), ends(:,3)]', 'Color', 'r')
axis equal

a_idx = arrayfun(@(s) s.chose == 0, skeleton_connections);
end_points = cat(1, skeleton_connections(a_idx).end_points);
starts = points(end_points(:, 1), :);
ends = points(end_points(:, 2), :);
line([starts(:,1), ends(:,1)]', [starts(:,2), ends(:,2)]', [starts(:,3), ends(:,3)]', 'Color', 'b')
axis equal

a_idx = arrayfun(@(s) s.chose == -1, skeleton_connections);
end_points = cat(1, skeleton_connections(a_idx).end_points);
starts = points(end_points(:, 1), :);
ends = points(end_points(:, 2), :);
line([starts(:,1), ends(:,1)]', [starts(:,2), ends(:,2)]', [starts(:,3), ends(:,3)]', 'Color', 'g')
axis equal



function D = points_to_lines(point1, vector1, points2)
    V = points2 - point1;
    N = cross(V, repmat(vector1, size(V, 1), 1), 2);
    D = sqrt(sum(N.^2, 2)) / norm(vector1);
end

function skeleton_connections = find_long_lines(points, skeleton_connections)
    a_idx = arrayfun(@(s) s.chose ~= 1, skeleton_connections);
    skeleton_tem = skeleton_connections(a_idx);
    num = length(skeleton_tem);
    
    node_points = reshape(cat(1, skeleton_tem.end_points), [num * 2, 1]);
    new_points = sort(unique(node_points));
    new_ps = points(new_points, :);
    num_points = size(points, 1);
    index_map = zeros(num_points, 1);
    index_map(new_points) = 1:length(new_points);
    
    fcrop = 0.05;
    dcrop = (max(points) - min(points)) * fcrop;
    bounding_box = [min(points), max(points)];
    p1 = cat(2, skeleton_tem.p1);
    p1p = points(p1, :);
    keep_a = abs(repmat(p1p, [1, 2]) - repmat(bounding_box, [length(p1), 1])) < repmat(dcrop, [1, 2]);
    p1 = p1(any(keep_a, 2));
    keep_a = keep_a(any(keep_a, 2), :);
    p1_new = index_map(p1);
    p1_new_sort = sort(p1_new);
    
    end_points = cat(1, skeleton_tem.new_end_points);
    starts = end_points(:, 1);
    ends = end_points(:, 2);
    weights = cat(1, skeleton_tem.distance);
    G = graph(starts, ends, weights);
    
    right_paths = struct();
    counts = 1;
    for i = 1: length(p1_new)
        face_group = find(keep_a(i, :) == 1);
        face_group = face_group(1);
        start_node = p1_new(i);
        end_nodes = p1_new(keep_a(:, face_group) ~= 1);
        paths = struct();
        for ii = 1: length(end_nodes)
            [path, dist] = shortestpath(G, start_node, end_nodes(ii));
            paths(ii).path = path;
            paths(ii).distance = dist;
            if (sum(ismember(path, p1_new_sort)) > 2) && (dist < 350)
                path = [];
            end
            if ~isempty(path)
                po1 = new_ps(path(1), :);
                v1 = new_ps(path(end), :) - po1;
                paths(ii).large_dist = max(points_to_lines(po1, v1, new_ps(path, :)));
            else
                paths(ii).large_dist = Inf;
            end
        end
        linear_path = find(arrayfun(@(s) s.large_dist < 50, paths));
        if ~isempty(linear_path)
            for iii = 1: length(linear_path)
                path1 = paths(linear_path(iii)).path;
                right_paths(counts).path = path1;
                right_paths(counts).distance = paths(linear_path(iii)).distance;
                path1 = new_points(path1);
                edges = [path1(1: end - 1), path1(2: end)];
                right_paths(counts).edges = edges;
                counts = counts + 1;
            end
        end
        clear paths;
    end
    
    p_num = length(right_paths);
    discarded = zeros(p_num);
    for i = 1: p_num
        path1 = right_paths(i).path;
        path1 = flip(path1);
        for j = i: p_num
            path2 = right_paths(j).path;
            if isequal(path1, path2)
                discarded(j) = 1;
            end
        end
    end
    right_paths(discarded == 1) = [];
    
    all_edges = unique(sort(cat(1, right_paths.edges), 2), 'rows');
    all_points = unique(cat(2, right_paths.path));
    all_points = sort(new_points(all_points));
    
    for i = 1: num
        ends = sort(skeleton_connections(i).end_points);
        if sum(ismember(ends, all_points)) >= 1
            skeleton_connections(i).chose = -1;
        else
            skeleton_connections(i).chose = 0;
        end
        if any(ismember(ends, all_edges, 'rows'))
            skeleton_connections(i).chose = 1;
        end
    end
end