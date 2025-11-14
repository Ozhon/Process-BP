% !!!!!!!!Attention!!!!!!
% The following three parts of the codes need to be run separately and 
% cannot be run simultaneously

%% Part 1
% Read "skeleton_connections_s1.mat" from a2_find_all_connections and
% "skeleton_connections_s1.mat" from a4_find_right_connections.m
% !!!!!!!!Run this part once!!!!
clear;
points = load('skeleton_coordinates_s1.mat').points;
skeleton_connections = load('skeleton_connections_s2.mat').skeleton_connections;

a_idx = find(arrayfun(@(s) s.chose == 1, skeleton_connections));
skeleton_connections_1 = skeleton_connections(a_idx);
num = length(skeleton_connections_1);

% Find the long linear disclinatiion lines in labelled edges. Save flagged
% edges as "skeleton_connections_s2.mat" and the disclination lines as
% "all_lines_1.mat"
[all_lines, skeleton_connections] = find_straight_lines(points, skeleton_connections, a_idx, skeleton_connections_1, num);
save('all_lines_1.mat', 'all_lines')
save('skeleton_connections_s2.mat', 'skeleton_connections')

%% Part 2
% Read "skeleton_connections_s1.mat" and "all_lines_1.mat" from Part 1
clear;
skeleton_connections = load('skeleton_connections_s2.mat').skeleton_connections;

% Find the residual disclination lines in the sub-graph. Save flagged
% edges as "skeleton_connections_s2.mat". Add new found disclination lines
% to "all_lines_1.mat".
% !!!!!!Run this part until all the blue lines are disappeared!!!!!
a_idx = find(arrayfun(@(s) s.chose == 1, skeleton_connections));
skeleton_connections_1 = skeleton_connections(a_idx);
num = length(skeleton_connections_1);
all_lines = load("all_lines_1.mat").all_lines;
[ex_lines, skeleton_connections] = find_straight_lines(points, skeleton_connections, a_idx, skeleton_connections_1, num);
all_lines = [all_lines, ex_lines];
save('all_lines_1.mat', 'all_lines')
save('skeleton_connections_s2.mat', 'skeleton_connections')

% Visualize the used, unused, and deleted edge connections as red, blue, 
% and green color, respectively
a_idx = arrayfun(@(s) s.chose == 2, skeleton_connections);
end_points = cat(1, skeleton_connections(a_idx).end_points);
draw_lines(points, end_points, 'r');

a_idx = arrayfun(@(s) s.chose == 1, skeleton_connections);
end_points = cat(1, skeleton_connections(a_idx).end_points);
draw_lines(points, end_points, 'b');

a_idx = arrayfun(@(s) s.chose == -1, skeleton_connections);
end_points = cat(1, skeleton_connections(a_idx).end_points);
draw_lines(points, end_points, 'g');

%% Part 3
% Read "skeleton_connections_s1.mat" and "all_lines_1.mat" from Part 2
% !!!!!!!!Run this part once!!!!
clear;
skeleton_connections = load('skeleton_connections_s2.mat').skeleton_connections;

% Find the residual disclination lines in deleted edges connection. Save 
% flagged edges as "skeleton_connections_s2.mat". Add new found 
% disclination lines to "all_lines_1.mat"
a_idx = find(arrayfun(@(s) s.chose == 1, skeleton_connections));
skeleton_connections_1 = skeleton_connections(a_idx);
num = length(skeleton_connections_1);
all_lines = load("all_lines_1.mat").all_lines;
a_idx = find(arrayfun(@(s) s.chose == -1, skeleton_connections));
skeleton_connections_1 = skeleton_connections(a_idx);
num = length(skeleton_connections_1);
ex_lines = struct();
count = 1;
for i = 1: num
    if skeleton_connections_1(i).distance > 120 && length(skeleton_connections_1(i).all_points) > 15
        ex_lines(count).path = skeleton_connections_1(i).end_points;
        ex_lines(count).edges = skeleton_connections_1(i).end_points;
        count = count + 1;
        skeleton_connections(a_idx(i)).chose = 2;
    end
end
all_lines = [all_lines, ex_lines];
save('all_lines_1.mat', 'all_lines')
save('skeleton_connections_s2.mat', 'skeleton_connections')


function D = points_to_lines(point1, vector1, points2)
    V = points2 - point1;
    N = cross(V, repmat(vector1, size(V, 1), 1), 2);
    D = sqrt(sum(N.^2, 2)) / norm(vector1);
end

function right_paths = find_linear_from_path(G_sub, nodes_sub, node_a, node_b, new_pos)
    nodes_new = 1: length(nodes_sub);
    nodes_new([node_a, node_b]) = [];
    nodes = [node_a, node_b];
    right_paths = struct();
    for i = 1: 2
        paths = struct();
        start_node = nodes(i);
        for ii = 1: length(nodes_new)
            end_node = nodes_new(ii);
            [path, dist] = shortestpath(G_sub, start_node, end_node);
            paths(ii).path = path;
            paths(ii).distance = dist;
            path = nodes_sub(path);
            if ~isempty(path)
                po1 = new_pos(path(1), :);
                v1 = new_pos(path(end), :) - po1;
                paths(ii).large_dist = max(points_to_lines(po1, v1, new_pos(path, :)));
            end
        end
        paths = paths(arrayfun(@(s) s.large_dist < 50, paths));
        dists = cat(1, paths.distance);
        [~, m_idx] = max(dists);
        right_paths(i).path = paths(m_idx).path;
        right_paths(i).distance = paths(m_idx).distance;
        clear paths;
    end
end

function [all_lines, skeleton_connections] = find_straight_lines(points, skeleton_connections, a_idx, skeleton_connections_1, num)
    end_points = reshape(cat(1, skeleton_connections_1.end_points), [num * 2, 1]);
    new_points = sort(unique(end_points));
    new_pos = points(new_points, :);
    num_points = size(points, 1);
    index_map = zeros(num_points, 1);
    index_map(new_points) = 1:length(new_points);
    for i = 1: num
        ends = skeleton_connections_1(i).end_points;
        skeleton_connections_1(i).new_end_points = index_map(ends)';
    end
    
    % 构建图，找到所有的联通区域 % 这里也可改成函数
    end_points = cat(1, skeleton_connections_1.new_end_points);
    starts = end_points(:, 1);
    ends = end_points(:, 2);
    weights = cat(1, skeleton_connections_1.distance);
    G = graph(starts, ends, weights);
    
    degree3 = find(degree(G) == 3);
    component_labels = conncomp(G);
    G_sub_n = max(component_labels);
    
    all_lines = struct();
    count = 1;
    for i = 1: G_sub_n
        nodes_sub = find(component_labels == i);  % 可以取出子图，上面的部分代码也可修改，所有的重新索引都可以重新写
        G_sub = subgraph(G, nodes_sub); % 注意又一次重新索引了 后续可以修改一下，边部的连接比较奇怪
        [~, node_a] = max(distances(G_sub, 1));
        [~, node_b] = max(distances(G_sub, node_a));
        [path, ~] = shortestpath(G_sub, node_a, node_b);
        path = nodes_sub(path);
        v1 = new_pos(nodes_sub(node_a), :) - new_pos(nodes_sub(node_b), :);
        dists = points_to_lines(new_pos(nodes_sub(node_a), :), v1, new_pos(path, :));
        if max(dists) > 50
            path_l = find_linear_from_path(G_sub, nodes_sub, node_a, node_b, new_pos);
            all_lines(count).path = nodes_sub(path_l(1).path);
            all_lines(count + 1).path = nodes_sub(path_l(2).path);
            count = count + 2;
            clear path_l;
        else
            all_lines(count).path = path;
            count = count + 1;
        end
    end
    
    for i = 1: length(all_lines)
        path1 = all_lines(i).path;
        path1 = new_points(path1);
        edges = [path1(1: end - 1), path1(2: end)];
        all_lines(i).edges = edges;
    end
    
    all_edges = unique(sort(cat(1, all_lines.edges), 2), 'rows');
    
    for i = 1: length(skeleton_connections_1)
        ends = sort(skeleton_connections_1(i).end_points);
        if any(ismember(ends, all_edges, 'rows'))
            skeleton_connections(a_idx(i)).chose = 2;
        end
    end
end

function [all_lines, skeleton_connections] = find_straight_lines_2(points, skeleton_connections, a_idx, skeleton_connections_1, num)
    end_points = reshape(cat(1, skeleton_connections_1.end_points), [num * 2, 1]);
    new_points = sort(unique(end_points));
    new_pos = points(new_points, :);
    num_points = size(points, 1);
    index_map = zeros(num_points, 1);
    index_map(new_points) = 1:length(new_points);
    for i = 1: num
        ends = skeleton_connections_1(i).end_points;
        skeleton_connections_1(i).new_end_points = index_map(ends)';
    end
    
    % 构建图，找到所有的联通区域 % 这里也可改成函数
    end_points = cat(1, skeleton_connections_1.new_end_points);
    starts = end_points(:, 1);
    ends = end_points(:, 2);
    weights = cat(1, skeleton_connections_1.distance);
    G = graph(starts, ends, weights);
    
    degree3 = find(degree(G) == 3);
    component_labels = conncomp(G);
    G_sub_n = max(component_labels);
    
    all_lines = struct();
    count = 1;
    for i = 1: G_sub_n
        nodes_sub = find(component_labels == i);  % 可以取出子图，上面的部分代码也可修改，所有的重新索引都可以重新写
        G_sub = subgraph(G, nodes_sub); % 注意又一次重新索引了 后续可以修改一下，边部的连接比较奇怪
        [~, node_a] = max(distances(G_sub, 1));
        [~, node_b] = max(distances(G_sub, node_a));
        [path, dist] = shortestpath(G_sub, node_a, node_b);
        if dist > 300
            path = nodes_sub(path);
            v1 = new_pos(nodes_sub(node_a), :) - new_pos(nodes_sub(node_b), :);
            dists = points_to_lines(new_pos(nodes_sub(node_a), :), v1, new_pos(path, :));
            if max(dists) > 50
                path_l = find_linear_from_path(G_sub, nodes_sub, node_a, node_b, new_pos);
                if path_l(1).distance > 250
                    all_lines(count).path = nodes_sub(path_l(1).path);
                    count = count + 1;
                end
                if path_l(2).distance > 250
                    all_lines(count).path = nodes_sub(path_l(2).path);
                    count = count + 1;
                end
                clear path_l;
            else
                all_lines(count).path = path;
                count = count + 1;
            end
        end
    end

    for i = 1: length(all_lines)
        path1 = all_lines(i).path;
        path1 = new_points(path1);
        edges = [path1(1: end - 1), path1(2: end)];
        all_lines(i).edges = edges;
    end

    all_edges = unique(sort(cat(1, all_lines.edges), 2), 'rows');
    
    for i = 1: length(skeleton_connections_1)
        ends = sort(skeleton_connections_1(i).end_points);
        if any(ismember(ends, all_edges, 'rows'))
            skeleton_connections(a_idx(i)).chose = 2;
        end
    end
end