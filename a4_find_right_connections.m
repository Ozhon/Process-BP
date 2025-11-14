% Read "skeleton_connections_s1.mat" from a2_find_all_connections and
% "skeleton_connections_s1.mat" from a3_simplify_skeleton_network.m
clear;

points = load('skeleton_coordinates_s1.mat').points;
skeleton_connections = load('skeleton_connections_s1.mat').skeleton_connections;
a_idx = arrayfun(@(s) s.chose == 0, skeleton_connections);
skeleton_connections_1 = skeleton_connections(a_idx);
num = length(skeleton_connections_1);

[new_points, slected_path] = extract_lines(points, skeleton_connections_1, num);

% Save flagged edges in the graph. Visualize the used, unused, and deleted 
% edge connections as red, blue, and green color, respectively.
% !!!!!Run this file until all the blue lines are disappeared.!!!!! 
if length(slected_path) > 1
    p_num = length(slected_path);
        discarded = [];
        for i = 1: p_num
            path1 = slected_path(i).path;
            path1 = flip(path1);
            for j = i: p_num
                path2 = slected_path(j).path;
                if isequal(path1, path2)
                    discarded = [discarded, j];
                end
            end
        end
        slected_path(discarded) = [];
    
    for i = 1: length(slected_path)
        path1 = slected_path(i).path;
        path1 = new_points(path1);
        edges = [path1(1: end - 1), path1(2: end)];
        slected_path(i).edges = edges;
    end
    
    all_edges = unique(sort(cat(1, slected_path.edges), 2), 'rows');
    
    all_points = unique(cat(2, slected_path.path));
    all_points = sort(new_points(all_points));
    
    for i = 1: length(skeleton_connections)
        ends = sort(skeleton_connections(i).end_points);
        if sum(ismember(ends, all_points)) >= 1
            skeleton_connections(i).chose = -1;
        end
        if any(ismember(ends, all_edges, 'rows'))
            skeleton_connections(i).chose = 1;
        end
        if length(skeleton_connections(i).p1) > 1 && skeleton_connections(i).distance < 100
            skeleton_connections(i).chose = -1;
        end
    end
    
    save('skeleton_connections_s1.mat', 'skeleton_connections');
    
    a_idx = arrayfun(@(s) s.chose == 1, skeleton_connections);
    end_points = cat(1, skeleton_connections(a_idx).end_points);
    draw_lines(points, end_points, 'r');
    
    a_idx = arrayfun(@(s) s.chose == 0, skeleton_connections);
    end_points = cat(1, skeleton_connections(a_idx).end_points);
    draw_lines(points, end_points, 'b');
else
    for i = 1: length(skeleton_connections)
        if skeleton_connections(i).chose == 0
            skeleton_connections(i).chose = -1;
        end
    end
    save('skeleton_connections_s2.mat', 'skeleton_connections');
    a_idx = arrayfun(@(s) s.chose == 1, skeleton_connections);
    end_points = cat(1, skeleton_connections(a_idx).end_points);
    draw_lines(points, end_points, 'r');
    
    a_idx = arrayfun(@(s) s.chose == -1, skeleton_connections);
    end_points = cat(1, skeleton_connections(a_idx).end_points);
    draw_lines(points, end_points, 'g');
end



function [new_points, slected_path] = extract_lines(points, skeleton_connections_1, num)
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
    
    end_points = cat(1, skeleton_connections_1.new_end_points);
    starts = end_points(:, 1);
    ends = end_points(:, 2);
    weights = cat(1, skeleton_connections_1.distance);
    G = graph(starts, ends, weights);
    
    component_labels = conncomp(G);
    G_sub_n = max(component_labels);
    new_ps = 1: length(new_points);
    slected_path = struct();
    counts = 1;
    for i = 1: G_sub_n
        nodes_sub = new_ps(component_labels == i);
        G_sub = subgraph(G, nodes_sub);
        [~, node_a] = max(distances(G_sub, 1));
        [maxdist, node_b] = max(distances(G_sub, node_a));
        if (maxdist > 100) && (length(nodes_sub) > 2)
            nodes_new = 1: length(nodes_sub);
            nodes_new([node_a, node_b]) = [];
            nodes = [node_a, node_b];
            for ii = 1: 2
                paths = struct();
                start_node = nodes(ii);
                for iii = 1: length(nodes_new)
                    end_node = nodes_new(iii);
                    [path, dist] = shortestpath(G_sub, start_node, end_node);
                    paths(iii).path = path;
                    paths(iii).distance = dist;
                    path = nodes_sub(path);
                    if ~isempty(path)
                        po1 = new_pos(path(1), :);
                        v1 = new_pos(path(end), :) - po1;
                        V = new_pos(path, :) - po1;
                        N = cross(V, repmat(v1, size(V, 1), 1), 2);
                        D = sqrt(sum(N.^2, 2)) / norm(v1);
                        paths(iii).large_dist = max(D);
                    else
                        paths(iii).large_dist = Inf;
                    end
                end
                keep_a = find(arrayfun(@(s) s.large_dist < 50, paths));
                if ~isempty(keep_a)
                    paths = paths(keep_a);
                    dists = cat(1, paths.distance);
                    [~, m_idx] = max(dists);
                    paths = paths(m_idx);
                    if paths(1).distance > 100
                        slected_path(counts).path = nodes_sub(paths.path);
                        slected_path(counts).distance = paths.distance;
                        slected_path(counts).large_dist = paths.large_dist;
                        counts = counts + 1;
                    end
                end
                clear paths;
            end
        end
    end
end