% Read data from a1_skeleton.m
clear;
filename = 'skeleton_coordinates.mat';
points = load(filename).coordinates;

connections = find_nearest_n(points);
points_n = length(connections);
point3 = find(arrayfun(@(s) length(s.nn) >= 3, connections));

[new_edges, discarded_edges, points, point3] = reconstruct_p3_connections(points, connections, point3, points_n);
all_edges = construct_edge_from_connections(connections, 1: points_n);
all_edges = setdiff(all_edges, discarded_edges, 'rows');
all_edges = [all_edges; new_edges];

conn_map = connect_matrix(all_edges, points_n);
nodes_degree_3 = find(sum(conn_map, 2) >= 3);
nodes_degree_1 = find(sum(conn_map, 2) == 1);
n3 = length(nodes_degree_3);
n1 = length(nodes_degree_1);
totol_n = 1;
skeleton_connections = struct();

for i = 1: n3
    point3_1 = nodes_degree_3(i);
    point3_1n = find(conn_map(point3_1, :) == 1);
    nnn = length(point3_1n);

    for ii = 1: nnn
        l_conns = bfs(conn_map, point3_1n(ii), 2);
        l_conns1 = intersect(l_conns, nodes_degree_1);
        l_conns3 = intersect(l_conns, nodes_degree_3);
        if length(l_conns3) + length(l_conns1) == 2
            skeleton_connections(totol_n).p1 = l_conns1';
            skeleton_connections(totol_n).p3 = l_conns3';
            skeleton_connections(totol_n).end_points = [l_conns3; l_conns1]';
            skeleton_connections(totol_n).all_points = l_conns;
            totol_n = totol_n + 1;
        end
    end
end

for i = 1: n1
    l_conns = bfs(conn_map, nodes_degree_1(i), 2);
    if l_conns >= 2
        l_conns1 = intersect(l_conns, nodes_degree_1);
        l_conns3 = intersect(l_conns, nodes_degree_3);
        if length(l_conns3) + length(l_conns1) == 2
            skeleton_connections(totol_n).p1 = l_conns1';
            skeleton_connections(totol_n).p3 = l_conns3';
            skeleton_connections(totol_n).end_points = [l_conns3; l_conns1]';
            skeleton_connections(totol_n).all_points = l_conns;
            totol_n = totol_n + 1;
        end
    end
end

keep_jug = sort(cat(1, skeleton_connections.end_points), 2);
[~, keep_idx, ~] = unique(keep_jug, 'rows');
skeleton_connections = skeleton_connections(keep_idx);
num = length(skeleton_connections);
end_points = reshape(cat(1, skeleton_connections.end_points), [num * 2, 1]);
new_points = sort(unique(end_points));
index_map = zeros(points_n, 1);
index_map(new_points) = 1:length(new_points);

for i = 1: num
    ends = skeleton_connections(i).end_points;
    skeleton_connections(i).new_end_points = sort(index_map(ends))';
    skeleton_connections(i).distance = norm(points(ends(1), :) - points(ends(2), :));
    skeleton_connections(i).chose = 0;
end

% Save the modified coordinates of vertices as skeleton_coordinates_s1.mat
% Save the edges of simplified graph as skeleton_connections.mat
save('skeleton_coordinates_s1.mat', "points");
save('skeleton_connections.mat', "skeleton_connections");


% Visualization of the edges of simplified graph
end_points = cat(1, skeleton_connections.end_points);
starts = points(end_points(:, 1), :);
ends = points(end_points(:, 2), :);
line([starts(:,1), ends(:,1)]', [starts(:,2), ends(:,2)]', [starts(:,3), ends(:,3)]')
axis equal



function connections = find_nearest_n(points)
    tree = KDTreeSearcher(points);
    k = 6;
    [idx, dist] = knnsearch(tree, points, 'K', k);
    idx = idx(:, 2: end);
    dist = dist(:, 2: end);
    keep = dist < 10;
    connections = struct();
    for i = 1: size(idx, 1)
        nn = idx(i, :);
        connections(i).nn = nn(keep(i, :));
    end
end

function edges = construct_edge_list(point, point_edges)
    edges = zeros(length(point_edges), 2);
    edges(:, 1) = point;
    edges(:, 2) = point_edges';
end

function edges = construct_edge_from_connections(connections, point_list)
    edges = [];
    for i = 1: length(point_list)
        nn = connections(point_list(i)).nn;
        if ~isempty(nn)
            edge3 = construct_edge_list(point_list(i), nn);
            edges = [edges; edge3];
        end
    end
    edges = sort(edges, 2);
    edges = unique(edges, 'rows');
end

function conn_map = connect_matrix(edge_list, num)
    conn_n = length(edge_list);    
    conn_map = zeros([num, num]);
    for i = 1: conn_n
        conn_map(edge_list(i, 1), edge_list(i, 2)) = 1;
        conn_map(edge_list(i, 2), edge_list(i, 1)) = 1;
    end
end

function [new_edges, p3_edge, points, point3] = reconstruct_p3_connections(points, connections, point3, points_n)
    p3_edge = construct_edge_from_connections(connections, point3);
    conn_map = connect_matrix(p3_edge, points_n);

    new_point3 = zeros(length(point3), 1);
    new_edges = [];
    for i = 1: length(point3)
        conns = bfs(conn_map, point3(i), 10);
        p3 = sort(intersect(conns, point3));
        p1 = sort(setdiff(conns, p3));
        new_edges = [new_edges; construct_edge_list(p3(1), p1)];
        new_point3(i) = p3(1);
        if length(conns) >= 2
            points(p3(1), :) = mean(points(conns, :));
        end
    end
    point3 = unique(new_point3, 'stable');
    new_edges = sort(new_edges, 2);
    new_edges = unique(new_edges, 'rows');
end

function visited_nodes = bfs(A, start_node, threhold)
    N = size(A, 1); 
    visited = false(N, 1);
    queue = start_node;
    visited(start_node) = true;

    while ~isempty(queue)
        node = queue(1);
        queue(1) = [];
        neighbors = find(A(node, :) == 1);
        if length(neighbors) <= threhold   
            for n = neighbors
                if ~visited(n)
                    queue(end + 1) = n; 
                    visited(n) = true;
                end
            end
        end
    end
    visited_nodes = find(visited);
end

function d = fitLine3D_LS(points)
    p0 = mean(points, 1);
    centered = points - p0;
    
    covariance = centered' * centered;
    [~, ~, V] = svd(covariance);
    d = V(:, 1)';
    d = d / norm(d);
end
