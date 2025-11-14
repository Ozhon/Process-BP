% Polish the final data

% Read "all_lines_last.mat" from a6_lines_get_points.m
clear;
ColorMatrix = [0.941, 0.749, 0.196;  % 黄色
       0.741, 0.294, 0.380;  % 红色
       0.384, 0.443, 0.580;  % 蓝色
       0.113, 0.474, 0.501];  % 绿色

all_lines = load("all_lines_last.mat").all_lines;
vectors = cat(1, all_lines.vector);

% Group the lines by DBSCSAN
[theta, phi] = turn_polar(vectors);
data = [theta, phi];
idx = dbscan(data, 0.3, 5);
gscatter(theta, phi, idx);

%% Integrate short lines
% Update the lines and idx
keep_a = idx ~= -1;
all_lines = all_lines(keep_a);
idx = idx(keep_a);

[all_lines, idx] = integrate_lines(all_lines, idx, 90);
all_lines = sort_points(all_lines, 30);
[all_lines.idx] = deal(idx);
%% Visualize the disclination lines
figure;
for i = 1: 4
    outlines = all_lines(idx == i);
    for ii = 1: length(outlines)
        l_points = outlines(ii).points;
        drawscatter3(l_points, 12, ColorMatrix(i, :));
    end
    axis equal;
end
%% Save the final disclination lines as "all_lines_last.mat"
save("all_lines_last.mat", "all_lines");


function [theta, phi] = turn_polar(points)
    [az, el, ~] = cart2sph(points(:, 1), points(:, 2), points(:, 3));
    theta = az + pi;
    phi = pi/2 - el;
end

function distance_matrix = construct_distance_matrix(lines1)
    n = length(lines1);
    distance_matrix = zeros([n, n]);
    for i = 1: n
        vector1 = lines1(i).vector;
        center1 = lines1(i).mean;
        for j = 1: n
            points2 = lines1(j).points;
            V = points2 - center1;
            N = cross(V, repmat(vector1, size(V, 1), 1), 2);
            N_norm = sqrt(sum(N.^2, 2));
            D = N_norm / norm(vector1);
            d = mean(D);
            if i == j
                distance_matrix(i, j) = 1000;
            else
                distance_matrix(i, j) = d;
            end
        end
    end
end

function [all_lines, idx] = integrate_lines(all_lines, idx, thr)
    for i = 1: 4
        lines_idx = find(idx == i);
        lines1 = all_lines(lines_idx);
        distance_matrix = construct_distance_matrix(lines1);
        conn_map = distance_matrix < thr;
        if any(any(conn_map))
            disp_idx = [];
            for ii = 1: length(conn_map)
                if ~ismember(ii, disp_idx)
                    inline_idx = find(conn_map(ii, :));
                    if ~isempty(inline_idx)
                        main_line = lines_idx(ii);
                        tree_lines = lines_idx(inline_idx);
                        main_points = all_lines(main_line).points;
                        tree_points = cat(1, all_lines(tree_lines).points);
                        full_points = [main_points; tree_points];
                        [bestline, inliners] = fitLineRANSAC(full_points, 20);
                        std_v = [1; 1; 1];
                        all_lines(main_line).points = full_points;
                        all_lines(main_line).mean = mean(full_points, 1);
                        all_lines(main_line).p0 = bestline(1: 3);
                        if dot(std_v, bestline(4: 6)) > 0
                            all_lines(main_line).vector = bestline(4: 6);
                        else
                            all_lines(main_line).vector = -bestline(4: 6);
                        end
                        all_lines(main_line).inliers = inliners;
                        disp_idx = [disp_idx, inline_idx];
                    end
                end
            end
            all_lines(lines_idx(disp_idx)) = [];
            idx(lines_idx(disp_idx)) = [];
        end
    end
end

function all_lines = sort_points(all_lines, ths)
    for i = 1: length(all_lines)
        if (size(all_lines(i).points, 1) - length(all_lines(i).inliers)) < ths
            line_points = all_lines(i).points(all_lines(i).inliers, :);
        else
            line_points = all_lines(i).points;
        end
        line_vector = all_lines(i).vector;
        line_mean = all_lines(i).mean;

        V = line_points - line_mean;
        t = V * line_vector';
        [~, p_idx] = sort(t);
        all_lines(i).points = line_points(p_idx, :);
    end
end
