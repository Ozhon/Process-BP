function [bestLine, inliers] = fitLineRANSAC(points, threshold)
    % 参数设置
    maxIterations = 2000;
    bestScore = 0;
    bestLine = [];
    inliers = [];
    
    if size(points,1) < 2
        return;
    end
    
    for iter = 1:maxIterations
        % 随机选择两个点
        sampleIds = randperm(size(points,1), 2);
        p1 = points(sampleIds(1),:);
        p2 = points(sampleIds(2),:);
        
        % 计算直线方向
        dir = p2 - p1;
        if norm(dir) < eps
            continue; % 跳过重合点
        end
        dir = dir / norm(dir);
        
        % 计算所有点到直线的距离
        distances = point_to_line_distance(points, p1, dir);
        
        % 统计内点
        currentInliers = find(distances < threshold);
        currentScore = length(currentInliers);
        
        % 更新最佳模型
        if currentScore > bestScore
            bestScore = currentScore;
            inliers = currentInliers;
            bestLine = [p1, dir]; % 存储点和方向
        end
    end
end

% 计算点到直线的距离
function distances = point_to_line_distance(points, linePoint, lineDir)
    % 向量AP
    AP = points - linePoint;
    
    % 计算垂直距离：|AP × dir| / |dir| (因为dir是单位向量，分母为1)
    crossProd = cross(AP, repmat(lineDir, size(points,1), 1));
    distances = sqrt(sum(crossProd.^2, 2));
end