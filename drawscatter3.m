function drawscatter3(data1, s, c)
    x = data1(:, 1);
    y = data1(:, 2);
    z = data1(:, 3);    
    scatter3(x, y, z, s, 'filled', 'MarkerFaceColor', c);
    hold on
    axis equal
end