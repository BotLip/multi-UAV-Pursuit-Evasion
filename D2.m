clc; clear;
%% INPUT 示例点集边界
P = [ 0.5    0;
      0      0.5;
     -0.5   -0.5;
     -0.2   -0.1;
     -0.1    0.1;
      0.1   -0.1;
      0.1    0.1 ];

% P = [1, 1; 2, 1; 3, 1; 2, 2];
% P = [0.500000000000000,0;0.0774601159010643,0.563135275354924;-0.445563486626429,-0.583870928091105;-0.101005050633883,-0.0858578643762691;-0.00136060761678562,0.0835601012694643;0.139971828289870,-0.191584100520549;0.178176130281010,0.162237088117667]
% P = [0.500000000000000,0;
%     0.405462267180715,0.434073411416459;
%     0.0224093910278814,-0.294904864720432;
%     0.443467170879759,-0.00807611844574879;
%     0.489206128922219,-0.0432829409898653;
%     0.597518786643932,-0.380845905278252;
%     0.674565231495186,0.234013816643841];
X_MAX = 3; X_MIN = -2;
Y_MAX = 3; Y_MIN = -2;
idx = 2 % evador 在P中的下标

sample_t = 0.01;
sample_f = 1 / sample_t;
T = 1;  % 计时

e_path = generateCircle([(X_MIN + X_MAX) / 2, (Y_MAX + Y_MIN) / 2], 2, sample_t);
e_step = 300;
P(idx, :) = e_path(e_step, :);

e_A = [];   % evador细胞面积

try
        voronoi(P(:, 1), P(:, 2));
    catch
        warning('Problem using delaunayTriangulation function. Points may lie on a line, Add a value of random');
        P_ = P + randi(size(P));
        voronoi(P_(:, 1), P_(:, 2));
end


%% run algorithm
while ~captured(P, idx)
    
    [miu_star, vncomb] = voronoiShrink(P, idx, X_MAX, X_MIN, Y_MAX, Y_MIN);
    
    %% show result
    hold on;
    axis([X_MIN, X_MAX, Y_MIN, Y_MAX]);
    plot_voronoi(vncomb, P, X_MAX, X_MIN, Y_MAX, Y_MIN);
    
    for i = 1:size(P, 1)
        % 在点坐标右上方偏移(0.1,0.1)处添加序号
        text(P(i,1)+0.05, P(i,2)+0.05, num2str(i), ...
            'FontSize', 8, ...
            'Color', 'red', ...
            'BackgroundColor', 'white', ...
            'HorizontalAlignment', 'left');
    end
    
    for i = 1 : size(vncomb, 2)
        if vncomb(i).neibor == 1
            x0 = vncomb(i).node(1); y0 = vncomb(i).node(2);
            nh = vncomb(i).nh .*0.1;
            nv = vncomb(i).nv.*0.1;
            u = vncomb(i).miustar.*0.1;
            quiver(x0, y0, nh(1), nh(2), 'LineWidth', 1, 'Color', 'r', 'MaxHeadSize', 0.1);
            quiver(x0, y0, nv(1), nv(2), 'LineWidth', 1, 'Color', 'r', 'MaxHeadSize', 0.1);
            quiver(x0, y0, u(1), u(2), 'LineWidth', 1, 'Color', 'b', 'MaxHeadSize', 0.1);
        end
    end
    hold off;
    
%     vertexs = [];
%     for i = 1 : size(vncomb, 2)
%         vertexs = [vertexs; vncomb(i).v1; vncomb(i).v2];
%     end
%     vertexs = unique(vertexs, 'rows');
%     e_A  = [e_A, polygonArea(vertexs)];
    
    drawnow;
    clf;
    % 执行一步控制
    for i = 1 : size(P, 1)
        if i ~= idx
            P(i, :) = P(i, :) + sample_t * miu_star{i};
        else
            e_step = e_step + 1;
            e_step = mod(e_step, size(e_path, 1));
            P(i, :) = e_path(e_step, :);
        end
    end
    T = T + 1;
end

% figure(2);
% plot(linspace(1, T), e_A);

%% Algorithm start
function [miu_star, vncomb] = voronoiShrink(P, idx, X_MAX, X_MIN, Y_MAX, Y_MIN)
    %P : 所有无人机的位置 
    %idx : evador在P的下标 XYMAXMIN是环境边界
    %输出最优控制u, 绘图需要vncomb 

    try
        dt = delaunayTriangulation(P);
    catch
        warning('Problem using delaunayTriangulation function. Points may lie on a line, Adding a value of random');
        P = P + randi(size(P));
        dt = delaunayTriangulation(P);
    end
    [V,r] = voronoiDiagram(dt);



    r{idx};
    % 控制顶点
    adjacentver = V(r{idx},:);
    % 获取所有三角形的邻接点
    adjacentIndices = [];
    for i = 1:size(dt.ConnectivityList, 1)
        % 检查每个三角形，看看第一个点是否在这个三角形的顶点中
        triangle = dt.ConnectivityList(i, :);

        if any(triangle == idx)  % 如果第一个点在三角形中
            % 找出三角形中与第一个点相邻的其他点
            adjacentIndices = [adjacentIndices, triangle(triangle ~= idx)];
        end
    end

    % 删除重复的下标
    adjacentIndices = unique(adjacentIndices)


    %% 对应点和边
    adjv_raw = adjacentver;
    adji = adjacentIndices;
    adjv = [];
    for i = 1 : size(adjv_raw, 1)
        b = any(ismember([inf,-inf], adjv_raw(i, :)));
        if ~b
            adjv = [adjv; adjv_raw(i, :)];
        end
    end

    %% 匹配和计算所有点和nv nh l L
    vncomb(1).neibor = 0;% 这个pursuer在范围边界内是否真相邻
    vncomb(1).v1 = [];  % cell边的第一个顶点
    vncomb(1).v2 = [];  % cell边的第二个顶点
    vncomb(1).i = 0;    % P中下标
    vncomb(1).k = 0;    % cell边斜率
    vncomb(1).node = [];% Pi
    vncomb(1).nh = [];  % nh
    vncomb(1).nv = [];  % nv
    vncomb(1).l = 0;    % l
    vncomb(1).L = 0;    % L
    vncomb(1).miustar = []; % 最优控制
    vncomb(1).epsilon = []; % e-p向量
    vncomb(1).epsilon_len = 0;  % 向量长度

    for i = 1 : size(adji, 2)
        vncomb(i).i = adji(i);
        vncomb(i).node = P(adji(i), :);
        Q1 = P(idx, :); Q2 = P(adji(i), :);
        vncomb(i).epsilon = Q1 - Q2;
        vncomb(i).epsilon_len = sqrt( sum(vncomb(i).epsilon.^2) );
        [ck, k] = findMidPerpendicularPoints(Q1, Q2, adjv);
        ck = adjv(ck, :);
        if (size(ck, 1) == 2)
            vncomb(i).v1 = ck(1, :);
            vncomb(i).v2 = ck(2, :);
            vncomb(i).neibor = 1;
        elseif (size(ck, 1) == 1)
            vncomb(i).v1 = ck(1, :);
            % 这里就应该给出至多两个无限长边与矩形范围相交的点了---------------------------------------
            % 如果这个v1已经在范围外了(正常情况，因为两细胞的交点可能因为斜率超出范围边界很远)就认为他俩不相邻
            % 不相邻就给neobor标false然后剩下的计算不用在做了
            % vncomb(i).v2 = [inf, inf];
            if (vncomb(i).v1(1) < X_MAX) && (vncomb(i).v1(1) > X_MIN) && (vncomb(i).v1(2) < Y_MAX) && (vncomb(i).v1(2) > Y_MIN)
                vncomb(i).neibor = 1;
                vncomb(i).v2 = findInfVertex(Q1, Q2, vncomb(i).v1, k, P(adji, :), X_MAX, X_MIN, Y_MAX, Y_MIN);
            end
            
        else
            error("more than two online vertex found, expected <=2 ");
        end
        
        if vncomb(i).neibor == 1
            vncomb(i).k = k;
            mid = (Q1 + Q2) / 2;
            [vncomb(i).nh, vncomb(i).nv, vncomb(i).l, vncomb(i).L] = ...
                findhvlL(mid, P(idx, :), P(adji(i), :), vncomb(i).epsilon, vncomb(i).v1, vncomb(i).v2, vncomb(i).k);
            alphah = - vncomb(i).L / 2;
            alphav = (vncomb(i).l.^2 - (vncomb(i).l - vncomb(i).L).^2) / vncomb(i).epsilon_len;
            vncomb(i).miustar = -(alphah ./ sqrt(alphah.^2 + alphav.^2) .* vncomb(i).nh +...
                alphav / sqrt(alphah.^2 + alphav.^2) .* vncomb(i).nv);
        end
        
    end

    % 邻接的pursuier就用miustar， 其他的简单朝向evador飞行
    miu_star = cell(size(P, 1), 1);
    for i = 1 : size(vncomb, 2)
        miu_star{vncomb(i).i} = vncomb(i).miustar;
    end
    for i = 1 : size(miu_star, 1)
        if isempty(miu_star{i})
            dir = P(idx, :) - P(i, :);
            dir = dir / sqrt(sum(dir.^2));
            miu_star{i} = dir;
        end
    end
    
end



%% uitils
function [result, k] = findMidPerpendicularPoints(Q1, Q2, points, epsilon)
    % 输入：Q1, Q2为线段端点，points为待检测点集，epsilon为容差（默认1e-6）
    % 输出：逻辑索引，标记符合条件的点 
    
    if nargin < 4 
        epsilon = 1e-6;
    end 
    
    % 计算中垂线参数 
    dx = Q2(1) - Q1(1);
    dy = Q2(2) - Q1(2);
    k = - 1 / (dy / dx);
    M = (Q1 + Q2) / 2;
    
    % 避免线段退化为点 
    if dx == 0 && dy == 0 
        error('线段端点重合，无法定义中垂线');
    end 
    
    % 计算方程误差 
    lhs = dx * points(:,1) + dy * points(:,2);
    rhs = dx * M(1) + dy * M(2);
    errors = abs(lhs - rhs);
    
    % 返回符合条件的逻辑索引 
    result = errors < epsilon;
end 

function dis = lineardis(Q1, Q2)
    if size(Q1, 1) ~= 1 || size(Q2, 1) ~= 1
        error(" error vertex input, wrong number of nodes recieved expected 2");
    end
    dis = sqrt( sum((Q1 - Q2) .^ 2) );
end

function isOnSegment = isPointOnSegment(Q1, Q2, P, tol)
    % 输入：Q1, Q2为线段端点，P为待测点，tol为容差（默认1e-10）
    % 输出：逻辑值，true表示P在线段Q1Q2上 
    
    if nargin < 4 
        tol = 1e-10;
    end 
    
    % 检查线段是否退化为点 
    if norm(Q2 - Q1) < tol 
        isOnSegment = (norm(P - Q1) < tol);
        return;
    end 
    
    % 共线性检验（叉积）
    crossProduct = (P(1) - Q1(1)) * (Q2(2) - Q1(2)) - (P(2) - Q1(2)) * (Q2(1) - Q1(1));
    isColinear = abs(crossProduct) < tol;
    
    % 坐标范围检验 
    withinX = (P(1) > min(Q1(1), Q2(1))) && (P(1) < max(Q1(1), Q2(1)));
    withinY = (P(2) > min(Q1(2), Q2(2))) && (P(2) < max(Q1(2), Q2(2)));
    
    isOnSegment = isColinear && withinX && withinY;
end 

% deprecated function
function position = checkPointPosition(points, testPoint)
    % points: 一个 n x 2 的矩阵，其中每一行是一个点的 (x, y) 坐标
    % testPoint: 一个 1 x 2 向量，表示待检测的点 (x, y)
    % 返回值 position: 一个字符串，表示 testPoint 的位置（'left', 'right', 'top', 'bottom', 'other'）

    % 获取点集的最小和最大 x, y 坐标
    minX = min(points(:,1));
    maxX = max(points(:,1));
    minY = min(points(:,2));
    maxY = max(points(:,2));

    % 判断 testPoint 是否在最左端、最右端、最上端或最下端
    if testPoint(1) == minX && testPoint(2) == minY
        position = [1, 0, 0, 1];   % 左下
    elseif testPoint(1) == minX && testPoint(2) == maxY
        position = [1, 0, 1, 0];      % 左上
    elseif testPoint(1) == maxX && testPoint(2) == minY
        position = [0, 1, 0, 1];  % 右下
    elseif testPoint(1) == maxX && testPoint(2) == maxY
        position = [0, 1, 1, 0];     % 右上
    else
        position = 'other';  % 如果不在四个端点上
        error(" error position ");
    end
end



function [nh, nv, l, L] = findhvlL(mid, e, p, epsilon, v1, v2, k)
    % mid: evador 和 pursuieri的中点， adjv:evador所有的细胞顶点，epsilon:evador 和 Pi的向量, i:P_i
    % v1, v2:Pi对应的细胞边的顶点， k:边的斜率
    l1 = lineardis(mid, v1);
    l2 = lineardis(mid, v2);
    L = lineardis(v1, v2);
    epsilon_len = sqrt( sum(epsilon(1)^2 + epsilon(2)^2) );
    nh = epsilon / epsilon_len;
    if isPointOnSegment(v1, v2, mid, 1e-3)
        nv = getOrthogonalUnitVector(nh);
%         if all(mid == v1) || all(mid == v2)
%             l
%         else
        if dot(nv, (v1 - mid)) > 0
            l = l2;
        elseif dot(nv, (v2 - mid)) > 0
            l = l1;
        end
%         end
    else    % mid不在L上，nv必须指向(v1 or v2 - mid)的相反方向 使得l存在 且l > L
        nv = getOrthogonalUnitVector(nh);   
        if dot(nv, (v1 - mid)) > 0 || dot(nv, (v2 - mid)) > 0
            nv = -nv;
        end
        l = max(l1, l2);
    end
end

function v2 = findInfVertex(e, p, v1, k, Pe, X_MAX, X_MIN, Y_MAX, Y_MIN)
    % 暂时
    % 出现开放细胞的情况是 adjv半包围 evador
    %  然后分情况讨论开放的边，最多就两条
    % 已知当前的点是开放细胞的一个开放点Pi，随机取adji里的任意一个不是这个开放点的点一定在evador和Pi直线的一侧，
    % (不然就相交而不是开放点了) 从 curvertex取另一侧作为射线找到与边界的交点作为Pi的v2
    % e是evador p是pursuer， v1是射线起点， k是射线斜率， Pe是所有的邻接细胞的节点
    dx = p(1) - e(1);
    dy = p(2) - e(2);
    k_ = dy / dx;
    b = p(2) - k_ * p(1);
    if isinf(k_)    % 防止inf出现NAN
        k_ = 1e9;
    elseif k_ == 0
        k_ = 1e-3;
    end
%     lhs = dx * p(1) + dy * p(2);
    other_node = [];
    for i = 1 : size(Pe, 1)
        if Pe(i, 1) ~= p(1) || Pe(i, 2) ~= p(2)
            other_node = Pe(i, :);
            break;
        end
    end
    if isempty(other_node)
        error("no other node found");
    end
%     rhs = dx * other_node(1) + dy * other_node(2);
    lhs = b + k_ * other_node(1);
    rhs = other_node(2);
    df = rhs - lhs;
    dirx = 0; diry = 0;
    if df > 0 && k_ > 0
        dirx = 1; diry = -1;
    elseif (df < 0 && k_ < 0)
        dirx = 1; diry = 1;
    elseif (df < 0 && k_ > 0)
        dirx = -1; diry = 1;
    elseif (df > 0 && k_ < 0)
        dirx = -1; diry = -1;
    end
    
    v2 = findRayIntersection(v1(1), v1(2), k, X_MIN, X_MAX, Y_MIN, Y_MAX, dirx, diry);
    
end


function intersectionPoints = findRayIntersection(x0, y0, slope, XMIN, XMAX, YMIN, YMAX, dirX, dirY)
    % x0, y0: 射线起点
    % slope: 斜率
    % XMIN, XMAX, YMIN, YMAX: 矩形的边界
    % direction: 射线的x轴方向，±1 表示
    % 返回值 intersectionPoints: 包含射线与矩形的交点

    % 射线方程： y = m * (x - x0) + y0 或 y = m * x + b
    b = y0 - slope * x0;  % 截距
    
    intersectionPoints = []; % 初始化交点
    
    if isinf(slope)
        if dirY == -1
            intersectionPoints = [intersectionPoints; x0, YMIN]; 
        elseif dirY == 1
            intersectionPoints = [intersectionPoints; x0, YMAX]; 
        else 
            error("error empty intersectionPoints");
        end
        return ;
    elseif slope == 0
        if dirX == -1
            intersectionPoints = [intersectionPoints; XMIN, y0];
        elseif dirX == 1
            intersectionPoints = [intersectionPoints; XMAX, y0;];
        else 
            error("error, empty intersectionPoints");
        end
        return ;
    end
    
    % 1. 计算射线与矩形的四个边的交点
    % 计算与左边界 X = XMIN 的交点
    yLeft = slope * XMIN + b;
    if yLeft >= YMIN && yLeft <= YMAX
        intersectionPoints = [intersectionPoints; XMIN, yLeft];
    end

    % 计算与右边界 X = XMAX 的交点
    yRight = slope * XMAX + b;
    if yRight >= YMIN && yRight <= YMAX
        intersectionPoints = [intersectionPoints; XMAX, yRight];
    end

    % 计算与下边界 Y = YMIN 的交点
    xBottom = (YMIN - b) / slope;
    if xBottom >= XMIN && xBottom <= XMAX
        intersectionPoints = [intersectionPoints; xBottom, YMIN];
    end

    % 计算与上边界 Y = YMAX 的交点
    xTop = (YMAX - b) / slope;
    if xTop >= XMIN && xTop <= XMAX
        intersectionPoints = [intersectionPoints; xTop, YMAX];
    end
   
    
    
    
    % 2. 根据射线的方向，筛选合适的交点
    if dirX == 1  % 射线沿正 x 方向
        % 筛选 x > x0 的交点
        intersectionPoints = intersectionPoints(intersectionPoints(:, 1) > x0, :);
    elseif dirX == -1  % 射线沿负 x 方向
        % 筛选 x < x0 的交点
        intersectionPoints = intersectionPoints(intersectionPoints(:, 1) < x0, :);
    end
    
    if dirY == 1
        intersectionPoints = intersectionPoints(intersectionPoints(:, 2) > y0, :);
    elseif dirY == -1
        intersectionPoints = intersectionPoints(intersectionPoints(:, 2) < y0, :);
    end
    
%     % 如果有多个交点，返回最小的一个交点 
%     if ~isempty(intersectionPoints)
%         % 按 x 坐标排序，取最小的那个交点
%         intersectionPoints = min(intersectionPoints, [], 1);
%     end
    intersectionPoints = unique(intersectionPoints, 'rows');  % 处理正好交在边界范围的顶角上
end


function orthogonalVector = getOrthogonalUnitVector(v)
    % v: 输入单位向量 (vx, vy)
    % 返回值 orthogonalVector: 与 v 正交的单位向量
    
    % 假设 v 已经是单位向量
    % 正交单位向量为 (-vy, vx)
    orthogonalVector = [-v(2), v(1)];
end

function plot_voronoi(vncomb, P, XMAX, XMIN, YMAX, YMIN)
    line([XMAX, XMAX], [YMAX, YMIN], 'Color', 'green', 'LineWidth', 1);
    line([XMIN, XMIN], [YMAX, YMIN], 'Color', 'green', 'LineWidth', 1);
    line([XMAX, XMIN], [YMAX, YMAX], 'Color', 'green', 'LineWidth', 1);
    line([XMAX, XMIN], [YMIN, YMIN], 'Color', 'green', 'LineWidth', 1);
    voronoi(P(:, 1), P(:, 2));
    scatter(P(:, 1), P(:, 2));
    for i = 1 : size(vncomb, 2)
        if vncomb(i).neibor == 1
            scatter(vncomb(i).node(1), vncomb(i).node(2));
            line([vncomb(i).v1(1), vncomb(i).v2(1)],[vncomb(i).v1(2), vncomb(i).v2(2)],'Color','red');
        end
    end
end

function over = captured(P, idx)
    over = false;
    for i = 1 : size(P, 1)
        if lineardis(P(i, :), P(idx, :)) < 0.1 && i ~= idx
            over = true;
        end
    end
end

function points = generateCircle(center, radius, sample_t)
    % center: 圆心坐标 [x, y]
    % radius: 圆的半径
    % sample_t: 离散点之间的最大距离
    % 返回值 points: 圆上离散点的坐标

    % 圆心坐标
    x0 = center(1);
    y0 = center(2);
    
    % 计算圆的周长
    circumference = 2 * pi * radius;
    
    % 计算点的数量，保证每个点之间的距离不超过 sample_t
    num_points = ceil(circumference / sample_t);
    
    % 角度增量
    theta_increment = 2 * pi / num_points;
    
    % 初始化存储圆上点的矩阵
    points = zeros(num_points, 2);
    
    % 生成圆上的离散点
    for i = 1:num_points
        theta = (i - 1) * theta_increment;  % 当前点的角度
        x = x0 + radius * cos(theta);       % 计算 x 坐标
        y = y0 + radius * sin(theta);       % 计算 y 坐标
        points(i, :) = [x, y];
    end
end

% deprecated function
function area = polygonArea(vertices)
    % vertices: 输入一个 N x 2 矩阵，表示多边形的顶点坐标
    
    % 将顶点按逆时针顺序排列（使用凸包算法）
    % 计算凸包的顶点顺序
    k = convhull(vertices(:, 1), vertices(:, 2));  % k 是凸包顶点的顺序
    sortedVertices = vertices(k, :);  % 排序后的顶点
    
    % 计算凸多边形的面积（应用格林公式）
    n = size(sortedVertices, 1);  % 顶点数
    area = 0;
    
    for i = 1:n
        j = mod(i, n) + 1;  % 循环连接最后一个点到第一个点
        area = area + (sortedVertices(i, 1) * sortedVertices(j, 2) - sortedVertices(j, 1) * sortedVertices(i, 2));
    end
    
    area = abs(area) / 2;  % 取绝对值并除以 2
end





