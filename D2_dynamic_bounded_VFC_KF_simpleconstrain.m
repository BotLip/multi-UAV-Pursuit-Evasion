clc; clear; % 不做严格的约束，再计算最优，只要求计算完最优之后限制一下方向即可，所有实际执行的不是最优的控制
%% INPUT 示例点集边界
P = [ 0.5    0;
      0      0.5;
     -0.5   -0.5;
     -0.2   -0.1;
     -0.1    0.1;
      0.1   -0.1;
      0.1    0.1 ];
% P = [ 3-0.2    -1.8;
%       0      0.5;
%      3-0.1   -1.8;
%      3-0.1   -1.6;
%      3-0.1    -1.7;
%       3-0.2   -1.7;
%       3-0.2    -1.6 ];
P = [1, 1; 2, 1; 3, 1; 2, 2];
% P = [0.500000000000000,0;0.0774601159010643,0.563135275354924;-0.445563486626429,-0.583870928091105;-0.101005050633883,-0.0858578643762691;-0.00136060761678562,0.0835601012694643;0.139971828289870,-0.191584100520549;0.178176130281010,0.162237088117667]
% P = [0.500000000000000,0;
%     0.405462267180715,0.434073411416459;
%     0.0224093910278814,-0.294904864720432;
%     0.443467170879759,-0.00807611844574879;
%     0.489206128922219,-0.0432829409898653;
%     0.597518786643932,-0.380845905278252;
%     0.674565231495186,0.234013816643841];
X_MAX = 10; X_MIN = -10;
Y_MAX = 10; Y_MIN = -10;

fov_rel = [-5, 5, -5, 5];

idx = 2 % evador 在P中的下标

sample_t = 0.01;
sample_f = 1 / sample_t;
T = 1;  % 计时

e_slow_scale = 0.5; % 假设evador的速度较慢

% e_path = generateCircle([(X_MIN + X_MAX) / 2, (Y_MAX + Y_MIN) / 2], 2, sample_t); % 让e走逆时针的圆
% e_path = ones(1000, 2); % 让e停在1，1这个位置
a = linspace(1, -100, 100*(sample_f+100));b = linspace(1, 1, 100*(sample_f+100));e_path = [a; b]; e_path = e_path'; % 让e以较慢的速度从(1，1)向左
e_step = 300;
P(idx, :) = e_path(e_step, :);

e_A = [];   % evador细胞面积

% 卡尔曼滤波器需要的历史观测
history.velocity = [0, 0; 0, 0];
history.position = [P(idx, :); P(idx, :);];
history.n = 10;

his_vel_plot = [0, 0; 0, 0; ];
his_e_pos_plot = [P(idx, :); P(idx, :);];
his_all_pos_plot = zeros(size(P, 1), 1, size(P, 2)); his_all_pos_plot(:, end, :) = P;
his_vel_pred_plot = [];
his_pos_pred_plot = [];

v_pre = NaN(size(P, 1), 2);
wmax = 100;

try
        voronoi(P(:, 1), P(:, 2));
catch
        warning('Problem using delaunayTriangulation function. Points may lie on a line, Add a value of random');
        P_ = P + randi(size(P));
        voronoi(P_(:, 1), P_(:, 2));
end


%% run algorithm
while ~captured(P, idx)
    
    e = P(idx, :);
%     e_vel = (e_path(T+e_step, :) - e_path(T-1+e_step, :)) * sample_f;
    
    fov = [fov_rel(1:2) + e(1), fov_rel(3:4) + e(2)];
    x_min_rel = fov_rel(1) + e(1); x_max_rel = fov_rel(2) + e(1); y_min_rel = fov_rel(3) + e(2); y_max_rel = fov_rel(4) + e(2);
    % 给evador分配控制e_vel 给pursuer分配miustar
    e_vel = virtualForceControl(P(idx, :), P([1:idx-1, idx+1:end], :), x_max_rel, x_min_rel, y_max_rel, y_min_rel);
    e_vel = e_vel * e_slow_scale;
    e_vel_pred = e_vel;
    % 可以假设直接知道e_vel，也可以用别的方法预测当前时刻evador的速度
    xk = KalmanFilter(sample_t, history.velocity, history.position);
    e_vel_pred = xk(3:4);
    [miu_star, vncomb] = voronoiShrink(P, idx, x_max_rel, x_min_rel, y_max_rel, y_min_rel, e_vel_pred, v_pre, sample_t, wmax, T);
    
    % show result
    hold on;
    
    axis(fov);
    plot_voronoi(vncomb, P, x_max_rel, x_min_rel, y_max_rel, y_min_rel);
    
    for i = 1:size(P, 1)
        % 在点坐标右上方偏移(0.1,0.1)处添加序号
        text(P(i,1)+0.05, P(i,2)+0.05, num2str(i), ...
            'FontSize', 8, ...
            'Color', 'red', ...
            'BackgroundColor', 'white', ...
            'HorizontalAlignment', 'left');
    end
    
    for i = 1 : size(P, 1)
        if i == idx
            quiver(P(idx, 1), P(idx, 2), e_vel(1)*0.2, e_vel(2)*0.2,'LineWidth', 1, 'Color', 'r', 'MaxHeadSize', 0.1);
        else
            x0 = P(i, 1); y0 = P(i, 2);
            u = miu_star{i}.*0.4;
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
            v_pre(i, :) = miu_star{i};
        else
            e_step = e_step + 1;
            e_step = mod(e_step, size(e_path, 1));
            P(i, :) = P(i, :) + e_vel * sample_t;
            v_pre(i, :) = e_vel;
        end
    end
    history.velocity(end + 1, :) = e_vel;   his_vel_plot(end+1, :) = e_vel;         his_pos_pred_plot(end+1, :) = xk(1:2);
    history.position(end + 1, :) = P(idx, :); his_e_pos_plot(end+1, :) = P(idx, :);   his_vel_pred_plot(end+1, :) = xk(3:4);
    refreshHistory(history);
    his_all_pos_plot(:, end+1, :) = P;
    
    T = T + 1;
end

%% show kalman filter
figure(3);
hold on;
x = 1 : T;

scatter(his_e_pos_plot(3:end, 1), his_e_pos_plot(3:end, 2));
scatter(his_pos_pred_plot(:, 1), his_pos_pred_plot(:, 2));
% his_vel_plot = his_pos_plot * sample_t; his_pos_pred_plot
% quiver(his_pos_pred_plot(3:end, 1))
legend("actual pos", "predict pos");
hold off;

figure(4);
diff = his_vel_pred_plot -  his_vel_plot(3:end, :);
diff_len = zeros(size(diff, 1), 1);
for i = 1 : size(diff, 1)
    diff_len(i) = norm(diff(i, :));
end

plot(linspace(0, T*sample_t, size(diff, 1)), diff_len);

figure(5);
hold on;
ph1 = []; ph2 = [];
for i = 1 : size(P, 1)
    if i == idx
        color = 'red';
        ph1 = plot3(his_all_pos_plot(i, :, 1), his_all_pos_plot(i, :, 2), linspace(0, T*sample_t, size(his_all_pos_plot, 2)), 'Color', color);
%         legend(ph, 'evader')
    else
        color = 'blue';
        ph2 = plot3(his_all_pos_plot(i, :, 1), his_all_pos_plot(i, :, 2), linspace(0, T*sample_t, size(his_all_pos_plot, 2)), 'Color', color);
%         legend(ph, 'pursuer')
    end
    xlabel("X");ylabel('Y');zlabel('t');
end
legend([ph1, ph2], ["evader", "pursuer"]);
hold off;



%% Algorithm start
function [miu_star, vncomb] = voronoiShrink(P, idx, X_MAX, X_MIN, Y_MAX, Y_MIN, e_vel, v_pre, sample_t, wmax, T)
    % P : 所有无人机的位置 
    % idx : evador在P的下标 XYMAXMIN是环境边界
    % 输出最优控制u, 绘图需要vncomb 

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
            alphav = (vncomb(i).l.^2 - (vncomb(i).l - vncomb(i).L).^2) / vncomb(i).epsilon_len / 2;
            vncomb(i).miustar = -(alphah ./ sqrt(alphah.^2 + alphav.^2) .* vncomb(i).nh +...
                alphav / sqrt(alphah.^2 + alphav.^2) .* vncomb(i).nv);
        end
        
    end

    % 邻接的pursuier就用miustar， 其他的简单朝向evador飞行
    miu_star = cell(size(P, 1), 1);
    for i = 1 : size(vncomb, 2)
        if vncomb(i).neibor
            ustar_p = vncomb(i).miustar;
            
            % ---------无约束的--------
            scale = -dot(ustar_p, e_vel) + sqrt(dot(ustar_p, e_vel)^2 + (1 - norm(e_vel)^2));
            miu = ustar_p * scale + e_vel;
            if T ~= 1
                % ---------有约束的--------
                miu = computeUltraUstar(ustar_p, v_pre(i, :), sample_t, wmax, 0.75);   % 无约束时的最优解不一定满足约束
            end
            miu_star{vncomb(i).i} = miu;
        end
    end
    for i = 1 : size(miu_star, 1)
        if idx ~= i && isempty(miu_star{i}) % 没有邻接的就直飞 / 根据1v1的策略飞
%             dir = P(idx, :) - P(i, :);
%             dir = dir / sqrt(sum(dir.^2));
%             miu_star{i} = dir;
            [miu_star{i}, ~, ~] = voronoiShrink_1v1([P(idx, :); P(i, :)], 1, X_MAX, X_MIN, Y_MAX, Y_MIN, e_vel, v_pre(i, :), sample_t, wmax, T);
        end
    end
    
end

%% 1v1 edition voronoi_shrink function
function [ustar_p, ustar_e, vertex] = voronoiShrink_1v1(P, idx, X_MAX, X_MIN, Y_MAX, Y_MIN, e_vel, v_pre, sample_t, wmax, T)
    % 输入 P 只包含逃避和追逐无人机， idx 为逃避无人机在P的下标， xymm是地图边界
    Q1 = P(idx, :);
    Q2 = P(3-idx, :);
    M = (Q1 + Q2) / 2;
    
    dx = Q1(1) - Q2(1); dy = Q1(2) - Q2(2);
    k_ = dy / dx;
    k = - 1 / k_;
    k_boo = k;
    if isinf(k_boo)
        k_boo = 1e9;
    elseif k_boo == 0
        k_boo = 1e-4;
    end

    if k_boo > 0
        dirx1 = 1; diry1 = 1;
        dirx2 = -1; diry2 = -1;
    else
        dirx1 = 1; diry1 = -1;
        dirx2 = -1; diry2 = 1;
    end

    x0 = M(1); y0 = M(2);

    vertex =[
        findRayIntersection(x0, y0, k, X_MIN, X_MAX, Y_MIN, Y_MAX, dirx1, diry1);
        findRayIntersection(x0, y0, k, X_MIN, X_MAX, Y_MIN, Y_MAX, dirx2, diry2)];

    l = 0;
    l1 = lineardis(M, vertex(1, :));
    l2 = lineardis(M, vertex(2, :));
    L = lineardis(vertex(1, :), vertex(2, :));
    cosi = Q1 - Q2;
    cosi_len = sqrt(sum(cosi.^2));
    nh = cosi / cosi_len;
    nv = getOrthogonalUnitVector(nh);

    if dot((vertex(1, :) - M), nv) > 0
        l = l2;
    else
        l = l1;
    end

    ah = - L / 2;
    av = (l^2 - (L - l)^2) / cosi_len / 2;
    ahv = sqrt(ah^2 + av^2);

    ustar_p = -( ah / ahv .* nh + av / ahv .* nv );
        % -----------无约束的-----------
        scale = -dot(ustar_p, e_vel) + sqrt(dot(ustar_p, e_vel)^2 + (1 - norm(e_vel)^2));
        ustar_p = scale * ustar_p + e_vel;
    if T ~= 1
        % -----------带约束的-----------
        ustar_p = computeUltraUstar(ustar_p, v_pre, sample_t, wmax, 0.75);
    end
    ustar_e = ( ah / ahv .* nh - av / ahv .* nv );
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

function norm = linearnorm(v)
    norm = v / sqrt(sum(v.^2));
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


function u = virtualForceControl(currentPos, Points, X_MAX, X_MIN, Y_MAX, Y_MIN)
    % 参数初始化 
    k_repel = 2.0;    % 无人机间斥力增益 
    k_boundary = 2; % 边界斥力增益 
    d_safe = 10;     % 安全距离（米）
    epsilon = 1e-5;   % 防止除零 
    
    % 初始化总控制力 
    u = [0, 0];
    
    % 计算其他无人机斥力 
    if ~isempty(Points)
        diff = currentPos - Points;                   % 坐标差向量 
        dist = sqrt(sum(diff.^2, 2)) + epsilon;       % 欧氏距离（避免零除）
        valid = dist < d_safe;                        % 筛选有效作用距离 
        if any(valid)
            dir_vec = diff(valid, :) ./ dist(valid);  % 单位方向向量 
            force = k_repel * (1./dist(valid) - 1/d_safe);
            u = u + sum(force .* dir_vec, 1);         % 累加斥力 
        end 
    end 
    
    % 计算边界斥力 
    th = 0.5;
    x = currentPos(1); y = currentPos(2);
    X_MID = (X_MAX + X_MIN) / 2; Y_MID = (Y_MAX + Y_MIN) / 2;
%     if abs(x - X_MAX) < th || abs(x - X_MIN) < th || abs(y-Y_MAX) < th || abs(y - Y_MIN) < th
    boundary_forces = zeros(4, 2);
    % 左右边界（X方向）
    boundary_forces(1, :) = [-k_boundary/exp(X_MID - x), 0];  % 右边界向左推
    boundary_forces(2, :) = [k_boundary/exp(x - X_MID), 0];  % 左边界向右推
    % 上下边界（Y方向）
    boundary_forces(3, :) = [0, -k_boundary/exp(Y_MID - y)];  % 上边界向下推
    boundary_forces(4, :) = [0, k_boundary/exp(y - Y_MID)];  % 下边界向上推
    u = u + sum(boundary_forces, 1);
%     end
    % 控制力归一化（可选）
    unorm = norm(u);
    if unorm > 1
        u = u / norm(u);
    end
end 


function xk = KalmanFilter(dt, velocity, position)
    % 位置数据：position = [x1, y1; x2, y2; ...]
    % 速度数据：velocity = [vx1, vy1; vx2, vy2; ...]
    % 加速度数据：acceleration = [ax1, ay1; ax2, ay2; ...]

    % 假设时间步长
    dt = 0.1;  % 每个采样时间间隔

    % 假设我们有 n 个历史数据点
    n = size(velocity, 1);

    % 初始化卡尔曼滤波器参数
    % 状态向量 [px, py, vx, vy]
    x = [position(1,1); position(1,2); velocity(1,1); velocity(1,2)];  % 初始状态
    P = eye(4) * 1000;  % 初始协方差矩阵，较大的初值表示初始的不确定性

    % 状态转移矩阵 F (考虑位置和速度变化)
    F = [1, 0, dt, 0;   % px = px + vx * dt
        0, 1, 0, dt;   % py = py + vy * dt
        0, 0, 1, 0;    % vx (速度保持不变)
        0, 0, 0, 1];   % vy (速度保持不变)

    % deprecated Buk_1
    % 控制矩阵 B (假设加速度或控制输入影响速度变化) 
%     B = [0.5 * dt^2, 0;    % 加速度对位置的影响
%         0, 0.5 * dt^2;    % 加速度对位置的影响
%         dt, 0;             % 加速度对速度的影响
%         0, dt];            % 加速度对速度的影响

    % 观测矩阵 H (我们可以直接观测速度)
    H = [1, 0, 0, 0;   % vx
        0, 1, 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1];  % vy

    % 过程噪声矩阵 Q
    Q = eye(4) * 0.01;

    % 观测噪声矩阵 R
    R = eye(4) * 0.01;

    % 用于存储卡尔曼滤波器的输出
    x_pred = zeros(n, 4);  % 存储每个时刻的预测速度

    % 卡尔曼滤波迭代
    for k = 2:n
        % 预测步骤
        x = F * x;  % 预测下一时刻的状态
        P = F * P * F' + Q;  % 更新协方差矩阵

        % 观测更新（假设我们通过传感器得到速度）
        z = [position(k, :), velocity(k, :)]';  % 当前时刻的实际速度

        % 卡尔曼增益
        K = P * H' / (H * P * H' + R);

        % 更新状态估计
        x = x + K * (z - H * x);

        % 更新协方差矩阵
        P = (eye(4) - K * H) * P;

        % 存储预测的速度
        x_pred(k, :) = x';  % 提取预测的速度分量
    end
    xk = x_pred(end, :);
end

function refreshHistory(his)
    if size(his.velocity,1) ~= size(his.position, 1)
        error("history error, size of velocities is not equal to positions");
    end
    if size(his.velocity, 1) > his.n
        his.velocity = his.velocity(2:end);
        his.position = his.position(2:end);
    end
end


%% 直接数学计算求ustar
% deprecated this shit always output NaN
% function [alpha, ultra_ustar] = computeAlpha(e_vel, ustar, v_pre, sample_t, wmax)
%     % 参数检查 
%     if norm(e_vel) < 0 || norm(e_vel) > 1 
%         error('e_vel必须满足0 < norm(e_vel) < 1');
%     end 
%     if isempty(ustar) || all(ustar == 0)
%         error('ustar不能为零向量');
%     end 
%     if sample_t <= 0 
%         error('sample_t必须为正数');
%     end 
%  
%     % 计算关键参数 
%     theta_max = wmax * sample_t;
%     cos_theta_max = cosd(theta_max);
%     e_dot_u = dot(e_vel, ustar);
%     e_dot_v = dot(e_vel, v_pre);
%     u_dot_v = dot(ustar, v_pre);
%     norm_e_sq = dot(e_vel, e_vel);
%     norm_u_sq = dot(ustar, ustar);
%     norm_v_sq = dot(v_pre, v_pre);
%     
%     % 1. 解范数约束区间 
%     % 上限约束：α²||u||² + 2α(e·u) + (||e||² -1) < 0 
%     a_upper = norm_u_sq;
%     b_upper = 2 * e_dot_u;
%     c_upper = norm_e_sq - 1;
%     roots_upper = roots([a_upper, b_upper, c_upper]);
%     I_upper = getInterval(roots_upper, a_upper, 'less');
%     
%     % 下限约束：α²||u||² + 2α(e·u) + (||e||² -0.5625) > 0 
%     a_lower = norm_u_sq;
%     b_lower = 2 * e_dot_u;
%     c_lower = norm_e_sq - 0.5625;
%     roots_lower = roots([a_lower, b_lower, c_lower]);
%     I_lower = getInterval(roots_lower, a_lower, 'greater');
%     
%     % 范数约束总区间 
%     I_norm = intersectIntervals(I_upper, I_lower);
%     
%     % 2. 解夹角约束区间 
%     % 构建二次不等式系数 
%     A = (u_dot_v)^2 - cos_theta_max^2 * norm_u_sq * norm_v_sq;
%     B = 2 * (e_dot_v * u_dot_v - cos_theta_max^2 * e_dot_u * norm_v_sq);
%     C = (e_dot_v)^2 - cos_theta_max^2 * norm_e_sq * norm_v_sq;
%     roots_angle = roots([A, B, C]);
%     I_angle = getInterval(roots_angle, A, 'greater');
%     
%     % 全局可行区间 
%     I_feasible = intersectIntervals(I_norm, I_angle);
%     
%     % 3. 选择最优α（优先非负，否则最小负值）
%     if isempty(I_feasible)
%         alpha = NaN;
%         ultra_ustar = NaN;
%         return;
%     end 
%     
%     candidates = [I_feasible(1), I_feasible(2), 0]; % 包含区间端点和零 
%     candidates = candidates(candidates >= I_feasible(1) & candidates <= I_feasible(2));
%     if isempty(candidates)
%         alpha = I_feasible(1); % 取左端点（最小负值）
%     else 
%         % 优先非负解 
%         non_neg = candidates(candidates >= 0);
%         if ~isempty(non_neg)
%             alpha = min(non_neg);
%         else 
%             alpha = max(candidates); % 最小绝对值负数 
%         end 
%     end 
%     
%     ultra_ustar = e_vel + alpha * ustar;
%     
%     % 验证范数约束 
%     if norm(ultra_ustar) < 0.75 || norm(ultra_ustar) > 1 
%         alpha = NaN;
%         ultra_ustar = NaN;
%     end 
% end 
%  
% % 辅助函数：计算二次不等式解区间 
% function interval = getInterval(roots, a, inequality)
%     real_roots = real(roots(abs(imag(roots)) < 1e-6));
%     real_roots = sort(real_roots);
%     if isempty(real_roots)
%         interval = [];
%         return;
%     end 
%     
%     if a > 0 
%         if strcmp(inequality, 'less')
%             % a > 0时，ax² + bx + c < 0 → 解在两个实根之间 
%             interval = [real_roots(1), real_roots(end)];
%         else 
%             % ax² + bx + c > 0 → 解在实根外 
%             interval = [-inf, real_roots(1), real_roots(end), inf];
%         end 
%     else 
%         if strcmp(inequality, 'less')
%             interval = [-inf, real_roots(1), real_roots(end), inf];
%         else 
%             interval = [real_roots(1), real_roots(end)];
%         end 
%     end 
% end 
%  
% % 辅助函数：区间交集 
% function result = intersectIntervals(I1, I2)
%     if isempty(I1) || isempty(I2)
%         result = [];
%         return;
%     end 
%     starts = max([I1(1), I2(1)]);
%     ends = min([I1(end), I2(end)]);
%     if starts > ends 
%         result = [];
%     else 
%         result = [starts, ends];
%     end 
% end 


%% 优化的方法求平滑ustar
% deprecated this shit always output unaccessable control
% function [alpha, ultra_ustar] = computeAlpha(e_vel, ustar, v_pre, sample_t, wmax)
%     % e_vel: 向量 0 < e_vel < 1
%     % ustar: 向量
%     % v_pre: 向量
%     % sample_t: 时间步长（标量）
%     % 返回值 alpha: 求得的系数
%     % 返回值 ultra_ustar: 最终的 ultra_ustar 向量
% 
%     % 计算夹角限制的余弦值
%     angle_limit = cosd(wmax * sample_t);
%     
%     % 初始猜测，尽量使得 α 为正
%     alpha_initial = 0;
%     
%     % 优化函数，目标是最大化 α，同时满足所有约束
%     objective = @(alpha) -alpha;  % 目标是最大化 α，所以最小化 -α
%     
%     % 初始 guess 和约束
%     options = optimset('Display', 'off');                % 不显示优化过程
%     alpha_bounds = [-10, 10]; % α 的边界范围
%     
%     % 使用 fmincon 求解优化问题
%     alpha = fmincon(objective, alpha_initial, [], [], [], [], alpha_bounds(1), alpha_bounds(2), ...
%                     @(alpha) constraints(alpha, e_vel, ustar, v_pre, angle_limit), options);
%     
%     % 计算最终的 ultra_ustar 向量
%     ultra_ustar = e_vel + alpha * ustar;
%     if norm(ultra_ustar) > 1 || norm(ultra_ustar) < 0.75
%         warning(" unaccessable control u length limit exceed");
%     end
%     cosuv = dot(ultra_ustar, v_pre) / (norm(ultra_ustar) * norm(v_pre));
%     if cosuv < angle_limit
%         warning(" unaccessable control u angle limit exceed ");
%     end
%     ultra_ustar = [];
% end
% 
% % 约束条件：夹角限制和模长限制
% function [c, ceq] = constraints(alpha, e_vel, ustar, v_pre, angle_limit)
%     ultra_ustar = e_vel + alpha * ustar;
%     norm_ultra_ustar = norm(ultra_ustar);
%     
%     % 计算夹角的余弦值
%     cosine_angle = dot(ultra_ustar, v_pre) / (norm_ultra_ustar * norm(v_pre));
%     
%     % 不等式约束
%     c = [
%         norm_ultra_ustar - 0.98;        % 模长大于 0.75
%         -norm_ultra_ustar + 0.77;   % 模长小于 1
%         -cosine_angle + angle_limit ; % 夹角限制
%     ];
%     
%     % 等式约束：没有
%     ceq = [];
% end
% 
% function ustar = computeSubOptustar(v_pre, dir, sample_t, wmax, lmax)
%     % 归一化输入向量 
%     v_pre = [v_pre, 0]; dir = [dir, 0];
%     v_pre_norm = v_pre / norm(v_pre);
%     dir_norm = dir / norm(dir);
%     
%     % 计算旋转轴n 
%     n = cross(v_pre_norm, dir_norm);
%     if norm(n) < 1e-6 
%         ustar = v_pre_norm; % 若v_pre与dir共线，直接返回v_pre 
%         return;
%     end 
%     n = n / norm(n);
%     
%     % 计算目标角度φ并验证可行性 
%     phi = wmax * sample_t;
%     theta = acosd(dot(v_pre_norm, dir_norm));
%     if phi > theta 
% %         warning('角度φ=%.2f°超过原始夹角θ=%.2f°，无法生成中间向量', phi, theta);
%         ustar = dir / norm(dir) * lmax;
%         ustar = ustar(1:2);
%         return;
%     end 
%     
%     % 构造旋转矩阵（绕n轴旋转φ角度）
%     R = @(angle) [cosd(angle) + n(1)^2*(1-cosd(angle)), ...
%                   n(1)*n(2)*(1-cosd(angle)) - n(3)*sind(angle), ...
%                   n(1)*n(3)*(1-cosd(angle)) + n(2)*sind(angle); ...
%                   n(2)*n(1)*(1-cosd(angle)) + n(3)*sind(angle), ...
%                   cosd(angle) + n(2)^2*(1-cosd(angle)), ...
%                   n(2)*n(3)*(1-cosd(angle)) - n(1)*sind(angle); ...
%                   n(3)*n(1)*(1-cosd(angle)) - n(2)*sind(angle), ...
%                   n(3)*n(2)*(1-cosd(angle)) + n(1)*sind(angle), ...
%                   cosd(angle) + n(3)^2*(1-cosd(angle))];
%     
%     % 生成ustar并归一化（可选）
%     ustar = R(phi) * v_pre_norm';
%     ustar = ustar / norm(ustar) * lmax; % 保持单位长度 
%     ustar = ustar(1:2)';
% end

%% 数学方法算

function ultra_ustar_p = computeUltraUstar(ustar_p, v_pre, sample_t, wmax, lmin)
    % 计算ustar_p的模长
    norm_u = norm(ustar_p);
    if norm_u == 0
        ultra_ustar_p = [0; 0];
        return;
    end
    
    % 计算当前方向
    theta = atan2(ustar_p(2), ustar_p(1));
    theta_v_pre = atan2(v_pre(2), v_pre(1));
    
    % 计算角度差
    delta_theta = theta - theta_v_pre;
    delta_theta = mod(delta_theta + pi, 2*pi) - pi;
    
    % 最大允许角度变化
    max_angle = sample_t * wmax / 180*pi;
    
    % 判断是否需要调整
    if abs(delta_theta) <= max_angle
        ultra_ustar_p = ustar_p;
        return;
    end
    
    % 调整角度差
    clamped_delta_theta = max(-max_angle, min(delta_theta, max_angle));
    
    % 新方向
    theta_new = theta_v_pre + clamped_delta_theta;
    
    % 构造新的向量
    ultra_ux = cos(theta_new);
    ultra_uy = sin(theta_new);
    
    ultra_ustar_p = [ultra_ux, ultra_uy] .* lmin;
end

