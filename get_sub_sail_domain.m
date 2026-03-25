function dom = get_sub_sail_domain(prob_fail)
    % 生成单个子帆面合法的状态集合
    if rand() < prob_fail
        dom = [randi([0, 1])]; % 发生卡死，固定在 0 或 1
    else
        dom = [0, 1];          % 健康，可取 0 或 1
    end
end