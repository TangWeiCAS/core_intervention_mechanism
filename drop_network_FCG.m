%% drop off
% 分网络drop，共7个不同的drop模式，针对干预组而非控制组
clear;clc;
cd C:\Users\tangw\Desktop\FCS\
load("vikO.mat");
%baseline
conn_matrix4 = zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_CG2.mat");
conn_matrix4(tril(true(400),-1))=mean(trildiag.CG2_harmonized,2);
conn_matrix4=conn_matrix4+conn_matrix4';
gcontrol = GradientMaps(kernel='na',approach='dm',n_components=10);
gcontrol = gcontrol.fit(conn_matrix4);%对400parcels进行拟合
labeling = load_parcellation('schaefer',400);
labeling = labeling.schaefer_400;
[surf_lh, surf_rh] = load_conte69('inflated');
CGgradients(:,1:2)=[gcontrol.gradients{1}(:,1) -1*gcontrol.gradients{1}(:,2)];%先在这里改正负号
h=plot_hemispheres(CGgradients(:,1:2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'Gradient 1','Gradient 2'});
colormap(h.handles.figure,[.4 .4 .4;vikO]);
for i=1:2
    for j=1:4
        h.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h.handles.colorbar(1,i).Ticks=[-0.2 0 0.2];
    h.handles.colorbar(1,i).TickLabels=[-0.2 0 0.2];
end
% interv ses-2
conn_matrix3 = zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_TG2.mat");
conn_matrix3(tril(true(400),-1))=mean(trildiag.TG2_harmonized,2);
conn_matrix3=conn_matrix3+conn_matrix3';

parcels=readtable("parcels.csv");%400parcels分别属于哪个网络
parcels=table2array(parcels(:,2));
network="FPCN";
cols_to_retain=find(parcels(:,1)~=network);
cols_to_remove=find(parcels(:,1)==network);
% 改成其他网络只需要修改名字即可
conn_matrix3=conn_matrix3(:,cols_to_retain);
conn_matrix3=conn_matrix3(cols_to_retain,:);

% 修改label进行绘图
% 找到label中等于要留下的labels中任何数字的索引
[ismember_idx, ~] = ismember(labeling, cols_to_remove);
% 提取出索引，这是所有不保留的DMN顶点的坐标信息
matched_idx = find(ismember_idx);
% 输出结果，凡是不保留的均归0
labeling(matched_idx(:,1),:)=0;

[ismember_idx2,~]=ismember(cols_to_retain,cols_to_retain);
matched_idx2=find(ismember_idx2);% 重新编号，剩下的编为1-3xx
for i=1:length(labeling(:,1))
    if labeling(i,1)~=0 % 还剩下的，就需要重新编号
        [change_idx,~]=find(cols_to_retain==labeling(i,1));% 原来这个位置的label排保留labels的第几位
        labeling(i,1)=matched_idx2(change_idx);% 改成那个位置的新编号
    else
    end
end

% 绘图
h = plot_hemispheres(labeling, {surf_lh,surf_rh});
colormap(h.handles.figure,[.4 .4 .4;vikO])%展示梯度所在脑区
gsub = GradientMaps(kernel='na',approach='dm',n_components=10);
gsub = gsub.fit(conn_matrix3);%对不重复的部分进行拟合
gradients(:,1:2)=[-1*gsub.gradients{1}(:,1) gsub.gradients{1}(:,2)];%先在这里改正负号
h_drop=plot_hemispheres(gradients(:,1:2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'Gradient 1','Gradient 2'});
colormap(h_drop.handles.figure,[.4 .4 .4;vikO]);
for i=1:2
    for j=1:4
        h_drop.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h_drop.handles.colorbar(1,i).Ticks=[-0.2 0 0.2];
    h_drop.handles.colorbar(1,i).TickLabels=[-0.2 0 0.2];
end

% 画解释方差
x = gsub.lambda{1};
x=100*x./sum(gsub.lambda{1});
x=x';
y = 1:length(x);
ind = zeros(1,length(x));
ind(1,x>=10)=1;
% 颜色定义
C1 = [0.7 0.7 0.7];
C2 = [217 116 43]./255;
C3 = [0.7 0.7 0.7];
% 棒棒
figure;
hold on;
GO = barh(y,x,0.05,'EdgeColor','none','FaceColor',C1);
xline(10,'LineStyle','--','Color',C2,'LineWidth',1.5);
% 圆圈
for i = 1:length(y)
    if ind(i) == 1
        scatter(x(i), y(i),100,C2,'filled')
    else
        scatter(x(i), y(i),100,C3,'filled')
    end
end
hXLabel = xlabel('Variance explained(%)','FontSize',14);
hYLabel = ylabel('number of components','FontSize',14);
hold off;
%%
% 相似性分析：哪个网络drop后与control的最像
% 将与训练组drop掉的parcel相同的位置归零，只提取想要的部分的FCG，便于计算相关性
control_gradients=CGgradients(cols_to_retain,1:2);

%检查梯度相似性
FCG_CG_ses2=control_gradients(:,2);
FCG_TG_ses2=gradients(:,2);
[sphere_lh, sphere_rh] = load_conte69('spheres');%这里是球面投影，不是surface
% 将数据在球面旋转打乱，看真实相关在随机相关里排多少
n_permutations = 1000; %1000次重取样
y_rand = spin_permutations({labeling(1:32492,1),labeling(32493:64984,1)}, ...
                  {sphere_lh,sphere_rh}, ...
                  n_permutations,'random_state',0);%旋转1000次
parcel_rotated=squeeze([y_rand{1}(:,1,:); y_rand{2}(:,1,:)]);
ExchangeLabel=zeros(64984,n_permutations);%没找到的归为0，在brainspace就是空皮质
FCG_TG_parcel_rotated=nan(length(gradients(:,1)),n_permutations);
for perm=1:n_permutations
    for i=1:length(gradients(:,1))
        OverlapElements = parcel_rotated(find(labeling(:,1)==i),perm);
        B = mode(OverlapElements);%找到重叠parcel的最多的label，重新赋值并计算相似性
        ExchangeLabel(find(parcel_rotated(:,perm)==B),perm)=i;%旋转后B与原来的label1覆盖最多，将B parcel赋值为label1，以此类推
        %目标是代替原位置的label参与相关性
        if B~=0
            FCG_TG_parcel_rotated(i,perm)=FCG_TG_ses2(B,1);
        end
    end
end
[r_origin, ~] = corr(FCG_CG_ses2,FCG_TG_ses2, ...
                'rows','pairwise','type','Pearson');%先把G12的轴确定好，不要出现两者相反的情况，正相关有意义
[r_rand,~] = corr(FCG_CG_ses2,FCG_TG_parcel_rotated, ...
            'rows','pairwise','type','Pearson');%以覆盖率最大的代表原始parcel，计算parcel级别相关性，拒绝过度效应
num=sum(logical(r_origin<r_rand));%有多少随机数比原始相关值大
prctile_rank = mean(r_origin > r_rand);
significant = prctile_rank >= 0.95;%是否显著
figure;
histogram(r_rand,FaceColor=[0.7,0.7,0.7],FaceAlpha=0.5);
xlim([0 1]);
ylim([0 200]);
hold on;  % 保持当前图形，以便在其上添加新的线条
% 添加一条线
xline(r_origin, 'k', 'LineWidth', 2);  % LineWidth设置线宽为2
hold off;  
title('Pearson corr between control and intervention group DCs');
xlabel('permutation r');
ylabel('Frequency');
% writematrix(gradients,'TG_drop_gradients.csv');
% writematrix(control_gradients,'baseline_drop_gradients.csv');
%% 相似性
distance=sqrt(sum((gradients(:,1:2) - control_gradients(:,1:2)).^2, 2));%二维距离
normal_distance=sum(distance,1)./length(cols_to_retain);
disp(['distance is ', num2str(normal_distance)]);