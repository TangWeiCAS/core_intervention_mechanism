clc;clear;
%%group based on median of DA FCSwithin (TG2-TG1)
diffdata=readmatrix("correlation_diff.csv");%17列SMN，18DA，19DMN

%高增强组和低改变组，高组应该翻转，低组不变
matrix_high=zeros(400,400);%单被试400*400矩阵
matrix_low=zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_TG2.mat");
DA_high=trildiag.TG2_harmonized(:,find(diffdata(:,18)==1));
DA_low=trildiag.TG2_harmonized(:,find(diffdata(:,18)==0));
matrix_high(tril(true(400),-1))=mean(DA_high,2);
matrix_low(tril(true(400),-1))=mean(DA_low,2);
matrix_high=matrix_high+matrix_high';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵
matrix_low=matrix_low+matrix_low';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵

%梯度对照组：CG2
matrix_CG2=zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_CG2.mat");
matrix_CG2(tril(true(400),-1))=mean(trildiag.CG2_harmonized,2);
matrix_CG2=matrix_CG2+matrix_CG2';

%G1-G3梯度绘制
load("vikO.mat")
labeling = load_parcellation('schaefer',400);
labeling = labeling.schaefer_400;
[surf_lh, surf_rh] = load_conte69('inflated');
gm = GradientMaps(kernel='na',approach='dm',n_components=10);
gm1 = gm.fit(matrix_high);%高
gm2 = gm.fit(matrix_low);%低
gm3 = gm.fit(matrix_CG2);
h1=plot_hemispheres([gm1.gradients{1}(:,1) gm2.gradients{1}(:,1) gm3.gradients{1}(:,1)],{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_1','low SMN_G_1','CG2_G_1'});
colormap(h1.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h1.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h1.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h1.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h2=plot_hemispheres([-1*gm1.gradients{1}(:,2) -1*gm2.gradients{1}(:,2) -1*gm3.gradients{1}(:,2)],{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_2','low SMN_G_2','CG2_G_2'});
colormap(h2.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h2.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h2.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h2.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h3=plot_hemispheres([gm1.gradients{1}(:,3) gm2.gradients{1}(:,3) gm3.gradients{1}(:,3)],{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_3','low SMN_G_3','CG2_G_3'});
colormap(h3.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h3.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h3.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h3.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
%%
% 1.1.2 基于Baum法的400 parcels spin test
%检验这差值分组G1和G2的相似性，可以说明两组G1/G2分布是否一致（没有大的改变）
%主要还是要区分high和low组
TG=GM2(:,2);
CG=GM2(:,3);
[sphere_lh, sphere_rh] = load_conte69('spheres');%这里是球面投影，不是surface
% 将数据在球面旋转打乱，看真实相关在随机相关里排多少
n_permutations = 1000; %1000次重取样
y_rand = spin_permutations({labeling(1:32492,1),labeling(32493:64984,1)}, ...
                  {sphere_lh,sphere_rh}, ...
                  n_permutations,'random_state',0);%旋转1000次
parcel_rotated=squeeze([y_rand{1}(:,1,:); y_rand{2}(:,1,:)]);
ExchangeLabel=zeros(64984,n_permutations);%没找到的归为0，在brainspace就是空皮质
TG_parcel_rotated=nan(400,n_permutations);
for perm=1:n_permutations
    for i=1:400
        OverlapElements = parcel_rotated(find(labeling(:,1)==i),perm);
        B = mode(OverlapElements);%找到重叠parcel的最多的label，重新赋值并计算相似性
        ExchangeLabel(find(parcel_rotated(:,perm)==B),perm)=i;%旋转后B与原来的label1覆盖最多，将B parcel赋值为label1，以此类推
        %目标是代替原位置的label参与相关性
        if B~=0
            TG_parcel_rotated(i,perm)=TG(B,1);
        end
    end
end
[r_origin, ~] = corr(CG,TG, ...
                'rows','pairwise','type','Pearson');
[r_rand,~] = corr(CG,TG_parcel_rotated, ...
            'rows','pairwise','type','Pearson');%以覆盖率最大的代表原始parcel，计算parcel级别相关性，拒绝过度效应
num=sum(r_origin<r_rand);%有多少随机数比原始相关值大
prctile_rank = mean(r_origin > r_rand);
significant = num/n_permutations;%是否显著
figure;
histogram(r_rand,FaceColor=[0.7,0.7,0.7],FaceAlpha=0.5);
xlim([0 1]);
ylim([0 200]);
hold on;  % 保持当前图形，以便在其上添加新的线条
% 添加一条线
xline(r_origin, 'k', 'LineWidth', 2);  % LineWidth设置线宽为2
hold off;  
title('Pearson corr between control and intervention group');
xlabel('permutation r');
ylabel('Frequency');
%%
clc;clear;
%%group based on DA FCSwithin (TG2)
TG2data=readmatrix("correlation_TG2.csv");%34列SMN，35DA，36DMN

%高增强组和低改变组，高组应该翻转，低组不变
matrix_high=zeros(400,400);%单被试400*400矩阵
matrix_low=zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_TG2.mat");
DA_high=trildiag.TG2_harmonized(:,find(TG2data(:,36)==1));
DA_low=trildiag.TG2_harmonized(:,find(TG2data(:,36)==0));
matrix_high(tril(true(400),-1))=mean(DA_high,2);
matrix_low(tril(true(400),-1))=mean(DA_low,2);
matrix_high=matrix_high+matrix_high';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵
matrix_low=matrix_low+matrix_low';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵

%梯度对照组：CG2
matrix_CG2=zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_CG2.mat");
matrix_CG2(tril(true(400),-1))=mean(trildiag.CG2_harmonized,2);
matrix_CG2=matrix_CG2+matrix_CG2';

%G1-G3梯度绘制
load("vikO.mat")
labeling = load_parcellation('schaefer',400);
labeling = labeling.schaefer_400;
[surf_lh, surf_rh] = load_conte69('inflated');
gm = GradientMaps(kernel='na',approach='dm',n_components=10);
gm1 = gm.fit(matrix_high);
gm2 = gm.fit(matrix_low);
gm3 = gm.fit(matrix_CG2);
h1=plot_hemispheres([gm1.gradients{1}(:,1) gm2.gradients{1}(:,1) gm3.gradients{1}(:,1)],{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high DA_G_1','low DA_G_1','CG2_G_1'});
colormap(h1.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h1.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h1.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h1.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h2=plot_hemispheres([gm1.gradients{1}(:,2) -1*gm2.gradients{1}(:,2) -1*gm3.gradients{1}(:,2)],{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high DA_G_2','low DA_G_2','CG2_G_2'});
colormap(h2.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h2.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h2.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h2.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h3=plot_hemispheres([-1*gm1.gradients{1}(:,3) gm2.gradients{1}(:,3) gm3.gradients{1}(:,3)],{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high DA_G_3','low DA_G_3','CG2_G_3'});
colormap(h3.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h3.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h3.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h3.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end

%%
% clc;clear;
%%group based on median of DA FCSwithin (TG2-TG1)
diffdata=readmatrix("correlation_diff.csv");%18列SMN，19DA，20DMN

%高增强组和低改变组，高组应该翻转，低组不变
matrix_high=zeros(400,400);%单被试400*400矩阵
matrix_low=zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_TG2.mat");
DA_high=trildiag.TG2_harmonized(:,find(diffdata(:,19)==1));
DA_low=trildiag.TG2_harmonized(:,find(diffdata(:,19)==0));
matrix_high(tril(true(400),-1))=mean(DA_high,2);
matrix_low(tril(true(400),-1))=mean(DA_low,2);
matrix_high=matrix_high+matrix_high';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵
matrix_low=matrix_low+matrix_low';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵

%梯度对照组：TG1
matrix_CG2=zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_TG1.mat");
matrix_CG2(tril(true(400),-1))=mean(trildiag.TG1_harmonized,2);
matrix_CG2=matrix_CG2+matrix_CG2';

%G1-G3梯度绘制
load("vikO.mat")
labeling = load_parcellation('schaefer',400);
labeling = labeling.schaefer_400;
[surf_lh, surf_rh] = load_conte69('inflated');
gm = GradientMaps(kernel='na',approach='dm',n_components=10);
gm1 = gm.fit(matrix_high);%高
gm2 = gm.fit(matrix_low);%低
gm3 = gm.fit(matrix_CG2);
GM1=[-1*gm1.gradients{1}(:,1), gm2.gradients{1}(:,1), gm3.gradients{1}(:,1)];
GM2=[-1*gm1.gradients{1}(:,2) gm2.gradients{1}(:,2) -1*gm3.gradients{1}(:,2)];
GM3=[-1*gm1.gradients{1}(:,3) gm2.gradients{1}(:,3) gm3.gradients{1}(:,3)];
h1=plot_hemispheres(GM1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high DAN_G_1','low DAN_G_1','TG1_G_1'});
colormap(h1.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h1.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h1.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h1.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h2=plot_hemispheres(GM2,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high DAN_G_2','low DAN_G_2','TG1_G_2'});
colormap(h2.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h2.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h2.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h2.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h3=plot_hemispheres(GM3,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high DAN_G_3','low DAN_G_3','TG1_G_3'});
colormap(h3.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h3.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h3.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h3.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
[fig,color]=gradient_surface([GM1(:,3) GM2(:,3)],{surf_lh,surf_rh},labeling,[]);
fig1=gradient_surface([GM1(:,1) GM2(:,1)],{surf_lh,surf_rh},labeling,color);
fig2=gradient_surface([GM1(:,2) GM2(:,2)],{surf_lh,surf_rh},labeling,color);
%% compared to baseine
%%
clc;clear;
%%group based on median of DA FCSwithin (TG2-TG1)
diffdata=readmatrix("correlation_diff.csv");%18列SMN，19DA，20DMN

%高增强组和低改变组，高组应该翻转，低组不变
matrix_high=zeros(400,400);%单被试400*400矩阵
matrix_low=zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_TG2.mat");
DA_high=trildiag.TG2_harmonized(:,find(diffdata(:,20)==1));
DA_low=trildiag.TG2_harmonized(:,find(diffdata(:,20)==0));
matrix_high(tril(true(400),-1))=mean(DA_high,2);
matrix_low(tril(true(400),-1))=mean(DA_low,2);
matrix_high=matrix_high+matrix_high';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵
matrix_low=matrix_low+matrix_low';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵

%梯度对照组：baseline
matrix_baseline=zeros(400,400);
trildiag=load("matrix_oldaverage.mat");
matrix_baseline=trildiag.FCzmatrix_average;

%G1-G3梯度绘制
load("vikO.mat")
labeling = load_parcellation('schaefer',400);
labeling = labeling.schaefer_400;
[surf_lh, surf_rh] = load_conte69('inflated');
gm = GradientMaps(kernel='na',approach='dm',n_components=10);
gm1 = gm.fit(matrix_high);%高
gm2 = gm.fit(matrix_low);%低
gm3 = gm.fit(matrix_baseline);
GM1=[gm1.gradients{1}(:,1), -1*gm2.gradients{1}(:,1), gm3.gradients{1}(:,1)];
GM2=[-1*gm1.gradients{1}(:,2) gm2.gradients{1}(:,2) gm3.gradients{1}(:,2)];
GM3=[gm1.gradients{1}(:,3) gm2.gradients{1}(:,3) -1*gm3.gradients{1}(:,3)];
h1=plot_hemispheres(GM1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_1','low SMN_G_1','baseline_G_1'});
colormap(h1.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h1.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h1.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h1.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h2=plot_hemispheres(GM2,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_2','low SMN_G_2','baseline_G_2'});
colormap(h2.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h2.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h2.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h2.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h3=plot_hemispheres(GM3,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_3','low SMN_G_3','baseline_G_3'});
colormap(h3.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h3.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h3.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h3.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
[fig,color]=gradient_surface([GM1(:,3) GM2(:,3)],{surf_lh,surf_rh},labeling,[]);
fig1=gradient_surface([GM1(:,1) GM2(:,1)],{surf_lh,surf_rh},labeling,color);
fig2=gradient_surface([GM1(:,2) GM2(:,2)],{surf_lh,surf_rh},labeling,color);
%%
clc;clear;
%%group based on median of DA FCSwithin TG2
diffdata=readmatrix("correlation_TG2.csv");%34列SMN，35DA，36DMN

%高增强组和低改变组，高组应该翻转，低组不变
matrix_high=zeros(400,400);%单被试400*400矩阵
matrix_low=zeros(400,400);%单被试400*400矩阵
trildiag=load("submatrix_TG2.mat");
DA_high=trildiag.TG2_harmonized(:,find(diffdata(:,36)==1));
DA_low=trildiag.TG2_harmonized(:,find(diffdata(:,36)==0));
matrix_high(tril(true(400),-1))=mean(DA_high,2);
matrix_low(tril(true(400),-1))=mean(DA_low,2);
matrix_high=matrix_high+matrix_high';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵
matrix_low=matrix_low+matrix_low';%去掉diag的下三角（对角归零）一列400*399/2重新恢复为400*400的对称矩阵

%梯度对照组：baseline
matrix_baseline=zeros(400,400);
trildiag=load("matrix_oldaverage.mat");
matrix_baseline=trildiag.FCzmatrix_average;

%G1-G3梯度绘制
load("vikO.mat")
labeling = load_parcellation('schaefer',400);
labeling = labeling.schaefer_400;
[surf_lh, surf_rh] = load_conte69('inflated');
gm = GradientMaps(kernel='na',approach='dm',n_components=10);
gm1 = gm.fit(matrix_high);%高
gm2 = gm.fit(matrix_low);%低
gm3 = gm.fit(matrix_baseline);
GM1=[gm1.gradients{1}(:,1), gm2.gradients{1}(:,1), gm3.gradients{1}(:,1)];
GM2=[gm1.gradients{1}(:,2) -1*gm2.gradients{1}(:,2) gm3.gradients{1}(:,2)];
GM3=[gm1.gradients{1}(:,3) -1*gm2.gradients{1}(:,3) -1*gm3.gradients{1}(:,3)];
h1=plot_hemispheres(GM1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_1','low SMN_G_1','baseline_G_1'});
colormap(h1.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h1.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h1.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h1.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h2=plot_hemispheres(GM2,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_2','low SMN_G_2','baseline_G_2'});
colormap(h2.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h2.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h2.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h2.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h3=plot_hemispheres(GM3,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'high SMN_G_3','low SMN_G_3','baseline_G_3'});
colormap(h3.handles.figure,[.4 .4 .4;vikO]);
for i=1:3
    for j=1:4
        h3.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h3.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h3.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
[fig,color]=gradient_surface([GM1(:,3) GM2(:,3)],{surf_lh,surf_rh},labeling,[]);
fig1=gradient_surface([GM1(:,1) GM2(:,1)],{surf_lh,surf_rh},labeling,color);
fig2=gradient_surface([GM1(:,2) GM2(:,2)],{surf_lh,surf_rh},labeling,color);