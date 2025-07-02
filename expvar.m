%%
clc;clear;

%CG1
matrix_CG1=zeros(400,400);%400*400矩阵
trildiag=load("submatrix_CG1.mat");
matrix_CG1(tril(true(400),-1))=mean(trildiag.CG1_harmonized,2);
matrix_CG1=matrix_CG1+matrix_CG1';

%CG2
matrix_CG2=zeros(400,400);%400*400矩阵
trildiag=load("submatrix_CG2.mat");
matrix_CG2(tril(true(400),-1))=mean(trildiag.CG2_harmonized,2);
matrix_CG2=matrix_CG2+matrix_CG2';

%G1-G2梯度绘制
colours=load("vikO.mat");
len=length(colours.vikO(:,1));
indices = linspace(0, len-1, 400);
% 对索引进行整数化处理，并使用线性内插法插值来计算新的颜色值
colors = interp1(1:len, colours.vikO, indices, 'linear', 'extrap');
colors(colors<0)=0;

labeling = load_parcellation('schaefer',400);
labeling = labeling.schaefer_400;
[surf_lh, surf_rh] = load_conte69('inflated');

gm = GradientMaps(kernel='na',approach='dm',n_components=10);
gm1 = gm.fit(matrix_CG1);
gm3 = gm.fit(matrix_CG2);
GM1=[gm1.gradients{1}(:,1), -1*gm1.gradients{1}(:,2), -1*gm1.gradients{1}(:,3)];
GM3=[gm3.gradients{1}(:,1), -1*gm3.gradients{1}(:,2), -1*gm3.gradients{1}(:,3)];
h1=plot_hemispheres(GM1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'CG1_G_1','CG1_G_2','CG1_G_3'});
colormap(h1.handles.figure,[.4 .4 .4;colors]);
for i=1:3
    for j=1:4
        h1.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h1.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h1.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end
h3=plot_hemispheres(GM3,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'CG2_G_1','CG2_G_2','CG2_G_3'});
colormap(h3.handles.figure,[.4 .4 .4;colors]);
for i=1:3
    for j=1:4
        h3.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h3.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h3.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end

[fig1,tempcolor]=gradient_surface([GM1(:,1) GM1(:,2)],{surf_lh,surf_rh},labeling,[]);
fig3=gradient_surface([GM3(:,1) GM3(:,2)],{surf_lh,surf_rh},labeling,tempcolor);

% 画解释方差
x = gm1.lambda{1};
%x = gm3.lambda{1};
x=100*x./sum(gm1.lambda{1});
%x=100*x./sum(gm3.lambda{1});
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
%% permutation test (exp. var.)
clc;clear;
trildiag_CG1 = load("submatrix_CG1.mat");
trildiag_CG2 = load("submatrix_CG2.mat");

all_data = [trildiag_CG1.CG1_harmonized, trildiag_CG2.CG2_harmonized];
[n_rows, n_cols] = size(all_data); % 79800×118
num_permutations = 1000; % 置换次数
perm_results = cell(num_permutations, 2); % 存储置换结果

for i = 1:num_permutations
    shuffled_indices = randperm(n_cols);
    
    group1_indices = shuffled_indices(1:59);
    group2_indices = shuffled_indices(60:118);

    matrix1 = all_data(:, group1_indices);
    matrix2 = all_data(:, group2_indices);
    
    perm_results{i, 1} = matrix1;
    perm_results{i, 2} = matrix2;
end
d_rand=zeros(2,num_permutations);
for i=1:1000
    example_matrix1 = perm_results{i, 1};
    example_matrix2 = perm_results{i, 2};

    matrix_CG1_perm = zeros(400,400);
    matrix_CG1_perm(tril(true(400),-1)) = mean(example_matrix1, 2);
    matrix_CG1_perm = matrix_CG1_perm + matrix_CG1_perm';

    matrix_CG2_perm = zeros(400,400);
    matrix_CG2_perm(tril(true(400),-1)) = mean(example_matrix2, 2);
    matrix_CG2_perm = matrix_CG2_perm + matrix_CG2_perm';

    gm = GradientMaps(kernel='na',approach='dm',n_components=10);
    gm1 = gm.fit(matrix_CG1_perm);
    gm2 = gm.fit(matrix_CG2_perm);
    var_CG1 = gm1.lambda{1};
    var_CG2 = gm2.lambda{1};
    var_CG1=100*var_CG1./sum(gm1.lambda{1});
    var_CG2=100*var_CG2./sum(gm2.lambda{1});

    d1=var_CG1(1)-var_CG2(2);
    d_rand(i,1)=d1;
    d2=var_CG2(1)-var_CG1(2);
    d_rand(i,2)=d2;
end

d1_origin = 0.6539;%primary gradient
d2_origin = 0.0749;%secondary gradient
% num=sum(d1_origin<d_rand(:,1));
% prctile_rank = mean(d1_origin > d_rand(:,1));
num=sum(d2_origin<d_rand(:,2));
prctile_rank = mean(d2_origin > d_rand(:,2));
significant = num/num_permutations;
figure;
% histogram(d_rand(:,1),FaceColor=[0.7,0.7,0.7],FaceAlpha=0.5);
histogram(d_rand(:,2),FaceColor=[0.7,0.7,0.7],FaceAlpha=0.5);
xlim([-2 2]);
ylim([0 150]);
hold on;
% xline(d1_origin, 'k', 'LineWidth', 2);
xline(d2_origin, 'k', 'LineWidth', 2);
xline(median(d_rand(:,2)),'LineStyle','--','Color',[217 116 43]./255,'LineWidth', 2);
hold off;  
title('permutation test');
xlabel('explained variance difference');
ylabel('Frequency');