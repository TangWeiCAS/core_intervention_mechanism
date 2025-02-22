%%
clc;clear;

%CG1
matrix_CG1=zeros(400,400);%400*400矩阵
trildiag=load("submatrix_CG1.mat");
matrix_CG1(tril(true(400),-1))=mean(trildiag.CG1_harmonized,2);
matrix_CG1=matrix_CG1+matrix_CG1';

%TG1
matrix_TG1=zeros(400,400);%400*400矩阵
trildiag=load("submatrix_TG1.mat");
matrix_TG1(tril(true(400),-1))=mean(trildiag.TG1_harmonized,2);
matrix_TG1=matrix_TG1+matrix_TG1';

%CG2
matrix_CG2=zeros(400,400);%400*400矩阵
trildiag=load("submatrix_CG2.mat");
matrix_CG2(tril(true(400),-1))=mean(trildiag.CG2_harmonized,2);
matrix_CG2=matrix_CG2+matrix_CG2';

%TG2
matrix_TG2=zeros(400,400);%400*400矩阵
trildiag=load("submatrix_TG2.mat");
matrix_TG2(tril(true(400),-1))=mean(trildiag.TG2_harmonized,2);
matrix_TG2=matrix_TG2+matrix_TG2';

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
gm2 = gm.fit(matrix_TG1);
gm3 = gm.fit(matrix_CG2);
gm4 = gm.fit(matrix_TG2);
GM1=[gm1.gradients{1}(:,1), -1*gm1.gradients{1}(:,2), -1*gm1.gradients{1}(:,3)];
GM2=[gm2.gradients{1}(:,1), -1*gm2.gradients{1}(:,2), gm2.gradients{1}(:,3)];
GM3=[gm3.gradients{1}(:,1), -1*gm3.gradients{1}(:,2), -1*gm3.gradients{1}(:,3)];
GM4=[gm4.gradients{1}(:,1), gm4.gradients{1}(:,2), -1*gm4.gradients{1}(:,3)];
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
h2=plot_hemispheres(GM2,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'TG1_G_1','TG1_G_2','TG1_G_3'});
colormap(h2.handles.figure,[.4 .4 .4;colors]);
for i=1:3
    for j=1:4
        h2.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h2.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h2.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
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
h4=plot_hemispheres(GM4,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'TG2_G_1','TG2_G_2','TG2_G_3'});
colormap(h4.handles.figure,[.4 .4 .4;colors]);
for i=1:3
    for j=1:4
        h4.handles.axes(i,j).CLim=[-0.2 0.2];
    end
    h4.handles.colorbar(1,i).Ticks=[-0.2 0.2];
    h4.handles.colorbar(1,i).TickLabels=[-0.2 0.2];
end

[fig2,tempcolor]=gradient_surface([GM2(:,1) GM2(:,2)],{surf_lh,surf_rh},labeling,[]);
fig1=gradient_surface([GM1(:,1) GM1(:,2)],{surf_lh,surf_rh},labeling,tempcolor);
fig3=gradient_surface([GM3(:,1) GM3(:,2)],{surf_lh,surf_rh},labeling,tempcolor);
fig4=gradient_surface([GM4(:,1) GM4(:,2)],{surf_lh,surf_rh},labeling,tempcolor);