%% step1.1 calculate the gradient of old group by using the group-averaged FC
load ('D:\data\LJ\results\final\networklabel_schaefer.mat');netwok=network(1:200,:);
labeling=load_parcellation('schaefer',400);labeling=labeling.schaefer_400;
sub_matrix = g_ls('D:\data\LJ\results\final\Training\0.5FD\FC\*\*\ses-1\*.mat');
for i=1:length(sub_matrix);
    [a,b,c]= fileparts(sub_matrix{i});
    subid(i)=str2double(regexp(b,'\d+','match'));
    submatrix=load(sub_matrix{i});submatrix=submatrix.r;
    trimatrix=tril(submatrix,-1);trimatrix=trimatrix(trimatrix~=0);
    sub_zmatrix=0.5.* reallog((1+trimatrix)./(1-trimatrix));
    sub_zmatrix(isnan(sub_zmatrix))=1e-10;
    subzmatrix(:,i)=sub_zmatrix;
end

load('D:\data\LJ\results\final\behav\old_baseline.mat');
age=Behav(:,5);sex=dummyvar(Behav(:,6));edu=Behav(:,8);fd=Behav(:,9);
% Neuorcombat
batch(1,1:31)=3;batch(1,32:72)=4;batch(1,73:104)=5;batch(1,105:190)=6;batch(1,191:217)=7;
mod=[age sex(:,1) edu fd];
FC_harmonized=combat(subzmatrix,batch,mod,1);FC_average=ZCX_no0mean(FC_harmonized,2);
FCzmatrix_average=zeros(400,400);
FCzmatrix_average(tril(true(400),-1))=FC_average;FCzmatrix_average=FCzmatrix_average+FCzmatrix_average';
SurfAvgInfL='D:\data\HCP_S1200_GroupAvg_v1\S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii';
SurfAvgInfR='D:\data\HCP_S1200_GroupAvg_v1\S1200.R.very_inflated_MSMAll.32k_fs_LR.surf.gii';
sL=gifti(SurfAvgInfL);  sR=gifti(SurfAvgInfR);
surf_l.coord=sL.vertices';surf_l.tri=sL.faces;
surf_r.coord=sR.vertices';surf_r.tri=sR.faces;
surf.coord=[surf_l.coord,surf_r.coord];surf.tri=[surf_l.tri;surf_r.tri+32492];
gm=GradientMaps('Kernel','na','approach','dm','alignment','none','n_components',10,'verbose',true);gm=gm.fit(FCzmatrix_average,'sparsity',90,'tolerance',1e1);
G_baseline=plot_hemispheres([gm.gradients{1,1}(:,1)*-1,gm.gradients{1,1}(:,2)*-1,gm.gradients{1,1}(:,3),gm.gradients{1,1}(:,4)*-1],{surf_l,surf_r},'parcellation',labeling,'labeltext',{'Gradient 1','Gradient 2','Gradient 3','Gradient 4'});
colormap(G_baseline.handles.figure,[.5 .5 .5;parula]);
scree_plot(gm.lambda{1});
gradient_in_euclidean([gm.gradients{1,1}(:,1)*-1,gm.gradients{1,1}(:,2)*-1],{surf_l,surf_r},labeling);

%% step 1.2 calculate the gradient of CG2/CG1 by using the group-averaged FC
load('D:\data\LJ\results\final\Training\0.5FD\FC\3CT\ID_C12_CT.mat');load('D:\data\LJ\results\final\Training\0.5FD\FC\3CT\ID_T12_CT.mat');
load('D:\data\LJ\results\final\Training\0.5FD\FC\4MI\ID_C12_MI.mat');load('D:\data\LJ\results\final\Training\0.5FD\FC\4MI\ID_T12_MI.mat');
load('D:\data\LJ\results\final\Training\0.5FD\FC\5Dance\ID_C12_Dance.mat');load('D:\data\LJ\results\final\Training\0.5FD\FC\5Dance\ID_T12_Dance.mat');
[ID_CTc,is_CTc,IS_CTc]=intersect(Behav(1:104,3),ID_C12_CT);[ID_CTt,is_CTt,IS_CTt]=intersect(Behav(1:104,3),ID_T12_CT);
[ID_MIc,is_MIc,IS_MIc]=intersect(Behav(1:104,3),ID_C12_MI);[ID_MIt,is_MIt,IS_MIt]=intersect(Behav(1:104,3),ID_T12_MI);
[ID_Dc,is_Dc,IS_Dc]=intersect(Behav(1:104,3),ID_C12_Dance);[ID_Dt,is_Dt,IS_Dt]=intersect(Behav(1:104,3),ID_T12_Dance);

CG12batch(1,1:12)=3;CG12batch(1,13:25)=4;CG12batch(1,26:38)=5;CG12batch(1,39:59)=6;
behav_CG2=[Behav(is_CTc,:);Behav(is_MIc,:);Behav(is_Dc,:);Behav(105:125,:)];
age_CG2=behav_CG2(:,5);sex_CG2=dummyvar(behav_CG2(:,6));edu_CG2=behav_CG2(:,8);fd_CG2=behav_CG2(:,10);fd_CG1=behav_CG2(:,9);
mod_CG2=[age_CG2 sex_CG2(:,1) edu_CG2 fd_CG2];mod_CG1=[age_CG2 sex_CG2(:,1) edu_CG2 fd_CG1];

CG2_matrix = g_ls('D:\data\LJ\results\final\Training\0.5FD\FC\*\CG\ses-2\*.mat');
for i=1:length(CG2_matrix);
CG2matrix=load(CG2_matrix{i});CG2matrix=CG2matrix.r;
CG2trimatrix=tril(CG2matrix,-1);CG2trimatrix=CG2trimatrix(CG2trimatrix~=0);
CG2_zmatrix=0.5.* reallog((1+CG2trimatrix)./(1-CG2trimatrix));
CG2_zmatrix(isnan(CG2_zmatrix))=1e-10;
CG2zmatrix(:,i)=CG2_zmatrix;
end

CG2_harmonized=combat(CG2zmatrix,CG12batch,mod_CG2,1);CG2_average=ZCX_no0mean(CG2_harmonized,2);
CG2zmatrix_average=zeros(400,400);
CG2zmatrix_average(tril(true(400),-1))=CG2_average;CG2zmatrix_average=CG2zmatrix_average+CG2zmatrix_average';
gm2=GradientMaps('Kernel','na','approach','dm','alignment','none','n_components',10,'verbose',true);gm2=gm2.fit(CG2zmatrix_average,'sparsity',90,'tolerance',1e1);
G_CG2=plot_hemispheres([gm2.gradients{1,1}(:,1)*-1,gm2.gradients{1,1}(:,2),gm2.gradients{1,1}(:,3)*-1],{surf_l,surf_r},'parcellation',labeling,'labeltext',{'Gradient 1','Gradient 2','Gradient 3'});
colormap(G_CG2.handles.figure,[.5 .5 .5;parula]);
scree_plot(gm2.lambda{1});
gradient_in_euclidean([gm2.gradients{1,1}(:,1)*-1,gm2.gradients{1,1}(:,2)],{surf_l,surf_r},labeling);

CG1zmatrix=[subzmatrix(:,is_CTc),subzmatrix(:,is_MIc),subzmatrix(:,is_Dc),subzmatrix(:,105:125)];
CG1_harmonized=combat(CG1zmatrix,CG12batch,mod_CG1,1);CG1_average=ZCX_no0mean(CG1_harmonized,2);
CG1zmatrix_average=zeros(400,400);
CG1zmatrix_average(tril(true(400),-1))=CG1_average;CG1zmatrix_average=CG1zmatrix_average+CG1zmatrix_average';
gm1=GradientMaps('Kernel','na','approach','dm','alignment','none','n_components',10,'verbose',true);gm1=gm1.fit(CG1zmatrix_average,'sparsity',90,'tolerance',1e1);
G_CG1=plot_hemispheres([gm1.gradients{1,1}(:,1)*-1,gm1.gradients{1,1}(:,2)*-1,gm1.gradients{1,1}(:,3)*-1],{surf_l,surf_r},'parcellation',labeling,'labeltext',{'Gradient 1','Gradient 2','Gradient 3'});
colormap(G_CG1.handles.figure,[.5 .5 .5;parula]);
scree_plot(gm1.lambda{1});
gradient_in_euclidean([gm1.gradients{1,1}(:,1)*-1,gm1.gradients{1,1}(:,2)],{surf_l,surf_r},labeling);
%change color
[fig1,color]=gradient_surface([gm1.gradients{1,1}(:,1)*-1,gm1.gradients{1,1}(:,2)],{surf_l,surf_r},labeling,[]);
fig2=gradient_surface([gm2.gradients{1,1}(:,1)*-1,gm2.gradients{1,1}(:,2)],{surf_l,surf_r},labeling,color);

%% step 1.3 calculate the gradient of TG2/TG1 by using the group-averaged FC
TG12batch(1,1:16)=3;TG12batch(1,17:32)=4;TG12batch(1,33:47)=5;TG12batch(1,48:112)=6;
behav_TG2=[Behav(is_CTt,:);Behav(is_MIt,:);Behav(is_Dt,:);Behav(126:190,:)];ID_TG2=behav_TG2(:,3);
age_TG2=behav_TG2(:,5);sex_TG2=dummyvar(behav_TG2(:,6));edu_TG2=behav_TG2(:,8);fd_TG2=behav_TG2(:,10);fd_TG1=behav_TG2(:,9);
mod_TG2=[age_TG2 sex_TG2(:,1) edu_TG2 fd_TG2];mod_TG1=[age_TG2 sex_TG2(:,1) edu_TG2 fd_TG1];

TG2_matrix = g_ls('D:\data\LJ\results\final\Training\0.5FD\FC\*\TG\ses-2\*.mat');
for i=1:length(TG2_matrix);
TG2matrix=load(TG2_matrix{i});TG2matrix=TG2matrix.r;
TG2trimatrix=tril(TG2matrix,-1);TG2trimatrix=TG2trimatrix(TG2trimatrix~=0);
TG2_zmatrix=0.5.* reallog((1+TG2trimatrix)./(1-TG2trimatrix));
TG2_zmatrix(isnan(TG2_zmatrix))=1e-10;
TG2zmatrix(:,i)=TG2_zmatrix;
end
TG2_harmonized=combat(TG2zmatrix,TG12batch,mod_TG2,1);TG2_average=ZCX_no0mean(TG2_harmonized,2);
TG2zmatrix_average=zeros(400,400);
TG2zmatrix_average(tril(true(400),-1))=TG2_average;TG2zmatrix_average=TG2zmatrix_average+TG2zmatrix_average';
gm=GradientMaps('Kernel','na','approach','dm','alignment','none','n_components',10,'verbose',true);gm=gm.fit(TG2zmatrix_average,'sparsity',90,'tolerance',1e1);
G_TG2=plot_hemispheres([gm.gradients{1,1}(:,1)*-1,gm.gradients{1,1}(:,2),gm.gradients{1,1}(:,3)*-1],{surf_l,surf_r},'parcellation',labeling,'labeltext',{'Gradient 1','Gradient 2','Gradient 3'});
colormap(G_TG2.handles.figure,[.5 .5 .5;parula]);
scree_plot(gm.lambda{1});
gradient_in_euclidean([gm.gradients{1,1}(:,1)*-1,gm.gradients{1,1}(:,2)],{surf_l,surf_r},labeling);

TG1zmatrix=[subzmatrix(:,is_CTt),subzmatrix(:,is_MIt),subzmatrix(:,is_Dt),subzmatrix(:,126:190)];
TG1_harmonized=combat(TG1zmatrix,TG12batch,mod_TG1,1);TG1_average=ZCX_no0mean(TG1_harmonized,2);
TG1zmatrix_average=zeros(400,400);
TG1zmatrix_average(tril(true(400),-1))=TG1_average;TG1zmatrix_average=TG1zmatrix_average+TG1zmatrix_average';
gm=GradientMaps('Kernel','na','approach','dm','alignment','none','n_components',10,'verbose',true);gm=gm.fit(TG1zmatrix_average,'sparsity',90,'tolerance',1e1);
G_CG1=plot_hemispheres([gm.gradients{1,1}(:,1)*-1,gm.gradients{1,1}(:,2)*-1,gm.gradients{1,1}(:,3)],{surf_l,surf_r},'parcellation',labeling,'labeltext',{'Gradient 1','Gradient 2','Gradient 3'});
colormap(G_CG1.handles.figure,[.5 .5 .5;parula]);
scree_plot(gm.lambda{1});
gradient_in_euclidean([gm.gradients{1,1}(:,1)*-1,gm.gradients{1,1}(:,2)],{surf_l,surf_r},labeling);

%% step 2 gradient alignment
Gref=GradientMaps('Kernel','na','approach','dm');Gref=Gref.fit(FCzmatrix_average);
Galign=GradientMaps('Kernel','na','approach','dm','alignment','pa');

for sub=1:342
    zmatrix=[CG1_harmonized,CG2_harmonized,TG1_harmonized,TG2_harmonized];
    zmatrix_subject=zeros(400,400);
    zmatrix_subject(tril(true(400),-1))=zmatrix(:,sub);zmatrix_subject=zmatrix_subject+zmatrix_subject';
    Galign=Galign.fit(zmatrix_subject,'reference',Gref.gradients{1});
end

%% step 3 FCS
FCS_old=getFCS_400(FCzmatrix_average);FCS_TG1=getFCS_400(TG1zmatrix_average);FCS_CG1=getFCS_400(CG1zmatrix_average);
FCS_TG2=getFCS_400(TG2zmatrix_average);FCS_CG2=getFCS_400(CG2zmatrix_average);
h_FCS=plot_hemispheres([[FCS_old.L,FCS_old.R]',[FCS_CG1.L,FCS_CG1.R]',[FCS_TG1.L,FCS_TG1.R]',[FCS_CG2.L,FCS_CG2.R]',[FCS_TG2.L,FCS_TG2.R]'],{surf_l,surf_r}, ...
             'parcellation', labeling, ...
             'labeltext',{'{old_F_C_S}','{CG1_F_C_S}', '{TG1_F_C_S}','{CG2_F_C_S}', '{TG2_F_C_S}'});

FCS_CG1sub=getFCS_subject(CG1_harmonized,59);FCS_CG2sub=getFCS_subject(CG2_harmonized,59);
FCS_TG1sub=getFCS_subject(TG1_harmonized,112);FCS_TG2sub=getFCS_subject(TG2_harmonized,112);


%% step 4 split TG2 according to FCS of SMN 
parcels=readtable("D:\data\LJ\results\final\Training\0.5FD\FC\parcels.csv");%400parcels分别属于哪个网络
label=table2array(parcels(:,3));TG2FCS=[FCS_TG2sub.L,FCS_TG2sub.R];
group=find(label==2);TG2FCS_SMN=TG2FCS(:,group);TG2mean_SMN=mean(TG2FCS_SMN,2);
SMN_low=TG2_harmonized.*(TG2mean_SMN(:,1)<median(TG2mean_SMN(:,1)))';
SMN_high=TG2_harmonized.*(TG2mean_SMN(:,1)>median(TG2mean_SMN(:,1)))';
low_average=ZCX_no0mean(SMN_low,2);high_average=ZCX_no0mean(SMN_high,2);
lowzmatrix_average=zeros(400,400);highzmatrix_average=zeros(400,400);
lowzmatrix_average(tril(true(400),-1))=low_average;lowzmatrix_average=lowzmatrix_average+lowzmatrix_average';
highzmatrix_average(tril(true(400),-1))=high_average;highzmatrix_average=highzmatrix_average+highzmatrix_average';
gm=GradientMaps('Kernel','na','approach','dm','alignment','none','n_components',10,'verbose',true);gm=gm.fit(lowzmatrix_average,'sparsity',90,'tolerance',1e1);
G_low=plot_hemispheres([gm.gradients{1,1}(:,1)*-1,gm.gradients{1,1}(:,2),gm.gradients{1,1}(:,3)*-1],{surf_l,surf_r},'parcellation',labeling,'labeltext',{'Gradient 1','Gradient 2','Gradient 3'});
colormap(G_low.handles.figure,[.5 .5 .5;parula]);
scree_plot(gm.lambda{1});
gradient_in_euclidean([gm.gradients{1,1}(:,1)*-1,gm.gradients{1,1}(:,2)*-1],{surf_l,surf_r},labeling);

% %% step 3 simiulation of lesion 
% parcels=readtable("D:\data\LJ\results\final\Training\0.5FD\FC\parcels.csv");%400parcels分别属于哪个网络
% parcels=table2array(parcels(:,3));
% 
% for value = 1:7
%     binaryVector = parcels ~= value;  % 如果不等于当前值为 1，否则为 0
%     disp(['When value = ', num2str(value), ':']);
%     disp(binaryVector);  % 显示生成的 0/1 向量
%     position=find(binaryVector);
%     TG2lesion=TG2zmatrix_average.*binaryVector.*binaryVector';
%     TG2lesion=TG2lesion(any(TG2lesion, 2),:);  % 只保留至少有一个非零元素的行
%     TG2lesion=TG2lesion(:,any(TG2lesion,1));  % 只保留至少有一个非零元素的列
%     gm=GradientMaps('Kernel','na','approach','dm','alignment','none','n_components',10,'verbose',true);gm_lesion=gm.fit(TG2lesion,'sparsity',90,'tolerance',1e1);
%     gradient=gm_lesion.gradients{1,1};Gradient=zeros(400,10);Gradient(position,:)=gradient;
%     G_TG2lesion=plot_hemispheres([Gradient(:,1)*-1,Gradient(:,2)*-1,Gradient(:,3)*-1],{surf_l,surf_r},'parcellation',labeling,'labeltext',{'Gradient 1','Gradient 2','Gradient 3'});
%     colormap(G_TG2lesion.handles.figure,[.5 .5 .5;parula]);
% end

%% step 5 relationship between behav and brain
group=find(label==2);FCS_SMN=[FCS_CG1sub.fcs;FCS_CG2sub.fcs;FCS_TG1sub.fcs;FCS_TG2sub.fcs];
subFCS_SMN=FCS_SMN(:,group);submean_SMN=mean(subFCS_SMN,2);
