function [h,C] = gradient_surface(gradients,surface,parcellation,color)
% 示例：原始参考图（例如老年人template）
% [fig1,color]=gradient_surface(gm1.gradients{1}(:,1:2),{surf_lh,surf_rh},labeling,[]);
% 参数4[]用来表示没有参考颜色文件，自动生成；color为自动生成的颜色图例
% 引用该颜色文件的图（例如相比模板，CG2翻转了，我要用2D散点证明这一点）
% fig2=gradient_surface(gm2.gradients{1}(:,1:2),{surf_lh,surf_rh},labeling,color);
% 参数4 color代表引用了作为基准的图片的调色板，生成的图就可以保证下方脑区绘制和原始图一模一样，但散点位置完全变了
% 将先确定散点色再确定脑区色，变成了先确定脑区色再确定散点色

if size(gradients,2) ~= 3 && size(gradients,2) ~= 2
    error('Input matrix must be numeric with two or three columns.');
end

h.figure = figure('Color','White','Units','Normalized','Position',[0 0 .9 .9]);
h.axes_scatter = axes('Position',[.38 .6 .3 .3]);

if ~isempty(color)
    C = color;
    if size(gradients,2) == 2
        h.scatter = scatter(gradients(:,1),gradients(:,2),200,C,'Marker','.');
    else
        h.scatter = scatter3(gradients(:,1),gradients(:,2),gradients(:,3),200,C,'Marker','.');
        zlabel('Gradient 3');
    end
else
    if size(gradients,2) == 2
        cart = (gradients - mean(gradients)) ./ max(abs(gradients(:)-mean(gradients(:))));
        [th,r] = cart2pol(cart(:,2),cart(:,1));
        r = r ./ 2 + .5;
        C = [cos(.75*(th + 0*pi)), cos(.75*th - .5*pi), cos(.75*th + .5*pi)];
        C = C .* r;
        C = C ./ max(C(:));
        C(C<0) = 0;
        h.scatter = scatter(gradients(:,1),gradients(:,2),200,C,'Marker','.');
    else
        C = (gradients - min(gradients)) ./ max(gradients - min(gradients));
        h.scatter = scatter3(gradients(:,1),gradients(:,2),gradients(:,3),200,C,'Marker','.');
        zlabel('Gradient 3');
    end
end

set(h.axes_scatter                              , ...
    'DataAspectRatio'       , [1 1 1]           , ...
    'PlotBoxAspectRatio'    , [1 1 1]           , ...
    'FontName'              , 'DroidSans'       , ...
    'FontSize'              , 14                );
xlabel('Gradient 1'); ylabel('Gradient 2');

if nargin > 1
    % If parcellated data is provided, bring it to the full mesh.
    if exist('parcellation','var')
        D_full = parcel2full((1:size(C,1))',parcellation);
    end
    
    % Deal with nans.
    D_full(all(isnan(D_full),2)) = .7; 
    
    % Make sure surfaces are in SurfStat format.
    S = convert_surface(surface);
    if ~iscell(S); S = {S}; end
    if numel(S) > 3; error('More than two surfaces are not accepted.'); end
    
    % Check correct dimensions of C_full and surfaces
    N1 = sum(cellfun(@(x)size(x.coord,2),S));
    N2 = size(D_full,1);
    if N1 ~= N2; error('Number of vertices on the surface and number of data points do not match.'); end
    
    % Split C_full over surfaces.
    D{1} = D_full(1:size(S{1}.coord,2));
    if numel(S) == 2
        D{2} = D_full(size(S{1}.coord,2)+1:end);
    end
    
    % Plot the surface(s).
    for ii = 1:numel(S)*2
        idx = ceil(ii/2);
        h.axes(ii) = axes('Position',[-.1+ii*.2 .2 .3 .3]);
        h.trisurf(ii) = trisurf(S{idx}.tri, ...
            S{idx}.coord(1,:), ...
            S{idx}.coord(2,:), ...
            S{idx}.coord(3,:), ...
            D{idx}, ...
            'EdgeColor', 'None');
        material dull; lighting phong;
        colormap(h.axes(ii),[.7 .7 .7; C])
    end
    set(h.axes                              , ...
        'Visible'           , 'off'         , ...
        'DataAspectRatio'   , [1 1 1]       , ...
        'PlotBoxAspectRatio', [1 1 1]       , ...
        'CLim'              , [0,size(C,1)] );
    % Add axes embelishments
    set(h.axes(1),'View',[-90 0]);
    set(h.axes(2),'View',[90 0]);
    if numel(h.axes) == 4
        set(h.axes(3),'View',[-90 0]);
        set(h.axes(4),'View',[90 0]);
    end

    % Add a camlight.
    for ii = 1:numel(h.axes)
        axes(h.axes(ii));
        h.camlight(ii) = camlight();
    end
end
end