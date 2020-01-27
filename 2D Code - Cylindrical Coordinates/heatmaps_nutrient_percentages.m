% plots nutrient percentages as a heatmap


%% algae
A = [0.1096 0.1288 0.0075 0.0006; ...
     0.1107 0.1342 0.0115 0.0010; ...
     0.1180 0.1377 0.0162 0.0015; ...
     0.1631 0.1331 0.0165 0.0015];

 dx = [1 10 100 1000];
 dz = [1 10 100 1000];
 figure(1)
 XTick = [ 1 2 3 4];
 YTick = [ 1 2 3 4];
  imagesc(XTick, YTick,A);
%   gca.XTickMode = 'manual';
%     gca_handle.YTickMode = 'manual';
% gca_handle.XTicklabels = dz;
% gca_handle.yticklabels = dx;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel',dz);
set(gca, 'YTick', [1 2 3 4], 'YTickLabel', dx) ;
  colormap(jet);
  colorbar
 title('Algal nutrient content');
 xlabel('dz');
 ylabel('dx');
  caxis([0 1]);
 
 %% Dissolved nutrients
 
 R = [0.8622 0.3632 0.1083 0.0441; ...
      0.8610 0.3772 0.1140 0.0651; ...
     0.8499 0.3855 0.1877 0.2103; ...
     0.7867 0.3762 0.2154 0.2809];

 dx = [1 10 100 1000];
 dz = [1 10 100 1000];
 figure(2)
 XTick = [ 1 2 3 4];
 YTick = [ 1 2 3 4];
  imagesc(XTick, YTick,R);
%   gca.XTickMode = 'manual';
%     gca_handle.YTickMode = 'manual';
% gca_handle.XTicklabels = dz;
% gca_handle.yticklabels = dx;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel',dz);
set(gca, 'YTick', [1 2 3 4], 'YTickLabel', dx) ;
  colormap(jet);
  colorbar
 title('Dissolved nutrient content');
 xlabel('dz');
 ylabel('dx');
  caxis([0 1]);
 
  %% sediment nutrients
 
 Rs = [0.0282 0.5080 0.8842 0.9553; ...
       0.0283 0.4886 0.8745 0.9339; ...
       0.0321 0.4768 0.7961 0.7883; ...
       0.0502 0.4907 0.7680 0.7176];

 dx = [1 10 100 1000];
 dz = [1 10 100 1000];
 figure(3)
 XTick = [ 1 2 3 4];
 YTick = [ 1 2 3 4];
  imagesc(XTick, YTick,Rs);
%   gca.XTickMode = 'manual';
%     gca_handle.YTickMode = 'manual';
% gca_handle.XTicklabels = dz;
% gca_handle.yticklabels = dx;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel',dz);
set(gca, 'YTick', [1 2 3 4], 'YTickLabel', dx) ;
  colormap(jet);
  colorbar
 title('Sediment nutrient content');
 xlabel('dz');
 ylabel('dx');
 caxis([0 1]);
