
clear all;
close all;
clc;
az =0;
el = 90;
figure(1)
ylabel('Water depth');
xlabel('Distance from lake center');
bottom = 0;
top = 100;

for y = 1:4
    for x=1:4
        % Loading data
        load("dx_"+ num2str(10^(x-1))+ "_dy_" + num2str(10^(y-1)) + "_res_50_50");
        
        % We plot the nutrient concentration of the alge, so that the
        % concentration values are directly comparable to the dissolved
        % algae.
        A = A.*p.q;
        
        %% Algae
        %subplot(6,3,6*(x-1)+y)
        if(y==1)
        subaxis(8, 4, 8*(x-1)+y, 'sh', 0.000, 'sv', 0.000, ...
            'PL', 0.00,  'PT', 0.015, 'PB', 0.001, ...
            'MT', 0.000, 'ML', 0.03, 'MR', 0.00);
        else 
            subaxis(8, 4, 8*(x-1)+y, 'sh', 0.000, 'sv', 0.000, ...
            'PL', 0.000,  'PT', 0.015, 'PB', 0.001, ...
            'MT', 0.000, 'ML', 0.00, 'MR', 0.00);
        end

        surf(p.X_vol,p.Y_vol,A,  'edgecolor','none', 'FaceAlpha', 1)
        view(az,el);
        set(gca, 'YDir','reverse', 'XDir', 'reverse')
        %set(gca, 'YDir','reverse')
        
        % interval for common color bar
        bottom = min(min(min(A)),min(min(Rd)));
        top  = max(max(max(A)),max(max(Rd)));
        
                % Common color bar
                colormap(jet);
%         caxis manual
% caxis([bottom top]);
     
        
        xticks(10:10:50);
        yticks(10:10:50);
        xlim([0 50])
        ylim([0 50])
        
        %%create a second axes.
        ax1 = gca; % the first axes
        ax1.FontSize = 7;
        ax1.FontWeight = 'bold';
        ax2 = axes('Position',ax1.Position,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none',...
            'Ylim',ax1.YLim,...
            'XLim',ax1.XLim,...
            'TickLength',[0 0],...
            'YTick', [ax1.YLim(1):5:ax1.YLim(2)], ...
            'XTick', [ax1.XLim(1):5:ax1.XLim(2)],  ...
            'YTickLabel', [],  ...
            'XTickLabel', []  );
        ax2.FontSize = 7;
        linkaxes([ax1 ax2],'xy')
        grid on
        temp = colorbar;
        set(temp,'YTick',[])
        set(gcf,'CurrentAxes',ax1); % return to ax1
        cbar = colorbar;
        cbar.FontWeight = 'bold';
        
        if y > 1
            set(gca,'YTickLabel',[]);
        end
        set(gca,'XTickLabel',[]);
        %caxis manual
        %caxis([bottom top]);
        
 
        %% Dissolved Nutrients
   if(y==1)
            subaxis(8, 4, 8*(x-1)+4+y, 'sh', 0.000, 'sv', 0.000, ...
                'PL', 0.00, 'PT', 0.001, 'PB', 0.015,...
                'MT', 0.00, 'ML', 0.03, 'MR', 0.00);
   else
                 subaxis(8, 4, 8*(x-1)+4+y, 'sh', 0.000, 'sv', 0.000, ...
                'PL', 0.000, 'PT', 0.001, 'PB', 0.015,...
                'MT', 0.000, 'ML', 0.00, 'MR', 0.00);
   end
        
        surf(p.X_vol,p.Y_vol,Rd,  'edgecolor','none', 'FaceAlpha', 1)
        alpha(1)
        view(az,el);
        set(gca, 'YDir','reverse', 'XDir', 'reverse')
        %set(gca, 'YDir','reverse')
        
        % Common color bar
        colormap(jet);
%         caxis manual
% caxis([bottom top]);
        
        xticks(10:10:50);
        yticks(10:10:50);
        xlim([0 50])
        ylim([0 50])
        
        %%create a second axes.
        ax1 = gca; % the first axes
        ax1.FontSize = 7;
        ax1.FontWeight = 'bold';
        ax2 = axes('Position',ax1.Position,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none',...
            'Ylim',ax1.YLim,...
            'XLim',ax1.XLim,...
            'TickLength',[0 0],...
            'YTick', [ax1.YLim(1):5:ax1.YLim(2)], ...
            'XTick', [ax1.XLim(1):5:ax1.XLim(2)],  ...
            'YTickLabel', [],  ...
            'XTickLabel', []  );
        
        
        linkaxes([ax1 ax2],'xy')
        grid on
        temp = colorbar;
        set(temp,'YTick',[])
        set(gcf,'CurrentAxes',ax1); % return to ax1
        cbar = colorbar;
        cbar.FontWeight = 'bold';
        % Removed axis unless on the bottom/left edge
        if y > 1
            set(gca,'YTickLabel',[]);
        end
        if x < 4
            set(gca,'XTickLabel',[]);
        end
        %         p.dx
        %         p.dy
    end
end

