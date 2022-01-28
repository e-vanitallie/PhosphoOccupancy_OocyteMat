function [fitting_plot,p_normal_V,data_p] = plot_SVD_out(fitting_plot,numP,dim_sp,sp,subplot_title,...
    FP_color_mat,NumConditions,InputMat,Vmat,Smat,normal_v,occv,impossible_flag)

% ploting code mainly written by Marc Presler
% slightly adjusted by E. Van Itallie 2020 
data_p = [];
planeColor = [0.5 0.5 0.5]; %colors(6,:);
fitting_plot.Visible = 'off';
%figure(fitting_plot)
subplot(dim_sp(1),dim_sp(2),sp)
hold on 

% instead of hard-coding this separately twice, extract and
% re-use subplot_title for 2 dims and 3 dim plots

%In 2 dimensions
if numP == 2

    % PC1_vector = quiver(0,0,Vmatrix(1,1),Vmatrix(2,1),Smatrix(1,1),'linewidth',1.5,'color',  planeColor);
    % PC1_vector_minus = quiver(0,0,-Vmatrix(1,1),-Vmatrix(2,1),Smatrix(1,1),'linewidth',1.5,'color',  planeColor);
    % PC1_vector = quiver(0,0,Vmat(1,1),Vmat(2,1),0.75*Smat(1,1),'linewidth',1,'color',  planeColor);
    % PC1_vector_minus = quiver(0,0,-Vmat(1,1),-Vmat(2,1),0.75*Smat(1,1),'linewidth',1,'color',  planeColor);
    
    PC1_vector = quiver(0,0,Vmat(1,1),Vmat(2,1),Smat(1,1),'linewidth',1,'color',  planeColor);
    PC1_vector_minus = quiver(0,0,-Vmat(1,1),-Vmat(2,1),Smat(1,1),'linewidth',1,'color',  planeColor);

    for v = 1:NumConditions
        data_p(v) = plot(InputMat(v,1), InputMat(v,2),...
            'o','MarkerEdgeColor','k','MarkerSize',10,'MarkerFaceColor',FP_color_mat(v,:));
    end
    
    p_normal_V = quiver(0,0,normal_v(1),normal_v(2),1,'linewidth',1.5,'color','k'); 
    
    set(gca,'fontsize',12)
    
    if (impossible_flag)
        title({subplot_title,...
            'Occ = INFEASIBLE'},...
            'fontsize',12,'color','k')
    else
        title({subplot_title,...
            ['Occ = [' num2str(occv(1)) '%, ' num2str(occv(2)) '%]']},...
            'fontsize',12,'color','k')
    end
    
    xlabel('Parent Form','fontsize',12)
    ylabel('Phospho Form','fontsize',12)
    
    axis equal
    set(gcf,'color','w')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    %In 3 dimensions
else
    
    hold on
    
    for v = 1:NumConditions
        data_p1(v) =  plot3(InputMat(v,1),InputMat(v,2),InputMat(v,3),...
            'o','MarkerEdgeColor','k','MarkerSize',10,'MarkerFaceColor',FP_color_mat(v,:));
    end
%     
%     plot3(normal_v(1)*range,normal_v(2)*range,normal_v(3)*range,'--','color','k','linewidth',1)
%     
%     S_scale = 0.75*Smat(1,1);
%     
    PC1_vector = quiver3(0,0,0,Vmat(1,1)*Smat(1,1),Vmat(2,1)*Smat(1,1),Vmat(3,1)*Smat(1,1),1,'linewidth',1.5,'color',  planeColor);
     PC1_vector_minus = quiver3(0,0,0,-Vmat(1,1)*Smat(1,1),-Vmat(2,1)*Smat(1,1),-Vmat(3,1)*Smat(1,1),1,'linewidth',1.5,'color',  planeColor);
%     
  PC2_vector = quiver3(0,0,0,Vmat(1,2)*Smat(2,2),Vmat(2,2)*Smat(2,2),Vmat(3,2)*Smat(2,2),1,'linewidth',1.5,'color',planeColor);
     PC2_vector_minus = quiver3(0,0,0,-Vmat(1,2)*Smat(2,2),-Vmat(2,2)*Smat(2,2),-Vmat(3,2)*Smat(2,2),1,'linewidth',1.5,'color',planeColor);
%     
     PC1_values = [PC1_vector.UData, PC1_vector.VData, PC1_vector.WData];
     PC2_values = [PC2_vector.UData, PC2_vector.VData, PC2_vector.WData];
%     
      X_grid = [PC1_values(1), PC2_values(1),;...
          -PC2_values(1), -PC1_values(1)];
      Y_grid = [PC1_values(2), PC2_values(2),;...
          -PC2_values(2), -PC1_values(2)];
     
          % ALTERNATIVE 3D plot option 
    t = 0.1; % transparency 
    quiv_D3 = sqrt(0.5*(Smat(1,1)+Smat(2,2)));
    
    plotcube(quiv_D3*[1 -1 1],[0 0 0],t,[0.5 0 0])
    plotcube(quiv_D3*[-1 1 1],[0 0 0],t,[0.5 0 0])
    
    plotcube(quiv_D3*[-1 1 -1],[0 0 0],t,[0 0 0.5])
    plotcube(quiv_D3*[1 -1 -1],[0 0 0],t,[0 0 0.5])
    
   % [x,y] = meshgrid(-1:0.1:1); % Generate x and y data
    
    z = -1/(normal_v(3))*(normal_v(1)*X_grid + normal_v(2)*Y_grid); % Solve for z data
    data_p2 = surf(X_grid,Y_grid,z,'FaceColor',planeColor,'FaceAlpha',0.4); %Plot the surface
    p_normal_V = quiver3(0,0,0,quiv_D3*normal_v(1),quiv_D3*normal_v(2),quiv_D3*normal_v(3),1,'-','color','k','linewidth',2);
    
    data_p = [data_p2 data_p1];
    
    h = get(gca,'DataAspectRatio');
    set(gca,'DataAspectRatio',[1 1 1])
    set(gcf,'color','w')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     legend([s,p3_norm],{'2D Plane','Normal Vector'},'Location','best');
    
    set(gca,'fontsize',10)
    ylabel('Phospho Form 1','fontsize',10)
    zlabel('Phospho Form 2','fontsize',10)
    xlabel('Parent Form','fontsize',10)
    
 
    if (impossible_flag)
        title({subplot_title,...
            'OV = INFEASIBLE'},'color','k')
    else
        title({subplot_title,...
            [' OV = [' num2str(occv(1)) '%, ' num2str(occv(2)) '%, ' num2str(occv(3)) '%]']},'color','k')
        
    end
    
    view(225,60) % is the default 3-D view.
    


%     
%     Z_grid = (-(normal_v(1).*X_grid + normal_v(2).*Y_grid))./(normal_v(3));
%     mesh(X_grid,Y_grid,Z_grid,...
%         'EdgeAlpha',0,...
%         'FaceColor',planeColor,...
%         'FaceAlpha',0.3)
%     hold off
%     grid on
%     
%     set(gca,'fontsize',8)
%     title(subplot_title,'fontsize',8)
%     zlabel('Parent Form','fontsize',8)
%     xlabel('Phospho Form 1','fontsize',8)
%     ylabel('Phospho Form 2','fontsize',8)
%     
%     h = get(gca,'DataAspectRatio');
%     set(gca,'DataAspectRatio',[1 1 1])
%     set(gcf,'color','w')
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     
%     upperaxis = max(max([X_grid; Y_grid;Z_grid]));
%     loweraxis = -upperaxis;
%     xlim([loweraxis, upperaxis])
%     ylim([loweraxis, upperaxis])
%     zlim([loweraxis, upperaxis])
%     campos([15 -4.5 3])
    

    
end