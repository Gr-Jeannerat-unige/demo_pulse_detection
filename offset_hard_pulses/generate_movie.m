clear all
%for main_ratio=[ 0.01 1 0.5 0.2 0.1],
for mainlooop=-3:6
    % for mainlooop=-3
    list_factor=mainlooop+0.01;
    list_an=[90 ];
    if (mainlooop==-1)
        list_factor=[0 1 2 3 4 5 6];
        
    end
    if mainlooop==-2 %not working does not know why
        list_factor=[-0.50:0.1:0.5];
        list_an=[90 180 270 360];
    end
    if mainlooop==-3 %not working does not know why
        list_factor=[-6:0.2:6];
        list_an=[90 ];
    end
    for main_ratio=list_an
        
        figure(mainlooop+4)
        Z = peaks;
        surf(Z)
        axis tight manual
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        
        switch mainlooop
            case -3
                writerObj = VideoWriter(['mov_-6_0.2_6_' num2str(main_ratio) 'deg.avi']);
            case -2
                writerObj = VideoWriter(['mov_-0.5_0.1_0.5_' num2str(main_ratio) 'deg.avi']);
                
            case -1
                writerObj = VideoWriter(['mov_0_1_6_' num2str(main_ratio) 'deg.avi']);
            otherwise
                writerObj = VideoWriter(['mov_' num2str(mainlooop) '_' num2str(main_ratio) 'deg.avi']);
        end
        
        open(writerObj);
        
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
        
        list_of_values= 0.00:1:main_ratio+1;
        loops = 360;
        loops = 1;
        
        F(size(list_of_values,2)+1) = struct('cdata',[],'colormap',[]);
        counter=1;
        for loopj = list_of_values
            % X = sin(loopj*pi/10)*Z;
            %surf(X,Z)
            clf
            
            anglet= loopj-1;
            % phase_a=cos(anglet/360*2*pi)+j*sin(anglet/360*2*pi);%as complex number
            %  phase_a=phase_a/(abs(phase_a));
            %  [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
            %  [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
            
            if main_ratio==270
                fig_gen_spheres(list_factor,loopj,[51 31]);%[152 19]
            else
                fig_gen_spheres(list_factor,loopj,[51 31])%[152 19]
                
            end
            draw_unit_circle([0;0;0]);
            set(gcf,'color','w');
            
            drawnow;
            tmp_frame = getframe;
            if counter==1
                si=size(tmp_frame.cdata);
            else
                tmp_frame.cdata=tmp_frame.cdata(1:si(1,1),1:si(1,2),1:si(1,3));
            end
            
            
            F(counter)=tmp_frame;
            writeVideo(writerObj,F(counter));
            counter=counter+1;
            
            
        end
        
        close(writerObj);
        %movie(F,2)
    end
    %     figure(111)
    %     clf
    %     plot(store_dis_in_hz(:,1),store_erro_in_deg(:,1),'b-');
    %     hold on
    %     plot(store_dis_in_hz(:,2),store_erro_in_deg(:,2),'r-');
    %     print('-depsc','-tiff','-r600',[ 'Phase_error_nearby_small_signals' num2str(main_ratio)  '.eps']);%here
    
    
end
