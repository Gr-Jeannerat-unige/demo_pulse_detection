clear all
%for main_ratio=[ 0.01 1 0.5 0.2 0.1],
%for mainlooop=0:6
for mainlooop=0:6
    % for mainlooop=-3
    clear stor_tr
    clear stor_tr_crude
    clear stor_t_crude
clear stor_t
    main_ratio=90;
    
    figure(1)
    %  Z = peaks;
    %  surf(Z)
    %  axis tight manual
    %  ax = gca;
    %  ax.NextPlot = 'replaceChildren';
    
    
    writerObj = VideoWriter(['mov_' num2str(mainlooop) '.avi']);
    %  writerObj = VideoWriter(['mov_' num2str(mainlooop) '.mp4']);
    writerObj.Quality = 75;%75 is default
    
    
    open(writerObj);
    
    axis tight
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    
    list_of_values= 0.00:1:main_ratio+1;
    loops = 360;
    loops = 1;
    count_sto=1;
    F(size(list_of_values,2)+1) = struct('cdata',[],'colormap',[]);
    counter_frames=1;
    
    pul_dur=10e-6;
    angle_pulse=90/180*pi;%90 deg
    ampli_hz=(angle_pulse/pul_dur)/(2*pi);
    disp(['pulse amplitude : ' num2str(ampli_hz) ' Hz'])
    loop_offset=0+000*1000;
    start_pt=[0 -1 0];
    if (mainlooop==0)
        tmax=0.3;
    else
        tmax=0.1;
        
    end
    pos_mag=[0 0 1];
    t_step=0.05;
    firstr1=0;
    tau=10;%(defaul no second pulse
    if mainlooop==4 tau=0.01;end
    if mainlooop==5 tau=0.03;end
    if mainlooop==6 tau=0.05;end
    
    for t = 0:0.000001:tmax
        % for t = 0:0.000001:0.05
        % X = sin(loopj*pi/10)*Z;
        %surf(X,Z)
        larmo=10.123122;
        
        increment_tilt=pi/100000;
        inc=0;
        
        %  tilt_angle=atan((ampli_hz/loop_offset));
        
        %  if tilt_angle<0, tilt_angle=tilt_angle+pi;end
        rf=0.0000001;
        % if mod(round(t*20),6)==0
        R2=0;
        R1=0;
        if (mainlooop==0)
            if t<(1/1000)
                rf=10;larmo=0000.00001;
            end
            if t>0.1
                R2=0.5;
            end
            if (t>0.18)
                R1=1;
            end
            
        else
            tim=400;
            rf_ref=20;
fa=1;
            tim=400*fa;
            rf_ref=20*fa;
            if t<(1/tim)
                if mainlooop==2
                    
                    rf=rf_ref/3;larmo=2;
                else
                    if mainlooop>=3
                        
                        rf=rf_ref*2;larmo=0.00000001;
                    else
                        rf=rf_ref;larmo=2;
                        
                    end
                end
            else
                % if t>0.1
                R2=2;
                % end
                % if t>0.18
                R1=1;
                %  end
                
            end
            if (t>=tau) && (t<=(tau+1/400))
                rf=20;larmo=0.01;
                firstr1=0;
            end
            
        end
        di=cross([rf 0 larmo ],pos_mag);
        %  di=di/norm(pos_mag);
        
        pos_mag=pos_mag+di*increment_tilt;
        pos_mag(1,1)=pos_mag(1,1)*exp(-R2*increment_tilt);
        pos_mag(1,2)=pos_mag(1,2)*exp(-R2*increment_tilt);
        if R1~=0
            if firstr1==0
                valref=pos_mag(1,3);diff=1-valref;
                firstr1=1;
            end
            diff=diff*((exp(-R1*increment_tilt)));
            pos_mag(1,3)=1-diff;
        end
        
        nu_eff=sqrt(loop_offset*loop_offset+ampli_hz*ampli_hz);
        
        %
        %      how_much_further=2;
        %     if size(factor,2)==1
        %         plot3(how_much_further*[0 sin(tilt_angle)],[0 0],how_much_further*[0 cos(tilt_angle)],'k--')
        %         text(how_much_further*[sin(tilt_angle)],[ 0],how_much_further*[ cos(tilt_angle)],'Beff')
        %
        %     end
        %     %%%%%%%%%%%%%%%%%
        %
        %
        
        
        
        %    figure('Units','inches', 'PaperPositionMode','auto', 'Position',[0 0 4 4]);
        
        
        
        
        
        %work on figure
        if mod(round(t*1e6),500)==0
            where_low=-3;
            disp(['time: ' num2str(t) '/' max(num2str(tmax)) ' ' num2str(rf)] )
            fig_gen_spheres(pos_mag)
            % store last point for full trajectory plot
            stor_tr_crude(count_sto,:)=pos_mag(1,:);
            stor_t_crude(count_sto,:)=t;
            count_sto=count_sto+1;
            % interpolation of
            [stor_tr, stor_t]=interp_to_smooth(stor_tr_crude,stor_t_crude);
            plot3(stor_tr(:,1),stor_tr(:,2),stor_tr(:,3),'k-','linewidth',1.25)
            %   plot3(stor_tr_crude(:,1),stor_tr_crude(:,2),stor_tr_crude(:,3),'o','linewidth',1.25)
            
            plot3( [0.5 0.5],[-1 1],[where_low where_low],'k-')
            plot3(-[0.5 0.5],[-1 1],[where_low where_low],'k-')
            plot3(-[0.0 0.0],[-1 1],[where_low where_low],'k:')
            plot3(+[1 1],[-1 1],[where_low where_low],'k:')
            plot3([-1 1],[-1 -1],[where_low where_low],'k-')
            plot3(-[1 1],[-1 1],[where_low where_low]+0.5,'k-')
            plot3(-[1 1],[-1 1],[where_low where_low]+1,'k:')
            plot3(-[1 1],[-1 1],[where_low where_low]+0,'k:')
            plot3(-[1 1],[-1 -1],[0 1 ]+where_low,'k-')
            
            %   list_t=2*(t/tmax)*(1:size(stor_tr,1))/size(stor_tr,1)-1;
            list_t=2*(t/tmax)*(stor_t)/max(abs(stor_t))-1;
            plot3( 0.5+0.5*stor_tr(:,1),list_t,stor_tr(:,3)*0+where_low,'g-','linewidth',1.25)
            
            
            plot3(-0.5+0.5*stor_tr(:,2),list_t,stor_tr(:,3)*0+where_low,'b-','linewidth',1.25)
            plot3(stor_tr(:,3)*0-1,list_t,0.5+0.5*stor_tr(:,3)+where_low,'r-','linewidth',1.25)
            
            % phase_a=cos(anglet/360*2*pi)+j*sin(anglet/360*2*pi);%as complex number
            %  phase_a=phase_a/(abs(phase_a));
            %  [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
            %  [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
            
            
            axis([ -1     1    -1     1    where_low     1])
            %plot3(store_traj(:,1),store_traj(:,2),store_traj(:,3),'k-','linewidth',1.5)
            %plot3(store_traj(:,1),store_traj(:,2),0*store_traj(:,3),'k-','linewidth',1)
            axis off
            
            set(gcf,'color','w');
            drawnow;
            tmp_frame = getframe;
            if counter_frames==1
                si=size(tmp_frame.cdata);
            else
                tmp_frame.cdata=tmp_frame.cdata(1:si(1,1),1:si(1,2),1:si(1,3));
            end
            
            
            F(counter_frames)=tmp_frame;
            writeVideo(writerObj,F(counter_frames));
            counter_frames=counter_frames+1;
        end
        
    end
    
    close(writerObj);
    %movie(F,2)
    
    %     figure(111)
    %     clf
    %     plot(store_dis_in_hz(:,1),store_erro_in_deg(:,1),'b-');
    %     hold on
    %     plot(store_dis_in_hz(:,2),store_erro_in_deg(:,2),'r-');
    %     print('-depsc','-tiff','-r600',[ 'Phase_error_nearby_small_signals' num2str(main_ratio)  '.eps']);%here
    
end

