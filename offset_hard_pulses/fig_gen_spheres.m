function fig_gen_spheres(factor,angle_deg,angle_view)
if nargin==0
    factor=1;
    angle_deg=90;
end

pul_dur=10e-6;
angle_pulse=90/180*pi;%deg
ampli_hz=(angle_pulse/pul_dur)/(2*pi);
%disp(['pulse amplitude : ' num2str(ampli_hz) ' Hz'])



cur_f=figure(1);clf;
inca=1;
if nargin<3
    view([51 31]);
else
    view(angle_view);
end
list_angl=(0:inca:360)*pi/180;
hold on
plot3(0*cos(list_angl),1*cos(list_angl),1*sin(list_angl),'k-','color',[1 1 1]*0.5)
plot3(1*cos(list_angl),0*cos(list_angl),1*sin(list_angl),'k-','color',[1 1 1]*0.5)
plot3(1*cos(list_angl),1*sin(list_angl),0*cos(list_angl),'k-','color',[1 1 1]*0.5)
plot3([0 0],[0 0],[-1 1],'k-','color',[1 1 1]*0.5)
plot3([0 0],[-1 1],[0 0],'k-','color',[1 1 1]*0.5)
plot3([-1 1],[0 0],[0 0],'k-','color',[1 1 1]*0.5)
axis('equal')
offsset_first_null=sqrt(15)/(4*pul_dur);
%disp(['offset_first_null : ' num2str(offsset_first_null) ' Hz ' ])
count_main=0;
start_pt=[0 -1 0];
start_pt_p=start_pt;
inc_store=1;
for loop_offset=factor*ampli_hz% not a loop anymore
    if count_main==0
        disp_on=1;
        count_main=count_main+1;
        
    else
        disp_on=0;
        
        count_main=count_main+1;
        if count_main==(5)
            count_main=0;
        end
    end
    nu_eff=sqrt(loop_offset*loop_offset+ampli_hz*ampli_hz);
    
    %plot field vector
    tilt_angle=atan((ampli_hz/loop_offset));
    if tilt_angle<0, tilt_angle=tilt_angle+pi;end
    % if disp_on
    % disp(['Offset : ' num2str(loop_offset) ' Hz w_eff=' num2str(nu_eff) ' Hz'])
    how_much_further=2;
    if size(factor,2)==1
        plot3(how_much_further*[0 sin(tilt_angle)],[0 0],how_much_further*[0 cos(tilt_angle)],'k--')
        text(how_much_further*[sin(tilt_angle)],[ 0],how_much_further*[ cos(tilt_angle)],'Beff')
        
    end
    %   plot3([sin(tilt_angle)],[ 0],[ cos(tilt_angle)],'ko','MarkerSize',5)
    
    % end
    pos_mag=[0 0 1];
    pos_mag_p=pos_mag;
    field=[0 0 1];
    increment_tilt=pi/100000;
    inc=0;
    inc_sto=1;
    for til_tim=0:increment_tilt:angle_deg/180*pi
        di=cross([sin((tilt_angle)) 0 cos((tilt_angle)) ],pos_mag);
        di=di/norm(di);
        
        pos_mag=pos_mag+di*increment_tilt;
        if inc==0
            %  if disp_on,
            %    plot3([pos_mag_p(1,1) pos_mag(1,1)],[pos_mag_p(1,2) pos_mag(1,2)],[pos_mag_p(1,3) pos_mag(1,3)],'r-')
            stor_tr(inc_sto,:)=pos_mag;inc_sto=inc_sto+1;
            %  end
            pos_mag_p=pos_mag;
            inc=inc+1;
            
        else
            inc=inc+1;
            if inc==1000
                inc=0;
            end
        end
        
    end
    if disp_on
        plot3(stor_tr(:,1),stor_tr(:,2),stor_tr(:,3),'k:','linewidth',1.25)
        
    end
    %  if disp_on,
    % plot3([pos_mag_p(1,1) ],[pos_mag_p(1,2) ],[ pos_mag(1,3)],'ko')
    % plot3([pos_mag_p(1,1) ],[pos_mag_p(1,2) ],[ pos_mag(1,3)],'k.','MarkerSize',12)
    plot3([0 pos_mag_p(1,1) ],[0 pos_mag_p(1,2) ],[0 pos_mag(1,3)])
    arrow([0 0 0 ],[pos_mag_p(1,1)  pos_mag_p(1,2)  pos_mag(1,3)],'linewidth',1.5)
    %  end
    drawnow
    start_pt_p=start_pt;
    start_pt=[pos_mag_p(1,1) pos_mag_p(1,2)  pos_mag(1,3)];
    store_traj(inc_store,:)=start_pt;inc_store=inc_store+1;
    % plot3([start_pt_p(1,1) start_pt(1,1)],[start_pt_p(1,2) start_pt(1,2)],[start_pt_p(1,3) start_pt(1,3)],'g-')
    
end
mi=min(min(factor));
ma=max(max(factor));
if size(factor,2)==1
    txt_ti=[' ' num2str(mi) ' x B1'];
else
    if mi==-ma
        txt_ti=[' +/- ' num2str(ma) ' x B1'];
        
    else
        txt_ti=[' (' num2str(mi) ':' num2str(ma) ') x B1'];
    end
end
text(0.4,0.4,-0.9,[txt_ti '  ' num2str(angle_deg) ' deg.'])
axis([ -1     1    -1     1    -1     1])
%plot3(store_traj(:,1),store_traj(:,2),store_traj(:,3),'k-','linewidth',1.5)
%plot3(store_traj(:,1),store_traj(:,2),0*store_traj(:,3),'k-','linewidth',1)
axis off

end
