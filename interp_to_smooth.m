function [traj_out, t2]=interp_to_smooth(traj_in,t)
%this function is used to smooth plot in angluar coordinates
max_angle=5;
if size(traj_in,1)>1
    pointer=1;
    [theta,rho,z] = cart2sph(traj_in(:,1),traj_in(:,2),traj_in(:,3));
    %    txt='';
    for loo=1:size(z,1)-1
        % check change in angle
        delta_theta=abs(theta(loo)-theta(loo+1))/pi*180;
        if abs(delta_theta)>180
            delta_theta=abs(360-abs(delta_theta));
        end
        delta_rho=abs(rho(loo)-rho(loo+1))/pi*180;
        if abs(delta_rho)>180
            delta_rho=abs(360-abs(delta_rho));
        end
        delta_largest=max([delta_theta delta_rho]);
        if delta_largest>max_angle
            size_int=round(delta_largest/max_angle+0.5);
            coef2=[0:1/size_int:1-1/size_int];
            coef1=1-coef2;
            theta_other=theta(loo+1);
            if (theta_other-theta(loo))>pi
                theta_other=theta_other-2*pi;
                % disp('fix 1')
            end
            if (theta_other-theta(loo))<-pi
                theta_other=theta_other+2*pi;
                %    disp('fix 2')
                
            end
            rho_other=rho(loo+1);
            if (rho_other-rho(loo))>pi
                rho_other=rho_other-2*pi;         %       disp('fix 3')
                
            end
            if (rho_other-rho(loo))<-pi
                rho_other=rho_other+2*pi;              %  disp('fix 4')
                
            end
            theta_int=theta(loo)*coef1+theta_other*coef2;
            rho_int=rho(loo)*coef1+rho_other*coef2;
            z_int=z(loo)*coef1+z(loo+1)*coef2;
            t_int=t(loo)*coef1+t(loo+1)*coef2;
            %    txt=[txt '1'];
        else
            
            size_int=1;
            theta_int=theta(loo);
            rho_int=rho(loo);
            z_int=z(loo);
            t_int=t(loo);
            
            %    txt=[txt '_'];
            
        end
        theta2(1,pointer:pointer+size_int-1)=theta_int;
        rho2(1,pointer:pointer+size_int-1)=rho_int;
        z2(1,pointer:pointer+size_int-1)=z_int;
        t2(1,pointer:pointer+size_int-1)=t_int;
        pointer=pointer+size_int;
    end
    %   disp([txt ' ' num2str(size(z,1)) ' ' num2str(size(z2,2)) ])
    theta2(pointer)=theta(size(z,1));
    rho2(pointer)=rho(size(z,1));
    z2(pointer)=z(size(z,1));
    t2(pointer)=t(size(z,1));
    
    [traj_out(:,1),traj_out(:,2),traj_out(:,3)] = sph2cart(theta2',rho2',z2');
    
else
    traj_out=traj_in;
    t2=t;
end

end