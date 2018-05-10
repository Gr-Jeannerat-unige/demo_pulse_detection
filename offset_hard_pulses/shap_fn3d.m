function [dist_in_hz erro_in_deg]=shap_fn3d(phi_an_orig,mainlooop,main_ratio)
if nargin<3,
    main_ratio=1;
end
if nargin<2,
    mainlooop=0;
end
if nargin<1,
    phi_an_orig=0;
end
phase_in=0;
dist_in_hz=0;
erro_in_deg=0;
add_text='';
sw=20;
dw=1/sw;

%t=0:dw:1000;
%t=0:dw:10000-dw;
th=0:dw/2:10000-dw/2;
t=th(1:2:end);
lb=1/2;%1/2
nu=[0 -8 -5 3 9.5];
%nu=[2.01321321 ];
amp=[1 -0.1 -0.1 -0.1 0.5];
height=[0.99 0.9 0.5 ];
height=[0.99 0.9 ];
heightfirst=[ 0.5 0.1 0.02  ];
heightfirstd=[ 0.1 0.02 ];

    %%sim td complex points
    % fid_sim= (main_ratio)*exp(j*2*pi*((t )*nu)).*exp(-t*2*pi*lb);
    %%seq td*2 real points
    fid_full=(0)*th;
    for loo_spi=1:size(nu,2),
        fid_full=fid_full+amp(1,loo_spi)*exp(j*2*pi*((th)*nu(1,loo_spi))).*exp(-th*2*pi*lb);
    end
    fid_sim=fid_full(1:2:end);% 1, 3, 5...
    fid_seq= 0*fid_full;
    %  figure(1111);clf;plot(real(phase_in));dfasdf

    %mamip to get seq data algebra (phase 90)
%     skip_detail=1;
%     if skip_detail==0,
%         for cursor=1:size(fid_full,2),
%             if mod(cursor-1,4)==0,
%                 fid_seq(cursor)=+real(fid_full(1,cursor));
%             end
%             if mod(cursor-1,4)==1,
%                 fid_seq(cursor)=+imag(fid_full(1,cursor));
%             end
%             if mod(cursor-1,4)==2,
%                 fid_seq(cursor)=-real(fid_full(1,cursor));
%             end
%             if mod(cursor-1,4)==3,
%                 fid_seq(cursor)=-imag(fid_full(1,cursor));
%             end
%         end
%     else
        %same as above, but using phase
        phi_an=(90-phi_an_orig)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
        phases_pt=0:size(fid_seq,2)-1;
        phases=phi_an*phases_pt;
        fid_seq=real(fid_full.*exp(-j*phases));
%     end

offset=0;
other_method=1;
if other_method==1,
    fid_sim(1)=fid_sim(1)/2;
    fid_seq(1)=fid_seq(1)/2;
    
else 
offset=sum(sum(abs(fid_sim)))*(1/2*pi)*(1/sw);%%% important !!
end
if mainlooop==0,
    spectrum=fftshift(fft((fid_sim)));
    title(['sim /  SHR'])

end
if mainlooop==1,
title(['seq / Redfield / TPPI'])

    spectrumtmp=(fft((fid_seq)));
    %keep second part
    spectrum=spectrumtmp(:,size(spectrumtmp,2)/2+1:end);
 %   spectrum=spectrum+fliplr(spectrumtmp(:,1:size(spectrumtmp,2)/2));
end
if mainlooop==2,
    title(['States-TPPI'])

  phi_an=2*(90-phi_an_orig)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
        phases_pt=0:size(fid_sim,2)-1;
        phases=phi_an*phases_pt;
        fid_sim=(fid_sim.*exp(-j*phases));
       spectrum=(fft((fid_sim)));

end
title(['tppi'])

spectrum=spectrum-offset;



top=max(max(abs(spectrum)));
spectrum=spectrum/top;
incsw=sw/size(t,2);

scale=-sw/2+incsw/2:incsw :sw/2-incsw/2;



where_projy=sw/2;
where_projx=-1;

% plot3(where_projy+0*scale,imag(spectrum),real(spectrum),'k-'); hold on
plot3(where_projy+0*scale,imag(spectrum),real(spectrum),'k-'); hold on
plot3(scale,1+0*imag(spectrum),real(spectrum),'k-'); hold on
plot3(scale,imag(spectrum),real(spectrum),'g-','LineWidth',2); hold on
[azg bzg]=max(abs(spectrum));
%text(scale(bzg),imag(spectrum(bzg)),real(spectrum(bzg)),[num2str( 180/pi*angle(spectrum(bzg)) ,'%.1f') 'deg.' add_text]);
text(10,imag(spectrum(bzg)),real(spectrum(bzg)),[num2str( 180/pi*angle(spectrum(bzg)) ,'%.1f') 'deg.' add_text]);
erro_in_deg=180/pi*angle(spectrum(bzg));
%plot3(scale,0*imag(spectrum),1*real(spectrum),'LineWidth',2); hold on
plot3(scale,1*imag(spectrum),-1+0*real(spectrum),'k-'); hold on
plot3(sw/2*[-1 1],[ 0 0],[ 0 0],'k-')
plot3(sw/2*[-1 1],1+[ 0 0],[ 0 0],'k-')
plot3(sw/2*[1 1],[ 1 -1],[ 0 0],'k-')
axis([ where_projy*[-1 1] -1 1 -1 1])

for loop=1:size(height,2),
    what=abs(spectrum);
    tmp=(what(1,1:size(what,2)-1)-height(1,loop)).*(what(1,2:size(what,2)-0)-height(1,loop));
    list_i= find(tmp < 0);
    [a b]=max(what);
    [a c]=max(abs(spectrum));
    list_i=[list_i b c];
    for i=list_i,
        %  plot3(scale(1,i)*[1 1],[ imag(spectrum(1,i)) 0],[ real(spectrum(1,i)) 0],'k:')
        
        av_al=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))   *0.5 + 0.5* angle(real(spectrum(1,i+1))+j*imag(spectrum(1,i+1)));
        delta=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))  - angle(real(spectrum(1,i+1))+j*imag(spectrum(1,i+1)));
        %text(where_projy*[1],[ imag(spectrum(1,i)) ],[ real(spectrum(1,i)) ],[num2str(360/(2*pi)*av_al,'%.1f') '(' num2str(360/(2*pi)*delta,'%.1f')  ')']);
        addtext=[num2str(height(1,loop)*100,'%.0f') '% '] ;
        signal_me=0;
        if i==c,
            addtext='max(abs) ';
            av_al=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))  ;
            signal_me=1;
        end
        if i==b,
            addtext='max(real) ';
            av_al=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))  ;
            
        end
        if i==b && i==c,
            addtext='max(real & abs) ';
            av_al=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))  ;
            signal_me=1;
        end
        if signal_me,
            %           plot3(where_projy*[1 1],[ imag(spectrum(1,i)) 0],[ real(spectrum(1,i)) 0],'k:','LineWidth',2)
        else
            %           plot3(where_projy*[1 1],[ imag(spectrum(1,i)) 0],[ real(spectrum(1,i)) 0],'k:')
        end
        
        %         text(where_projy*[1],[ imag(spectrum(1,i)) ],[ real(spectrum(1,i)) ],[addtext num2str(360/(2*pi)*av_al,'%.1f') ]);
        %          text(scale(1,i)*[1],[ imag(spectrum(1,i)) ],[ real(spectrum(1,i)) ],[addtext num2str(360/(2*pi)*av_al,'%.1f') ]);
    end
end
plot3([sw/2 sw/2],[0 0 ],[0 1 ],'k','LineWidth',2);
plot3([sw/2 sw/2 -sw/2],[0 0 0],[1 -1 -1],'k');
plot3([-sw/2 -sw/2],[-1 1 ],-1+[0 0 ],'k');
plot3([-sw/2 sw/2 sw/2],[1 1 1],-1+[0 0 2],'k');
plot3([-sw/2 sw/2 sw/2],[1 1 -1],-1+[0 0 0],'k');
view([-24,22])
%view([-90,0])
axis([ where_projy*[-1 1] -1 1 -1 1])
if abs(dist_in_hz)>0.88 && abs(dist_in_hz)<0.89,
    print('-depsc','-tiff','-r600',[ 'Phase_error_near0.885' num2str(main_ratio)  '.eps']);%here
end
