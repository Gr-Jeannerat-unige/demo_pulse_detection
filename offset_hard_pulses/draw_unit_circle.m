function draw_unit_circle(offset)
if nargin==0
    offset=[0;0;0];
end
%figure(1)
[az,el] = view;
inca=2*pi/1000;
list_an=0:inca:2*pi;

%% generate main circle
%figure(2);clf;

l=[1 0 0 ;0 1 0];
p=zeros(3,size(list_an,2));
p2=zeros(3,size(list_an,2));
for lo=1:3
    p(lo,:)=l(1,lo)*cos(list_an')+l(2,lo)*sin(list_an');
end
%plot3(p(1,:),p(2,:),p(3,:),'k-');hold on

%first rotation
now=90-el;
%now=-az;
l1=[1 0 0 ;0 0 0;0 0  1]*cos(now/180*pi);
l2=[0 0 1 ;0 0 0;-1 0 0]*sin(now/180*pi);
l3=[0 0 0 ;0 1 0; 0 0 0];
for lo1=1:3
    for lo2=1:3
      %  disp([num2str(l1(lo1,lo2)) ' '  num2str(l2(lo1,lo2)) ' '  num2str(l3(lo1,lo2)) ]);
        p2(lo1,:)=p2(lo1,:)+p(lo2,:)*l1(lo1,lo2)+p(lo2,:)*l2(lo1,lo2)+p(lo2,:)*l3(lo1,lo2);
    end
end

%plot3(p2(1,:)+offset(1,1),p2(2,:)+offset(2,1),p2(3,:)+offset(3,1),'k:');hold on

p=p2;
p2=p2*0;
%second rotation

now=90-az;
%now=-az;
l1=[1 0 0 ;0  1 0;0 0 0]*cos(now/180*pi);
l2=[0 1 0 ;-1 0 0;0 0 0]*sin(now/180*pi);
l3=[0 0 0 ;0 0 0; 0 0 1];
for lo1=1:3
    for lo2=1:3
      %  disp([num2str(l1(lo1,lo2)) ' '  num2str(l2(lo1,lo2)) ' '  num2str(l3(lo1,lo2)) ]);
        p2(lo1,:)=p2(lo1,:)+p(lo2,:)*l1(lo1,lo2)+p(lo2,:)*l2(lo1,lo2)+p(lo2,:)*l3(lo1,lo2);
    end
end
plot3(p2(1,:)+offset(1,1),p2(2,:)+offset(2,1),p2(3,:)+offset(3,1),'k-');hold on
end