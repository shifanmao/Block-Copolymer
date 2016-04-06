clear;
NQ=3;

q=1;
qv=zeros(NQ,3);

if NQ==2 %Lamellar
    qv(1,:)=q*[1,0,0];
    qv(2,:)=-qv(1,:);
elseif NQ==3 %cylinder
%     qv(1,:)=q*[1,0,0];
%     qv(2,:)=q*[-.5,sqrt(3)/2,0];
%     qv(3,:)=q*[-.5,-sqrt(3)/2,0];

    qv(1,:)=[1,0,0];
    qv(2,:)=(1/2)*[-1,+sqrt(3),0];
    qv(3,:)=(1/2)*[-1,-sqrt(3),0];
    qv(4,:)=-qv(1,:);
    qv(5,:)=-qv(2,:);
    qv(6,:)=-qv(3,:);
NQ=6;

elseif NQ==6 %BCC
    qv(1,:)=q*[1,0,0];
    qv(2,:)=q*[0,1,0];
    qv(3,:)=q*[0,0,1];
    qv(4,:)=-qv(1,:);
    qv(5,:)=-qv(2,:);
    qv(6,:)=-qv(3,:);
elseif NQ==50 %random
    for ii=1:2:NQ
        randv=rand(1,3)-0.5;
        qv(ii,:)=q*randv/norm(randv);
        qv(ii+1,:)=-qv(ii,:);
    end
end

NX=20;
x=linspace(0,10,NX);
y=x;
z=x;
rho=zeros(NX,NX,NX);
for ix=1:NX
    ix
    for iy=1:NX
        for iz=1:NX
            r=[x(ix),y(iy),z(iz)];
            rhov=0;
            for ii=1:NQ
                rhov=rhov+exp(1i*dot(r,qv(ii,:)));
            end
            rho(ix,iy,iz)=real(rhov);
        end
    end
end

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1)
set(gca,'fontsize',15);

NISO=5;
isovalues=linspace(min(min(min(rho))),max(max(max(rho))),NISO)
for ii=1:NISO
    isovalue=isovalues(ii);
    col = (ii-1)/(NISO-1);
    surf1=isosurface(x,y,z,rho,isovalue);
    p1 = patch(surf1);
    isonormals(x,y,z,rho,p1);
    set(p1,'FaceColor',[col 0 1-col],'EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface
end

% isovalue=2;
% surf1=isosurface(x,y,z,rho,isovalue);
% p1 = patch(surf1);
% isonormals(x,y,z,rho,p1);
% set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface
% 
% isovalue=1;
% surf1=isosurface(x,y,z,rho,isovalue);
% p1 = patch(surf1);
% isonormals(x,y,z,rho,p1);
% set(p1,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface
% 
% isovalue=0;
% surf1=isosurface(x,y,z,rho,isovalue);
% p1 = patch(surf1);
% isonormals(x,y,z,rho,p1);
% set(p1,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface

daspect([1,1,1])
view([20,10,10]); axis tight
camlight; lighting gouraud

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
set(gca,'fontsize',18)
title('From side')

subplot(1,2,2)
set(gca,'fontsize',15);
% 
% isovalue=rho*2;
% surf1=isosurface(x,y,z,rho,isovalue);
% p1 = patch(surf1);
% isonormals(x,y,z,rho,p1);
NISO=5;
isovalues=linspace(min(min(min(rho))),max(max(max(rho))),NISO)
for ii=1:NISO
    isovalue=isovalues(ii);
    col = (ii-1)/(NISO-1);
    surf1=isosurface(x,y,z,rho,isovalue);
    p1 = patch(surf1);
    isonormals(x,y,z,rho,p1);
    set(p1,'FaceColor',[col 0 1-col],'EdgeColor','none','FaceAlpha',0.8); % set the color, mesh and transparency level of the surface
end

set(p1,'FaceColor','blue','EdgeColor','none','FaceAlpha',1.); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
% view(2); axis tight
view([88,0]); axis tight
camlight; lighting gouraud
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
set(gca,'fontsize',18)
title('From top')