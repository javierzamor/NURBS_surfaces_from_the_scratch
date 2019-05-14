function NURBS_from_the_scratch()                 

    set(0,'DefaultFigureWindowStyle','docked')
    clear;
    fclose all;
    close all;
    clc;
    format short;
    
    sonidito=audioread('sound5.wav');
    player=audioplayer(sonidito,44100);
   
    play(player);
    for i=1:7
        [M,U,V,w,nx,ny]=cargar_superficie(i,1);
        c=nurbs(M,U,V,w,10);
        ver(c,M,U,V,w,nx,ny,i+1,player);  
        if(isplaying(player)==0);play(player);end
    end
    
    esfera(player);

end

   
  
function c_todo=nurbs(M,U,V,w,resolucion)
    
    Nx=2*size(M,1)-size(U,2)-1;
    Ny=2*size(M,2)-size(V,2)-1;
    nx=size(U,2)-size(M,1)+1;
    ny=size(V,2)-size(M,2)+1;
    
    c_todo=zeros(0,resolucion+1,size(M,3));
    for j=1:Ny
        for i=1:Nx
            ax=U(i+nx-1);    
            bx=U(i+nx);
            ay=V(j+ny-1);    
            by=V(j+ny);
            if (ax~=bx & ay~=by)
                u=ax:(bx-ax)/resolucion:bx;
                v=ay:(by-ay)/resolucion:by;
                c=super_r_boor(M,U,V,w,u,v);
                c_todo=[c_todo;c];
            end
        end
    end
end
    
function [y]=super_r_boor(M,U,V,w,u,v)
        
    M1(:,:,1)=w;
    for cont=1:size(M,3)
        M1(:,:,cont+1)=M(:,:,cont).*w;
    end
    c=superboor(M1,U,V,u,v);
    for cont=1:(size(c,3)-1)
        y(:,:,cont)=c(:,:,cont+1)./c(:,:,1);
    end
end
            
function [y]=superboor(M,U,V,u,v)
    
    c=zeros(0,size(u,2),size(M,3));
    for i=1:size(M,2)
        MM=squeeze(M(:,i,:));
        c_tmp=multideboor(MM,U,u);
        A=reshape(c_tmp,1,size(c_tmp,1),size(c_tmp,2));
        c=[c;A];    
    end
       
    d=zeros(0,size(v,2),size(M,3));
    for j=1:size(c,2)
        cc=squeeze(c(:,j,:));
        d_tmp=multideboor(cc,V,v);
        A=reshape(d_tmp,1,size(d_tmp,1),size(d_tmp,2));
        d=[d;A];    
    end
    y=d;
end

       
function [y]=multideboor(M,U,t)

    N=2*size(M,1)-size(U,2)-1;
    n=size(U,2)-size(M,1)+1;
    pivote=(t(1)+t(size(t,2)))/2; %para mantener la compatibilidad con que t sea un único valor
    
    contador=1;
    while (pivote>U(contador))
        contador=contador+1;
    end
    u_inicio=contador-n;
    m_inicio=contador-n;
    
    UU=U(u_inicio);
    for cont=(u_inicio+1):1:(u_inicio+2*n-1)
        UU=[UU U(cont)];
    end
    MM=M(m_inicio,:);
    for cont=(m_inicio+1):1:(m_inicio+n)
        MM=[MM;M(cont,:)];
    end
    y=deboor(MM,UU,t);
end

function [y]=deboor(M,U,t)             

    dimens=size(M,2);
    n=size(M,1)-1;
    a=ones(1,dimens);    
    t=(a(:)*t)';
    dim_t=size(t,1);
    d.flag=zeros(1,1+n+n*n);
    d.tmp={[zeros(dim_t,dimens)]};
    for contador=1:(1+n+n*n)
        d.tmp={[zeros(dim_t,dimens)] d.tmp{:}};
    end
    [d]=D(M,U,t,n,1,d);
    y=d.tmp{n*n+n+1};
end

function [d]=D(M,U,t,r,i,d) 
    
    dim_t=size(t,1);
    n=size(M,1)-1;
    if r==0
        a=ones(1,dim_t);    
        d.tmp{i}=a(:)*M(i,:);
        d.flag(i)=1;
    else
        if (d.flag((r-1)*(n+1)+i))==0
            [d]=D(M,U,t,r-1,i,d);
        end
        if (d.flag((r-1)*(n+1)+i+1))==0
            [d]=D(M,U,t,r-1,i+1,d);
        end
        a=d.tmp{(r-1)*(n+1)+i};
        b=d.tmp{(r-1)*(n+1)+i+1};
        u=t;
        t=(u-U(i+r-1))./(U(i+n)-U(i+r-1));
        d.tmp{r*(n+1)+i}=(1-t).*a+t.*b;
        d.flag(r*(n+1)+i)=1;
    end
end

function [U]=construir_nudos(M,n)
   
    N=size(M,1)-n;
    U=zeros(1,2*n+N-1);
    for i=1:(size(M,1)-1)
        distn(i)=sqrt(sum((M(i+1,:)-M(i,:)).^2));
    end
    if (sum(distn)==0)
        distn=distn+N/2;
    end
    for i=1:n
        U(i)=0;
    end
    for i=1:N-1
        U(i+n)=((i/(N))*distn(i+1)+sum(distn(1:i)))*(N)/sum(distn);
    end
    for i=(N+n):(2*n+N-1)
        U(i)=N;
    end
end


function ver(c,M,U,V,w,nx,ny,i,player)
    
    if(isplaying(player)==0);play(player);end
    Nx=size(M,1)-nx;
    Ny=size(M,2)-ny;
    U
    V
    w
    disp('nx, Nx, ny, Ny =');
    disp([nx Nx ny Ny]);
    
    plot3(c(:,:,1),c(:,:,2),c(:,:,3),'.','LineWidth',2,'color','black');
    axis equal;
    hold on;
    axis off;
    plot3(M(:,:,1),M(:,:,2),M(:,:,3),'s','LineWidth',2,'color','blue');
    plot3(c(:,:,1),c(:,:,2),c(:,:,3),c(:,:,1)',c(:,:,2)',c(:,:,3)','LineWidth',1,'color','red');
    
    hold off;
    for cont1=0:11
        for cont2=0:1:30
            az=cont1*30+cont2;
            el=-90+az/2;
            view(az,el)
            drawnow
            if(isplaying(player)==0);play(player);end
        end
    end
    view(360,90)
    drawnow()
end

   
function esfera(player)
    if(isplaying(player)==0);play(player);end
    hold off;
    clc;
    disp('esfera');
    
    if(isplaying(player)==0);play(player);end
    [M,U,V,w,nx,ny]=cargar_superficie(7,0.8);
    e=nurbs(M,U,V,w,20);
    Nx=size(M,1)-nx;
    Ny=size(M,2)-ny;
    U
    V
    w
    disp('nx, Nx, ny, Ny =');
    disp([nx Nx ny Ny]);
   
    axis equal;
    axis off;
    view(360,90);
    hold on;
    drawnow;
    hand6=plot3(e(:,:,1),e(:,:,2),e(:,:,3),'.','Markersize',1,'color','black');


    drawnow;
    if(isplaying(player)==0);play(player);end
 
    for az=0:0.5:360
        el=-90+az/2;
        view(az,el)
        drawnow
        if(isplaying(player)==0);play(player);end
    end
    for az=361:0.5:540
        el=-90+az/2;
        view(az,el)
        drawnow
        if(isplaying(player)==0);play(player);end
    end
%     cont=1;
    for cont=1:-0.001:0.1
        zlim([-cont cont]);
        xlim([-cont cont]);
        ylim([-cont cont]);
        drawnow;
        if(isplaying(player)==0);play(player);end
    end
    
    drawnow;
    disp('');

end

function [M,U,V,w,nx,ny]=cargar_superficie(i,parametro)

    switch i
    case 1 % 1 superficie
        clc;
        disp('superficie de Bezier');
        nx=3;
        ny=3;
        M(:,:,1)=[-15 -15 -15 -15;-5 -5 -5 -5;5 5 5 5;15 15 15 15];
        M(:,:,2)=[0 5 5 0;5 5 5 5;5 5 5 5;0 5 5 0];
        M(:,:,3)=[15 5 -5 -15;15 5 -5 -15;15 5 -5 -15;15 5 -5 -15];
        w=ones(size(M,1),size(M,2));
        U=construir_nudos(squeeze(M(:,1,:)),nx);
        V=construir_nudos(squeeze(M(1,:,:)),ny);
    case 2 % 2 juntas
        clc;
        disp('dos superficies unidas, splines');
        nx=3;
        ny=3;
        M(:,:,1)=[-15 -15 -15 -15;-5 -5 -5 -5;5 5 5 5;15 15 15 15;25 25 25 25;35 35 35 35;45 45 45 45];
        M(:,:,2)=[0 5 5 0;5 5 5 5;5 5 5 5;0 5 5 0;5 5 5 5;5 5 5 5;0 5 5 0];
        M(:,:,3)=[15 5 -5 -15;15 5 -5 -15;15 5 -5 -15;15 5 -5 -15;15 5 -5 -15;15 5 -5 -15;15 5 -5 -15];
        w=ones(size(M,1),size(M,2));
        U=construir_nudos(squeeze(M(:,1,:)),nx);
        V=construir_nudos(squeeze(M(1,:,:)),ny);
    case 3 % cilindro cuadrado
        clc;
        disp('cilindro cuadrado, sin asignar pesos');
        nx=3;
        ny=2;
        M(:,:,3)=[0 1 2;0 1 2;0 1 2;0 1 2;0 1 2;0 1 2;0 1 2;0 1 2;0 1 2;0 1 2;0 1 2];
        M(:,:,2)=[-1 -1 -1;-1 -1 -1;0 0 0;1 1 1;1 1 1;1 1 1;0 0 0;-1 -1 -1;-1 -1 -1;-1 -1 -1;0 0 0];
        M(:,:,1)=[0 0 0;1 1 1;1 1 1;1 1 1;0 0 0;-1 -1 -1;-1 -1 -1;-1 -1 -1;0 0 0;1 1 1;1 1 1];
        w=ones(size(M,1),size(M,2));
        U=[0 1 2 3 4 5 6 7 8 9 10 11 12];
        V=construir_nudos(squeeze(M(1,:,:)),ny);
    case 4 % cilindro 
        clc;
        disp('cilindro constuido como una superficie de spline continuos');
        nx=3;
        ny=3;
        M(:,:,3)=[0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5;0 1 2 3 4 5];
        M(:,:,2)=[-1 -1 -1 -1 -1 -1;-1 -1 -1 -1 -1 -1;0 0 0 0 0 0;1 1 1 1 1 1;1 1 1 1 1 1;1 1 1 1 1 1;0 0 0 0 0 0;-1 -1 -1 -1 -1 -1;-1 -1 -1 -1 -1 -1;-1 -1 -1 -1 -1 -1;0 0 0 0 0 0];
        M(:,:,1)=[0 0 0 0 0 0;1 1 1 1 1 1;1 1 1 1 1 1;1 1 1 1 1 1;0 0 0 0 0 0;-1 -1 -1 -1 -1 -1;-1 -1 -1 -1 -1 -1;-1 -1 -1 -1 -1 -1;0 0 0 0 0 0;1 1 1 1 1 1;1 1 1 1 1 1];
        w=ones(size(M,1),size(M,2));
        w(2,:)=sqrt(2)/4;
        w(4,:)=sqrt(2)/4;
        w(6,:)=sqrt(2)/4;
        w(8,:)=sqrt(2)/4;
        w(10,:)=sqrt(2)/4;
        U=[0 1 2 3 4 5 6 7 8 9 10 11 12];
        V=construir_nudos(squeeze(M(1,:,:)),ny);
    case 5 % otro cilindro, otros w y U 
        clc;
        disp('cilindro constuido mediante una concatenación de splines');
        nx=3;
        ny=2;
        M(:,:,3)=[-1 0 1;-1 0 1;-1 0 1;-1 0 1;-1 0 1;-1 0 1;-1 0 1;-1 0 1;-1 0 1];
        M(:,:,2)=[-1 -1 -1;-1 -1 -1;0 0 0;1 1 1;1 1 1;1 1 1;0 0 0;-1 -1 -1;-1 -1 -1];
        M(:,:,1)=[0 0 0;1 1 1;1 1 1;1 1 1;0 0 0;-1 -1 -1;-1 -1 -1;-1 -1 -1;0 0 0];
        w=ones(size(M,1),size(M,2));
        w(2,:)=sqrt(2)/2;
        w(4,:)=sqrt(2)/2;
        w(6,:)=sqrt(2)/2;
        w(8,:)=sqrt(2)/2;
        U=[0 0 1 1 2 2 3 3 4 4];
        V=construir_nudos(squeeze(M(1,:,:)),ny);
    case 6 % hemiesfera 
        clc;
        disp('hemiesfera');
        nx=3;
        ny=3;
        M(:,:,3)=[0 0 0 0 0;0 1 1 1 0;0 1 1 1 0;0 1 1 1 0;0 0 0 0 0];
        M(:,:,1)=[0 0 0 0 0;1 1 0 -1 -1;1 1 0 -1 -1;1 1 0 -1 -1;0 0 0 0 0];
        M(:,:,2)=[-1 -1 -1 -1 -1;-1 -1 -1 -1 -1;0 0 0 0 0;1 1 1 1 1;1 1 1 1 1];
        a=ones(1,size(M,2));
        b=[1 sqrt(2)/2 1 sqrt(2)/2 1]';
        w=b*a;
        w(3,2:2:5)=sqrt(2)/2;
        U=[0 0 1 1 2 2];
        V=[0 0 1 1 2 2];
    case 7 % esfera 
        R=parametro; %el radio
        nx=2;
        ny=2;
        M(:,:,1)=[0 0 0 0 0 0 0 0 0;0 1 1 1 0 -1 -1 -1 0;0 1 1 1 0 -1 -1 -1 0;0 1 1 1 0 -1 -1 -1 0;0 0 0 0 0 0 0 0 0];
        M(:,:,3)=[0 0 0 0 0 0 0 0 0;1 1 0 -1 -1 -1 0 1 1;1 1 0 -1 -1 -1 0 1 1;1 1 0 -1 -1 -1 0 1 1;0 0 0 0 0 0 0 0 0];
        M(:,:,2)=[-1 -1 -1 -1 -1 -1 -1 -1 -1;-1 -1 -1 -1 -1 -1 -1 -1 -1;0 0 0 0 0 0 0 0 0;1 1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1 1];
        M=R*M;
        a=ones(1,size(M,2));
        b=[1 sqrt(2)/2 1 sqrt(2)/2 1]';
        w=b*a;
        w(3,2:2:9)=sqrt(2)/2;
        U=[0 0 1 1 2 2];
        V=[0 0 1 1 2 2 3 3 4 4];
    end
end

