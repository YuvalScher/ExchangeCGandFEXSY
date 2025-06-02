%A 3D Monte Carlo simulation of a two exchanging population and a third non-exchanging compartment (An extension of the FOV box):
%An "intracellular" spherical population (restricted compartment), and an
%extracullar compartment, suffering from turtuosity as function of the packing
%density. Several modes of packing are available.
%A FEXSY experiment is simulated. Filter is simulated as well.
tic
clear all

%User defined tissue model parameters:

dt = 10^-3;                                %Millisecond, simulation time step
spinsnum = 10^6;                           %Number of spins

                                           %Containing volume
FOVx = 5.5*4; FOVx_ExtraComp=FOVx+2;                             %Micrometers. In case of closed-packing, these numbers will be changed to the nearest values that can exactly accomodate the packing
FOVy = 5.5*4;
FOVz = 5.5*4;

a = 5;                                     %Spherical compartment diameter in micrometers
R = a/2;

De = 2;                                    %Extracellular diffusion coef. micrometer^2/millisec
Di = 2;                                    %Intracellular diffusion coef. micrometer^2/millisec

ki = 1/550;                                %The "real" exchange rate from intra to extra
kappa = (R*ki)/3;                          %Micrometer/millisecond. Permeability.


spheresnum = 43;   %Number of spheres to be generated, given random packing (otherwise determined by the volume of the FOV and of single sphere)

%Few packing arrangements are available (determined by the value of the parameter 'packing'):
% 0 =  Random packing 
% 1 =  The spheres centers lie on lines parallel to the axes. The number of spheres must equal (FOVx*FOVy*FOVz/a^3). The code fixes num of spheres automatically!
% 2 =  HCP (hexagonal closed-packing) lattice. The number of spheres must equal (FOVx*FOVy*FOVz/a^3 - FOVz/2a *(FOVy/a)). The code fixes num of spheres automatically!
% 3 =  FCC (cubic closed-packing) lattice. The number of spheres must equal a complicated expression. The code fixes num of spheres automatically!
% For HCP and FCC (AFTER GENERATING THE SPHERES), FOVz is then changed to fit num of layers, alpha is recalculated, and the center of all spheres is moved in the z-axis direction to achieve z-axis symmetry.

packing=0;

% When spination is 0 (see below), "scale" determines the size of the
% rectangular cuboid in which the spins are initially distributed, as
% compared to the FOV.

scale=1;
if scale>1 || scale<0
    error(1)
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  NMR Experiment Parameters  %%%%%%%%%%%%%%%%

%%%%%%%%% FEXSY %%%%%%%%%%%%%%%%%%
g1_user = 150;                           % filter block gradient pulse strengh g/cm
g2_user = linspace(5,150,30);            % measurement block gradient pulse strengh g/cm
d=2;                                     % gradient pulse duration ms
Tm=[8, 25, 50, 100, 150, 200, 250, 300, 400, 500, 750, 1000,1250]';    %ms
% Tm=[8, 25, 50,100,150,200,250,300];          %ms

tau=13;                  %ms
%%%%
g1 = g1_user*(10^(-8));  g2 = g2_user*(10^(-8));   %T/micrometer 
gyro= 42.57747 * 10^3;                             %gyromagentic ratio of protons 1/(ms T)
q1 = 2 * pi * gyro * g1 * d;          %q vector
q2 = 2 * pi * gyro * g2 * d;          %q vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AUTOMATIC FROM HERE

if packing ~= 0   %Fix FOV size such that packing 1 is exactly fitting in the box. The extra volume for packing 2-3 will be erased later
[FOVx, FOVy, FOVz, spheresnum] = UpdateSpheresNum(a,FOVx,FOVy,FOVz);    
end
                          
N=length(Tm)+1;    

%Intra- and extra-cellular diffusion coefficients. Gaussian 3D trajectory during dt.
%The direction will be chosen randomly in each iteration (see MovingSpins)
dre=sqrt(6*De*dt);
dri=sqrt(6*Di*dt);

%Spheres will be generated according to 'packing'.
%For Packing=0: Random (centers of) spheres position inside the given 3D FOV.
%Spheres are randomly generated one by one, while overlap between spheres is forbiden and results in recalculation of the position


[FOVx, FOVy, FOVz, FOV_volume, alpha, spheresnum, xcenter, ycenter, zcenter, rcenter] = CreateSpheres(a,FOVx,FOVy,FOVz,spheresnum,packing); %Create Spheres according the Packing
    
[xs,ys,zs] = sphere;  %For ploting 
for e=1:spheresnum    %Plot spheres with color
    figure (2)
    surf(R*xs+xcenter(e),R*ys+ycenter(e),R*zs+zcenter(e))
    hold on
end
xlim([-FOVx/2 FOVx_ExtraComp/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])

%Spination determines whether INITIALLY:
% 0 = The spins are randomly distributed within a rectangular cuboid which is the
% FOV scaled down by "scale" (0<scale<1);
% 1 = Spins are distributed equally between the spheres, and uniformly
% inside each sphere. No spins are initially in the extra-cellular space
% 2 = Spins are distributed only in a single sphere, closest as possible to the origin and uniformly inside that sphere. 


%%%%%%%%%%%%%%%% Prior to the filter, spination starts at 0  %%%%%%%%%%%%%%%%%

spination=0;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if spination==1
    sspss=floor(spinsnum/spheresnum);
    spinsnum=sspss*spheresnum;
elseif spination==2
    sspss=spinsnum;
elseif spination==0
    sspss=NaN;
end

ck=zeros(spinsnum,1);
for i=1:spinsnum
   ck(i) = 1 + floor((i-1)/sspss) ;    
end    
  

distsphere=zeros(spheresnum,1);   %Distance of the origin from ALL sphere centers
for e=1:spheresnum  
    distsphere(e) = sqrt((xcenter(e)).^2 + (ycenter(e)).^2 + (zcenter(e)).^2);
end
[~,central] = min(distsphere);     %Find the min distance, i.e., closest sphere to the origin.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%Begin Parfor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0=zeros(spinsnum,1);y0=zeros(spinsnum,1);z0=zeros(spinsnum,1);indicator0=zeros(spinsnum,1); indicator=zeros(spinsnum,1);  sps=zeros(spheresnum,1);xfin=zeros(spinsnum,1);yfin=zeros(spinsnum,1);zfin=zeros(spinsnum,1);

spinsnumi0=0;spinsnume0=0;spinsnumi=0;spinsnume=0;
E=zeros(length(q2),N); E1=zeros(length(q2),1); E2=zeros(length(q2),N-1); 

parfor s=1:spinsnum
    %%%
    [x0(s),y0(s),z0(s),indicator0(s),spheretrackeri,sps] = CreateSpins_tri(FOVx,FOVx_ExtraComp,FOVy,FOVz,R,spheresnum,xcenter,ycenter,zcenter,central,spination,ck,s,scale); %Create Spins according the Spination
    %%%

    x = x0(s); y = y0(s); z = z0(s); indicator(s)=indicator0(s); 
     loopE=zeros(length(q2),N);loopE1=zeros(length(q2),N); 
    %%%%%%%%%%%%%%%%%%%%%%%%%Finally: ITERATIONS
   
    %%%%%%%%filter block
    
   i=0;       
    while dt*i < tau+dt    

        [sps,x,y,z,indicator(s),spheretrackeri] = MovingSpins_tri(kappa,R,dre,dri,dt,De,Di,FOVx,FOVx_ExtraComp,FOVy,FOVz,xcenter,ycenter,zcenter,sps,spheresnum,x,y,z,indicator(s),spheretrackeri); 
    
        i=i+1
    end

    for j=2:N
         loopE1(:,j) = exp( - 1i * q1 * (x-x0(s)) );         
    end

    for j=1:length(q2)
         loopE1(j,1) = 1;         
    end
     
    %%%%%%%%mixing time and measurement block
  
    k=0;  counter=1;   
    while dt*k < max(Tm)+dt   %collect rms data at points DELTA

        [sps,x,y,z,indicator(s),spheretrackeri] = MovingSpins_tri(kappa,R,dre,dri,dt,De,Di,FOVx,FOVx_ExtraComp,FOVy,FOVz,xcenter,ycenter,zcenter,sps,spheresnum,x,y,z,indicator(s),spheretrackeri); 

          if dt*k == Tm(counter)    %collect rms data at points DELTA
             f=0;   

                     xP=x; yP=y; zP=z; indicatorP=indicator(s); spheretrackeriP=spheretrackeri;

                    while dt*f < tau+dt    %collect rms data at points DELTA

                       [sps,x,y,z,indicator(s),spheretrackeri] = MovingSpins_tri(kappa,R,dre,dri,dt,De,Di,FOVx,FOVx_ExtraComp,FOVy,FOVz,xcenter,ycenter,zcenter,sps,spheresnum,x,y,z,indicator(s),spheretrackeri); 

                        f=f+1
                    end

                    for j=1:length(q2)
                        loopE(j,counter+1) = exp(-1i * q2(j) *(x-xP)) * loopE1(j,counter+1) ;                
                    end    

                    if counter==N-1
                       for j=1:length(q2)  
                        loopE(j,1) = exp(-1i * q2(j)*(x-xP)) * loopE1(j,1) ;     
                       end 
                    end    

               counter=counter+1;                                                      
               x=xP; y=yP; z=zP;  indicator(s)=indicatorP;  spheretrackeri=spheretrackeriP;
          end



        k=k+1
    end
   E=E+loopE; 
   xfin(s)=x;yfin(s)=y;zfin(s)=z;
   s
end


E=E/spinsnum;

Enorm=zeros(length(q2),N);
for count = 1:N
    Enorm(:,count) = abs(E(:,count))/abs(E(1,count));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%End Parfor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b=zeros(length(q2),N);
  for i=1:N                  
 b(:,i) = q2.^2 .* tau;  
  end

ini=find(indicator0==1);
inen=find(indicator0==2);

ine = find(x0(inen)<=FOVx/2); 
inn = find(x0(inen)>FOVx/2);

for e=1:spheresnum      %Plot spheres with no color, when intial spins location will be genereated, they will be added to this figure
    figure (1)
    surf(R*xs+xcenter(e),R*ys+ycenter(e),R*zs+zcenter(e),'FaceColor', 'none')
    hold on
end
xlim([-FOVx/2 FOVx_ExtraComp/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])
plot3(x0(ini),y0(ini),z0(ini),'.red')
hold on
bg0=find(abs(x0(ine))>FOVx/2);
plot3(x0(ine(bg0)),y0(ine(bg0)),z0(ine(bg0)),'.green')
gb0=find(abs(x0(ine))<FOVx/2);
plot3(x0(ine(gb0)),y0(ine(gb0)),z0(ine(gb0)),'.blue')

plane_x_position = FOVx/ 2;
y_range = linspace(-FOVy/2,FOVy/2, 100);
z_range = linspace(-FOVy/2, FOVz/2, 100);
[Y, Z] = meshgrid(y_range, z_range);
X = plane_x_position * ones(size(Y));

surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none','FaceColor','black');






fe0=length(ine)/spinsnum;  
fn0=length(inn)/spinsnum;  
%%%%%%%%%%%%
%%%%%%%%%%%%
for e=1:spheresnum
    figure (3)
    surf(R*xs+xcenter(e),R*ys+ycenter(e),R*zs+zcenter(e),'FaceColor', 'none')
    hold on
end
xlim([-FOVx/2 FOVx_ExtraComp/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])

fni=find(indicator==1);
fnen=find(indicator==2);

fne = find(xfin(fnen)<=FOVx/2); 
fnn = find(xfin(fnen)>FOVx/2);

figure(3)
plot3(xfin(fni),yfin(fni),zfin(fni),'.red')
hold on
bg=find(abs(xfin(fne))>FOVx/2);
plot3(xfin(fne(bg)),yfin(fne(bg)),zfin(fne(bg)),'.green')
gb=find(abs(xfin(fne))<FOVx/2);
plot3(xfin(fne(gb)),yfin(fne(gb)),zfin(fne(gb)),'.blue')

plane_x_position = FOVx/ 2;
y_range = linspace(-FOVy/2,FOVy/2, 100);
z_range = linspace(-FOVy/2, FOVz/2, 100);
[Y, Z] = meshgrid(y_range, z_range);
X = plane_x_position * ones(size(Y));

surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none','FaceColor','black');

fef=length(fne)/spinsnum; 
fnf=length(fnn)/spinsnum; 
%%%%%%%%%%%%%%%%%%%%%%%%


% First, a plot of ALL experiments.

for i=2:N 
figure(4)
plot(b(:,i),log(Enorm(:,i)),'.-','MarkerSize',20,'LineWidth',2)
hold on
end
plot(b(:,1),log(Enorm(:,1)),'.-k','MarkerSize',20,'LineWidth',3)


for i = 1:N-1
    legendCell(i) = strcat(string(num2cell(Tm(i))),'ms');
end
legendCell(N)=('No Filter');
legend(legendCell,'Location','eastoutside')
ylabel('$ln(E)$','Interpreter','latex')
xlabel('$\tilde{q}^2\Delta \hspace{1ex} [s/m^2]$','Interpreter','latex')

ax=gca;
ax.FontSize=14;
ax.FontWeight='bold';

% title('x-FEXSY') %% USER!!


% The next step is finding the diffusion coefficient of each experiment
% when b---> 0 . The limit of this condition is left to the user to define
% (the var 'finish)

ft = fittype({'x','1'},'coefficients',{'a','b'}); % Define ft as linear fit
coeffnames(ft);

start=1; 
finish=8;  % USER!!
sece=8;
tri=4;

yL = get(gca,'YLim');
line([b(sece) b(sece)],yL,'Color','black','LineStyle',':','LineWidth',2);

D=zeros(N,1);
Rsquare=zeros(N,1);
const=zeros(N,1);


%Fitting:
for i=1:N

EE(:,i)=log(E(start:finish, i));
B(:,i)=b(start:finish, i);
 
    
[fitobject,gof] =fit(B(:,i),EE(:,i),ft);

D(i,1)=-fitobject.a;
const(i,1)=fitobject.b;
Rsquare(i,1)=gof.rsquare;

end
 
for i=1:(N-1)
Dnorm(i,1)=D(i+1,1)/D(1,1);     % Normalize  D Values with respect to 'no filter' experiment (this exp is no longer needed)
end

%Fitting the D values to a Growth exponential function. plotting the result
%t_res is the apparent residual time. 1/t_res=k is the in-exchange cnst.
 fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[inf,inf,1],...
               'StartPoint',[1 300 1]);
ft2 = fittype('c-a*exp(-x/b)','coefficient',{'a','b','c'},'options',fo);
coeffnames(ft2);
[fitobject,gof] =fit(Tm,Dnorm(:,1),ft2);


t_res=fitobject.b;
AXR=1/(t_res*0.001);
a=fitobject.a;
c=fitobject.c;
Rsquare2=gof.rsquare;


% %%%%%%%%5
%  fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,0],...
%                'Upper',[inf,inf,1],...
%                'StartPoint',[1, 300, 0.5]);
% ft3 = fittype('1-s-a*exp(-x/b)','coefficient',{'a','b','s'},'options',fo);
% coeffnames(ft2);
% [fitobject,gof] =fit(Tm,Dnorm(:,1),ft3);
% 
% 
% t_res3=fitobject.b;
% AXR=1/(t_res*0.001);
% a3=fitobject.a;
% s3=fitobject.s;
% Rsquare2_3=gof.rsquare;
% 
% %%%%%%%%5



%the USER can decide truncating the fit and plot at a certain point:
finish2=N-1; %If Finish2 is set to N-1, no truncation will be performed
TM=Tm(1:finish2);
DNORM=Dnorm(1:finish2);


figure (5)
plot(TM,DNORM,'kO')
hold on
Fitra=fitobject.c-fitobject.a*exp(-TM/fitobject.b);
plot(TM,Fitra,'k-')

ylabel('Normalized ADC')
xlabel('Mixing Time [ms]')

ax=gca;
ax.FontSize=14;
ax.FontWeight='bold';


%%% For cftool
LOOP=30;
for  i=1:length(Tm)
EEE((LOOP*(i-1)+1):(i*LOOP))=abs(Enorm(:,i+1));
BBB((LOOP*(i-1)+1):(i*LOOP))=b(:,i+1);
QQQ((LOOP*(i-1)+1):(i*LOOP))=q2(1,:)/2/pi;
  for j=0:(LOOP-1)
  TTT(j+(LOOP*(i-1)+1))= Tm(i);
  end
end


t_res
toc
save('workspace3')