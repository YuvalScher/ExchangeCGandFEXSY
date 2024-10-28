%A 3D Monte Carlo simulation of a two exchanging population: 
%An "intracellular" spherical population (restricted compartment), and an
%extracullar compartment, suffering from turtuosity as function of the packing
%density. Several modes of packing are available.
%A CG-PFG experiment is simulated.

clear all

%User defined tissue model parameters:

dt = 10^-3;                                 %Millisecond, simulation time step
spinsnum = 10^4;                            %Number of spins

                                            %Containing volume
FOVx = 5.5*2;                               %Micrometers. In case of closed-packing, these numbers will be changed to the nearest values that can exactly accomodate the packing
FOVy = 5.5*2;
FOVz = 5.5*2;

a = 5;                                      %Spherical compartment diameter in micrometers
R = a/2;

De = 2;                                     %Extracellular diffusion coef. micrometer^2/millisec
Di = 2;                                     %Intracellular diffusion coef. micrometer^2/millisec

ki = 0;                                 %The "real" exchange rate from intra to extra
kappa = (R*ki)/3;                           %Micrometer/millisecond. Permeability.


spheresnum = 1;   %Number of spheres to be generated, given random packing (otherwise determined by the volume of the FOV and of single sphere)


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

scale=0.5;
if scale>1 || scale<0
    error(1)
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  NMR Experiment Parameters  %%%%%%%%%%%%%%%%

%%%%%%%%% CG %%%%%%%%%%%%%%%%%%
g_user = linspace(1,500,60);     % Measurement block gradient pulse strengh g/cm
d=2;                             % Gradient pulse duration, ms
tau=100;                         % ms
%%%%
g = g_user*(10^(-8));            %T/micrometer 
gyro= 42.57747 * 10^3;           %gyromagentic ratio of protons, 1/(ms T)
q = 2 * pi * gyro * g * d;       %q vector, rad/micrometer


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
                          
%Intra- and extra-cellular diffusion coefficients. Gaussian 3D trajectory during dt.
%The direction will be chosen randomly in each iteration (see MovingSpins)
dre=sqrt(6*De*dt);
dri=sqrt(6*Di*dt);

%Spheres will be generated according to 'packing'.
%For Packing=0: Random (centers of) spheres position inside the given 3D FOV.
%Spheres are randomly generated one by one, while overlap between spheres is forbiden and results in recalculation of the position

[FOVx, FOVy, FOVz, FOV_volume, alpha, spheresnum, xcenter, ycenter, zcenter, rcenter] = CreateSpheres(a,FOVx,FOVy,FOVz,spheresnum,packing); %Create Spheres according the Packing

[xs,ys,zs] = sphere;    %For ploting 
for e=1:spheresnum      %Plot spheres with no color, when the intial spins locations will be genereated, they will be added to this figure
    figure (1)
    surf(R*xs+xcenter(e),R*ys+ycenter(e),R*zs+zcenter(e),'FaceColor', 'none')
    hold on
end
xlim([-FOVx/2 FOVx/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])


for e=1:spheresnum    %Plot spheres with color
    figure (2)
    surf(R*xs+xcenter(e),R*ys+ycenter(e),R*zs+zcenter(e))
    hold on
end
xlim([-FOVx/2 FOVx/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])

%Spination determines whether INITIALLY:
% 0 = The spins are randomly distributed within a rectangular cuboid which is the
% FOV scaled down by "scale" (0<scale<1);
% 1 = Spins are distributed equally between the spheres, and uniformly
% inside each sphere. No spins are initially in the extra-cellular space
% 2 = Spins are distributed only in a single sphere, closest as possible to the origin and uniformly inside that sphere. 

%%%%%%%%%%%%%%%%%%%%% For Diffraction check

spination=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if spination==1   %For later, when distributing the spins, see CreateSpins
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

distsphere=zeros(spheresnum,1);        %Distance of the origin from ALL sphere centers
for e=1:spheresnum           
    distsphere(e) = sqrt((xcenter(e)).^2 + (ycenter(e)).^2 + (zcenter(e)).^2);
end
[~,central] = min(distsphere);          %Find the min distance, i.e., closest sphere to the origin.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%Begin Parfor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=zeros(spinsnum,1); y0=zeros(spinsnum,1); z0=zeros(spinsnum,1); indicator0=zeros(spinsnum,1); indicator=zeros(spinsnum,1); sps=zeros(spheresnum,1);
spinsnumi0=0; spinsnume0=0; spinsnumi=0; spinsnume=0;
E=zeros(length(q),1);

parfor s=1:spinsnum

      %%%
      [x0(s),y0(s),z0(s),indicator0(s),spheretrackeri,sps] = CreateSpins(FOVx,FOVy,FOVz,R,spheresnum,xcenter,ycenter,zcenter,central,spination,ck,s,scale); %Create Spins according to Spination
      %%%
    
      x = x0(s); y = y0(s); z = z0(s); indicator(s)=indicator0(s); 
    
      %%%%%%%%%%%%%%%%%%%%%%%%%Finally: ITERATIONS
      counter=1; loopE=zeros(length(q),1);
      
        i=0;       
        while dt*i < tau+dt     
            
            [sps,x,y,z,indicator(s),spheretrackeri] = MovingSpins(kappa,R,dre,dri,dt,De,Di,FOVx,FOVy,FOVz,xcenter,ycenter,zcenter,sps,spheresnum,x,y,z,indicator(s),spheretrackeri); 

        i=i+1
        end      
      
      
      
      for j=1:length(q)
         loopE(j) =    exp( - 1i * q(j)* (x-x0(s)) )   ;           
      end    

       E=E+loopE;
       xfin(s)=x;yfin(s)=y;zfin(s)=z;
       s
end

E=E/spinsnum; 


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%End Parfor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ini=find(indicator0==1);   %Calcluate inital spins fractions in the intra-and extra-cellular domains
ine=find(indicator0==2);
fei=length(ine)/spinsnum;  

figure(1)                  %Plot spins initial locations
plot3(x0(ini),y0(ini),z0(ini),'.red')
hold on
plot3(x0(ine),y0(ine),z0(ine),'.blue')


for e=1:spheresnum    %Plot spins final locations
    figure (3)
    surf(R*xs+xcenter(e),R*ys+ycenter(e),R*zs+zcenter(e),'FaceColor', 'none')
    hold on
end
xlim([-FOVx/2 FOVx/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])


%Calcluate final spins fractions in the intra-and extra-cellular domains
fni=find(indicator==1);
fne=find(indicator==2);
fef=length(fne)/spinsnum;  

figure(3)   %Plot spins final locations
plot3(xfin(fni),yfin(fni),zfin(fni),'.red')
hold on
plot3(xfin(fne),yfin(fne),zfin(fne),'.blue')


figure(4)

plot(q*R/2/pi,log((E)),'.k','MarkerSize',20,'LineWidth',3)

ylabel('$E$','Interpreter','latex')
xlabel('$\tilde{q}R$','Interpreter','latex')

hold on


analytical =  (abs(3*(q*R .* cos(q*R) - sin(q*R))./((q*R).^3))).^2;

plot(q*R/2/pi,log(analytical),'-b','MarkerSize',20,'LineWidth',3)

