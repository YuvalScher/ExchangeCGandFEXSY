%A 3D Monte Carlo simulation of a two exchanging population: 
%An "intracellular" spherical population (restricted compartment), and an
%extracullar compartment, suffering from turtuosity as function of the packing
%density. Several modes of packing are available.
%A CG-PFG experiment is simulated.

clear all

%User defined tissue model parameters:

dt = 10^-3;                                 %Millisecond, simulation time step
spinsnum = 10^5;                            %Number of spins

                                            %Containing volume
FOVx = 5.5*4;                               %Micrometers. In case of closed-packing, these numbers will be changed to the nearest values that can exactly accomodate the packing
FOVy = 5.5*4;
FOVz = 5.5*4;

a = 5;                                      %Spherical compartment diameter in micrometers
R = a/2;

De = 2;                                     %Extracellular diffusion coef. micrometer^2/millisec
Di = 2;                                     %Intracellular diffusion coef. micrometer^2/millisec

ki = 1/550;                                 %The "real" exchange rate from intra to extra
kappa = (R*ki)/3;                           %Micrometer/millisecond. Permeability.


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

scale=0.5;
if scale>1 || scale<0
    error(1)
end    

exp_time = 600;                                  %millisecond, largest diffusion time
NDELTA = 1000;                                   %number of diffusion times to be sampled

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  NMR Experiment Parameters  %%%%%%%%%%%%%%%%

%%%%%%%%% CG %%%%%%%%%%%%%%%%%%
g_user = 100;                    % gradient pulse strengh, g/cm
d=2;                             % gradient pulse duration, ms
%%%%
g = g_user*(10^(-8));            %T/micrometer 
gyro= 42.57747 * 10^3;           %gyromagentic ratio of protons, 1/(ms T)
q = 2 * pi * gyro * g * d;       %q vector, rad/micrometer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                          
DELTA=linspace(0,exp_time,NDELTA+1).';

for fu=1:NDELTA                  %Making DELTA multiples of dt by slightly changing them (supposed to be negligible)                
    dev=DELTA(fu)/dt;           
    dev=round(dev);
    DELTA(fu)=dev*dt;
end
exp_time=DELTA(NDELTA+1);
N=exp_time/dt;                    %Number of iterations

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

%%%%%%%%%%%%%%%%%%%%% For CG -- spination starts at 0  %%%%%%%%%%%%%%%%%%

spination=0;

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
E=zeros(NDELTA+1,1); sdis=zeros(NDELTA+1,1);

parfor s=1:spinsnum

      %%%
      [x0(s),y0(s),z0(s),indicator0(s),spheretrackeri,sps] = CreateSpins(FOVx,FOVy,FOVz,R,spheresnum,xcenter,ycenter,zcenter,central,spination,ck,s,scale); %Create Spins according to Spination
      %%%
    
      x = x0(s); y = y0(s); z = z0(s); indicator(s)=indicator0(s); 
    
      %%%%%%%%%%%%%%%%%%%%%%%%%Finally: ITERATIONS
      counter=1; loopsdis=zeros(length(DELTA),1); loopE=zeros(length(DELTA),1);
      for i=0:N
    
                 if dt*i == DELTA(counter)    %collect data at points DELTA to create the attenuation curve later
    
                            loopsdis(counter,1)  =  (x-x0(s))^2;
        
    %                      
                            loopE(counter,1) =  exp(  -1i * q *  (x-x0(s)));
    
                            counter=counter+1;
    
    
                 end
            %Move the spin
            [sps,x,y,z,indicator(s),spheretrackeri] = MovingSpins(kappa,R,dre,dri,dt,De,Di,FOVx,FOVy,FOVz,xcenter,ycenter,zcenter,sps,spheresnum,x,y,z,indicator(s),spheretrackeri); 
    
      end
   sdis=sdis+loopsdis;
      
   E=E+loopE;
   xfin(s)=x;yfin(s)=y;zfin(s)=z;
   s
end

msdis=sdis/spinsnum; E=E/spinsnum; 

rmsdis=sqrt(msdis);



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




figure(3)
plot(sqrt(DELTA(:)),rmsdis(:),'LineWidth',2)
hold on
    
ylabel('rms displacement [\mum]')
xlabel('\surdDiffusion time [ms^1^/^2]')

ax=gca;
ax.FontSize=14;
ax.FontWeight='bold';

for e=1:spheresnum    %Plot spins final locations
    figure (4)
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

figure(4)   %Plot spins final locations
plot3(xfin(fni),yfin(fni),zfin(fni),'.red')
hold on
plot3(xfin(fne),yfin(fne),zfin(fne),'.blue')


 
figure(5)  %Attenuation curve
plot(DELTA(:),log(abs(E(:))),'LineWidth',2)
hold on
    
ylabel('ln(E)')
xlabel('Diffusion time [ms]')
 
save('DELTA_CG_Simulation_g100.mat',"DELTA")
save('E_CG_Simulation_g100.mat',"E") 
 
eenorm=log(abs(E));Ng = length(g); b = (q.^2).*(DELTA) ;


figure(6) %Attenuation curve
for i=1:Ng 
plot(DELTA(1:length(DELTA)),eenorm(1:length(DELTA),i),'o')
hold on
end                     

ylabel('ln(E)')
xlabel('\Delta [ms]')

ax=gca;
ax.FontSize=14;
ax.FontWeight='bold';

% title('x-CG')
legend show

for i = 1:Ng
 legendCell(i) = strcat('q = ',string(num2cell(q(i))),' rad/\mu m');
end
legend(legendCell);

ft = fittype({'x','1'},'coefficients',{'a','b'}); % Define ft as linear fit
coeffnames(ft);

NumOfPoints=length(DELTA); %In each experiment
minus=500;
start=NumOfPoints-minus; % user! Define to get the slope in the long diffusion time regime
finish=NumOfPoints-100; % user!

Tauilim=zeros(Ng,1);
for i=1:Ng
[fitobject,gof] = fit(DELTA(start:finish),eenorm((start:finish),i),ft);
Tauilim(i) = - 1/fitobject.a ;
Fitra=fitobject.b+fitobject.a*DELTA;
plot(DELTA(1:length(DELTA)),Fitra(1:length(DELTA)),'black--')
end


%and now for the bi-exponential fit:

figure(7)
for i=1:Ng 
plot(DELTA,E(:,i),'o')
hold on
end                     

ylabel('E')
xlabel('\Delta [s]')

ax=gca;
ax.FontSize=14;
ax.FontWeight='bold';

% title('x-CG') %% USER!!
legend show
 for i = 1:Ng
 legendCell(i) = strcat('q = ',string(num2cell(q(i))),' rad/m');
 end
legend(legendCell);

%%%
 fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[1,inf,inf],...
               'StartPoint',[0.5 1 1]);
ft2 = fittype('fA*exp(-dA*x*1e-1)+(1-fA)*exp(-dB*x*1e-0)','coefficient',{'fA','dA','dB'},'options',fo);
coeffnames(ft2);

for i=1:Ng
[fitobject,gof] =fit(b(:,i),E(:,i),ft2);

fA(i)=fitobject.fA;
dA(i)=fitobject.dA*1e-1;
dB(i)=fitobject.dB*1e-0;

Rsquare2(i)=gof.rsquare;

Fitra2(:,i)=fA(i)*exp(-dA(i)*b(:,i))+(1-fA(i))*exp(-dB(i)*b(:,i));
plot(DELTA(1:length(DELTA)),Fitra2(:,i),'-black')

ax=gca;
ax.FontSize=14;
ax.FontWeight='bold';

end
 % Aproximation from the bi-exponential fit, in the limit of slow exchange
for i=1:Ng
tau_in_approx(i)=1/(dA(i)*q(i)^2);
end
 
Tauilim
tau_in_approx