
clear all

%User defined tissue model parameters:

dt = 10^-1;                                   %Millisecond, simulation time step
spinsnum = 10^7;                           %Number of spins

FOVx = 5.5*3;                              %Micrometers. In case of closed packing these numbers will be changed to nearest values that can exactly accomodate the packing
FOVy = 5.5*3;
FOVz = 5.5*3;

a = 5;                                    %Spheric compartment diameter in micrometers
R = a/2;

De = 1;                                    %Extracellular diffusion coef. micrometer^2/millisec
Di = 1;                                    %Intracellular diffusion coef. micrometer^2/millisec

ki = 0;                                 %The "real" exchange rate from intra to extra
kappa = (R*ki)/3;                             %Micrometer/millisecond. Permeability.

spheresnum = 1;   

packing=0;

%User defined experiment parameters:
%%%%%%%%% FEXSY %%%%%%%%%%%%%%%%%%                        % filter block gradient pulse strengh g/cm
g_user = linspace(1,500,60);            % measurement block gradient pulse strengh g/cm
d=2;                     % gradient pulse duration ms
% Tm=[8, 25, 50, 100, 150, 200, 250, 300, 400, 500, 750, 1000];          %ms

tau=100;                  %ms
%%%%
g = g_user*(10^(-8));    %T/micrometer 
gyro= 42.57747 * 10^3;   %gyromagentic ratio of protons 1/(ms T)
q = 2 * pi * gyro * g * d;          %q vector


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AUTOMATIC FROM HERE

if packing ~= 0   %fix FOV size such that packing 1 is exactly fitting in the box. The extra volume for packing 2-3 will be erased later
[FOVx, FOVy, FOVz, spheresnum] = UpdateSpheresNum(a,FOVx,FOVy,FOVz);    
end

                         
dre=sqrt(6*De*dt);
dri=sqrt(6*Di*dt);

[x,y,z] = sphere; 

[FOVx, FOVy, FOVz, FOV_volume, alpha, spheresnum, xcenter, ycenter, zcenter, rcenter] = CreateSpheres(a,FOVx,FOVy,FOVz,spheresnum,packing); %Create Spheres according the Packing

for e=1:spheresnum      
    figure (1)
    surf(R*x+xcenter(e),R*y+ycenter(e),R*z+zcenter(e),'FaceColor', 'none')
    hold on
end
xlim([-FOVx/2 FOVx/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])


for e=1:spheresnum    
    figure (2)
    surf(R*x+xcenter(e),R*y+ycenter(e),R*z+zcenter(e))
    hold on
end
xlim([-FOVx/2 FOVx/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])

global x0i; global y0i; global z0i; global r0i; global xi; global yi; global zi; global ri; global xe; global ye; global ze; global re; global x0e; global y0e; global z0e; global r0e; global spheretrackeri;
global dephasei; global dephasee;


spination=2;

%%%%
[sspss,sps,spinsnum,spinsnumi,spinsnume] = CreateSpins(FOVx,FOVy,FOVz,R,spheresnum,xcenter,ycenter,zcenter,spinsnum,spination); %Create Spins according the Spination
%%%
dephasei=zeros(spinsnumi,1); dephasee = zeros(spinsnume,1);



figure(1)
plot3(x0i,y0i,z0i,'.red',x0e,y0e,z0e,'.blue')
hold on



    %more spins counters :(
    
% fi=zeros(NDELTA,1);
% fe=zeros(NDELTA,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%Finally: ITERATIONS
i=0;       
while dt*i < tau+dt     %collect rms data at points DELTA
% 
%                 xrmsdis_i  =  (xi-x0i).^2;
%                 %   xrmsdis_i =  (xi-x0i).^2 + (yi-y0i).^2 + (zi-z0i).^2 ;
% 
%                 xrmsdis_e =  (xe-x0e).^2;
%                 %   xrmsdis_e =  (xe-x0e).^2 + (ye-y0e).^2 + (ze-z0e).^2 ;
% 
%                 xmsdis(counter,1) = ( sum(xrmsdis_i) +  sum(xrmsdis_e) )/ spinsnum ;
%                 rxmsdis(counter,1) = sqrt(xmsdis(counter));
% 
%                 E(counter,1) = (  sum(exp(  1i * q *  (xi-x0i) ) )   + sum(exp(  1i * q *  (xe-x0e) ) )   )   / spinsnum ;
% 
%                 fe(counter,1)=spinsnume/(spinsnume+spinsnumi);
%                 fi(counter,1)=spinsnumi/(spinsnume+spinsnumi);
[sps,spinsnumi,spinsnume] = MovingSpinsFEXSY(kappa,R,dre,dri,dt,De,Di,FOVx,FOVy,FOVz,xcenter,ycenter,zcenter,sps,spinsnumi,spinsnume,spheresnum); 

i=i+1
end



for j=1:length(q)
  E(j) = (  sum( exp( - 1i * q(j)* (xi-x0i) )  )   +  sum( exp( -1i * q(j)* (xe-x0e) )  ) )/ spinsnum ;           
end 


Enorm = E(:)/E(1);  

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% 
% figure(3)
% plot(sqrt(DELTA(:)),rxmsdis(:),'LineWidth',2)
% hold on
%     
% ylabel('rms displacement [\mum]')
% xlabel('\surdDiffusion time [ms^1^/^2]')
% 
% ax=gca;
% ax.FontSize=14;
% ax.FontWeight='bold';

% legend show
% for c = 1:length(kappa)
%     legendCell(c) = strcat('\kappa =',string(num2cell(kappa(c))),' \mum/ms');
% end
% legend(legendCell);
% 
% tau(1)=0; %= infinity
% for c=2:length(kappa)
%     tau(c) = R/(2*kappa(c));
% end
% for c=1:length(tau)
%     figure(6)
%     plot(sqrt(DELTA+dt),rmsdis(:,c),'LineWidth',2)
%     hold on
% end
% ylabel('rms displacement [\mum]')
% % xlabel('\surdDiffusion time [ms^1^/^2]')
% 
% ax=gca;
% ax.FontSize=14;
% ax.FontWeight='bold';
% 
% legend show
% legendCell(1) = strcat('\tau_{i} = \infty');
% for c = 2:length(tau)
%     legendCell(c) = strcat('\tau_{i} = ',string(num2cell(tau(c))),' ms');
% end
% legend(legendCell);

for e=1:spheresnum
    figure (4)
    surf(R*x+xcenter(e),R*y+ycenter(e),R*z+zcenter(e),'FaceColor', 'none')
    hold on
end
xlim([-FOVx/2 FOVx/2])
ylim([-FOVy/2 FOVy/2])
zlim([-FOVz/2 FOVz/2])
for e=1:spheresnum
    figure(4)
    plot3(xi,yi,zi,'.red')
    hold on
end
figure(4)
plot3(xe,ye,ze,'.blue')
 

figure(5)

plot(q*R/2/pi,log((Enorm)),'.k','MarkerSize',20,'LineWidth',3)

ylabel('$E$','Interpreter','latex')
xlabel('$\tilde{q}R$','Interpreter','latex')

hold on


analytical =  (abs(3*(q*R .* cos(q*R) - sin(q*R))./((q*R).^3))).^2;

plot(q*R/2/pi,log(analytical),'-b','MarkerSize',20,'LineWidth',3)

save('checkq.mat',"q")
save('checkE.mat',"E") 
 