function [sps,x,y,z,indicator,spheretrackeri] = MovingSpins_tri(kappa,R,dre,dri,dt,De,Di,FOVx,FOVx_ExtraComp,FOVy,FOVz,xcenter,ycenter,zcenter,sps,spheresnum,x,y,z,indicator,spheretrackeri)
skip=0;
if abs(x)<=FOVx/2
        %infitisimal gaussian displacement in the random direction
        phi = 2*pi*rand();            %random direction for spins (uniform distribution of a sphere SURFACE)
        theta = acos(2*rand()-1);
        
%         dsi=zeros(spinsnumi,1);Pi=zeros(spinsnumi,1);
%         dse=zeros(spinsnume,1);Pe=zeros(spinsnume,1);
        if indicator==2
            dx=dre*sin(theta).*cos(phi);             %projections
            dy=dre*sin(theta).*sin(phi);
            dz=dre*cos(theta);
        elseif indicator==1
            dx=dri*sin(theta).*cos(phi);             %projections
            dy=dri*sin(theta).*sin(phi);
            dz=dri*cos(theta);
        end
%         rnxe=randn([spinsnume 1]); 
%         rnye=randn([spinsnume 1]);
%         rnze=randn([spinsnume 1]);        
%         
%         dxe=sqrt(2*De*dt)*rnxe; 
%         dye=sqrt(2*De*dt)*rnye;
%         dze=sqrt(2*De*dt)*rnze;
% 
%        xi=x0i; yi=y0i; zi=z0i; ri=r0i;
%        xe=x0e; ye=y0e; ze=z0e; re=r0e;
%Intra --> Extra. Dealing with the probability of membrane crossing and making sure that all spins stay inside the defined FOV after crossing
                     
if indicator==1          
   
    if  sqrt(  (x+dx-xcenter(spheretrackeri))^2 + (y+dy-ycenter(spheretrackeri))^2 +(z+dz-zcenter(spheretrackeri))^2  ) < R     
    x=x+dx;y=y+dy;z=z+dz;
    else   %if a spin crossed from INTRA to EXTRA
                
                dsi = R - sqrt(  (x-xcenter(spheretrackeri))^2 + (y-ycenter(spheretrackeri))^2 +(z-zcenter(spheretrackeri))^2  );
                dti = (dsi^2)/(6*Di);   %the time that takes the particle that is about to cross to reach the membrane. remember that dri=sqrt(6*Di*dt);
                dte = dt-dti;                 %the time that left for diffusion outside the cell


                Pi = 2*dsi*kappa/Di;         %crossing probability
%                   Pi = 2/3*dti*kappa/Di;
%                     for i=1:length(k)
%                         if Pi(i)<0 || Pi(i)>1
%                              print('error line 74')
%                              error(1)
%                         end
%                     end    

                dice = rand();
                     
                if dice < Pi    %crossing case   
                        
                      if De~=Di
                           dx = sqrt(6*Di*dti + 6*De*dte) * sin(theta) * cos(phi) ; %if de=di, then dxi=dxicross. otherwise, this is not a futile calulation!
                           dy = sqrt(6*Di*dti + 6*De*dte) * sin(theta) * sin(phi) ;
                           dz = sqrt(6*Di*dti + 6*De*dte) * cos(theta) ;
                      end
                           

                     if (   abs( x + dx ) < FOVx/2   &&   abs( y + dy ) < FOVy/2  &&   abs( z + dz ) <FOVz/2  ) 
                            
                         x=x+dx;y=y+dy;z=z+dz;
                         distfromsphere=zeros(spheresnum,1);
                            for e=1:spheresnum   %for each spin, distance from ALL sphere centers
                                distfromsphere(e) = sqrt((x-xcenter(e)).^2 + (y-ycenter(e)).^2 + (z-zcenter(e)).^2);
                            end

                            [~,closestsphere] = min(distfromsphere);
                            
                            
                            
%                             mindistfromsphere = min(distfromsphere).'; %find the min distance for each spin
%                             closestsphere = find(temp_distfromsphere==mindistfromsphere);  %the index of the closest sphere


                            if distfromsphere(closestsphere) > R
                                       moon=1;
                            elseif distfromsphere(closestsphere) < R
                                       moon=2;
                            end

                                    
                            if moon==1
                                    %%%count
%                                     spinsnume=spinsnume+1; spinsnumi=spinsnumi-1;
                                    indicator=2;
                                    sps(spheretrackeri)=sps(spheretrackeri)-1;
                                    %Add the spin to Extra
                                    skip=1;                                  

    %                                     for e=1:spheresnum   %for each spin, distance from ALL sphere centers
    %                                         distfromsphere(spinsnume,e) = sqrt((x(spinsnume)-xcenter(e)).^2 + (y(spinsnume)-ycenter(e)).^2 + (z(spinsnume)-zcenter(e)).^2);
    %                                     end
    %                                     mindistfromsphere(spinsnume,1) = min(distfromsphere(spinsnume),[],2).'; %find the min distance for each spin
                            elseif moon==2
                                sps(spheretrackeri)=sps(spheretrackeri)-1;
                                spheretrackeri=closestsphere;
                                sps(spheretrackeri)=sps(spheretrackeri)+1;   
                            end    

                    else %reflection                          
%                                     dxi(k(j))=(sqrt(6*Di*dti(j)) - sqrt(6*Di*dte(j))) * sin(thetai(k(j))) * cos(phii(k(j))) ;  
%                                     dyi(k(j))=(sqrt(6*Di*dti(j)) - sqrt(6*Di*dte(j))) * sin(thetai(k(j))) * sin(phii(k(j))) ;
%                                     dzi(k(j))=(sqrt(6*Di*dti(j)) - sqrt(6*Di*dte(j))) * cos(thetai(k(j))) ;     

                                    dx=0;  dy=0;  dz=0;

                                    x=x+dx;y=y+dy;z=z+dz;


                     end    

                else %reflection                          
%                             dxi(k(j))=(sqrt(6*Di*dti(j)) - sqrt(6*Di*dte(j))) * sin(thetai(k(j))) * cos(phii(k(j))) ;  
%                             dyi(k(j))=(sqrt(6*Di*dti(j)) - sqrt(6*Di*dte(j))) * sin(thetai(k(j))) * sin(phii(k(j))) ;
%                             dzi(k(j))=(sqrt(6*Di*dti(j)) - sqrt(6*Di*dte(j))) * cos(thetai(k(j))) ;     
                                    dx=0;  dy=0;  dz=0;

                                    x=x+dx;y=y+dy;z=z+dz;
                end    
                     
                     
    end  
                
        
 end
       
         

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %Extra

if indicator==2 && skip==0
  
     %if a spin plots to get out of the defined FOV -- reflect it back in
            if  (  (abs(x+dx) > FOVx/2) || (abs(y+dy) > FOVy/2)  ||  (abs(z+dz) > FOVz/2)        )            
                dx=0; dy=0; dz=0;
            end
      
%             for j=1:length(f)                  
%                     if xe(f(j))+dxe(f(j)) > FOVx/2
%                     xe(f(j)) = FOVx/2 - (xe(f(j))+dxe(f(j))-FOVx/2);                        
%                     elseif xe(f(j))+dxe(f(j)) < -FOVx/2
%                     xe(f(j)) = -FOVx/2 + (-FOVx/2-xe(f(j))-dxe(f(j)));                          
%                     end 
%                     
%                     if ye(f(j))+dye(f(j)) > FOVy/2
%                     ye(f(j)) = FOVy/2 - (ye(f(j))+dye(f(j))-FOVy/2);     
%                     elseif ye(f(j))+dye(f(j)) < -FOVy/2
%                     ye(f(j)) = -FOVy/2 + (-FOVy/2-ye(f(j))-dye(f(j)));                          
%                     end 
%                     
%                     if ze(f(j))+dze(f(j)) > FOVz/2
%                     ze(f(j)) = FOVz/2 - (ze(f(j))+dze(f(j))-FOVz/2);     
%                     elseif ze(f(j))+dze(f(j)) < -FOVz/2
%                     ze(f(j)) = -FOVz/2 + (-FOVz/2-ze(f(j))-dze(f(j)));                          
%                     end                  
%             end
%                 dse(h) =  re(h) - R;
%                 
%                 Pe(h) = 2*dse(h)*kappa/De;                        %crossing probability
%                 dte = (dse(h).^2)./(6.*Di);   %the time that takes the particle that is about to cross to reach the membrane. remember that dri=sqrt(6*Di*dt);
%                 dti = dt-dte;                 %the time that left for diffusion outside the cell          
                  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %Extra --> Intra. Dealing with the probability of membrane crossing.
       
            %complicated sequence in order to find NEAREST sphere center, so later it will be asked if the spin crosses this sphere membrane or not:

            
           for e=1:spheresnum   %for each spin, distance from ALL sphere centers
                    distfromsphere(e) = sqrt((x+dx-xcenter(e)).^2 + (y+dy-ycenter(e)).^2 + (z+dz-zcenter(e)).^2);
           end

           [~,closestsphere] = min(distfromsphere);           
                                                        
                           
            if  sqrt(  (x+dx-xcenter(closestsphere))^2 + (y+dy-ycenter(closestsphere))^2 +(z+dz-zcenter(closestsphere))^2  )     >  R
                  x=x+dx;y=y+dy;z=z+dz;
            else
                dse =    sqrt(  (x-xcenter(closestsphere))^2 + (y-ycenter(closestsphere))^2 +(z-zcenter(closestsphere))^2  )  - R;
                      
                    
                dte = (dse^2)/(6*De);   %the time that takes the particle that is about to cross to reach the membrane. remember that dri=sqrt(6*Di*dt);
                dti = dt-dte;                 %the time that left for diffusion outside the cell          
                Pe = 2*dse*kappa/De;                        %crossing probability


                        if Pe<0 || Pe>1
                             print('error line 292')
                              error(1)
                        end
                     
               
                dice = rand();
                if dice < Pe                                  %crossing case:
                      indicator=1;
%                       spinsnumi=spinsnumi+1; spinsnumi=spinsnumi-1;



                       if De~=Di
                          dx = sqrt(6*De*dte + 6*Di*dti) * sin(theta) * cos(phi); %if de=di, then dxi=dxicross. otherwise, this is not a futile calulation!
                          dy = sqrt(6*De*dte + 6*Di*dti) * sin(theta) * sin(phi);
                          dz = sqrt(6*De*dte + 6*Di*dti) * cos(theta);
                       end
             

                      x=x+dx;y=y+dy;z=z+dz;                                                 
                            
%                             for e=1:spheresnum   %for the spin, distance from ALL sphere centers
%                                  edistfromsphere(e,1) = sqrt(( xe(h(g)) + dxcross -xcenter(e)).^2 + (ye(h(g)) + dycross-ycenter(e)).^2 + (ze(h(g)) + dzcross-zcenter(e)).^2);
%                             end
% 
%                             emindistfromsphere = min(edistfromsphere).'; %find the min distance for the spin
% 
%                                                 
%                             etemp_distfromsphere(:,1) =  edistfromsphere(:,1) - emindistfromsphere*ones(spheresnum,1); %this way, the coloumn where the min distance is 0!
%                             eclosestsphere = find(etemp_distfromsphere(:)==0);  %the number of the coloumn with 0 in it is the number of the sphere!
%                                                     
%                             if  eclosestsphere==closestsphere(h(g))        
%                              spheretrackeri(spinsnumi,1)=closestsphere(h(g));  
%                              sps(closestsphere(h(g)))=sps(closestsphere(h(g))) + 1;
%                             elseif eclosestsphere~=closestsphere(h(g))  
%                              spheretrackeri(spinsnumi,1)=eclosestsphere;  
%                              sps(eclosestsphere)=sps(eclosestsphere) + 1;
%                             end    
%                             
       
                      spheretrackeri=closestsphere;  
                      sps(closestsphere)=sps(closestsphere) + 1;

              
              elseif  dice > Pe           %reflection   
%                              dxe(h(g)) = (sqrt(6*De*dte(g)) - sqrt(6*De*dti(g)))* sin(thetae(h(g))) * cos(phie(h(g))); 
%                              dye(h(g)) = (sqrt(6*De*dte(g)) - sqrt(6*De*dti(g)))* sin(thetae(h(g))) * sin(phie(h(g)));
%                              dze(h(g)) = (sqrt(6*De*dte(g)) - sqrt(6*De*dti(g)))* cos(thetae(h(g)));           
%                              if (  (abs(xe(h(g))+dxe(h(g))) > FOVx/2) || (abs(ye(h(g))+dye(h(g))) > FOVy/2)  ||  (abs(ze(h(g))+dze(h(g))) > FOVz/2)     )              
                             dx=0; dy=0; dz=0;   
%                              end
                             
                             x=x+dx;y=y+dy;z=z+dz;
%                                  if (~ismember(h(g),f)) 
%                                         dxe(h(g)) = (sqrt(6*De*dte(g)) - sqrt(6*De*dti(g)))* sin(thetae(h(g))) * cos(phie(h(g))); 
%                                         dye(h(g)) = (sqrt(6*De*dte(g)) - sqrt(6*De*dti(g)))* sin(thetae(h(g))) * sin(phie(h(g)));
%                                         dze(h(g)) = (sqrt(6*De*dte(g)) - sqrt(6*De*dti(g)))* cos(thetae(h(g)));                
%                                 else
%                                         xe(h(g)) = xe(h(g))+sqrt(2*De*dte(g));              %first reach FOV wall
%                                         ye(h(g)) = ye(h(g))+sqrt(2*De*dte(g));
%                                         ze(h(g)) = ze(h(g))+sqrt(2*De*dte(g));    
%                                         re(h(g))= sqrt( xe(h(g)).^2 + ye(h(g)).^2 + ze(h(g)).^2 );
%                                         dse(h) =  re(h) - R;
%                                 end
                end     
            end
            
                                       
end
elseif abs(x)>FOVx/2              
        phi = 2*pi*rand();            %random direction for spins (uniform distribution of a sphere SURFACE)
        theta = acos(2*rand()-1);

        dx=dri*sin(theta).*cos(phi);             %projections
        dy=dri*sin(theta).*sin(phi);
        dz=dri*cos(theta);

       if ( FOVx < abs( x + dx ) && abs( x + dx ) < FOVx_ExtraComp   &&  abs( y + dy ) < FOVx_ExtraComp  &&   abs( z + dz ) <FOVx_ExtraComp  ) 
                         
          x=x+dx;y=y+dy;z=z+dz;
       
       end

end       