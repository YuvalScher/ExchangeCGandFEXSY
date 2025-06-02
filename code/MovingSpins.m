function [sps,x,y,z,indicator,spheretrackeri] = MovingSpins(kappa,R,dre,dri,dt,De,Di,FOVx,FOVy,FOVz,xcenter,ycenter,zcenter,sps,spheresnum,x,y,z,indicator,spheretrackeri)
skip=0;
        %Infinitesimal Gaussian displacement in a random direction
        phi = 2*pi*rand();            %Random direction for spins (uniform distribution of a sphere SURFACE)
        theta = acos(2*rand()-1);
        
        if indicator==2
            dx=dre*sin(theta).*cos(phi);             %projections
            dy=dre*sin(theta).*sin(phi);
            dz=dre*cos(theta);
        elseif indicator==1
            dx=dri*sin(theta).*cos(phi);             %projections
            dy=dri*sin(theta).*sin(phi);
            dz=dri*cos(theta);
        end

if indicator==1  %If the spin is in the intra-cellular domain        
   
    if  sqrt(  (x+dx-xcenter(spheretrackeri))^2 + (y+dy-ycenter(spheretrackeri))^2 +(z+dz-zcenter(spheretrackeri))^2  ) < R     
    x=x+dx;y=y+dy;z=z+dz;
    else   %If the spin crossed from INTRA to EXTRA
                
                dsi = R - sqrt(  (x-xcenter(spheretrackeri))^2 + (y-ycenter(spheretrackeri))^2 +(z-zcenter(spheretrackeri))^2  );
                dti = (dsi^2)/(6*Di);   %The time that takes the particle, that is about to cross, to reach the membrane. Remember that dri=sqrt(6*Di*dt);
                dte = dt-dti;           %The time that left for diffusion outside the cell


                Pi = 2*dsi*kappa/Di;    %Crossing probability

                        if Pi<0 || Pi>1
                             print('error line 74')
                             error(1)
                        end

                dice = rand();
                     
                if dice < Pi             %Crossing case   
                        
                      if De~=Di
                           dx = sqrt(6*Di*dti + 6*De*dte) * sin(theta) * cos(phi) ; %If de=di, then dxi=dxicross. Otherwise, this is not a futile approximation
                           dy = sqrt(6*Di*dti + 6*De*dte) * sin(theta) * sin(phi) ;
                           dz = sqrt(6*Di*dti + 6*De*dte) * cos(theta) ;
                      end
                           

                     if (   abs( x + dx ) < FOVx/2   &&   abs( y + dy ) < FOVy/2  &&   abs( z + dz ) <FOVz/2  ) 
                            
                         x=x+dx;y=y+dy;z=z+dz;
                         distfromsphere=zeros(spheresnum,1);
                            for e=1:spheresnum   %For each spin, the distance from ALL sphere centers
                                distfromsphere(e) = sqrt((x-xcenter(e)).^2 + (y-ycenter(e)).^2 + (z-zcenter(e)).^2);
                            end

                            [~,closestsphere] = min(distfromsphere);
                            

                            if distfromsphere(closestsphere) > R
                                       moon=1;
                            elseif distfromsphere(closestsphere) < R
                                       moon=2;
                            end
                                    
                            if moon==1
                                %Add the spin to Extra
                                    indicator=2;
                                    sps(spheretrackeri)=sps(spheretrackeri)-1;
                                    skip=1;                                  
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
%If the spin is in the extra-cellular domain (if it just crossed in this iteration than skip=1 and we skip this part)       

if indicator==2 && skip==0
  
     %If a spin plots to get out of the defined FOV -- reflect it back in
            if  (  (abs(x+dx) > FOVx/2) || (abs(y+dy) > FOVy/2)  ||  (abs(z+dz) > FOVz/2)        )            
                dx=0; dy=0; dz=0;
            end
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %Extra --> Intra. Dealing with the probability of membrane crossing.
       
            %Find the NEAREST sphere center, so later it will be asked if the spin crosses this sphere membrane or not:
      
           for e=1:spheresnum   %for each spin, distance from ALL sphere centers
                    distfromsphere(e) = sqrt((x+dx-xcenter(e)).^2 + (y+dy-ycenter(e)).^2 + (z+dz-zcenter(e)).^2);
           end

           [~,closestsphere] = min(distfromsphere);           
                                                        
                           
            if  sqrt(  (x+dx-xcenter(closestsphere))^2 + (y+dy-ycenter(closestsphere))^2 +(z+dz-zcenter(closestsphere))^2  )     >  R
                  x=x+dx;y=y+dy;z=z+dz;
            else
                dse =    sqrt(  (x-xcenter(closestsphere))^2 + (y-ycenter(closestsphere))^2 +(z-zcenter(closestsphere))^2  )  - R;
                      
                    
                dte = (dse^2)/(6*De);   %The time that takes the particle, that is about to cross, to reach the membrane.Remember that dri=sqrt(6*Di*dt)
                dti = dt-dte;           %The time that left for diffusion outside the cell          
                Pe = 2*dse*kappa/De;    %Crossing probability


                        if Pe<0 || Pe>1
                             print('error line 292')
                              error(1)
                        end
                     
               
                dice = rand();
                if dice < Pe                                  %Crossing case:
                      indicator=1;

                       if De~=Di
                          dx = sqrt(6*De*dte + 6*Di*dti) * sin(theta) * cos(phi); %if de=di, then dxi=dxicross. otherwise, this is not a futile calulation!
                          dy = sqrt(6*De*dte + 6*Di*dti) * sin(theta) * sin(phi);
                          dz = sqrt(6*De*dte + 6*Di*dti) * cos(theta);
                       end
             

                      x=x+dx;y=y+dy;z=z+dz;                                                                                             
       
                      spheretrackeri=closestsphere;  
                      sps(closestsphere)=sps(closestsphere) + 1;

              
              elseif  dice > Pe           %reflection   
%                              dx = (sqrt(6*De*dte) - sqrt(6*De*dti))* sin(theta) * cos(phi); 
%                              dy = (sqrt(6*De*dte) - sqrt(6*De*dti))* sin(theta) * sin(phi);
%                              dz = (sqrt(6*De*dte) - sqrt(6*De*dti))* cos(theta);           
%                              if (  (abs(x+dx) > FOVx/2) || (abs(ye+dye) > FOVy/2)  ||  (abs(ze+dze) > FOVz/2)     )              
                             dx=0; dy=0; dz=0;   
%                              end
                             
                             x=x+dx;y=y+dy;z=z+dz;

                end     
            end
            
                                       
end
                

        