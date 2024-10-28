function[x0,y0,z0,indicator0,spheretrackeri,sps]=CreateSpins(FOVx,FOVy,FOVz,R,spheresnum,xcenter,ycenter,zcenter,central,spination,ck,s,side)
   
    sps=zeros(spheresnum,1);

    if spination==0       
       x0 = side*FOVx*(0.5 - rand()); %For theta, sin and cos are opossite here, since theta is from -pi/2 to pi/2 and not from 0 to pi as usual.
       y0 = side*FOVy*(0.5 - rand());
       z0 = side*FOVz*(0.5 - rand());
       distfromsphere=zeros(spheresnum,1);  
 
       %Find NEAREST sphere center          
                for e=1:spheresnum   %For each spin, distance from ALL sphere centers
                    distfromsphere(e) = sqrt((x0-xcenter(e)).^2 + (y0-ycenter(e)).^2 + (z0-zcenter(e)).^2);
                end

                [~,closestsphere] = min(distfromsphere);            
                    %Is the spin in the closest sphere or out of it?
                    if distfromsphere(closestsphere) > R    
                       indicator0=2;
                       spheretrackeri=NaN;
                    elseif distfromsphere(closestsphere) < R
                       
                       indicator0=1;                                                                  
                       spheretrackeri=closestsphere;
                    end
                       
    elseif spination==1
    
     %  Initial location
        sps(ck(s))=sps(ck(s))+1;  
                  
        phii = 2*pi*rand();         %for each sphere, sps spins are generated randomly in an imagined sphere centered at the origin, then moved to the real sphere using 'x,y,z-center'
        thetai = asin(2*rand()-1);
        temp_ri=R*((rand()).^(1/3));
        
        x0 =  temp_ri.*cos(thetai).*cos(phii) + xcenter(ck(s)) ; %for theta, sin and cos are opossite here, since theta is from -pi/2 to pi/2 and not from 0 to pi as usual.
        y0 =  temp_ri.*cos(thetai).*sin(phii) + ycenter(ck(s)) ;
        z0 =  temp_ri.*sin(thetai)            + zcenter(ck(s)) ;
        spheretrackeri = ck(s); 
    
        indicator0=1;       
   
    elseif spination==2
                      
        phii = 2*pi*rand();      %for each sphere, sps spins are generated randomly in an imagined sphere centered at the origin, then moved to the real sphere using 'x,y,z-center'
        thetai = asin(2*rand()-1);
        temp_ri=R*((rand()).^(1/3));
        
        x0 =  temp_ri.*cos(thetai).*cos(phii) + xcenter(central) ; %for theta, sin and cos are opossite here, since theta is from -pi/2 to pi/2 and not from 0 to pi as usual.
        y0 =  temp_ri.*cos(thetai).*sin(phii) + ycenter(central) ;
        z0 =  temp_ri.*sin(thetai)            + zcenter(central) ;
        spheretrackeri = central;         
        
        indicator0=1;
       
        sps(central)=sps(central)+1;
    
    end  
    

    
    
    