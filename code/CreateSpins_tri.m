function[x0,y0,z0,indicator0,spheretrackeri,sps]=CreateSpins_tri(FOVx,FOVx_ExtraComp,FOVy,FOVz,R,spheresnum,xcenter,ycenter,zcenter,central,spination,ck,s,side)
   
    sps=zeros(spheresnum,1);

    if spination==0       
       x0 = side*FOVx_ExtraComp*(0.5 - rand())+(FOVx_ExtraComp-FOVx)/2; %for theta, sin and cos are opossite here, since theta is from -pi/2 to pi/2 and not from 0 to pi as usual.
       y0 = side*FOVy*(0.5 - rand());
       z0 = side*FOVz*(0.5 - rand());
       distfromsphere=zeros(spheresnum,1);  
%        xi=NaN(spinsnum,1);yi=NaN(spinsnum,1);zi=NaN(spinsnum,1);ri=NaN(spinsnum,1);    
              %complicated sequence in order to find NEAREST sphere center, so it
       %will be asked if the spin is inside or outside and at which sphere
            
                for e=1:spheresnum   %for each spin, distance from ALL sphere centers
                    distfromsphere(e) = sqrt((x0-xcenter(e)).^2 + (y0-ycenter(e)).^2 + (z0-zcenter(e)).^2);
                end

                [~,closestsphere] = min(distfromsphere);
%                 mindistfromsphere = min(distfromsphere).'; %find the min distance for each spin
%                 closestsphere = find(temp_distfromsphere==mindistfromsphere);  %the index of the closest sphere
                               
               

                    if distfromsphere(closestsphere) > R
                       
                       indicator0=2;
                       spheretrackeri=NaN;
                    elseif distfromsphere(closestsphere) < R
                       
                       indicator0=1;                                                                  
                       spheretrackeri=closestsphere;
                    end
                       
    elseif spination==1
        %num of spins PER SPHERE. It's a vector, since number of spins in each sphere is allowed to change indepently.
       
     %initial locations , note the use of a matrix for computational ease, although the number of spins in each sphere might be different later on, so for a sphere with less spins - the excess cells will be ignored using the data in the sps vector.
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
       
%         xi=NaN(spinsnum,1);yi=NaN(spinsnum,1);zi=NaN(spinsnum,1);ri=NaN(spinsnum,1);    
%         spheretracker=NaN(1,spinsnum);
       %complicated sequence in order to find NEAREST sphere center, so it
       %will be asked if the spin is inside or outside and at which sphere
                  
                
        phii = 2*pi*rand();         %for each sphere, sps spins are generated randomly in an imagined sphere centered at the origin, then moved to the real sphere using 'x,y,z-center'
        thetai = asin(2*rand()-1);
        temp_ri=R*((rand()).^(1/3));
        
        x0 =  temp_ri.*cos(thetai).*cos(phii) + xcenter(central) ; %for theta, sin and cos are opossite here, since theta is from -pi/2 to pi/2 and not from 0 to pi as usual.
        y0 =  temp_ri.*cos(thetai).*sin(phii) + ycenter(central) ;
        z0 =  temp_ri.*sin(thetai)            + zcenter(central) ;
        spheretrackeri = central;         
        
        indicator0=1;
       
        sps(central)=sps(central)+1;
    
    end  
    

    
    
    