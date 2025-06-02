function[FOVx, FOVy, FOVz, FOV_volume, alpha, spheresnum, xcenter, ycenter, zcenter, rcenter]=CreateSpheres(a,FOVx,FOVy,FOVz,spheresnum,packing)

   R=a/2;  sphere_volume = (4/3) * pi * R^3; FOV_volume = FOVx*FOVy*FOVz;

    if packing==0 %then RANDOM PACKING
      
           xcenter(1)=0;
           ycenter(1)=0;
           zcenter(1)=0;
           rcenter(1)=0;
            for e=2:spheresnum
                    logical=true;
                    while(logical)          %While there is no ovelap:
                        logical=false;
                        xcenter(e) = FOVx * (rand() - 0.5);    %Random uniform distributiom in carteasian coordinates.
                        ycenter(e) = FOVy * (rand() - 0.5);
                        zcenter(e) = FOVz * (rand() - 0.5);
                        rcenter(e) = sqrt(xcenter(e)^2 + ycenter(e)^2 + zcenter(e)^2);

                        for l=1:e
                            distcenter(l) = sqrt (    ( xcenter(e)- xcenter(l) )^2 + ( ycenter(e)- ycenter(l) )^2 + ( zcenter(e)- zcenter(l) )^2     );   %chekcing for overlap
                            if ((distcenter(l) < a) && (l~=e) || ((abs(xcenter(e))+R)>FOVx/2) || ((abs(ycenter(e))+R)>FOVy/2) || ((abs(zcenter(e))+R)>FOVz/2))
                                logical=true;
                            end
                        end
                    end             
            end
    elseif packing==1  %Layers are exactly on top of each other
        counter=1;
        for ez=1:(FOVz/a)
            for ey=1:(FOVy/a)
                for ex=1:(FOVx/a)
                    xcenter(counter) = -FOVx/2 + R + a*(ex-1);
                    ycenter(counter) = -FOVy/2 + R + a*(ey-1);
                    zcenter(counter) = -FOVz/2 + R + a*(ez-1);
                    rcenter(counter) = sqrt(xcenter(counter)^2 + ycenter(counter)^2 + zcenter(counter)^2);
                    counter=counter+1;
                end
            end
        end
    elseif packing==2  % HCP lattice
        counter=0;
        for ez=0:(FOVz/a)
            for ey=0:(FOVy/a)
                    for ex=0:(FOVx/a)
                        counter=counter+1;
                        xcenter(counter) = -FOVx/2 - 0.5*R   + a * ex + R * mod(ey+ez,2);
                        ycenter(counter) = -FOVy/2 + sqrt(3) * R * ey + R / 3 * (mod(ez,2));
                        zcenter(counter) = -FOVz/2 + 2 / 3 * sqrt(6) * R * ez;
                        rcenter(counter) = sqrt(xcenter(counter)^2 + ycenter(counter)^2 + zcenter(counter)^2);                      
                    end     
            end
        end
        spheresnum=counter;
        ycenter(:) = ycenter(:)+(FOVy/2-(max(ycenter)))/2;
        zcenter(:) = zcenter(:)+(FOVz/2-(max(zcenter)))/2;
        
        FOVx = (max(abs(xcenter))+R) * 2;
        FOVy = (max(abs(ycenter))+R) * 2;
        FOVz = (max(abs(zcenter))+R) * 2;

    elseif packing==3  % FCC lattice
        counter=0;
        alpha = (sphere_volume/FOV_volume) * spheresnum;
        for ez=1:(FOVz/a)
            if mod(ez,3)==2
                for ey=1:(FOVy/a)-2
                    if mod(ez,3)==0
                        for ex=1:(FOVx/a)
                            counter=counter+1;
                            xcenter(counter) = -FOVx/2 + R+a*(ex-1);
                            ycenter(counter) = -FOVy/2 + sqrt(3)*R+a*(ey-1)-(R*(sqrt(3)-1));
                            zcenter(counter) = -FOVz/2 + 2/3*sqrt(6)*R*(ez-1)+R;
                            rcenter(counter) = sqrt(xcenter(counter)^2 + ycenter(counter)^2 + zcenter(counter)^2);
                            
                        end
                    elseif mod(ez,3)==1
                        for ex=1:(FOVx/a)-1
                            counter=counter+1;
                            xcenter(counter) = -FOVx/2 + R+a*(ex-1) + R;
                            ycenter(counter) = -FOVy/2 + sqrt(3)*R+a*(ey-1)-(R*(sqrt(3)-1));
                            zcenter(counter) = -FOVz/2 +  2/3*sqrt(6)*R*(ez-1)+R;
                            rcenter(counter) = sqrt(xcenter(counter)^2 + ycenter(counter)^2 + zcenter(counter)^2);                          
                        end
                    else
                        for ex=1:(FOVx/a)-1
                            counter=counter+1;
                            xcenter(counter) = -FOVx/2 + R+a*(ex-1) + R;
                            ycenter(counter) = -FOVy/2 + sqrt(3)*R+a*(ey-1)-(R*(sqrt(3)-1))+R;
                            zcenter(counter) = -FOVz/2 +  2/3*sqrt(6)*R*(ez-1)+R;
                            rcenter(counter) = sqrt(xcenter(counter)^2 + ycenter(counter)^2 + zcenter(counter)^2);                            
                        end
                    end
                end
            else
                for ey=1:(FOVy/a)-1
                    if mod(ez,3)==0
                        for ex=1:(FOVx/a)
                            counter=counter+1;
                            xcenter(counter) = -FOVx/2 + R+a*(ex-1);
                            ycenter(counter) = -FOVy/2 + sqrt(3)*R+a*(ey-1)-(R*(sqrt(3)-1));
                            zcenter(counter) = -FOVz/2 + 2/3*sqrt(6)*R*(ez-1)+R;
                            rcenter(counter) = sqrt(xcenter(counter)^2 + ycenter(counter)^2 + zcenter(counter)^2);                          
                        end
                    elseif mod(ez,3)==1
                        for ex=1:(FOVx/a)-1
                            counter=counter+1;
                            xcenter(counter) = -FOVx/2 + R+a*(ex-1) + R;
                            ycenter(counter) = -FOVy/2 + sqrt(3)*R+a*(ey-1)-(R*(sqrt(3)-1));
                            zcenter(counter) = -FOVz/2 +  2/3*sqrt(6)*R*(ez-1)+R;
                            rcenter(counter) = sqrt(xcenter(counter)^2 + ycenter(counter)^2 + zcenter(counter)^2);     
                        end
                    else
                        for ex=1:(FOVx/a)-1
                            counter=counter+1;
                            xcenter(counter) = -FOVx/2 + R+a*(ex-1) + R;
                            ycenter(counter) = -FOVy/2 + sqrt(3)*R+a*(ey-1)-(R*(sqrt(3)-1))-R;
                            zcenter(counter) = -FOVz/2 +  2/3*sqrt(6)*R*(ez-1)+R;
                            rcenter(counter) = sqrt(xcenter(counter)^2 + ycenter(counter)^2 + zcenter(counter)^2); 
                        end
                    end
                end
            end

        end

        spheresnum=counter;
        ycenter(:) = ycenter(:)+(FOVy/2-(R+max(ycenter)))/2;
        zcenter(:) = zcenter(:)+(FOVz/2-(R+max(zcenter)))/2;
        
        FOVx = (max(abs(xcenter))+R) * 2;
        FOVy = (max(abs(ycenter))+R) * 2;
        FOVz = (max(abs(zcenter))+R) * 2;
    end
    
FOV_volume = FOVx*FOVy*FOVz;
alpha = (sphere_volume/FOV_volume) * spheresnum;
xcenter=xcenter';ycenter=ycenter'; zcenter=zcenter'; rcenter=rcenter';