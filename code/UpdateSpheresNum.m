function[FOVx, FOVy, FOVz, spheresnum]=UpdateSpheresNum(a,FOVx,FOVy,FOVz)

    FOVx = FOVx - mod(FOVx, a);
    FOVy = FOVy - mod(FOVy, a);
    FOVz = FOVz - mod(FOVz, a);
    
    R=a/2;
    sphere_volume = (4/3) * pi * R^3;
    FOV_volume = FOVx * FOVy * FOVz;

    spheresnum=(FOVx)*(FOVy)*(FOVz)/(a^3);

    alpha = (sphere_volume / FOV_volume) * spheresnum;
    
