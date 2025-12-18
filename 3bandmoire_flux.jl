using ITensors
function threebandmoire_flux(eA,t1,t2,t3,phi,Nx,Ny,V1,V2,flux)
    os = OpSum()
    # A site occupation
    for i = 1:Nx
        for j = 1:Ny
            n = (i-1)*Ny+j
            if j%6 == 2 || j%6 == 5
                os += eA,"N",n
            end
        end
    end
    # C to A (A to C) hopping
    for i = 1:Nx
        for j = 1:Ny
            n = (i-1)*Ny+j
            if j%6 == 1 && j != 1
                if i != Nx
                    os += t2,"Cdag",(n+Ny+1),"C",n
                    os += t2,"Cdag",n,"C",(n+Ny+1)
                    os += V1,"N",n,"N",(n+Ny+1)
                    # phi = pi
                end
                os += -t2*exp(1im*phi),"Cdag",(n+1),"C",n
                os += -t2*exp(-1im*phi),"Cdag",n,"C",(n+1)
                os += V1,"N",n,"N",(n+1)
                #phi = pi/3
                os += -t2*exp(-1im*phi),"Cdag",(n-2),"C",n
                os += -t2*exp(1im*phi),"Cdag",n,"C",(n-2)
                os += V1,"N",n,"N",(n-2)
                #phi = -pi/3
            end
            if j%6 == 1 && j == 1
                if i != Nx
                    os += t2,"Cdag",(n+Ny+1),"C",n
                    os += t2,"Cdag",n,"C",(n+Ny+1)
                    os += V1,"N",n,"N",(n+Ny+1)
                    # phi = pi
                end
                os += -t2*exp(1im*phi),"Cdag",(n+1),"C",n
                os += -t2*exp(-1im*phi),"Cdag",n,"C",(n+1)
                os += V1,"N",n,"N",(n+1)
                #phi = pi/3
                os += -t2*exp(-1im*phi)*exp(1im*flux),"Cdag",(n+Ny-2),"C",n
                os += -t2*exp(1im*phi)*exp(-1im*flux),"Cdag",n,"C",(n+Ny-2)
                os += V1,"N",n,"N",(n+Ny-2)
                #phi = -pi/3
            end
            if j%6 == 4
                if i != 1
                    os += -t2*exp(1im*phi),"Cdag",(n-Ny+1),"C",n
                    os += -t2*exp(-1im*phi),"Cdag",n,"C",(n-Ny+1)
                    os += V1,"N",n,"N",(n-Ny+1)
                #phi = pi/3
                end
                os += -t2*exp(-1im*phi),"Cdag",(n-2),"C",n
                os += -t2*exp(1im*phi),"Cdag",n,"C",(n-2)
                os += V1,"N",n,"N",(n-2)
                #phi = -pi/3
                os += t2,"Cdag",(n+1),"C",n
                os += t2,"Cdag",n,"C",(n+1)
                os += V1,"N",n,"N",(n+1)
                #phi = pi
            end
        end
    end
    # B to A (A to B) hopping
    for i = 1:Nx
        for j = 1:Ny
            n = (i-1)*Ny+j
            if j%6 == 3 
                if i != Nx
                    os += -t1*exp(-1im*2*phi),"Cdag",(n+Ny-1),"C",n
                    os += -t1*exp(1im*2*phi),"Cdag",n,"C",(n+Ny-1)
                    os += V1,"N",n,"N",(n+Ny-1)
                #phi = -2*pi/3
                end
                os += -t1,"Cdag",(n-1),"C",n
                os += -t1,"Cdag",n,"C",(n-1)
                os += V1,"N",n,"N",(n-1)
                #phi = 0
                os += -t1*exp(1im*2*phi),"Cdag",(n+2),"C",n
                os += -t1*exp(-1im*2*phi),"Cdag",n,"C",(n+2)
                os += V1,"N",n,"N",(n+2)
                #phi = 2*pi/3
            end
            if j%6 == 0 && j != Ny
                if i != 1
                    os += -t1,"Cdag",(n-Ny-1),"C",n
                    os += -t1,"Cdag",n,"C",(n-Ny-1)
                    os += V1,"N",n,"N",(n-Ny-1)
                #phi = 0
                end
                os += -t1*exp(-1im*2*phi),"Cdag",(n-1),"C",n
                os += -t1*exp(1im*2*phi),"Cdag",n,"C",(n-1)
                os += V1,"N",n,"N",(n-1)
                #phi = -2*pi/3
                os += -t1*exp(1im*2*phi),"Cdag",(n+2),"C",n
                os += -t1*exp(-1im*2*phi),"Cdag",n,"C",(n+2)
                os += V1,"N",n,"N",(n+2)
                #phi = 2*pi/3
            end
            if j%6 == 0 && j == Ny
                if i != 1
                    os += -t1,"Cdag",(n-Ny-1),"C",n
                    os += -t1,"Cdag",n,"C",(n-Ny-1)
                    os += V1,"N",n,"N",(n-Ny-1)
                #phi = 0
                end
                os += -t1*exp(-1im*2*phi),"Cdag",(n-1),"C",n
                os += -t1*exp(1im*2*phi),"Cdag",n,"C",(n-1)
                os += V1,"N",n,"N",(n-1)
                #phi = -2*pi/3
                os += -t1*exp(1im*2*phi)*exp(-1im*flux),"Cdag",(n+2-Ny),"C",n
                os += -t1*exp(-1im*2*phi)*exp(1im*flux),"Cdag",n,"C",(n+2-Ny)
                os += V1,"N",n,"N",(n+2-Ny)
                #phi = 2*pi/3
            end
        end
    end
    # B to C (C to B) hopping
    for i = 1:Nx
        for j = 1:Ny
            n = (i-1)*Ny+j
            if j%6 == 3 
                if i != Nx
                    os += -t3,"Cdag",(n+Ny+1),"C",n
                    os += -t3,"Cdag",n,"C",(n+Ny+1)
                    os += V2,"N",n,"N",(n+Ny+1)
                end
                os += -t3,"Cdag",(n+1),"C",n
                os += -t3,"Cdag",n,"C",(n+1)
                os += V2,"N",n,"n",(n+1)
                os += -t3,"Cdag",(n-2),"C",n
                os += -t3,"Cdag",n,"C",(n-2)
                os += V2,"N",n,"N",(n-2)
            end
            if j%6 == 0 && j != Ny
                if i != 1
                    os += -t3,"Cdag",(n-Ny+1),"C",n
                    os += -t3,"Cdag",n,"C",(n-Ny+1)
                    os += V2,"N",n,"N",(n-Ny+1)
                end
                os += -t3,"Cdag",(n-2),"C",n
                os += -t3,"Cdag",n,"C",(n-2)
                os += V2,"N",n,"N",(n-2)
                os += -t3,"Cdag",(n+1),"C",n
                os += -t3,"Cdag",n,"C",(n+1)
                os += V2,"N",n,"N",(n+1)
            end
            if j%6 == 0 && j == Ny
                if i != 1
                    os += -t3*exp(-1im*flux),"Cdag",(n-2*Ny+1),"C",n
                    os += -t3*exp(1im*flux),"Cdag",n,"C",(n-2*Ny+1)
                    os += V2,"N",n,"N",(n-2*Ny+1)
                end
                os += -t3,"Cdag",(n-2),"C",n
                os += -t3,"Cdag",n,"C",(n-2)
                os += V2,"N",n,"N",(n-2)
                os += -t3*exp(-1im*flux),"Cdag",(n-Ny+1),"C",n
                os += -t3*exp(1im*flux),"Cdag",n,"C",(n-Ny+1)
                os += V2,"N",n,"N",(n-Ny+1)
            end
        end
    end
    return os
end



