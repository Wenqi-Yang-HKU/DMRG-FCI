using MKL
using ITensors
using ITensors.HDF5
using KrylovKit
include("./3bandmoire_flux.jl")
function dmrg_(H,psi0,nsweep,maxdim,cutoff)
    if length(maxdim) < nsweep
        for i = length(maxdim)+1:nsweep
            push!(maxdim,maxdim[end])
        end
    end
    PH = ProjMPO(H) 
    psi = copy(psi0)
    N = length(psi)
    if !isortho(psi) || orthoCenter(psi) != 1
        psi = ITensors.orthogonalize!(PH,psi,1)
    end
    PH = position!(PH,psi,1)
    energy = 0.0
    for sw in 1:nsweep
        maxtruncerr = 0.0
        for (b,ha) in sweepnext(N)
            PH = position!(PH,psi,b)
            phi = psi[b]*psi[b+1]
            vals,vecs = eigsolve(PH,phi,1,:SR;ishermitian=true,tol=1e-14,krylovdim=3,maxiter=1)
            energy = vals[1]
            phi = vecs[1]
            ortho = ha ==1 ? "left" : "right"
            spec = replacebond!(PH,psi,b,phi;maxdim=maxdim[sw],cutoff=cutoff,ortho=ortho,normalize=true)
            maxtruncerr = max(maxtruncerr,spec.truncerr)
        end
        io = open("result.txt","a+")
        write(io,"After sweep $sw energy=$energy maxlinkdim = $(maxlinkdim(psi)) maxerr = $(maxtruncerr)\n")
        close(io)
    end
    return energy,psi
end

let 
    t1 = 1
    t2 = 1
    t3 = 0
    eA = 10
    phi = pi/3
    Nx = 24
    Ny = 12
    V1 = 1
    V2 = 1
    nsweep_initial = 20
    nsweeps = 6
    #maxdim = Array{Int}(undef,nsweep_initial)
    #for i = 0:nsweep_initial-1
    #    if round(256*2^(0.25*(i÷2))) <= 1024
    #        maxdim[i+1] = round(256*2^(0.25*(i÷2)))
    #    else
    #        maxdim[i+1] = 1024
    #    end
    #end
    cutoff = 1E-8
    f = h5open("psi0.h5","r")
    psi0 = read(f,"psi",MPS);
    normalize!(psi0)
    close(f) 
    sites = siteinds(psi0)
    os = threebandmoire_flux(eA,t1,t2,t3,phi,Nx,Ny,V1,V2,0)
    H = MPO(os,sites)
    energy,psi = dmrg_(H,psi0,nsweep_initial,[256],cutoff)
    occupation_number = expect(psi,"N")
    left_charge = sum(occupation_number[1:(Nx*Ny÷2)])
    io = open("result.txt","a+")
    write(io,"flux = 0 left edge charge = $left_charge\n")
    close(io)
    f = h5open("flux_0.h5","w")
    write(f,"psi",psi)
    close(f)
    for flux in π/8:π/8:6*π
        os = threebandmoire_flux(eA,t1,t2,t3,phi,Nx,Ny,V1,V2,flux)
        H = MPO(os,sites)
        energy,psi = dmrg_(H,psi,nsweeps,[256],cutoff)
        occupation_number = expect(psi,"N")
        left_charge = sum(occupation_number[1:(Nx*Ny÷2)])
        io = open("result.txt","a+")
        write(io,"flux = $flux left edge charge = $left_charge\n")
        close(io)
	    if round(flux/(pi/8)) in 16*(-10:10)
            f = h5open("flux_$(round(flux/(π/8))).h5","w")
            write(f,"psi",psi)
            close(f)
    	end
    end
end
