program structure_factor

    implicit none

    integer, parameter          :: dp  = selected_real_kind(15,307)
    real(kind=dp), parameter    :: pi = 4.0_dp*atan(1.0_dp)

    integer                     :: j1, j2, nx, ny, nz, ki, rk, np, nxi, nyi
    real(kind=dp)               :: vol, boxl, kr, ckr, skr, dk, kmag, kmax
    real(kind=dp)               :: beta, k(3), kx, ky, kz

    real(kind=dp), allocatable  :: r(:,:), sk(:), s(:)
    integer, allocatable        :: sc(:)

!   npart: number of particles
    integer, parameter          :: npart = 1000
!   ndump: number of snapshots
    integer, parameter          :: ndump = 25
!   ap: aspect ratio
    real(kind=dp), parameter    :: ap=0.2, alpha=1.0_dp+1.5_dp*ap-0.5_dp*ap**3
!   rho: density/ packing fraction
    real(kind=dp), parameter    :: rho = 0.5_dp
!   packingt: true if using packing fraction
    logical, parameter          :: packingt = .false.
!   dumbellt: true if using dumbbells
    logical, parameter          :: dumbellt = .false.
!   tetrahedrat: true if using tetrahedra centres
    logical, parameter          :: tetrahedrat = .false.
!   nk: maximum value of one component of the wavevector (i.e., kmax = [nk,nk,nk])
    integer, parameter          :: nk = 60, ntot = 2*nk+1, ntot2 = ntot*ntot
!   kres: resolution of wavevectors for plotting purposes 
    integer, parameter          :: kres = 180
!   path to file containing coordinates of the particles
    open (unit = 33, file = '../pos.dat', status = 'old')

    if(.not. tetrahedrat) then 
        np = npart
        allocate(r(3,npart))
        r  = 0.0_dp
    endif

    allocate( sk(0:ntot**3), s(0:kres), sc(0:kres) )
    
    vol = real(npart,dp) / rho
    if(dumbellt) vol = vol*alpha
    if(packingt) vol = vol*(pi/6.0_dp)
    boxl = vol**(1.0_dp/3.0_dp)
    
    sk   = 0.0_dp
    s    = 0.0_dp
    sc   = 0
    beta = 2.0_dp*pi / boxl

    kmax = beta*sqrt(3.0_dp)*nk
    dk   = kmax / real(kres,dp)

    print *, boxl, kmax, dk

    do j1 = 1, ndump

        print *, "On Snapshot:", j1

        if(tetrahedrat) then 
            read(33,*) np
            allocate(r(3,np))
            r  = 0.0_dp
        endif
        
        do j2 = 1, np
            read(33,*) r(:,j2)
        enddo

        do nx = -nk, nk
            kx  = beta*real(nx,dp)
            nxi = (nx+nk)*ntot2 + nk
            do ny = -nk, nk
                ky  = beta*real(ny,dp)
                nyi = nxi + (ny+nk)*ntot
                do nz = -nk, nk
                    kz = beta*real(nz,dp)
                    k  = [kx,ky,kz]
                    ckr = 0.0_dp
                    skr = 0.0_dp
                    do j2 = 1, np
                        kr = dot_product(k,r(:,j2))
                        ckr = ckr + cos(kr)
                        skr = skr + sin(kr)
                    enddo
                    ki     = nyi + nz
                    sk(ki) = sk(ki) + ckr**2 + skr**2
                enddo
            enddo
        enddo

        if(tetrahedrat) then 
            sk = sk/real(np,dp)
            deallocate(r)
        endif
    
    enddo

    close(33)

    if(.not. tetrahedrat) sk = sk / ( real(npart,dp)*real(ndump,dp) )

    open(unit = 34, file ='sk.dat', status = 'replace', access = 'append')

    do nx = -nk, nk
        kx  = beta*real(nx,dp)
        nxi = (nx+nk)*ntot2 + nk
        do ny = -nk, nk
            ky  = beta*real(ny,dp)
            nyi = nxi + (ny+nk)*ntot
            do nz = -nk, nk
                kz     = beta*real(nz,dp)
                k      = [kx,ky,kz] 
                kmag   = norm2(k)
                rk     = int(anint(kmag/dk))
                ki     = nyi + nz
                s(rk)  = s(rk) + sk(ki)
                sc(rk) = sc(rk) + 1
            enddo
        enddo
    enddo

    do j1 = 1, kres
        kmag = dk*real(j1,dp)
        if(sc(j1)>0) then
            write(34,*) kmag, s(j1)/real(sc(j1),dp)
        else
            write(34,*) kmag, 0.0
        endif
    enddo

end program
