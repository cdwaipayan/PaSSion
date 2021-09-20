program histogram_reweighting
    implicit none
    integer, parameter :: sp = selected_real_kind(6, 37), dp = selected_real_kind(15, 307), qp = selected_real_kind(33, 4931)
    
    integer                    :: resolution, ntemps, npress, j1, j2, j3, j4
    
    real(kind=qp)              :: line(3), nc, v_min, v_max, e_min, e_max, b_min, b_max, p_min, p_max
    real(kind=qp)              :: ec, vc
    real(kind=qp), allocatable :: volume(:), energy(:), counts(:), temps(:), press(:)
    
    real(kind=qp)              :: thresh, gg, oi, c1, c3, c4, e, v, b, p, e1, e2, gi
    real(kind=qp), allocatable :: GDiff(:)
    real(kind=qp), allocatable :: S(:), G(:), GOld(:)

    integer                    :: npnew, vlm_count, vlm_res
    real(kind=qp)              :: ti, bi, cnst, counter
    real(kind=qp), allocatable :: pi(:), tmp(:), Pr(:), vols(:)

    character(len=100)         :: arg, fmt1, string

    logical                    :: readin

    ! nc         = 3200000.0_qp
    ! resolution = 22673
    nc         = 2400000.0_qp
    resolution = 22132
    ntemps     = 7
    npress     = 9
    npnew      = 16
    vlm_res    = 6
    thresh     = 1.e-2
    readin     = .true.

    allocate(volume(resolution), energy(resolution), counts(resolution), temps(ntemps), press(npress))
    
    open (unit = 33, file = '../merged_histogram.dat', status = 'old')
    
    do j1 = 1, resolution
        read(33,*) line
        volume(j1) = line(1)
        energy(j1) = line(2)
        counts(j1) = log(line(3)/nc)
    enddo

    close(unit=33)

    temps  = [0.105_qp,0.106_qp,0.107_qp,0.108_qp,0.109_qp,0.11_qp,0.112_qp]
    press  = [0.005_qp,0.006_qp,0.0065_qp,0.007_qp,0.0075_qp,0.008_qp,0.0085_qp,0.009_qp,0.01_qp]

    v_min = minval(volume);       v_max = maxval(volume)
    e_min = minval(energy);       e_max = maxval(energy)
    p_min = minval(press);        p_max = maxval(press)
    b_min = 1.0_qp/maxval(temps); b_max = 1.0_qp/minval(temps)

    ec = (e_max+e_min)/2.0_qp
    vc = (v_max+v_min)/2.0_qp

    allocate(S(resolution), G(ntemps*npress), GOld(ntemps*npress), GDiff(ntemps*npress))
    
    S     = 0.0_qp
    G     = 0.0_qp
    GOld  = 0.0_qp
    GDiff = 0.0_dp
    gg    = 1.0_dp

    if(readin) then
        open (unit = 35, file = '../../gibbs.dat', status = 'old')
        do j1 = 1, ntemps*npress
            read(35,*) line
            G(j1) = line(3)
        enddo
        close(unit=35)
    endif

    do while(gg > thresh)
        GOld = G
        do j1 = 1, resolution
            e  = energy(j1)
            v  = volume(j1)
            oi = 0.0_qp
            c1 = -b_max*(e+p_max*v)
            do j2 = 1, ntemps
                b = 1.0_qp / temps(j2)
                do j3 = 1, npress
                    j4 = (j2-1)*npress + j3
                    p  = press(j3)
                    e1 = -b*(e+p*v) - c1
                    e2 = b*G(j4)
                    oi = oi + exp(e1+e2)
                enddo
            enddo
            S(j1) = counts(j1) - (c1+log(oi))
        enddo

        c4 = (maxval(S)+minval(S))/2.0_qp
        do j1 = 1, ntemps
            b = 1.0_qp / temps(j1)
            do j2 = 1, npress
                j4 = (j1-1)*npress + j2
                gi = 0.0_qp
                c3 = -b*(ec + p*vc)
                p  = press(j2) 
                do j3 = 1, resolution
                    e = energy(j3)
                    v = volume(j3)
                    gi = gi + exp( (S(j3)-c4) + (-b*(e+p*v)-c3) )
                enddo
                G(j4) = -temps(j1)*( c3+c4+log(gi) )
            enddo
        enddo

        GDiff = G-GOld
        gg    = norm2(GDiff)
        print *, real(gg,sp)
    enddo

    open (unit = 36, file = './gibbs.dat', status='replace', access='append')
    do j1 = 1, ntemps
        do j2 = 1, npress
            j4 = (j1-1)*npress + j2
            write(36,*) real(temps(j1),dp), real(press(j2),dp), G(j4)
        enddo
    enddo
    close(unit=36)

    allocate( pi(npnew), tmp(resolution) )

!   https://gcc.gnu.org/onlinedocs/gfortran/GET_005fCOMMAND_005fARGUMENT.html
    call get_command_argument(1, arg) ! result is stored in arg, see 
    if(len_trim(arg) == 0) then ! no argument given
        stop "No input temperature provided"
    else
        string = trim(arg)
        read(string,'(f10.4)')ti
        bi = 1.0_qp / real(ti,qp)
    endif

    pi = [0.005_qp, 0.006_qp, 0.007_qp, 0.0075_qp, 0.0078_qp, 0.0079_qp, 0.008_qp, 0.008025_qp, & 
          0.00805_qp, 0.0081_qp, 0.0082_qp, 0.0083_qp, 0.0084_qp, 0.0086_qp, 0.009_qp, 0.01_qp]

    vlm_count = int(((v_max - v_min - 1) / vlm_res)) + 1
    allocate(vols(vlm_count), Pr(vlm_count))

    do j1 = 1, vlm_count
        vols(j1) = v_min + real(j1,qp)*vlm_res
    enddo

    do j1 = 1, npnew
        tmp = 0.0_qp
        do j2 = 1, resolution
            tmp(j2) = S(j2) - bi*(energy(j2)+pi(j1)*volume(j2))
        enddo

    !   Calculate probabilities for energy, volume pairs.
    !   Density of states shifted by an arbitrary constant not accounted for when 
    !   performing self-consistent equations. Shift is required in order to allow
    !   for the computation of exponential.
        cnst = minval(tmp)
        tmp  = exp( tmp-cnst )
        tmp  = tmp / sum(tmp)

        write(fmt1,'(I3)') j1
        string = 'dos_'//trim(adjustl(fmt1))//'.dat'
        open (unit=34, file=string, status='replace', access='append')
        do j2 = 1, resolution
            write(34,*) real(volume(j2),dp), real(energy(j2),dp), tmp(j2)
        enddo
        close(unit = 34)

    !   Calculate probability distributions for volume
        Pr       = 0.0_qp
        do j2 = 1, vlm_count
            counter = 0.0_qp
            do j3 = 1, resolution
                if(vols(j2)==volume(j3)) then
                    counter = counter + tmp(j3)
                endif
            enddo
            Pr(j2) = counter
        enddo

        write(fmt1,'(I3)') j1
        string = 'rho_'//trim(adjustl(fmt1))//'.dat'
        open (unit=36, file=string, status='replace', access='append')
        do j2 = 1, vlm_count
            write(36,*) 1000.0_dp/real(vols(j2),dp), Pr(j2)
        enddo
        close(unit = 36)

    enddo

end program