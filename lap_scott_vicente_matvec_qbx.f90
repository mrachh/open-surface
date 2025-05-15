      implicit real *8 (a-h,o-z)
      real *8 rmid
      integer npars(3)
      real *8, allocatable :: ptinfo(:,:), srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: ptcoefs(:,:)
      integer, allocatable :: ixyzs(:), iptype(:), norders(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:), qwts(:)
      real *8, allocatable :: cms(:,:), rads(:)
      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: centers(:,:)
      real *8, allocatable :: wlege(:)
      complex *16, allocatable :: fints(:)
      complex *16, allocatable :: amat(:,:)
      real *8 pols(100)
      real *8 xj(200), wj(200)
      real *8 xl(200), wl(200)
      real *8 vtmp(3)
      real *8 umattmp, vmattmp
      complex *16 fcoefs(20)
      integer, allocatable :: nfars(:), ixyzso(:)
      real *8, allocatable :: srcover(:,:), wtsover(:)
      real *8, allocatable :: uvs_jac(:,:), wts_jac(:)
      real *8, allocatable :: pmat(:,:)
      real *8, allocatable :: charges(:)
      real *8, allocatable :: pvals(:)
      complex *16, allocatable :: taexp(:,:,:)
      complex *16 z, ztmp, ima

      complex *16 zk
      character *100 fname, ftarg, dirname 
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      external h3d_slp_disk, h3d_slp

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4.0d0
      norder_jac = 120
      nuse = 60

      dirname = 'scott_vicente/'
      ep = 0.0d0

      if(ep.eq.0) then
        ftarg = trim(dirname)//'eps_0.csv'
      elseif(abs(ep-1.0d-1).le.1.0d-14) then
        ftarg = trim(dirname)//'eps_1e-1.csv'
      elseif(abs(ep-1.0d-2).le.1.0d-14) then
        ftarg = trim(dirname)//'eps_1e-2.csv'
      elseif(abs(ep-1.0d-3).le.1.0d-14) then
        ftarg = trim(dirname)//'eps_1e-3.csv'
      endif

      ntarg = 1280
      allocate(targs(3,ntarg), centers(3,ntarg))

      open(unit=37, file=trim(ftarg))
      do i=1,ntarg
        read(37,*) targs(1,i), targs(2,i)
        targs(3,i) = 0
      enddo

      call prin2('targs=*',targs,24)
      
  
      npow = 5
      
      rmid = 0.3d0
      npars(1) = 6
      npars(2) = 6
      npars(3) = 10
      
      iort = 1
      
      norder = 4
      iptype0 = 11
      call mesh_circle_pts_mem(rmid, npars, iort, norder, iptype0, &
        npatches, npts)
      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)
      
      allocate(ptinfo(6,npts), srcvals(12,npts), srccoefs(9,npts))
      allocate(ixyzs(npatches+1), iptype(npatches), norders(npatches))
      call mesh_circle_pts(rmid, npars, iort, norder, iptype0, &
        npatches, npts, ixyzs, ptinfo)
      do i=1,npatches
        iptype(i) = iptype0
        norders(i) = norder
      enddo
      call prinf('iptype=*',iptype,20)
      call prinf('norders=*',norders,20)
      allocate(ptcoefs(6,npts))
      call surf_vals_to_coefs(6, npatches, norders, ixyzs, iptype, npts, &
        ptinfo, ptcoefs)

      do i=1,npts
        x = ptinfo(1,i)
        y = ptinfo(2,i)
        r = sqrt(x**2 + y**2)
        thet = atan2(y,x)
        cnt = cos((n+1)*thet)
        snt = sin((n+1)*thet)
        rn1 = r**(n+1)
        srcvals(1,i) = ptinfo(1,i) + ep*rn1*cnt
        srcvals(2,i) = ptinfo(2,i) + ep*rn1*snt
        srcvals(3,i) = 0

        drdx = x/r
        drdy = y/r

        dthetdx = -y/r**2
        dthetdy = x/r**2

        df1dx = 1 + ep*(n+1)*rn1*(cnt*x + y*snt)/r**2
        df1dy = ep*(n+1)*rn1*(cnt*y - x*snt)/r**2
        
        df2dx = ep*(n+1)*rn1*(snt*x - y*cnt)/r**2
        df2dy = 1 + ep*(n+1)*rn1*(snt*y + x*cnt)/r**2
                   
        dxdu = ptinfo(3,i)
        dydu = ptinfo(4,i)
        
        dxdv = ptinfo(5,i)
        dydv = ptinfo(6,i)

        srcvals(4,i) = df1dx*dxdu + df1dy*dydu 
        srcvals(5,i) = df2dx*dxdu + df2dy*dydu 
        srcvals(6,i) = 0

        srcvals(7,i) = df1dx*dxdv + df1dy*dydv 
        srcvals(8,i) = df2dx*dxdv + df2dy*dydv 
        srcvals(9,i) = 0

        srcvals(10,i) = 0
        srcvals(11,i) = 0
        srcvals(12,i) = 1
        
      enddo
      call surf_vals_to_coefs(9, npatches, norders, ixyzs, iptype, npts, &
        srcvals(1:9,1:npts), srccoefs) 

      call prin2('srcvals=*',srcvals,24)
      call prin2('srccoefs=*',srccoefs,24)

!
!
!
      nptuse = 2
      allocate(ipatch_id(npts), uvs_src(2,npts))
      allocate(fints(nptuse))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_src)



      if(norder_jac.eq.20) then
        open(unit=32, file='jacpts20.bin', access='stream', &
          form='unformatted')
        read(32) xj(1:20)
        read(32) wj(1:20)
        close(32)
      elseif(norder_jac.eq.30) then
        open(unit=32, file='jacpts30.bin', access='stream', &
          form='unformatted')
        read(32) xj(1:30)
        read(32) wj(1:30)
        close(32)
      elseif(norder_jac.eq.60) then
        open(unit=32, file='jacpts60.bin', access='stream', &
          form='unformatted')
        read(32) xj(1:60)
        read(32) wj(1:60)
        close(32)
      elseif(norder_jac.eq.120) then
        open(unit=32, file='jacpts120.bin', access='stream', &
          form='unformatted')
        read(32) xj(1:120)
        read(32) wj(1:120)
        close(32)
      endif


      nfcoef = 20
      open(unit=32, file='fcoefs.bin', access='stream', &
        form='unformatted')
      read(32) fcoefs 
      close(32)


!
!  generate oversampled surface
!
!
      allocate(nfars(npatches), ixyzso(npatches+1))
      do i=1,npatches
        nfars(i) = norder_jac - 1 
        ixyzso(i) = norder_jac*norder_jac*(i-1) + 1
      enddo
      nptso = norder_jac*norder_jac*npatches
      ixyzso(npatches+1) = nptso + 1
      allocate(srcover(12,nptso), wtsover(nptso))

      call prinf('nfars=*', nfars, 20)

      call oversample_geom(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, nfars, ixyzso, nptso, srcover)
      
      call get_qwts(npatches, nfars, ixyzso, iptype, nptso, srcover, &
        wtsover)

      allocate(charges(nptso))
      
      do i=1,nptso
        z = srcover(1,i) + ima*srcover(2,i)
        ztmp = z*exp(-ep*(z**n))
        dd = abs(ztmp)**2 
        charges(i) = 1.0d0/sqrt(1.0d0 - dd)*wtsover(i)
      enddo
!
!
!  now replace patches on the periphery with the tensor 
!  product jacobi-legendre patches
!
!
      itype = 1
      call legeexps(itype, norder_jac, xl, utmp, vtmp, wl)
      npol_jac = norder_jac*norder_jac
      allocate(uvs_jac(2,npol_jac))
      allocate(wts_jac(npol_jac))
      call prin2('xl=*',xl,20)
      call prin2('wl=*',wl,20)

      rsum = 0.0d0
      do i=1,norder_jac
        rsum = rsum + wj(i)
      enddo
      call prin2('sum of jacobi weights=*',rsum,1)
      erra = rsum - sqrt(8.0d0)
      call prin2('error in jacobi weights=*',erra,1)

      
      nptmp = (norder+1)*(norder+1)
      allocate(pmat(nptmp, npol_jac))

      do i=1,norder_jac
        do j=1,norder_jac
          ipt = (i-1)*norder_jac + j
          uvs_jac(1,ipt) = xl(j)
          uvs_jac(2,ipt) = xj(i)

          wts_jac(ipt) = wl(j)*wj(i)
          call get_basis_pols(uvs_jac(1,ipt), norder, nptmp, iptype0, &
             pmat(1,ipt))
        enddo
      enddo

      do i=1,4*npars(3)
        ipatch = npars(1)*npars(1) + npars(2)*i

        do j = 1,npol_jac
          jpt = ixyzso(ipatch) + j - 1
          srcover(1:12,jpt) = 0
          do l=1,nptmp
            ipt = ixyzs(ipatch) + l - 1
            srcover(1:9,jpt) = srcover(1:9,jpt) + &
               srccoefs(1:9,ipt)*pmat(l,j)
          enddo
          call cross_prod3d(srcover(4,jpt), srcover(7,jpt), vtmp)
          wtsover(jpt) = sqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)* &
            wts_jac(j)
          z = srcover(1,jpt) + ima*srcover(2,jpt)
          ztmp = z*exp(-ep*(z**n))
          dd = abs(ztmp)**2 
          srcover(12,jpt) = iort/abs(iort)
          charges(jpt) = 1.0d0/sqrt(1.0d0 - dd)* & 
             wtsover(jpt)*sqrt(1-uvs_jac(2,j))
          if(1.0d0-dd.le.0) then
            print *, "here0"
            stop

          endif
        enddo
      enddo


      print *, ""
      print *, ""
      print *, "====================="
      print *, ""
      print *, "Qbx testing begins now"
      print *, ""
      print *, "====================="

      
      print *, "ntarg=",ntarg
      do i=1,ntarg
        d1 = 0.25d0
        z = targs(1,i) + ima*targs(2,i)
        ztmp = z*exp(-ep*(z**n))
        dd = abs(ztmp)**2 
        dd = (1.0d0 - sqrt(dd))/2
        if(dd.le.d1) d1 = dd
        centers(1,i) = targs(1,i)
        centers(2,i) = targs(2,i) 
        centers(3,i) = d1
      enddo
      call prin2('centers=*',centers,24)
      call prin2('targs=*',targs,24)


      allocate(pvals(ntarg))
      nlege = 100
      lw7 = (nlege+1)**2*4
      allocate(wlege(lw7))
      call ylgndrfwini(nlege, wlege, lw7, lused7)

      
      allocate(taexp(0:nuse,-nuse:nuse,npts))
      rscales = 1.0d0
      charges = charges/4/pi

      fname = trim(dirname)//'results_lap_scott_vicente.dat'
      open(unit=35,file=trim(fname), access='append')
!!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,30
        print *, "itarg=",i
        pvals(i) = 0
        
        if(abs(centers(3,i)).ge.0.001d0) then
          taexp(0:nuse,-nuse:nuse,i) = 0
          call l3dformtac(1, rscales, srcover(1:3,1:nptso), &
            charges, nptso, centers(1,i), nuse, taexp(0,-nuse,i), &
            wlege, nlege)

          
          call l3dtaevalp(1, rscales, centers(1,i), &
            taexp(0,-nuse,i), nuse, targs(1,i), 1, pvals(i), &
            wlege, nlege)
        endif

      enddo
!!!$OMP END PARALLEL DO
      do i=1,30
        write(35,*) nuse, norder_jac, targs(1,i), targs(2,i), &
             centers(3,i), pvals(i)
      enddo
      close(35)
      
      
       



      return
      end
