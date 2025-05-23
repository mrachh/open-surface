      implicit real *8 (a-h,o-z)
      integer norders_jac(4)

      norders_jac(1) = 20
      norders_jac(2) = 30
      norders_jac(3) = 60
      norders_jac(4) = 120
      do i=1,5
        nuse = 5*i
        do j=4,4
          call compute_error(norders_jac(j), nuse)
        enddo
      enddo
      
      
      stop
      end

      subroutine compute_error(norder_jac, nuse)
      implicit real *8 (a-h,o-z)
      real *8 rmid
      integer npars(3)
      real *8, allocatable :: ptinfo(:,:), srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: ptcoefs(:,:)
      integer, allocatable :: ixyzs(:), iptype(:), norders(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:), qwts(:)
      real *8, allocatable :: cms(:,:), rads(:)
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
      complex *16, allocatable :: charges(:)
      complex *16, allocatable :: pvals(:), pex(:)
      complex *16, allocatable :: taexp(:,:,:)

      complex *16 zk
      character *100 fname 
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      external h3d_slp_disk, h3d_slp
      
      call prini(6,13)
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
        srcvals(1,i) = ptinfo(1,i)
        srcvals(2,i) = ptinfo(2,i)
        srcvals(3,i) = 0

        srcvals(4,i) = ptinfo(3,i)
        srcvals(5,i) = ptinfo(4,i)
        srcvals(6,i) = 0

        srcvals(7,i) = ptinfo(5,i)
        srcvals(8,i) = ptinfo(6,i)
        srcvals(9,i) = 0

        srcvals(10,i) = 0
        srcvals(11,i) = 0
        srcvals(12,i) = 1
        
        srccoefs(1,i) = ptcoefs(1,i)
        srccoefs(2,i) = ptcoefs(2,i)
        srccoefs(3,i) = 0

        srccoefs(4,i) = ptcoefs(3,i)
        srccoefs(5,i) = ptcoefs(4,i)
        srccoefs(6,i) = 0

        srccoefs(7,i) = ptcoefs(5,i)
        srccoefs(8,i) = ptcoefs(6,i)
        srccoefs(9,i) = 0
      enddo

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
      
      rsum = 0
      do i=1,nptso
        rsum = rsum + wtsover(i)
        dd = srcover(1,i)**2 + srcover(2,i)**2
        charges(i) = 1.0d0/sqrt(1.0d0 - dd)*wtsover(i)
      enddo
      call prin2('area of circle=*',rsum,1)
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
          dd = srcover(1,jpt)**2 + srcover(2,jpt)**2
          srcover(12,jpt) = iort/abs(iort)
          charges(jpt) = 1.0d0/sqrt(1.0d0 - dd)* & 
             wtsover(jpt)*sqrt(1-uvs_jac(2,j))
          if(1.0d0-dd.le.0) then
            print *, "here0"
            stop

          endif
        enddo
      enddo


      rsum = 0
      do i=1,nptso
        rsum = rsum + real(charges(i))
      enddo


      done = 1.0d0
      pi = atan(done)*4.0d0

      call prin2('rsum=*',rsum,1)
      call prin2('error in integral=*',abs(rsum-2*pi),1)

      rsum2 = 0
      ipatch = npars(1)*npars(1) + npars(2)*1
      call prinf('ipatch=*',ipatch,1)
      do j=ixyzso(ipatch),ixyzso(ipatch+1)-1
        rsum2 = rsum2 + real(charges(j))
      enddo
      call prin2('rsum2=*',rsum2,1)

      print *, ""
      print *, ""
      print *, "====================="
      print *, ""
      print *, "Qbx testing begins now"
      print *, ""
      print *, "====================="
 
      allocate(cms(3,npatches), rads(npatches))
      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, cms, rads)

      allocate(centers(3,npts))
      do i=1,npatches
        do j = ixyzs(i), ixyzs(i+1)-1
          d1 = rads(i)
          dd = (1.0d0 - sqrt(srcvals(1,j)**2 + srcvals(2,j)**2))/2
          if(dd.le.d1) d1 = dd
          centers(1,j) = srcvals(1,j)
          centers(2,j) = srcvals(2,j)
          centers(3,j) = d1
        enddo
      enddo

      allocate(pvals(npts), pex(npts))
      nlege = 100
      lw7 = (nlege+1)**2*4
      allocate(wlege(lw7))
      call ylgndrfwini(nlege, wlege, lw7, lused7)


      zk = 2.1d0
      
      allocate(taexp(0:nuse,-nuse:nuse,npts))
      rscales = 1.0d0
      charges = charges/4/pi

      fname = './data/results.dat'
      open(unit=35,file=trim(fname), access='append')

      do i=1,npts
        x = srcvals(1,i)
        y = srcvals(2,i)
        rr = sqrt(x**2 + y**2)
        ruse = 2*rr - 1

        thet = atan2(y,x)
        pvals(i) = 0
        pex(i) = 0
        
        if(abs(thet).le.0.01d0.and.rr.gt.0.95d0) then
          taexp = 0
          call h3dformtac(1, zk, rscales, srcover(1:3,1:nptso), &
            charges, nptso, centers(1,i), nuse, taexp(0,-nuse,i), &
            wlege, nlege)
          
          call h3dtaevalp(1, zk, rscales, centers(1,i), &
            taexp(0,-nuse,i), nuse, srcvals(1,i), 1, pvals(i), &
            wlege, nlege)

          call legepols(ruse, nfcoef-1, pols)
          do j=1,nfcoef
            pex(i) = pex(i) + fcoefs(j)*pols(j)
          enddo

          write(35,*) nuse, norder_jac, 1.0d0-rr, centers(3,i), abs(pvals(i)-pex(i))
          print *, nuse, norder_jac, abs(pvals(i)-pex(i))
        endif

      enddo
      close(35)
      
      
       



      return
      end
