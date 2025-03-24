      implicit real *8 (a-h,o-z)
      real *8 rmid
      integer npars(3)
      real *8, allocatable :: ptinfo(:,:), srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: ptcoefs(:,:)
      integer, allocatable :: ixyzs(:), iptype(:), norders(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:), qwts(:)
      complex *16, allocatable :: fints(:)
      complex *16, allocatable :: amat(:,:)
      real *8 xj20(20), wj20(20)
      real *8 xl20(20), wl20(20)
      real *8 xj30(30), wj30(30)
      real *8 vtmp(3)
      real *8 umattmp, vmattmp
      complex *16 fcoefs(20)
      integer, allocatable :: nfars(:), ixyzso(:)
      real *8, allocatable :: srcover(:,:), wtsover(:)
      real *8, allocatable :: uvs_jac(:,:), wts_jac(:)
      real *8, allocatable :: pmat(:,:)
      complex *16, allocatable :: charges(:)

      complex *16 zk
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      external h3d_slp_disk, h3d_slp
      
      call prini(6,13)
      rmid = 0.3d0
      npars(1) = 3
      npars(2) = 3
      npars(3) = 5
      
      iort = 1
      

      norder = 9
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


      open(unit=32, file='jacpts20.bin', access='stream', &
        form='unformatted')
      read(32) xj20
      read(32) wj20
      close(32)

      call prin2('xj20=*',xj20,20)
      call prin2('wj20=*',wj20,20)

      
      open(unit=32, file='jacpts30.bin', access='stream', &
        form='unformatted')
      read(32) xj30
      read(32) wj30
      close(32)

      call prin2('xj30=*',xj30,30)
      call prin2('wj30=*',wj30,30)

      open(unit=32, file='fcoefs.bin', access='stream', &
        form='unformatted')
      read(32) fcoefs 
      close(32)


!
!  generate oversampled surface
!
!
      allocate(nfars(npatches), ixyzso(npatches+1))
      norder_jac = 20 
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
      call legeexps(itype, norder_jac, xl20, utmp, vtmp, wl20)
      npol_jac = norder_jac*norder_jac
      allocate(uvs_jac(2,npol_jac))
      allocate(wts_jac(npol_jac))
      call prin2('xl20=*',xl20,20)
      call prin2('wl20=*',wl20,20)

      rsum = 0.0d0
      do i=1,norder_jac
        rsum = rsum + wj20(i)
      enddo
      call prin2('sum of jacobi weights=*',rsum,1)
      erra = rsum - sqrt(8.0d0)
      call prin2('error in jacobi weights=*',erra,1)

      
      nptmp = (norder+1)*(norder+1)
      allocate(pmat(nptmp, npol_jac))

      do i=1,norder_jac
        do j=1,norder_jac
          ipt = (i-1)*norder_jac + j
          uvs_jac(1,ipt) = xl20(j)
          uvs_jac(2,ipt) = xj20(i)

          wts_jac(ipt) = wl20(j)*wj20(i)
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
        enddo
      enddo

      call prin2('uvs_jac=*',uvs_jac(2,1:npol_jac),npol_jac)

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
    


      



      stop
      end
