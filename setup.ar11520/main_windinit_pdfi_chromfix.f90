Program main
       Use problemsize
       Use grid
       Use blk
       Use blkdata
       Use input_namelists
       implicit none
!
       integer i, j, k, jn_in, kn_in
       real(8), allocatable :: b1(:,:)
       real(8), allocatable :: bp1(:,:,:)
       real(8), allocatable :: ph(:,:,:)
       real(8), allocatable :: tfield(:,:,:)
       real(8) delx
       real(8) anq,amq,cnq,cmq,err,referr,chk,rhserr,rhsnorm,qnorm
       real(8) g_const,rgas,kboltz,mproton,mpovme,gamma,muconst, &
               unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time, &
               msol,rsol,gacc_s,gacc,d0,rho_c,temp_c,rstar, &
               as2_c, hp_c, as2, ovhp0
       real(8), allocatable :: fint(:)
       real(8) tovt0
!
       real(8) r1, r2, pr1, pr2
!
!--------------------------------------------------------
!      read in input parameters and initialization
!--------------------------------------------------------
!
       write(6,*) 'read input parameters'
       call flush(6)
       call main_input
!
       call initialize_problemsize
!
       call initialize_grid
!
!---------------------------------------------------------
!     initialize blk
!---------------------------------------------------------
!
       allocate(ar(inmax-5))
       allocate(br(inmax-5))
       allocate(cr(inmax-5))
       allocate(ath(jnmax-5))
       allocate(bth(jnmax-5))
       allocate(cth(jnmax-5))
       allocate(bthk(jnmax-5))
!
       np_blk=1
       mp_blk=1
       m_blk=inmax-5
       n_blk=jnmax-5
       idimy_blk=inmax-5
       kk=dlog(dble(n_blk))/dlog(dble(2))+1
       ll=2**(kk+1)
       nww=(kk-2)*ll+kk+5+max(2*n_blk,6*m_blk)
       nrlpts=kunif*2
       write(6,*) 'nww= ', nww
       write(6,*) 'in,jn,kn,m_blk,n_blk,nrlpts', &
         in,jn,kn,m_blk,n_blk,nrlpts
       call flush(6)
!
       allocate(an(n_blk))
       allocate(bn(n_blk))
       allocate(cn(n_blk))
       allocate(am(m_blk))
       allocate(bm(m_blk))
       allocate(cm(m_blk))
       allocate(yy(m_blk,n_blk))
       allocate(ww(nww))
       allocate(rldat(nrlpts))
       allocate(wsave(2*nrlpts+15))
       allocate(qq(in-5,jn-5,nrlpts))
!
       do i=is,in-3
        cr(i-is+1)=dx1ai(i) &
          *g2a(i+1)*g31a(i+1)*dx1bi(i+1)
        ar(i-is+1)=dx1ai(i) &
          *g2a(i)*g31a(i)*dx1bi(i)
        br(i-is+1)=-ar(i-is+1)-cr(i-is+1)
       enddo
!      
       delx=(x3a(kep1)-x3a(ks))/dble(nrlpts/2)
       write(6,*) 'delx=',delx
       do j=js,jn-3
        cth(j-js+1)=g32bi(j)*dx2ai(j) &
          *g32a(j+1)*dx2bi(j+1)
        ath(j-js+1)=g32bi(j)*dx2ai(j) &
          *g32a(j)*dx2bi(j)
        bth(j-js+1)=-cth(j-js+1)-ath(j-js+1)
        bthk(j-js+1)=1.D0/delx/delx &
          *g32bi(j)*g32bi(j)
       enddo
!
       call rffti(nrlpts,wsave)
!
!---------------------------------------------------------
!      start computing the potential field
!---------------------------------------------------------
!
       allocate(b1(jn,kn))
       allocate(bp1(in,jn,kn))
       allocate(ph(in,jn,kn))
       allocate(tfield(in,jn,kn))
!
         open(unit=13, file=trim(lbfile), &
           form='unformatted',access='stream')
         read(13) jn_in,kn_in
         if(jn_in .ne. jn &
           .or. kn_in .ne. kn) then
           write(6,*) 'jn_in kn_in= ', &
             jn_in, kn_in, 'but jn kn = ', jn,kn
           stop
         endif
         read(13) ((b1(j,k),j=js,je),k=ks,ke)
         close(13)
         do j=js,je
           b1(j,ksm1)=b1(j,ks)
           b1(j,ksm2)=b1(j,ks+1)
           b1(j,kep1)=b1(j,ke)
           b1(j,kep2)=b1(j,ke-1)
           b1(j,kep3)=b1(j,ke-2)
         enddo
!
         write(6,*) 'x3b(ksm2),x3b(kep2)', x3b(ksm2),x3b(kep2)
         call flush(6)
         do j=js,je
           i=is
           do k=1,nrlpts/2
             rldat(k)=-b1(j,k-1+ks)*dx1b(is)*ar(1)
           enddo
           do k=nrlpts/2+1,nrlpts
             rldat(k)=rldat(nrlpts/2-(k-nrlpts/2-1))
           enddo
           call rfftf(nrlpts,rldat,wsave)
           do k=1,nrlpts
             qq(i-is+1,j-js+1,k)=rldat(k)/dble(nrlpts)
           enddo
           do i=is+1,ie
             do k=1,nrlpts
               qq(i-is+1,j-js+1,k)=0.D0
             enddo
           enddo
         enddo
         write(6,*) 'computed rhs'
         call flush(6)
!
         do k=1,nrlpts
!
           kint=k/2
           phk=2.D0*pi*dble(kint)/dble(nrlpts)

           do i=1,m_blk
             am(i)=ar(i)
             bm(i)=br(i)
             cm(i)=cr(i)
           enddo
           am(1)=0.D0
           bm(1)=bm(1)+ar(1)
           cm(m_blk)=0.D0
           bm(m_blk)=bm(m_blk)-cr(m_blk)
!
           do j=1,n_blk
             an(j)=ath(j)
             bn(j)=bth(j) &
               -2.D0*(1.D0-dcos(phk))*bthk(j)
             cn(j)=cth(j)
           enddo
           an(1)=0.D0
           bn(1)=bn(1)+ath(1)
           cn(n_blk)=0.D0
           bn(n_blk)=bn(n_blk)+cth(n_blk)
           iflg_blk=0
           call blktri(iflg_blk,np_blk,n_blk,an,bn,cn,mp_blk, &
             m_blk,am,bm,cm,idimy_blk,yy,ierr_blk,ww)
!
           do j=js,jn-3
             do i=is,in-3
               yy(i-is+1,j-js+1)=qq(i-is+1,j-js+1,k)
             enddo
           enddo
           iflg_blk=1
           call blktri(iflg_blk,np_blk,n_blk,an,bn,cn,mp_blk, &
             m_blk,am,bm,cm,idimy_blk,yy,ierr_blk,ww)
!
!          checking
!
           err=0.D0
           referr=0.D0
           do j=1,n_blk
             do i=1,m_blk
               amq=0.D0
               cmq=0.D0
               anq=0.D0
               cnq=0.D0
               if (i .gt. 1) amq=yy(i-1,j)
               if (i .lt. m_blk) cmq=yy(i+1,j)
               if (j .gt. 1) anq=yy(i,j-1)
               if (j .lt. n_blk) cnq=yy(i,j+1)
               chk=an(j)*anq+am(i)*amq+(bn(j)+bm(i))*yy(i,j) &
                      +cn(j)*cnq+cm(i)*cmq
               err=dmax1(err,abs(chk-qq(i,j,k)))
               referr=dmax1(referr,abs(qq(i,j,k)))
             enddo
           enddo
           write(6,*) 'k, kint, ierr_blk=',k,kint,ierr_blk
           write(6,*) 'err,referr',err,referr
           call flush(6)
!
           do j=js,jn-3
             do i=is,in-3
               qq(i-is+1,j-js+1,k)=yy(i-is+1,j-js+1)
             enddo
           enddo
!
         enddo
!
         do j=js,je
           do i=is,ie
             do k=1,nrlpts
               rldat(k)=qq(i-is+1,j-js+1,k)
             enddo
             call rfftb(nrlpts,rldat,wsave)
             do k=ks,ke
               ph(i,j,k)=rldat(k-ks+1)
             enddo
           enddo
         enddo
!
!        boundaries of ph
!
         do k=ks,ke
           do j=js,je
             ph(ism1,j,k)=ph(is,j,k)+dx1b(is)*b1(j,k)
             ph(iep1,j,k)=-ph(ie,j,k)
           enddo
         enddo
         do k=ks,ke
           do i=ism1,iep1
             ph(i,jsm1,k)=ph(i,js,k)
             ph(i,jep1,k)=ph(i,je,k)
           enddo
         enddo
         do j=jsm1,jep1
           do i=ism1,iep1
             ph(i,j,ksm1)=ph(i,j,ks)
             ph(i,j,kep1)=ph(i,j,ke)
           enddo
         enddo
!
!        checking
!
         rhserr=0.d0
         rhsnorm=0.d0
         do k=ks,ke
           do j=js,je
             do i=is,ie
               qq(i-is+1,j-js+1,k-ks+1) = ath(j-js+1)*ph(i,j-1,k) &
                 +ar(i-is+1)*ph(i-1,j,k)+(bth(j-js+1)+br(i-is+1)) &
                 *ph(i,j,k)+cr(i-is+1)*ph(i+1,j,k) &
                 +cth(j-js+1)*ph(i,j+1,k) &
                 +g32bi(j)*g32bi(j)*dx3ai(k) &
                 *(dx3bi(k+1)*(ph(i,j,k+1)-ph(i,j,k)) &
                  -dx3bi(k)*(ph(i,j,k)-ph(i,j,k-1)))
               qnorm = dabs(ath(j-js+1)*ph(i,j-1,k)) &
                 +dabs(ar(i-is+1)*ph(i-1,j,k)) &
                 +dabs((bth(j-js+1)+br(i-is+1)) &
                 *ph(i,j,k))+dabs(cr(i-is+1)*ph(i+1,j,k)) &
                 +dabs(cth(j-js+1)*ph(i,j+1,k)) &
                 +dabs(g32bi(j)*g32bi(j)*dx3ai(k) &
                 *(dx3bi(k+1)*(ph(i,j,k+1)-ph(i,j,k)) &
                  -dx3bi(k)*(ph(i,j,k)-ph(i,j,k-1))))
               rhserr=rhserr+dabs(qq(i-is+1,j-js+1,k-ks+1))
               rhsnorm=rhsnorm+qnorm
             enddo
           enddo
         enddo
         write(6,*) 'rhserr,rhsnorm= ',rhserr,rhsnorm
!        open(unit=13, file=trim(qqfile), &
!          form='unformatted',access='stream')
!        write(13) in,jn,kn
!        do k=1,ke-ks+1
!        write(13) ((qq(i,j,k),i=1,ie-is+1),j=1,je-js+1)
!        enddo
!        close(13)
!
!        write out phfile
!
         time=0.D0
         unit_b=0.20000000000000D+02
         write(6,*) 'write ph file'
         call flush(6)
         open(unit=13, file=trim(phfile), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm1,kep1
         write(13) ((ph(i,j,k),i=ism1,iep1),j=jsm1,jep1)
         enddo
         close(13)
!
!        compute bp1,bp2,bp3 from ph
!
         do k=ksm1,kep1
           do j=jsm1,jep1
             do i=is,iep1
               bp1(i,j,k)=-(ph(i,j,k)-ph(i-1,j,k))*dx1bi(i)
             enddo
           enddo
         enddo
         do k=ksm1,kep1
           do j=jsm1,jep1
             bp1(ism2,j,k)=bp1(is,j,k)
             bp1(ism1,j,k)=bp1(is,j,k)
             bp1(iep2,j,k)=bp1(iep1,j,k)
             bp1(iep3,j,k)=bp1(iep2,j,k)
           enddo
         enddo
         do j=jsm1,jep1
           do i=ism2,iep3
             bp1(i,j,ksm2)=bp1(i,j,ks+1)
             bp1(i,j,kep2)=bp1(i,j,ke-1)
           enddo
         enddo
         do k=ksm2,kep2
           do i=ism2,iep3
             bp1(i,jsm2,k)=bp1(i,js+1,k)
             bp1(i,jep2,k)=bp1(i,je-1,k)
           enddo
         enddo
!
         write(6,*) 'write bp1'
         call flush(6)
         open(unit=13, file=trim(b1file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
!        write(13) ((bp1(i,j,k)/unit_b,i=ism2,iep3),j=jsm2,jep2)
         write(13) ((bp1(i,j,k)/unit_b,i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
!
         do k=ksm1,kep1
           do j=js,jep1
             do i=ism1,iep1
               bp1(i,j,k)=-(ph(i,j,k)-ph(i,j-1,k)) &
                         *g2bi(i)*dx2bi(j)
             enddo
           enddo
         enddo
         do k=ksm1,kep1
           do j=js,jep1
             bp1(ism2,j,k)=bp1(ism1,j,k)
             bp1(iep2,j,k)=bp1(iep1,j,k)
           enddo
         enddo
         do k=ksm1,kep1
           do i=ism2,iep2
             bp1(i,jsm1,k)=-bp1(i,js+1,k)
             bp1(i,jsm2,k)=-bp1(i,js+2,k)
             bp1(i,jep2,k)=-bp1(i,je,k)
             bp1(i,jep3,k)=-bp1(i,je-1,k)
           enddo
         enddo
         do j=jsm2,jep3
           do i=ism2,iep2
             bp1(i,j,ksm2)=bp1(i,j,ks+1)
             bp1(i,j,kep2)=bp1(i,j,ke-1)
           enddo
         enddo
!
         write(6,*) 'write bp2'
         open(unit=13, file=trim(b2file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
!        write(13) ((bp1(i,j,k)/unit_b,i=ism2,iep2),j=jsm2,jep3)
         write(13) ((bp1(i,j,k)/unit_b,i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
!
         do k=ks,kep1
           do j=jsm1,jep1
             do i=ism1,iep1
               bp1(i,j,k)=-(ph(i,j,k)-ph(i,j,k-1)) &
                         *g31bi(i)*g32bi(j)*dx3bi(k)
             enddo
           enddo
         enddo
         do k=ks,kep1
           do j=jsm1,jep1
             bp1(ism2,j,k)=bp1(ism1,j,k)
             bp1(iep2,j,k)=bp1(iep1,j,k)
           enddo
         enddo
         do j=jsm1,jep1
           do i=ism2,iep2
             bp1(i,j,ksm1)=-bp1(i,j,ks+1)
             bp1(i,j,ksm2)=-bp1(i,j,ks+2)
             bp1(i,j,kep2)=-bp1(i,j,ke)
             bp1(i,j,kep3)=-bp1(i,j,ke-1)
           enddo
         enddo
         do k=ksm2,kep3
           do i=ism2,iep2
             bp1(i,jsm2,k)=bp1(i,js+1,k)*g32b(js+1)*g32bi(jsm2)
             bp1(i,jep2,k)=bp1(i,je-1,k)*g32b(je-1)*g32bi(jep2)
           enddo
         enddo
!
         write(6,*) 'write bp3'
         open(unit=13, file=trim(b3file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
!        do k=ksm2,kep3
         do k=ksm2,kep2
         write(13) ((bp1(i,j,k)/unit_b,i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
!
!        Physical parameters and units specified in input_param
!
         g_const = 6.672D-8
         rgas = 8.3D7
         kboltz = 1.381D-16
         mproton = 1.673D-24
         mpovme = 1836.15267245D0
         muconst = 0.5D0
         gamma=1.666666667D0
         unit_rho=0.83650000000000D-15
         unit_len=0.69600000000000D+11
         unit_temp=0.10000000000000D+07
         unit_v=unit_b/dsqrt(4.D0*pi*unit_rho)
         unit_time=unit_len/unit_v

         msol=1.99D33
         rsol=6.96D10
         gacc_s=g_const*msol/rsol/rsol
         gacc=gacc_s/(unit_v*unit_v/unit_len)*(rsol/unit_len)**2
         temp_c=2.D4
         rho_c=1.D12*mproton
         rstar=kboltz/(muconst*mproton)
         as2_c=rstar*temp_c
         hp_c=as2_c/gacc_s

         as2=as2_c/unit_v/unit_v
         ovhp0=unit_len/hp_c
         d0=rho_c/unit_rho

         allocate(fint(inmax))
         fint(is)=0.D0
         fint(ism1)=fint(is) &
           -ovhp0*(x1a(is)/x1b(is))**2 &
           /tovt0(x1b(is)/x1a(is))*dx1a(is)*0.5D0 &
           -ovhp0*(x1a(is)/x1b(ism1))**2 &
           /tovt0(x1b(ism1)/x1a(is))*dx1a(ism1)*0.5D0
         fint(ism2)=fint(ism1) &
           -ovhp0*(x1a(is)/x1b(ism1))**2 &
           /tovt0(x1b(ism1)/x1a(is))*dx1a(ism1)*0.5D0 &
           -ovhp0*(x1a(is)/x1b(ism2))**2 &
           /tovt0(x1b(ism2)/x1a(is))*dx1a(ism2)*0.5D0
!
!      output physical constants and units
!
         open(unit=18,file='physparams.dat', &
           form='unformatted',access='stream')
         write(18) g_const
         write(18) rgas
         write(18) kboltz
         write(18) mproton
         write(18) pi
         write(18) gamma
         write(18) unit_rho
         write(18) unit_len
         write(18) unit_b
         write(18) unit_temp
         write(18) unit_v
         write(18) unit_time
         write(18) msol
         write(18) rsol
         write(18) temp_c
         write(18) rho_c
         write(18) muconst
         close(18)
!
!        initialize e field
!
       do i=is+1,inmax-1
         fint(i)=fint(i-1) &
          +ovhp0*(x1a(is)/x1b(i-1))**2 &
          /tovt0(x1b(i-1)/x1a(is))*dx1a(i-1)*0.5D0 &
          +ovhp0*(x1a(is)/x1b(i))**2 &
          /tovt0(x1b(i)/x1a(is))*dx1a(i)*0.5D0
       enddo
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             bp1(i,j,k) = 1.D0/(gamma-1.D0) &
               *as2*d0*dexp(-fint(i))
           enddo
         enddo
       enddo
!
!      initialize d field
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             ph(i,j,k)=d0/tovt0(x1b(i)/x1a(is)) &
                      *dexp(-fint(i))
             tfield(i,j,k)=bp1(i,j,k)*(gamma-1.D0) &
                /ph(i,j,k)*unit_v**2/(rstar*unit_temp)
           enddo
         enddo
       enddo
!
!
!      lower outer pressure
!
       r1=x1b(inmax-3-10)
       r2=x1b(inmax-3)
       pr1 = as2*d0*dexp(-fint(inmax-3-10))
       pr2=  as2*d0*dexp(-fint(inmax-3))*0.2D0
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,ie
             if(x1b(i) .gt. r1) then
               bp1(i,j,k) = 1.D0/(gamma-1.D0) &
                 *(pr1+(x1b(i)-r1)/(r2-r1)*(pr2-pr1))
               ph(i,j,k) = bp1(i,j,k)*(gamma-1.D0) &
                 /tfield(i,j,k)*unit_v**2/(rstar*unit_temp)
             endif
           enddo
         enddo
       enddo

         do k=ksm2,kep2
           do j=jsm2,jep2
             bp1(iep1,j,k)=bp1(ie,j,k)
             bp1(iep2,j,k)=bp1(ie,j,k)
             ph(iep1,j,k)=ph(ie,j,k)
             ph(iep2,j,k)=ph(ie,j,k)
           enddo
         enddo
         open(unit=13, file=trim(efile), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
         write(13) ((bp1(i,j,k),i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)

         open(unit=13, file=trim(dfile), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
         write(13) ((ph(i,j,k),i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
!
!        initialize v1,v2,v3
!
         do k=ksm2,kep2
           do j=jsm2,jep2
             do i=ism2,iep2
               bp1(i,j,k)=0.D0
             enddo
           enddo
         enddo
         open(unit=13, file=trim(v1file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
         write(13) ((bp1(i,j,k),i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
         open(unit=13, file=trim(v2file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
         write(13) ((bp1(i,j,k),i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
         open(unit=13, file=trim(v3file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
         write(13) ((bp1(i,j,k),i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
!
!        initialize q flux file
!
         do k=ksm2,kep2
           do j=jsm2,jep2
             do i=ism2,iep2
               bp1(i,j,k)=0.D0
             enddo
           enddo
         enddo
         open(unit=13, file=trim(q1file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
         write(13) ((bp1(i,j,k),i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
         open(unit=13, file=trim(q2file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
         write(13) ((bp1(i,j,k),i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
         open(unit=13, file=trim(q3file), &
           form='unformatted',access='stream')
         write(13) time
         write(13) in,jn,kn
         do k=ksm2,kep2
         write(13) ((bp1(i,j,k),i=ism2,iep2),j=jsm2,jep2)
         enddo
         close(13)
!
!        output grid
!
         open(unit=13, file=trim(gridfile_out), &
           form='unformatted',access='stream')
         write(13) in,jn,kn
         write(13) (x1a(i),i=ism2,iep3)
         write(13) (x1b(i),i=ism2,iep2)
         write(13) (x2a(j),j=jsm2,jep3)
         write(13) (x2b(j),j=jsm2,jep2)
         write(13) (x3a(k),k=ksm2,kep3)
         write(13) (x3b(k),k=ksm2,kep2)

         write(13) (dx1a(i),i=ism2,iep2)
         write(13) (dx1b(i),i=ism2,iep3)
         write(13) (dx2a(j),j=jsm2,jep2)
         write(13) (dx2b(j),j=jsm2,jep3)
         write(13) (dx3a(k),k=ksm2,kep2)
         write(13) (dx3b(k),k=ksm2,kep3)

         write(13) (g2a(i),i=ism2,iep3)
         write(13) (g2b(i),i=ism2,iep2)
         write(13) (g31a(i),i=ism2,iep3)
         write(13) (g31b(i),i=ism2,iep2)
         write(13) (g32a(j),j=jsm2,jep3)
         write(13) (g32b(j),j=jsm2,jep2)

         write(13) (dg2bd1(i),i=ism2,iep3)
         write(13) (dg2ad1(i),i=ism2,iep2)
         write(13) (dg31bd1(i),i=ism2,iep3)
         write(13) (dg31ad1(i),i=ism2,iep2)
         write(13) (dg32bd2(j),j=jsm2,jep3)
         write(13) (dg32ad2(j),j=jsm2,jep2)

         write(13) (dvl1a(i),i=ism2,iep2)
         write(13) (dvl1b(i),i=ism2,iep3)
         write(13) (dvl2a(j),j=jsm2,jep2)
         write(13) (dvl2b(j),j=jsm2,jep3)
         write(13) (dvl3a(k),k=ksm2,kep2)
         write(13) (dvl3b(k),k=ksm2,kep3)

         write(13) (dx1ai(i),i=ism2,iep2)
         write(13) (dx1bi(i),i=ism2,iep3)
         write(13) (dx2ai(j),j=jsm2,jep2)
         write(13) (dx2bi(j),j=jsm2,jep3)
         write(13) (dx3ai(k),k=ksm2,kep2)
         write(13) (dx3bi(k),k=ksm2,kep3)

         write(13) (g2ai(i),i=ism2,iep3)
         write(13) (g2bi(i),i=ism2,iep2)
         write(13) (g31ai(i),i=ism2,iep3)
         write(13) (g31bi(i),i=ism2,iep2)
         write(13) (g32ai(j),j=jsm2,jep3)
         write(13) (g32bi(j),j=jsm2,jep2)

         write(13) (dvl1ai(i),i=ism2,iep2)
         write(13) (dvl1bi(i),i=ism2,iep3)
         write(13) (dvl2ai(j),j=jsm2,jep2)
         write(13) (dvl2bi(j),j=jsm2,jep3)
         write(13) (dvl3ai(k),k=ksm2,kep2)
         write(13) (dvl3bi(k),k=ksm2,kep3)
         close(13)
!
       stop
end Program main

Function tovt0(rovrs)
  Use problemsize
  Use grid
  implicit none

  real(8) :: tovt0
  real(8), intent(in) :: rovrs
  real(8), parameter :: rbovrs=1.D0+3.D8/6.96D10, &
                        rtovrs=1.D0+1.D9/6.96D10, &
                        tjump=74.D0, tbot=1.D0
!---------------------------------------------------------------------------
       if (rovrs .ge. rtovrs) then
         tovt0=tbot+tjump
       else
         if (rovrs .ge. rbovrs) then
           tovt0=tbot+(rovrs-rbovrs)/(rtovrs-rbovrs)*tjump
         else
           tovt0=tbot
         endif
       endif

end Function tovt0
