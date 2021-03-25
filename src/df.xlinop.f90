  subroutine dfxlinop(j,n19)
! --- calculations of more complicated line profiles, i.e. for hydrogen lines, and Stark profile lines.

    use types
    use dfcomm
    use hprof
    
    implicit none

    include 'common.rhoxbl'
    include 'common.tempbl' 
    include 'common.turbpr'
    include 'common.stateb'
    include 'common.depart'
    include 'common.opsblk'
    include 'common.ionsbl'

      integer, intent(in) :: j, n19
!---- Local variables -------------------------------------------!      

      integer n, k,  iline,  itemp1
      integer(kind=4) nlast,  ifstrong, minred, maxred, ibuff, maxblue, i, ib
!these are read as integers, but processed later as real

      integer(kind=4)  nblo, nbup, nelion, type, ncon, nelionx, nbuff
      real(kind=4) gammar, gammas, gammaw
      real(kind=4) bshore,  g,  cgf, ashore
      real(kind=4) gf

      equivalence (gammas,ashore),(gammaw,bshore)
      equivalence (gf,g,cgf),(type,nlast)!,(gammar,gaunt)


      real(kind=4) nstark,ndopp,nmerge,inglis
      real(kind=4) kappa,kapmin,kappa0,kapcen
      real(kind=4) kappa0red,kappared,kappa0blue,kappablue
      real(kind=4) oldkappa,newkappa,dkappa
      real(kind=8) dopph(maxd)


      real(kind=dp)  wave,wcon,wmerge,wshift,contx(26,17),wtail
      real(kind=dp)  wl,emergeh(maxd)
      real(kind=dp)  bluecut,wlplus1,wlplus2,redcut,wlminus1,wlminus2,vacair
      real(kind=dp)  ehyd(100),conth(15),alphahyd(99)
!     default implicit rule
      real(kind=dp) frelin, freq, epsil,  edgeblue, adamp, dopwl, vvoigt 
      integer(kind=4)  nbuff1, nbuff2, nbuff3, ixwl 
      real(kind=dp) dnbuff, xsectg, tail
!---------- External ----------------------------------------------!
      real(kind=dp) voigt !
   

      logical whilep, anotherp

  conth = (/ 109678.764,  27419.659,  12186.462,  6854.871,  4387.113,  &  !1.00
     &         3046.604,   2238.320,   1713.711,  1354.044,  1096.776,  & 
     &          906.426,    761.650,    648.980,   559.579,   487.456/) 


  contx= 0.0d0
  contx(1:10,1) = (/109678.764,  27419.659,  12186.462,  6854.871,  4387.113,  &   !1.00
     &         3046.604,   2238.320,   1713.711,  1354.044,  1096.776/)
  contx(1:6,2) =(/198310.760,  38454.691,  32033.214,  29223.753,27175.760,15073.868/) !2.00
  contx(1:6,3) =(/438908.850, 109726.529,  48766.491,  27430.925,17555.715,12191.437/) !2.01
  contx(1:10,4) =(/90883.840,  90867.420,  90840.420,  90820.420,90804.000,  &   !6.00
     &        90777.000,  80691.180,  80627.760,  69235.820,69172.400/)
  contx(1:4,6)= (/61671.020,  39820.615,  39800.556,  39759.842/) !12.00
  contx(1:2,8) =(/48278.370,  48166.309/) !13.00
  contx(1:10,10) = (/66035.000,  65957.885,  65811.843,  65747.550,65670.435,  &   !14.00
     &        65524.393,59736.150,59448.700,50640.630,50553.180/) 

   itemp1 = 0
!
!


  if(itemp.ne.itemp1) then 

    ehyd(1)=0.
    ehyd(2)=82259.105
    ehyd(3)=97492.302
    ehyd(4)=102823.893
    ehyd(5)=105291.651
    ehyd(6)=106632.160
    ehyd(7)=107440.444
    ehyd(8)=107965.051
 
    do  n=9,100
     ehyd(n)=109678.764d0-109677.576d0/n**2
    end do  

    do n=1,99
      alphahyd(n)=1.d7/(ehyd(n+1)-ehyd(n))
    end do 

    do k=1,nrhox
!     empirical
      inglis=1600./xne(k)**(2./15.)
      nmerge=inglis-1.5
      emerge(k)=109737.312/nmerge**2
      emergeh(k)=109677.576/nmerge**2
    end do 
      itemp1=itemp
  end if 
 
      rewind 19
!
!--- open big loop over all ntle lines ---!
  do iline= 1, n19
  


    read(19)wl,gf,nblo,nbup,nelion,type,ncon,nelionx,gammar,gammas,gammaw,nbuff
 

 
!--- Note that He lines do not work yet! -----
!!  ---- the following types are considered:
!
!     IF(TYPE.EQ.2)    coronal lines not implemented
!
!     IF(TYPE.LT.-2) OR IF(TYPE.EQ.0) or IF(TYPE.EQ.3) -> normal line
!
!      IF(TYPE.EQ.-1)GO TO 600 Hydrogen line
!      IF(TYPE.EQ.-2)GO TO 600 Hydrogen line 
!
!      IF(TYPE.EQ.1)GO TO 700 autoionizing line 
!--------------------------------------------------------------------

    if ((type .eq. 0 ) .or. (type .lt.-2) .or. (type .eq. 3)) then  

!------     normal line ------------------------

      kappa0=cgf*xnfdop(nelion)
      if(kappa0 .ge. continuum(nbuff)) then 
      ifstrong=0
!     since cont=cont/1000, cont*500 is 1/2 the actual continuum.  lines are
!     strong if the line center is 1/2 the continuum.
      if( kappa0 .gt. 500.*continuum(nbuff)) ifstrong=1

      mlines=mlines+1
      wcon=0.
      wtail=0.


      if(ncon.gt.0)then
        wcon=1.e7/(contx(ncon,nelionx)-emerge(j))
        wtail=1.e7/(contx(ncon,nelionx)-emerge(j)-500.)
      endif

      adamp=(gammar+gammas*xne(j)+gammaw*txnxn(j))/dopple(nelion)
      kapcen=kappa0*voigt(0.0d0,adamp)
      dopwl=dopple(nelion)*wl

        if(wl .le. wlend) then 

!------     red wing
          minred=max(1,nbuff)
          maxred=min(nbuff+1000,length)
          if(ifstrong .eq. 1) maxred=length
          wave=wbegin*ratio**(minred-1)

          ibuff = minred -1 
          whilep = .true.

          do while ((ibuff .lt. maxred) .and. whilep )      
            ibuff = ibuff +1
            if(wave.ge.wcon) then 
              vvoigt=abs(wave-wl)/dopwl
              kappa=kappa0*voigt(vvoigt,adamp)
              if(wave.lt.wtail)kappa=kappa*(wave-wcon)/(wtail-wcon)
              buffer(ibuff)=buffer(ibuff)+kappa
           
              if(kappa.lt.continuum(ibuff)) then 
                whilep = .false.
              else
                wave=wave*ratio
              end if
            end if 
            
        end do 

       end if 
       whilep = .false.

       if (wl .le. wlend) then 
         if(minred.ne.1) then 
          if(wl.ge.wlbeg) then 
            whilep = .true.
          end if
         end if 
       else
          whilep = .true.
       end if

       if (whilep) then 
 
!------     blue wing
         ibuff=min(length+1,nbuff)
         maxblue=ibuff-1
         if(ifstrong.eq.0) maxblue=min(maxblue,1000)
    
         wave=wbegin*ratio**(ibuff-1)
         i = 0
         
         do while ((i .lt. maxblue) .and. (whilep)) 
            i= i+1
            ibuff=ibuff-1
            wave=wave/ratio
  
            if(wave .ge. wcon) then 
              vvoigt=abs(wave-wl)/dopwl
              kappa=kappa0*voigt(vvoigt,adamp)
              if(wave.lt.wtail)kappa=kappa*(wave-wcon)/(wtail-wcon)
              buffer(ibuff)=buffer(ibuff)+kappa

              if(kappa.lt.continuum(ibuff)) then 
                whilep = .false.
              end if 
            end if 
           
         end do !end do while 

       end if 
      end if !top if

    else if ((type .eq. -1) .or. (type .eq. -2) ) then   

!-----  hydrogen line
    
      kappa0=cgf*xnfdop(1)
      if(kappa0.ge. continuum(nbuff)) then   


        dopph(j)=dopple(1)
!--     deuterium

        if(type.eq.-2)dopph(j)=dopph(j)/1.4142
 
        mlines=mlines+1


!      if(ncon.eq.0)go to 620
!      if(nbup.eq.nblo+1)go to 620
!      if(nbup.eq.nblo+2)go to 630
!------------------------------------------------------------------------------
!----- three conditions 1) non of them are true, 2) either 1st or 2nd is true, 3) 3rd is true: 

      if((ncon .ne. 0) .and. (nbup .ne. nblo+1) .and. (nbup .ne. nblo+2)) then

 
        wshift=1.0d7/(conth(ncon)-109677.576d0/81.d0**2)
        wmerge=1.0d7/(conth(ncon)-emergeh(j))
  
        if(wmerge.lt.0.)wmerge=wshift+wshift
   
        wcon=max(wshift,wmerge)
        wtail=1.0d7/(1.0d7/wcon-500.0d0)
        wcon=min(wshift+wshift,wcon)
       
        if(wtail .lt. 0.) wtail=wcon+wcon
        wtail=min(wcon+wcon,wtail)

       
        whilep = .false.
 
        if(wl  .le. wlend) then 
!----------------------------------------------------------------------
!     computes the first 20 points, then every 20th point and log interpolates
!     the 19 missing points
!
!     red wing
!----------------------------------------------------------------------

          if(wlbeg .gt. alphahyd(nblo)) then 
            whilep = .true.
          else 
            
            redcut=1.d7/(109678.764d0-109677.576d0/(nbup-0.8d0)**2-ehyd(nblo))
            wlminus1=1.d7/(ehyd(nbup-1)-ehyd(nblo))
            wlminus2=1.d7/(ehyd(nbup-2)-ehyd(nblo))
            kappa0red=kappa0*hfnm(nblo,nbup-2)/hfnm(nblo,nbup)/ & 
     &                (ehyd(nbup-2)-ehyd(nblo))*(ehyd(nbup)-ehyd(nblo))

            minred=max(1,nbuff)
            wave=wbegin*ratio**(minred-1)
            ib=0

            ibuff = minred-1 
            whilep = .true.
 
            do while ((ibuff .lt. length) .and. whilep ) 
        
              ibuff = ibuff + 1     
              if(wave.lt.wcon) then    ! 1
              ! nothing 
              else
                if(wave.gt.wlminus1) then ! 2 
                  whilep = .false. 
                else

                  ib=ib+1
                  if((ib.gt.21) .and. (mod(ib,20) .ne. 1)) then !3 
                    wave = wave*ratio    
                  else 

                    kappa=kappa0*hprof4(nblo,nbup,j,wl,wave-wl,dopph)

                    if(wave.lt.wtail)kappa=kappa*(wave-wcon)/(wtail-wcon)

                    if(wave.gt.redcut)then
                      kappared=kappa0red*hprof4m(nblo,nbup-2,j,wl,wave-wlminus2,dopph)
                      if(wave.lt.wtail)kappared=kappared*(wave-wcon)/(wtail-wcon)
                      if(kappared.ge.kappa) whilep = .false. 
                    endif

                    if (whilep) then !4 
                      buffer(ibuff)=buffer(ibuff)+kappa
                      newkappa=log(max(kappa,1.e-20))
         
                      if(ib.gt.21)then
                         dkappa=(newkappa-oldkappa)/20.
                         do  k=1,19
                           buffer(ibuff-k)=buffer(ibuff-k)+exp(newkappa-dble(k)*dkappa)
                         end do 
                      endif

                      oldkappa=newkappa
         
                      if(kappa.lt.continuum(ibuff)) then 
                        whilep = .false.                   
                      end if 
                    end if ! 4 

                    if (whilep) wave = wave*ratio 

                  end if !3 
                end if ! 2
              end if! 1
 
            end do ! while do loop 

 
            whilep = .false. 
           end if 

        end if 
!------- this is to decide if the blue wing has to be done :       
        if (wl .le. wlend ) then 
          if (minred .ne. 1) then 
            if(wl.ge.wlbeg) whilep = .true.
          end if 
        else
          whilep = .true.
   
        end if



        if (whilep) then 
!     blue wing
    
          bluecut=1.d7/(109678.764d0-109677.576d0/(nbup+0.8d0)**2- & 
     &             ehyd(nblo))
          wlplus1=1.d7/(ehyd(nbup+1)-ehyd(nblo))
          wlplus2=1.d7/(ehyd(nbup+2)-ehyd(nblo))
          kappa0blue=kappa0*hfnm(nblo,nbup+2)/hfnm(nblo,nbup)/     &
     &            (ehyd(nbup+2)-ehyd(nblo))*(ehyd(nbup)-ehyd(nblo))
          ibuff=min(length+1,nbuff)
          maxblue=ibuff-1
          wave=wbegin*ratio**(ibuff-1)

          i = 0
          do while ((i .lt. maxblue) .and. (whilep)) !614
            i = i + 1
            ibuff=ibuff-1
            wave=wave/ratio

           if((wave.lt.wcon) .or. wave .lt. wlplus1) then ! 1 
             whilep = .false.

           else 

             if(( i .gt. 21) .and. (mod(i,20) .ne. 1))  then ! 2 
             
             else 

               kappa=kappa0*hprof4(nblo,nbup,j,wl,wave-wl,dopph)
               if(wave.lt.wtail)kappa=kappa*(wave-wcon)/(wtail-wcon)
    
               if(wave .lt. bluecut)then
                 kappablue=kappa0blue*hprof4p(nblo,nbup+2,j,wl,wave-wlplus2,dopph)
                 if(wave.lt.wtail)kappablue=kappablue*(wave-wcon)/(wtail-wcon)
                 if(kappablue.gt.kappa) whilep = .false. 
               
               endif
     
               if (whilep ) then 
                 buffer(ibuff)=buffer(ibuff)+kappa
                 newkappa=log(max(kappa,1.e-20))

                 if(i.gt.21)then
                   dkappa=(newkappa-oldkappa)/20.
        
                  do  k=1,19
                   buffer(ibuff+k)=buffer(ibuff+k)+exp(newkappa-dble(k)*dkappa)
                  end do 
                 endif

               oldkappa=newkappa
             
               if(kappa.lt.continuum(ibuff)) whilep = .false. 
               
               end if 
             end if !2 
            end if !  1
          end do !do while loop 614

        end if 

!--- second contitions 

      else if ((ncon .eq. 0) .or. (nbup .eq. nblo+1) ) then


!-------------------------------------------------------- 
!     alpha presumed isolated lines and beta blue wings
!---------------------------------------------------------
       if(wl .le. wlend) then 
!----     red wing-------!
         minred=max(1,nbuff)
         ib=0
         wave=wbegin*ratio**(minred-1)
     
         whilep = .true.
         ibuff = minred-1 

         do while ((ibuff .lt. length) .and. whilep)
          ibuff = ibuff +1
          ib=ib+1
          if((ib.gt.21 ) .and. (mod(ib,20) .ne.1)) then
            wave = wave*ratio
          else 

          kappa=kappa0*hprof4(nblo,nbup,j,wl,wave-wl,dopph)
          buffer(ibuff)=buffer(ibuff)+kappa
          newkappa=log(max(kappa,1.e-20))

          if(ib.gt.21)then
            dkappa=(newkappa-oldkappa)/20.
            do  k=1,19
             buffer(ibuff-k)=buffer(ibuff-k)+exp(newkappa-dble(k)*dkappa)
            end do  
          endif

          oldkappa=newkappa
       
            if(kappa.lt.continuum(ibuff)) then 
               whilep = .false. 
            else 
               wave = wave*ratio 
            end if
          end if 
         end do 
       end if 



       whilep = .false.

       if (wl .le. wlend) then
         if(minred.ne.1) then
          if(wl.ge.wlbeg) then
            whilep = .true.
          end if
         end if
       else
          whilep = .true.
       end if


        if (whilep ) then   
!-------- alpha or beta blue wing
          ibuff=min(length+1,nbuff)

          maxblue=ibuff-1
          wave=wbegin*ratio**(ibuff-1)

          i = 0

          do while ((i .lt. maxblue) .and. whilep) 
            i = i +1 
            ibuff=ibuff-1
            wave=wave/ratio
        
            if((i.gt.21) .and. (mod(i,20).ne. 1)) then

            else            

              kappa=kappa0*hprof4(nblo,nbup,j,wl,wave-wl,dopph)
              buffer(ibuff)=buffer(ibuff)+kappa
              newkappa=log(max(kappa,1.e-20))
              if(i.gt.21)then
                dkappa=(newkappa-oldkappa)/20.
                do k=1,19
                  buffer(ibuff+k)=buffer(ibuff+k)+exp(newkappa-dble(k)*dkappa)
                end do
              endif

              oldkappa=newkappa

              if(kappa.lt.continuum(ibuff)) then 
               whilep = .false. 
              end if  
            end if 
          end do 
        end if 

!--- third condition:

      else if (nbup .eq. nblo+2) then 
       anotherp = .false. 
!-----------------------------------------------------------
!     beta lines red wing
    
      if( wl .le. wlend) then 
! 
!-----  red wing

        minred=max(1,nbuff)
        ib=0
        wave=wbegin*ratio**(minred-1)

        ibuff = minred -1 
        whilep = .true.

        do while ((ibuff .lt. length) .and. whilep) ! 631 
          ibuff = ibuff +1 

         if(wave.gt.alphahyd(nblo)) then ! (1) 
            anotherp = .true.  
              
           if(nbuff.ge.1) then !(2) 
             ibuff=nbuff
             maxblue=ibuff-1
             wave=wbegin*ratio**(ibuff-1)

             i = 0
             do while ((i.lt.maxblue) .and. whilep  ) 
                i = i +1
                ibuff=ibuff-1
                wave=wave/ratio

                if((i.gt.21) .and. (mod(i,20).ne.1)) then 
                   whilep = .false.        
                else
    
                 kappa=kappa0*hprof4(nblo,nbup,j,wl,wave-wl,dopph)
                 buffer(ibuff)=buffer(ibuff)+kappa
                 newkappa=log(max(kappa,1.e-20))

                 if(i.gt.21)then
                  dkappa=(newkappa-oldkappa)/20.
                  do  k=1,19
                      buffer(ibuff+k)=buffer(ibuff+k)+ & 
     &                              exp(newkappa-dble(k)*dkappa)
                  end do 
                 endif

                 oldkappa=newkappa
                 if(kappa.lt.continuum(ibuff)) then 
                  whilep = .false.
                 end if 
                end if
              end do
           else 
             whilep = .false.  
           end if ! (2) 
!
           whilep = .false.  
! else from (1)
         else 
           ib=ib+1

           if((ib.gt.21) .and. mod(ib,20) .ne. 1) then !(3)  
              wave = wave*ratio 
           else 
    
             kappa=kappa0*hprof4(nblo,nbup,j,wl,wave-wl,dopph)
             buffer(ibuff)=buffer(ibuff)+kappa
             newkappa=log(max(kappa,1.e-20))

             if(ib.gt.21)then
               dkappa=(newkappa-oldkappa)/20.
               do k=1,19
                  buffer(ibuff-k)=buffer(ibuff-k)+exp(newkappa-dble(k)*dkappa)
               end do  
             endif

             oldkappa=newkappa
       
             if (kappa .lt. continuum(ibuff)) then  
               whilep = .false. 
             else 
               wave=wave*ratio
             end if 

           end if ! (close (3))  
         end if ! close(1)  
        end do
!--------------- 
       end if 

       if (wl .gt. wlend  ) then 
        whilep = .true.
       else
        if ((minred .ne. 1 ) .and. (wl .ge. wlbeg)) whilep = .true.
       end if 
       
       if (anotherp) whilep = .false.

       if (whilep ) then      
!-------- (alpha or) beta blue wing
          ibuff=min(length+1,nbuff)

          maxblue=ibuff-1
          wave=wbegin*ratio**dble(ibuff-1)

          i = 0

          do while ((i .lt. maxblue) .and. whilep)
            i = i + 1
            ibuff=ibuff-1
            wave=wave/ratio

            if((i.gt.21) .and. (mod(i,20).ne. 1)) then
            !  i = i +1
            else

              kappa=kappa0*hprof4(nblo,nbup,j,wl,wave-wl,dopph)
              buffer(ibuff)=buffer(ibuff)+kappa
              newkappa=log(max(kappa,1.e-20))

              if(i.gt.21)then
                dkappa=(newkappa-oldkappa)/20.
                do k=1,19
                  buffer(ibuff+k)=buffer(ibuff+k)+exp(newkappa-dble(k)*dkappa)
                end do
              endif

              oldkappa=newkappa

              if(kappa.lt.continuum(ibuff)) then
               whilep = .false.
              end if 
            end if
          end do

      end if 


      
      end if ! three conditions

      end if ! -- top if 


!--------------------------------------------------------------!
     else if (type .eq. 1) then 
!-------------------------------
!     autoionizing line
!-------------------------------    

       kappa0=bshore*g*xnfpel(nelion)
       if(kappa0.ge.continuum(nbuff)) then 

         mlines=mlines+1
         frelin=2.99792458e17/wl

         if(wl .le. wlend) then  
!          red wing

           minred=max(1,nbuff)
           freq=2.99792458e17/(wbegin*ratio**(minred-1))

           ibuff = minred-1 
           whilep = .true.
           do while (whilep .and. (ibuff .lt. length) )  
              ibuff = ibuff +1 

              epsil=2.*(freq-frelin)/gammar
              kappa=kappa0*(ashore*epsil+bshore)/(epsil**2+1.)/bshore
              buffer(ibuff)=buffer(ibuff)+kappa
              freq=freq/ratio

              if (kappa .lt. continuum(ibuff)) then 
                whilep = .false.
              end if 

           end do 
         end if 

         whilep= .false. 
         
         if (wl .le. wlend) then 
           if (nbuff .ne. 1 ) then  
             whilep = .true. 
           end if 
         else 
           whilep = .true.
         end if  

         if (whilep ) then  
!          blue wing
           ibuff=min(length+1,nbuff)
           maxblue=ibuff-1
           freq=2.99792458e17/(wbegin*ratio**(ibuff-1))
           
           i = 0 
           whilep = .true. 
           do while ((i .lt. maxblue) .and. whilep) 
             i = i +1
             ibuff=ibuff-1
             freq=freq*ratio
             epsil=2.*(freq-frelin)/gammar
             kappa=kappa0*(ashore*epsil+bshore)/(epsil**2+1.)/bshore
             buffer(ibuff)=buffer(ibuff)+kappa
            if(kappa.lt.continuum(ibuff)) then 
             whilep = .false.
            end if 
           end do  


         end if
       end if

!----------------------------------------------------------------!
    else if (type .eq. 2 ) then 
!------------------------------
!   coronal line
!   500 go to 900
!----------------------------------------------------------------! 
    else  
!----------------------------------
!     merged continuum
!     edge wavelengths are in vacuum
!-----------------------------------
      wshift=1.d7/(1.d7/wl-109737.312d0/nlast**2)
      wmerge=1.d7/(1.d7/wl-emerge(j))

      if(nelion.eq.1)then
        wshift=1.d7/(1.d7/wl-109677.576d0/nlast**2)
        wmerge=1.d7/(1.d7/wl-emergeh(j))
      endif
  
      if(wmerge .lt. 0.0d0) wmerge=wshift+wshift
 
      wmerge=max(wmerge,wshift)
      wmerge=min(wshift+wshift,wmerge)
      wtail=1.d7/(1.d7/wmerge-500.)

      if(wtail .lt. 0.0d0) wtail=wmerge+wmerge

      wtail=min(wmerge+wmerge,wtail)
      ixwl= log(wl)/ratiolg
      edgeblue=exp(ixwl*ratiolg)

      if(edgeblue .gt. wl) ixwl=ixwl-1

      nbuff1=ixwl+1-ixwlbeg+1
      ixwl= log(wmerge)/ratiolg+.5
      nbuff2=ixwl-ixwlbeg+1
      ixwl= log(wtail)/ratiolg+.5
      nbuff3=ixwl-ixwlbeg+1

      if(nbuff1.le.length) then 
       
        if(nbuff3.ge.1) then 
      
          dnbuff=dble(nbuff3-nbuff2) 
          nbuff1=max(nbuff1,1)
          xsectg=dble(gf) 
          kappa=xsectg*xnfpel(nelion)
          tail=1.
        
         do ibuff=nbuff1, min(nbuff3,length)
          if(ibuff.gt.nbuff2) tail=(nbuff3-ibuff)/dnbuff
          buffer(ibuff) = buffer(ibuff)+kappa*tail
         end do
        end if
      end if 
    end if 
 
!
   end do  ! 900 loop

  end 




  double precision function voigt(v,a)
!
!
      use types
      use dfcomm

      implicit none
      real(kind=dp) v, a
      real(kind=dp) aa, vv, u, aau, vvu 
      integer iv, ia     
!     fast voigt


      if(a .ge. 0.20d0) then 

        if((a .gt. 1.4d0) .or. ((a+v) .gt. 3.2d0)) then 

        aa=a*a
        vv=v*v
        u=(aa+vv)*1.4142d0 
        voigt=a*.79788d0/u
        if(a.gt. 100.0d0)return
        aau=aa/u
        vvu=vv/u
        voigt=voigt* & 
     &  ((((aau-10.*vvu)*aau+5.*vvu*vvu)/u+vvu-aau/3.)*3./u+1.)

    
        return
        else 

          iv=v*200.+1.5
          ia=a*200.+1.5
          voigt=atab(ia,iv)

        end if 

      else

        if(v .le. 10.) then 
          iv=v*200.+1.5
          voigt=(h2tab(iv)*a+h1tab(iv))*a+h0tab(iv)
        else
          voigt=0.5642d0*a/(v**2) 
        end if 

      end if

      return
  end


