      double precision function pdg2pdf_timed(ih,ipdg,ibeam,x,xmu)
c        function
         double precision pdg2pdf
         external pdg2pdf

c        argument

         integer ih, ipdg, ibeam
         DOUBLE  PRECISION x,xmu

c timing statistics
         include "timing_variables.inc"

         call cpu_time(tbefore)
         pdg2pdf_timed = pdg2pdf(ih,ipdg,ibeam,x,xmu)
         call cpu_time(tAfter)
         tPDF = tPDF + (tAfter-tBefore)
         return

      end

      double precision function pdg2pdf(ih,ipdg,ibeam,x,xmu)
c***************************************************************************
c     Based on pdf.f, wrapper for calling the pdf of MCFM
c     ih is now signed <0 for antiparticles
c     if ih<0 does not have a dedicated pdf, then the one for ih>0 will be called
c     and the sign of ipdg flipped accordingly.
c
c     ibeam is the beam identity 1/2
c      if set to -1/-2 it meand that ipdg should not be flipped even if ih<0
c      usefull for re-weighting
c***************************************************************************
      implicit none
c
c     Arguments
c
      DOUBLE  PRECISION x,xmu
      INTEGER IH,ipdg, ibeam
C
C     Include
C
      include 'pdf.inc'
C      
      integer i,j,ihlast(20),ipart,iporg,ireuse,imemlast(20),iset,imem
     &     ,i_replace,ii,ipartlast(20),isetlast(20)
      double precision xlast(20),xmulast(20),pdflast(20)
      save ihlast,xlast,xmulast,pdflast,imemlast,ipartlast
      data ihlast/20*-99/
      data ipartlast/20*-99/
      data xlast/20*-99d9/
      data xmulast/20*-99d9/
      data pdflast/20*-99d9/
      data imemlast/20*-99/
      data isetlast/20*-99/
      data i_replace/20/

      if (ih.eq.0) then
c     Lepton collisions (no PDF). 
         pdg2pdf=1d0
         return
      endif

c     Make sure we have a reasonable Bjorken x. Note that even though
c     x=0 is not reasonable, we prefer to simply return pdg2pdf=0
c     instead of stopping the code, as this might accidentally happen.
      if (x.eq.0d0) then
         pdg2pdf=0d0
         return
      elseif (x.lt.0d0 .or. x.gt.1d0) then
         write (*,*) 'PDF not supported for Bjorken x ', x
         stop 1
      endif
      if (ibeam.gt.0)then
         ipart=sign(1,ih)*ipdg
      else
         ipart=ipdg
      endif
      
      if(iabs(ipart).eq.21) then
         ipart=0
      else if(iabs(ipart).eq.22) then
         ipart=7
      else if(iabs(ipart).eq.7) then
         ipart=7
      else if(iabs(ipart).gt.7)then
c     This will be called for any PDG code, but we only support up to 7
C         write(*,*) 'PDF not supported for pdg ',ipdg
C         write(*,*) 'For lepton colliders, please set the lpp* '//
C     $    'variables to 0 in the run_card'  
C         open(unit=26,file='../../../error',status='unknown')
C         write(26,*) 'Error: PDF not supported for pdg ',ipdg
C         stop 1
         pdg2pdf=0d0
         return
      endif

c     Determine the iset used in lhapdf
      call getnset(iset)
c     Determine the member of the set (function of lhapdf)
      call getnmem(iset,imem)

      iporg=ipart
      ireuse = 0
      ii=i_replace
      do i=1,20
c     Check if result can be reused since any of last twenty
c     calls. Start checking with the last call and move back in time
         if (ih.eq.ihlast(ii)) then
            if (ipart.eq.ipartlast(ii)) then
               if (x.eq.xlast(ii)) then
                  if (xmu.eq.xmulast(ii)) then
                     if (imem.eq.imemlast(ii)) then
                        if (iset.eq.isetlast(ii)) then
                           ireuse = ii
                           exit
                        endif
                     endif
                  endif
               endif
            endif
         endif
         ii=ii-1
         if (ii.eq.0) ii=ii+20
      enddo

c     Reuse previous result, if possible
      if (ireuse.gt.0) then
         if (pdflast(ireuse).ne.-99d9) then
            pdg2pdf=pdflast(ireuse)
            return 
         endif
      endif

c Calculated a new value: replace the value computed longest ago
      i_replace=mod(i_replace,20)+1

c     Call lhapdf and give the current values to the arrays that should
c     be saved
      call evolvepartm(iset,ipart,x,xmu,pdg2pdf)

c     Ansatz                                   
      pdg2pdf=pdg2pdf/x                                                             
c      pdg2pdf=(1d0-x)/x

c$$$c     Turning off quarks                                                                        
c$$$      if(ipart.ne.0) then
c$$$	 pdg2pdf=0d0
c$$$      endif
c
c     Turning off gluons                                  
      if(ipart.eq.0) then
	 pdg2pdf=0d0
      endif

      pdflast(i_replace)=pdg2pdf
      xlast(i_replace)=x
      xmulast(i_replace)=xmu
      ihlast(i_replace)=ih
      imemlast(i_replace)=imem
      isetlast(i_replace)=iset


      
      return
      end

      double precision function get_ee_expo()
      ! return the exponent used in the
      ! importance-sampling transformation to sample
      ! the Bjorken x's
      implicit none
      stop 21
      return
      end

      double precision function compute_eepdf(x,omx_ee, xmu, n_ee, id, idbeam)
      implicit none
      double precision x, xmu, omx_ee(*)
      integer n_ee, id, idbeam
      stop 21
      return
      end

      double precision function ee_comp_prod(comp1, comp2)
      ! compute the scalar product for the two array
      ! of eepdf components
      implicit none
      double precision comp1(*), comp2(*)
      stop 21
      return
      end
      
      
