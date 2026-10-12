!WRF:MODEL_LAYER:DYNAMICS
!

!  SMALL_STEP code for the geometric height coordinate model
!
!---------------------------------------------------------------------------

MODULE module_small_step_em

   USE module_configure
   USE module_model_constants

CONTAINS

!----------------------------------------------------------------------

SUBROUTINE small_step_prep( u_1, u_2, v_1, v_2, w_1, w_2, &
                            t_1, t_2, ph_1, ph_2,         &
                            mub, mu_1, mu_2,              &
                            muu, muus, muv, muvs,         &
                            mut, muts, mudf,              &
                            u_save, v_save, w_save,       &
                            t_save, ph_save, mu_save,     &
                            ww, ww_save,  change_cs, mub_rk1, rebalance_rk,        & !PZ added move
                            c2a, pb, p, alt,              &
                            msfux, msfuy, msfvx,          &
                            msfvx_inv,                    &
                            msfvy, msftx, msfty,          &
                            rdx, rdy,                     &
                            rk_step,                      &
                            ids,ide, jds,jde, kds,kde,    & 
                            ims,ime, jms,jme, kms,kme,    &
                            its,ite, jts,jte, kts,kte    )

  IMPLICIT NONE  ! religion first

! declarations for the stuff coming in

  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  INTEGER,      INTENT(IN   )    :: rk_step

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT) :: u_1,   &
                                                              v_1,   &
                                                              w_1,   &
                                                              t_1,   &
                                                              ph_1

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(  OUT) :: u_save,   &
                                                              v_save,   &
                                                              w_save,   &
                                                              t_save,   &
                                                              ph_save

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT) :: u_2,   &
                                                              v_2,   &
                                                              w_2,   &
                                                              t_2,   &
                                                              ph_2

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(  OUT) :: c2a, &
                                                               ww_save

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::  pb,  &
                                                                p,   &
                                                                alt, &
                                                                ww

  REAL, DIMENSION(ims:ime, jms:jme)         , INTENT(INOUT) :: mu_1,mu_2

  REAL, DIMENSION(ims:ime, jms:jme)         , INTENT(IN   ) :: mub,  &
                                                               muu,  &
                                                               muv,  &
                                                               mut,  &
                                                               msfux,&
                                                               msfuy,&
                                                               msfvx,&
                                                               msfvx_inv,&
                                                               msfvy,&
                                                               msftx,&
                                                               msfty,&
                                                               mub_rk1 !PZ added

  REAL, DIMENSION(ims:ime, jms:jme)         , INTENT(  OUT) :: muus, &
                                                               muvs, &
                                                               muts, &
                                                               mudf

  REAL, DIMENSION(ims:ime, jms:jme)         , INTENT(  OUT) :: mu_save

  REAL, INTENT(IN) :: rdx,rdy
  LOGICAL, INTENT(IN):: change_cs, rebalance_rk 

! local variables

  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  INTEGER :: i_endu, j_endv


!<DESCRIPTION>
!
!  small_step_prep prepares the prognostic variables for the small timestep.
!  This includes switching time-levels in the arrays and computing coupled 
!  perturbation variables for the small timestep 
!  (i.e. mu*u" = mu(t)*u(t)-mu(*)*u(*); mu*u" is advanced during the small 
!  timesteps
!
!</DESCRIPTION>

    i_start = its
    i_end   = min(ite,ide-1)
    j_start = jts
    j_end   = min(jte,jde-1)
    k_start = kts
    k_end = min(kte,kde-1)

    i_endu = ite
    j_endv = jte

    !  if this is the first RK step, reset *_1 to *_2
    !  (we are replacing the t-dt fields with the time t fields)

    IF ((rk_step == 1) ) THEN

      DO j=j_start, j_end
      DO i=i_start, i_end
        mu_1(i,j)=mu_2(i,j)
        ww_save(i,kde,j) = 0.
        ww_save(i,1,j) = 0.
        mudf(i,j) = 0.  !  initialize external mode div damp to zero
      ENDDO
      ENDDO

      DO j=j_start, j_end
      DO k=k_start, k_end
      DO i=i_start, i_endu
        u_1(i,k,j) = u_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO

      DO j=j_start, j_endv
      DO k=k_start, k_end
      DO i=i_start, i_end
        v_1(i,k,j) = v_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO

      DO j=j_start, j_end
      DO k=k_start, k_end
      DO i=i_start, i_end
        t_1(i,k,j) = t_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO

      DO j=j_start, j_end
      DO k=k_start, min(kde,kte)
      DO i=i_start, i_end
        w_1(i,k,j)  = w_2(i,k,j)
        ph_1(i,k,j) = ph_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO

    DO j=j_start, j_end
      DO i=i_start, i_end
        muts(i,j)=mub(i,j)+mu_2(i,j)
      ENDDO
      DO i=i_start, i_endu
!  rk_step==1, WCS fix for tiling
!        muus(i,j)=0.5*(mub(i,j)+mu_2(i,j)+mub(i-1,j)+mu_2(i-1,j))
        muus(i,j) = muu(i,j)
      ENDDO
    ENDDO

    DO j=j_start, j_endv
    DO i=i_start, i_end
!  rk_step==1, WCS fix for tiling
!      muvs(i,j)=0.5*(mub(i,j)+mu_2(i,j)+mub(i,j-1)+mu_2(i,j-1))
        muvs(i,j) = muv(i,j)
    ENDDO
    ENDDO

    DO j=j_start, j_end
      DO i=i_start, i_end
        mu_save(i,j)=mu_2(i,j)
        mu_2(i,j)=0.
!       mu_2(i,j)=mu_2(i,j)-mu_2(i,j)
      ENDDO
    ENDDO

    ELSE  !rk_step >1

    DO j=j_start, j_end
      DO i=i_start, i_end
        if (rebalance_rk) then !PZ added
          muts(i,j) = mub_rk1(i,j)+mu_1(i,j)
        else
          muts(i,j)=mub(i,j)+mu_1(i,j)
        endif 
        if (i==1 .and. j==1) write(*,*) 'in small prep, muts(rk step>1=)', muts(1,1), muts(2,1)  !PZ added
      ENDDO
      DO i=i_start, i_endu
        if (rebalance_rk) then !PZ added
          muus(i,j) = 0.5*(mub_rk1(i,j)+mu_1(i,j) + mub_rk1(i-1,j)+mu_1(i-1,j))
        else
          muus(i,j)=0.5*(mub(i,j)+mu_1(i,j)+mub(i-1,j)+mu_1(i-1,j))
        endif
      ENDDO
    ENDDO

    DO j=j_start, j_endv
    DO i=i_start, i_end
      muvs(i,j)=0.5*(mub(i,j)+mu_1(i,j)+mub(i,j-1)+mu_1(i,j-1))
    ENDDO
    ENDDO

    DO j=j_start, j_end
      DO i=i_start, i_end
        mu_save(i,j)=mu_2(i,j)
        mu_2(i,j)=mu_1(i,j)-mu_2(i,j)
      ENDDO
    ENDDO


    END IF

    ! set up the small timestep variables

      DO j=j_start, j_end
      DO i=i_start, i_end
        ww_save(i,kde,j) = 0.
        ww_save(i,1,j) = 0.
      ENDDO
      ENDDO

    DO j=j_start, j_end
    DO k=k_start, k_end
    DO i=i_start, i_end
      if (change_cs) then   !PZ added for moving bottom
          c2a(i,k,j) = cpovcv*(pb(i,k,j)+p(i,k,j))/alt(i,k,j)          

          if (i==2 .and. k==1 .and. j==1) then
            write(*,*) 'In small_prep: Actual c2a=', cpovcv*(pb(i,k,j)+p(i,k,j))/alt(i,k,j), 'actual cs=', sqrt(cpovcv*(pb(i,k,j)+p(i,k,j))*alt(i,k,j)), 'c2a=', c2a(i,k,j)
          endif
      else
          c2a(i,k,j) = cpovcv*(pb(i,k,j)+p(i,k,j))/alt(i,k,j)
      ENDIF
    ENDDO
    ENDDO
    ENDDO

    DO j=j_start, j_end
    DO k=k_start, k_end
    DO i=i_start, i_endu
      u_save(i,k,j) = u_2(i,k,j)
      ! u coupled with my
      u_2(i,k,j) = (muus(i,j)*u_1(i,k,j)-muu(i,j)*u_2(i,k,j))/msfuy(i,j)
    ENDDO
    ENDDO
    ENDDO

    DO j=j_start, j_endv
    DO k=k_start, k_end
    DO i=i_start, i_end
      v_save(i,k,j) = v_2(i,k,j)
      ! v coupled with mx
!      v_2(i,k,j) = (muvs(i,j)*v_1(i,k,j)-muv(i,j)*v_2(i,k,j))/msfvx(i,j)
      v_2(i,k,j) = (muvs(i,j)*v_1(i,k,j)-muv(i,j)*v_2(i,k,j))*msfvx_inv(i,j)
    ENDDO
    ENDDO
    ENDDO

    DO j=j_start, j_end
    DO k=k_start, k_end
    DO i=i_start, i_end
      t_save(i,k,j) = t_2(i,k,j)
      t_2(i,k,j) = muts(i,j)*t_1(i,k,j)-mut(i,j)*t_2(i,k,j)
    ENDDO
    ENDDO
    ENDDO

    DO j=j_start, j_end
!    DO k=k_start, min(kde,kte)
    DO k=k_start, kde
    DO i=i_start, i_end
      w_save(i,k,j) = w_2(i,k,j)
      !PZ added for moving bottom
      !IF (j==1 .and. i<=3 .and. k==1) THEN
        !write(*,*) 'I=',i, 'w_save', w_save(i,k,j) 
      !ENDIF

      ! w coupled with my
      w_2(i,k,j)  = (muts(i,j)* w_1(i,k,j)-mut(i,j)* w_2(i,k,j))/msfty(i,j)
      ph_save(i,k,j) = ph_2(i,k,j)
      ph_2(i,k,j) = ph_1(i,k,j)-ph_2(i,k,j)
    ENDDO
    ENDDO
    ENDDO

      DO j=j_start, j_end
!      DO k=k_start, min(kde,kte)
      DO k=k_start, kde
      DO i=i_start, i_end
        ww_save(i,k,j) = ww(i,k,j)
      ENDDO
      ENDDO
      ENDDO

END SUBROUTINE small_step_prep

!-------------------------------------------------------------------------


SUBROUTINE small_step_finish( u_2, u_1, v_2, v_1, w_2, w_1,    &
                              t_2, t_1, ph_2, ph_1, ww, ww1,   &
                              mu_2, mu_1,                      &
                              mut, muts, muu, muus, muv, muvs, &
                              u_save, v_save, w_save,          &
                              t_save, ph_save, mu_save,        &
                              msfux, msfuy, msfvx, msfvy,      &
                              msftx, msfty,                    &
                              h_diabatic,                      &
                              number_of_small_timesteps,dts,   &
                              rk_step, rk_order,               &
                              ids,ide, jds,jde, kds,kde,       &
                              ims,ime, jms,jme, kms,kme,       &
                              its,ite, jts,jte, kts,kte,       &
                              amp, omega, kw, t_rk, tau, rdx, rdy, ht, &
                              cf1, cf2, cf3, move, ramp, test_wave       ) !PZ added for movingbottom


  IMPLICIT NONE  ! religion first

!  stuff passed in

  INTEGER,                  INTENT(IN   ) :: ids,ide, jds,jde, kds,kde
  INTEGER,                  INTENT(IN   ) :: ims,ime, jms,jme, kms,kme
  INTEGER,                  INTENT(IN   ) :: its,ite, jts,jte, kts,kte
  INTEGER,                  INTENT(IN   ) :: number_of_small_timesteps
  INTEGER,                  INTENT(IN   ) :: rk_step, rk_order
  REAL,                     INTENT(IN   ) :: dts

  REAL,   DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN   ) :: u_1, &
                                                                 v_1, &
                                                                 w_1, &
                                                                 t_1, &
                                                                 ww1, &
                                                                 ph_1

  REAL,   DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: u_2, &
                                                                 v_2, &
                                                                 w_2, &
                                                                 t_2, &
                                                                 ww,  &
                                                                 ph_2

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(IN   ) :: u_save,   &
                                                              v_save,   &
                                                              w_save,   &
                                                              t_save,   &
                                                              ph_save,  &
                                                             h_diabatic

  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: muus, muvs
  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: mu_2, mu_1
  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: mut, muts, &
                                                        muu, muv, mu_save
  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(IN   ) :: msfux, msfuy, &
                                                        msfvx, msfvy, &
                                                        msftx, msfty


! local stuff

  INTEGER         :: i,j,k
  INTEGER :: i_start, i_end, j_start, j_end, i_endu, j_endv
 
! PZ added for moving wave
  REAL, INTENT(IN   ) :: amp, omega, kw, t_rk, tau, rdx, rdy, cf1, cf2, cf3 
  REAL:: x_dist
  LOGICAL, INTENT(IN   ) :: move, ramp
  INTEGER, INTENT(IN   ) :: test_wave
  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(IN   ) :: ht

!<DESCRIPTION>
!
!  small_step_finish reconstructs the full uncoupled prognostic variables
!  from the coupled perturbation variables used in the small timesteps.
!
!</DESCRIPTION>

  i_start = its
  i_end   = min(ite,ide-1)
  j_start = jts
  j_end   = min(jte,jde-1)

  i_endu = ite
  j_endv = jte

!    addition of time level t back into variables

  DO j = j_start, j_endv
  DO k = kds, kde-1
  DO i = i_start, i_end
    ! v coupled with mx
    v_2(i,k,j) = (msfvx(i,j)*v_2(i,k,j) + v_save(i,k,j)*muv(i,j))/muvs(i,j)
  ENDDO
  ENDDO
  ENDDO

  DO j = j_start, j_end
  DO k = kds, kde-1
  DO i = i_start, i_endu
    ! u coupled with my
    u_2(i,k,j) = (msfuy(i,j)*u_2(i,k,j) + u_save(i,k,j)*muu(i,j))/muus(i,j)
  ENDDO
  ENDDO
  ENDDO

  DO j = j_start, j_end
  DO k = kds, kde
  DO i = i_start, i_end
    ! w coupled with my

    w_2(i,k,j) = (msfty(i,j)*w_2(i,k,j) + w_save(i,k,j)*mut(i,j))/muts(i,j)
    ph_2(i,k,j) = ph_2(i,k,j) + ph_save(i,k,j)


    ww(i,k,j) = ww(i,k,j) + ww1(i,k,j)
  ENDDO
  ENDDO
  ENDDO
  
!  IF (move) THEN
!    DO j = j_start, j_end
!    DO i = i_start, i_end
!      x_dist = (i-0.5)/rdx
!      SELECT CASE (test_wave)
!      CASE (2)
!        IF (ramp) THEN
!          w_2(i,1,j) =  -( amp*omega*cos(kw*x_dist -omega*t_rk)*tanh(2.0*t_rk/tau) - 2*amp/tau*sin(kw*x_dist - omega*t_rk)*(1 - tanh(2*t_rk/tau)**2)) + msfty(i,j)*.5*rdy*(                        &
!                      (ht(i,j+1)-ht(i,j  ))             &
!        *(cf1*v_2(i,1,j+1)+cf2*v_2(i,2,j+1)+cf3*v_2(i,3,j+1))    &
!                      +(ht(i,j  )-ht(i,j-1))             &
!        *(cf1*v_2(i,1,j  )+cf2*v_2(i,2,j  )+cf3*v_2(i,3,j  ))  ) &
!              +msftx(i,j)*.5*rdx*(                        &
!                       (ht(i+1,j)-ht(i,j  ))             &
!        *(cf1*u_2(i+1,1,j)+cf2*u_2(i+1,2,j)+cf3*u_2(i+1,3,j))    &
!                      +(ht(i,j  )-ht(i-1,j))             &
!        *(cf1*u_2(i  ,1,j)+cf2*u_2(i  ,2,j)+cf3*u_2(i  ,3,j)) ) 
!          ph_2(i,1,j) =  amp*sin(kw*x_dist -omega*t_rk)*tanh(2.0*t_rk/tau)*g - 0
!      ENDIF
!    END SELECT
!   END DO
!   END DO
!     
!  ENDIF

#ifdef REVERT
  DO j = j_start, j_end
  DO k = kds, kde-1
  DO i = i_start, i_end
    t_2(i,k,j) = (t_2(i,k,j) + t_save(i,k,j)*mut(i,j))/muts(i,j)
  ENDDO
  ENDDO
  ENDDO
#else
  IF ( rk_step < rk_order ) THEN
    DO j = j_start, j_end
    DO k = kds, kde-1
    DO i = i_start, i_end
      t_2(i,k,j) = (t_2(i,k,j) + t_save(i,k,j)*mut(i,j))/muts(i,j)
    ENDDO
    ENDDO
    ENDDO
  ELSE

    DO j = j_start, j_end
    DO k = kds, kde-1
    DO i = i_start, i_end
      t_2(i,k,j) = (t_2(i,k,j) - dts*number_of_small_timesteps*mut(i,j)*h_diabatic(i,k,j) &
                    + t_save(i,k,j)*mut(i,j))/muts(i,j)
    ENDDO
    ENDDO
    ENDDO
  ENDIF
#endif

  DO j = j_start, j_end
  DO i = i_start, i_end
    mu_2(i,j) = mu_2(i,j) + mu_save(i,j)
  ENDDO
  ENDDO

END SUBROUTINE small_step_finish

!-----------------------------------------------------------------------

SUBROUTINE calc_p_rho( al, p, ph,                    &
                       alt, t_2, t_1, c2a, pm1,      &
                       mu, muts, znu, t0,            &
                       rdnw, dnw, smdiv,             &
                       non_hydrostatic, step,        &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its,ite, jts,jte, kts,kte, &
                       write_acous    )

  IMPLICIT NONE  ! religion first

! declarations for the stuff coming in

  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  INTEGER,      INTENT(IN   )    :: step

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(  OUT) :: al,   &
                                                               p

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(IN   ) :: alt,   &
                                                              t_2,   &
                                                              t_1,   &
                                                              c2a

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT) :: ph, pm1

  REAL, DIMENSION(ims:ime, jms:jme)         , INTENT(IN   ) :: mu,   &
                                                               muts

  REAL, DIMENSION(kms:kme)         , INTENT(IN   ) :: dnw,  &
                                                      rdnw, &
                                                      znu

  REAL,                                       INTENT(IN   ) :: t0, smdiv

  LOGICAL, INTENT(IN   )  :: non_hydrostatic

! local variables

  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  REAL    :: ptmp

  CHARACTER (LEN=26) :: profile_format !PZ added
  LOGICAL, INTENT(IN) :: write_acous

!<DESCRIPTION>
!
!  For the nonhydrostatic option,
!  calc_p_rho computes the perturbation inverse density and 
!  perturbation pressure from the hydrostatic relation and the 
!  linearized equation of state, respectively.
!
!  For the hydrostatic option,
!  calc_p_rho computes the perturbation pressure, perturbation density,
!  and perturbation geopotential
!  from the vertical coordinate definition, linearized equation of state 
!  and the hydrostatic relation, respectively.
!
!  forward weighting of the pressure (divergence damping) is also
!  computed here.
!
!</DESCRIPTION>

   i_start = its
   i_end   = min(ite,ide-1)
   j_start = jts
   j_end   = min(jte,jde-1)
   k_start = kts
   k_end = min(kte,kde-1)

   IF (non_hydrostatic) THEN
     DO j=j_start, j_end
     DO k=k_start, k_end
     DO i=i_start, i_end

!  al computation is all dry, so ok with moisture
 
      al(i,k,j)=-1./muts(i,j)*(alt(i,k,j)*mu(i,j)  &
             +rdnw(k)*(ph(i,k+1,j)-ph(i,k,j)))

!  this is temporally linearized p, no moisture correction needed

      p(i,k,j)=c2a(i,k,j)*(alt(i,k,j)*(t_2(i,k,j)-mu(i,j)*t_1(i,k,j))  &
                       /(muts(i,j)*(t0+t_1(i,k,j)))-al (i,k,j))
      if (i<=10 .and. j==1 .and. k==1) write(*,*) 'after calc_p_rho, al=', al (i,k,j),'dphdz', ph(i,k+1,j)-ph(i,k,j), 'mu',mu(i,j),'muts', muts(i,j),'ph k2', ph(i,k+1,j),'ph k1', ph(i,k,j),'p',p(i,k,j)  
      !if (k==k_end) p(i,k,j) = 0 !PZ modified
     ENDDO
     ENDDO
     ENDDO

   ELSE  ! hydrostatic calculation

       DO j=j_start, j_end
       DO k=k_start, k_end
       DO i=i_start, i_end
         p(i,k,j)=mu(i,j)*znu(k)
         al(i,k,j)=alt(i,k,j)*(t_2(i,k,j)-mu(i,j)*t_1(i,k,j))            &
                      /(muts(i,j)*(t0+t_1(i,k,j)))-p(i,k,j)/c2a(i,k,j)
         ph(i,k+1,j)=ph(i,k,j)-dnw(k)*(muts(i,j)*al (i,k,j)              &
                          +mu(i,j)*alt(i,k,j))
       ENDDO
       ENDDO
       ENDDO

   END IF


!  divergence damping setup
 
     IF (step == 0) then   ! we're initializing small timesteps
       DO j=j_start, j_end
       DO k=k_start, k_end
       DO i=i_start, i_end
         pm1(i,k,j)=p(i,k,j)
       ENDDO
       ENDDO
       ENDDO 
     ELSE                     ! we're in the small timesteps 
       DO j=j_start, j_end    ! and adding div damping component
       DO k=k_start, k_end
       DO i=i_start, i_end
         ptmp = p(i,k,j)
         p(i,k,j) = p(i,k,j) + smdiv*(p(i,k,j)-pm1(i,k,j))
         pm1(i,k,j) = ptmp
       ENDDO
       ENDDO
       ENDDO
     END IF

    if (write_acous) then
    OPEN(UNIT=100, FILE=TRIM('acoustic_p_pp001'), STATUS='unknown',POSITION='append', FORM='formatted')

    profile_format = '(000(f14.10,1x))'
    WRITE(profile_format(2:4),'(I3.3)') k_end
    WRITE(UNIT=100,FMT=profile_format) p(1,1:k_end, 1) 
    CLOSE(UNIT=100)

    OPEN(UNIT=100, FILE=TRIM('acoustic_p_pp010'), STATUS='unknown',POSITION='append',FORM='formatted')
    WRITE(UNIT=100,FMT=profile_format) p(10,1:k_end, 1) 
    CLOSE(UNIT=100)
    endif
END SUBROUTINE calc_p_rho

!----------------------------------------------------------------------

SUBROUTINE calc_coef_w( a,alpha,gamma,              &
                        mut, cqw,                   &
                        rdn, rdnw, c2a,             &
                        dts, g, epssm, top_lid,     &
                        ids,ide, jds,jde, kds,kde,  & ! domain dims
                        ims,ime, jms,jme, kms,kme,  & ! memory dims
                        its,ite, jts,jte, kts,kte,set_w )  ! tile   dims
                                                   
      IMPLICIT NONE  ! religion first

!  passed in through the call

  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  LOGICAL,      INTENT(IN   )    :: top_lid

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN   ) :: c2a,  &
                                                               cqw

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: alpha, &
                                                               gamma, &
                                                               a

  REAL, DIMENSION(ims:ime, jms:jme),          INTENT(IN   ) :: mut

  REAL, DIMENSION(kms:kme),                   INTENT(IN   ) :: rdn,   &
                                                               rdnw

  REAL,                                       INTENT(IN   ) :: epssm, &
                                                               dts,   &
                                                               g

!  Local stack data.

  REAL, DIMENSION(ims:ime)                         :: cof
  REAL  :: b, c

  INTEGER :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end
  INTEGER :: ij, ijp, ijm, lid_flag
 
  logical,  INTENT(IN   ) :: set_w !PZ added

!<DESCRIPTION>
!
!  calc_coef_w calculates the coefficients needed for the 
!  implicit solution of the vertical momentum and geopotential equations.
!  This requires solution of a tri-diagonal equation.
!
!</DESCRIPTION>

      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
      k_start = kts
      k_end   = kte-1

      lid_flag=1
      IF(top_lid)lid_flag=0

     outer_j_loop:  DO j = j_start, j_end

        DO i = i_start, i_end
          cof(i)  = (.5*dts*g*(1.+epssm)/mut(i,j))**2
          a(i, 2 ,j) = 0
          a(i,kde,j) =-2.*cof(i)*rdnw(kde-1)**2*c2a(i,kde-1,j)*lid_flag
          gamma(i,1 ,j) = 0.
        ENDDO

        DO k=3,kde-1
        DO i=i_start, i_end
          a(i,k,j) =   -cqw(i,k,j)*cof(i)*rdn(k)* rdnw(k-1)*c2a(i,k-1,j)    
          
        ENDDO
        ENDDO


        DO k=2,kde-1
        DO i=i_start, i_end
          b = 1.+cqw(i,k,j)*cof(i)*rdn(k)*(rdnw(k  )*c2a(i,k,j  )  &
                                       +rdnw(k-1)*c2a(i,k-1,j))
          c =   -cqw(i,k,j)*cof(i)*rdn(k)*rdnw(k  )*c2a(i,k,j  )
          alpha(i,k,j) = 1./(b-a(i,k,j)*gamma(i,k-1,j))
          gamma(i,k,j) = c*alpha(i,k,j)
        ENDDO
        ENDDO

        DO i=i_start, i_end
          b = 1.+2.*cof(i)*rdnw(kde-1)**2*c2a(i,kde-1,j)
          c = 0.
          alpha(i,kde,j) = 1./(b-a(i,kde,j)*gamma(i,kde-1,j))
          gamma(i,kde,j) = c*alpha(i,kde,j)
        ENDDO

        !PZ added
        IF (set_w) then
        Do i= i_start, i_end
            a(i,2,j) = -cof(i)*rdn(2)*rdnw(1)*c2a(i,1,j)
        ENDDO
        ENDIF
      ENDDO outer_j_loop

END SUBROUTINE calc_coef_w

!----------------------------------------------------------------------

SUBROUTINE advance_uv ( u, ru_tend, v, rv_tend,        &
                        p, pb,                         &
                        ph, php, alt, al, mu,          &
                        muu, cqu, muv, cqv, mudf,      &
                        msfux, msfuy, msfvx,           &
                        msfvx_inv, msfvy,              &
                        rdx, rdy, dts,                 &
                        cf1, cf2, cf3, fnm, fnp,       &
                        ht, amp, kw, H, omega, t_old, move, write_acous, iteration, test_wave, u_save, phb, muus, & !PZ for moving bottom
                        emdiv,                         &
                        rdnw, config_flags, spec_zone, &
                        non_hydrostatic, top_lid,      &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



      IMPLICIT NONE  ! religion first

! stuff coming in

      TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      LOGICAL, INTENT(IN   ) :: non_hydrostatic, top_lid

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      INTEGER,      INTENT(IN   )    :: spec_zone

      REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),  &
            INTENT(INOUT) ::                          &
                                                  u,  &
                                                  v,  &
                                                  u_save, phb !PZ added

      REAL, DIMENSION( ims:ime , kms:kme, jms:jme ), &
            INTENT(IN   ) ::                          &
                                             ru_tend, &
                                             rv_tend, &
                                             ph,      &
                                             php,     &
                                             p,       &
                                             pb,      &
                                             alt,     &
                                             al,      &
                                             cqu,     &
                                             cqv


      REAL, DIMENSION( ims:ime , jms:jme ),    INTENT(IN   ) :: muu,  &
                                                                muv,  &
                                                                mu,   &
                                                                mudf, &
                                                                muus, ht !PZ added for moving bottom


      REAL, DIMENSION( kms:kme ),              INTENT(IN   ) :: fnm,    &
                                                                fnp ,   &
                                                                rdnw

      REAL, DIMENSION( ims:ime , jms:jme ),    INTENT(IN   ) :: msfux,  &
                                                                msfuy,  &
                                                                msfvx,  &
                                                                msfvy,  &
                                                                msfvx_inv

      REAL,                                    INTENT(IN   ) :: rdx,    &
                                                                rdy,    &
                                                                dts,    &
                                                                cf1,    &
                                                                cf2,    &
                                                                cf3,    &
                                                                emdiv
    

!  Local 3d array from the stack (note tile size)

      REAL, DIMENSION (its:ite, kts:kte) :: dpn, dpxy
      REAL, DIMENSION (its:ite) :: mudf_xy
      REAL                      :: dx, dy

      INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
      INTEGER :: i_endu, j_endv, k_endw
      INTEGER :: i_start_up, i_end_up, j_start_up, j_end_up
      INTEGER :: i_start_vp, i_end_vp, j_start_vp, j_end_vp

      INTEGER :: i_start_u_tend, i_end_u_tend, j_start_v_tend, j_end_v_tend

! PZ added for moving bottom
      INTEGER, INTENT(IN   ) :: iteration, test_wave 
      REAL, INTENT(IN ) :: amp, omega, kw, H, t_old
      LOGICAL, INTENT(IN ) :: move, write_acous
      REAL:: x_dist, z
      CHARACTER (LEN=26) :: profile_format

!<DESCRIPTION>
!
!  advance_uv advances the explicit perturbation horizontal momentum 
!  equations (u,v) by adding in the large-timestep tendency along with 
!  the small timestep pressure gradient tendency.
!
!</DESCRIPTION>

!  now, the real work.
!  set the loop bounds taking into account boundary conditions.

    IF( config_flags%nested .or. config_flags%specified ) THEN
      i_start = max( its,ids+spec_zone )
      i_end   = min( ite,ide-spec_zone-1 )
      j_start = max( jts,jds+spec_zone )
      j_end   = min( jte,jde-spec_zone-1 )
      k_start = kts
      k_end   = min( kte, kde-1 )

      i_endu = min( ite,ide-spec_zone )
      j_endv = min( jte,jde-spec_zone )
      k_endw = k_end

      IF( config_flags%periodic_x) THEN
        i_start = its
        i_end   = min(ite,ide-1)
        i_endu = ite
      ENDIF
    ELSE
      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
      k_start = kts
      k_end   = kte-1

      i_endu = ite
      j_endv = jte
      k_endw = k_end
    ENDIF

      i_start_up = i_start
      i_end_up   = i_endu
      j_start_up = j_start
      j_end_up   = j_end

      i_start_vp = i_start
      i_end_vp   = i_end
      j_start_vp = j_start
      j_end_vp   = j_endv

      IF ( (config_flags%open_xs   .or.     &
            config_flags%symmetric_xs   )   &
            .and. (its == ids) )            &
                 i_start_up = i_start_up + 1

      IF ( (config_flags%open_xe    .or.  &
            config_flags%symmetric_xe   ) &
             .and. (ite == ide) )         &
                 i_end_up   = i_end_up - 1

      IF ( (config_flags%open_ys    .or.   &
            config_flags%symmetric_ys  .or.   &
            config_flags%polar   )  &
                     .and. (jts == jds) )  &
                 j_start_vp = j_start_vp + 1

      IF ( (config_flags%open_ye     .or. &
            config_flags%symmetric_ye  .or.   &
            config_flags%polar   )  &
            .and. (jte == jde) )          &
                 j_end_vp   = j_end_vp - 1

      i_start_u_tend = i_start
      i_end_u_tend   = i_endu
      j_start_v_tend = j_start
      j_end_v_tend   = j_endv

      IF ( config_flags%symmetric_xs .and. (its == ids) ) &
                     i_start_u_tend = i_start_u_tend+1
      IF ( config_flags%symmetric_xe .and. (ite == ide) ) &
                     i_end_u_tend = i_end_u_tend-1
      IF ( config_flags%symmetric_ys .and. (jts == jds) ) &
                     j_start_v_tend = j_start_v_tend+1
      IF ( config_flags%symmetric_ye .and. (jte == jde) ) &
                     j_end_v_tend = j_end_v_tend-1

   dx = 1./rdx
   dy = 1./rdy

!  start real calculations.
!  first, u

  u_outer_j_loop: DO j = j_start, j_end

   DO k = k_start, k_end
   DO i = i_start_u_tend, i_end_u_tend
     u(i,k,j) = u(i,k,j) + dts*ru_tend(i,k,j)
   ENDDO
   ENDDO
   !PZ added
   write(*,*) 'in advance_uv, dts*ru_tend(k1)=', dts*ru_tend(5,1,1),'ru_tend(k2)=',dts*ru_tend(5,2,1)

   DO i = i_start_up, i_end_up
     mudf_xy(i)= -emdiv*dx*(mudf(i,j)-mudf(i-1,j))/msfuy(i,j)
   ENDDO

   DO k = k_start, k_end
   DO i = i_start_up, i_end_up

!     Comments on map scale factors:   
!     x pressure gradient: ADT eqn 44, penultimate term on RHS
!        = -(mx/my)*(mu/rho)*partial dp/dx
!     [i.e., first rho->mu; 2nd still rho; alpha=1/rho]
!     Klemp et al. splits into 2 terms:
!        mu alpha partial dp/dx + partial dp/dnu * partial dphi/dx 
!     then into 4 terms:
!        mu alpha partial dp'/dx + nu mu alpha' partial dmubar/dx +
!      + mu partial dphi/dx + partial dphi'/dx * (partial dp'/dnu - mu')
!
!     first 3 terms:
!     ph, alt, p, al, pb not coupled
!     since we want tendency to fit ADT eqn 44 (coupled) we need to
!     multiply by (mx/my):
!
     dpxy(i,k)= (msfux(i,j)/msfuy(i,j))*.5*rdx*muu(i,j)*(         &
       ((ph (i,k+1,j)-ph (i-1,k+1,j))+(ph (i,k,j)-ph (i-1,k,j)))  &
      +(alt(i,k  ,j)+alt(i-1,k  ,j))*(p  (i,k,j)-p  (i-1,k,j))  &
      +(al (i,k  ,j)+al (i-1,k  ,j))*(pb (i,k,j)-pb (i-1,k,j)) ) 

   ENDDO
   ENDDO
   !PZ added
   if (j==1) write(*,*) 'dpxy(k1)=', dpxy(5,1), 'dpxy(k2)=', dpxy(5,2)

   IF (non_hydrostatic) THEN

     DO i = i_start_up, i_end_up
       x_dist = (i-0.5)*dx
!       IF (move) THEN  ! PZ added for moving bottom
!           SELECT CASE (test_wave)
!             CASE (1)
!                 dpn(i,1) = 0.5*(-1/alt(i,1,j)*amp*g*sin(kw*x_dist)*cos(omega*(t_old + dts*iteration)) - 1/alt(i,1,j)*g*(phb(i,1,j)+ph(i,1,j))/g &
!                  -1/alt(i-1,1,j)*amp*g*sin(kw*(x_dist-dx))*cos(omega*(t_old + dts*iteration)) - 1/alt(i-1,1,j)*g*(phb(i-1,1,j)+ph(i-1,1,j))/g)
!!               
!                 if (i==1 .and. j==1) write(*, *) 'ana p=', dpn(i,1), 'ori p', .5*( cf1*(p(i,1,j)+p(i-1,1,j)) +cf2*(p(i,2,j)+p(i-1,2,j))+cf3*(p(i,3,j)+p(i-1,3,j)) )

!             CASE(2)
!               x_dist = (i-1)*dx
!               dpn(i,1) =  .5*( cf1*(p(i,1,j)+p(i-1,1,j))  &
!                      +cf2*(p(i,2,j)+p(i-1,2,j))  &
!                      +cf3*(p(i,3,j)+p(i-1,3,j)) ) !(-1/alt(i,1,j)*amp*omega**2/kw*cosh(kw*(ht(i,j)-H))/sinh(kw*H) * sin(kw*x_dist-omega*(t_old + dts*iteration)) + (-1/alt(i-1,1,j)*amp*omega**2/kw*cosh(kw*(ht(i-1,j)/g-H))/sinh(kw*H) * sin(kw*(x_dist+dx)-omega*(t_old + dts*iteration)) ) )*0.5
!               write(*,*) 'i=', i, 'dpn=', dpn(i,1),'k=',k, 'dpn ori=',  .5*(cf1*(p(i,1,j)+p(i-1,1,j)) +cf2*(p(i,2,j)+p(i-1,2,j)) + cf3*(p(i,3,j)+p(i-1,3,j))), 'p(1)=', 0.5*(p(i,1,j)+p(i-1,1,j)), 'p(2)=', 0.5*(p(i,2,j)+p(i-1,2,j)),'p(3)=', 0.5*(p(i,3,j)+p(i-1,3,j)) 
!
!             CASE (6)
!               x_dist = (i-0.5)*dx
!               dpn(i,1) = 0.5*(-1/alt(i,1,j)*amp*g*cosh(kw*((phb(i,1,j)+ph(i,1,j))/g-H))/sinh(kw*H) * cos(kw*x_dist)*sin(omega*(t_old + dts*iteration)) - 1/alt(i,1,j)*g*(phb(i,1,j)+ph(i,1,j))/g  &
!                        + (-1/alt(i-1,1,j)*amp*g*cosh(kw*((phb(i-1,1,j)+ph(i-1,1,j))/g-H))/sinh(kw*H)*cos(kw*(x_dist-dx) )*sin(omega*(t_old + dts*iteration)) - 1/alt(i-1,1,j)*g*(phb(i-1,1,j)+ph(i-1,1,j))/g) )
!           END SELECT
!       ELSE
         dpn(i,1) = .5*( cf1*(p(i,1,j)+p(i-1,1,j))  &
                      +cf2*(p(i,2,j)+p(i-1,2,j))  &
                      +cf3*(p(i,3,j)+p(i-1,3,j)) )
!       ENDIF
       dpn(i,kde) = 0.
     ENDDO
     IF (top_lid) THEN
       DO i = i_start_up, i_end_up
         dpn(i,kde) =.5*( cf1*(p(i-1,kde-1,j)+p(i,kde-1,j))   &
                         +cf2*(p(i-1,kde-2,j)+p(i,kde-2,j))   &
                         +cf3*(p(i-1,kde-3,j)+p(i,kde-3,j))  )
       ENDDO
     ENDIF

     DO k = k_start+1, k_end
       DO i = i_start_up, i_end_up
         dpn(i,k) = .5*( fnm(k)*(p(i,k  ,j)+p(i-1,k  ,j))   &
                        +fnp(k)*(p(i,k-1,j)+p(i-1,k-1,j)) )
       ENDDO
     ENDDO

!    Comments on map scale factors:
!    4th term:
!    php, dpn, mu not coupled
!    since we want tendency to fit ADT eqn 44 (coupled) we need to
!    multiply by (mx/my):

     DO k = k_start, k_end
       DO i = i_start_up, i_end_up
         dpxy(i,k)=dpxy(i,k) + (msfux(i,j)/msfuy(i,j))*rdx*(php(i,k,j)-php(i-1,k,j))*        &
           (rdnw(k)*(dpn(i,k+1)-dpn(i,k))-.5*(mu(i-1,j)+mu(i,j)))
       ENDDO
     ENDDO

 
   END IF


   DO k = k_start, k_end
     DO i = i_start_up, i_end_up
       u(i,k,j)=u(i,k,j)-dts*cqu(i,k,j)*dpxy(i,k)+mudf_xy(i) 
       if (i<4 .and. k<=3 .and. j==1) write(*,*) 'i=',i, 'k=',k, 'upp=', u(i,k,j)
     ENDDO
   ENDDO

   ENDDO u_outer_j_loop

    !PZ added
   if (write_acous) then
   OPEN(UNIT=100, FILE=TRIM('acoustic_u001'), STATUS='unknown',POSITION='append', FORM='formatted')
    profile_format = '(000(f14.10,1x))'
    WRITE(profile_format(2:4),'(I3.3)') k_end
    WRITE(UNIT=100,FMT=profile_format) (u(1,1:k_end, 1) + muu(1,1)*u_save(1,1:k_end,1))/muus(1,1)
    CLOSE(UNIT=100)

    OPEN(UNIT=100, FILE=TRIM('acoustic_u002'),STATUS='unknown',POSITION='append',FORM='formatted')
    WRITE(UNIT=100,FMT=profile_format) (u(2,1:k_end, 1) +muu(2,1)*u_save(2,1:k_end,1))/muus(2,1)
    CLOSE(UNIT=100)
   end if
 

! now v

  v_outer_j_loop: DO j = j_start_v_tend, j_end_v_tend


   DO k = k_start, k_end
   DO i = i_start, i_end
     v(i,k,j) = v(i,k,j) + dts*rv_tend(i,k,j)
   ENDDO
   ENDDO

   DO i = i_start, i_end
     mudf_xy(i)= -emdiv*dy*(mudf(i,j)-mudf(i,j-1))*msfvx_inv(i,j)
   ENDDO

     IF (     ( j >= j_start_vp)  &
       .and.( j <= j_end_vp  ) )  THEN

         DO k = k_start, k_end
	 DO i = i_start, i_end

!      Comments on map scale factors:
!      y pressure gradient: ADT eqn 45, penultimate term on RHS
!         = -(my/mx)*(mu/rho)*partial dp/dy
!      [i.e., first rho->mu; 2nd still rho; alpha=1/rho]
!      Klemp et al. splits into 2 terms:
!         mu alpha partial dp/dy + partial dp/dnu * partial dphi/dy 
!      then into 4 terms:
!         mu alpha partial dp'/dy + nu mu alpha' partial dmubar/dy +
!       + mu partial dphi/dy + partial dphi'/dy * (partial dp'/dnu - mu')
!
!      first 3 terms:
!      ph, alt, p, al, pb not coupled
!      since we want tendency to fit ADT eqn 45 (coupled) we need to
!      multiply by (my/mx):
!      mudf_xy is NOT a map scale factor coupling
!      it is some sort of divergence damping

       dpxy(i,k)= (msfvy(i,j)/msfvx(i,j))*.5*rdy*muv(i,j)*(       &
        ((ph(i,k+1,j)-ph(i,k+1,j-1))+(ph (i,k,j)-ph (i,k,j-1)))   &
        +(alt(i,k  ,j)+alt(i,k  ,j-1))*(p  (i,k,j)-p  (i,k,j-1))  &
        +(al (i,k  ,j)+al (i,k  ,j-1))*(pb (i,k,j)-pb (i,k,j-1)) ) 
	    ENDDO
	    ENDDO


     IF (non_hydrostatic) THEN

       DO i = i_start, i_end     
!         if (move) then  ! PZ added for moving bottom
!           x_dist = (i-0.5)*dx
!           SELECT CASE (test_wave)
!             CASE (1)
!                 dpn(i,1) = 0.5*(-1/alt(i,1,j-1)*amp*g*sin(kw*x_dist)*cos(omega*(t_old + dts*iteration)) - 1/alt(i,1,j-1)*g*(phb(i,1,j-1)+ph(i,1,j-1))/g &
!                  -1/alt(i,1,j)*amp*g*sin(kw*x_dist)*cos(omega*(t_old +dts*iteration)) - 1/alt(i,1,j)*g*(phb(i,1,j)+ph(i,1,j))/g)
!                 dpn(i,1) = .5*( cf1*(p(i,1,j)+p(i-1,1,j))  &
!                      +cf2*(p(i,2,j)+p(i-1,2,j))  &
!                      +cf3*(p(i,3,j)+p(i-1,3,j)) )
!             CASE(2)
!               x_dist = (i-0.5)*dx
!               dpn(i,1) = .5*( cf1*(p(i,1,j)+p(i,1,j-1))  &
!                        +cf2*(p(i,2,j)+p(i,2,j-1))  &
!                        +cf3*(p(i,3,j)+p(i,3,j-1)) ) !-1/alt(i,1,j)*amp*omega**2/kw*cosh(kw*(php(i,1,j)/g -H))/sinh(kw*H) * sin(kw*x_dist-omega*(t_old + dts*iteration)) 
!             CASE(6)
!               x_dist = (i-0.5)*dx
!               dpn(i,1) = 0.5*(-1/alt(i,1,j)*amp*g*cosh(kw*((phb(i,1,j)+ph(i,1,j))/g-H))/sinh(kw*H) * cos(kw*x_dist)*sin(omega*(t_old + dts*iteration)) - 1/alt(i,1,j)*g*(phb(i,1,j)+ph(i,1,j))/g + (-1/alt(i,1,j-1)*amp*g*cosh(kw*((phb(i,1,j-1)+ph(i,1,j-1))/g-H))/sinh(kw*H) * cos(kw*x_dist)*sin(omega*(t_old + dts*iteration)) - 1/alt(i,1,j-1)*g*(phb(i,1,j-1)+ph(i,1,j-1))/g) )
!           END SELECT
!         else
           dpn(i,1) = .5*( cf1*(p(i,1,j)+p(i,1,j-1))  &
                        +cf2*(p(i,2,j)+p(i,2,j-1))  &
                        +cf3*(p(i,3,j)+p(i,3,j-1)) )
!         endif
         dpn(i,kde) = 0.
       ENDDO
     IF (top_lid) THEN
       DO i = i_start, i_end
         dpn(i,kde) =.5*( cf1*(p(i,kde-1,j-1)+p(i,kde-1,j))   &
                         +cf2*(p(i,kde-2,j-1)+p(i,kde-2,j))   &
                         +cf3*(p(i,kde-3,j-1)+p(i,kde-3,j))  )
       ENDDO
     ENDIF

       DO k = k_start+1, k_end
         DO i = i_start, i_end
           dpn(i,k) = .5*( fnm(k)*(p(i,k  ,j)+p(i,k  ,j-1))  &
                          +fnp(k)*(p(i,k-1,j)+p(i,k-1,j-1)) )
         ENDDO
       ENDDO

!      Comments on map scale factors:
!      4th term:
!      php, dpn, mu not coupled
!      since we want tendency to fit ADT eqn 45 (coupled) we need to
!      multiply by (my/mx):

       DO k = k_start, k_end
         DO i = i_start, i_end
           dpxy(i,k)=dpxy(i,k) + (msfvy(i,j)/msfvx(i,j))*rdy*(php(i,k,j)-php(i,k,j-1))*    &
             (rdnw(k)*(dpn(i,k+1)-dpn(i,k)) -.5*(mu(i,j-1)+mu(i,j)))
         ENDDO
       ENDDO


     END IF


     DO k = k_start, k_end
       DO i = i_start, i_end
         v(i,k,j)=v(i,k,j)-dts*cqv(i,k,j)*dpxy(i,k)+mudf_xy(i) 
       ENDDO
     ENDDO
   END IF

  ENDDO  v_outer_j_loop

  ! The check for j_start_vp and j_end_vp makes sure that the edges in v
  ! are not updated.  Let's add a polar boundary condition check here for
  ! safety to ensure that v at the poles never gets to be non-zero in the
  ! small time steps.
  IF (config_flags%polar) THEN
     IF (jts == jds) THEN
        DO k = k_start, k_end
        DO i = i_start, i_end
           v(i,k,jds) = 0.
        ENDDO
        ENDDO
     END IF
     IF (jte == jde) THEN
        DO k = k_start, k_end
        DO i = i_start, i_end
           v(i,k,jde) = 0.
        ENDDO
        ENDDO
     END IF
  END IF

  !PZ added
  write(*,*) 'END of advance_uv, u=', u(1,1,1), u(2,1,1), u(3,1,1), u(1,2,1)


END SUBROUTINE advance_uv

!---------------------------------------------------------------------

SUBROUTINE advance_mu_t( ww, ww_1, u, u_1, v, v_1,            &
                         mu, mut, muave, muts, muu, muv,      &
                         mudf, uam, vam, wwam, t, t_1,        &
                         t_ave, ft, mu_tend,                  &
                         rdx, rdy, dts, epssm,                &
                         dnw, fnm, fnp, rdnw,                 &
                         msfux, msfuy, msfvx, msfvx_inv,      &
                         msfvy, msftx, msfty,                 &
                         step, config_flags,                  &
                         ids, ide, jds, jde, kds, kde,        &
                         ims, ime, jms, jme, kms, kme,        &
                         its, ite, jts, jte, kts, kte,        &
                         mu_save,write_acous, dmub_dt  )  !PZ added

  IMPLICIT NONE  ! religion first

! stuff coming in

  TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  INTEGER,      INTENT(IN   )    :: step

  REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),   &
            INTENT(IN   ) ::                       &
                                              u,   &
                                              v,   &
                                              u_1, &
                                              v_1, &
                                              t_1, &
                                              ft  

  REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),      &
            INTENT(INOUT) ::                          &
                                              ww,     &
                                              ww_1,   &
                                              t,      &
                                              t_ave,  &
                                              uam,    &
                                              vam,    &
                                              wwam
                                              
  REAL, DIMENSION( ims:ime , jms:jme ),    INTENT(IN   ) :: muu,  &
                                                            muv,  &
                                                            mut,  &
                                                            msfux,&
                                                            msfuy,&
                                                            msfvx,&
                                                            msfvx_inv,&
                                                            msfvy,&
                                                            msftx,&
                                                            msfty,&
                                                            mu_tend, &
                                                            mu_save, dmub_dt !PZ added

  REAL, DIMENSION( ims:ime , jms:jme ),    INTENT(  OUT) :: muave, &
                                                            muts,  &
                                                            mudf

  REAL, DIMENSION( ims:ime , jms:jme ),    INTENT(INOUT) :: mu

  REAL, DIMENSION( kms:kme ),              INTENT(IN   ) :: fnm,    &
                                                            fnp,    &
                                                            dnw,    &
                                                            rdnw


  REAL,                                    INTENT(IN   ) :: rdx,    &
                                                            rdy,    &
                                                            dts,    &
                                                            epssm

!  Local arrays from the stack (note tile size)

  REAL, DIMENSION (its:ite, kts:kte) :: wdtn, dvdxi
  REAL, DIMENSION (its:ite) :: dmdt

  INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
  INTEGER :: i_endu, j_endv
  REAL    :: acc

  CHARACTER (LEN=26) :: profile_format  !PZ added
  LOGICAL, INTENT(IN   ) :: write_acous

!<DESCRIPTION>
!
!  advance_mu_t advances the explicit perturbation theta equation and the mass
!  conservation equation.  In addition, the small timestep omega is updated,
!  and some quantities needed in other places are squirrelled away.
!
!</DESCRIPTION>

!  now, the real work.
!  set the loop bounds taking into account boundary conditions.

  i_start = its
  i_end   = min(ite,ide-1)
  j_start = jts
  j_end   = min(jte,jde-1)
  k_start = kts
  k_end   = kte-1
  IF ( .NOT. config_flags%periodic_x )THEN
     IF ( config_flags%specified .or. config_flags%nested ) then
       i_start = max(its,ids+1)
       i_end   = min(ite,ide-2)
     ENDIF
  ENDIF
  IF ( config_flags%specified .or. config_flags%nested ) then
     j_start = max(jts,jds+1)
     j_end   = min(jte,jde-2)
  ENDIF

  i_endu = ite
  j_endv = jte

!        CALCULATION OF WW (dETA/dt)
   DO j = j_start, j_end

     DO i=i_start, i_end
            dmdt(i) = 0.
     ENDDO
!  NOTE:  mu is not coupled with the map scale factor.
!         ww (omega) IS coupled with the map scale factor.
!         Being coupled with the map scale factor means 
!           multiplication by (1/msft) in this case.

!  Comments on map scale factors
!  ADT eqn 47: 
!  partial drho/dt = -mx*my[partial d/dx(rho u/my) + partial d/dy(rho v/mx)]
!                    -partial d/dz(rho w)
!  with rho -> mu, dividing by my, and with partial d/dnu(rho nu/my [=ww])
!  as the final term (because we're looking for d_nu_/dt)
!
!  begin by integrating with respect to nu from bottom to top
!  BCs are ww=0 at both
!  final term gives 0
!  first term gives Integral([1/my]partial d mu/dt) over total column = dm/dt
!  RHS remaining is Integral(-mx[partial d/dx(mu u/my) + 
!                                partial d/dy(mu v/mx)]) over column
!  lines below find RHS terms at each level then set dmdt = sum over all levels
!
!  [don't divide the below by msfty until find ww, since dmdt is used in
!   the meantime]

     DO k=k_start, k_end
     DO i=i_start, i_end
         dvdxi(i,k) = msftx(i,j)*msfty(i,j)*(                                      &
                     rdy*( (v(i,k,j+1)+muv(i,j+1)*v_1(i,k,j+1)*msfvx_inv(i,j+1))   &
                          -(v(i,k,j  )+muv(i,j  )*v_1(i,k,j  )*msfvx_inv(i,j  )) ) &
                    +rdx*( (u(i+1,k,j)+muu(i+1,j)*u_1(i+1,k,j)/msfuy(i+1,j))       &
                          -(u(i,k,j  )+muu(i  ,j)*u_1(i,k,j  )/msfuy(i  ,j)) ))
        dmdt(i)    = dmdt(i) + dnw(k)*dvdxi(i,k) 
     ENDDO
     ENDDO
     DO i=i_start, i_end
       muave(i,j) = mu(i,j)
       dmdt(i) = dmdt(i) - dmub_dt(i,j) !PZ added
       mu(i,j) = mu(i,j)+dts*(dmdt(i)+mu_tend(i,j) ) 
       mudf(i,j) = (dmdt(i)+mu_tend(i,j)) ! save tendency for div damp filter
       muts(i,j) = mut(i,j)+mu(i,j)
       muave(i,j) =.5*((1.+epssm)*mu(i,j)+(1.-epssm)*muave(i,j))
       
       if (j==1 .and. i<=5) write(*,*) 'mu_tend=', mu_tend(i,j),'dmdt=', dmdt(i), 'mu=', mu(i,j), 'dmub_dt=', dmub_dt(i,j) !PZ added
     ENDDO

     DO k=2,k_end
     DO i=i_start, i_end
       ww(i,k,j)=ww(i,k-1,j)-dnw(k-1)*(dmdt(i) +dvdxi(i,k-1)+mu_tend(i,j))/msfty(i,j) 
       if (i<=3 .and. j==1 .and. k<=5) write(*,*) 'k=', k, 'ww=', ww(i,k,j)
     ENDDO
     END DO

!  NOTE:  ww_1 (large timestep ww) is already coupled with the 
!         map scale factor

     DO k=1,k_end
     DO i=i_start, i_end
       ww(i,k,j)=ww(i,k,j)-ww_1(i,k,j)
     END DO
     END DO

    DO i=i_start, i_end
       !PZ added for moving bottom    
       IF (j==1 .and. step==1) THEN
         write(*,*) 'End of advance_mut'
         write(*,*) 'mu=', mu(1,1), mu(2,1), mu(3,1) 
       ENDIF
    ENDDO


   ENDDO

   write(*,*) 'i_end-i_start+1=', i_end-i_start+1
   !******* PZ added to print out mu at each acoustic step
    if (write_acous) then
    OPEN(UNIT=101, FILE=TRIM('acoustic_mu'), STATUS='unknown',POSITION='append', FORM='formatted')
    profile_format = '(000(f14.10,1x))'
    WRITE(profile_format(2:4),'(I3.3)') (i_end-i_start+1)

    WRITE(UNIT=101,FMT=profile_format) mu(:,1) + mu_save(:,1)
    CLOSE(UNIT=101)
    endif
    !*******

! CALCULATION OF THETA

! NOTE: theta'' is not coupled with the map-scale factor, 
!       while the theta'' tendency is coupled (i.e., mult by 1/msft)

! Comments on map scale factors
! BUT NOTE THAT both are mass coupled
! in flux form equations (Klemp et al.) Theta = mu*theta
!
! scalar eqn: partial d/dt(rho q/my) = -mx[partial d/dx(q rho u/my) + 
!                                          partial d/dy(q rho v/mx)]
!                                      - partial d/dz(q rho w/my)
! with rho -> mu, and with partial d/dnu(q rho nu/my) as the final term
!
! adding previous tendency contribution which was map scale factor coupled
! (had been divided by msfty)
! need to uncouple before updating uncoupled Theta (by adding)

   DO j=j_start, j_end
     DO k=1,k_end
     DO i=i_start, i_end
       t_ave(i,k,j) = t(i,k,j)
       t   (i,k,j) = t(i,k,j) + msfty(i,j)*dts*ft(i,k,j)
     END DO
     END DO
   ENDDO   

   DO j=j_start, j_end

     DO i=i_start, i_end
       wdtn(i,1  )=0.
       wdtn(i,kde)=0.
     ENDDO

     DO k=2,k_end
     DO i=i_start, i_end
        ! for scalar eqn RHS term 3
        wdtn(i,k)= ww(i,k,j)*(fnm(k)*t_1(i,k  ,j)+fnp(k)*t_1(i,k-1,j))
     ENDDO
     ENDDO

! scalar eqn, RHS terms 1, 2 and 3
! multiply by msfty to uncouple result for Theta from map scale factor

     DO k=1,k_end
     DO i=i_start, i_end
       ! multiplication by msfty uncouples result for Theta
       t(i,k,j) = t(i,k,j) - dts*msfty(i,j)*(              &
                          ! multiplication by mx needed for RHS terms 1 & 2
                          msftx(i,j)*(                     &
               .5*rdy*                                     &
              ( v(i,k,j+1)*(t_1(i,k,j+1)+t_1(i,k, j ))     &
               -v(i,k,j  )*(t_1(i,k, j )+t_1(i,k,j-1)) )   &
             + .5*rdx*                                     &
              ( u(i+1,k,j)*(t_1(i+1,k,j)+t_1(i  ,k,j))     &
               -u(i  ,k,j)*(t_1(i  ,k,j)+t_1(i-1,k,j)) ) ) &
             + rdnw(k)*( wdtn(i,k+1)-wdtn(i,k) ) )       
     ENDDO
     ENDDO

   ENDDO
   

END SUBROUTINE advance_mu_t
          


!------------------------------------------------------------

SUBROUTINE advance_w( w, rw_tend, ww, w_save, u, v,  &
                      mu1, mut, muave, muts,      &
                      t_2ave, t_2, t_1,           &
                      ph, ph_1, phb, ph_tend,     &
                      ht, c2a, cqw, alt, alb,     &
                      a, alpha, gamma,            &
                      rdx, rdy, dts, t0, epssm,   &
                      u_save, v_save, omega, iteration, t_old,t_rk, move, write_acous,set_w, add_phi_adv, cfn, cfn1, & !PZ for moving bottom
                      ramp, tau, amp,  k_wave, H, test_wave, muu, muv,dphb_dt, num_sound_steps,phb_rk1, & !PZ for moving bottom
                      dnw, fnm, fnp, rdnw, rdn,   &
                      cf1, cf2, cf3, msftx, msfty,&
                      config_flags, top_lid,      &
                      ids,ide, jds,jde, kds,kde,  & ! domain dims
                      ims,ime, jms,jme, kms,kme,  & ! memory dims
                      its,ite, jts,jte, kts,kte  )  ! tile   dims

! We have used msfty for msft_inv but have not thought through w equation
! pieces properly yet, so we will have to hope that it is okay
! We think we have found a slight error in surface w calculation

  IMPLICIT NONE ! religion first
      
! stuff coming in


  TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  LOGICAL,      INTENT(IN   )    :: top_lid

      REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), &
            INTENT(INOUT) ::                          &
                                             t_2ave,  &
                                             w,       &
                                             ph,      &
                                             phb_rk1  !PZ added

      REAL, DIMENSION(  ims:ime , kms:kme, jms:jme ), &
            INTENT(IN   ) ::                          &
                                             rw_tend, &
                                             ww,      &
                                             w_save,  &
                                             u_save, dphb_dt,v_save,  & !PZ added
                                             u,       &
                                             v,       &
                                             t_2,     &
                                             t_1,     &
                                             ph_1,    &
                                             phb,     &
                                             ph_tend, &
                                             alpha,   &
                                             gamma,   &
                                             a,       &
                                             c2a,     &
                                             cqw,     &
                                             alb,     &
                                             alt     

      REAL, DIMENSION( ims:ime , jms:jme ), &
            INTENT(IN   )  ::               &
                                   mu1,     &
                                   mut,     &
                                   muave,   &
                                   muts,    &
                                   ht,      &
                                   msftx,   &
                                   msfty,   &
                                   muu,muv      !PZ added   


      REAL, DIMENSION( kms:kme ),  INTENT(IN   )  :: fnp,     &
                                                     fnm,     &
                                                     rdnw,    &
                                                     rdn,     &
                                                     dnw

      REAL,   INTENT(IN   )  :: rdx,     &
                                rdy,     &
                                dts,     &
                                cf1,     &
                                cf2,     &
                                cf3,     &
                                t0,      &
                                epssm,   &
                                tau,amp, k_wave, omega, t_old, H, cfn, cfn1 !PZ added for moving bottom 
     REAL:: x_dist, t_rk !PZ added for moving bottom  
      LOGICAL,      INTENT(IN   ):: move, ramp, write_acous, set_w, add_phi_adv !PZ added for moving bottom
      INTEGER, INTENT(IN   ) :: iteration, test_wave, num_sound_steps !PZ added for moving bottom
      CHARACTER (LEN=26) :: profile_format

!  Stack based 3d data, tile size.

      REAL, DIMENSION( its:ite ) :: mut_inv, msft_inv
      REAL, DIMENSION( its:ite, kts:kte ) :: rhs, wdwn
      INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
      REAL, DIMENSION (kts:kte) :: dampwt
      real :: htop,hbot,hdepth,hk
      real    :: pi,dampmag

!<DESCRIPTION>
!
!  advance_w advances the implicit w and geopotential equations.
!
!</DESCRIPTION>

!  set loop limits.
!  Currently set for periodic boundary conditions

      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
      k_start = kts
      k_end   = kte-1
      IF ( .NOT. config_flags%periodic_x )THEN
         IF ( config_flags%specified .or. config_flags%nested ) then
           i_start = max(its,ids+1)
           i_end   = min(ite,ide-2)
         ENDIF
      ENDIF
      IF ( config_flags%specified .or. config_flags%nested ) then
         j_start = max(jts,jds+1)
         j_end   = min(jte,jde-2)
      ENDIF

   pi = 4.*atan(1.)
   dampmag = dts*config_flags%dampcoef
   hdepth=config_flags%zdamp

! calculation of phi and w equations

!   Comments on map scale factors:
!   phi equation is:
!    partial d/dt(rho phi/my) = -mx partial d/dx(phi rho u/my) 
!                               -mx partial d/dy(phi rho v/mx)
!                               - partial d/dz(phi rho w/my) + rho g w/my
!   as with scalar equation, use uncoupled value (here phi) to find the
!   coupled tendency (rho phi/my)
!   here as usual rho -> ~'mu'
!
!   w eqn [divided by my] is:
!     partial d/dt(rho w/my) = -mx partial d/dx(w rho u/my)
!                              -mx partial d/dy(v rho v/mx)
!                              - partial d/dz(w rho w/my)
!                              +rho[(u*u+v*v)/a + 2 u omega cos(lat) -
!                                   (1/rho) partial dp/dz - g + Fz]/my
!   here as usual rho -> ~'mu'
!
!  'u,v,w' sent here must be coupled variables (= rho w/my etc.) by their usage


    DO i=i_start, i_end
      rhs(i,1) = 0.
    ENDDO

    write(*,*) 'in advance_w, t_rk=', t_rk, 'iter=', iteration,'ts=', t_old +dts*iteration  !PZ added
  j_loop_w:  DO j = j_start, j_end
    DO i=i_start, i_end
       mut_inv(i) = 1./mut(i,j)
       msft_inv(i) = 1./msfty(i,j)
    ENDDO

    DO k=1, k_end
    DO i=i_start, i_end
      t_2ave(i,k,j)=.5*((1.+epssm)*t_2(i,k,j)       &
                    +(1.-epssm)*t_2ave(i,k,j))
      t_2ave(i,k,j)=(t_2ave(i,k,j) + muave(i,j)*t0) &
                    /(muts(i,j)*(t0+t_1(i,k,j)))
      wdwn(i,k+1)=.5*(ww(i,k+1,j)+ww(i,k,j))*rdnw(k)    &
           *(ph_1(i,k+1,j)-ph_1(i,k,j)+phb(i,k+1,j)-phb(i,k,j)) 
      rhs(i,k+1) = dts*(ph_tend(i,k+1,j) + .5*g*(1.-epssm)*w(i,k+1,j))

      !PZ added
      if (add_phi_adv) then
        if (k==k_end) then
          rhs(i, k+1) = rhs(i,k+1) - dts* (0.5*rdx/msfty(i,j))*(                             &
                 ( (cfn*u(i+1,k,j)+cfn1*u(i+1,k-1,j))  &
                  +(cfn*u(i  ,k,j)+cfn1*u(i,k-1,j)) )*(1./60.)*(  &
                   45.*(ph(i+1,k-1,j)-ph(i-1,k-1,j))   &
                   -9.*(ph(i+2,k-1,j)-ph(i-2,k-1,j))   &
                      +(ph(i+3,k-1,j)-ph(i-3,k-1,j))   &
                  +45.*(phb(i+1,k-1,j)-phb(i-1,k-1,j)) &
                   -9.*(phb(i+2,k-1,j)-phb(i-2,k-1,j)) &
                      +(phb(i+3,k-1,j)-phb(i-3,k-1,j))  )     )
        else 
          rhs(i,k+1) = rhs(i,k+1) + dts*(-0.25*rdx/msfty(i,j))*( &
                 ( u(i+1,k+1,j)+u(i+1,k,j) &
                  +(u(i,k+1,j  )+u(i,k,j  )) )* (1./60.)*( &
                   45.*(ph(i+1,k+1,j)-ph(i-1,k+1,j)) &
                   -9.*(ph(i+2,k+1,j)-ph(i-2,k+1,j)) &
                      +(ph(i+3,k+1,j)-ph(i-3,k+1,j)) &
                  +45.*(phb(i+1,k+1,j)-phb(i-1,k+1,j)) &
                   -9.*(phb(i+2,k+1,j)-phb(i-2,k+1,j)) &
                      +(phb(i+3,k+1,j)-phb(i-3,k+1,j))  )   )      
        endif
      endif

      !PZ added
      !IF (iteration ==1 .and. i<=5 .and. k==1 .and. j==1) THEN
      !  write(*,*) 'ph_tend=',ph_tend(i,k+1,j), 'rhs=', rhs(i,k+1)
      !  write(*,*) 'wdwn', wdwn(i,k+1), 'ww', .5*(ww(i,k+1,j)+ww(i,k,j)), 'ph_1(i,2,1)', ph_1(i,k+1,j), 'ph_1(i,1,1)', ph_1(i,k,j),'phb(i,2,1)', phb(i,k+1,j),'phb(i,1,1)', phb(i,k,j), 'msfty=', msfty(i,j) 
      !ENDIF

    ENDDO
    ENDDO

!   building up RHS of phi equation
!   ph_tend contains terms 1 and 2; now adding 3 and 4 in stages:
!   here rhs = delta t [ph_tend + ~g*w/2 - ~ww * partial d phi/dz]
    DO k=2,k_end
    DO i=i_start, i_end
       rhs(i,k) = rhs(i,k)-dts*( fnm(k)*wdwn(i,k+1)  &
                                +fnp(k)*wdwn(i,k  ) ) 
    ENDDO
    ENDDO

!  NOTE:  phi'' is not coupled with the map-scale factor  (1/m), 
!         but it's tendency is, so must multiply by msft here

!   Comments on map scale factors:
!   building up RHS of phi equation
!   ph_tend contains terms 1 and 2; now adding 3 and 4 in stages:
!   here rhs = ph_previous + (msft/mu)*[rhs_previous - ~ww * delta t * 
!                                                      partial d phi/dz]
!            = ph_previous + (msft*delta t/mu)*[ph_tend + ~g*w/2 - ~ww *
!                                                         partial d phi/dz]
   
    !PZ added
    !IF (iteration == num_sound_steps) phb_rk1 = ph

    DO k=2,k_end+1
    DO i=i_start, i_end
      rhs(i,k) = ph(i,k,j) + msfty(i,j)*rhs(i,k)*mut_inv(i)
      if(top_lid .and. k.eq.k_end+1)rhs(i,k)=0.
      
      !PZ added for moving bottom: specify BC of ph
      IF (k==2 .and. move ) THEN
        x_dist = (i-0.5)/rdx
        SELECT CASE (test_wave)
        CASE (1)
          ! rhs(i,k-1) = ph(i,k-1,j) + dts*(-amp*omega*sin(k_wave*x_dist)*sin(omega*(t_old + dts*iteration) ))*g
           rhs(i,k-1) = amp*sin(k_wave*x_dist)*cos(omega*(t_old +dts*iteration) )*g - amp*sin(k_wave*x_dist)*g - ph_1(i, k-1, j) 
        CASE (2)
          IF (set_w) THEN
              rhs(i, k-1) = ph(i,1,j) + dts*1/muts(i,j)*(g*(1-epssm)*0.5*w(i,1,j) &
                      -.5*rdx*(                        &
                    (phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)) )  &
           *muu(i+1,j)*(cf1*u_save(i+1,1,j)+cf2*u_save(i+1,2,j)+cf3*u_save(i+1,3,j)) &
                   +(phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j)) ) &
           *muu(i,j)*(cf1*u_save(i  ,1,j)+cf2*u_save(i,2,j)+cf3*u_save(i,3,j))) + g*mut(i,j)*w_save(i,1,j) & !-muts(i,j)*g*(amp*sin(k_wave*x_dist -omega*(t_old + dts*iteration))*tanh(2.0*(t_old + dts*iteration)/tau) - amp*sin(k_wave*x_dist -omega*(t_old + dts*(iteration-1)))*tanh(2.0*(t_old + dts*(iteration-1))/tau))/dts )
                      -.5*rdy*(                        &
                    (phb(i,1,j+1)+ph_1(i,1,j+1) - (phb(i,1,j)+ph_1(i,1,j)) ) &
           *muv(i,j+1)*(cf1*v_save(i,1,j+1)+cf2*v_save(i,2,j+1)+cf3*v_save(i,3,j+1)) &
                    +(phb(i,1,j)+ph_1(i,1,j) - (phb(i,1,j-1)+ph_1(i,1,j-1))) &
           *muv(i,j)*(cf1*v_save(i,1,j  )+cf2*v_save(i,2,j  )+cf3*v_save(i,3,j  ))  ) )
           
              IF (add_phi_adv) THEN
                rhs(i, k-1) = rhs(i,k-1) - dts*1/muts(i,j)*.5*rdx*(  &
                (phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)) )*(cf1*u(i +1 ,1,j)+cf2*u(i+1,2,j)+cf3*u(i+1,3,j)) &
                + (phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j)) ) &
                  *(cf1*u(i,1,j)+cf2*u(i,2,j)+cf3*u(i,3,j)) )  
              ENDIF

          ELSE !not setting w
            IF (ramp) THEN
              IF (iteration == num_sound_steps) then
              rhs(i,k-1) =  amp*sin(k_wave*x_dist -omega*(t_old + dts*iteration))*tanh(2.0*(t_old + dts*iteration)/tau)*g - ph(i,k-1,j) !- ph_1(i, k-1, j) 
              ENDIF
            ELSE 
              rhs(i,k-1) = amp*sin(k_wave*x_dist -omega*(t_old + dts*iteration))*g - amp*sin(k_wave*x_dist)*g - ph_1(i, k-1, j)
            END IF
          ENDIF

        CASE (3)
          rhs(i,k-1) = ph(i,k-1,j) + dts*g*amp*exp(-2*config_flags%khdif*k_wave**2*(t_old + dts*iteration) )*(-2*config_flags%khdif*k_wave**2*sin(k_wave*x_dist +omega*(t_old + dts*iteration)) + omega*cos(k_wave*x_dist +omega*(t_old + dts*iteration)))
        CASE (4)
          rhs(i,k-1) = ph(i,k-1,j) + dts*g*amp*exp(-2*config_flags%khdif*k_wave**2*(t_old + dts*iteration))*(-2*config_flags%khdif*k_wave**2*sin(k_wave*x_dist))
        CASE (5)  !Flat plate moving up
          IF (set_w) THEN
            rhs(i, k-1) = ph(i,1,j) + dts*1/muts(i,j)*(g*(1-epssm)*0.5*w(i,1,j)&
                      -.5*rdx*(                        &
                    (phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)) )  &
           *muu(i+1,j)*(cf1*u_save(i+1,1,j)+cf2*u_save(i+1,2,j)+cf3*u_save(i+1,3,j))&
                   +(phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j)) ) &
           *muu(i,j)*(cf1*u_save(i  ,1,j)+cf2*u_save(i,2,j)+cf3*u_save(i,3,j))) + g*mut(i,j)*w_save(i,1,j) - muts(i,j)*g*amp) 
          ELSE
            rhs(i,k-1) = ( amp*(t_old+dts*iteration))*g - ph_1(i, k-1, j)
          ENDIF

        CASE (6)  !standing wave slowly moving up
          rhs(i,k-1) = (amp*cos(k_wave*x_dist)*sin(omega*(t_old + dts*iteration) ))*g  - ph_1(i, k-1, j) !amp*omega*cos(k_wave*x_dist)*sin(omega*t_rk ))*g
          !rhs(i,k-1) = ph(i,k-1,j) + dts*(amp*omega*cos(k_wave*x_dist)*cos(omega*(t_old + dts*iteration) ))*g
        END SELECT
      !ELSE IF (k==2 .and. move ) THEN
      !  rhs(i,k-1) = ph(i,k-1,j)
      ENDIF
    ENDDO
    ENDDO

!  lower boundary condition on w

!   Comments on map scale factors:
!   Chain rule: if Z=Z(X,Y) [true at the surface] then
!      dZ/dt = dZ/dX * dX/dt + dZ/dY * dY/dt, U=dX/dt, V=dY/dt
!   using capitals to denote actual values
!   In mapped values, u=U, v=V, z=Z, 1/dX=mx/dx, 1/dY=my/dy
!      w = dz/dt = mx u dz/dx + my v dz/dy
!   [where dz/dx is just the surface height change between x
!    gridpoints, and dz/dy is the change between y gridpoints]
!   [cf1, cf2 and cf3 do vertical weighting of u or v values nearest
!    the surface]

    !KAL add no-slip boundary condition
    IF (config_flags%no_slip_zs) THEN
       DO i=i_start,i_end
          w(i,1,j)=0.0
       ENDDO
    ELSE
     IF (.not. move) then  !PZ added
       DO i=i_start,i_end
!          w(i,1,j)=                                           &
!                   msfty(i,j)*.5*rdy*(                        &
!                            (ht(i,j+1)-ht(i,j  ))             &
!           *(cf1*v(i,1,j+1)+cf2*v(i,2,j+1)+cf3*v(i,3,j+1))    &
!                           +(ht(i,j  )-ht(i,j-1))             &
!           *(cf1*v(i,1,j  )+cf2*v(i,2,j  )+cf3*v(i,3,j  ))  ) &
!                  +msftx(i,j)*.5*rdx*(                        &
!                            (ht(i+1,j)-ht(i,j  ))             &
!           *(cf1*u(i+1,1,j)+cf2*u(i+1,2,j)+cf3*u(i+1,3,j))    &
!                          +(ht(i,j  )-ht(i-1,j))             &
!          *(cf1*u(i  ,1,j)+cf2*u(i  ,2,j)+cf3*u(i  ,3,j))  ) 
           w(i,1,j)=                                           &!PZ modified
                   msfty(i,j)*.5*rdy*(                        &
                            (ht(i,j+1)-ht(i,j  ))             &
           *(cf1*v(i,1,j+1)+cf2*v(i,2,j+1)+cf3*v(i,3,j+1))    &
                           +(ht(i,j  )-ht(i,j-1))             &
           *(cf1*v(i,1,j  )+cf2*v(i,2,j  )+cf3*v(i,3,j  ))  ) &
                  +msftx(i,j)*.5*rdx*(                        &
                            (ht(i+1,j)-ht(i,j  ))             &
           *(cf1*u(i+1,1,j)+cf2*u(i+1,2,j)+cf3*u(i+1,3,j))    &
                          +(ht(i,j  )-ht(i-1,j))             &
          *(cf1*u(i  ,1,j)+cf2*u(i  ,2,j)+cf3*u(i  ,3,j))  )
           if (i<=5 .and. j==1) write(*,*) 'not move, w=',w(i,1,j) 
       ENDDO

     ELSE ! moving bottom
       DO i=i_start,i_end
          x_dist = (i-0.5)/rdx
          !w(i,1,j)= w(i,1,j) - muts(i,j)*w_wave(i, j)/ msfty(i,j)
          SELECT CASE (test_wave)
          CASE (1)
            w(i,1,j)=  muts(i,j)* (-amp*omega*sin(k_wave*x_dist)*sin(omega*(t_old + dts*iteration) )) / msfty(i,j) - mut(i,j)*w_save(i,1,j)/msfty(i,j) + msfty(i,j)*.5*rdy*(                        &
                            (ht(i,j+1)-ht(i,j  ))             &
           *(cf1*v(i,1,j+1)+cf2*v(i,2,j+1)+cf3*v(i,3,j+1))    &
                           +(ht(i,j  )-ht(i,j-1))             &
           *(cf1*v(i,1,j  )+cf2*v(i,2,j  )+cf3*v(i,3,j  ))  ) &
                  +msftx(i,j)*.5*rdx*(                        &
                            (ht(i+1,j)-ht(i,j  ))             &
           *(cf1*u(i+1,1,j)+cf2*u(i+1,2,j)+cf3*u(i+1,3,j))    &
                          +(ht(i,j  )-ht(i-1,j))             &
          *(cf1*u(i  ,1,j)+cf2*u(i  ,2,j)+cf3*u(i  ,3,j))  ) 

          CASE (2)
            IF (ramp) THEN
         
            IF (set_w) THEN
              w(i,1,j)=  (muts(i,j)*(amp*sin(k_wave*x_dist -omega*(t_old + dts*iteration))*tanh(2.0*(t_old + dts*iteration)/tau) - amp*sin(k_wave*x_dist -omega*(t_old + dts*(iteration-1)))*tanh(2.0*(t_old + dts*(iteration-1))/tau))/dts/ msfty(i,j)  & !(muts(i,j)* (-(amp*omega*cos(k_wave*x_dist -omega*(t_old + dts*iteration))*tanh(2.0*(t_old + dts*iteration)/tau) - 2*amp/tau*sin(k_wave*x_dist - omega*(t_old + dts*iteration))*(1 - tanh(2*(t_old + dts*iteration)/tau)**2) ) ) / msfty(i,j)    & 
         
                  - mut(i,j)*w_save(i,1,j) &
                  +.5*rdx*(                        &
                    1/g*(phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)) ) &
           *muu(i+1,j)*(cf1*u_save(i+1,1,j)+cf2*u_save(i+1,2,j)+cf3*u_save(i+1,3,j)) &
                   +1/g*(phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j)) ) &
           *muu(i,j)*(cf1*u_save(i  ,1,j)+cf2*u_save(i  ,2,j)+cf3*u_save(i,3,j))  ) &
                  +.5*rdy*(                        &
                    1/g*(phb(i,1,j+1)+ph_1(i,1,j+1) - (phb(i,1,j)+ph_1(i,1,j)) ) &
           *muv(i,j+1)*(cf1*v_save(i,1,j+1)+cf2*v_save(i,2,j+1)+cf3*v_save(i,3,j+1))    &
                    +1/g*(phb(i,1,j)+ph_1(i,1,j) - (phb(i,1,j-1)+ph_1(i,1,j-1))) &
           *muv(i,j)*(cf1*v_save(i,1,j  )+cf2*v_save(i,2,j  )+cf3*v_save(i,3,j  ))  ) &
           -(1-epssm)*0.5*w(i,1,j) )*2/(1+epssm)
              !IF (iteration == num_sound_steps) then
                !if (i==1 .and. j==1) write(*,*) 'num sound steps=', num_sound_steps
                !w(i,1,j) = w(i,1,j) + 2/(1+epssm)*muts(i,j)*(amp*sin(k_wave*x_dist -omega*(t_old + dts*iteration))*tanh(2.0*(t_old + dts*iteration)/tau) - amp*sin(k_wave*x_dist-omega*(t_old+dts*(iteration-1)))*tanh(2.0*(t_old + dts*(iteration-1))/tau))/dts/ msfty(i,j)
    !          ENDIF
                         
            ELSE
              w(i,1,j)= 1/msfty(i,j)*(muts(i,j)*((rhs(i,1) - ph(i,1,j))/dts/g) & 
                  - mut(i,j)*w_save(i,1,j) &
                  +.5*rdx*(                        &
                    1/g*(phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)) )  &
           *muu(i+1,j)*(cf1*u_save(i+1,1,j)+cf2*u_save(i+1,2,j)+cf3*u_save(i+1,3,j))    &
                   +1/g*(phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j)) ) &
           *muu(i,j)*(cf1*u_save(i  ,1,j)+cf2*u_save(i  ,2,j)+cf3*u_save(i  ,3,j))  ) &
           -(1-epssm)*0.5*w(i,1,j) )*2/(1+epssm)  
            ENDIF !end if set_w

          ELSE  ! not doing ramping
            IF (set_w) THEN
              w(i,1,j)=  (muts(i,j)* ( -amp*omega*cos(k_wave*x_dist -omega*(t_old + dts*iteration) ))/ msfty(i,j) - mut(i,j)*w_save(i,1,j)/msfty(i,j) + msfty(i,j)*.5*rdy*(                        &
                            (ht(i,j+1)-ht(i,j  ))             &
           *(cf1*v(i,1,j+1)+cf2*v(i,2,j+1)+cf3*v(i,3,j+1))    &
                           +(ht(i,j  )-ht(i,j-1))             &
           *(cf1*v(i,1,j  )+cf2*v(i,2,j  )+cf3*v(i,3,j  ))  ) &
  +.5*rdx*(                        &
                    1/g*(phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)) )&
           *muu(i+1,j)*(cf1*u_save(i+1,1,j)+cf2*u_save(i+1,2,j)+cf3*u_save(i+1,3,j)) &
                   +1/g*(phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j)) ) &
           *muu(i,j)*(cf1*u_save(i  ,1,j)+cf2*u_save(i  ,2,j)+cf3*u_save(i,3,j))) &
           -(1-epssm)*0.5*w(i,1,j) )*2/(1+epssm)
            ELSE
              w(i,1,j)= 1/msfty(i,j)*(muts(i,j)*((rhs(i,1) - ph(i,1,j))/dts/g) &
                  - mut(i,j)*w_save(i,1,j) &
                  +.5*rdx*(                        &
                    1/g*(phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)) ) &
           *muu(i+1,j)*(cf1*u_save(i+1,1,j)+cf2*u_save(i+1,2,j)+cf3*u_save(i+1,3,j)) &
                   +1/g*(phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j)) ) &
           *muu(i,j)*(cf1*u_save(i  ,1,j)+cf2*u_save(i  ,2,j)+cf3*u_save(i,3,j))  ) &
           -(1-epssm)*0.5*w(i,1,j) )*2/(1+epssm)
            ENDIF !end if set_w
          ENDIF ! end if ramping
             
              if (i==1 .and. j==1) write(*,*) 'ori u_surf=', cf1*u(i+1,1,j)+cf2*u(i+1,2,j)+cf3*u(i+1,3,j), 'ana u_surf=', muts(i,j)*amp*omega*sin(k_wave*x_dist -omega*(t_old + dts*iteration) ) - mut(i,j)*u_save(i,1,j)
   
          CASE (3)
            w(i,1,j)= muts(i,j)* (amp*omega*exp(-2*config_flags%khdif*k_wave**2*(t_old+ dts*iteration) )*cos(k_wave*x_dist + omega*(t_old+ dts*iteration) ) &
            -2*config_flags%khdif*k_wave**2*amp*exp(-2*config_flags%khdif*k_wave**2*(t_old+ dts*iteration) )*sin(k_wave*x_dist + omega*(t_old+ dts*iteration)) ) / msfty(i,j)
          
          CASE (4)
            w(i,1,j)= w(i,1,j) + muts(i,j)* (-2*config_flags%khdif*k_wave**2*amp*exp(-2*config_flags%khdif*k_wave**2*(t_old+ dts*iteration) )*sin(k_wave*x_dist) )/ msfty(i,j) 
          CASE (5)
            if (set_w) then
              if (i==1 .and. j==1 ) write (*,*) 'before set w, w=', w(i,1,j), 'muts=', muts(i,j), 'mut', mut(i,j), 'udhdx=', .5*rdx*(                        &
                    1/g*(phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)))&
           *muu(i+1,j)*(cf1*u_save(i+1,1,j)+cf2*u_save(i+1,2,j)+cf3*u_save(i+1,3,j))&
                   +1/g*(phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j)))&
           *muu(i,j)*(cf1*u_save(i  ,1,j)+cf2*u_save(i,2,j)+cf3*u_save(i,3,j))) 

              w(i,1,j) = (muts(i,j)*amp - mut(i,j)*w_save(i,1,j) &
                  +.5*rdx*(                        &
                    1/g*(phb(i+1,1,j)+ph_1(i+1,1,j) - (phb(i,1,j)+ph_1(i,1,j)) )&
           *muu(i+1,j)*(cf1*u_save(i+1,1,j)+cf2*u_save(i+1,2,j)+cf3*u_save(i+1,3,j)) &
                   +1/g*(phb(i,1,j)+ph_1(i,1,j) - (phb(i-1,1,j)+ph_1(i-1,1,j))) &
           *muu(i,j)*(cf1*u_save(i  ,1,j)+cf2*u_save(i  ,2,j)+cf3*u_save(i,3,j))) &
           -(1-epssm)*0.5*w(i,1,j) )*2/(1+epssm)
              if (i==1 .and. j==1 ) write (*,*) 'set w to:', w(i,1,j)
            else
            w(i,1,j) = muts(i,j)*amp/ msfty(i,j) - muts(i,j)*w_save(i,1,j)/msfty(i,j)
            endif

          CASE (6)
            w(i,1,j) = muts(i,j)*(amp*omega *cos(k_wave*x_dist)*cos(omega*(t_old + dts*iteration) )) / msfty(i,j) - mut(i,j)*w_save(i,1,j)/msfty(i,j) 
          END SELECT

       ENDDO

     ENDIF !PZ: end of if not move
    ENDIF ! end of if no-slip
 

! Jammed 3 doubly nested loops over k/i into 1 for slight improvement
! in efficiency.  No change in results (bit-for-bit). JM 20040514
! (left a blank line where the other two k/i-loops were)
!
!   above surface, begin by adding delta t * previous (coupled) w tendency
!    if (j==1 .and. iteration ==1) then
!      write(*,*) 'before updating w, rw_tend=', rw_tend(1,2,1)
!      write(*,*) 'rhs(k+1)-rhs(k)=', rhs(1,3  )-rhs(1,2), rhs(1,2)-rhs(1,1)
!      write(*,*) 'rhs(2)=', rhs(1,2),'rhs(1)=', rhs(1,1)
!      write(*,*) 't2ave=', t_2ave(1,1,1), t_2ave(2,1,1), t_2ave(1,2,1)
!      write(*,*) 'muave=', muave(1,1), muave(2,1)
!      write(*,*) 'last term in w eq=', dts*g*msft_inv(i)*rdn(2)*(c2a(i,2  ,j)*alt(i,2  ,j)*t_2ave(i,2  ,j)            &
!             -c2a(i,k,j)*alt(i,1,j)*t_2ave(i,1,j))
!      write(*,*) 'rhs term=', -cqw(1,2,1)*.5*dts*g*mut_inv(1)*rdn(2)*c2a(1,1,j)*rdnw(1)*((1.+epssm)*(rhs(1,2    )-rhs(1,1  ))), -cqw(1,3,1)*.5*dts*g*mut_inv(1)*rdn(3)*c2a(1,2,j)*rdnw(2)*((1.+epssm)*(rhs(1,3    )-rhs(1,2  )))
!    endif

    DO k=2,k_end
      DO i=i_start, i_end
        w(i,k,j)=w(i,k,j)+dts*rw_tend(i,k,j)                       &
                 + msft_inv(i)*cqw(i,k,j)*(                        &
            +.5*dts*g*mut_inv(i)*rdn(k)*                           &
             (c2a(i,k  ,j)*rdnw(k  )                               &
        *((1.+epssm)*(rhs(i,k+1  )-rhs(i,k    ))                   &
         +(1.-epssm)*(ph(i,k+1,j)-ph(i,k  ,j)))                    &
             -c2a(i,k-1,j)*rdnw(k-1)                               &
        *((1.+epssm)*(rhs(i,k    )-rhs(i,k-1  ))                   &
         +(1.-epssm)*(ph(i,k  ,j)-ph(i,k-1,j)))))                  &

                +dts*g*msft_inv(i)*(rdn(k)*                        &
             (c2a(i,k  ,j)*alt(i,k  ,j)*t_2ave(i,k  ,j)            &
             -c2a(i,k-1,j)*alt(i,k-1,j)*t_2ave(i,k-1,j))           &
              - muave(i,j))
       !PZ added
       !if ((k<=3) .and. j==1 .and. i==1) THEN
       !  WRITE(*,*) 'k=',k, 'rhs term:', rdnw(k)*(1.+epssm)*(rhs(i,k+1  )-rhs(i,k)), &
       !  'ph term:', rdnw(k)*(1.-epssm)*(ph(i,k+1,j)-ph(i,k  ,j)), &
!         'alt term',dts*g*msft_inv(i)*rdn(k)*                        &
!             (c2a(i,k  ,j)*alt(i,k  ,j)*t_2ave(i,k  ,j)            &
!             -c2a(i,k-1,j)*alt(i,k-1,j)*t_2ave(i,k-1,j)), &   
!         '2nd rhs:', rdn(k)*                           &
!             (c2a(i,k  ,j)*rdnw(k  )                               &
!        *(1.+epssm)*(rhs(i,k+1  )-rhs(i,k    ))                   &
!             -c2a(i,k-1,j)*rdnw(k-1)                               &
!        *(1.+epssm)*(rhs(i,k    )-rhs(i,k-1  ))) , &
!         '2nd ph:', rdn(k)*                           &
!             (c2a(i,k  ,j)*rdnw(k  )                               &
!        *(1.+epssm)*(ph(i,k+1, j ) - ph(i,k,j    ))                   &
!             -c2a(i,k-1,j)*rdnw(k-1)                               &
!        *(1.+epssm)*(ph(i,k,j)-ph(i,k-1, j )))
!       ENDIF

      ENDDO
    ENDDO

    if (j==1 .and. iteration ==1 ) write(*,*) 'after w first loop, w(1,2,1)=', w(1,2,1), w(1,3,1)   !PZ added

    K=k_end+1

    DO i=i_start, i_end
       w(i,k,j)=w(i,k,j)+dts*rw_tend(i,k,j)                         &
           +msft_inv(i)*(                                        &
         -.5*dts*g*mut_inv(i)*rdnw(k-1)**2*2.*c2a(i,k-1,j)            &
             *((1.+epssm)*(rhs(i,k  )-rhs(i,k-1  ))                 &
              +(1.-epssm)*(ph(i,k,j)-ph(i,k-1,j)))                  &
         -dts*g*(2.*rdnw(k-1)*                                      &
                   c2a(i,k-1,j)*alt(i,k-1,j)*t_2ave(i,k-1,j)        &
              + muave(i,j)) )
       if(top_lid)w(i,k,j) = 0.
    ENDDO

    DO k=2,k_end+1
    DO i=i_start, i_end
       w(i,k,j)=(w(i,k,j)-a(i,k,j)*w(i,k-1,j))*alpha(i,k,j)
    ENDDO
    ENDDO

    DO k=k_end,2,-1
    DO i=i_start, i_end
       w (i,k,j)=w (i,k,j)-gamma(i,k,j)*w(i,k+1,j)
    ENDDO
    ENDDO
     
     !k = k_end+1
     !if (move) then  !PZ added
     !    w(i,k,j) = (rhs(i,k  ) - ph(i,k,j))/g/dts*muts(i,j)/ msfty(i,j) - mut(i,j)*w_save(i,k,j)/msfty(i,j)
     !endif
 
    IF (config_flags%damp_opt .eq. 3) THEN
      DO k=k_end+1,2,-1
      DO i=i_start, i_end
          htop=(ph_1(i,k_end+1,j)+phb(i,k_end+1,j))/g
          hk=(ph_1(i,k,j)+phb(i,k,j))/g
          hbot=htop-hdepth
          dampwt(k) = 0.
          if(hk .ge. hbot)then
            dampwt(k) = dampmag*sin(0.5*pi*(hk-hbot)/hdepth)*sin(0.5*pi*(hk-hbot)/hdepth)
          endif
          w(i,k,j) = (w(i,k,j) - dampwt(k)*mut(i,j)*w_save(i,k,j))/(1.+dampwt(k))
      ENDDO
      ENDDO
    ENDIF

    if (move .and. j==1) then
      write(*,*) 't=', t_old+ dts*iteration 
      write(*,*) 'w1=',w(1,1,1),muts(1,1), mut(1,1)*w_save(1,1,1) 
      write(*,*) 'w2=', w(1,2,1), mut(1,1)*w_save(1,2,1)
      write(*,*) 'w3=', w(1,3,1), mut(1,1)*w_save(1,3,1)
    endif

    DO k=k_end+1,2,-1
    DO i=i_start, i_end
       ph(i,k,j) = rhs(i,k)+msfty(i,j)*.5*dts*g*(1.+epssm) &
                      *w(i,k,j)/muts(i,j)
   
       !PZ added
       IF (k==2 .and. move) THEN
           IF (set_w) THEN
             ph(i,k-1,j) = rhs(i,k-1) +msfty(i,j)*.5*dts*g*(1.+epssm) &
                      *w(i,k-1,j)/muts(i,j)
           ELSE
             ph(i,k-1,j) = rhs(i,k-1) 
           END IF
       ENDIF
    ENDDO
    ENDDO
 
    !PZ added
    !IF (move .and. j==1) then
    !    Do i = 1, 3
    !      write(*,*) 'i=', i,'ph=', ph(i,1,1)/g, ph(i,2,1)/g, ph(i,3,1)/g
    !    ENDDO
    !ENDIF
      
      
    

  ENDDO j_loop_w
  
  if (write_acous) then
  IF (its==1) THEN
    OPEN(UNIT=100, FILE=TRIM('acoustic_w001'), STATUS='unknown',POSITION='append', FORM='formatted')
    profile_format = '(000(f14.10,1x))'
    WRITE(profile_format(2:4),'(I3.3)') kde
    WRITE(UNIT=100,FMT=profile_format) (w(1,1:kde, 1) + mut(1,1)*w_save(1,1:kde,1))/muts(1,1)   
    CLOSE(UNIT=100)
  
    OPEN(UNIT=100, FILE=TRIM('acoustic_w010'), STATUS='unknown',POSITION='append',FORM='formatted')
    WRITE(UNIT=100,FMT=profile_format) (w(10,1:kde, 1) + mut(10,1)*w_save(10,1:kde,1))/muts(10,1)
    CLOSE(UNIT=100)

    OPEN(UNIT=100, FILE=TRIM('acoustic_ph001'), STATUS='unknown',POSITION='append', FORM='formatted')
    WRITE(UNIT=100,FMT=profile_format) ph(1,1:kde, 1) + ph_1(1,1:kde,1)
    CLOSE(UNIT=100)

    OPEN(UNIT=100, FILE=TRIM('acoustic_ph010'), STATUS='unknown',POSITION='append',FORM='formatted')
    WRITE(UNIT=100,FMT=profile_format) ph(10,1:kde, 1) + ph_1(10,1:kde,1)
    CLOSE(UNIT=100)

  ENDIF
  endif
END SUBROUTINE advance_w

!---------------------------------------------------------------------

SUBROUTINE sumflux ( ru, rv, ww,                             &
                     u_lin, v_lin, ww_lin,                   &
                     muu, muv,                               &
                     ru_m, rv_m, ww_m, epssm,                &
                     msfux, msfuy, msfvx, msfvx_inv, msfvy,  &
                     iteration , number_of_small_timesteps,  &
                     ids,ide, jds,jde, kds,kde,              &
                     ims,ime, jms,jme, kms,kme,              &
                     its,ite, jts,jte, kts,kte              )


  IMPLICIT NONE  ! religion first

! declarations for the stuff coming in

  INTEGER,      INTENT(IN   )    :: number_of_small_timesteps
  INTEGER,      INTENT(IN   )    :: iteration
  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),  INTENT(IN   ) :: ru, &
                                                                rv, &
                                                                ww, &
                                                                u_lin,  &
                                                                v_lin,  &
                                                                ww_lin


  REAL, DIMENSION(ims:ime, kms:kme, jms:jme) , INTENT(INOUT) :: ru_m, &
                                                                rv_m, &
                                                                ww_m
  REAL, DIMENSION(ims:ime, jms:jme) , INTENT(IN   ) :: muu, muv,      &
                                                       msfux, msfuy,  &
                                                       msfvx, msfvy, msfvx_inv

  INTEGER :: mini, minj, mink


  REAL, INTENT(IN   )  ::  epssm
  INTEGER   :: i,j,k


!<DESCRIPTION>
!
!  update the small-timestep time-averaged mass fluxes;  these
!  are needed for consistent mass-conserving scalar advection.
!
!</DESCRIPTION>

    IF (iteration == 1 )THEN
      DO  j = jts, jte
      DO  k = kts, kte
      DO  i = its, ite
        ru_m(i,k,j)  = 0.
        rv_m(i,k,j)  = 0.
        ww_m(i,k,j)  = 0.
      ENDDO
      ENDDO
      ENDDO
    ENDIF

  mini = min(ide-1,ite)
  minj = min(jde-1,jte)
  mink = min(kde-1,kte)


    DO  j = jts, minj
    DO  k = kts, mink
    DO  i = its, mini
      ru_m(i,k,j)  = ru_m(i,k,j) + ru(i,k,j)
      rv_m(i,k,j)  = rv_m(i,k,j) + rv(i,k,j)
      ww_m(i,k,j)  = ww_m(i,k,j) + ww(i,k,j)
    ENDDO
    ENDDO
    ENDDO
 
    IF (ite .GT. mini) THEN
      DO  j = jts, minj
      DO  k = kts, mink
      DO  i = mini+1, ite
        ru_m(i,k,j)  = ru_m(i,k,j) + ru(i,k,j)
      ENDDO
      ENDDO
      ENDDO
    END IF
    IF (jte .GT. minj) THEN
      DO  j = minj+1, jte
      DO  k = kts, mink
      DO  i = its, mini
	rv_m(i,k,j)  = rv_m(i,k,j) + rv(i,k,j)
      ENDDO
      ENDDO
      ENDDO
    END IF
    IF ( kte .GT. mink) THEN
      DO  j = jts, minj
      DO  k = mink+1, kte
      DO  i = its, mini
        ww_m(i,k,j)  = ww_m(i,k,j) + ww(i,k,j)
      ENDDO
      ENDDO
      ENDDO
    END IF

  IF (iteration == number_of_small_timesteps) THEN

    DO  j = jts, minj
    DO  k = kts, mink
    DO  i = its, mini
      ru_m(i,k,j)  = ru_m(i,k,j) / number_of_small_timesteps   &
                     + muu(i,j)*u_lin(i,k,j)/msfuy(i,j)
      rv_m(i,k,j)  = rv_m(i,k,j) / number_of_small_timesteps   &
                     + muv(i,j)*v_lin(i,k,j)*msfvx_inv(i,j) 
      ww_m(i,k,j)  = ww_m(i,k,j) / number_of_small_timesteps   &
                     + ww_lin(i,k,j) 
    ENDDO
    ENDDO
    ENDDO


    IF (ite .GT. mini) THEN
      DO  j = jts, minj
      DO  k = kts, mink
      DO  i = mini+1, ite
        ru_m(i,k,j)  = ru_m(i,k,j) / number_of_small_timesteps   &
                     + muu(i,j)*u_lin(i,k,j)/msfuy(i,j)
      ENDDO
      ENDDO
      ENDDO
    END IF
    IF (jte .GT. minj) THEN
      DO  j = minj+1, jte
      DO  k = kts, mink
      DO  i = its, mini
	rv_m(i,k,j)  = rv_m(i,k,j) / number_of_small_timesteps   &
                     + muv(i,j)*v_lin(i,k,j)*msfvx_inv(i,j)
      ENDDO
      ENDDO
      ENDDO
    END IF
    IF ( kte .GT. mink) THEN
      DO  j = jts, minj
      DO  k = mink+1, kte
      DO  i = its, mini
        ww_m(i,k,j)  = ww_m(i,k,j) / number_of_small_timesteps   &
                     + ww_lin(i,k,j)
      ENDDO
      ENDDO
      ENDDO
    END IF

  ENDIF


END SUBROUTINE sumflux 

!---------------------------------------------------------------------

  SUBROUTINE init_module_small_step
  END SUBROUTINE init_module_small_step

END MODULE module_small_step_em
