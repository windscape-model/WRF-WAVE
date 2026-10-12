!WRF:MEDIATION_LAYER:SOLVER

MODULE module_after_all_rk_steps

CONTAINS

   !  This subroutine is called once per domain per time step.  It is outside
   !  of and after the end of the Runge-Kutta time steps, after the calls to 
   !  the explicit moisture driver, and after the polar filtering calls.  The
   !  variables in here are all up-to-date with the end of this current time 
   !  step.


   SUBROUTINE after_all_rk_steps ( grid, config_flags,                  &
                                   moist, chem, tracer, scalar,         &
                                   th_phy, pi_phy, p_phy, rho_phy,      & 
                                   p8w, t8w, dz8w,                      &
                                   curr_secs2,                          &
                                   diag_flag,                           &
                                   ids,  ide,  jds,  jde,  kds,  kde,   &
                                   ims,  ime,  jms,  jme,  kms,  kme,   &
                                   ips,  ipe,  jps,  jpe,  kps,  kpe,   &
                                   imsx, imex, jmsx, jmex, kmsx, kmex,  &
                                   ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                                   imsy, imey, jmsy, jmey, kmsy, kmey,  &
                                   ipsy, ipey, jpsy, jpey, kpsy, kpey   )


      !=============================================================
      !  USE Association for Generic WRF Infrastructure
      !=============================================================

      !  Pick up the number of members for each of the 4d arrays - for declaration purposes.

      USE module_state_description, ONLY: num_moist, num_chem, num_tracer, num_scalar

      !  This gives us the type definition for grid (domain)

      USE module_domain, ONLY : domain

      !  All of the information from the namelist is in config_flags.  The
      !  type declaration for this puppy must be available.

      USE module_configure, ONLY : grid_config_rec_type

#ifdef DM_PARALLEL
      !  Ensure some of the fancy diagnostics variables that need to
      !  talk to other patches can do so.

      USE module_dm, ONLY : &
                  local_communicator, mytask, ntasks, ntasks_x, ntasks_y                   &
                 ,local_communicator_periodic, wrf_dm_maxval

      USE module_comm_dm, ONLY : &
                  halo_em_phys_w_sub
#endif

      !=============================================================
      !  USE Association for the Diagnostic Packages
      !=============================================================
      
      USE module_diagnostics_driver, ONLY : diagnostics_driver


      IMPLICIT NONE


      !=============================================================
      !  Subroutine Arguments
      !=============================================================

      !  Arguments passed in.  All of the diagnostics are part of the grid structure, so
      !  even though we are not changing any of the fundamental variables, we are computing
      !  the diagnostics.  Therefore grid is INOUT.

      TYPE ( domain ), INTENT(INOUT) :: grid

      !  We are not changing any of the namelist settings.

      TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags

      !  The 4d arrays are input only, no mods to them.

      REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_moist ) , INTENT(IN) :: moist
      REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_chem  ) , INTENT(IN) :: chem
      REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer) , INTENT(IN) :: tracer
      REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar) , INTENT(IN) :: scalar

      !  A few handy 3d arrays computed for the physics scheme: pressure (Pa) and
      !  temperature (K), on both half (_phy) and full levels.

      REAL , DIMENSION(ims:ime,kms:kme,jms:jme)            , INTENT(IN) :: th_phy  , &
                                                                           p_phy   , &
                                                                           pi_phy  , &
                                                                           rho_phy , &
                                                                           dz8w    , &
                                                                           p8w     , &
                                                                           t8w

      !  Time (s) since the beginning of the restart.

      REAL :: curr_secs2

      !  Is this to be a history output time?  If so, compute the diagnostics.

      LOGICAL :: diag_flag

      !  The sundry dimensions required to keep a model running smoothly:
      !     The first letter:
      !        i: refers to the nominally west east direction, the inner-most (fastest)
      !           incrementing index
      !        j: refers to the nominally south north direction, the outer-most (slowest)
      !           incrementing index
      !        k: refers to the vertical direction form bottom to top, the second dimension
      !           in all 3d arrays
      !     The second letter: 
      !        d: refers to the domain size, the geophysical extent of the entire domain,
      !           not used in dimensions or looping, used to determine when we are close to
      !           the edge of the boundary
      !        m: refers to the memory size size, all 2d and 3d arrays from the Registry
      !           (passed into here via the grid structure or the I1 variables [such as
      !           p_phy, for example]) use these values for dimensioning
      !        p: refers to the patch size, the extent over which computational loops run

      INTEGER , INTENT(IN) :: ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe

      !  Hopefully unnecessary, these are the filtered dimensions.

      INTEGER , INTENT(IN) :: imsx,imex,jmsx,jmex,kmsx,kmex,    &
                              ipsx,ipex,jpsx,jpex,kpsx,kpex,    &
                              imsy,imey,jmsy,jmey,kmsy,kmey,    &
                              ipsy,ipey,jpsy,jpey,kpsy,kpey

#ifdef DM_PARALLEL
      !=============================================================
      !  Include patch communications
      !=============================================================
#     include "HALO_EM_PHYS_W.inc"
#endif

      !=============================================================
      !  Start of executable code
      !=============================================================

      CALL wrf_debug ( 100 , '--> TOP OF AFTER ALL RK STEPS' ) 
      CALL wrf_debug ( 100 , '--> CALLING DIAGNOSTICS DRIVER' )


      CALL diagnostics_driver ( grid, config_flags,               &
                                moist, chem, tracer, scalar,         &
                                th_phy, pi_phy, p_phy, rho_phy,      & 
                                p8w, t8w, dz8w,                      &
                                curr_secs2,                          &
                                diag_flag,                           &
                                ids,  ide,  jds,  jde,  kds,  kde,   &
                                ims,  ime,  jms,  jme,  kms,  kme,   &
                                ips,  ipe,  jps,  jpe,  kps,  kpe,   &
                                imsx, imex, jmsx, jmex, kmsx, kmex,  &
                                ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                                imsy, imey, jmsy, jmey, kmsy, kmey,  &
                                ipsy, ipey, jpsy, jpey, kpsy, kpey   )


   END SUBROUTINE after_all_rk_steps

SUBROUTINE temporal_spatial_averaging( &
      grid, config_flags, called_from_start, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      ips, ipe, jps, jpe, kps, kpe)

      use mpi
      USE module_configure
      USE module_state_description, only: p_m13
      USE module_domain
      !USE module_domain, ONLY : domain, get_ijk_from_grid, get_ijk_from_subgrid &
             !    ,domain_get_current_time, domain_get_start_time &
             !    ,domain_get_sim_start_time, domain_clock_get,is_alarm_tstep
      use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
      use, intrinsic :: iso_fortran_env, only: real32

      IMPLICIT NONE

      TYPE(domain), INTENT(INOUT) :: grid

      !  We are not changing any of the namelist settings.
      TYPE(grid_config_rec_type), INTENT(IN) :: config_flags

      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             ips, ipe, jps, jpe, kps, kpe

      
      DOUBLE PRECISION, DIMENSION(9, kms:kme) :: var_avg_xy !(u,uu,v,vv,w,ww,uw,vw,mag)
      REAL :: u_vol, v_vol, w_vol, uc2_xy, vc2_xy, wc2_xy, tke, taub13, vol_z,dz , dhdx, p_drag, p_drag_full ! PZ added, vol_z=sum_{i,j,k}dz(i,j,k) 
      integer, DIMENSION(ids:ide) :: idx_shift
      logical, intent(in) :: called_from_start

      integer :: i, j, k, comm, IERR
      TYPE(WRFU_Time) :: time_avg_start_date, lcl_currtime
      TYPE(WRFU_TimeInterval) :: avg_start_to_currtime, curr_time_step
      INTEGER :: Sn, Sd, rc, time_avg_sample_old, sample_count, avg_start !in seconds
      real :: nan

      DOUBLE PRECISION :: u, v, w, uu, vv, ww, uw, vw, u_c, v_c, w_c, umag,tau13, time_avg_x_shift
      CHARACTER*30 :: msg

      var_avg_xy = 0.
      taub13 = 0.
      tke = 0.
      u_vol = 0.
      v_vol = 0.
      w_vol = 0.
      uc2_xy = 0.
      vc2_xy = 0.
      wc2_xy = 0.
      vol_z = 0.
      dz = 0.
      p_drag = 0.
      p_drag_full = 0.

      nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)


      do j = jps, jpe
      do k = kps, kpe
      do i = ips, ipe
         if (j < jde) then
            var_avg_xy(1, k) = var_avg_xy(1, k) + grid%u_2(i, k, j)/((ide - ids + 1)*(jde - jds))
            var_avg_xy(2, k) = var_avg_xy(2, k) + grid%u_2(i, k, j)*grid%u_2(i, k, j)/((ide - ids + 1)*(jde - jds))
            
         endif

         if (i < ide) then
            var_avg_xy(3, k) = var_avg_xy(3, k) + grid%v_2(i, k, j)/((ide - ids)*(jde - jds + 1))
            var_avg_xy(4, k) = var_avg_xy(4, k) + grid%v_2(i, k, j)*grid%v_2(i, k, j)/((ide - ids)*(jde - jds + 1))
         endif

         if (i < ide .and. j < jde) then
            var_avg_xy(5, k) = var_avg_xy(5, k) + grid%w_2(i, k, j)/((ide - ids)*(jde - jds))
            var_avg_xy(6, k) = var_avg_xy(6, k) + grid%w_2(i, k, j)*grid%w_2(i, k, j)/((ide - ids)*(jde - jds))
         endif

         if (i < ide .and. j < jde .and. k < kde) then
            u_c = (grid%u_2(i, k, j) + grid%u_2(i+1, k, j))/2
            v_c = (grid%v_2(i, k, j) + grid%v_2(i, k, j+1))/2
            w_c = (grid%w_2(i, k, j) + grid%w_2(i, k+1, j))/2

            var_avg_xy(7, k) = var_avg_xy(7, k) + u_c*w_c/((ide - ids)*(jde- jds))
            var_avg_xy(8, k) = var_avg_xy(8, k) + v_c*w_c/((ide - ids)*(jde- jds))
            var_avg_xy(9, k) = var_avg_xy(9, k) +sqrt(u_c*u_c+v_c*v_c)/((ide - ids)*(jde - jds))
            ! PZ added to compute tke
            dz = (grid%phb(i,k+1,j) + grid%ph_2(i,k+1,j) - (grid%phb(i,k,j) + grid%ph_2(i,k,j)) )/9.81
            vol_z = vol_z + dz
            u_vol = u_vol + u_c*dz 
            uc2_xy = uc2_xy + u_c*u_c*dz 
            vc2_xy = vc2_xy + v_c*v_c*dz 
            wc2_xy = wc2_xy + w_c*w_c*dz 
            v_vol = v_vol + v_c*dz 
            w_vol = w_vol + w_c*dz
         endif

         if (i < ide .and. j < jde .and. k == kps) then
             taub13 = taub13 + grid%nba_mij(i,kps,j,P_m13)/((ide - ids)*(jde- jds))
             ! PZ added to compute pressure drag
             dhdx = (grid%phb(i+1,2,j)+grid%ph_2(i+1,2,j) +grid%phb(i+1,1,j) +grid%ph_2(i+1,1,j) - (grid%phb(i-1,2,j)+grid%ph_2(i-1,2,j) + grid%phb(i-1,1,j) + grid%ph_2(i-1,1,j)) )/4/9.81/config_flags%dx
             p_drag = p_drag + grid%p(i,k,j)*dhdx*grid%alt(i,k,j)/((ide - ids)*(jde- jds))
             p_drag_full = p_drag_full + (grid%p(i,k,j)+grid%pb(i,k,j))*dhdx*grid%alt(i,k,j)/((ide - ids)*(jde- jds))
         end if
!get the xy-averaged bottom stress 
      enddo
      enddo
      enddo



#ifdef DM_PARALLEL
      call wrf_get_dm_communicator(comm)

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, var_avg_xy, product(shape(var_avg_xy)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, u_vol, product(shape(u_vol)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, v_vol, product(shape(v_vol)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, w_vol, product(shape(w_vol)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, uc2_xy, product(shape(uc2_xy)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, vc2_xy, product(shape(vc2_xy)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, wc2_xy, product(shape(wc2_xy)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, vol_z, product(shape(vol_z)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      !CALL MPI_ALLREDUCE(MPI_IN_PLACE, tke, product(shape(tke)) &
      !&                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, taub13, product(shape(taub13)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, p_drag, product(shape(p_drag)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, p_drag_full, product(shape(p_drag_full)) &
      &                  , MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
#endif

      !PZ added
      ! compute tke_xyavg using var_avg_xy
      tke = 0.5*(uc2_xy/vol_z - (u_vol/vol_z)**2 + vc2_xy/vol_z - (v_vol/vol_z)**2 + wc2_xy/vol_z - (w_vol/vol_z)**2 )      

      grid%u_avg_xy(kps:kpe) = var_avg_xy(1, kps:kpe)
      grid%uu_avg_xy(kps:kpe) = var_avg_xy(2, kps:kpe)
      grid%v_avg_xy(kps:kpe) = var_avg_xy(3, kps:kpe)
      grid%vv_avg_xy(kps:kpe) = var_avg_xy(4, kps:kpe)
      grid%w_avg_xy(kps:kpe) = var_avg_xy(5, kps:kpe)
      grid%ww_avg_xy(kps:kpe) = var_avg_xy(6, kps:kpe)
      grid%uw_avg_xy(kps:kpe) = var_avg_xy(7, kps:kpe)
      grid%vw_avg_xy(kps:kpe) = var_avg_xy(8, kps:kpe)
      grid%umag_avg_xy(kps:kpe) = var_avg_xy(9, kps:kpe)
      grid%u_vol = u_vol/vol_z
      grid%tke = tke
      grid%taub13 = taub13
      grid%p_drag = p_drag
      grid%p_drag_full = p_drag_full
      !if (ips==1 .and. jps==1) print *,"u_avg_y(1,3) ", grid%u_avg_y(1,3,1)

      CALL wrf_atotime(config_flags%time_avg_start_date, time_avg_start_date)
      CALL domain_clock_get(grid, current_time=lcl_currtime, &
                            current_timestr_frac=msg, time_step=curr_time_step)

      lcl_currtime = lcl_currtime + curr_time_step
      avg_start_to_currtime = lcl_currtime - time_avg_start_date
      ! print *, " lcl_currtime ", msg, "",
      ! trim(config_flags%time_avg_start_date)
      ! print *, " time_diff ", avg_start_to_currtime%basetime%S, &
      !       config_flags%time_avg_interval, avg_start_to_currtime%basetime%Sn

      !Because we only enter this subroutine at whole second time step, so there
      !is no need to check fractional second
      !Does TIME_AVG_START_DATE arrive?
      if (avg_start_to_currtime%basetime%S <= 0) then
         !if TIME_AVG_START_DATE doesn't arrive, set avg value to NAN.
         grid%time_avg_sample = 0
         grid%u_avg_t_dp = nan
         grid%v_avg_t_dp = nan
         grid%w_avg_t_dp = nan
         grid%uu_avg_t_dp = nan
         grid%vv_avg_t_dp = nan
         grid%ww_avg_t_dp = nan
         grid%uw_avg_t_dp = nan
         grid%vw_avg_t_dp = nan
         grid%umag_avg_t_dp = nan
      elseif (called_from_start) then
         !if TIME_AVG_START_DATE arrives,
         !but CALLED_FROM_START, read the avg value from last run
         grid%u_avg_t_dp = grid%u_avg_t
         grid%v_avg_t_dp = grid%v_avg_t
         grid%w_avg_t_dp = grid%w_avg_t
         grid%uu_avg_t_dp = grid%uu_avg_t
         grid%vv_avg_t_dp = grid%vv_avg_t
         grid%ww_avg_t_dp = grid%ww_avg_t
         grid%uw_avg_t_dp = grid%uw_avg_t
         grid%vw_avg_t_dp = grid%vw_avg_t
         grid%umag_avg_t_dp = grid%umag_avg_t
      else !if (avg_start_to_currtime%basetime%Sn==0) then
        !only use data at whole senconds to compute the average 
        do k = kps, kpe
            if (config_flags%avg_in_xyt) then
                u = var_avg_xy(1, k)
                uu = var_avg_xy(2, k)
                v = var_avg_xy(3, k)
                vv = var_avg_xy(4, k)
                w = var_avg_xy(5, k)
                ww = var_avg_xy(6, k)
                uw = var_avg_xy(7, k)
                vw = var_avg_xy(8, k)
                umag = var_avg_xy(9, k)
                
                if (grid%time_avg_sample == 0) then
                  grid%u_avg_t_dp(k) = 0.
                  grid%v_avg_t_dp(k) = 0.
                  grid%w_avg_t_dp(k) = 0.
                  grid%uu_avg_t_dp(k) = 0.
                  grid%vv_avg_t_dp(k) = 0.
                  grid%ww_avg_t_dp(k) = 0.
                  grid%uw_avg_t_dp(k) = 0.
                  grid%vw_avg_t_dp(k) = 0.
                  grid%umag_avg_t_dp(k) = 0.
                end if

                grid%u_avg_t_dp(k) = grid%u_avg_t_dp(k) + u/config_flags%time_avg_interval
                grid%v_avg_t_dp(k) = grid%v_avg_t_dp(k) + v/config_flags%time_avg_interval
                grid%w_avg_t_dp(k) = grid%w_avg_t_dp(k) +w/config_flags%time_avg_interval
                grid%uu_avg_t_dp(k) = grid%uu_avg_t_dp(k) +uu/config_flags%time_avg_interval
                grid%vv_avg_t_dp(k) = grid%vv_avg_t_dp(k) +vv/config_flags%time_avg_interval
                grid%ww_avg_t_dp(k) = grid%ww_avg_t_dp(k) +ww/config_flags%time_avg_interval
                grid%uw_avg_t_dp(k) = grid%uw_avg_t_dp(k) +uw/config_flags%time_avg_interval
                grid%vw_avg_t_dp(k) = grid%vw_avg_t_dp(k) +vw/config_flags%time_avg_interval
                grid%umag_avg_t_dp(k) = grid%umag_avg_t_dp(k) +umag/config_flags%time_avg_interval
        
            end if ! if avg_in_xyt is true
        end do

        grid%time_avg_sample = mod(avg_start_to_currtime%basetime%S, config_flags%time_avg_interval)
      endif ! end if avg_start_to_currtime%basetime%S <= 0

!------Comment out WZ's original version, which computes running averages using
!only data at the average interval---------------------------------------------
!      elseif (config_flags%time_avg_interval <= 0 .or. &
!              mod(avg_start_to_currtime%basetime%S, config_flags%time_avg_interval) == 0) then
!         !if TIME_AVG_START_DATE arrives, and NOT CALLED_FROM_START
!         !is AVG_START_TO_CURRTIME multiple of TIME_AVG_INTEVAL?
!         time_avg_sample_old = grid%time_avg_sample
!         grid%time_avg_sample = grid%time_avg_sample + 1
!
!         !print *, " time_avg_sample, time_avg_x_shift", grid%time_avg_sample
!
!         if (time_avg_sample_old == 0) then
!            grid%u_avg_t_dp = 0.
!            grid%v_avg_t_dp = 0.
!            grid%w_avg_t_dp = 0.
!            grid%uu_avg_t_dp = 0.
!            grid%vv_avg_t_dp = 0.
!            grid%ww_avg_t_dp = 0.
!            grid%uw_avg_t_dp = 0.
!            grid%vw_avg_t_dp = 0.
!            grid%umag_avg_t_dp = 0.
!         endif
!
!         do k = kps, kpe
!            if (config_flags%avg_in_xyt) then
!                u = var_avg_xy(1, k)
!                uu = var_avg_xy(2, k)
!                v = var_avg_xy(3, k)
!                vv = var_avg_xy(4, k)
!                w = var_avg_xy(5, k)
!                ww = var_avg_xy(6, k)
!                uw = var_avg_xy(7, k)
!                vw = var_avg_xy(8, k)
!                umag = var_avg_xy(9, k)
!!            else
!!               if (j < jde) then
!!                  u = dble(grid%u_2(i, k, j))
!!                  uu = dble(grid%u_2(i, k, j))*dble(grid%u_2(i, k, j))
!!               endif
!!               if (i < ide) then
!!                  v = dble(grid%v_2(i, k, j))
!!                  vv = dble(grid%v_2(i, k, j))*dble(grid%v_2(i, k, j))
!!               endif
!!               if (i < ide .and. j < jde) then
!!                  w = dble(grid%w_2(i, k, j))
!!                  ww = dble(grid%w_2(i, k, j))*dble(grid%w_2(i, k, j))
!!               endif
!!               if (i < ide .and. j < jde .and. k < kde) then
!!                  u_c = (dble(grid%u_2(i, k, j)) + dble(grid%u_2(i+1, k,   j
!!)))/2
!!                  v_c = (dble(grid%v_2(i, k, j)) + dble(grid%v_2(i,   k,
!!j+1)))/2
!!                  w_c = (dble(grid%w_2(i, k, j)) + dble(grid%w_2(i,   k+1, j
!!)))/2
!!                  uv = u_c*v_c
!!                  uw = u_c*w_c
!!                  vw = v_c*w_c
!!                  umag = sqrt(u_c*u_c+v_c*v_c)
!!               endif
!                grid%u_avg_t_dp(k) = (grid%u_avg_t_dp(k)*time_avg_sample_old + u)/grid%time_avg_sample
!                grid%uu_avg_t_dp(k) = (grid%uu_avg_t_dp(k)*time_avg_sample_old + uu)/grid%time_avg_sample
!                grid%v_avg_t_dp(k) = (grid%v_avg_t_dp(k)*time_avg_sample_old + v)/grid%time_avg_sample
!                grid%vv_avg_t_dp(k) = (grid%vv_avg_t_dp(k)*time_avg_sample_old + vv)/grid%time_avg_sample
!                grid%w_avg_t_dp(k) = (grid%w_avg_t_dp(k)*time_avg_sample_old + w)/grid%time_avg_sample
!                grid%ww_avg_t_dp(k) = (grid%ww_avg_t_dp(k)*time_avg_sample_old + ww)/grid%time_avg_sample
!                grid%uw_avg_t_dp(k) = (grid%uw_avg_t_dp(k)*time_avg_sample_old + uw)/grid%time_avg_sample
!                grid%vw_avg_t_dp(k) = (grid%vw_avg_t_dp(k)*time_avg_sample_old + vw)/grid%time_avg_sample
!                grid%umag_avg_t_dp(k) = (grid%umag_avg_t_dp(k)*time_avg_sample_old + umag)/grid%time_avg_sample
!            endif
!
!         enddo
!      endif
!------End of WZ's original version
      grid%u_avg_t = grid%u_avg_t_dp
      grid%v_avg_t = grid%v_avg_t_dp
      grid%w_avg_t = grid%w_avg_t_dp
      grid%uu_avg_t = grid%uu_avg_t_dp
      grid%vv_avg_t = grid%vv_avg_t_dp
      grid%ww_avg_t = grid%ww_avg_t_dp
      grid%uw_avg_t = grid%uw_avg_t_dp
      grid%vw_avg_t = grid%vw_avg_t_dp
      grid%umag_avg_t = grid%umag_avg_t_dp

      return

   END SUBROUTINE temporal_spatial_averaging

END MODULE module_after_all_rk_steps
