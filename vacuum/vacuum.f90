!> Common parameters and variables for the JOREK free boundary extension.
!!
!! @see vacuum_response, vacuum_equilibrium
module vacuum
  use mod_parameters, only: var_zj, var_psi
  use phys_module, only: rst_hdf5_version, freeb_change_indices
  
  implicit none
  
  !> @name General parameters
  logical, parameter  :: vacuum_debug          = .false. !< Enable additional output and tests
  logical, parameter  :: vacuum_decouple_modes = .false. !< Option to switch off 3D wall mode coupling
  integer             :: n_dof_bnd                       !< Total number of boundary dofs per harmonic
  integer             :: n_dof_starwall                  !< Total number of boundary dofs in STARWALL response
  logical             :: CARIDDI_mode                    !< CARIDDI or STARWALL
  logical             :: halo_freeb                      !< Switch on halo current coupling with walls.
  logical             :: halo_imp                        !< Establish implicit time stepping for halos
  logical             :: vacuum_min                      !< Mode to minimalize memory consumption
  
  !> @name Resistive wall only
  real*8              :: wall_resistivity_fact           !< Scaling factor for the wall and coil resistivities specified in STARWALL
  real*8              :: wall_resistivity                !< Resistivity of the external wall
  logical             :: wall_curr_initialized = .false. !< Have the wall currents been initialized?
  integer             :: n_wall_curr                     !< Number of wall current potentials.
  integer             :: n_wall_curr_i                   !< Halo coupling: number of internal current potentials (null space)
  integer             :: n_wall_curr_e                   !< Halo coupling: number of external current potentials (range space)
  integer             :: n_shell_curr_e
  real*8, allocatable :: wall_curr(:)                    !< Wall current potentials (\f$Y_k\f$).
  real*8, allocatable :: dwall_curr(:)                   !< Change of wall current potentials (\f$\delta Y_k\f$).
  real*8, allocatable :: wall_curr_i(:)                  !< Internal wall current potentials (\f$Yi_k\f$).
  real*8, allocatable :: dwall_curr_i(:)                 !< Change of internal wall current potentials (\f$\delta Yi_k\f$).
  real*8, allocatable :: wall_curr_e(:)                  !< External wall current potentials (\f$Yi_k\f$).
  real*8, allocatable :: dwall_curr_e(:)                 !< Change of external wall current potentials (\f$\delta Yi_k\f$).
  real*8, allocatable :: old_dpsibnd_vec(:)              !< Previous delta Psi values required for time-stepping with zeta/=0
  real*8, allocatable :: old_dIelec_w(:)                 !< Previous dIelec_w
  real*8, allocatable :: old_dIelec_J(:)                 !< Previous dIelec_J
  real*8, allocatable :: dIelec_w(:)                     !< Delta Wall electrode current for halo coupling
  real*8, allocatable :: dIelec_J(:)                     !< Delta Shell electrode current for halo coupling
  real*8, allocatable :: Ielec_w(:)                     !< Wall electrode current for halo coupling
  real*8, allocatable :: Ielec_J(:)                     !< Shell electrode current for halo coupling
  real*8, allocatable :: Ielec_w_old(:)                 !< Wall electrode current for halo coupling, previous step
  real*8, allocatable :: Ielec_J_old(:)                 !< Shell electrode current for halo coupling, previous step
  ! Matrices for offline solution of the wall current
  real*8, allocatable  :: psibnd_mat(:,:) , Ielecw_mat(:,:), IelecJ_mat(:,:), Ielecw_phys_mat_past(:,:), IelecJ_phys_mat_past(:,:),  Ielecw_phys_mat_past_full(:,:), IelecJ_phys_mat_past_full(:,:), wall_curr_i_mat(:,:), coil_curr_i_mat(:,:)
  real*8, allocatable  :: Ielecw_phys_mat_pres(:,:), IelecJ_phys_mat_pres(:,:),  Ielecw_phys_mat_pres_full(:,:), IelecJ_phys_mat_pres_full(:,:)

  !> Data type for response matrices distributed over the MPI tasks
  type :: t_distrib_mat
    real*8, allocatable :: loc_mat(:,:)                  !< Local chunk of the matrix
    logical             :: distrib                       !< Is the matrix distributed?
    logical             :: row_wise                      !< Is the matrix distributed rowwise (otherwise columnwise)?
    integer             :: ind_start                     !< Minimum row/column index of local chunk.
    integer             :: ind_end                       !< Maximum row/column index of local chunk.
    integer             :: step                          !< "chunk size" of each MPI task. Note that the last MPI task has a larger chunk size than this!
    integer             :: dim(2)                        !< global matrix dimensions
  end type t_distrib_mat

  !> @name JOREK vacuum response matrices
  !! Response matrices derived from STARWALL response (w=wall, p=plasma)
  real*8, allocatable  :: response_d_b(:)                !< \f$\hat{B}\f$ in the documentation
  real*8, ALLOCATABLE  :: response_d_b_i(:)
  real*8, allocatable  :: response_d_c(:)                !< \f$\hat{C}\f$ in the documentation
  real*8, allocatable  :: response_d_c_i(:)               
  real*8, allocatable  :: response_m_e(:,:)              !< \f$\hat{E}\f$ in the documentation
  real*8, allocatable  :: response_m_e_i(:,:)            !< \f$\hat{E_i}\f$ in the documentation
  real*8, allocatable  :: response_m_h(:,:)              !< \f$\hat{H}\f$ in the documentation
  real*8, allocatable  :: response_m_j(:,:)              !< \f$\hat{J}\f$ in the documentation
  real*8, allocatable  :: response_m_j_i(:,:)            !< \f$\hat{J_i}\f$ in the documentation
  real*8, allocatable  :: response_m_eh(:,:)             !< \f$\hat{E_H}\f$ in the documentation
  real*8, allocatable  :: response_m_zh(:,:)              !< \f$\hat{Z_H}\f$ in the documentation
  real*8, allocatable  :: response_m_ph(:,:)             !< \f$\hat{P_H}\f$ in the documentation
  real*8, allocatable  :: response_m_qh(:,:)              !< \f$\hat{Q_H}\f$ in the documentation
  real*8, allocatable  :: response_m_xh(:,:)             !< \f$\hat{X_H}\f$ in the documentation
  real*8, allocatable  :: response_m_vh(:,:)              !< \f$\hat{V_H}\f$ in the documentation
  real*8, allocatable  :: response_m_eq(:,:)             !< Response matrix for vacuum_equil
  ! Static pre-jacobian matrices in the halo time stepping
  real*8,allocatable :: response_m_jewj_i(:,:), response_m_wtotww_i(:,:), response_m_t_h(:,:)

  !> @name Equilibrium coil contributions
  integer             :: n_coils                         !< number of poloidal field coils in coil_field.dat
  integer             :: n_pf_coils                      !< number of poloidal field coils
  logical             :: starwall_equil_coils            !< specify wheter the equilibrium PF coils will be given by STARWALL or not
  logical             :: find_pf_coil_currents           !< search for optimal pf_coil currents to build a free-bnd equil? [[jorek-starwall-faqs|fbnd_eq_FAQs]]
  real*8, allocatable :: I_coils(:)                      !< coil currents
  real*8, allocatable :: Y_coils0(:)                     !< imposed STARWALL coil currents source
  real*8, ALLOCATABLE :: Yi_coils0(:)                    !< imposed STARWALL coil currents in null space
  real*8              :: vertical_FB                     !< a variable for the feedback control of the plasma's vertical position
  real*8              :: radial_FB                       !< a variable for the feedback control of the plasma's radial position   (during equilibrium)
  real*8, allocatable :: bext_tan(:,:)                   !< external tangential field
  real*8, allocatable :: bext_nor(:,:)                   !< external normal field
  real*8, allocatable :: bext_psi(:,:)                   !< external poloidal flux      
  real*8              :: psi_offset_freeb                !< Allows to shift the value of psi by a global constant for freeb_equil (improves convergence)
  
  !> @name Equilibrium parameters for feedback
  real*8              :: current_ref                     !< Target total plasma current Ip for the feedback (FB) [[jorek-starwall-faqs|fbnd_eq_FAQs]]
  real*8              :: FB_Ip_position                  !< Amplification factor for Ip feedback (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  real*8              :: FB_Ip_integral                  !< Amplification factor for Ip feedback (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  real*8              :: Z_axis_ref                      !< Target magnetic axis vertical position (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  real*8              :: R_axis_ref                      !< Optional target magnetic axis radial position (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  real*8              :: FB_Zaxis_position               !< Amplification factor for Zaxis feedback (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  real*8              :: FB_Zaxis_derivative             !< Amplification factor for Zaxis feedback (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  real*8              :: FB_Zaxis_integral               !< Amplification factor for Zaxis feedback (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  real*8              :: cte_current_FB_fact             !< Constant factor that scales FF'& T profiles before freebnd GS iterations (switches off current FB) 
  integer             :: start_VFB                       !< Iteration for starting vertical feedback (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  integer             :: n_feedback_current              !< Feedback will be performed each n_... iterations (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  integer             :: n_feedback_vertical             !< Feedback will be performed each n_... iterations (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])
  integer             :: n_iter_freeb                    !< Number of iterations for freeboundary equilibirum (see [[jorek-starwall-faqs|fbnd_eq_FAQs]])

  
  ! ### various variables, some need to be removed
  real*8, allocatable :: R_coils(:), Z_coils(:)          ! ### old
  real*8, allocatable :: dR_coils(:), dZ_coils(:)        ! ### old
  real*8, allocatable :: coil_voltages(:)                !< Coil voltages (not ready yet)
  real*8              :: current_FB_fact  = 1.d0         !< Factor used for current feedback during the freeboundary equilibrium
  real*8, allocatable :: diag_coil_curr(:,:), pf_coil_curr(:,:), rmp_coil_curr(:,:)
  integer,parameter   :: COIL_NAME_LEN=64                !< Max name length per coil (same as in starwall)
  real*8, allocatable :: net_tor_wall_curr(:)            !< Time trace of net toroidal wall current (live data)
  character(len=COIL_NAME_LEN), allocatable:: diag_coil_name(:), pf_coil_name(:), rmp_coil_name(:)  !< Names for coils needed for live_data
  
  type :: t_starwall_response
    integer :: file_version           = 9999
    logical :: initialized            = .false.
    integer :: n_bnd                  = -1
    integer :: nd_bez                 = -1
    integer :: n_pts_real_w           = -1
    integer :: n_pts_real_J           = -1
    integer :: n_electrode_w           = -1
    integer :: n_electrode_J           = -1
    integer :: ncoil                  = -1
    integer :: npot_w                 = -1
    integer :: n_w                    = -1
    integer :: n_nw                   = -1
    integer :: ntri_w                 = -1
    integer :: n_tor                  = -1
    integer :: n_tor0                 = -1

    ! --- Additional STARWALL input parameters (version >= 5)
    integer :: nv                     = -1  !< Number of toroidal points of the control surface
    integer :: n_points               = -1  !< Number of triangles per JOREK boundary element (control surface)
    integer :: iwall                  = -1  !< 1 if the wall is represented with Fourier harmonics 
    integer :: nwu                    = -1  !< Number of poloidal grid points of the wall
    integer :: nwv                    = -1  !< Number of toroidal grid points of the wall
    integer :: mn_w                   = -1  !< Number of Fourier harmonics to represent the wall contour
    integer :: MAX_MN_W               = -1  !< MAX number of mn_w
    integer, allocatable :: m_w(:)          !< Wall contour: poloidal wall harmonics
    integer, allocatable :: n_w_fourier(:)  !< Wall contour: toroidal wall harmonics
    real*8,  allocatable :: rc_w(:)         !< Wall contour: Fourier cosine coefficients for R
    real*8,  allocatable :: rs_w(:)         !< Wall contour: Fourier sine   coefficients for R
    real*8,  allocatable :: zc_w(:)         !< Wall contour: Fourier cosine coefficients for Z
    real*8,  allocatable :: zs_w(:)         !< Wall contour: Fourier sine   coefficients for Z

    integer :: ntri_c                 = 0  !< Number of coil triangles
    integer :: n_pol_coils            = 0  !< Number of poloidal field coils
    integer :: n_rmp_coils            = 0  !< Number of RMP coils
    integer :: n_voltage_coils        = 0  !< Number of "voltage coils"
    integer :: n_diag_coils           = 0  !< Number of diagnostic coils
    integer :: ind_start_coils        = -1 !< starting index of coils in wall matrix
    integer :: ind_start_pol_coils    = -1 !< Index of first poloildal field coil
    integer :: ind_start_rmp_coils    = -1 !< Index of first RMP coil
    integer :: ind_start_voltage_coils= -1 !< Index of first voltage coil
    integer :: ind_start_diag_coils   = -1 !< Index of first diagnostic coil
    integer, allocatable :: jtri_c(:)      !< Number of triangles per coil
    real*8,  allocatable :: x_coil(:,:)      !< x-position of coil triangle nodes
    real*8,  allocatable :: y_coil(:,:)      !< y-position of coil triangle nodes
    real*8,  allocatable :: z_coil(:,:)      !< z-position of coil triangle nodes
    real*8,  allocatable :: phi_coil(:,:)    !< "Potential" at coil triangle nodes
    real*8,  allocatable :: eta_thin_coil(:) !< Thin wall resistivity of coil triangles
    real*8,  allocatable :: coil_resist(:)   !< Resistance of each coil
    character(len=COIL_NAME_LEN),  allocatable :: coil_name(:)   !< Name of each coil
    real*8  :: eta_thin_w             = 1. !< Thin wall resistivity of wall triangles
    integer, allocatable :: i_tor(:)
    real*8,  allocatable :: d_yy(:)
    real*8,  allocatable :: dii_yy(:)
    type(t_distrib_mat)  :: a_ye
    type(t_distrib_mat)  :: a_ey
    type(t_distrib_mat)  :: a_ee
    type(t_distrib_mat)  :: a_id
    type(t_distrib_mat)  :: s_ww
    type(t_distrib_mat)  :: s_ww_inv
    type(t_distrib_mat)  :: s_ww_ii_inv
   
    !!!! Matrices for halo coupling  !!!!
    type(t_distrib_mat)  :: aeq_el
    type(t_distrib_mat)  :: a_yie
    type(t_distrib_mat)  :: d_yieqe
    type(t_distrib_mat)  :: e_yiwe
    type(t_distrib_mat)  :: g_yiwe
    type(t_distrib_mat)  :: a_eyi
    type(t_distrib_mat)  :: a_eye
    type(t_distrib_mat)  :: fwtfw_inv
    type(t_distrib_mat)  :: w_h  ! Matrix M_eye * inv(FwtFw) for time-stepping
    real*8, allocatable  :: ielsphi_real_w(:,:)
    real*8, allocatable  :: ielsphi_real_J(:,:)
    ! Electrode mapping Matrices in real space
    real*8, allocatable  :: Jmat_real_w(:,:)
    real*8, allocatable  :: Jmat_real_J(:,:)
    ! Electrode mapping Matrices in Fourier space
    real*8, allocatable  :: Jmat_w(:,:) 
    real*8, allocatable  :: Jmat_J(:,:)
    ! type(t_distrib_mat)  :: Jmat_w 
    ! type(t_distrib_mat)  :: Jmat_J

    real*8,  allocatable :: xyzpot_w(:,:)
    real*8,  allocatable :: phi0_w(:,:)   !< Toroidal net wall current potentials
    integer, allocatable :: jpot_w(:,:)
  end type t_starwall_response
  
  type(t_starwall_response) :: sr              !< STARWALL response
  
  !> Coil current specification.
  type :: t_coil_curr_input
    real*8             :: current   = 0.d0  !< Current of the coil in Ampere*Turns
    real*8             :: pert      = 0.d0  !< Pert. of coil current in Ampere*Turns to speed-up VDE.
    real*8             :: pert_start_time  = 1.d33   !< Starting time of pert. of coil current in JOREK_time.
    real*8             :: pert_growth_time = 1.d-12  !< Ramp-up time of pert. of coil current in JOREK_time.
    character(len=256) :: curr_file = 'none'!< Ascii file with coil current time trace.
    real*8             :: time_shift    = 0.d0  !< Shift time of time trace.
    real*8             :: time_scale    = 1.d0  !< Scale time of time trace.
    real*8             :: curr_scale    = 1.d0  !< Scale amplitude of time trace.
    character(len=512) :: curr_expr = 'none'!< Analyt. expression for time trace (Python).
    real*8             :: max_time  = 0.d0  !< Evaluate analytical expression up to this time.
    integer            :: len       = 0     !< Evaluate analytical expression on this number of time points.
  end type t_coil_curr_input
  type :: t_coil_curr_time_trace
    integer            :: len      !< Number of points in numerical time trace.
    real*8, allocatable:: time(:)  !< time-values of numerical time trace
    real*8, allocatable:: curr(:)  !< current-values of numerical time trace
  end type t_coil_curr_time_trace
  integer, parameter              :: MAX_COILS = 299
  type(t_coil_curr_input), target :: diag_coils(MAX_COILS)    ! see [[jorek-starwall-faqs|jorek_starwall_FAQs]]
  type(t_coil_curr_input), target :: rmp_coils(MAX_COILS)     ! see [[jorek-starwall-faqs|jorek_starwall_FAQs]]
  type(t_coil_curr_input), target :: voltage_coils(MAX_COILS) ! not ready yet (see [[jorek-starwall-faqs|jorek_starwall_FAQs]])
  type(t_coil_curr_input), target :: pf_coils(MAX_COILS)      ! see [[jorek-starwall-faqs|jorek_starwall_FAQs]]
  type(t_coil_curr_time_trace)    :: coil_curr_time_trace(4*MAX_COILS)
  real*8 :: vert_FB_amp(MAX_COILS) = 0.d0 !< Tune direction and magnitude of vert feedback for each poloidal field coil ([[jorek-starwall-faqs|eq_FAQs]])
  real*8 :: rad_FB_amp(MAX_COILS) = 0.d0  !< Tune direction and magnitude of vert feedback for each poloidal field coil ([[jorek-starwall-faqs|eq_FAQs]])
  
  ! --- Parameters for the feedback on the vertical position during timestepping (VFB), see ([[active_controller_model_for_vertical_stabilization|documentation]])
  character(len=256)  :: vert_pos_file = 'none'
  !> Time trace of axis position to match
  type :: t_Z_axis_ref_ts     
    integer                :: len = 0      !< Number of points in numerical time trace
    real*8, allocatable    :: time(:)      !< time-values of numerical time trace
    real*8, allocatable    :: position(:)  !< evolution of vertical axis position over time
  end type t_Z_axis_ref_ts
  real*8                        :: start_VFB_ts                  !< start time of active VFB during simulation ([JOREK units])
  real*8                        :: vert_FB_amp_ts(MAX_COILS)     !< Amplitude and sign of vert feedback for each coil ([[jorek-starwall-faqs|eq_FAQs]])
  real*8                        :: I_coils_max(MAX_COILS)        !< Current limit of each coil ([Ampere])
  real*8                        :: vert_FB_gain(3)               !< Gain parameters for vertical feedback controller
  real*8                        :: vert_FB_tact                  !< Time interval between two controller actions ([JOREK units])
  real*8                        :: dZ_axis_integral              !< Integrated values of Z_axis-Z_reference for controller
  real*8, allocatable           :: vert_FB_response(:,:)         !< Controller response (PID gain * err) and target axis
  type(t_Z_axis_ref_ts), target :: Z_axis_ref_ts                 !< Time trace of axis target position
  
  
  
  contains
  
  
  
  !> Set coil_curr_time_trace based on user input.
  subroutine set_coil_curr_time_trace()
    use profiles, only: readProf
    
    integer :: i, j, k, l, i_start_coil, i_end_coil, err
    character(len=60) :: s, filename
    real*8 :: r
    class(t_coil_curr_input), pointer :: coil_curr_input
    
    !--- Make sure that the main coil current vector (I_coils) is properly allocated
    if (sr%ncoil > 0) then
      n_coils = sr%ncoil
      if (.not. allocated(I_coils)) then 
        allocate(I_coils(n_coils))
        I_coils(:)                =  0.d0
        i_start_coil = sr%ind_start_pol_coils
        i_end_coil   = i_start_coil + sr%n_pol_coils - 1
        I_coils(i_start_coil:i_end_coil) =  pf_coils(1:sr%n_pol_coils)%current
        
        i_start_coil = sr%ind_start_rmp_coils
        i_end_coil   = i_start_coil + sr%n_rmp_coils - 1
        I_coils(i_start_coil:i_end_coil) =  rmp_coils(1:sr%n_rmp_coils)%current 
        n_coils                   =  sr%ncoil
          
      endif
    endif
        
    do i = 1, sr%ncoil
            
      ! --- Which coil are we processing?
      if ( (i>=sr%ind_start_pol_coils) .and. (i<sr%ind_start_pol_coils+sr%n_pol_coils) ) then
        j = i - sr%ind_start_pol_coils + 1
        coil_curr_input => pf_coils(j)
      else if ( (i>=sr%ind_start_rmp_coils) .and. (i<sr%ind_start_rmp_coils+sr%n_rmp_coils) ) then
        j = i - sr%ind_start_rmp_coils + 1
        coil_curr_input => rmp_coils(j)
      else if ( (i>=sr%ind_start_diag_coils) .and. (i<sr%ind_start_diag_coils+sr%n_diag_coils) ) then
        j = i - sr%ind_start_diag_coils + 1
        coil_curr_input => diag_coils(j)
      else if ( (i>=sr%ind_start_voltage_coils) .and. (i<sr%ind_start_voltage_coils+sr%n_voltage_coils) ) then
        j = i - sr%ind_start_voltage_coils + 1
        coil_curr_input => voltage_coils(j)
      end if
      
      ! --- Set up numerical coil current time trace based on the input type...
      if ( coil_curr_input%curr_file /= 'none' ) then ! ... numerical input by ascii file
        
        call readProf(coil_curr_time_trace(i)%time, coil_curr_time_trace(i)%curr, &
          coil_curr_time_trace(i)%len, coil_curr_input%curr_file)

      else if ( coil_curr_input%curr_expr /= 'none' ) then ! ... analytical Python expression
        
        ! --- Python script
        call random_seed()
        err = 1
        do while ( err /= 0 )
          call random_number(r)
          l = r * 99999999
          write(s,*) l
          filename='./jorek_curr_expr_'//trim(adjustl(s))//'.py'
          open(42, file=trim(filename), status='new', iostat=err)
        end do
        111 format(2a)
        112 format(a,i16)
        113 format(a,es25.16)
        write(42,111) 'from math import *'
        write(42,111) 'def f(t):'
        write(42,111) '  return ', trim(coil_curr_input%curr_expr)
        write(42,112) 'len=', coil_curr_input%len
        write(42,113) 'tmin=', - coil_curr_input%time_shift / coil_curr_input%time_scale - 1.d-12
        write(42,113) 'tmax=', ( coil_curr_input%max_time - coil_curr_input%time_shift ) / coil_curr_input%time_scale + 1.d-12
        write(42,111) 'for x in range(1,len):'
        write(42,111) '  t=tmin+(x-1)/float(len-1)*(tmax-tmin)'
        write(42,111) '  s = "%25.16e"%t'
        write(42,111) '  s += "%25.16e"%f(t)'
        write(42,111) '  print(s)'
        close(42)
        
        ! --- Call Python
        call system('python ./jorek_curr_expr_'//trim(adjustl(s))//'.py > ./jorek_curr_expr_'//trim(adjustl(s))//'.dat')
        
        ! --- Read the result
        call readProf(coil_curr_time_trace(i)%time, coil_curr_time_trace(i)%curr, &
          coil_curr_time_trace(i)%len, './jorek_curr_expr_'//trim(adjustl(s))//'.dat')

        ! --- Delete temporary files
        call system('rm ./jorek_curr_expr_'//trim(adjustl(s))//'.py ./jorek_curr_expr_'//trim(adjustl(s))//'.dat')
        
      else ! ... just a value plus a perturbation
        
        if (allocated(coil_curr_time_trace(i)%time)) deallocate(coil_curr_time_trace(i)%time)
        if (allocated(coil_curr_time_trace(i)%curr)) deallocate(coil_curr_time_trace(i)%curr)
        
        allocate(coil_curr_time_trace(i)%time(4) )
        allocate(coil_curr_time_trace(i)%curr(4) )
        coil_curr_time_trace(i)%len = 4
        
        coil_curr_time_trace(i)%time(1)   = -1.d50
        coil_curr_time_trace(i)%time(2)   = coil_curr_input%pert_start_time
        coil_curr_time_trace(i)%time(3)   = coil_curr_input%pert_start_time + coil_curr_input%pert_growth_time
        coil_curr_time_trace(i)%time(4)   = 1.d50
        
        coil_curr_time_trace(i)%curr(1:2) = coil_curr_input%current
        coil_curr_time_trace(i)%curr(3:4) = coil_curr_input%current + coil_curr_input%pert
        
      end if
      
      coil_curr_time_trace(i)%time = (coil_curr_time_trace(i)%time * coil_curr_input%time_scale) &
        + coil_curr_input%time_shift
      coil_curr_time_trace(i)%curr = coil_curr_time_trace(i)%curr * coil_curr_input%curr_scale
      
      if ( coil_curr_time_trace(i)%time(1) > 0.d0 ) then
        write(*,*) 'ERROR: A coil current time trace does not start at time 0. Maybe you have time_shift wrong?', i
        stop
      end if
      
    end do
    
  end subroutine set_coil_curr_time_trace
  
  
  
  !> Check whether the user has supplied the right coil current time trace input.
  subroutine check_coil_curr_time_trace_input(coils_number)
  
    integer, intent(in) :: coils_number
    
    integer :: i, j
    logical :: changed_by_user
    class(t_coil_curr_input), pointer :: coil_curr_input
    
    do i = 1, 4
      do j = 1, MAX_COILS
        
        if ( i == 1 ) then
          coil_curr_input => pf_coils(j)
        else if ( i == 2 ) then
          coil_curr_input => rmp_coils(j)
        else if ( i == 3 ) then
          coil_curr_input => diag_coils(j)
        else if ( i == 4 ) then
          coil_curr_input => voltage_coils(j)
        end if
        
        changed_by_user = ( (coil_curr_input%current    /= 0.d0  ) .or. (coil_curr_input%pert       /= 0.d0)  &
                       .or. (coil_curr_input%curr_file  /= 'none') .or. (coil_curr_input%time_shift /= 0.d0)  &
                       .or. (coil_curr_input%time_scale /= 1.d0  ) .or. (coil_curr_input%curr_scale /= 1.d0)  &
                       .or. (coil_curr_input%curr_expr  /= 'none') .or. (coil_curr_input%max_time   /= 0.d0)  &
                       .or. (coil_curr_input%len        /= 0     ) .or. (vert_FB_amp(j)             /= 0.d0)  )

        if ( (j > coils_number) .and. changed_by_user ) then
          write(*,*) 'ERROR: Coil current input has been provided for a coil not existing.', i, j
          stop
        end if
        
        if ( (coil_curr_input%curr_expr /= 'none') .and. ( (coil_curr_input%len <= 0) .or. (coil_curr_input%max_time <= 0.d0) ) ) then
          write(*,*) 'ERROR: When prescribing a coil current via %curr_expr, you need to specify also %len and %max_time.'
          stop
        end if
        
      end do
    end do
  end subroutine check_coil_curr_time_trace_input
  
  
  
  
  !> Read the prescribed time evolution profile of Z_axis from a file, if prsent. Otherwise use
  !! either the input value Z_axis_ref or the equilibrium value
  subroutine read_Z_axis_profile()

    use profiles, only: readProf
    use equil_info, only: ES
    if (vert_pos_file /= 'none') then
      call readProf(Z_axis_ref_ts%time, Z_axis_ref_ts%position, Z_axis_ref_ts%len, vert_pos_file)
    else
      if (Z_axis_ref > 1.d10) Z_axis_ref = ES%Z_axis

      if (allocated(Z_axis_ref_ts%time))     deallocate(Z_axis_ref_ts%time)
      if (allocated(Z_axis_ref_ts%position)) deallocate(Z_axis_ref_ts%position)
      allocate(Z_axis_ref_ts%time    (2))
      allocate(Z_axis_ref_ts%position(2))
      Z_axis_ref_ts%len          =  2
      Z_axis_ref_ts%time     (1) = -1.d12
      Z_axis_ref_ts%time     (2) =  1.d12
      Z_axis_ref_ts%position (1) =  Z_axis_ref 
      Z_axis_ref_ts%position (2) =  Z_axis_ref
    endif
    call check_Z_axis_profile() 
  end subroutine read_Z_axis_profile
  
  
  
  !> Basic checks that the prescribed Z_axis profile provided makes sense
  subroutine check_Z_axis_profile()
    if (sum(abs(vert_FB_amp_ts(1:n_pf_coils)))>1.d-6) then
      if  (maxval(abs(Z_axis_ref_ts%position(:))) > 1.d10) then
        write(*,*) 'ERROR: target Z_axis beyond Machine limits'
        stop
      else if (minval(I_coils_max(1:n_coils)) .lt. 0) then
        write(*,*) 'ERROR: The maximum value of the coil cannot be smaller than 0.'
        stop
      else if (Z_axis_ref_ts%time(1)>0.d0) then        
        write(*,*) 'ERROR: The Z_axis time trace does not start at time 0. Check your input file'
        stop
      else if (Z_axis_ref_ts%len .lt. 2) then        
        write(*,*) 'ERROR: The length of the profile for the axis target position must be larger than 1'
        stop
      endif
    endif
  end subroutine check_Z_axis_profile
  
  
  
  !> Preset freeboundary related input parameters to reasonable default values.
  subroutine vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(out)   :: freeboundary_equil, freeboundary, resistive_wall
    
    ! --- Preset namelist input parameters.
    CARIDDI_mode         = .false.            !< True if CARIDDI input file
    halo_freeb           = .false.
    halo_imp             = .false.
    vacuum_min           = .false.
    freeboundary_equil   = .false.
    starwall_equil_coils = .false.
    find_pf_coil_currents= .false.
    freeboundary         = .false.
    resistive_wall       = .true.
    wall_resistivity     = 0.d0
    wall_resistivity_fact= 1.d0
        
    current_ref          = 1.d22
    cte_current_FB_fact  = -1d99
    FB_Ip_position       = 0.2d0
    FB_Ip_integral       = 0.01d0
    n_feedback_current   = 2
        
    Z_axis_ref           = 1.d22
    R_axis_ref           = -99.d0
    FB_Zaxis_position    = 1.d0
    FB_Zaxis_derivative  = 0.d0
    FB_Zaxis_integral    = 0.d0
    n_feedback_vertical  = 1
    start_VFB            = 10
    
    n_iter_freeb         = 900
    
    psi_offset_freeb     = 0.d0

  ! ---- Parameters for vertical feedback (VFB)
    start_VFB_ts          = 0.d0
    vert_FB_amp_ts        = 0.d0   ! amplification factor (of PF coil)
    vert_FB_gain(:)       = 0.d0   ! Proportional, derivative, integral gain of VFB controller
    vert_FB_tact          = 1.d-9  ! Tact of VFB controller
    I_coils_max           = 1.d99  ! Maximum absolute value for coils
    dZ_axis_integral      = 0.d0   ! Integrated error of the Z-axis
  end subroutine vacuum_preset
  
  
  
  !> Initialize vacuum, ensure consistency of input parameters.
  subroutine vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(inout) :: freeboundary_equil, freeboundary, resistive_wall

    integer  :: i
    
    ! --- Make input parameters consistent.
    freeboundary   = freeboundary .or. freeboundary_equil
    resistive_wall = freeboundary .and. resistive_wall
        
    ! --- Initialize some variables.
    sr%n_bnd  = 0
    sr%nd_bez = 0
    sr%n_pts_real_w = 0
    sr%n_pts_real_J = 0
    sr%n_electrode_w = 0
    sr%n_electrode_J = 0
    sr%ncoil  = 0
    sr%npot_w = 0
    sr%n_w    = 0
    sr%ntri_w = 0
    sr%n_tor  = 0
    sr%n_tor0 = 0

    do i=1, MAX_COILS
      if ( (my_id == 0) .and. (pf_coils(i)%pert > 0) .and. ( pf_coils(i)%pert_start_time>1.d30 ) ) then
        write(*,*) 'WARNING: Poloidal field coil perturbation pf_coils%pert has been set by the user, but will not be applied since pert_start_time was not set to a reasonable value.'
     end if
   end do
    
  end subroutine vacuum_init
  
  
  
  !> Allows to decide if free- or fixed-boundary conditions apply for a certain combination of
  !! toroidal harmonic i_tor and variable i_var. This is important, as the STARWALL response
  !! must not be provided for all toroidal modes and fixed boundary conditions remain implemented
  !! for some of the JOREK variables.
  logical function is_freebound(i_tor, i_var)
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_tor
    integer, intent(in) :: i_var
    
    ! --- Local variables
    integer :: i
    
    ! --- Free boundary conditions only if STARWALL response provided for this toroidal harmonic
    is_freebound = .false.
    do i = 1, sr%n_tor
      is_freebound = is_freebound .or. ( i_tor == sr%i_tor(i) )
    end do
    
    ! --- Free boundary conditions only for certain variables
    is_freebound = is_freebound .and. ( (i_var == var_zj) .or. (i_var == var_psi) )
    
  end function is_freebound
  
  
  
  !> Import some vacuum-related data from the restart file
  !!
  !! @todo Does not work currently if variable freeboundary is changed between export and import!
  subroutine import_restart_vacuum(file_handle, freeboundary, resistive_wall)
    
    use phys_module, only: t_start
    
    ! --- Routine parameters
    integer, intent(in) :: file_handle
    logical, intent(in) :: freeboundary
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    logical :: resistive_wall_rst
    integer :: ierr
    
    if ( freeboundary ) then
      
      read(file_handle, iostat=ierr) resistive_wall_rst
      if ( ierr /= 0 ) then
        write(*,*) 'WARNING: Restarting a simulation with freeboundary=.t. which was run with'
        write(*,*) '  freeboundary=.f. so far.'
        return
      else if ( resistive_wall .neqv. resistive_wall_rst ) then
        write(*,*) 'ERROR: It is currently not possible to restart a JOREK simulation with a'
        write(*,*) '  modified setting for resistive_wall.'
        stop
      end if
      
      if ( resistive_wall ) then

        if (halo_freeb .AND. CARIDDI_mode) then

          read(file_handle) n_wall_curr_i, n_wall_curr_e, n_shell_curr_e, n_wall_curr_e,n_dof_starwall

          if ( allocated(wall_curr_i) ) deallocate(wall_curr_i)
          allocate( wall_curr_i(n_wall_curr_i) )
          read(file_handle) wall_curr_i(:)
        
          if ( allocated(dwall_curr_i) ) deallocate(dwall_curr_i)
          allocate( dwall_curr_i(n_wall_curr_i) )
          read(file_handle) dwall_curr_i(:)

          if ( allocated(wall_curr_e) ) deallocate(wall_curr_e)
          allocate( wall_curr_e(n_wall_curr_e) )
          read(file_handle) wall_curr_e(:)
        
          if ( allocated(dwall_curr_e) ) deallocate(dwall_curr_e)
          allocate( dwall_curr_e(n_wall_curr_e) )
          read(file_handle) dwall_curr_e(:)

          if ( allocated(Ielec_J) ) deallocate(Ielec_J)
          allocate( Ielec_J(n_shell_curr_e) )
          read(file_handle) Ielec_J(:)

          if ( allocated(Ielec_w) ) deallocate(Ielec_w)
          allocate( Ielec_w(n_wall_curr_e) )
          read(file_handle) Ielec_w(:)

          if ( allocated(Ielec_J_old) ) deallocate(Ielec_J_old)
          allocate( Ielec_J_old(n_shell_curr_e) )
          read(file_handle) Ielec_J_old(:)

          if ( allocated(Ielec_w_old) ) deallocate(Ielec_w_old)
          allocate( Ielec_w_old(n_wall_curr_e) )
          read(file_handle) Ielec_w_old(:)

          if ( allocated(dIelec_J) ) deallocate(dIelec_J)
          allocate( dIelec_J(n_shell_curr_e) )
          read(file_handle) dIelec_J(:)

          if ( allocated(dIelec_w) ) deallocate(dIelec_w)
          allocate( dIelec_w(n_wall_curr_e) )
          read(file_handle) dIelec_w(:)
        else

          read(file_handle) n_wall_curr, n_dof_starwall
          
          if ( allocated(wall_curr) ) deallocate(wall_curr)
          allocate( wall_curr(n_wall_curr) )
          read(file_handle) wall_curr(:)
          
          if ( allocated(dwall_curr) ) deallocate(dwall_curr)
          allocate( dwall_curr(n_wall_curr) )
          read(file_handle) dwall_curr(:)
        
        endif
          
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
        old_dpsibnd_vec = 0.d0 !###

        if (halo_freeb) then 
          if ( allocated(old_dIelec_w) ) deallocate(old_dIelec_w)
          allocate( old_dIelec_w(n_wall_curr_e) )
          old_dIelec_w = 0.d0 
          if ( allocated(old_dIelec_J) ) deallocate(old_dIelec_J)
          allocate( old_dIelec_J(n_shell_curr_e) )
          old_dIelec_J = 0.d0 
        endif
        
        if ( vacuum_debug .and. resistive_wall ) then
          write(*,*) 'DEBUG: Checksums'
          write(*,*) 'wall_curr', sum(abs(wall_curr))
          write(*,*) 'dwall_curr', sum(abs(dwall_curr))
          write(*,*) 'END: Checksums'
        end if
        
        wall_curr_initialized = .true.
        
      end if
      
      read(file_handle) current_FB_fact
      read(file_handle) dZ_axis_integral
      read(file_handle) n_coils
      if ( n_coils /= 0 ) then
        if ( allocated(I_coils) ) deallocate(I_coils)
        allocate( I_coils(n_coils) )
        read(file_handle) I_coils(:)
      end if
      
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
    
  end subroutine import_restart_vacuum


  !> Import some vacuum-related data from the HDF5 restart file
  !!
  !! @todo Does not work currently if variable freeboundary is changed between export and import!
  subroutine import_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)
    
    use phys_module, only: t_start, nstep, index_start
#ifdef USE_HDF5
    use hdf5
    use hdf5_io_module
    
    ! --- Routine parameters
    integer(HID_T), intent(in) :: file_id
#else
    integer,        intent(in) :: file_id
#endif
    logical,        intent(in) :: freeboundary
    logical,        intent(in) :: resistive_wall
    
#ifdef USE_HDF5
    ! --- Local variables
    logical   :: freeboundary_rst, resistive_wall_rst
    integer   :: n_diag_coil = 0, n_pf_coil = 0, n_rmp_coil = 0
    character :: t_freeboundary, t_resistive_wall
    real*8, allocatable :: t_diag_coil_curr(:,:), t_pf_coil_curr(:,:), t_rmp_coil_curr(:,:)
    character(len=COIL_NAME_LEN), allocatable :: t_diag_coil_name(:), t_pf_coil_name(:), t_rmp_coil_name(:)
    real*8, allocatable :: t_net_tor_wall_curr(:)   
    real*8, allocatable :: t_vert_FB_response(:,:)
 
    call HDF5_char_reading(file_id,t_freeboundary,"freeboundary")
    freeboundary_rst = (t_freeboundary == "T")
    
    if ( freeboundary ) then
      
      if ( .not. freeboundary_rst ) then
        write(*,*) 'WARNING: Restarting a simulation with freeboundary=.t. which was run with freeboundary=.f. so far.'
        if ( .not. freeb_change_indices ) then
          write(*,*) 'WARNING: Your free boundary simulation might be parallelized badly if you re-start a fixed'
          write(*,*) 'WARNING:   boundary simulation with free boundary, unless you set freeb_change_indices=.t.'
          write(*,*) 'WARNING:   from the very beginning of the simulation (grid construction)'
        end if
      end if
      
      if ( freeboundary_rst ) then
        call HDF5_char_reading(file_id,t_resistive_wall,"resistive_wall")
        resistive_wall_rst = (t_resistive_wall == "T")
      else
        resistive_wall_rst = .false.
      end if
      
      if ( resistive_wall .and. resistive_wall_rst ) then
        
        if (halo_freeb .and. CARIDDI_mode) then
          call HDF5_integer_reading(file_id,n_wall_curr_i,"n_wall_curr_i")
          call HDF5_integer_reading(file_id,n_wall_curr_e,"n_wall_curr_e")
          call HDF5_integer_reading(file_id, n_shell_curr_e, "n_shell_curr_e")
        else
         call HDF5_integer_reading(file_id,n_wall_curr,"n_wall_curr")
       endif
        
       call HDF5_integer_reading(file_id,n_dof_starwall,"n_dof_starwall")

       if (halo_freeb .AND. CARIDDI_mode) then

           if ( allocated(wall_curr_i) ) deallocate(wall_curr_i)
           allocate( wall_curr_i(n_wall_curr_i) )
           call HDF5_array1D_reading(file_id,wall_curr_i,"wall_curr_i")
          
           if ( allocated(dwall_curr_i) ) deallocate(dwall_curr_i)
           allocate( dwall_curr_i(n_wall_curr_i) )
           call HDF5_array1D_reading(file_id,dwall_curr_i,"dwall_curr_i")

           if ( allocated(wall_curr_e) ) deallocate(wall_curr_e)
           allocate( wall_curr_e(n_wall_curr_e) )
           call HDF5_array1D_reading(file_id,wall_curr_e,"wall_curr_e")
          
           if ( allocated(dwall_curr_e) ) deallocate(dwall_curr_e)
           allocate( dwall_curr_e(n_wall_curr_e) )
           call HDF5_array1D_reading(file_id,dwall_curr_e,"dwall_curr_e")

           if ( allocated(Ielec_J) ) deallocate(Ielec_J)
           allocate( Ielec_J(n_shell_curr_e) )
           call HDF5_array1D_reading(file_id,Ielec_J,"Ielec_J")

           if ( allocated(Ielec_w) ) deallocate(Ielec_w)
           allocate( Ielec_w(n_wall_curr_e) )
           call HDF5_array1D_reading(file_id,Ielec_w,"Ielec_w")

           if ( allocated(Ielec_J_old) ) deallocate(Ielec_J_old)
           allocate( Ielec_J_old(n_shell_curr_e) )
           call HDF5_array1D_reading(file_id,Ielec_J_old,"Ielec_J_old")

           if ( allocated(Ielec_w_old) ) deallocate(Ielec_w_old)
           allocate( Ielec_w_old(n_wall_curr_e) )
           call HDF5_array1D_reading(file_id,Ielec_w_old,"Ielec_w_old")

           if ( allocated(dIelec_J) ) deallocate(dIelec_J)
           allocate( dIelec_J(n_shell_curr_e) )
           call HDF5_array1D_reading(file_id,dIelec_J,"dIelec_J")

           if ( allocated(dIelec_w) ) deallocate(dIelec_w)
           allocate( dIelec_w(n_wall_curr_e) )
           call HDF5_array1D_reading(file_id,dIelec_w,"dIelec_w")

       else  

          if ( allocated(wall_curr) ) deallocate(wall_curr)
          allocate( wall_curr(n_wall_curr) )
          call HDF5_array1D_reading(file_id,wall_curr,"wall_curr")
           
          if ( allocated(dwall_curr) ) deallocate(dwall_curr)
          allocate( dwall_curr(n_wall_curr) )
          call HDF5_array1D_reading(file_id,dwall_curr,"dwall_curr")

       endif
        
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )

        old_dpsibnd_vec(:) = 0.d0
        if (halo_freeb) then 
          if ( allocated(old_dIelec_w) ) deallocate(old_dIelec_w)
          allocate( old_dIelec_w(n_wall_curr_e) )
          old_dIelec_w = 0.d0 
          if ( allocated(old_dIelec_J) ) deallocate(old_dIelec_J)
          allocate( old_dIelec_J(n_shell_curr_e) )
          old_dIelec_J = 0.d0 
        endif
        call HDF5_array1D_reading(file_id,old_dpsibnd_vec,"old_dpsibnd_vec")
        if (halo_freeb) then 
          call HDF5_array1D_reading(file_id,old_dIelec_J,"old_dIelec_J")
          call HDF5_array1D_reading(file_id,old_dIelec_w,"old_dIelec_w")
        endif
        call HDF5_integer_reading(file_id,n_coils,"n_coils")
        if ( index_start > 1 ) then

          if ( allocated(diag_coil_curr) ) deallocate(diag_coil_curr)
          if ( allocated(diag_coil_name) ) deallocate(diag_coil_name)
          call HDF5_integer_reading(file_id,n_diag_coil,"n_diag_coil")
          if (n_diag_coil > 0 ) then ! if sr%n_diag_coil > 0 the array will be allocated in evolve_wall_currents

            allocate( t_diag_coil_curr(index_start,n_diag_coil) )
            call HDF5_array2D_reading(file_id,t_diag_coil_curr,"diag_coil_curr")
            allocate( diag_coil_curr(index_start+nstep,n_diag_coil) )
            diag_coil_curr = 0.d0
            diag_coil_curr(1:index_start,:) = t_diag_coil_curr(1:index_start,:)
            deallocate(t_diag_coil_curr)

            allocate( t_diag_coil_name(n_diag_coil) )
            call HDF5_array1D_reading_char(file_id,t_diag_coil_name,"diag_coil_name")
            allocate( diag_coil_name(n_diag_coil) )
            diag_coil_name = transfer(t_diag_coil_name,diag_coil_name)
            deallocate(t_diag_coil_name)
            
          endif


          if ( allocated(pf_coil_curr) ) deallocate(pf_coil_curr)
          if ( allocated(pf_coil_name) ) deallocate(pf_coil_name)
          call HDF5_integer_reading(file_id,n_pf_coil,"n_pf_coil")
          if (n_pf_coil > 0 ) then    ! if sr%n_pol_coil > 0 the array will be allocated in evolve_wall_currents

            allocate( t_pf_coil_curr(index_start,n_pf_coil) )
            call HDF5_array2D_reading(file_id,t_pf_coil_curr,"pf_coil_curr")
            allocate( pf_coil_curr(index_start+nstep,n_pf_coil) )
            pf_coil_curr  = 0.d0
            pf_coil_curr(1:index_start,:) = t_pf_coil_curr(1:index_start,:)
            deallocate(t_pf_coil_curr)

            allocate( t_pf_coil_name(n_pf_coil) )
            call HDF5_array1D_reading_char(file_id,t_pf_coil_name,"pf_coil_name")
            allocate( pf_coil_name(n_pf_coil) )
            pf_coil_name = transfer(t_pf_coil_name,pf_coil_name)
            deallocate(t_pf_coil_name)

          endif

          if ( allocated(rmp_coil_curr) ) deallocate(rmp_coil_curr)
          if ( allocated(rmp_coil_name) ) deallocate(rmp_coil_name)
          call HDF5_integer_reading(file_id,n_rmp_coil,"n_rmp_coil")

          if (n_rmp_coil > 0 ) then! if sr%n_rmp_coil > 0 the array will be allocated in evolve_wall_currents

            allocate( t_rmp_coil_curr(index_start,n_rmp_coil) )
            call HDF5_array2D_reading(file_id,t_rmp_coil_curr,"rmp_coil_curr")
            allocate( rmp_coil_curr(index_start+nstep,n_rmp_coil) )
            rmp_coil_curr = 0.d0
            rmp_coil_curr(1:index_start,:) = t_rmp_coil_curr(1:index_start,:)
            deallocate(t_rmp_coil_curr)

            allocate( t_rmp_coil_name(n_rmp_coil) )
            call HDF5_array1D_reading_char(file_id,t_rmp_coil_name,"rmp_coil_name")
            allocate( rmp_coil_name(n_rmp_coil) )
            rmp_coil_name = transfer(t_rmp_coil_name,rmp_coil_name)
            deallocate(t_rmp_coil_name)

          endif

          if ( allocated(net_tor_wall_curr) ) deallocate(net_tor_wall_curr)
          allocate( t_net_tor_wall_curr(index_start) )
          call HDF5_array1D_reading(file_id,t_net_tor_wall_curr,"net_tor_wall_curr")
          allocate( net_tor_wall_curr(index_start+nstep))
          net_tor_wall_curr = 0.d0
          net_tor_wall_curr(1:index_start) = t_net_tor_wall_curr(1:index_start)
          deallocate(t_net_tor_wall_curr)

          
          if ( n_coils > 1  ) then
            if ( allocated(vert_FB_response)) deallocate(vert_FB_response)
            allocate( t_vert_FB_response(index_start,4) )
            call HDF5_array2D_reading(file_id,t_vert_FB_response,"vert_FB_response")
            allocate( vert_FB_response(index_start+nstep,4) )
            vert_FB_response = 0.d0
            vert_FB_response(1:index_start,:) = t_vert_FB_response(1:index_start,:)
            deallocate(t_vert_FB_response)
          endif
          
        end if
        
        if ( vacuum_debug .and. resistive_wall ) then
          write(*,*) 'DEBUG: Checksums'
          write(*,*) 'wall_curr', sum(abs(wall_curr))
          write(*,*) 'dwall_curr', sum(abs(dwall_curr))
          write(*,*) 'END: Checksums'
        end if
        
        wall_curr_initialized = .true.
        
      else if ( resistive_wall ) then
        
        write(*,*) 'WARNING: Continuing a JOREK simulation with resistive_wall=.f.'
        write(*,*) '   now with resistive_wall=.t. Wall currents will be initialized to zero.'
        wall_curr_initialized = .false.
        
      end if
      
      call HDF5_real_reading(file_id,current_FB_fact,'current_FB_fact')
      call HDF5_real_reading(file_id,dZ_axis_integral,'dZ_axis_integral')

      if ( n_coils /= 0 ) then
        if ( allocated(I_coils) ) deallocate(I_coils)
        allocate( I_coils(n_coils) )
        call HDF5_array1D_reading(file_id,I_coils,"I_coils")
      end if
      
    end if
     
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if

#endif

  end subroutine import_HDF5_restart_vacuum
  
  
  !> Export some vacuum-related data to the restart file
  subroutine export_restart_vacuum(file_handle, freeboundary, resistive_wall)
    
    ! --- Routine parameters
    integer, intent(in) :: file_handle
    logical, intent(in) :: freeboundary
    logical, intent(in) :: resistive_wall
    
    if ( freeboundary ) then
      
      write(file_handle) resistive_wall
      if ( resistive_wall ) then
        
        if (halo_freeb .AND. CARIDDI_mode) then

          if ( (.not. allocated(wall_curr_i)) .or. (.not. allocated(dwall_curr_i)) .or.                    &
          (.not. allocated(wall_curr_e)) .or. (.not. allocated(dwall_curr_e))      .or.  &
          (.not. allocated(old_dpsibnd_vec)) .or. (.not. allocated(old_dIelec_J)) .or. &
          (.not. allocated(old_dIelec_w)) .or. (.not. allocated(Ielec_J)) .or. &
          (.not. allocated(Ielec_w)) .or. (.not. allocated(Ielec_J_old))  .or. (.not. allocated(Ielec_w_old)) .or. & 
          (.not. allocated(dIelec_J)) .or. (.not. allocated(dIelec_w)) ) then
            write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
            stop
          end if


          write(file_handle) n_wall_curr_i, n_wall_curr_e, n_shell_curr_e, n_dof_starwall
        
          write(file_handle) wall_curr_i(:)
        
          write(file_handle) dwall_curr_i(:)
          write(file_handle) wall_curr_e(:)
        
          write(file_handle) dwall_curr_e(:)

          write(file_handle) Ielec_J(:)
          write(file_handle) Ielec_w(:)
          write(file_handle) Ielec_J_old(:)
          write(file_handle) Ielec_w_old(:)
          write(file_handle) dIelec_J(:)
          write(file_handle) dIelec_w(:)

        else

          if ( (.not. allocated(wall_curr)) .or. (.not. allocated(dwall_curr)) .or.                    &
          (.not. allocated(old_dpsibnd_vec)) ) then
            write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
            stop
          end if
        
          write(file_handle) n_wall_curr, n_dof_starwall
        
          write(file_handle) wall_curr(:)
        
          write(file_handle) dwall_curr(:)
        
        endif
        
        

        
      end if
      
      write(file_handle) current_FB_fact
      write(file_handle) dZ_axis_integral
      
      if ( (n_coils/=0) .and. (.not. allocated(I_coils)) ) then
        write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: I_coils not allocated.'
        stop
      end if
    
      write(file_handle) n_coils
      if ( n_coils /= 0 ) write(file_handle) I_coils(:)
      
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
    
  end subroutine export_restart_vacuum
  
  
  !> Export some vacuum-related data to the HDF5 restart file
  subroutine export_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)
    
    use phys_module, only: index_now
#ifdef USE_HDF5
    use hdf5
    use hdf5_io_module

    ! --- Routine parameters
    integer(HID_T), intent(in) :: file_id
#else
    integer,        intent(in) :: file_id
#endif
    logical,        intent(in) :: freeboundary
    logical,        intent(in) :: resistive_wall
    
#ifdef USE_HDF5
    ! --- Local variables
    character           :: t_freeboundary, t_resistive_wall
    real*8, allocatable :: t_diag_coil_curr(:,:), t_pf_coil_curr(:,:), t_rmp_coil_curr(:,:)
    real*8, allocatable :: t_net_tor_wall_curr(:)
    real*8, allocatable :: t_vert_FB_response(:,:)

    t_freeboundary = "F"
    if (freeboundary) t_freeboundary = "T"
    call HDF5_char_saving(file_id,t_freeboundary,"freeboundary"//char(0))
    
    if ( freeboundary ) then
    
      t_resistive_wall = "F"
      if (resistive_wall) t_resistive_wall = "T"
      call HDF5_char_saving(file_id,t_resistive_wall,"resistive_wall"//char(0))

      if ( resistive_wall ) then

        if (halo_freeb .AND. CARIDDI_mode) then
          if ( (.not. allocated(wall_curr_i)) .or. (.not. allocated(dwall_curr_i)) .or. &
                (.not. allocated(wall_curr_e)) .or. (.not. allocated(dwall_curr_e)) .or. &
           (.not. allocated(old_dpsibnd_vec)) .or. (.not. allocated(old_dIelec_w)) &
           .or. (.not. allocated(old_dIelec_J)) .or. (.not. allocated(Ielec_J)) .or. &
           (.not. allocated(Ielec_J_old)) .or. (.not. allocated(Ielec_w_old)) .or. &
          (.not. allocated(Ielec_w)) .or. (.not. allocated(dIelec_J)) .or. (.not. allocated(dIelec_w))) then
          write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
          stop
          end if
        else  
          if ( (.not. allocated(wall_curr)) .or. (.not. allocated(dwall_curr)) .or. &
            (.not. allocated(old_dpsibnd_vec)) ) then
            write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
            stop
          end if
        endif

        if (halo_freeb .AND. CARIDDI_mode) then
          call HDF5_integer_saving(file_id,n_wall_curr_i,"n_wall_curr_i"//char(0))
          call HDF5_integer_saving(file_id,n_wall_curr_e,"n_wall_curr_e"//char(0))
          call HDF5_integer_saving(file_id,n_shell_curr_e,"n_shell_curr_e"//char(0))
        else 
         call HDF5_integer_saving(file_id,n_wall_curr,"n_wall_curr"//char(0))
        endif

        call HDF5_integer_saving(file_id,n_dof_starwall,"n_dof_starwall"//char(0))

        if (halo_freeb .AND. CARIDDI_mode) then    
          call HDF5_array1D_saving(file_id,wall_curr_i,n_wall_curr_i,"wall_curr_i"//char(0))
          call HDF5_array1D_saving(file_id,dwall_curr_i,n_wall_curr_i,"dwall_curr_i"//char(0))
          call HDF5_array1D_saving(file_id,wall_curr_e,n_wall_curr_e,"wall_curr_e"//char(0))
          call HDF5_array1D_saving(file_id,dwall_curr_e,n_wall_curr_e,"dwall_curr_e"//char(0))
          call HDF5_array1D_saving(file_id,Ielec_J,n_shell_curr_e,"Ielec_J"//char(0))
          call HDF5_array1D_saving(file_id,Ielec_w,n_wall_curr_e,"Ielec_w"//char(0))
          call HDF5_array1D_saving(file_id,Ielec_J_old,n_shell_curr_e,"Ielec_J_old"//char(0))
          call HDF5_array1D_saving(file_id,Ielec_w_old,n_wall_curr_e,"Ielec_w_old"//char(0))
          call HDF5_array1D_saving(file_id,dIelec_J,n_shell_curr_e,"dIelec_J"//char(0))
          call HDF5_array1D_saving(file_id,dIelec_w,n_wall_curr_e,"dIelec_w"//char(0))
        else
          call HDF5_array1D_saving(file_id,wall_curr,n_wall_curr,"wall_curr"//char(0))
          call HDF5_array1D_saving(file_id,dwall_curr,n_wall_curr,"dwall_curr"//char(0))
        endif

        call HDF5_integer_saving(file_id,sr%n_diag_coils,"n_diag_coil"//char(0))
        call HDF5_integer_saving(file_id,sr%n_pol_coils,"n_pf_coil"//char(0))
        call HDF5_integer_saving(file_id,sr%n_rmp_coils,"n_rmp_coil"//char(0))
        if(sr%n_diag_coils .gt. 0) call HDF5_array1D_saving_char(file_id,diag_coil_name,sr%n_diag_coils,"diag_coil_name"//char(0))
        if(sr%n_pol_coils  .gt. 0) call HDF5_array1D_saving_char(file_id,pf_coil_name,  sr%n_pol_coils,"pf_coil_name"//char(0))
        if(sr%n_rmp_coils  .gt. 0) call HDF5_array1D_saving_char(file_id,rmp_coil_name, sr%n_rmp_coils,"rmp_coil_name"//char(0))
        

        if ( index_now > 0) then        

          allocate(t_net_tor_wall_curr(index_now))
          t_net_tor_wall_curr(1:index_now) = net_tor_wall_curr(1:index_now)
          call HDF5_array1D_saving(file_id,t_net_tor_wall_curr,index_now,"net_tor_wall_curr"//char(0))
          deallocate(t_net_tor_wall_curr)

          if ( sr%n_diag_coils > 0 ) then
            allocate(t_diag_coil_curr(index_now,sr%n_diag_coils))
            t_diag_coil_curr(1:index_now,:) = diag_coil_curr(1:index_now,:)
            call HDF5_array2D_saving(file_id,t_diag_coil_curr,index_now,sr%n_diag_coils,"diag_coil_curr"//char(0))
            deallocate(t_diag_coil_curr)
          end if

          if ( sr%n_pol_coils > 0 ) then
            allocate(t_pf_coil_curr(index_now,sr%n_pol_coils))
            t_pf_coil_curr(1:index_now,:) = pf_coil_curr(1:index_now,:)
            call HDF5_array2D_saving(file_id,t_pf_coil_curr,index_now,sr%n_pol_coils,"pf_coil_curr"//char(0))
            deallocate(t_pf_coil_curr)
          end if

          if ( sr%n_rmp_coils > 0 ) then
            allocate(t_rmp_coil_curr(index_now,sr%n_rmp_coils))
            t_rmp_coil_curr(1:index_now,:) = rmp_coil_curr(1:index_now,:)
            call HDF5_array2D_saving(file_id,t_rmp_coil_curr,index_now,sr%n_rmp_coils,"rmp_coil_curr"//char(0))
            deallocate(t_rmp_coil_curr)
          end if

          if ( sr%ncoil > 0 ) then
            allocate(t_vert_FB_response(index_now,4))
            t_vert_FB_response(1:index_now,:) = vert_FB_response(1:index_now,:)
            call HDF5_array2D_saving(file_id,t_vert_FB_response,index_now,4,"vert_FB_response"//char(0))
            deallocate(t_vert_FB_response)
          endif

        end if !--- index now

     end if !--- resistive wall

      call HDF5_array1D_saving(file_id,old_dpsibnd_vec,n_dof_starwall,"old_dpsibnd_vec"//char(0))
      if (halo_freeb .AND. CARIDDI_mode) then  
        call HDF5_array1D_saving(file_id,old_dIelec_J,n_shell_curr_e,"old_dIelec_J"//char(0))
        call HDF5_array1D_saving(file_id,old_dIelec_w,n_wall_curr_e,"old_dIelec_w"//char(0))
      endif  
      call HDF5_integer_saving(file_id,n_coils,"n_coils"//char(0))      
      call HDF5_real_saving(file_id,current_FB_fact,'current_FB_fact'//char(0))
      call HDF5_real_saving(file_id,dZ_axis_integral,'dZ_axis_integral'//char(0))
      if ( (n_coils/=0) .and. (.not. allocated(I_coils)) )  then
        write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: I_coils not allocated.'
        stop
      end if

      if ( n_coils /= 0 ) call HDF5_array1D_saving(file_id,I_coils,n_coils,"I_coils"//char(0))
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
#endif
  
  end subroutine export_HDF5_restart_vacuum
  
  
  !> Broadcast vacuum information between MPI processes
  subroutine broadcast_vacuum(my_id, resistive_wall)
    use mpi_mod
    implicit none
    
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    integer :: ierr, sz_diag(2), sz_pol(2), sz_rmp(2), sz_net, sz_VFB(2)

    if ( my_id == 0 ) then
      write(*,*) '*************************************'
      write(*,*) '*        Broadcast_vacuum           *'
    end if

    call MPI_BCAST(n_dof_starwall,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    sz_net     = 0
    sz_diag(:) = 0
    sz_pol(:)  = 0
    sz_rmp(:)  = 0
    sz_VFB(:)  = 0

    if ( resistive_wall ) then
      if (halo_freeb .and. CARIDDI_mode) then
        if ( my_id == 0 ) then
          write(*,*) '*************************************'
          write(*,*) '*        Broadcasting n_halo           *'
        end if
        call MPI_BCAST(n_wall_curr_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(n_wall_curr_e,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(n_shell_curr_e,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if ( my_id == 0 ) then
          write(*,*) '*************************************'
          write(*,*) '*     Done broadcasting n_halo      *'
        end if       
      else
        call MPI_BCAST(n_wall_curr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      endif
               if ( my_id == 0 ) then
          write(*,*) '*************************************'
          write(*,*) '*    Broadcasting n_coils           *'
        end if 
      call MPI_BCAST(n_coils,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      if (halo_freeb .and. CARIDDI_mode) then
                     if ( my_id == 0 ) then
          write(*,*) '*************************************'
          write(*,*) '*    n_wall_curr_i+n_wall_curr_e           *'
        end if 
        if (n_wall_curr_i + n_wall_curr_e == 0) RETURN
      else 
        if ( n_wall_curr == 0 ) return
      endif
                           if ( my_id == 0 ) then
          write(*,*) '*************************************'
          write(*,*) '*    DONE n_wall_curr_i+n_wall_curr_e           *'
        end if 
      if ( my_id == 0 ) then
        if ( allocated(net_tor_wall_curr) ) then
          sz_net = size(net_tor_wall_curr)
        end if
        if ( allocated(diag_coil_curr) ) then
          sz_diag(:) = (/ size(diag_coil_curr,1), size(diag_coil_curr,2) /)
        end if
        if ( allocated(pf_coil_curr) ) then
          sz_pol(:) = (/ size(pf_coil_curr,1), size(pf_coil_curr,2) /)
        end if
        if ( allocated(rmp_coil_curr) ) then
          sz_rmp(:) = (/ size(rmp_coil_curr,1), size(rmp_coil_curr,2) /)
        end if
        if ( allocated(vert_FB_response) ) then
          sz_VFB(:) = (/ size(vert_FB_response,1), size(vert_FB_response,2) /)
        end if
      end if
      call MPI_BCAST( sz_net,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sz_diag,  2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( sz_pol,  2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( sz_rmp,  2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( sz_VFB,  2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                                 if ( my_id == 0 ) then
          write(*,*) '*************************************'
          write(*,*) '*    About to enter my_id=/0           *'
        end if 
      if ( my_id /= 0 ) then

        if (halo_freeb .AND. CARIDDI_mode) then
          if ( my_id == 1 ) then
            write(*,*) '*************************************'
            write(*,*) '* Allocating halo-related variables *'
          end if 
          if ( allocated(wall_curr_i) ) deallocate(wall_curr_i)
          allocate( wall_curr_i(n_wall_curr_i) )
          
          if ( allocated(dwall_curr_i) ) deallocate(dwall_curr_i)
          allocate( dwall_curr_i(n_wall_curr_i) )

          if ( allocated(wall_curr_e) ) deallocate(wall_curr_e)
          allocate( wall_curr_e(n_wall_curr_e) )
          
          if ( allocated(dwall_curr_e) ) deallocate(dwall_curr_e)
          allocate( dwall_curr_e(n_wall_curr_e) )

          if (allocated(Ielec_J)) deallocate(Ielec_J)
          allocate(Ielec_J(n_shell_curr_e))

          if (allocated(Ielec_w)) deallocate(Ielec_w)
          allocate(Ielec_w(n_wall_curr_e))

          if (allocated(Ielec_J_old)) deallocate(Ielec_J_old)
          allocate(Ielec_J_old(n_shell_curr_e))

          if (allocated(Ielec_w_old)) deallocate(Ielec_w_old)
          allocate(Ielec_w_old(n_wall_curr_e))

          if (allocated(dIelec_J)) deallocate(dIelec_J)
          allocate(dIelec_J(n_shell_curr_e))

          if (allocated(dIelec_w)) deallocate(dIelec_w)
          allocate(dIelec_w(n_wall_curr_e))  

                    if ( my_id == 1 ) then
            write(*,*) '*************************************'
            write(*,*) '* DONE Allocating halo-related variables *'
          end if 
        else
          if ( allocated(wall_curr) ) deallocate(wall_curr)
          allocate( wall_curr(n_wall_curr) )
          
          if ( allocated(dwall_curr) ) deallocate(dwall_curr)
          allocate( dwall_curr(n_wall_curr) )
        endif


        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
        if (halo_freeb .AND. CARIDDI_mode) then  
                            if ( my_id == 1 ) then
            write(*,*) '*****************************************'
            write(*,*) '* Allocating old_dIelec_J, old_dIelec_w *'
          end if   
          if ( allocated(old_dIelec_J) ) deallocate(old_dIelec_J)
          allocate( old_dIelec_J(n_shell_curr_e) )
          if ( allocated(old_dIelec_w) ) deallocate(old_dIelec_w)
          allocate( old_dIelec_w(n_wall_curr_e) )
                               if ( my_id == 1 ) then
            write(*,*) '*****************************************'
            write(*,*) '* DONE Allocating old_dIelec_J, old_dIelec_w *'
          end if

        endif
                                if ( my_id == 1 ) then
            write(*,*) '*****************************************'
            write(*,*) '* Allocating net_tor_wall_curr *'
          end if
        

        if ( allocated(net_tor_wall_curr) ) deallocate(net_tor_wall_curr)
        if ( allocated(diag_coil_curr) )    deallocate(diag_coil_curr)
        if ( allocated(pf_coil_curr  ) )    deallocate(pf_coil_curr)
        if ( allocated(rmp_coil_curr ) )    deallocate(rmp_coil_curr)
        if ( allocated(vert_FB_response) )  deallocate(vert_FB_response)
        if (         sz_net  > 0 ) allocate( net_tor_wall_curr(sz_net) )
        if ( minval(sz_diag) > 0 ) then
          if ( allocated(diag_coil_name) )    deallocate(diag_coil_name)
          allocate( diag_coil_curr(sz_diag(1), sz_diag(2)) )
          allocate( diag_coil_name(sz_diag(2)) )
        endif
        if ( minval(sz_pol)  > 0 ) then
          if ( allocated(pf_coil_name) )      deallocate(pf_coil_name)
          allocate(   pf_coil_curr( sz_pol(1),  sz_pol(2)) )
          allocate(   pf_coil_name( sz_pol(2)) )
        endif
        if ( minval(sz_rmp)  > 0 ) then
          if ( allocated(rmp_coil_name) )     deallocate(rmp_coil_name)
          allocate(  rmp_coil_curr( sz_rmp(1),  sz_rmp(2)) )
          allocate(  rmp_coil_name( sz_rmp(2)) )
        endif
        if ( minval(sz_VFB) > 0 ) then
          allocate( vert_FB_response(sz_VFB(1), sz_VFB(2)) )
        endif
        if ( n_coils>0 ) then
          if (.not. allocated(I_coils)) then
            allocate( I_coils(n_coils))
            I_coils=0.d0
          end if
        end if

                                 if ( my_id == 0 ) then
            write(*,*) '*****************************************'
            write(*,*) '* DONE Allocating net_tor_wall_curr part my_id = 0 *'
          end if

                                 if ( my_id == 1 ) then
            write(*,*) '*****************************************'
            write(*,*) '* DONE Allocating net_tor_wall_curr part my_id =1 *'
          end if

      end if

                                                  if ( my_id == 0 ) then
            write(*,*) '*****************************************'
            write(*,*) '* BEFORE wall_curr_i, wall_curr_e my_id=0 *'
          end if

                                      if ( my_id == 1 ) then
            write(*,*) '*****************************************'
            write(*,*) '* BEFORE wall_curr_i, wall_curr_e my_id=1*'
          end if

      if (halo_freeb .AND. CARIDDI_mode) then
                if ( my_id == 1 ) then
            write(*,*) '*************************************'
            write(*,*) '* Broadcasting wall_curr_i, wall_curr_e my_id=1 *'
          end if 
                       if ( my_id == 0 ) then
            write(*,*) '*************************************'
            write(*,*) '* Broadcasting wall_curr_i, wall_curr_e my_id=0 *'
          end if 
        call MPI_BCAST(wall_curr_i,n_wall_curr_i,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dwall_curr_i,n_wall_curr_i,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wall_curr_e,n_wall_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dwall_curr_e,n_wall_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Ielec_J,n_shell_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Ielec_w,n_wall_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Ielec_J_old,n_shell_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Ielec_w_old,n_wall_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dIelec_J,n_shell_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dIelec_w,n_wall_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


                   if ( my_id == 1 ) then
            write(*,*) '*************************************'
            write(*,*) '* DONE Broadcasting wall_curr_i, wall_curr_e my_id=1 *'
          end if 
                                 if ( my_id == 0 ) then
            write(*,*) '*************************************'
            write(*,*) '* DONE Broadcasting wall_curr_i, wall_curr_e my_id=0 *'
          end if 
      else
        call MPI_BCAST(wall_curr,n_wall_curr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dwall_curr,n_wall_curr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      endif

      call MPI_BCAST(old_dpsibnd_vec,n_dof_starwall,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
      if (halo_freeb .and. CARIDDI_mode) then 
          if ( my_id == 1 ) then
            write(*,*) '*************************************'
            write(*,*) '* Broadcasting old_dIelec_J,old_dIelec_w my_id=1 *'
            write(*,*) 'Printing n_shell_curr_e for rank 1 =' , n_shell_curr_e
            write(*,*) 'Printing n_wall_curr_e for rank 1 =' , n_wall_curr_e
          end if 
          if ( my_id == 0 ) then
            write(*,*) '*************************************'
            write(*,*) '* Broadcasting old_dIelec_J,old_dIelec_w my_id=0 *'
            write(*,*) 'Printing n_shell_curr_e for rank 0 =' , n_shell_curr_e
            write(*,*) 'Printing n_wall_curr_e for rank 0 =' , n_wall_curr_e

          end if
        call MPI_BCAST(old_dIelec_J,n_shell_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
        call MPI_BCAST(old_dIelec_w,n_wall_curr_e,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
                  if ( my_id == 1 ) then
            write(*,*) '*************************************'
            write(*,*) '* DONE Broadcasting old_dIelec_J,old_dIelec_w my_id=1 *'
          end if 
          if ( my_id == 0 ) then
            write(*,*) '*************************************'
            write(*,*) '* DONE Broadcasting old_dIelec_J,old_dIelec_w my_id=0 *'
          end if
      endif
      call MPI_BCAST(wall_curr_initialized,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      if (         sz_net  > 0 ) call MPI_BCAST(net_tor_wall_curr, sz_net               ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if ( minval(sz_diag) > 0 ) call MPI_BCAST(diag_coil_curr,    sz_diag(1)*sz_diag(2),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if ( minval( sz_pol) > 0 ) call MPI_BCAST(  pf_coil_curr,    sz_pol(1)*sz_pol(2)  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)      
      if ( minval( sz_rmp) > 0 ) call MPI_BCAST( rmp_coil_curr,    sz_rmp(1)*sz_rmp(2)  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if ( minval(sz_diag) > 0 ) call MPI_BCAST(diag_coil_name, sz_diag(2)*COIL_NAME_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      if ( minval( sz_pol) > 0 ) call MPI_BCAST(  pf_coil_name,  sz_pol(2)*COIL_NAME_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      if ( minval( sz_rmp) > 0 ) call MPI_BCAST( rmp_coil_name,  sz_rmp(2)*COIL_NAME_LEN,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      if ( minval( sz_VFB) > 0 ) call MPI_BCAST(vert_FB_response,  sz_VFB(1)*sz_VFB(2)  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if ( n_coils > 0 ) call MPI_BCAST(I_coils,  n_coils,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    end if
    
    call MPI_BCAST(current_FB_fact,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dZ_axis_integral,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                  if ( my_id == 1 ) then
            write(*,*) '*************************************'
            write(*,*) '* DONE broadcast_vacuum my_id 1 *'
          end if 
                          if ( my_id == 0 ) then
            write(*,*) '*************************************'
            write(*,*) '* DONE broadcast_vacuum my_id 0 *'
          end if 
  end subroutine broadcast_vacuum



  !> Allocate the memory for a matrix distributed over the MPI tasks.
  subroutine alloc_distr(my_id, matrix, dim, row_wise)

    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer,              intent(in)     :: my_id
    type(t_distrib_mat),  intent(inout)  :: matrix
    integer,              intent(in)     :: dim(2)  !< Global dimension
    logical,              intent(in)     :: row_wise !< if  .true. - rowwise; .false. - columnwise    

    ! --- Local variables
    integer :: loc_size, ntasks, ierr

    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

    if (allocated(matrix%loc_mat)) call dealloc_distr(matrix)
    if (row_wise) then

      matrix%step      = dim(1)/ntasks
      loc_size         = matrix%step
      if (my_id==ntasks-1) loc_size = matrix%step + dim(1) - ntasks * matrix%step

      matrix%row_wise  = .true.
      matrix%ind_start = my_id*matrix%step+1
      matrix%ind_end   = my_id*matrix%step+loc_size

      allocate( matrix%loc_mat(loc_size, dim(2)) )
    else
      matrix%step      = dim(2)/ntasks
      loc_size         = matrix%step
      if (my_id==ntasks-1) loc_size = matrix%step + dim(2) - ntasks * matrix%step

      matrix%row_wise  = .false.
      matrix%ind_start = my_id*matrix%step+1
      matrix%ind_end   = my_id*matrix%step+loc_size

      allocate( matrix%loc_mat(dim(1), loc_size) )
    end if

    matrix%distrib    = .true.
    matrix%dim(1)     = dim(1) ! Save global matrix dimensions
    matrix%dim(2)     = dim(2)

  end subroutine alloc_distr



  !> Deallocate the memory of a distributed matrix and set all parameters to default values.
  subroutine dealloc_distr(matrix)
    implicit none
    ! --- Routine parameters
    type(t_distrib_mat),  intent(inout)  :: matrix

    if ( allocated(matrix%loc_mat) ) deallocate(matrix%loc_mat)
    
    ! --- Set all parameters to default values.
    matrix%dim(:)    = 0
    matrix%step      = 0
    matrix%ind_start = 0
    matrix%ind_end   = 0
    matrix%distrib   = .false.
    matrix%row_wise  = .false.

  end subroutine dealloc_distr



  !> Print information about a distributed matrix for debugging purposes.
  subroutine print_distr(my_id, matrix_name, matrix)

    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer,              intent(in)     :: my_id
    character(len=*),     intent(in)     :: matrix_name
    type(t_distrib_mat),  intent(inout)  :: matrix
  
    ! --- Local variables
    integer ::  ierr
    
    800 format(1x,3a)
    801 format(3x,a,2i12)
    802 format(3x,a,l12)
    803 format(5x,2i12,a,i6)
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)    
    if ( my_id == 0 ) then
      write(*,800) 'Distributed matrix ', trim(matrix_name), ':'
      write(*,801) 'dim:      ', matrix%dim
      write(*,802) 'row_wise: ', matrix%row_wise
      write(*,802) 'distrib:  ', matrix%distrib
      write(*,801) 'step:     ', matrix%step
      write(*,801) 'ind_start, ind_end:'
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    write(*,803) matrix%ind_start, matrix%ind_end, '@my_id=', my_id
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if ( my_id == 0 ) write(*,801) 'local dimensions:'
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    write(*,803) size(matrix%loc_mat,1), size(matrix%loc_mat,2) , '@my_id=', my_id
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  end subroutine print_distr

  !> Compute J_normal in the mapping points written by CARIDDI, then get I_electrode on shell and wall
  subroutine get_Ielec_halo(node_list, element_list, bnd_elm_list, bnd_node_list, mapping, my_id, i_tor_min, i_tor_max)
 
  use constants
  use basis_at_gaussian
  use mod_basisfunctions
  use mod_parameters, only: n_degrees
  use data_structure
  use phys_module
  use gauss

  implicit none
  ! Input, output variables
  type (type_node_list),        intent(in)    :: node_list
  type (type_element_list),     intent(in)    :: element_list
  type (type_bnd_element_list), intent(in)    :: bnd_elm_list   
  type (type_bnd_node_list),    intent(in)    :: bnd_node_list
  character(len=*),             intent(in)    :: mapping
  integer,                      intent(in)    :: my_id, i_tor_min, i_tor_max

  ! Local variables
  type(type_node)       :: node_k
  type(type_element)    :: elm_k
  type(type_bnd_element):: bndelem  
  real*8, allocatable   :: R(:), Z(:), R_s(:), Z_s(:), R_global(:), Z_global(:), phi_global(:)
  real*8                :: R_c, Z_c
  real*8, allocatable   :: zj(:), psi_s(:), rho(:), psi_s_mat(:,:), zj_mat(:,:)
  real*8, allocatable   :: T0(:), rho_s(:), T0_s(:), P0_s(:)
  real*8, allocatable   :: fourier_psi_s_0(:,:), fourier_zj_0(:,:), psi_s_0(:), psi_s_0_mat(:,:), zj_0(:), zj_0_mat(:,:)
  real*8, allocatable   :: J_normal_temp(:), J_normal(:), J_normal2(:), J_normal_dim(:,:),J_normal_rescale(:), J_normal_mat(:,:), Jn_back2real(:), Jn_CAR(:)
  real*8                :: phi_v, size_i(2), phi_c
  integer               :: k_vertex, k_dof, k_node, k_dir, m_bndelem, m_bndelem_arr(1), i ,j, index, i_vertex, i_node, i_start, i_node_bnd, i_dof
  integer               :: ms,nctrl, nsub_ctrl, npol, mp, in, l_tor, f_tor, i_resp_1, i_resp_2, m_pt, itor, npts_elm, npts_CAR, i_CAR_glob
  real*8                :: k_size,vec_inside(2), grad_t(2), sign_out, phi(n_plane), G(4, n_degrees), s, fact_mu0
  type(type_bnd_element):: bndelem_m
  real*8, allocatable   :: Jn_bez(:) , Jn_bez_mat(:,:), Bez2Real(:,:)                   !< output matrix with Bezier
  logical               :: reconvert_to_real, is_c, full_Jn, debug_Jn, check_solenoid, remove_axi_halo, enforce_kirch
  integer               :: mode_n(sr%n_tor), k_mode, n_modes, el, iv, idx,npol_2D, mt
  integer, allocatable  :: mode_u(:)
  real*8, allocatable   :: iel_real(:), s_real(:), phi_real(:), s_elm(:), phi_elm(:)
  integer, allocatable  :: iel_real_int(:), ielm_indexes(:)  
  real*8, allocatable   :: H1_CAR(:,:,:), H1_CAR_s(:,:,:), HZ_CAR(:,:)
  real*8, allocatable   :: I_elec_J_complete(:), I_elec_w_complete(:), i_elec_tor_J(:), i_elec_tor_w(:), i_elec_J_new(:), i_elec_w_new(:)
  real*8, allocatable   :: Jn_four_coeff(:,:), Jn_four_vec(:)
  complex*16, allocatable:: Jn_coeff_c(:,:), Jn_four_c_vec(:), i_halo_int_c(:), jphi_c(:,:), i_plasma_int_c(:)
  type(type_element)    :: element 
  type (type_node)      :: nodes(n_vertex_max)
  real*8                :: zj0(element_list%n_elements,n_gauss, n_gauss, n_plane ), xjac, wgauss_copy(n_gauss), wst
  real*8                :: x_g(element_list%n_elements,n_gauss, n_gauss), x_s(element_list%n_elements,n_gauss, n_gauss), x_t(element_list%n_elements,n_gauss, n_gauss), y_s(element_list%n_elements,n_gauss, n_gauss), y_t(element_list%n_elements,n_gauss, n_gauss)
  real*8                :: jphi(element_list%n_elements*n_gauss*n_gauss*n_plane), jphi_dim(element_list%n_elements*n_gauss*n_gauss*n_plane,1)
  complex*16            :: i_cplx
  integer               :: n_tor_elec_J, n_tor_elec_w, n_pol_elec_J, n_pol_elec_w
  if (my_id .eq. 0) write(*,*) 'Inside get_Ielec_halo'

  
  !!!!!!  TESTING !!!!!!!!!
  ! Reconvert for debugging?
  reconvert_to_real = .false.
  ! Compute full current or implement harmonic beat frequency mitigation (ONLY GAUSS)
  full_Jn = .true.  !!!!!!!!!! LEAVE TO TRUE FOR CONSISTENCY !!!!!!!!!!!!!!!!
  ! Debugging
  debug_Jn = .true.
  ! Check solenoidality
  check_solenoid = .false.
  ! Remove axisymmetry from electrode currents
  remove_axi_halo = .false.
  if (my_id .eq. 0) then 
    if (remove_axi_halo) write(*,*) 'Removing axisymmetric components from halo currents'
  endif
  ! Enforce Kirchhoff law to numerical precision
  enforce_kirch = .false.
  if (my_id .eq. 0) then 
    if (enforce_kirch) write(*,*) 'Electrodes: Enforcing Kirchhoff law to numerical precision'
  endif
  
  ! FACT_MU0
  fact_mu0 = 1.d0/mu_zero

  ! Define electrode currents (complete and full rank)

  !!!!!!! ================================================================================================
  !!!!!!!========= R.S. > PROBLEM: This repeated deallocation/allocation will slow things down! Fix it! ==========
  !!!!!! =================================================================================================
  if (allocated(I_elec_J_complete)) deallocate(I_elec_J_complete)
  allocate(I_elec_J_complete(sr%n_electrode_J))
  I_elec_J_complete = 0.d0

  if (allocated(I_elec_w_complete)) deallocate(I_elec_w_complete)
  allocate(I_elec_w_complete(sr%n_electrode_w))
  I_elec_w_complete = 0.d0

  if (allocated(Ielec_J)) deallocate(Ielec_J)
  allocate(Ielec_J(sr%n_electrode_J-1))
  Ielec_J = 0.d0

  if (allocated(Ielec_w)) deallocate(Ielec_w)
  allocate(Ielec_w(sr%n_electrode_w-1))
  Ielec_w = 0.d0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  select case (mapping)
    
    case ('equidistant')
      nsub_ctrl = 3
      nctrl = bnd_elm_list%n_bnd_elements * nsub_ctrl * n_plane

    case ('CARIDDI') 

      if (ALLOCATED(iel_real)) deallocate(iel_real)
      if (ALLOCATED(s_real)) deallocate(s_real)
      if (ALLOCATED(phi_real)) deallocate(phi_real)
      if (ALLOCATED(iel_real_int)) deallocate(iel_real_int)  
      allocate(iel_real(sr%n_pts_real_J), s_real(sr%n_pts_real_J), phi_real(sr%n_pts_real_J), iel_real_int(sr%n_pts_real_J))

      iel_real = sr%ielsphi_real_J(:,1)
      iel_real_int = int(iel_real)
      s_real   = sr%ielsphi_real_J(:,2)
      phi_real = sr%ielsphi_real_J(:,3)
      npts_CAR = size(s_real)
      

    case ('gauss') 
      nsub_ctrl = n_gauss
      if (my_id.eq.0) write(*,*) 'n_gauss = ', n_gauss
      nctrl = bnd_elm_list%n_bnd_elements * nsub_ctrl * n_plane
      npol  = nctrl / n_plane

  end select
  
  ! Match the actual toroidal mode numbers to the free boundary modes list 
  mode_n = mode_number(sr%i_tor) ! N.B. n_tor = number of JOREK modes, differs from sr%n_tor which is the number of freeb modes
 
  if (ALLOCATED(zj)) deallocate(zj)
  if (ALLOCATED(psi_s)) deallocate(psi_s)
  if (ALLOCATED(zj_mat)) deallocate(zj_mat)
  if (ALLOCATED(psi_s_mat)) deallocate(psi_s_mat)
  if (ALLOCATED(rho)) deallocate(rho)
  if (ALLOCATED(T0)) deallocate(T0)
  if (ALLOCATED(rho_s)) deallocate(rhO_s)
  if (ALLOCATED(T0_s)) deallocate(T0_s)
  if (ALLOCATED(P0_s)) deallocate(P0_s)
  if (ALLOCATED(fourier_psi_s_0)) deallocate(fourier_psi_s_0)
  if (ALLOCATED(fourier_zj_0)) deallocate(fourier_zj_0)
  if (ALLOCATED(psi_s_0)) deallocate(psi_s_0)
  if (ALLOCATED(zj_0)) deallocate(zj_0)
  if (ALLOCATED(psi_s_0_mat)) deallocate(psi_s_0_mat)
  if (ALLOCATED(zj_0_mat)) deallocate(zj_0_mat)   
  if (ALLOCATED(J_normal_temp)) deallocate(J_normal_temp)
  if (ALLOCATED(J_normal)) deallocate(J_normal)
  if (ALLOCATED(J_normal2)) deallocate(J_normal2)
  if (ALLOCATED(J_normal_mat)) deallocate(J_normal_mat)
  if (ALLOCATED(Jn_bez)) deallocate(Jn_bez)
  if (ALLOCATED(Jn_bez_mat)) deallocate(Jn_bez_mat)
  if (ALLOCATED(Bez2Real)) deallocate(Bez2Real)
  if (ALLOCATED(Jn_back2real)) deallocate(Jn_back2real)
  if (ALLOCATED(R)) deallocate(R)
  if (ALLOCATED(R_s)) deallocate(R_s)
  if (ALLOCATED(Z)) deallocate(Z)
  if (ALLOCATED(Z_s)) deallocate(Z_s)
  if (ALLOCATED(R_global)) deallocate(R_global)
  if (ALLOCATED(Z_global)) deallocate(Z_global)
  if (ALLOCATED(phi_global)) deallocate(phi_global)

  select case (mapping)
    case ('equidistant') 
      allocate(zj(nsub_ctrl*n_plane), psi_s(nsub_ctrl*n_plane), psi_s_mat(nsub_ctrl*n_plane,1), &
      zj_mat(nsub_ctrl*n_plane,1), rho(nsub_ctrl*n_plane), T0(nsub_ctrl*n_plane), &
      rho_s(nsub_ctrl*n_plane) , T0_s(nsub_ctrl*n_plane) , P0_s(nsub_ctrl*n_plane), fourier_psi_s_0(nsub_ctrl,1), &
      fourier_zj_0(nsub_ctrl,1),  psi_s_0(nsub_ctrl*n_plane)  , zj_0(nsub_ctrl*n_plane), &
      psi_s_0_mat(nsub_ctrl*n_plane,1), zj_0_mat(nsub_ctrl*n_plane,1), Jn_back2real(nctrl), &
      R(nsub_ctrl), Z(nsub_ctrl), R_s(nsub_ctrl), Z_s(nsub_ctrl), R_global(nctrl), Z_global(nctrl), phi_global(nctrl), &
      J_normal_temp(nsub_ctrl*n_plane) , J_normal(nctrl), J_normal2(nctrl), J_normal_rescale(nctrl),J_normal_mat(nctrl,1), Jn_bez(n_dof_bnd*sr%n_tor) , Jn_bez_mat(n_dof_bnd*sr%n_tor,1), Bez2Real(nctrl, n_dof_bnd*sr%n_tor)   )

      
      J_normal = 0.d0 ; J_normal2=0.d0; J_normal_rescale=0.d0; Jn_bez = 0.d0; Jn_bez_mat = 0.d0; J_normal_mat=0.d0 ; Bez2Real = 0.d0; Jn_back2real=0.d0
      R_global = 0.d0; Z_global = 0.d0; phi_global = 0.d0

    case ('gauss') 
      allocate(zj(nsub_ctrl*n_plane), psi_s(nsub_ctrl*n_plane), psi_s_mat(nsub_ctrl*n_plane,1), &
      zj_mat(nsub_ctrl*n_plane,1), rho(nsub_ctrl*n_plane), T0(nsub_ctrl*n_plane), &
      rho_s(nsub_ctrl*n_plane) , T0_s(nsub_ctrl*n_plane) , P0_s(nsub_ctrl*n_plane), fourier_psi_s_0(nsub_ctrl,1), &
      fourier_zj_0(nsub_ctrl,1),  psi_s_0(nsub_ctrl*n_plane)  , zj_0(nsub_ctrl*n_plane), &
      psi_s_0_mat(nsub_ctrl*n_plane,1), zj_0_mat(nsub_ctrl*n_plane,1), Jn_back2real(nctrl), &
      R(nsub_ctrl), Z(nsub_ctrl), R_s(nsub_ctrl), Z_s(nsub_ctrl), R_global(nctrl), Z_global(nctrl), phi_global(nctrl), &
      J_normal_temp(nsub_ctrl*n_plane) , J_normal(nctrl), J_normal2(nctrl), J_normal_rescale(nctrl), J_normal_mat(nctrl,1), Jn_bez(n_dof_bnd*sr%n_tor) , Jn_bez_mat(n_dof_bnd*sr%n_tor,1), Bez2Real(nctrl, n_dof_bnd*sr%n_tor)   )

      
      J_normal = 0.d0 ; J_normal2 = 0.d0; J_normal_rescale=0.d0; Jn_bez = 0.d0; Jn_bez_mat = 0.d0; J_normal_mat=0.d0 ; Bez2Real = 0.d0; Jn_back2real=0.d0
      R_global = 0.d0; Z_global = 0.d0; phi_global = 0.d0

  end select
   
  !--- go through the boundary elements
  do m_bndelem = 1, bnd_elm_list%n_bnd_elements
    
    bndelem = bnd_elm_list%bnd_element(m_bndelem) ! Type bnd_element
    elm_k   = element_list%element(bndelem%element) ! Type element
    
    select case (mapping)
    !--- calculate values at gaussian points on the element
      case('gauss')
        R     = 0.d0;      Z = 0.d0;   R_s = 0.d0;   Z_s = 0.d0      
        zj    = 0.d0;      psi_s = 0.d0;   rho = 0.d0;   T0 = 0.d0 ; 
        rho_s = 0.d0;      T0_s = 0.d0;  P0_s = 0.d0 ; 
        fourier_psi_s_0 = 0.d0 ; fourier_zj_0 = 0.d0; psi_s_0 = 0.d0; zj_0 = 0.d0 ; J_normal_temp = 0.d0
        psi_s_mat = 0.d0; zj_mat = 0.d0; psi_s_0_mat = 0.d0 ; zj_0_mat = 0.d0
    
      case('CARIDDI')
        ! In the case of CARIDDI, we have to first know how many points per element we have for this particular element
        ! Find indexes in the ielm vector
        m_bndelem_arr(1) = m_bndelem
        call find_loc_indexes(iel_real_int, m_bndelem_arr, ielm_indexes)
        npts_elm = size(ielm_indexes)
        ! Proceed with allocation
        if (allocated(s_elm)) deallocate(s_elm)
        if (allocated(phi_elm)) deallocate(phi_elm)
        if (allocated(H1_CAR)) deallocate(H1_CAR)
        if (allocated(H1_CAR_s)) deallocate(H1_CAR_s)
        if (allocated(R)) deallocate(R)
        if (allocated(Z)) deallocate(Z)
        if (allocated(R_s)) deallocate(R_s)
        if (allocated(Z_s)) deallocate(Z_s)
        if (allocated(HZ_CAR)) deallocate(HZ_CAR)
        if (allocated(psi_s)) deallocate(psi_s)
        if (allocated(zj)) deallocate(zj)
        if (allocated(rho)) deallocate(rho)
        if (allocated(rho_s)) deallocate(rho_s)
        if (allocated(T0)) deallocate(T0)
        if (allocated(T0_s)) deallocate(T0_s)
        if (allocated(P0_s)) deallocate(P0_s)


        allocate(s_elm(npts_elm), phi_elm(npts_elm), H1_CAR(2,n_degrees_1d, npts_elm), H1_CAR_s(2, n_degrees_1d, npts_elm), &
                  R(npts_elm), Z(npts_elm), R_s(npts_elm), Z_s(npts_elm), HZ_CAR(sr%n_tor, npts_elm), &
                  psi_s(npts_elm), zj(npts_elm), rho(npts_elm), rho_s(npts_elm), T0(npts_elm), T0_s(npts_elm), P0_s(npts_elm))

        s_elm = 0.d0 ; phi_elm = 0.d0 ; H1_CAR = 0.d0 ; H1_CAR_s = 0.d0; R = 0.d0; Z = 0.d0; R_s = 0.d0; Z_s = 0.d0
        HZ_CAR = 0.d0 ; psi_s = 0.d0 ; zj = 0.d0; rho=0.d0; rho_s=0.d0; T0=0.d0; T0_s=0.d0; P0_s=0.d0
    
        ! Find s and phi for the element
        s_elm    = s_real(ielm_indexes)
        phi_elm  = phi_real(ielm_indexes)
    
    end select 

    do k_vertex = 1, 2
      do k_dof = 1, 2

        k_node   = bndelem%vertex(k_vertex)
        k_dir    = bndelem%direction(k_vertex,k_dof)
        k_size   = bndelem%size(k_vertex,k_dof)
        node_k   = node_list%node(k_node)

        select case (mapping)
          case ('gauss')     
            R  (:)   = R  (:)  + node_k%x(1,k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
            Z  (:)   = Z  (:)  + node_k%x(1,k_dir,2) * k_size * H1  (k_vertex,k_dof,:)
            R_s(:)   = R_s(:)  + node_k%x(1,k_dir,1) * k_size * H1_s(k_vertex,k_dof,:)
            Z_s(:)   = Z_s(:)  + node_k%x(1,k_dir,2) * k_size * H1_s(k_vertex,k_dof,:)

            do mp=1,n_plane
              do in=1,n_tor  ! USING JOREK HARMONICS!!!!! AND NOT STARWALL
                psi_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) = psi_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + node_k%values(in,k_dir,1) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,mp)
                zj   ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)   = zj   ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + node_k%values(in,k_dir,3) * k_size*H1(k_vertex,k_dof,:) * HZ(in,mp)
                if (jorek_model > 200 ) then
                  rho  ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) = rho  ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + node_k%values(in,k_dir,5) * k_size*  H1(k_vertex,k_dof,:) * HZ(in,mp)
                  rho_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) = rho_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + node_k%values(in,k_dir,5) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,mp)
                  T0   ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) = T0   ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + node_k%values(in,k_dir,6) * k_size*  H1(k_vertex,k_dof,:) * HZ(in,mp)
                  T0_s ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) = T0_s ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + node_k%values(in,k_dir,6) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,mp)
                else
                  rho = 1.d0; rho_s = 0.d0; T0 = 0.d0; T0_s = 0.d0 
                endif        
              enddo
            enddo
            
          case('CARIDDI')
        
            do i = 1, npts_elm
              s = s_elm(i)
              call basisfunctions1(s, H1_CAR(1:2,1:n_degrees_1d,i), H1_CAR_s(1:2,1:n_degrees_1d,i))
            end do
            ! All these arrays have the size of ielm_indexes. 
            R  (:)   = R  (:)  + node_k%x(1,k_dir,1) * k_size * H1_CAR  (k_vertex,k_dof,:)
            Z  (:)   = Z  (:)  + node_k%x(1,k_dir,2) * k_size * H1_CAR  (k_vertex,k_dof,:)
            R_s(:)   = R_s(:)  + node_k%x(1,k_dir,1) * k_size * H1_CAR_s(k_vertex,k_dof,:)
            Z_s(:)   = Z_s(:)  + node_k%x(1,k_dir,2) * k_size * H1_CAR_s(k_vertex,k_dof,:)
          
            ! Get toroidal basis functions
            ! - note that the basis for physics variables and coordinates have opposite signs
            do j = 1,npts_elm     
              phi_c = -phi_elm(j)
              if ( mode_n(1) .eq. 0) then ! First mode is n=0
                do i=1,sr%n_tor,2
                  if (sr%i_tor(i) .eq. 1) then 
                    HZ_CAR(i,j)    = 1.d0
                  else ! If we are at 3,5,7...
                    HZ_CAR(i-1,j) = cos(mode_n(i-1)  *phi_c)
                    HZ_CAR(i,j)= sin(mode_n(i)*phi_c)
                  end if
                end do
              else ! n_tor is even because we are without n=0
                do i=2,sr%n_tor,2
                  HZ_CAR(i-1,j) = cos(mode_n(i-1)  *phi_c)
                  HZ_CAR(i,j)= sin(mode_n(i)*phi_c)
                end do
              end if
            end do
           
            ! Values are placed sequentially in these local arrays for psi, zj, rho, T0
            ! The order is the one the ielm_indexes are placed at!
            do i= 1,npts_elm
              do in = 1,sr%n_tor
                psi_s(i) =  psi_s(i) + node_k%values(sr%i_tor(in),k_dir,1) * k_size*H1_CAR_s(k_vertex,k_dof,i) * HZ_CAR(in,i)
                zj(i)    =  zj(i) + node_k%values(sr%i_tor(in),k_dir,3) * k_size*H1_CAR(k_vertex,k_dof,i) * HZ_CAR(in,i)
                if (jorek_model > 200 ) then
                  rho  (i) = rho  (i) + node_k%values(sr%i_tor(in),k_dir,5) * k_size*  H1_CAR(k_vertex,k_dof,i) * HZ_CAR(in,i)
                  rho_s(i) = rho_s(i) + node_k%values(sr%i_tor(in),k_dir,5) * k_size*H1_CAR_s(k_vertex,k_dof,i) * HZ_CAR(in,i)
                  T0   (i) = T0   (i) + node_k%values(sr%i_tor(in),k_dir,6) * k_size*  H1_CAR(k_vertex,k_dof,i) * HZ_CAR(in,i)
                  T0_s (i) = T0_s (i) + node_k%values(sr%i_tor(in),k_dir,6) * k_size*H1_CAR_s(k_vertex,k_dof,i) * HZ_CAR(in,i)
                else
                  rho = 1.d0; rho_s = 0.d0; T0 = 0.d0; T0_s = 0.d0 
                endif 
              end do
            end do
            
        end select
      enddo
    enddo
    
    !--- Find out correct sign of the normal (it has to point outwards the domain)
    !---------------------------------------------------------------------------------- 
    ! --- Calculate an inside point on the element to calculate the
    ! direction of bnd normals
    call basisfunctions(xgauss(2),xgauss(2), G)  
    R_c = 0.d0 ;  Z_c = 0.d0 
    do i = 1, n_vertex_max
      do j = 1, n_degrees
        node_k = node_list%node(elm_k%vertex(i)) 
        R_c    = R_c + node_k%x(1,j,1) * elm_k%size(i,j) * G(i,j)
        Z_c    = Z_c + node_k%x(1,j,2) * elm_k%size(i,j) * G(i,j)
      enddo
    enddo  
    vec_inside = (/ R_c - R(2), Z_c - Z(2) /)       ! vector pointing towards the domain
    grad_t     = (/ -Z_s(2) , R_s(2) /)     ! gradient of the coordinate t (normal to the boundary here)
    ! grad_t could point either outwards or inwards depending on the orientation of the bnd_element!
    sign_out   = -1.d0 * sign( 1.d0, ( vec_inside(1)*grad_t(1) + vec_inside(2)*grad_t(2) ) )  
    !--------------------------------------------------------------------------------
    select case (mapping)
      case ('gauss')
        do mp=1, n_plane
          P0_s ((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)     = T0((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) * rho_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + &
                                                        T0_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) * rho((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)
          phi   (mp)      = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
        enddo

        if (full_Jn) then

          ! --- Computation of full normal current

          ! do mp = 1, n_plane

          !   J_normal_temp((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) = ((zj((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) * &
          !                                                 psi_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)) + &
          !                                                 + R(:)**2.d0 * P0_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)) &
          !                                                 /sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out
          ! enddo

          ! do mp = 1,n_plane ! Fill up the extended Jnorm by first going along the poloidal angle and then proceeding to next poloidal plane (along phi)

          !   J_normal( (m_bndelem-1)*nsub_ctrl+1+(mp-1)*npol : m_bndelem*nsub_ctrl+(mp-1)*npol ) = J_normal_temp((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)                          

          ! enddo

          do mp = 1,n_plane ! Fill up the extended Jnorm by first going along the poloidal angle and then proceeding to next poloidal plane (along phi)

            J_normal( (m_bndelem-1)*nsub_ctrl+1+(mp-1)*npol : m_bndelem*nsub_ctrl+(mp-1)*npol ) = ((zj((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) * &
                                                                                                psi_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)) + &
                                                                                                + R(:)**2.d0 * P0_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)) &
                                                                                                /sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out         
          enddo  


        else
        ! HARMONICS SEPARATION (OLD)
          ! ! ---- Implementation of Jn = 1/Bphi ( Bn(n=0) * Jphi + Jphi(n=0) * Bn )

          ! ! ! Extract n=0 harmonics for field and toroidal current
          ! psi_s_mat(:,1) = psi_s(:)
          ! zj_mat(:,1)    = zj(:)
          
          ! call calc_fourier_coeff_single(psi_s_mat,  n_plane, nsub_ctrl, 1, 1, phi, fourier_psi_s_0, my_id  )
          ! call calc_fourier_coeff_single( zj_mat, n_plane, nsub_ctrl , 1, 1, phi,  fourier_zj_0, my_id    )
      
          ! ! Re-convert the n=0 to real space
          ! call fourier2realspace(fourier_psi_s_0, 1, nsub_ctrl*n_plane, nsub_ctrl, n_plane, psi_s_0_mat)
          ! call fourier2realspace(fourier_zj_0, 1, nsub_ctrl*n_plane,   nsub_ctrl, n_plane,    zj_0_mat)

          ! psi_s_0(:) = psi_s_0_mat(:,1) ; zj_0(:) = zj_0_mat(:,1)

          ! ! Output variables are arrays in real space, on the whole toroidal "slice" identified by the nsub_ctrl poloidal points

          ! ! Loop over poloidal direction first, then toroidal
          ! do mp=1,n_plane
          !   J_normal_temp((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) = (zj_0((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) * psi_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + &
          !                                               psi_s_0((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)*zj((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)  &
          !                                               - zj_0((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) * psi_s_0((mp-1)*nsub_ctrl+1:mp*nsub_ctrl) + &
          !                                               R(:)**2.d0 * P0_s((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)) /sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out
          ! enddo

          ! do mp = 1,n_plane
          !   J_normal( (m_bndelem-1)*nsub_ctrl+1+(mp-1)*npol : m_bndelem*nsub_ctrl+(mp-1)*npol ) = J_normal_temp((mp-1)*nsub_ctrl+1:mp*nsub_ctrl)                                            
          ! enddo

        endif

        do mp = 1,n_plane
          R_global( (m_bndelem-1)*nsub_ctrl+1+(mp-1)*npol : m_bndelem*nsub_ctrl+(mp-1)*npol ) = R(:)                                            
          Z_global( (m_bndelem-1)*nsub_ctrl+1+(mp-1)*npol : m_bndelem*nsub_ctrl+(mp-1)*npol ) = Z(:)
          phi_global((m_bndelem-1)*nsub_ctrl+1+(mp-1)*npol : m_bndelem*nsub_ctrl+(mp-1)*npol) =  phi(mp)
        enddo

      case ('CARIDDI')
        do i = 1,npts_elm
          P0_s (i)  = T0(i) * rho_s(i) + T0_s(i) * rho(i)
        end do
        
        ! Compute normal current
        do i = 1,npts_elm
          i_CAR_glob = ielm_indexes(i)
          Jn_CAR(i_CAR_glob) = ((zj(i) * psi_s(i)) + R(i)**2.d0 * P0_s(i)) / &
                               sqrt(R_s(i)**2.d0 + Z_s(i)**2.d0) / (F0*R(i)) * sign_out
        end do

          Jn_CAR = fact_mu0 * Jn_CAR
    end select 

  enddo ! end bnd elements loop

  if (debug_Jn) then
    if (my_id == 0) then
      OPEN(UNIT=12, FILE="Jn_gauss_BEFORE_RESCALE.dat", ACTION="write", STATUS="replace")
      DO i=1,nctrl
        WRITE(12,*) J_normal(i)
      END DO
      CLOSE(12)

    endif
  endif
  ! --- RESCALE BY 1/mu_0
  J_normal_rescale(:) = fact_mu0 * J_normal(:)

  if (debug_Jn) then
    if (my_id == 0) then
      OPEN(UNIT=12, FILE="Jn_gauss_AFTER_RESCALE.dat", ACTION="write", STATUS="replace")
      DO i=1,nctrl
        WRITE(12,*) J_normal_rescale(i)
      END DO
      CLOSE(12)
    endif
  endif

  ! CHECK FOR HARMONIC-BY-HARMONIC SOLENOIDALITY
  if (check_solenoid) then

    if (.not. ALLOCATED(Jn_four_coeff)) allocate(Jn_four_coeff(npol,1))
    if (.not. ALLOCATED(Jn_coeff_c)) allocate(Jn_coeff_c(npol,1))
    if (.not. ALLOCATED(J_normal_dim)) allocate(J_normal_dim(nctrl,1))
    if (.not. ALLOCATED(Jn_four_c_vec)) allocate(Jn_four_c_vec(npol))   

    ! Separate Fourier harmonics
    J_normal_dim(:,1) = J_normal_rescale 
    Jn_four_coeff(:,:) = 0.d0
    Jn_coeff_c(:,:) = (0.d0,0.d0)
    Jn_four_c_vec(:)   = (0.d0, 0.d0)

    ! Produce array with the actual mode numbers and no repetitions (mode_n)
    if (.not. allocated(mode_u)) allocate(mode_u(sr%n_tor))
    k_mode=1
    mode_u(k_mode) = mode_n(1)
    do i = 2 , size(mode_n)
      if (mode_n(i) .ne. mode_n(i-1) ) then 
        k_mode = k_mode+1
        mode_u(k_mode) = mode_n(i)
      endif 
    enddo 
    mode_u = mode_u(1:k_mode)
    n_modes=k_mode

    if (.not. ALLOCATED(I_halo_int_c)) allocate(I_halo_int_C(n_modes))
    I_halo_int_c(:)      = (0.d0,0.d0)

    do i = 1, n_modes
      write(*,*) 'size(mode_u) check =' , size(mode_u)
      call calc_fourier_coeff_complex(J_normal_dim, n_plane, npol, mode_u(i), 1, phi, Jn_coeff_c, my_id )

      ! Integrate coefficient along poloidal length
      Jn_four_c_vec(:) = Jn_coeff_c(:,1)

      do m_bndelem=1,bnd_elm_list%n_bnd_elements
        
        bndelem = bnd_elm_list%bnd_element(m_bndelem) ! Type bnd_element

        ! Calculate radial coordinate
        do k_vertex = 1, 2
          do k_dof = 1, 2
            k_node   = bndelem%vertex(k_vertex)
            k_dir    = bndelem%direction(k_vertex,k_dof)
            k_size   = bndelem%size(k_vertex,k_dof)
            node_k   = node_list%node(k_node)
            R  (:)   = R  (:)  + node_k%x(1,k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
          enddo 
        enddo

        do ms=1,nsub_ctrl
          ! Should the jacobian still be canceled out because it compensates with dpsi/ds derivative?
          I_halo_int_c(i) = I_halo_int_c(i) + wgauss(ms) * R(ms) * Jn_four_c_vec((m_bndelem-1)*nsub_ctrl+ms)
        enddo

      enddo !bndelem loop
    enddo ! modes loop

    ! Take Jphi in real space and convert it to Fourier
    do el=1,element_list%n_elements
      element = element_list%element(el)

      do iv = 1, n_vertex_max
        i_node     = element%vertex(iv)
        call make_deep_copy_node(node_list%node(i_node), nodes(iv))
      enddo
      
      ! Fill in values for R,Z and zj
      zj0(:,:,:,:) = 0.d0
      x_g(:,:,:)   = 0.d0 
      x_s(:,:,:)   = 0.d0 
      x_t(:,:,:)   = 0.d0 
      y_s(:,:,:)   = 0.d0 
      y_t(:,:,:)   = 0.d0

      do i=1,n_vertex_max
        do j=1,n_degrees
          do ms=1, n_gauss
            do mt=1, n_gauss
              do mp=1,n_plane
                !do in=1,n_coord_tor
                  x_g(el,ms,mt) = x_g(el,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H(i,j,ms,mt) !* HZ_coord(in,mp)
                  x_s(el,ms,mt) = x_s(el,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_s(i,j,ms,mt) !* HZ_coord(in,mp)
                  x_t(el,ms,mt) = x_t(el,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_t(i,j,ms,mt) !* HZ_coord(in,mp)
                  y_s(el,ms,mt) = y_s(el,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_s(i,j,ms,mt) !* HZ_coord(in,mp)
                  y_t(el,ms,mt) = y_t(el,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_t(i,j,ms,mt) !* HZ_coord(in,mp)

                !enddo
                do in=1,n_tor
                  zj0(el,ms,mt,mp) = zj0(el,ms,mt,mp) + nodes(i)%values(in,j,var_zj) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                enddo 
              enddo 
            enddo 
          enddo 
        enddo 
      enddo
    enddo 
    
    ! Rearrange zj0 so that it's a vector scanning points first poloidally and then toroidally
    jphi(:) = 0.d0
    jphi_dim(:,:) = 0.d0

    npol_2D = element_list%n_elements * n_gauss * n_gauss ! ARE n_gauss THE TOTAL GAUSS POINTS OR THE GAUSS POINTS PER DIRECTION?
    do mp=1,n_plane
      do el =1,element_list%n_elements
        do ms=1,n_gauss 
          do mt=1,n_gauss
            idx = (mp-1)*npol_2D + (el-1)*n_gauss*n_gauss + (ms-1)*n_gauss + mt
            jphi(idx) = zj0(el,ms,mt,mp) !/ x_g(el,ms,mt) It simplifies with the later multiplication?
          enddo 
        enddo 
      enddo 
    enddo

    jphi_dim(:,1) = jphi(:)
    if (.not. allocated(jphi_c)) allocate(jphi_c(npol_2D,1))
    jphi_c = (0.d0,0.d0)
    if (.not. ALLOCATED(I_plasma_int_c)) allocate(I_plasma_int_C(n_modes))
    I_plasma_int_c(:)      = (0.d0,0.d0)

    open(unit=12, file="jphi_dim.dat", action="write", status="replace")
    do j = 1, npol_2D*n_plane
      write(12,*) jphi_dim(j,1)
    enddo 
    close(12)

    do i=1,n_modes

      call calc_fourier_coeff_complex(jphi_dim, n_plane, npol_2D, mode_u(i), 1, phi, JPHI_c, my_id )
      if (my_id .eq. 0) then
      open(unit=12, file="jphi_c.dat", action="write", status="replace")
      do j = 1, npol_2D
        write(12,*) real(jphi_c(j,1)) , aimag(jphi_c(j,1))
      enddo 
      close(12)
      endif
      ! For each mode, perform the surface integral over the 2D domain 
      do el=1,element_list%n_elements
        do ms=1,n_gauss 
          do mt=1,n_gauss 
            xjac = x_s(el,ms,mt)*y_t(el,ms,mt) - x_t(el,ms,mt)*y_s(el,ms,mt)
            !write(*,*) 'xjac = ', xjac
            idx = (el-1)*n_gauss*n_gauss + (ms-1)*n_gauss + mt
            wgauss_copy = wgauss
            wst  = wgauss_copy(ms)*wgauss_copy(mt)
            !write(*,*) 'wst = ', wst

            I_plasma_int_c(i) = I_plasma_int_c(i) + wst * xjac * jphi_c(idx,1)
          enddo 
        enddo 
      enddo
      i_cplx = (0.d0,1.d0)
      I_plasma_int_c(i) = i_cplx*i*I_plasma_int_c(i)
    enddo
    
    ! PRINT THE TWO ARRAYS
    if (my_id .eq. 0) then
      open(unit=12, file="I_halo_int_c.dat", action="write", status="replace")
      do i = 1, n_modes
        write(12,*) real(I_halo_int_c(i)) , aimag(i_halo_int_c(i))
      enddo 
      close(12)
    endif

    if (my_id .eq. 0) then
      open(unit=12, file="I_jphi_int_c.dat", action="write", status="replace")
      do i = 1, n_modes
        write(12,*) real(I_plasma_int_c(i)) , aimag(i_plasma_int_c(i))
      enddo 
      close(12)
    endif



  endif ! FIN SOLENOIDALITY TEST


  select case (mapping)

    case ('CARIDDI')
      
      ! ---------  ELECTRODE CURRENTS ----------------
      
      ! Shell
      do i =1,sr%n_electrode_J
        I_elec_J_complete(i) = sum(sr%Jmat_real_J(i,:)*Jn_CAR(:))
      end do

      ! Wall
      do i =1,sr%n_electrode_w
        I_elec_w_complete(i) = sum(sr%Jmat_real_w(i,:)*Jn_CAR(:))
      end do

      Ielec_J = I_elec_J_complete(1:sr%n_electrode_J-1)
      Ielec_w = I_elec_w_complete(1:sr%n_electrode_w-1)
      
      if (debug_Jn) then
        ! Debugging
        OPEN(UNIT=12, FILE="Jn_CAR.dat", ACTION="write", STATUS="replace")
        DO i=1,npts_CAR
          WRITE(12,*) Jn_CAR(i)
        END DO
        CLOSE(12)

        OPEN(UNIT=12, FILE="I_elec_J_complete.dat", ACTION="write", STATUS="replace")
        DO i=1,sr%n_electrode_J
          WRITE(12,*) I_elec_J_complete(i)
        END DO
        CLOSE(12)

        OPEN(UNIT=12, FILE="I_elec_w_complete.dat", ACTION="write", STATUS="replace")
        DO i=1,sr%n_electrode_w
          WRITE(12,*) I_elec_w_complete(i)
        END DO
        CLOSE(12)
      end if


    case ('gauss')

      J_normal_mat(:,1) = J_normal_rescale(:)
  
      ! Transform the assembled array in Bezier
      call fourier( J_normal_mat, nctrl, n_plane, nsub_ctrl, 1, phi , mapping, Jn_bez_mat, my_id , bnd_node_list, bnd_elm_list)

      Jn_bez(:) = Jn_bez_mat(:,1)

      !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      !===================  N.B. We use the CARIDDI-computed matrices in JOREK space to map to the electrode surfaces ===========
      ! =========================================================================================================================
      
      ! Shell - reference normal orientation
      do i =1,sr%n_electrode_J
        I_elec_J_complete(i) = sum(sr%Jmat_J(i,:)*Jn_bez(:))
      end do

      ! Wall - opposite normal orientation
      do i =1,sr%n_electrode_w
        I_elec_w_complete(i) = -sum(sr%Jmat_w(i,:)*Jn_bez(:))
      end do


      if (debug_Jn) then 
        if (my_id.eq.0) then 
          OPEN(UNIT=12, FILE="I_elec_J_complete.dat", ACTION="write", STATUS="replace")
          DO i=1,sr%n_electrode_J
            WRITE(12,*) I_elec_J_complete(i)
          END DO
          CLOSE(12)

          OPEN(UNIT=12, FILE="I_elec_w_complete.dat", ACTION="write", STATUS="replace")
          DO i=1,sr%n_electrode_w
            WRITE(12,*) I_elec_w_complete(i)
          END DO
          CLOSE(12)
        endif 
      endif
      
      ! Enforce solenoidality of electrodes currents
      if (enforce_kirch) then 
        I_elec_J_complete(:) = I_elec_J_complete(:) - sum(I_elec_J_complete) / sr%n_electrode_J
        I_elec_w_complete(:) = I_elec_w_complete(:) - sum(I_elec_w_complete) / sr%n_electrode_w
        
        if (my_id.eq.0) then 
          if (debug_Jn) then 
            OPEN(UNIT=12, FILE="I_elec_J_complete_Kirchhoff.dat", ACTION="write", STATUS="replace")
            DO i=1,sr%n_electrode_J
              WRITE(12,*) I_elec_J_complete(i)
            END DO
            CLOSE(12)

            OPEN(UNIT=12, FILE="I_elec_w_complete_Kirchhoff.dat", ACTION="write", STATUS="replace")
            DO i=1,sr%n_electrode_w
              WRITE(12,*) I_elec_w_complete(i)
            END DO
            CLOSE(12)
          endif
        endif
      endif
      
      ! CLEANUP AXISYMMETRIC COMPONENT FROM HALO CURRENTS
      if (remove_axi_halo) then ! For now, only valid for an all-plastered shell, where n_electrodes=n_bnd_elements
        
        ! Take chunks of size n_tor
        n_tor_elec_J = 40
        n_tor_elec_w = 40
        n_pol_elec_J = sr%n_electrode_J/n_tor_elec_J 
        n_pol_elec_w = sr%n_electrode_w/n_tor_elec_w

        if (.not. allocated(i_elec_tor_J)) allocate(i_elec_tor_J(n_tor_elec_J))
        if (.not. allocated(i_elec_tor_w)) allocate(i_elec_tor_w(n_tor_elec_w))
        if (.not. allocated(I_elec_J_new)) allocate(I_elec_J_new(sr%n_electrode_J))
        i_elec_J_new(:) = 0.d0
        if (.not. allocated(I_elec_w_new)) allocate(I_elec_w_new(sr%n_electrode_w))
        i_elec_w_new(:) = 0.d0
        
        ! SHELL
        ! Identify a poloidal location
        do i = 1,n_pol_elec_J
          ! Initialize array of toroidal currents to 0
          i_elec_tor_j(:) = 0.d0

          do j =1,n_tor_elec_J 
            i_elec_tor_J(j) = I_elec_J_complete((j-1) * n_pol_elec_J + i) 
          enddo 

          ! Subtract by mean value
          i_elec_tor_j = i_elec_tor_J - (  sum(I_elec_tor_J) / real(n_tor_elec_J) )

          ! Put the values back in the main array 
          do j=1,n_tor_elec_J
            I_elec_J_new ((j-1) * n_pol_elec_J + i) = i_elec_tor_J(j)
          enddo 

        enddo  
        
        ! WALL
        do i = 1,n_pol_elec_w
          ! Initialize array of toroidal currents to 0
          i_elec_tor_w(:) = 0.d0

          do j =1,n_tor_elec_w 
            i_elec_tor_w(j) = I_elec_w_complete((j-1) * n_pol_elec_w + i) 
          enddo 

          ! Subtract by mean value
          i_elec_tor_w = i_elec_tor_w - (  sum(I_elec_tor_w) / real(n_tor_elec_w) )

          ! Put the values back in the main array 
          do j=1,n_tor_elec_w
            I_elec_w_new ((j-1) * n_pol_elec_w + i) = i_elec_tor_w(j)
          enddo 

        enddo  

        I_elec_J_complete = I_elec_J_new
        I_elec_w_complete = I_elec_w_new

        if (my_id.eq.0) then 
          if (debug_Jn) then 
            OPEN(UNIT=12, FILE="I_elec_J_complete_WO_AXI.dat", ACTION="write", STATUS="replace")
            DO i=1,sr%n_electrode_J
              WRITE(12,*) I_elec_J_complete(i)
            END DO
            CLOSE(12)

            OPEN(UNIT=12, FILE="I_elec_w_complete_WO_AXI.dat", ACTION="write", STATUS="replace")
            DO i=1,sr%n_electrode_w
              WRITE(12,*) I_elec_w_complete(i)
            END DO
            CLOSE(12)
          endif
        endif

      endif

      Ielec_J = I_elec_J_complete(1:sr%n_electrode_J-1)
      Ielec_w = I_elec_w_complete(1:sr%n_electrode_w-1)

      if (debug_Jn) then

        !Debugging
        if (my_id == 0) then
          OPEN(UNIT=12, FILE="Jn_gauss_from_Ielec_halo.dat", ACTION="write", STATUS="replace")
          DO i=1,nctrl
            WRITE(12,*) J_normal(i)
          END DO
          CLOSE(12)
        endif

        if (my_id == 0) then
          OPEN(UNIT=12, FILE="Jn_bez_from_Ielec_halo.dat", ACTION="write", STATUS="replace")
          DO i=1,n_dof_bnd*sr%n_tor
            WRITE(12,*) Jn_bez(i)
          END DO
          CLOSE(12)
        endif
        

      endif

      if (reconvert_to_real) then

        L_PHI:do mp = 1, n_plane

          ! phi value
          phi_v = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)

          L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements  ! Elements loop
            bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
            L_MP: do m_pt = 1, nsub_ctrl  ! Poloidal points in an element
              index = npol*(mp-1)+(m_bndelem-1)*nsub_ctrl + m_pt
        
              select case (mapping)
                case ('equidistant') 
                  s = float(m_pt-1)/float(nsub_ctrl)           
                  !call basisfunctions1(s, H1, H1_s, H1_ss)
                end select
                !step_dist       = sr%a_phiye%step ! for distributed matrices       
                L_IV: do i_vertex = 1, 2 ! (Loop over nodes in element m_bndelem)
                  i_start = 1

                  L_LS:  do l_tor = 1, sr%n_tor ! STARWALL harmonics
                    f_tor    = sr%i_tor(l_tor) ! Index for STARWALL harmonics
                    !if ( (f_tor < i_tor_min) .or. (f_tor > i_tor_max) ) cycle 
                    i_node      = bndelem_m%vertex(i_vertex)
                    i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)
                    size_i = bndelem_m%size(i_vertex,:)
                    !write(*,*) ' H1 :', H1(i_vertex,i_dof)

                    i_resp_1   = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor &
                           + bnd_node_list%bnd_node(i_node_bnd)%n_dof*(i_start-1) &
                           + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) &
                           - bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1

                    i_resp_2   = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor &
                           + bnd_node_list%bnd_node(i_node_bnd)%n_dof*(i_start-1) &
                           + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(2) &
                           - bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1

                    Bez2Real(index, i_resp_1 ) = Bez2Real (index, i_resp_1) + H1(i_vertex,1, m_pt) * size_i(1) &
                                       * cos(mode_n(l_tor)*phi_v)
                    Bez2Real(index, i_resp_2 ) = Bez2Real (index, i_resp_2) + H1(i_vertex,2, m_pt) * size_i(2) &
                                       * cos(mode_n(l_tor)*phi_v)

                    if (mode_n(l_tor) .ne. 0) then
                      i_resp_1   = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor &
                           + bnd_node_list%bnd_node(i_node_bnd)%n_dof*(i_start) &
                           + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) &
                           - bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1

                      i_resp_2   = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor &
                           + bnd_node_list%bnd_node(i_node_bnd)%n_dof*(i_start) &
                           + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(2) &
                           - bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1


                      Bez2Real(index, i_resp_1 ) = Bez2Real (index, i_resp_1) + H1(i_vertex,1, m_pt) * size_i(1) &
                                         * sin(mode_n(l_tor)*phi_v)
                      Bez2Real(index, i_resp_2 ) = Bez2Real (index, i_resp_2) + H1(i_vertex,2, m_pt) * size_i(2) &
                                         * sin(mode_n(l_tor)*phi_v)
                    endif

                    if (mode_n(l_tor) .eq. 0) then
                      i_start = i_start + 1
                    else
                      i_start = i_start + 2
                    end if

                    ! if (i_resp>=sr%a_ey_pol%ind_start .AND. i_resp<=sr%a_ey_pol%ind_end) then
                    ! is_c = is_c_JOREK(l_tor)
                    !  itor = itor_fourier(l_tor)
              
                    ! if (is_c) then
                    ! else              
                    ! endif

                  enddo L_LS
                enddo L_IV
              enddo L_MP
            enddo L_MB
          enddo L_PHI

        do i = 1,nctrl
          Jn_back2real(i) = sum( Bez2Real(i,:) * Jn_bez(:) )
        enddo

        if (my_id == 0) then
          OPEN(UNIT=12, FILE="Jn_back2real.dat", ACTION="write", STATUS="replace")
          DO i=1,nctrl
            WRITE(12,*) Jn_back2real(i)
          END DO
          CLOSE(12)
        endif

      endif ! reconvert2real
  end select
  if (my_id .eq. 0) write(*,*) 'Finished get_Ielec_halo'
  end subroutine get_Ielec_halo
  
  subroutine halo_time_stepping(my_id, node_list, element_list, bnd_elm_list, bnd_node_list, response_R_zj, response_R_psis, debug, update_electrode_current, dzjbnd_vec, dpsibnd_vec)
  
  use data_structure
  use gauss
  use mod_settings
  use phys_module, only:  PI, n_period, F0, index_now, index_start, nstep_n
  use basis_at_gaussian, only: H1, H1_s, HZ
  use constants, only: MU_ZERO
  use mod_basisfunctions

  integer,                      intent(in)    :: my_id                !< MPI process ID
  type (type_node_list),        intent(in)    :: node_list
  type (type_element_list),     intent(in)    :: element_list
  type (type_bnd_element_list), intent(in)    :: bnd_elm_list   
  type (type_bnd_node_list),    intent(in)    :: bnd_node_list
  real*8,                       intent(out)   :: response_R_zj(:,:)
  real*8,                       intent(out)   :: response_R_psis(:,:)
  logical,                      intent(in)    :: debug
  logical,                      intent(in)    :: update_electrode_current              ! Update electrode currents?
  real*8, OPTIONAL,             intent(in)    :: dzjbnd_vec(:)
  real*8, OPTIONAL,             intent(in)    :: dpsibnd_vec(:)

  !> LOCAL VARIABLES
  real*8              :: deltaJn(sr%nd_bez), deltaIelec_full(sr%n_electrode_w), Ielec_w_pre_full(sr%n_electrode_w), Ielec_J_pre_full(sr%n_electrode_J), Ielec_w_pre(sr%n_electrode_w-1), Ielec_J_pre(sr%n_electrode_J-1)
  real*8              :: Ielec_w_now_full(sr%n_electrode_w), Ielec_J_now_full(sr%n_electrode_J), Ielec_w_now(sr%n_electrode_w-1), Ielec_J_now(sr%n_electrode_J-1)
  integer             :: Npts_tot, npol,jj,mp, ms, i, j
  real*8, allocatable :: fact_R_zj(:), fact_R_psis(:), R_zj(:,:), R_psis(:,:), J_normal_past(:), J_normal_pres(:), jn_bez_past(:), jn_bez_pres(:), test_sig(:), bez_sig(:), phi_t
  real*8, allocatable :: identity_mat(:,:), real2bez(:,:), fact_R_zj_point, fact_R_psis_point
  real*8   :: psi_s(n_gauss*n_plane), zj(n_gauss*n_plane), R(n_gauss), Z(n_gauss), R_s(n_gauss), Z_s(n_gauss),vec_inside(2), grad_t(2), G(4, n_degrees)
  real*8   :: psi_s_array_pres(n_gauss*n_plane) , zj_array_pres(n_gauss*n_plane)
  real*8   :: phi_vec(n_plane), R_c, Z_c
  type(type_bnd_element) :: bndelem_s, bndelem_e, bndelem_l
  type(type_element)     :: elem_s    , elem_l
  type(type_node)        :: node_k, node_s, node_l
  integer             :: s_bndelem, s_vertex, s_node, s_dof, s_dir, s_tor, k, n_steps
  integer             :: l_bndelem, l_vertex, l_node, l_dof, l_dir, l_tor
  real*8   :: s_size, ss_size, j_size, fact_mu0, sign_out, e_size, psi_prec, zj_prec
  real*8   :: l_size, psi_pres, zj_pres
  integer             :: index_start_realsp, index_end_realsp, real_entry_index, s_resp, s_harm
  integer             :: e_bndelem, e_vertex, e_node, e_node_bnd, s_node_bnd, e_dof, e_index, e_dir, e_tor, e_harm, e_resp


  ! Writing transformation matrix real-space -> Bezier
  Npts_tot = n_gauss * n_plane * bnd_elm_list%n_bnd_elements

  if (.not. ALLOCATED(R_psis)) allocate(R_psis(npts_tot, sr%nd_bez))
  if (.not. ALLOCATED(R_zj)) allocate(R_zj(npts_tot, sr%nd_bez ))
  if (.not. ALLOCATED(fact_R_psis)) allocate(fact_R_psis(npts_tot))
  if (.not. ALLOCATED(fact_R_zj)) allocate(fact_R_zj(npts_tot))
  if (.not. ALLOCATED(J_normal_past)) allocate(J_normal_past(npts_tot))
  if (.not. ALLOCATED(J_normal_pres)) allocate(J_normal_pres(npts_tot))
  if (.not. ALLOCATED(jn_bez_past)) allocate(jn_bez_past(sr%nd_bez))
  if (.not. ALLOCATED(jn_bez_pres)) allocate(jn_bez_pres(sr%nd_bez))
  if (.not. ALLOCATED(test_sig)) ALLOCATE(test_sig(npts_tot))
  if (.not. ALLOCATED(bez_sig)) ALLOCATE(bez_sig(sr%nd_bez))
      
  npol     = bnd_elm_list%n_bnd_elements * n_gauss
  if (my_id.eq. 0) write(*,*) 'get_delta_Ielec: NPOL =', npol
  if (.not. allocated(identity_mat)) allocate(identity_mat(Npts_tot, Npts_tot))

  identity_mat = 0.d0
  do jj = 1 , Npts_tot
    identity_mat(jj,jj) = 1.d0
  end do
  phi_vec = 0.d0
  do mp = 1, n_plane
    phi_vec   (mp)      = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
  enddo
  ! N.B. Real2Bez WILL BE ALLOCATED BASED ON THE NUMBER OF FREE BOUNDARY MODES IN STARWALL!    
  call fourier(identity_mat, Npts_tot, n_plane, n_gauss, Npts_tot, phi_vec, 'gauss', Real2Bez, my_id, bnd_node_list, bnd_elm_list)
  
  if (update_electrode_current) then 

    ! Calculate real space coefficients for implicit time stepping at previous time step
    fact_R_psis = 0.d0
    fact_R_zj   = 0.d0 
    J_normal_past   = 0.d0
    J_normal_pres   = 0.d0
    jn_bez_past      = 0.d0
    jn_bez_pres      = 0.d0
    fact_mu0 = 1.d0/mu_zero
    R_zj = 0.d0 
    R_psis = 0.d0

    do s_bndelem = 1,bnd_elm_list%n_bnd_elements

      bndelem_s = bnd_elm_list%bnd_element(s_bndelem)
      elem_s    = element_list%element(bndelem_s%element)
      psi_s = 0.d0
      zj    = 0.d0
      R     = 0.d0
      Z     = 0.d0
      R_s   = 0.d0
      Z_s   = 0.d0

      do s_vertex = 1, 2 ! (loop over nodes in element s_bndelem)
        
        do s_dof = 1, 2 ! (loop over node dofs)
          s_node      = bndelem_s%vertex(s_vertex)
          node_s      = node_list%node(s_node)
          s_dir       = bndelem_s%direction(s_vertex,s_dof)
          s_size      = bndelem_s%size(s_vertex,s_dof)
          
          R  (:)   = R  (:)  + node_s%x(1,s_dir,1) * s_size * H1  (s_vertex,s_dof,:)
          Z  (:)   = Z  (:)  + node_s%x(1,s_dir,2) * s_size * H1  (s_vertex,s_dof,:)
          R_s(:)   = R_s(:)  + node_s%x(1,s_dir,1) * s_size * H1_s(s_vertex,s_dof,:)
          Z_s(:)   = Z_s(:)  + node_s%x(1,s_dir,2) * s_size * H1_s(s_vertex,s_dof,:)

          do mp=1,n_plane
            do s_tor = 1, n_tor ! (loop over toroidal harmonics)
              ! Get values at previous t step
              psi_prec = node_s%values(s_tor,s_dir,var_psi) - node_s%deltas(s_tor,s_dir,var_psi)
              zj_prec =  node_s%values(s_tor,s_dir,var_zj)  - node_s%deltas(s_tor,s_dir,var_zj)
              ! Chunks of size n_gauss * n_plane
              psi_s((mp-1)*n_gauss+1:mp*n_gauss) =    psi_s((mp-1)*n_gauss+1:mp*n_gauss) + &
                                                      psi_prec * s_size*H1_s(s_vertex,s_dof,:) * HZ(s_tor,mp)

              zj((mp-1)*n_gauss+1:mp*n_gauss)    =    zj((mp-1)*n_gauss+1:mp*n_gauss) + &
                                                      zj_prec * s_size*H1(s_vertex,s_dof,:) * HZ(s_tor,mp)
            enddo            
          enddo

        enddo
      enddo

      !---------------------- BOUNDARY ELEMENT ORIENTATION-------------------------- 
      ! --- Calculate an inside point on the element to then calculate the
      ! direction of bnd normals
      vec_inside(2) = 0.d0
      grad_t(2)     = 0.d0
      sign_out   = 0.d0

      call basisfunctions(xgauss(2),xgauss(2), G)  
      R_c = 0.d0 ;  Z_c = 0.d0 
      do i = 1, n_vertex_max
        do j = 1, n_degrees
          node_k = node_list%node(elem_s%vertex(i)) 
          R_c    = R_c + node_k%x(1,j,1) * elem_s%size(i,j) * G(i,j)
          Z_c    = Z_c + node_k%x(1,j,2) * elem_s%size(i,j) * G(i,j)
        enddo
      enddo  

      vec_inside = (/ R_c - R(2), Z_c - Z(2) /)       ! vector pointing towards the domain
      !--- Find out correct sign of the normal (it has to point outwards the domain)
      grad_t     = (/ -Z_s(2) , R_s(2) /)     ! gradient of the coordinate t (normal to the boundary here) 
      ! grad_t could point either outwards or inwards depending on the orientation of the bnd_element!
      sign_out   = -1.d0 * sign( 1.d0, ( vec_inside(1)*grad_t(1) + vec_inside(2)*grad_t(2) ) )  ! (+1 if normal outwards)
      
      do s_vertex = 1, 2 ! (loop over nodes in element s_bndelem)
        s_node_bnd  = bndelem_s%bnd_vertex(s_vertex)
        do s_dof = 1, 2 ! (loop over node dofs)
            ss_size = bndelem_s%size(s_vertex,s_dof)
          do ms = 1,n_gauss
            do mp=1,n_plane

              ! Active loops:
              ! s_bndelem
              ! s_vertex 
              ! s_dof
              ! ms
              ! mp

              ! We are sitting on a localized Gauss point in the whole mesh
              fact_R_zj_point   = &
                      psi_s ((mp-1)*n_gauss+ms) / sqrt(R_s(ms)**2.d0 + Z_s(ms)**2.d0) / (F0*R(ms)) * sign_out   

              fact_R_psis_point   = &
                  zj( (mp-1)*n_gauss+ms ) / sqrt(R_s(ms)**2.d0 + Z_s(ms)**2.d0) / (F0*R(ms)) * sign_out             

              do s_harm = 1, sr%n_tor
                s_tor  = sr%i_tor(s_harm)

                ! Real space index
                real_entry_index = n_gauss*(s_bndelem-1) + ms + (mp-1)*npol
                ! Bezier space index
                s_resp   = (bnd_node_list%bnd_node(s_node_bnd)%index_starwall(1) - 1)*sr%n_tor0 &
                              + bnd_node_list%bnd_node(s_node_bnd)%n_dof*(s_harm-1)              &
                              + bnd_node_list%bnd_node(s_node_bnd)%index_starwall(s_dof)             &
                              - bnd_node_list%bnd_node(s_node_bnd)%index_starwall(1) + 1
              
                R_zj( real_entry_index, s_resp ) = &
                        fact_mu0*fact_R_zj_point &
                      * ss_size * H1(s_vertex,s_dof,ms) * HZ(s_tor,mp)                  

                R_psis( real_entry_index, s_resp ) = &
                        fact_mu0*fact_R_psis_point &
                      * ss_size * H1_s(s_vertex,s_dof,ms) * HZ(s_tor,mp)

              enddo  ! s_harm
            enddo  ! mp
          enddo  ! ms
        enddo  ! s_dof
      enddo  ! s_vertex
      
      do mp=1,n_plane 

        index_start_realsp =  (s_bndelem-1)*n_gauss+1+(mp-1)*npol
        index_end_realsp   =   s_bndelem*n_gauss+(mp-1)*npol
        ! Compute the reference J_normal at the previous step
        J_normal_past (index_start_realsp : index_end_realsp) = &
            ( ( zj((mp-1)*n_gauss+1:mp*n_gauss) * psi_s((mp-1)*n_gauss+1:mp*n_gauss)) ) &
                    / sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out

      enddo
    enddo   ! s_bndelem

    ! ============ Loop for calculation of Jnormal in the present (n) =================

    do l_bndelem = 1,bnd_elm_list%n_bnd_elements

      bndelem_l = bnd_elm_list%bnd_element(l_bndelem)
      elem_l    = element_list%element(bndelem_l%element)
      psi_s_array_pres = 0.d0
      zj_array_pres    = 0.d0
      R     = 0.d0
      Z     = 0.d0
      R_s   = 0.d0
      Z_s   = 0.d0

      do l_vertex = 1, 2 ! (loop over nodes in element l_bndelem)
        do l_dof = 1, 2 ! (loop over node dofs)

          l_node      = bndelem_l%vertex(l_vertex)
          node_l      = node_list%node(l_node)
          l_dir       = bndelem_l%direction(l_vertex,l_dof)
          l_size      = bndelem_l%size(l_vertex,l_dof)

          R  (:)   = R  (:)  + node_l%x(1,l_dir,1) * l_size * H1  (l_vertex,l_dof,:)
          Z  (:)   = Z  (:)  + node_l%x(1,l_dir,2) * l_size * H1  (l_vertex,l_dof,:)
          R_s(:)   = R_s(:)  + node_l%x(1,l_dir,1) * l_size * H1_s(l_vertex,l_dof,:)
          Z_s(:)   = Z_s(:)  + node_l%x(1,l_dir,2) * l_size * H1_s(l_vertex,l_dof,:)

          do mp=1,n_plane
            do l_tor = 1, n_tor ! (loop over toroidal harmonics)
              ! Get values at previous t step
              psi_pres = node_l%values(l_tor,l_dir,var_psi) 
              zj_pres =  node_l%values(l_tor,l_dir,var_zj)  
              ! Chunks of size n_gauss * n_plane
              psi_s_array_pres((mp-1)*n_gauss+1:mp*n_gauss) =    psi_s_array_pres((mp-1)*n_gauss+1:mp*n_gauss) + &
                                                      psi_pres * l_size*H1_s(l_vertex,l_dof,:) * HZ(l_tor,mp)

              zj_array_pres((mp-1)*n_gauss+1:mp*n_gauss)    =    zj_array_pres((mp-1)*n_gauss+1:mp*n_gauss) + &
                                                      zj_pres * l_size*H1(l_vertex,l_dof,:) * HZ(l_tor,mp)
            enddo            
          enddo

        enddo
      enddo
      
      !---------------------- BOUNDARY ELEMENT ORIENTATION-------------------------- 
      ! --- Calculate an inside point on the element to then calculate the
      ! direction of bnd normals
      vec_inside(2) = 0.d0
      grad_t(2)     = 0.d0
      sign_out   = 0.d0

      call basisfunctions(xgauss(2),xgauss(2), G)  
      R_c = 0.d0 ;  Z_c = 0.d0 
      do i = 1, n_vertex_max
        do j = 1, n_degrees
          node_k = node_list%node(elem_l%vertex(i)) 
          R_c    = R_c + node_k%x(1,j,1) * elem_l%size(i,j) * G(i,j)
          Z_c    = Z_c + node_k%x(1,j,2) * elem_l%size(i,j) * G(i,j)
        enddo
      enddo  

      vec_inside = (/ R_c - R(2), Z_c - Z(2) /)       ! vector pointing towards the domain
      !--- Find out correct sign of the normal (it has to point outwards the domain)
      grad_t     = (/ -Z_s(2) , R_s(2) /)     ! gradient of the coordinate t (normal to the boundary here) 
      ! grad_t could point either outwards or inwards depending on the orientation of the bnd_element!
      sign_out   = -1.d0 * sign( 1.d0, ( vec_inside(1)*grad_t(1) + vec_inside(2)*grad_t(2) ) )  ! (+1 if normal outwards)

      do mp=1,n_plane 

        index_start_realsp =  (l_bndelem-1)*n_gauss+1+(mp-1)*npol
        index_end_realsp   =   l_bndelem*n_gauss+(mp-1)*npol
        ! Compute the reference J_normal
        J_normal_pres (index_start_realsp : index_end_realsp) = &
            ( ( zj_array_pres((mp-1)*n_gauss+1:mp*n_gauss) * psi_s_array_pres((mp-1)*n_gauss+1:mp*n_gauss)) ) &
                    / sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out

      enddo
    enddo   ! l_bndelem
    
    ! Convert the matrices to Bezier space
    response_R_zj = 0.d0 
    response_R_psis = 0.d0
    response_R_zj = MATMUL(real2bez,R_zj) 
    response_R_psis = MATMUL(real2bez, R_psis )
    
    ! RESCALE J_normal for SI units! 
    J_normal_past = J_normal_past * fact_mu0
    J_normal_pres = J_normal_pres * fact_mu0
    ! Calculate J_normal in Bezier
    do jj = 1,sr%nd_bez
      Jn_bez_past(jj) = sum(real2bez(jj,:) * J_normal_past(:))
      Jn_bez_pres(jj) = sum(real2bez(jj,:) * J_normal_pres(:))
    enddo
    
    Ielec_w_pre_full(:) = 0.d0 
    Ielec_J_pre_full(:) = 0.d0
    Ielec_w_pre(:) = 0.d0 
    Ielec_J_pre(:) = 0.d0 

    ! Compute I_elec from physics (at the previous time step!)
    do i =1,sr%n_electrode_w
      Ielec_w_pre_full(i) = -sum(sr%Jmat_w(i,:)*Jn_bez_past(:))
      Ielec_w_now_full(i) = -sum(sr%Jmat_w(i,:)*Jn_bez_pres(:))
    end do
    Ielec_w_pre(:) = Ielec_w_pre_full(1:sr%n_electrode_w-1)
    Ielec_w_now(:) = Ielec_w_now_full(1:sr%n_electrode_w-1)

    do i =1,sr%n_electrode_J
      Ielec_J_pre_full(i) = sum(sr%Jmat_J(i,:)*Jn_bez_past(:))
      Ielec_J_now_full(i) = sum(sr%Jmat_J(i,:)*Jn_bez_pres(:))
    end do
    Ielec_J_pre(:) = Ielec_J_pre_full(1:sr%n_electrode_J-1)
    Ielec_J_now(:) = Ielec_J_now_full(1:sr%n_electrode_J-1)

    n_steps = sum(nstep_n)
    
    

    !!!!!! ALLOCATE ELECTRODE CURRENT MATRICES FOR POSTPROC !!!!!!!!!!!!!!!!!!!!!!

    !!! Ielecw_phys_mat
    if (.not. allocated (Ielecw_phys_mat_past)) then 
        write(*,*) 'Allocating Ielecw_phys_mat_past with', sr%n_electrode_w-1 , n_steps
        allocate( Ielecw_phys_mat_past(sr%n_electrode_w-1 , n_steps) )
        Ielecw_phys_mat_past(:,:) = 0.d0
    end if  
    if (.not. allocated (Ielecw_phys_mat_pres)) then 
        write(*,*) 'Allocating Ielecw_phys_mat_pres with', sr%n_electrode_w-1 , n_steps
        allocate( Ielecw_phys_mat_pres(sr%n_electrode_w-1 , n_steps) )
        Ielecw_phys_mat_pres(:,:) = 0.d0
    end if 
    !!! IelecJ_phys_mat
    if (.not. allocated (IelecJ_phys_mat_past)) then 
        write(*,*) 'Allocating Ielecj_phys_mat_past with', sr%n_electrode_j-1 , n_steps
        allocate( IelecJ_phys_mat_past(sr%n_electrode_j-1 , n_steps) )
        IelecJ_phys_mat_past(:,:) = 0.d0
    end if  
    if (.not. allocated (IelecJ_phys_mat_pres)) then 
        write(*,*) 'Allocating Ielecj_phys_mat_pres with', sr%n_electrode_j-1 , n_steps
        allocate( IelecJ_phys_mat_pres(sr%n_electrode_j-1 , n_steps) )
        IelecJ_phys_mat_pres(:,:) = 0.d0
    end if   

    !!! Ielecw_phys_mat_full 
    if (.not. allocated (Ielecw_phys_mat_past_full)) then 
      write(*,*) 'Allocating Ielecw_phys_mat_past_full with', sr%n_electrode_w, n_steps
      allocate( Ielecw_phys_mat_past_full(sr%n_electrode_w , n_steps) )
      Ielecw_phys_mat_past_full(:,:) = 0.d0
    end if

    if (.not. allocated (Ielecw_phys_mat_pres_full)) then 
      write(*,*) 'Allocating Ielecw_phys_mat_pres_full with', sr%n_electrode_w, n_steps
      allocate( Ielecw_phys_mat_pres_full(sr%n_electrode_w , n_steps) )
      Ielecw_phys_mat_pres_full(:,:) = 0.d0
    end if

    !!! IelecJ_phys_mat_full 
    if (.not. allocated (IelecJ_phys_mat_past_full)) then 
        write(*,*) 'Allocating Ielecj_phys_mat_past_full with', sr%n_electrode_j, n_steps
        allocate( IelecJ_phys_mat_past_full(sr%n_electrode_j , n_steps) )
        IelecJ_phys_mat_past_full(:,:) = 0.d0
    end if  
    if (.not. allocated (IelecJ_phys_mat_pres_full)) then 
        write(*,*) 'Allocating Ielecj_phys_mat_pres_full with', sr%n_electrode_j, n_steps
        allocate( IelecJ_phys_mat_pres_full(sr%n_electrode_j , n_steps) )
        IelecJ_phys_mat_pres_full(:,:) = 0.d0
    end if  

    if (index_now >= 2) then
      ! I_elec_pre are calculated at the previous time step
      Ielecw_phys_mat_past(:,index_now-1) = Ielec_w_pre(:) 
      IelecJ_phys_mat_past(:,index_now-1) = Ielec_J_pre(:)
      Ielecw_phys_mat_past_full(:,index_now-1) = Ielec_w_pre_full(:)
      IelecJ_phys_mat_past_full(:,index_now-1) = Ielec_J_pre_full(:)
    endif
    ! I_elec_now are calculated at current time step
    Ielecw_phys_mat_pres(:,index_now) = Ielec_w_now(:) 
    IelecJ_phys_mat_pres(:,index_now) = Ielec_J_now(:)
    Ielecw_phys_mat_pres_full(:,index_now) = Ielec_w_now_full(:)
    IelecJ_phys_mat_pres_full(:,index_now) = Ielec_J_now_full(:)

    if (my_id .eq. 0) then 

      open(unit=10, file="Ielecw_phys_mat_past.txt", status="replace")
      do i = 1, n_wall_curr_e
        write(10,'(*(ES24.16,1X))') (Ielecw_phys_mat_past(i,j), j = 1, n_steps)
      end do
      close(10)

      open(unit=10, file="IelecJ_phys_mat_past.txt", status="replace")
      do i = 1, n_shell_curr_e
        write(10,'(*(ES24.16,1X))') (IelecJ_phys_mat_past(i,j), j = 1, n_steps)
      end do
      close(10)

      open(unit=10, file="Ielecw_phys_mat_pres.txt", status="replace")
      do i = 1, n_wall_curr_e
        write(10,'(*(ES24.16,1X))') (Ielecw_phys_mat_pres(i,j), j = 1, n_steps)
      end do
      close(10)

      open(unit=10, file="IelecJ_phys_mat_pres.txt", status="replace")
      do i = 1, n_shell_curr_e
        write(10,'(*(ES24.16,1X))') (IelecJ_phys_mat_pres(i,j), j = 1, n_steps)
      end do
      close(10)

      open(unit=10, file="Ielecw_phys_mat_past_full.txt", status="replace")
      do i = 1, sr%n_electrode_w
        write(10,'(*(ES24.16,1X))') (Ielecw_phys_mat_past_full(i,j), j = 1, n_steps)
      end do
      close(10)

      open(unit=10, file="IelecJ_phys_mat_past_full.txt", status="replace")
      do i = 1, sr%n_electrode_J
        write(10,'(*(ES24.16,1X))') (IelecJ_phys_mat_past_full(i,j), j = 1, n_steps)
      end do
      close(10)

      open(unit=10, file="Ielecw_phys_mat_pres_full.txt", status="replace")
      do i = 1, sr%n_electrode_w
        write(10,'(*(ES24.16,1X))') (Ielecw_phys_mat_pres_full(i,j), j = 1, n_steps)
      end do
      close(10)

      open(unit=10, file="IelecJ_phys_mat_pres_full.txt", status="replace")
      do i = 1, sr%n_electrode_J
        write(10,'(*(ES24.16,1X))') (IelecJ_phys_mat_pres_full(i,j), j = 1, n_steps)
      end do
      close(10)

    endif

    if (debug) then 
      if (my_id .eq. 0) write(*,*) 'HALO TIME STEPPING: Writing J_normal_halo_time_stepping.dat'
      open(unit=12, file="J_normal_halo_time_stepping.dat", status="replace", action="write", form="formatted")
      do i = 1,npts_tot
        write(12,*) J_normal_pres(i)
      end do
      close(12)
      if (my_id .eq. 0) write(*,*) 'Done writing J_normal.dat'   
      
    
      if (my_id .eq. 0) write(*,*) 'HALO TIME STEPPING: Writing J_bez_halo_time_stepping.dat'
      open(unit=12, file="J_bez_halo_time_stepping.dat", status="replace", action="write", form="formatted")
      do i = 1,n_dof_starwall
        write(12,*) jn_bez_pres(i)
      end do
      close(12)
      if (my_id .eq. 0) write(*,*) 'Done writing Jn_bez.dat' 

    endif
    
    ! --- COMPUTE DELTAS

    if (present(dzjbnd_vec) .and. present(dpsibnd_vec)) then 

      ! if (index_now <= 1) then 
      !   ! Ielec_w(:) = Ielec_w_pre(:) 
      !   ! Ielec_J(:) = Ielec_J_pre(:)
      !   dIelec_J = 0.d0 
      !   dIelec_w = 0.d0 
      ! else 
      deltaJn(:) = 0.d0
      ! Get delta Jn
      do i =1 ,sr%nd_bez
        deltaJn(i) = sum(response_R_zj(i,:) * dzjbnd_vec(:)) + sum(response_R_psis(i,:) * dpsibnd_vec(:))
      enddo
      deltaIelec_full(:) = 0.d0
      ! Get delta I_elec wall at step (n-1) by using conversion matrix from CARIDDI
      do i =1,sr%n_electrode_w
      deltaIelec_full(i) = -sum(sr%Jmat_w(i,:)*deltaJn(:))
      end do
      dIelec_w(:) = deltaIelec_full(1:sr%n_electrode_w-1)

      deltaIelec_full(:) = 0.d0
      ! Get delta I_elec shell at step (n-1) by using conversion matrix
      do i =1,sr%n_electrode_J
      deltaIelec_full(i) = sum(sr%Jmat_J(i,:)*deltaJn(:))
      end do
      dIelec_J(:) = deltaIelec_full(1:sr%n_electrode_J-1)
      !endif
      
      if (index_now .eq. 2) then 

        OPEN(UNIT=12, FILE="I_elec_w_pre_2_stepping.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) Ielec_w(k)
        END DO

        OPEN(UNIT=12, FILE="delta_I_elec_w_1_stepping.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) dIelec_w(k)
        END DO

        OPEN(UNIT=12, FILE="I_elec_w_now_1.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) Ielec_w_pre(k)
        END DO

        OPEN(UNIT=12, FILE="I_elec_J_now_1.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_J-1
          WRITE(12,*) Ielec_J_pre(k)
        END DO

      endif 

      if (index_now .eq. 3) then 

        OPEN(UNIT=12, FILE="I_elec_w_pre_3_stepping.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) Ielec_w(k)
        END DO

        OPEN(UNIT=12, FILE="delta_I_elec_w_2_stepping.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) dIelec_w(k)
        END DO

        OPEN(UNIT=12, FILE="I_elec_w_now_2.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) Ielec_w_pre(k)
        END DO

        OPEN(UNIT=12, FILE="I_elec_J_now_2.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_J-1
          WRITE(12,*) Ielec_J_pre(k)
        END DO

      endif 
      
      if (index_now .eq. 4) then 
        OPEN(UNIT=12, FILE="I_elec_w_now_3.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) Ielec_w_pre(k)
        END DO

        OPEN(UNIT=12, FILE="I_elec_J_now_3.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_J-1
          WRITE(12,*) Ielec_J_pre(k)
        END DO
      endif

      ! ---- UPDATE THE PRESENT ELECTRODE CURRENTS
      ! --- Initialization from the reference
      if (index_now .eq. 1) then
        Ielec_w(:) = Ielec_w_now(:)
        Ielec_J(:) = Ielec_J_now(:)
      endif

      Ielec_w_old(:) = Ielec_w(:) 
      Ielec_J_old(:) = Ielec_J(:)
      ! WALL
      Ielec_w(:) = Ielec_w(:) + dIelec_w(:)
      ! SHELL
      Ielec_J(:) = Ielec_J(:) + dIelec_J(:)

      if (index_now .eq. 2) then 
        OPEN(UNIT=12, FILE="I_elec_w_post_2_stepping.dat", ACTION="write", STATUS="replace" )
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) Ielec_w(k)
        END DO
      endif

      if (index_now .eq. 3) then 
        OPEN(UNIT=12, FILE="I_elec_w_post_3_stepping.dat", ACTION="write", STATUS="replace" )
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) Ielec_w(k)
        END DO
      endif

      !  instant n : I(n-1) --> I(n)  ALL DELTAS ARE DELTA(n-1)

      ! M_CAR * sum(response_R_zj(i,:) * dzjbnd_vec(:)) + sum(response_R_psis(i,:) * dpsibnd_vec(:)) =

      ! Ielec_w (t2) - Ielec_w (t1)

    endif

  else  ! Field BCs: We do not want electrode currents to evolve! Just calculate the jacobian matrices in the present 

    ! Calculate real space coefficients for implicit time stepping at (n-1)
    fact_R_psis = 0.d0
    fact_R_zj   = 0.d0 
    fact_mu0 = 1.d0/mu_zero
    R_zj = 0.d0 
    R_psis = 0.d0

    do s_bndelem = 1,bnd_elm_list%n_bnd_elements

      bndelem_s = bnd_elm_list%bnd_element(s_bndelem)
      elem_s    = element_list%element(bndelem_s%element)
      psi_s = 0.d0
      zj    = 0.d0
      R     = 0.d0
      Z     = 0.d0
      R_s   = 0.d0
      Z_s   = 0.d0

      do s_vertex = 1, 2 ! (loop over nodes in element s_bndelem)
        
        do s_dof = 1, 2 ! (loop over node dofs)
          s_node      = bndelem_s%vertex(s_vertex)
          node_s      = node_list%node(s_node)
          s_dir       = bndelem_s%direction(s_vertex,s_dof)
          s_size      = bndelem_s%size(s_vertex,s_dof)
          
          R  (:)   = R  (:)  + node_s%x(1,s_dir,1) * s_size * H1  (s_vertex,s_dof,:)
          Z  (:)   = Z  (:)  + node_s%x(1,s_dir,2) * s_size * H1  (s_vertex,s_dof,:)
          R_s(:)   = R_s(:)  + node_s%x(1,s_dir,1) * s_size * H1_s(s_vertex,s_dof,:)
          Z_s(:)   = Z_s(:)  + node_s%x(1,s_dir,2) * s_size * H1_s(s_vertex,s_dof,:)

          do mp=1,n_plane
            do s_tor = 1, n_tor ! (loop over toroidal harmonics)
              ! Get values at previous t step
              psi_pres = node_s%values(s_tor,s_dir,var_psi)
              zj_pres =  node_s%values(s_tor,s_dir,var_zj) 
              ! Chunks of size n_gauss * n_plane
              psi_s((mp-1)*n_gauss+1:mp*n_gauss) =    psi_s((mp-1)*n_gauss+1:mp*n_gauss) + &
                                                      psi_pres * s_size*H1_s(s_vertex,s_dof,:) * HZ(s_tor,mp)

              zj((mp-1)*n_gauss+1:mp*n_gauss)    =    zj((mp-1)*n_gauss+1:mp*n_gauss) + &
                                                      zj_pres * s_size*H1(s_vertex,s_dof,:) * HZ(s_tor,mp)
            enddo            
          enddo

        enddo
      enddo

      !---------------------- BOUNDARY ELEMENT ORIENTATION-------------------------- 
      ! --- Calculate an inside point on the element to then calculate the
      ! direction of bnd normals
      vec_inside(2) = 0.d0
      grad_t(2)     = 0.d0
      sign_out   = 0.d0

      call basisfunctions(xgauss(2),xgauss(2), G)  
      R_c = 0.d0 ;  Z_c = 0.d0 
      do i = 1, n_vertex_max
        do j = 1, n_degrees
          node_k = node_list%node(elem_s%vertex(i)) 
          R_c    = R_c + node_k%x(1,j,1) * elem_s%size(i,j) * G(i,j)
          Z_c    = Z_c + node_k%x(1,j,2) * elem_s%size(i,j) * G(i,j)
        enddo
      enddo  

      vec_inside = (/ R_c - R(2), Z_c - Z(2) /)       ! vector pointing towards the domain
      !--- Find out correct sign of the normal (it has to point outwards the domain)
      grad_t     = (/ -Z_s(2) , R_s(2) /)     ! gradient of the coordinate t (normal to the boundary here) 
      ! grad_t could point either outwards or inwards depending on the orientation of the bnd_element!
      sign_out   = -1.d0 * sign( 1.d0, ( vec_inside(1)*grad_t(1) + vec_inside(2)*grad_t(2) ) )  ! (+1 if normal outwards)
      
      do s_vertex = 1, 2 ! (loop over nodes in element s_bndelem)
        s_node_bnd  = bndelem_s%bnd_vertex(s_vertex)
        do s_dof = 1, 2 ! (loop over node dofs)
            ss_size = bndelem_s%size(s_vertex,s_dof)
          do ms = 1,n_gauss
            do mp=1,n_plane

              ! Active loops:
              ! s_bndelem
              ! s_vertex 
              ! s_dof
              ! ms
              ! mp

              ! We are sitting on a localized Gauss point in the whole mesh
              fact_R_zj_point   = &
                      psi_s ((mp-1)*n_gauss+ms) / sqrt(R_s(ms)**2.d0 + Z_s(ms)**2.d0) / (F0*R(ms)) * sign_out   

              fact_R_psis_point   = &
                  zj( (mp-1)*n_gauss+ms ) / sqrt(R_s(ms)**2.d0 + Z_s(ms)**2.d0) / (F0*R(ms)) * sign_out             

              do s_harm = 1, sr%n_tor

                s_tor  = sr%i_tor(s_harm)
                ! Real space index
                real_entry_index = n_gauss*(s_bndelem-1) + ms + (mp-1)*npol
                ! Bezier space index
                s_resp   = (bnd_node_list%bnd_node(s_node_bnd)%index_starwall(1) - 1)*sr%n_tor0 &
                              + bnd_node_list%bnd_node(s_node_bnd)%n_dof*(s_harm-1)              &
                              + bnd_node_list%bnd_node(s_node_bnd)%index_starwall(s_dof)             &
                              - bnd_node_list%bnd_node(s_node_bnd)%index_starwall(1) + 1
              
                R_zj( real_entry_index, s_resp ) = &
                        fact_mu0*fact_R_zj_point &
                      * ss_size * H1(s_vertex,s_dof,ms) * HZ(s_tor,mp)                  

                R_psis( real_entry_index, s_resp ) = &
                        fact_mu0*fact_R_psis_point &
                      * ss_size * H1_s(s_vertex,s_dof,ms) * HZ(s_tor,mp)

              enddo  ! s_harm
            enddo  ! mp
          enddo  ! ms
        enddo  ! s_dof
      enddo  ! s_vertex
    enddo
    
    ! Convert the matrices to Bezier space
    response_R_zj = 0.d0 
    response_R_psis = 0.d0
    response_R_zj = MATMUL(real2bez,R_zj) 
    response_R_psis = MATMUL(real2bez, R_psis )
  
  endif


  end subroutine halo_time_stepping

 
  subroutine fourier(A, nctrl, ntor_f, nsub_ctrl, ndim, phi_ctrl,mapping, global ,my_id, node_list, element_list)
  use data_structure
    real*8 :: A(nctrl,ndim)                               !< Function to fourier decompose
    real*8, allocatable  :: global(:,:)                 !< output matrix with Bezier
    type (type_bnd_node_list),        intent(in)    :: node_list
    type (type_bnd_element_list),     intent(in)    :: element_list
    integer, intent(in)  :: nctrl, ndim, ntor_f, nsub_ctrl !< number of poloidal points, dim, toroidal points
    real*8, intent(in) :: phi_ctrl(:) 
    character(len=*),                 intent(in)    :: mapping
    character*2 :: rank
    integer :: i, j, inode, i_v, i_order, iv
    integer :: ind, ind0, ind_glob
    integer :: my_id
    integer :: n_dof_tot, npol
    real*8, allocatable, dimension(:,:) :: fourier_coeff, bez

   ! n_dof_bnd = node_list%n_dof_bnd
    n_dof_tot = n_dof_bnd*sr%n_tor
    npol = nctrl/ntor_f
    
    allocate(fourier_coeff(npol,ndim), bez(n_dof_bnd, ndim))
    if (allocated(global)) deallocate(global)
    allocate(global(n_dof_tot,ndim))
    if ( my_id .eq. 0 ) write(*,*) 'transform into realspace and bezier'

    do i = 1, sr%n_tor
       call calc_fourier_coeff_single(A,  ntor_f, npol, i, ndim, phi_ctrl, fourier_coeff,my_id)
       call realspace2bezier(bez, fourier_coeff, ndim, nsub_ctrl, mapping, node_list, element_list, my_id)
       do inode = 1, node_list%n_bnd_nodes
          do i_order=1,2
             ind  = global_JOREK_index(node_list,inode, i_order, i, sr%n_tor)
             ind0 = JOREK_index_1mode(node_list,inode, i_order)
             global(ind,:) = bez(ind0,:)
          end do
       end do
    end do
    deallocate(fourier_coeff, bez)
    if ( my_id .eq. 0 ) write(*,*) 'end of fourier routine'

  end subroutine fourier

  !< Calculate toroidal Fourier coefficients of function
  subroutine calc_fourier_coeff_single(f, ntor_f, npol, itor_J ,ndim, phi_ctrl, fourier, my_id)

    real*8,  allocatable  ::  fourier(:,:)      !< Fourier coefficient, f(npol, ndim)
    real*8,  intent(in)   ::  f(:,:)            !< function to decompose
                                                !< size( npol*ntor_f:ndim )
                                                !< Ordering first in poloidal direction then in toroidal
    integer, intent(in)   ::  npol,  ntor_f     !< # of poloidal field points, toroidal points
    integer, intent(in)   ::  itor_J            !< harmonic in JOREK convention
    integer, intent(in)   ::  ndim              !< second dimension of input function
    real*8,  intent(in)   ::  phi_ctrl(:)       !< phi of control surface (in realspace convention)
                                                !< dimension = ntor_f
    integer, intent(in)   ::  my_id           !< MPI rank
    
    !<----------------Local variables ----------------------------
    integer :: itor      !< number of harmonic in realspace convention
    integer :: jphi, k,i, i_start, i_end, mode_n(sr%n_tor)
    real*8  :: phi       
    logical :: is_c      !< T = cos, F = sin
    
    if (.not. allocated(fourier)) allocate(fourier(npol,ndim))

    !< Get the right mode number

    mode_n = mode_number(sr%i_tor)
    ! if (my_id .eq. 0) then
    ! open(unit=12, file="mode_n.dat", action='write', status='replace')
    ! do i = 1, sr%n_tor
    !   write(12,*) mode_n(i)
    ! enddo 
    ! close(12)
    ! endif

    !< calculate real harmonics from JOREK input
    !if (my_id .eq. 0 ) write(*,*) 'itor_J = ', itor_J
    itor = itor_fourier(itor_J)
    !if (my_id .eq. 0 ) write(*,*) 'itor = ', itor 
    is_c = is_c_JOREK(itor_J)

    fourier=0.d0
    
    !write(*,*) ' Fourier subroutine: itor_J = ', itor_J, ', itor = ', itor, ', is_c = ', is_c, ', nharm(itor) = ', nharm(itor)
 
    ! Calculate Fourier series coefficients
    do jphi = 1, ntor_f   ! toroidal slices 
       phi     =  phi_ctrl(jphi)
       i_start = (jphi-1)*npol + 1
       i_end   =  jphi   *npol
       if(is_c) then
          fourier(:,:)   = fourier(:,:) + 2*cos(phi * mode_n(itor) )/ntor_f * f(i_start:i_end,:)
       else
          ! minus from the JOREK convention, phi has to be given in realspace convention
          fourier(:,:)   = fourier(:,:) + 2*sin(phi * mode_n(itor) )/ntor_f * f(i_start:i_end,:)
      end if
   end do

   if (mode_n(itor)==0) fourier=fourier/2.d0

  end subroutine calc_fourier_coeff_single

  subroutine calc_fourier_coeff_complex(f, ntor_f, npol, i_mode ,ndim, phi_ctrl, fourier_c, my_id)

  complex*16,allocatable::  fourier_c(:,:)    !< Exp. Fourier coefficient, f(npol, ndim)
  real*8,  intent(in)   ::  f(:,:)            !< function to decompose
                                              !< size( npol*ntor_f:ndim )
                                              !< Ordering first in poloidal direction then in toroidal
  integer, intent(in)   ::  npol,  ntor_f     !< # of poloidal field points, toroidal points
  integer, intent(in)   ::  i_mode            !< mode number
  integer, intent(in)   ::  ndim              !< second dimension of input function
  real*8,  intent(in)   ::  phi_ctrl(:)       !< phi of control surface (in realspace convention)
                                              !< dimension = ntor_f
  integer, intent(in)   ::  my_id             !< MPI rank
  
  !<----------------Local variables ----------------------------
  integer :: itor      !< number of harmonic in realspace convention
  integer :: jphi, k,i, i_start, i_end, mode_n(sr%n_tor)
  real*8  :: phi       
  real*8,allocatable :: fourier_a(:,:), fourier_b(:,:)
  complex*16 :: i_cplx, factor_compl
  COMPLEX*16, allocatable :: fourier_c_ref(:,:)
  
  if (.not. allocated(fourier_c)) allocate(fourier_c(npol,ndim))
  if (.not. allocated(fourier_c_ref)) allocate(fourier_c_ref(npol,ndim))
  if (.not. allocated(fourier_a)) allocate(fourier_a(npol,ndim))
  if (.not. allocated(fourier_b)) allocate(fourier_b(npol,ndim))

  fourier_c=(0.d0,0.d0)
  fourier_c_ref = (0.d0,0.d0)
  i_cplx   =(0.d0,1.d0)
  fourier_a = 0.d0
  fourier_b = 0.d0

  ! Calculate Fourier series coefficients
  do jphi = 1, ntor_f   ! toroidal slices 
    phi     =  phi_ctrl(jphi)
    i_start = (jphi-1)*npol + 1
    i_end   =  jphi   *npol


    factor_compl = exp(-i_cplx * i_mode * phi  ) / ntor_f
    fourier_c(:,:) = fourier_c(:,:) + factor_compl * f(i_start:i_end, :)

    ! For check, compare with the trigonometric series
    fourier_a(:,:) = fourier_a(:,:) + 2.d0*cos(phi * i_mode )/ntor_f * f(i_start:i_end,:)
    fourier_b(:,:) = fourier_b(:,:) - 2*sin(phi * i_mode )/ntor_f * f(i_start:i_end,:)

  end do
  if (i_mode .eq. 0) then
    fourier_a(:,:) = fourier_a(:,:) / 2.d0
    fourier_b(:,:) = 0.d0
  endif
  ! Print arrays
  ! if (i_mode .eq. 0) then 
  !   write(*,*) 'size fourier_a=', size(fourier_a,DIM=1)
  !   write(*,*) 'size fourier_c=', size(fourier_c,DIM=1)

  !   ! a0=c0
  !   if (my_id .eq. 0) then 
  !     open(unit=12, file='a0.dat', action='write', status='replace')
  !     do i = 1,size(fourier_a,DIM=1)
  !       write(12,*) fourier_a(i,1) 
  !     enddo 
  !   endif 

  !   if (my_id .eq. 0) then 
  !     open(unit=12, file='c0.dat', action='write', status='replace')
  !     do i = 1,size(fourier_c,DIM=1)
  !       write(12,*) real(fourier_c(i,1) )
  !     enddo 
  !   endif

  ! else  

  !   ! Calculate comparison ref.
  !   fourier_c_ref = (fourier_a -  (-i_cplx * fourier_b))/2.d0

  !   if (my_id .eq. 0) then 
  !     open(unit=12, file='cn_ref.dat', action='write', status='replace')
  !     do i = 1,size(fourier_c_ref,DIM=1)
  !       write(12,*) real(fourier_c_ref(i,1) ) , AIMAG(fourier_c_ref(i,1))
  !     enddo 
  !   endif

  !   if (my_id .eq. 0) then 
  !     open(unit=12, file='cn.dat', action='write', status='replace')
  !     do i = 1,size(fourier_c,DIM=1)
  !       write(12,*) real(fourier_c(i,1) ) , AIMAG(fourier_c(i,1))
  !     enddo 
  !   endif

  ! endif
    

end subroutine calc_fourier_coeff_complex



  subroutine realspace2bezier(bez, f, ndim, nsub, mapping, node_list, element_list, my_id)
  use data_structure
  use mod_basisfunctions
  use gauss,          only: n_gauss, xgauss, wgauss
    integer, intent(in)  :: nsub        !< Number of evaluation points to be used for quadrature
    integer, intent(in)  :: ndim        !< Second Dimension of function f
    real*8 , intent(in)  :: f(:,:)      !< function (Bfield) to be projected on Bezier basis (nsub*nelements:ndim)
    character(len=*),             intent(in)    :: mapping  
    integer, intent(in)   ::  my_id           !< MPI rank
    real*8 , intent(out), allocatable :: bez(:,:)    !< Bezier coefficients
    type (type_bnd_node_list),        intent(in)    :: node_list
    type (type_bnd_element_list),     intent(in)    :: element_list

    
    real*8, dimension(4,4)                     :: LSQ            !> Coefficients of LHS matrix (s. wiki)
    real*8, dimension(n_dof_bnd,n_dof_bnd)     :: aa, inv        !> matrix of LHS and its inverse
    ! real*8, dimension(node_list%n_dof_bnd,   ndim)       :: integral       !> integrated function f of RHS
    real*8, allocatable       :: integral(:,:)       !> integrated function f of RHS
    real*8, dimension(2,2)                     :: H, Hs, Hss   !> basis functions
    real*8, allocatable, dimension(:,:,:)              :: H_g            !> arrays of basis functions at gaussian points
    real*8, allocatable, dimension(:)                  :: x, w           !> Evaluation points and Gaussian weights
    real*8                                     :: s, sz(2,2)     !> coordinate in element, size of element
    type(type_bnd_element)                     :: elm            !> element
    integer                                    :: ielm, i1, i2, i3, i4, i, j, k, ii, dd, gg
    integer                                    :: nd, info
    integer                                    :: n_points_poloidal
    integer, allocatable, dimension(:)         :: subseq
    
    if (allocated(bez)) deallocate(bez)
    allocate(bez(n_dof_bnd,ndim))
    allocate(integral(n_dof_bnd,ndim))
    
    LSQ = reshape([  52., 22., 18.,  13., &          !>(https://www.jorek.eu/wiki/doku.php?id=real_space2bezier.f90)       
                 22., 12., 13.,   9., &
                 18., 13., 52.,  22., &
                 13.,  9., 22.,  12.] / 140., shape(LSQ))
    
    !< Vector used for the trapezoidal rule
    n_points_poloidal = element_list%n_bnd_elements * nsub
    allocate(subseq(n_points_poloidal))
    subseq(1:n_points_poloidal-1) = (/ (ii, ii = 1,n_points_poloidal-1,1) /)
    subseq(n_points_poloidal) = 0
    subseq = subseq +1
    
    select case (mapping)
     case ('equidistant')
        if (my_id .eq. 0) write(*,*) 'realspace2bezier: Using equistant points'
        if (ALLOCATED(x)) deallocate(x)
        if (ALLOCATED(w)) deallocate(w)
        if (ALLOCATED(H_g)) deallocate(H_g)
        allocate(x(nsub+1))
        allocate(w(nsub))
        allocate(H_g(nsub+1,2,2))
        w = 1.d0 / nsub
      do j=1,nsub+1
        x(j) = (j-1.d0)/nsub
        call  basisfunctions1((x(j)), H, Hs, Hss)
        H_g(j,:,:) = H                        ! basis function at Gaussian points       
      end do
    case ('gauss')
        if (my_id .eq. 0) write(*,*) 'realspace2bezier: Using Gauss points'
        if (ALLOCATED(x)) deallocate(x)
        if (ALLOCATED(w)) deallocate(w)
        if (ALLOCATED(H_g)) deallocate(H_g)
        allocate(x(nsub))
        allocate(w(nsub))
        allocate(H_g(nsub,2,2))
        x = Xgauss
        w = Wgauss
        do j=1,nsub
          call  basisfunctions1((x(j)), H, Hs, Hss)
          H_g(j,:,:) = H                        ! basis function at Gaussian points
        enddo
    end select

    aa = 0.d0
    integral = 0.d0
    do ielm=1,element_list%n_bnd_elements  ! N_bnd: number of boundary elements from JOREK
      elm = element_list%bnd_element(ielm)
      i1 = node_list%bnd_node(elm%bnd_vertex(1))%index_starwall(1)
      i2 = node_list%bnd_node(elm%bnd_vertex(1))%index_starwall(2)
      i3 = node_list%bnd_node(elm%bnd_vertex(2))%index_starwall(1)
      i4 = node_list%bnd_node(elm%bnd_vertex(2))%index_starwall(2)
      sz = elm%size
      
  
      ! ------- Calculate analytical integral of Bezier coefficients ------------------
      aa(i1,i1) = aa(i1,i1) + LSQ(1,1) * sz(1,1) * sz(1,1)
      aa(i1,i2) = aa(i1,i2) + LSQ(1,2) * sz(1,2) * sz(1,1)
      aa(i1,i3) = aa(i1,i3) + LSQ(1,3) * sz(2,1) * sz(1,1)
      aa(i1,i4) = aa(i1,i4) + LSQ(1,4) * sz(2,2) * sz(1,1)

      aa(i2,i1) = aa(i2,i1) + LSQ(2,1) * sz(1,1) * sz(1,2)
      aa(i2,i2) = aa(i2,i2) + LSQ(2,2) * sz(1,2) * sz(1,2)
      aa(i2,i3) = aa(i2,i3) + LSQ(2,3) * sz(2,1) * sz(1,2)
      aa(i2,i4) = aa(i2,i4) + LSQ(2,4) * sz(2,2) * sz(1,2)

      aa(i3,i1) = aa(i3,i1) + LSQ(3,1) * sz(1,1) * sz(2,1)
      aa(i3,i2) = aa(i3,i2) + LSQ(3,2) * sz(1,2) * sz(2,1)
      aa(i3,i3) = aa(i3,i3) + LSQ(3,3) * sz(2,1) * sz(2,1)
      aa(i3,i4) = aa(i3,i4) + LSQ(3,4) * sz(2,2) * sz(2,1)

      aa(i4,i1) = aa(i4,i1) + LSQ(4,1) * sz(1,1) * sz(2,2)
      aa(i4,i2) = aa(i4,i2) + LSQ(4,2) * sz(1,2) * sz(2,2)
      aa(i4,i3) = aa(i4,i3) + LSQ(4,3) * sz(2,1) * sz(2,2)
      aa(i4,i4) = aa(i4,i4) + LSQ(4,4) * sz(2,2) * sz(2,2)


      !--------------- Interpolate function using Gauss points ------------------------
      ! write(*,*) 'In element loop'
      ! write(6,*) 'In element loop'
    select case (mapping)
    case ('equidistant')
     do k = 1, nsub
        j = k + nsub*(ielm-1)
        integral(i1, :) =  integral(i1, :) + (H_g(k,1,1)*f(j,:) + H_g(k+1,1,1)*f(subseq(j),:))/2 * sz(1,1)* w(k)
        integral(i2, :) =  integral(i2, :) + (H_g(k,1,2)*f(j,:) + H_g(k+1,1,2)*f(subseq(j),:))/2 * sz(1,2)* w(k) 
        integral(i3, :) =  integral(i3, :) + (H_g(k,2,1)*f(j,:) + H_g(k+1,2,1)*f(subseq(j),:))/2 * sz(2,1)* w(k) 
        integral(i4, :) =  integral(i4, :) + (H_g(k,2,2)*f(j,:) + H_g(k+1,2,2)*f(subseq(j),:))/2 * sz(2,2)* w(k) 
      enddo
      ! write(*,*) 'In element loop, after nsub loop'
      ! write(6,*) 'In element loop, after nsub loop'
    case ('gauss')
      do k = 1, nsub
        j = k + nsub*(ielm-1)
        integral(i1, :) =  integral(i1, :) + H_g(k,1,1)*f(j,:) * sz(1,1)* w(k)
        integral(i2, :) =  integral(i2, :) + H_g(k,1,2)*f(j,:)  * sz(1,2)* w(k) 
        integral(i3, :) =  integral(i3, :) + H_g(k,2,1)*f(j,:)  * sz(2,1)* w(k) 
        integral(i4, :) =  integral(i4, :) + H_g(k,2,2)*f(j,:)  * sz(2,2)* w(k) 
      enddo
    end select

  enddo
    ! write(*,*) 'After element loop'
    ! write(6,*) 'After element loop'

    !     do dd = 1,n_dof_bnd
    !   do gg = 1,n_dof_bnd
    !     write(*,*) 'realspace2bez debug', aa(dd,gg)
    !   enddo
    ! enddo
    inv = 0.
    do  i = 1, n_dof_bnd
      inv(i,i) = 1.d0
    enddo
    nd  = n_dof_bnd

    call dpotrf('U',nd,aa,nd,info)
    call dpotrs('U',nd,nd,aa,nd,inv,nd,info)


    ! ----------------------------------------------------------------------
    bez = 0.

    do i=1, n_dof_bnd
      do k=1, n_dof_bnd
        bez(i,:) = bez(i,:) + inv(i,k)*integral(k,:)
      enddo
    enddo
    deallocate(integral)
  end subroutine realspace2bezier


  ! --- Routine to determine if which fourier mode corresponds to JOREK
  integer function itor_fourier(itor_J)
    integer, intent(in) :: itor_J      ! --- mode in JOREK counting
    integer             :: mode_n(sr%n_tor),i
    
    mode_n = mode_number(sr%i_tor)

    if (mode_n(1) .eq. 0 ) then
      itor_fourier = itor_J/2 + 1   
    else
      itor_fourier = (itor_J+1)/2
    end if


  end function itor_fourier

    ! --- Routine to determine if JOREK mode is sin or cos
  logical function is_c_JOREK(itor_J)
    integer, intent(in) :: itor_J      ! --- mode in JOREK counting
    integer             :: mode_n(sr%n_tor),i

    mode_n = mode_number(sr%i_tor)

    if (mode_n(1) .eq. 0 ) then
      is_c_JOREK = mod(itor_J,2)==0
      if (itor_J==1) is_c_JOREK = .true.
    else
      is_c_JOREK = mod(itor_J+1,2) == 0
    end if
  end function is_c_JOREK
  
  ! Called just for n=0
  subroutine fourier2realspace(fourier, n_modes, nctrl, npol, ntor_f, real_func)
  use constants
  real*8 :: fourier(npol, n_modes)
  integer :: npol, ntor_f, dim, nctrl, n_modes
  real*8 :: phi, real_func(nctrl, n_modes)
  integer :: i, j, k
  integer :: mode_n(sr%n_tor)
  mode_n = mode_number(sr%i_tor)
    
  do i = 1, npol
    do k = 1, ntor_f
      phi = float(k-1)/ntor_f*2*pi
      do j = 1, n_modes
        if (is_c_JOREK(j)) then
          real_func ((k-1)*npol+1:k*npol,j) = real_func((k-1)*npol+1:k*npol,j) + fourier(:,j) * cos(phi*mode_n(j))
        else
         real_func ((k-1)*npol+1:k*npol,j) = real_func((k-1)*npol+1:k*npol,j) - fourier(:,j) * sin(phi*mode_n(j))
        end if
      end do     
    end do
      
  end do
  end subroutine fourier2realspace

    ! --- Routine to calculate the global JOREK index
  pure integer function global_JOREK_index(node_list,i_v, i_dof, i_start, n_modes)
    use data_structure
    type (type_bnd_node_list),        intent(in)    :: node_list
    integer, intent(in) :: i_v      ! --- Index of the vertex in the node_list
    integer, intent(in) :: i_dof    ! --- Local dof of the node 
    integer, intent(in) :: i_start  ! --- Mode counter
    integer, intent(in) :: n_modes  ! --- Number of modes in JOREK (counting sine and cosine separately)

    global_JOREK_index = (  node_list%bnd_node(i_v)%index_starwall(1)-1)*n_modes        &   
        + node_list%bnd_node(i_v)%n_dof*(i_start-1)                &
        - node_list%bnd_node(i_v)%index_starwall(1)                   & 
        + node_list%bnd_node(i_v)%index_starwall(i_dof) + 1 

  end function global_JOREK_index
  ! --- Routine to calculate the global JOREK index
  pure integer function JOREK_index_1mode(node_list,i_v, i_dof)
    use data_structure
    type (type_bnd_node_list),        intent(in)    :: node_list
    integer, intent(in) :: i_v      ! --- Index of the vertex in the node_list
    integer, intent(in) :: i_dof    ! --- Local dof of the node

    JOREK_index_1mode  =  node_list%bnd_node(i_v)%index_starwall(i_dof)

  end function JOREK_index_1mode


! --- This function adapts the index in the sr%i_tor array (which, for e.g. n_tor in the STARWALL input equal to 0,1,2
! will contain indexes 1,2,3,4,5 where 1 is the axisymmetric one and the others are cosine and sine for each harmonic)
! -------> and CONVERTS it into the actual mode number n
  pure function mode_number(JOREK_harmonics) result (mode_n)
  implicit none
  integer, DIMENSION(:), intent(in) :: JOREK_harmonics
  integer, DIMENSION(size(JOREK_harmonics)) :: mode_n
  integer             :: i
  
  do i = 1, sr%n_tor

    if (JOREK_harmonics(i) .eq. 1 ) mode_n(i) = 0
    
    if (mod(JOREK_harmonics(i),2) .eq. 0) then
      mode_n(i) = JOREK_harmonics(i)/2  ! Harmonic number is even
    else
      mode_n(i) = (JOREK_harmonics(i) -1)/2
    endif
  enddo

  end function mode_number

  subroutine find_loc_indexes(array, sub_array, loc_indexes) 
  ! Finds the indexes of the second array's elements occurring in the first array
  implicit none
  integer,               intent(in) :: array(:)
  integer,               intent(in) :: sub_array(:)
  integer, ALLOCATABLE,  intent(out):: loc_indexes(:)
  integer                           :: i

  ! Initialize result with all possible starting indices
  loc_indexes = [(i, i = 1, size(array)-size(sub_array)+1)]
  ! Sequentially keep only those indices that satisfy a condition
  do i = 1, size(sub_array)
    ! Condition for keeping: the value of array(results + i - 1) must be equal to the value of sub_array(i)
    ! Use PACK to only keep results that satisfy this condition
    loc_indexes = PACK(loc_indexes, array(loc_indexes + i - 1) == sub_array(i))
    if (size(loc_indexes) == 0) then
      exit
    end if
  end do
  end subroutine find_loc_indexes

 end module vacuum
