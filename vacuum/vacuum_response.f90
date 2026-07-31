!> Implements the interplay of the plasma with a conducting wall.
!!
!! The plasma-wall interaction is characterized by vacuum response matrices which are calculated by
!! the STARWALL code and imported into JOREK by the routine read_starwall_response(). The matrices
!! allow to express the magnetic field parallel to the interface (boundary of the JOREK domain) in
!! terms of the poloidal flux at the interface and the wall currents. The vacuum
!! response enters into the boundary integral of the current equation which vanishes for fixed
!! boundary conditions. The boundary integral is implemented in the routine
!! vacuum_boundary_integral().
!!
!! @note The variable s in a boundary element may correspond to s or t of the respective
!! 2D element depending on element orientation.
module vacuum_response
  
  use vacuum
  
  implicit none
  
  
  
  contains
  
  
  
  !> Is the filehandle associated with a formatted file?
  logical function is_formatted(filehandle)
    implicit none
    ! --- Routine parameters
    integer,          intent(in) :: filehandle
    ! --- Local variables
    character(len=64) :: format_type
    
    inquire(filehandle,form=format_type)
    is_formatted = ( trim(format_type) == 'FORMATTED' )
  end function is_formatted
  
  
  
  !> Get the vacuum response for an ideal or resistive wall.
  subroutine get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list,                    &
    freeboundary_equil, resistive_wall)
  
    use mpi_mod
    use constants
    use mod_parameters, only: n_tor
    use data_structure, only: type_node_list, type_bnd_element_list, type_bnd_node_list
    use phys_module,    only: central_mass, central_density

    implicit none    
    
    integer,                     intent(in) :: my_id              !< MPI proc ID
    type(type_node_list),        intent(in) :: node_list          !< List of boundary nodes
    type(type_bnd_element_list), intent(in) :: bnd_elm_list       !< List of boundary elements
    type(type_bnd_node_list),    intent(in) :: bnd_node_list      !< List of boundary nodes
    logical,                     intent(in) :: freeboundary_equil !< Use free boundary equilibrium?
    logical,                     intent(in) :: resistive_wall     !< Resistive or ideal wall?

    integer :: i,j, ierr, dim
    logical :: exists

    n_dof_bnd = 0

    ! --- Determine total number of boundary degrees of freedom per harmonic (skipping duplicates).
    do i=1, bnd_node_list%n_bnd_nodes
      exists = .false.
      do j=1, i-1
        if (bnd_node_list%bnd_node(i)%index_jorek .eq. bnd_node_list%bnd_node(j)%index_jorek) then
          exists  = .true.
          exit
        end if
      end do
      if (.not. exists) then
        n_dof_bnd = n_dof_bnd + bnd_node_list%bnd_node(i)%n_dof ! Number of boundary degrees of freedom per harmonic
      end if
    end do

    if (my_id == 0) write(*,'(A,i5)') '   total number of degrees of freedom on the boundary : ',n_dof_bnd
    
    ! --- Write out the boundary information for STARWALL.
    if (my_id == 0) call export_boundary(node_list, bnd_elm_list, bnd_node_list)
    
    if ( vacuum_debug  ) call log_starwall_response(my_id, sr) 
    if (my_id == 0) write(*,*) 'Calling read_SW_response'
    call read_starwall_response(my_id, sr,'starwall-response.dat',bnd_elm_list%n_bnd_elements)
    if (my_id == 0) write(*,*) 'Calling broadcast SW response'
    call broadcast_starwall_response(my_id, sr)
    if (my_id == 0) write(*,*) 'Finished broadcasting SW response'
    if( sr%ncoil .gt. 0 ) call distribute_coil_names

    ! --- Set the "wall resistivity" to be used inside JOREK (actually it is the normalized thin wall resistivity)
    if ( (sr%file_version == 1)  .and. (my_id == 0) ) then
      write(*,*) 'Remark: STARWALL response file_version==1 means that wall_resistivity is specified in the JOREK namelist file.'
      write(*,*) '        Thus, the input parameter wall_resistivity_fact is ignored (WARNING).'
    else
      if (my_id ==0) then
          write(*,*) 'Remark: STARWALL response file_version>=2 means that eta_thin_w is specified in the STARWALL input.'
          write(*,*) '        The JOREK variable wall_resistivity is automatically calculated from it.'
          write(*,*) '        Thus, the input parameter wall_resistivity is ignored (WARNING).'
      end if
      wall_resistivity = wall_resistivity_fact * sr%eta_thin_w * &
        sqrt( central_density * 1.d20 * central_mass * mass_proton / mu_zero )
    end if
   
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call log_starwall_response(my_id, sr)
    call init_vacuum_response(my_id,  freeboundary_equil)
              
  end subroutine get_vacuum_response



  !> Read an integer parameter from the STARWALL response file.
  integer function read_intparam(filehandle, parameter_name)

    implicit none

    ! --- Routine parameters
    integer,          intent(in) :: filehandle
    character(len=*), intent(in) :: parameter_name

    ! --- Local variables
    character(len=12) :: marker
    character(len=24) :: name
    integer           :: ierr

    if ( is_formatted(filehandle) ) then
      read(filehandle,'(A12,A24,I12)',iostat=ierr) marker, name, read_intparam
    else
      read(filehandle, iostat=ierr) marker, name, read_intparam
    end if

    if ( (ierr /= 0) .or. (trim(adjustl(marker)) /= '#@intparam') .or. (trim(adjustl(name)) /= trim(parameter_name)) ) then
      write(*,*) 'ERROR: Could not read parameter "', trim(parameter_name) ,'"from STARWALL response.'
      stop
    end if

    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"=",I12)',iostat=ierr) name,read_intparam

  end function read_intparam



  !> Read an array from the STARWALL respone file
  subroutine read_array(filehandle, array_name, dim, int1d, int2d, float1d, float2d)

    implicit none

    ! --- Routine parameters
    integer, intent(in) :: filehandle
    character(len=*), intent(in) :: array_name
    integer, intent(in) :: dim(2)
    integer, allocatable, optional, intent(inout)   :: int1d(:)
    integer, allocatable, optional, intent(inout)   :: int2d(:,:)
    real*8,  allocatable, optional, intent(inout)   :: float1d(:)
    real*8,  allocatable, optional, intent(inout)   :: float2d(:,:)

    ! --- Local variables
    character(len=12) :: marker
    integer           :: nd, d(2), ierr
    character(len=24) :: name, datatype, requested_type
    logical           :: error

    if ( present(int1d) .or. present(int2d) ) then
      requested_type = 'int'
    else
      requested_type = 'float'
    end if

    if ( is_formatted(filehandle) ) then
      read(filehandle,'(A12,A24,I12,A24,2I12)',iostat=ierr) marker, name, nd, datatype, d
    else
      read(filehandle,iostat=ierr) marker, name, nd, datatype, d
    end if
    marker   = adjustl(marker)
    name     = adjustl(name)
    datatype = adjustl(datatype)

    error = ( ierr /= 0 ) .or. ( trim(marker) /= '#@array' ) .or. ( trim(name) /= trim(array_name) )                     &
      .or. ( dim(1) /= d(1) ) .or. ( dim(2) /= d(2) ) .or. ( trim(datatype) /= trim(requested_type) )

    if ( error ) then
      write(*,*) 'ERROR: Could not read array ', trim(array_name), ' from STARWALL response.'
      stop
    end if

    if ( present(int1d) ) then

      if ( allocated(int1d) ) deallocate( int1d )
      allocate( int1d(dim(1)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) int1d(:)
      else
        read(filehandle) int1d(:)
      end if

    else if ( present(int2d) ) then

      if ( allocated(int2d) ) deallocate( int2d )
      allocate( int2d(dim(1),dim(2)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) int2d(:,:)
      else
        read(filehandle) int2d(:,:)
      end if

    else if ( present(float1d) ) then

      if ( allocated(float1d) ) deallocate( float1d )
      allocate( float1d(dim(1)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) float1d(:)
      else
        read(filehandle) float1d(:)
      end if

    else if ( present(float2d) ) then
      
      if ( allocated(float2d) ) deallocate( float2d )
      allocate( float2d(dim(1),dim(2)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,'(4ES24.16)') float2d(:,:)
      else
        read(filehandle) float2d(:,:)
      end if

    end if

    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"> type ",A," size ",2I7)')name, trim(datatype), d(1:nd)

  end subroutine read_array

  

  !> Read an integer parameter from a STARWALL response file (MPI I/O).
  integer function read_intparam_parallel(filehandle, parameter_name,disp)

    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer,                       intent(in)    :: filehandle
    character(len=*),              intent(in)    :: parameter_name
    integer(kind=MPI_OFFSET_KIND), intent(inout) :: disp !< Present location in the file

    ! --- Local variables
    integer, dimension (MPI_STATUS_SIZE) :: status
    character(len=12) :: marker
    character(len=24) :: name
    integer           :: err

    call MPI_FILE_READ(filehandle, marker,int(sizeof(marker),4),MPI_CHARACTER,status,err)
    call MPI_FILE_READ(filehandle, name,  int(sizeof(name),4)  ,MPI_CHARACTER,status,err)
    call MPI_FILE_READ(filehandle, read_intparam_parallel, 1,MPI_INTEGER,status,err)
    disp = disp + sizeof(marker) + sizeof(name) + sizeof(10)! any integer 

    if ( (err /= 0) .or. (trim(adjustl(marker)) /= '#@intparam') .or. (trim(adjustl(name)) /= trim(parameter_name)) ) then
      write(*,*) 'ERROR: Could not read parameter "', trim(parameter_name) ,'"from STARWALL response.'
      stop
    end if
    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"=",I12)',iostat=err) name, read_intparam_parallel

  end function read_intparam_parallel  
 

 
  
  !> Read an array from the STARWALL respone file (MPI I/O).
  subroutine read_array_not_distr(filehandle, array_name, dim, disp, int1d, int2d, float1d, float2d, char1d, deb)

    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer,                        intent(in)      :: filehandle
    character(len=*),               intent(in)      :: array_name
    integer,                        intent(in)      :: dim(2)
    integer(kind=MPI_OFFSET_KIND),  intent(inout)   :: disp !< Present location in the file
    integer, allocatable, optional, intent(inout)   :: int1d(:)
    integer, allocatable, optional, intent(inout)   :: int2d(:,:)
    real*8,  allocatable, optional, intent(inout)   :: float1d(:)
    real*8,  allocatable, optional, intent(inout)   :: float2d(:,:)
    logical,                        intent(in)      :: deb
    character(len=*), allocatable, optional, intent(inout)   :: char1d(:)

    ! --- Local variables
    integer, dimension (MPI_STATUS_SIZE) :: status
    character(len=12) :: marker
    integer           :: nd, d(2), ierr, err, str_len
    character(len=24) :: name, datatype, requested_type
    logical           :: error
    character, allocatable :: tmp_char1d(:)
    
    if ( present(int1d) .or. present(int2d) ) then
      requested_type = 'int'
    elseif (present(char1d))then
      requested_type = 'char'
    else
      requested_type = 'float'
    end if

    call MPI_FILE_READ(filehandle,  marker,   int(sizeof(marker),4),    MPI_CHARACTER,status,ierr)
    if (deb) write(*,*) 'ierr after marker', ierr
    call MPI_FILE_READ(filehandle,  name,     int(sizeof(name),4),      MPI_CHARACTER,status,ierr)
    if (deb) write(*,*) 'ierr after name', ierr
    call MPI_FILE_READ(filehandle,  nd,       1,                        MPI_INTEGER  ,status,ierr)
    if (deb) write(*,*) 'ierr after nd', ierr
    call MPI_FILE_READ(filehandle,  datatype, int(sizeof(datatype),4),  MPI_CHARACTER,status,ierr)
    if (deb) write(*,*) 'ierr after datatype', ierr
    call MPI_FILE_READ(filehandle,  d,        2,                        MPI_INTEGER  ,status,ierr)
    if (deb) write(*,*) 'ierr after d', ierr
    disp = disp + sizeof(marker) + sizeof(name) + sizeof(nd) + sizeof(datatype) + sizeof(d) 

    marker   = adjustl(marker)
    name     = adjustl(name)
    datatype = adjustl(datatype)

    if (deb) then 
      
      write(*,*) '( ierr /= 0 )', (ierr /= 0 )
      write(*,*) '( trim(marker) /= #@array )', ( trim(marker) /= '#@array' )
      write(*,*) 'marker is', trim(marker) 
      write(*,*) '( trim(name) /= trim(array_name) )' , ( trim(name) /= trim(array_name) )
      write(*,*) 'name is', trim(name) 
      write(*,*) '( dim(1) /= d(1) )', ( dim(1) /= d(1) )
      write(*,*) 'd(1) is', d(1)
      write(*,*) '( dim(2) /= d(2) )', ( dim(2) /= d(2) )
      write(*,*) 'd(2) is', d(2)
      write(*,*) '( trim(datatype) /= trim(requested_type) )', ( trim(datatype) /= trim(requested_type) )
      write(*,*) 'datatype is', trim(datatype)
    endif

    error = ( ierr /= 0 ) .or. ( trim(marker) /= '#@array' ) .or. ( trim(name) /= trim(array_name) ) &
      .or. ( dim(1) /= d(1) ) .or. ( dim(2) /= d(2) ) .or. ( trim(datatype) /= trim(requested_type) )

    if ( error ) then
      write(*,*) 'ERROR: Could not read array ', trim(array_name), ' from STARWALL response.'
      stop
    end if
 
    if ( present(int1d) ) then

      if ( allocated(int1d) ) deallocate( int1d ); allocate( int1d(dim(1)) )
      call MPI_FILE_READ(filehandle, int1d, dim(1), MPI_INTEGER  ,status,ierr)
      disp = disp + sizeof(int1d)

    else if ( present(int2d) ) then

      if ( allocated(int2d) ) deallocate( int2d ); allocate( int2d(dim(1),dim(2)) )
      call MPI_FILE_READ(filehandle, int2d, dim(1)*dim(2), MPI_INTEGER  ,status,ierr)
      disp = disp + sizeof(int2d)

    else if ( present(float1d) ) then

      if ( allocated(float1d) ) deallocate( float1d ); allocate( float1d(dim(1)) )
      call MPI_FILE_READ(filehandle, float1d, dim(1), MPI_DOUBLE_PRECISION  ,status,ierr)
      disp = disp + sizeof(float1d)

    else if ( present(float2d) ) then

      if ( allocated(float2d) ) deallocate( float2d ); allocate( float2d(dim(1),dim(2)) )
      call MPI_FILE_READ(filehandle, float2d, dim(1)*dim(2), MPI_DOUBLE_PRECISION  ,status,ierr)
      disp = disp + sizeof(float2d)

    else if (present(char1d)) then
      
      if (allocated(char1d)) deallocate(char1d)
      allocate(char1d(1))
      str_len = dim(1) * sizeof(char1d(1))
      if ( allocated(tmp_char1d) ) deallocate( tmp_char1d ); allocate( tmp_char1d(str_len))
      if ( allocated(char1d) ) deallocate( char1d ); allocate( char1d(dim(1)))
      call MPI_FILE_READ(filehandle, tmp_char1d, str_len, MPI_CHARACTER  ,status,ierr)
      char1d = transfer(tmp_char1d, char1d)
      deallocate(tmp_char1d)
      disp = disp + str_len
    end if

    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"> type ",A," size ",2I7)')name, trim(datatype), d(1:nd)

  end subroutine read_array_not_distr
  
  
  
  !> Read an array from the STARWALL respone file
  subroutine read_array_par(filehandle, array_name, dim, disp, my_id, float2d, row_wise, loc_start, loc_len, loc_row, deb)

    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer,                        intent(in)      :: filehandle
    character(len=*),               intent(in)      :: array_name
    integer,                        intent(in)      :: dim(2)
    integer(kind=MPI_OFFSET_KIND),  intent(inout)   :: disp !< Present location in the file
    integer,                        intent(in)      :: my_id
    type(t_distrib_mat),            intent(inout)   :: float2d
    logical,                        intent(in)      :: row_wise ! if  .true. -  rowwise reading; .false. - columnwise
    integer, optional,              intent(in)      :: loc_start    ! start index of partial read in direction loc_row
    integer, optional,              intent(in)      :: loc_len      ! length partial read in direction loc_row
    logical, optional,              intent(in)      :: loc_row      ! direction of partial read (has to be opposite of row_wise)
    logical, optional,              intent(in)      :: deb


    ! --- Local variables
    integer, dimension (MPI_STATUS_SIZE) :: status
    character(len=12) :: marker
    integer           :: nd,d(2),loc_sizes(2),loc_starts(2),ierr,err,ntasks
    character(len=24) :: name, datatype
    logical           :: error
    integer           :: my_subarray
    integer           :: num_read_elements, n_step_read, step_read
    integer(KIND=8)   :: local_num_elements, num_read_elements_const, i, i_ind, j_ind
    integer           :: read_start, read_len

    if  (.not. ( (present(loc_start) .eqv.  present(loc_len))  .and.  (present(loc_start) .eqv. present(loc_row))) ) then
      write(*,*) 'If you want to read a matrix partially you have to specify loc_start, loc_len, loc_row'
      stop
    end if

    if (present(loc_start)) then
      if (loc_row .eqv. row_wise) then
        write(*,*) 'partial read not implemented for this combination, loc_row, row_wise', loc_row, row_wise, array_name
        stop
      end if
      read_len=loc_len
      read_start =loc_start -1
    else
      read_start = 0
      if (row_wise) then
        read_len   =  dim(2)
      else
        read_len   =  dim(1)
      end if
    end if

    if (my_id==0) then
      call MPI_FILE_READ(filehandle, marker,   int(sizeof(marker), 4),   MPI_CHARACTER, status, ierr)
      call MPI_FILE_READ(filehandle, name,     int(sizeof(name), 4),     MPI_CHARACTER, status, ierr)
      call MPI_FILE_READ(filehandle, nd,       1,                        MPI_INTEGER,   status, ierr)
      call MPI_FILE_READ(filehandle, datatype, int(sizeof(datatype), 4), MPI_CHARACTER, status, ierr)
      call MPI_FILE_READ(filehandle, d,        2,                        MPI_INTEGER,   status, ierr)

      disp = disp + sizeof(marker) + sizeof(name) + sizeof(nd) + sizeof(datatype) + sizeof(d)
      marker   = adjustl(marker)
      name     = adjustl(name)
      datatype = adjustl(datatype)

      ! if (deb) then 
        
      !   write(*,*) '( ierr /= 0 )', (ierr /= 0 )
      !   write(*,*) '( trim(marker) /= #@array )', ( trim(marker) /= '#@array' )
      !   write(*,*) 'marker is', trim(marker) 
      !   write(*,*) '( trim(name) /= trim(array_name) )' , ( trim(name) /= trim(array_name) )
      !   write(*,*) 'name is', trim(name) 
      !   write(*,*) '( dim(1) /= d(1) )', ( dim(1) /= d(1) )
      !   write(*,*) 'd(1) is', d(1)
      !   write(*,*) '( dim(2) /= d(2) )', ( dim(2) /= d(2) )
      !   write(*,*) 'd(2) is', d(2)
      !   write(*,*) '( trim(datatype) /= trim(requested_type) )', ( trim(datatype) /= 'float' )
      !   write(*,*) 'datatype is', trim(datatype)
      ! endif

      error = ( ierr /= 0 ) .or. ( trim(marker) /= '#@array' ) .or. (trim(name) /= trim(array_name) ) &
        .or. ( dim(1) /= d(1) ) .or. ( dim(2) /= d(2) ) .or. (trim(datatype) /= 'float' )

      if ( error ) then
        write(*,*) 'ERROR: Could not read array ', trim(array_name), ' from STARWALL response.'
        stop
      end if
    end if 
          
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)  
    call MPI_BCAST(disp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if(row_wise) then ! If we read matrix rowwise
      float2d%step       = dim(1) / ntasks
      loc_starts(:)      = (/ my_id*float2d%step, read_start /)
      loc_sizes(:)       = (/ float2d%step, read_len /)
      if (my_id==ntasks-1) loc_sizes(1) = float2d%step + dim(1)-ntasks*float2d%step
      float2d%row_wise   = .true.
    else             ! If we read matrix columnwise
      float2d%step       = dim(2) / ntasks
      loc_starts(:)      = (/ read_start, my_id*float2d%step /)
      loc_sizes(:)       = (/ read_len, float2d%step /)
      if (my_id==ntasks-1) loc_sizes(2) = float2d%step + dim(2)-ntasks*float2d%step
      float2d%row_wise   = .false.
    endif

    local_num_elements = int(loc_sizes(1),8) * int(loc_sizes(2),8)
    if (row_wise) then
      call  alloc_distr(my_id, float2d, (/dim(1) , read_len /) , row_wise)
    else
      call  alloc_distr(my_id, float2d, (/read_len, dim(2)/) , row_wise)
    end if
    
    call MPI_TYPE_CREATE_SUBARRAY(2,dim,loc_sizes,loc_starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,my_subarray,ierr)
    call MPI_Type_commit(my_subarray,ierr)

    call MPI_File_set_view(filehandle, disp, MPI_DOUBLE_PRECISION, my_subarray, "native", MPI_INFO_NULL, ierr)

    ! -----------------------------------
    ! The following is a work around for an Intel MPI limitation not allowing to read more than 2 GB
    ! of data per MPI task in each call corresponding to 268435455 double precision values.
    num_read_elements_const = 250000000

    if (num_read_elements_const>local_num_elements) then
      num_read_elements = local_num_elements
      n_step_read       = 1
    else
      n_step_read       = local_num_elements/num_read_elements_const + 1
      num_read_elements = num_read_elements_const
    end if

    call MPI_AllREDUCE(MPI_IN_PLACE,n_step_read,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

    i_ind=1
    j_ind=1

    do i = 1, n_step_read
      call MPI_File_read_all(filehandle, float2d%loc_mat(i_ind, j_ind),num_read_elements, &
        MPI_DOUBLE_PRECISION, status, ierr)

      ! Calculation of starting indices of each slice in local array     
      j_ind = int(num_read_elements,8)*i/loc_sizes(1) + 1
      i_ind = int(num_read_elements,8)*i-int((j_ind-1),8)*int(loc_sizes(1),8) + 1

      ! If number of reading step is not eaqual for each MPI task.         
      if (num_read_elements_const*int((i+1),8)>local_num_elements) then
        num_read_elements = local_num_elements-i*num_read_elements_const
        if (num_read_elements <= 0) then
          num_read_elements = 0
          j_ind = 1
          i_ind = 1
        endif
      endif
    end do
    ! End of workaround for Intel MPI limitation
    ! -----------------------------------
      
    !The following line can replace the workaround in the future
    !call MPI_File_read_all(filehandle, float2d%loc_mat, loc_sizes(1)*loc_sizes(2), MPI_DOUBLE_PRECISION,status, ierr)
    disp = disp + sizeof(1.d0)*dim(1)*dim(2)

    ! Return to ordinary file view
    call MPI_File_set_view(filehandle, disp, MPI_BYTE, MPI_BYTE,"native",MPI_INFO_NULL, ierr)

    if ( vacuum_debug  .and. (my_id==0) ) write(*,'(3x,"Read: ",A24,"> type ",A," size",2I7)')name,&
      trim(datatype), d(1:nd)

  end subroutine read_array_par
  
  !===========================================================================================================
  
  !> Read the STARWALL response matrices from a single file.
  !!
  !! file_version 1: Original
  !! file_version 2: Includes eta_thin_w
  !! file_version 3: Includes additional coil information
  !! file_version 4: Includes coil names                                    
  !! file_version 5: Includes additional information on wall resolution, wall net potentials, control surface
  !! file_version 6: Coil and wall currents sign reversed to follow JOREK coordinate system
  subroutine read_starwall_response(my_id, sr, filename, n_bnd)

    use constants
    use mod_parameters, only: n_tor, n_period
    use phys_module, only : resistive_wall
    use mpi_mod

    implicit none

    ! --- Routine parameters
    integer,                   intent(in)    :: my_id
    type(t_starwall_response), intent(inout) :: sr
    character(len=*),          intent(in)    :: filename
    integer,                   intent(in)    :: n_bnd !< Number of boundary elements for consistency check

    ! --- Local variables
    integer            :: filehandle
    character(len=512) :: comment, tmp_int2str
    integer            :: i, j, i_starw, n, is_sin, err, i_tmp, ier
    real*8             :: r_tmp
    real*8, allocatable:: tmp(:)
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer(kind=MPI_OFFSET_KIND)        :: disp !< Present location in the file
    real*8             :: test_sum
    integer            :: loc_sizes(2),ntasks

    disp = 0
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, err)

    ! --- Open file
    if (my_id ==0 ) write(*,*) 'Trying to open response file...'

    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'starwall-response.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, filehandle, err)
    if ( (err == 0) .and. (my_id == 0) ) then
      call MPI_FILE_READ(filehandle, i_tmp, 1, MPI_INTEGER,          status, err)
      call MPI_FILE_READ(filehandle, r_tmp, 1, MPI_DOUBLE_PRECISION, status, err)
      disp = disp + sizeof(42) + sizeof(42.d0)
      if ( (i_tmp/=42) .or. (r_tmp/=42.d0) ) err=-42         
    end if
    call MPI_BCAST(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier )
    if (err/=0) call MPI_FILE_CLOSE(filehandle, err)

    if ( (err /= 0) .and. (my_id == 0) ) then
      write(*,*) '  ... failed.'
      write(*,*) 'ERROR: STARWALL response file (',trim(filename),') could not be opened.'
      stop
    end if
    if (my_id ==0) write(*,*) '  ... succeeded.'

    ! --- Read data from STARWALL response file
    if (my_id == 0) then

      call MPI_FILE_READ(filehandle, comment,int(sizeof(comment),4),MPI_CHARACTER,status,err)
      disp = disp + sizeof(comment)

      sr%file_version = read_intparam_parallel(filehandle, 'file_version', disp)
      if ( sr%file_version > 6 ) then
        write(*,*) 'ERROR: STARWALL response file version ', sr%file_version, ' is not supported.'
        stop
      end if

      if ( sr%file_version < 6 ) then
        write(*,*) 'WARNING: You are using an old STARWALL file version and the wall and coil currents    '
        write(*,*) '         sign do not follow the JOREK phi direction (positive means -phi direction)   '
      end if


      sr%n_bnd  = read_intparam_parallel(filehandle, 'n_bnd' , disp)
      if ( n_bnd /= sr%n_bnd ) then
        write(*,*) 'ERROR: The number of boundary elements in the STARWALL response file is different from your grid.'
        stop
      end if
      
      sr%nd_bez = read_intparam_parallel(filehandle, 'nd_bez', disp)

      if (CARIDDI_mode .and. halo_freeb) then ! .and. add version 7 if clause
        sr%n_pts_real_w = read_intparam_parallel(filehandle, 'n_pts_real_w', disp)
        sr%n_pts_real_J = read_intparam_parallel(filehandle, 'n_pts_real_J', disp)
        sr%n_electrode_w = read_intparam_parallel(filehandle, 'n_electrode_w', disp)
        sr%n_electrode_J = read_intparam_parallel(filehandle, 'n_electrode_J', disp)
      else 
        sr%n_pts_real_w = 0
        sr%n_pts_real_J = 0
        sr%n_electrode_w = 0
        sr%n_electrode_J = 0       
      end if

      if ( sr%file_version >= 5 ) then
        sr%nv       = read_intparam_parallel(filehandle, 'nv',        disp)
        sr%n_points = read_intparam_parallel(filehandle, 'n_points',  disp)
      endif

      sr%ncoil  = read_intparam_parallel(filehandle, 'ncoil' , disp)
      sr%npot_w = read_intparam_parallel(filehandle, 'npot_w', disp)
      sr%n_w    = read_intparam_parallel(filehandle, 'n_w'   , disp)
      if (CARIDDI_mode .and. halo_freeb) then
        ! Read null space dimension for wall
        sr%n_nw   =  read_intparam_parallel(filehandle, 'n_nw'   , disp)
      else
        sr%n_nw = 0
      endif
      sr%ntri_w = read_intparam_parallel(filehandle, 'ntri_w', disp)

      if ( sr%file_version >= 5 .and. (.not. CARIDDI_mode)) then
        sr%iwall    = read_intparam_parallel(filehandle, 'iwall',    disp)
        sr%nwu      = read_intparam_parallel(filehandle, 'nwu',      disp)
        sr%nwv      = read_intparam_parallel(filehandle, 'nwv' ,     disp)
        sr%mn_w     = read_intparam_parallel(filehandle, 'mn_w',     disp)
        sr%max_mn_w = read_intparam_parallel(filehandle, 'MAX_MN_W', disp)

        call read_array_not_distr(filehandle, 'm_w',         (/sr%max_mn_w,0/),  disp,  int1d=sr%m_w, deb=.false. )
        call read_array_not_distr(filehandle, 'n_w_fourier', (/sr%max_mn_w,0/),  disp,  int1d=sr%n_w_fourier, deb=.false.)
        call read_array_not_distr(filehandle, 'rc_w',        (/sr%max_mn_w,0/),  disp,  float1d=sr%rc_w, deb=.false.)
        call read_array_not_distr(filehandle, 'rs_w',        (/sr%max_mn_w,0/),  disp,  float1d=sr%rs_w, deb=.false.)
        call read_array_not_distr(filehandle, 'zc_w',        (/sr%max_mn_w,0/),  disp,  float1d=sr%zc_w, deb=.false.)
        call read_array_not_distr(filehandle, 'zs_w',        (/sr%max_mn_w,0/),  disp,  float1d=sr%zs_w, deb=.false.)
      endif

      sr%n_tor  = read_intparam_parallel(filehandle, 'n_tor' , disp)
      sr%n_tor0 = sr%n_tor

      call read_array_not_distr(filehandle, 'i_tor', (/sr%n_tor,0/), disp, int1d=sr%i_tor, deb=.false.)

      if ( sr%file_version >= 3) then
         if (.not. CARIDDI_mode)   sr%ntri_c                 = read_intparam_parallel(filehandle, 'ntri_c',          disp)
         sr%n_pol_coils            = read_intparam_parallel(filehandle, 'n_pol_coils',     disp)
         sr%n_rmp_coils            = read_intparam_parallel(filehandle, 'n_rmp_coils',     disp)
         sr%n_voltage_coils        = read_intparam_parallel(filehandle, 'n_voltage_coils', disp)
         sr%n_diag_coils           = read_intparam_parallel(filehandle, 'n_diag_coils',    disp)

        if (sr%n_voltage_coils > 0 ) then
          write(*,*) 'ERROR: voltage_coils not yet implemented.'
          stop
        end if
      
        if ( sr%ncoil /= sr%n_pol_coils + sr%n_rmp_coils + sr%n_voltage_coils + sr%n_diag_coils) then
          write(*,*) 'ERROR: STARWALL response is inconsistent: ncoil does not match sum.'
          stop
        end if
        if (CARIDDI_mode) then
          sr%ind_start_coils = read_intparam_parallel(filehandle, 'ind_start_coils',     disp)
        else
           sr%ind_start_coils = 1
        end if
        sr%ind_start_pol_coils     = read_intparam_parallel(filehandle, 'ind_start_pol_coils',     disp)
        sr%ind_start_rmp_coils     = read_intparam_parallel(filehandle, 'ind_start_rmp_coils',     disp)
        sr%ind_start_voltage_coils = read_intparam_parallel(filehandle, 'ind_start_voltage_coils', disp)
        sr%ind_start_diag_coils    = read_intparam_parallel(filehandle, 'ind_start_diag_coils',    disp)

        if ( sr%ncoil > 0 .and. .not. CARIDDI_mode) then
          call read_array_not_distr(filehandle, 'jtri_c',        (/sr%ncoil,0/),  disp,  int1d=sr%jtri_c, deb=.false.)
          call read_array_not_distr(filehandle, 'x_coil',        (/sr%ntri_c,3/), disp,  float2d=sr%x_coil, deb=.false.)
          call read_array_not_distr(filehandle, 'y_coil',        (/sr%ntri_c,3/), disp,  float2d=sr%y_coil, deb=.false.)
          call read_array_not_distr(filehandle, 'z_coil',        (/sr%ntri_c,3/), disp,  float2d=sr%z_coil, deb=.false.)
          call read_array_not_distr(filehandle, 'phi_coil',      (/sr%ntri_c,3/), disp,  float2d=sr%phi_coil, deb=.false.)
          call read_array_not_distr(filehandle, 'eta_thin_coil', (/sr%ntri_c,0/), disp,  float1d=sr%eta_thin_coil, deb=.false.)
          call read_array_not_distr(filehandle, 'coil_resist',   (/sr%ncoil,0/),  disp,  float1d=sr%coil_resist, deb=.false.)
          if (sr%file_version > 3)  then
            call read_array_not_distr(filehandle, 'coil_name',     (/sr%ncoil,0/),  disp,  char1d=sr%coil_name, deb=.false.)
           endif
        end if
        
        if ( (sr%ncoil .gt. 0) .and. sr%file_version .le. 3 .or. CARIDDI_mode) then
         write(*,*) "Coil names not yet supported by starwall (ver <=3). Generic names used."
         if (allocated(sr%coil_name)) deallocate(sr%coil_name)
         allocate(sr%coil_name(sr%ncoil))

         do i=1, sr%n_diag_coils
           write(tmp_int2str,'(I3.3)'), i
           sr%coil_name(sr%ind_start_diag_coils + i - 1) = trim("Dia_"//trim(tmp_int2str))
         end do

         do i=1, sr%n_pol_coils
           write(tmp_int2str,'(I3.3)'), i
           sr%coil_name(sr%ind_start_pol_coils + i - 1) = trim("PF_"//trim(tmp_int2str))
         end do
         
         do i=1, sr%n_rmp_coils
           write(tmp_int2str,'(I3.3)'), i
           sr%coil_name(sr%ind_start_rmp_coils + i - 1) = trim("RMP_"//trim(tmp_int2str))
         end do
       endif
     end if

      ! --- eta_thin_w is only part of the STARWALL response file since file_version 2
      if ( sr%file_version >= 2 ) then
        allocate(tmp(1))
        call read_array_not_distr(filehandle, 'eta_thin_w', (/1,0/), disp, float1d=tmp, deb=.false.)
        sr%eta_thin_w = tmp(1)
        deallocate(tmp)
      else
        sr%eta_thin_w = 0.
      end if

      if (CARIDDI_mode .and. halo_freeb) then
        call read_array_not_distr(filehandle, 'yy_ii', (/sr%n_nw,0 /), disp, float1d=sr%dii_yy , deb=.false.  )
      else
        call read_array_not_distr(filehandle, 'yy', (/sr%n_w,0/), disp, float1d=sr%d_yy, deb=.false.)
      endif
    end if

    if(my_id==0) call memory_prediction(sr%n_w, sr%nd_bez, ntasks)

    call MPI_BCAST(sr%n_w,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
    call MPI_BCAST(sr%n_nw,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
    call MPI_BCAST(sr%nd_bez, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
    call MPI_BCAST(sr%n_pts_real_w, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)  
    call MPI_BCAST(sr%n_pts_real_J, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
    call MPI_BCAST(sr%n_electrode_w, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)  
    call MPI_BCAST(sr%n_electrode_J, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
    call MPI_bcast(sr%ind_start_coils,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, err)
    call MPI_bcast(sr%ncoil,                   1, MPI_INTEGER,  0, MPI_COMM_WORLD, err)

    ! Last parameter defines how matrix should be distributed: .true. - rowwise; .false. - columnwise
    if (CARIDDI_mode .and. halo_freeb) then
      call read_array_par(filehandle, 'yie', (/sr%n_nw, sr%nd_bez/), disp, my_id, sr%a_yie, .false. ) 
      call read_array_par(filehandle, 'd_yieqe', (/sr%n_nw, sr%n_electrode_J-1/), disp, my_id, sr%d_yieqe, .true. ) 
      call read_array_par(filehandle, 'e_yiwe', (/sr%n_nw, sr%n_electrode_w-1/), disp, my_id, sr%e_yiwe, .true. ) 
      call read_array_par(filehandle, 'g_yiwe', (/sr%n_nw, sr%n_electrode_w-1/), disp, my_id, sr%g_yiwe, .true. ) 
    else
      call read_array_par    (filehandle, 'ye',       (/sr%n_w,sr%nd_bez/),    disp, my_id,  sr%a_ye,     .false.) 
    endif

    if (CARIDDI_mode .and. halo_freeb) then
      call read_array_par  (filehandle, 'eyi',    (/sr%nd_bez,sr%n_nw/),    disp, my_id,  sr%a_eyi,     .true. )
      call read_array_par  (filehandle, 'eye',    (/sr%nd_bez,sr%n_electrode_w-1/),    disp, my_id,  sr%a_eye,  .true. )
    else
      call read_array_par    (filehandle, 'ey',       (/sr%nd_bez,sr%n_w/),    disp, my_id,  sr%a_ey,     .true. )
    endif

    call read_array_par    (filehandle, 'ee',       (/sr%nd_bez,sr%nd_bez/), disp, my_id,  sr%a_ee,     .true. )
    
    ! =========================================  Aeq_el =============================================================
    if (CARIDDI_mode .and. halo_freeb) then
      call read_array_par(filehandle, 'eq_el',  (/sr%nd_bez,sr%n_electrode_J-1/), disp, my_id,  sr%aeq_el,     .true. )
    endif
    ! ============================== S_ww_inv, S_ww_ii_inv ===========================================================
    if (.not. vacuum_min) then
      if (.not. CARIDDI_mode) call read_array_par    (filehandle, 's_ww',     (/sr%n_w,sr%n_w/),       disp, my_id,  sr%s_ww,     .false.)
      if (CARIDDI_mode .and. halo_freeb) then
        call read_array_par(filehandle, 's_ww_ii_inv', (/sr%n_nw,sr%n_w/),       disp, my_id,  sr%s_ww_ii_inv, .true., deb=.true.)
      else
        call read_array_par    (filehandle, 's_ww_inv', (/sr%n_w,sr%n_w/),       disp, my_id,  sr%s_ww_inv, .true.)
      endif
    else
      if (.not. CARIDDI_mode) call read_array_par    (filehandle, 's_ww',     (/sr%n_w,sr%n_w/),       disp, my_id,  sr%s_ww,     .false. , &
          loc_start=sr%ind_start_coils, loc_len=sr%ncoil, loc_row=.true. )
      if (CARIDDI_mode .and. halo_freeb) then
        call read_array_par    (filehandle, 's_ww_ii_inv', (/sr%n_nw,sr%n_w/), disp, my_id,  sr%s_ww_ii_inv, .true., &
          loc_start=sr%ind_start_coils, loc_len=sr%ncoil, loc_row=.false. )
      else
        call read_array_par    (filehandle, 's_ww_inv', (/sr%n_w,sr%n_w/),       disp, my_id,  sr%s_ww_inv, .true., &
          loc_start=sr%ind_start_coils, loc_len=sr%ncoil, loc_row=.false. )
      endif
      sr%ind_start_coils = 1
    end if
    if (halo_freeb) then
      if (CARIDDI_mode)  call read_array_par(filehandle, 'FwtFw_inv', (/sr%n_electrode_w-1, sr%n_electrode_w-1/), disp, my_id, sr%fwtfw_inv, .true.)
      ! M_eye * (Fw^T Fw)^-1
      if (CARIDDI_mode)  call read_array_par(filehandle, 'W_h', (/sr%nd_bez, sr%n_electrode_w-1/) , disp, my_id, sr%w_h,  .true.  )
      if (my_id == 0 .and. CARIDDI_mode) then ! use read_array_not_distr ?
        call read_array_not_distr(filehandle, 'ielsphi_w', (/sr%n_pts_real_w,3/), disp, float2d=sr%ielsphi_real_w, deb=.false.)
        call read_array_not_distr(filehandle, 'ielsphi_J', (/sr%n_pts_real_J,3/), disp, float2d=sr%ielsphi_real_J, deb=.false.)
        ! write(*,*) 'Reading Jmat_real_w'
        ! call read_array_not_distr(filehandle, 'Jmat_real_w', (/sr%n_electrode_w,sr%n_pts_real_w/), disp, float2d=sr%Jmat_real_w, deb=.true.)
        ! write(*,*) 'Reading Jmat_real_J'
        ! call read_array_not_distr(filehandle, 'Jmat_real_J', (/sr%n_electrode_J,sr%n_pts_real_J/), disp, float2d=sr%Jmat_real_J, deb=.true.) 
        call read_array_not_distr(filehandle, 'Jmat_w', (/sr%n_electrode_w,sr%nd_bez /), disp, float2d=sr%Jmat_w, deb=.false.)
        call read_array_not_distr(filehandle, 'Jmat_J', (/sr%n_electrode_J,sr%nd_bez /), disp, float2d=sr%Jmat_J, deb=.false.)
      endif   

      ! READ JMAT AS DISTRIBUTED 

      ! if (CARIDDI_mode) then
      !   call read_array_par(filehandle, 'Jmat_w', (/sr%n_electrode_w,sr%nd_bez/),       disp, my_id,  sr%Jmat_w, .false., deb=.true.)
      !   call read_array_par(filehandle, 'Jmat_J', (/sr%n_electrode_J,sr%nd_bez/),       disp, my_id,  sr%Jmat_J, .false., deb=.true.)
      ! end if   
    endif

    if(my_id == 0 .and. .not. CARIDDI_mode) then
     call read_array_not_distr(filehandle, 'xyzpot_w', (/sr%npot_w,3/), disp, float2d=sr%xyzpot_w, deb=.false.)
     call read_array_not_distr(filehandle, 'jpot_w',   (/sr%ntri_w,3/), disp, int2d=sr%jpot_w, deb=.false.)
     if ( sr%file_version >= 5 ) then
       call read_array_not_distr(filehandle, 'phi0_w',   (/sr%ntri_w,3/), disp, float2d=sr%phi0_w, deb=.false.)
     endif
    end if

    call MPI_FILE_CLOSE(filehandle, err)

    if ( vacuum_debug .and. (my_id==0) ) write(*,*) 'Finished reading vacuum response.'

    ! --- Import normalization
    sr%a_ee%loc_mat(:,:) = sr%a_ee%loc_mat(:,:) * 2.d0*PI
    if (CARIDDI_mode .and. halo_freeb) then
      ! Change of sign to make it consistent with eddy coupling (check documentation!) 
      sr%a_yie%loc_mat(:,:)= -sr%a_yie%loc_mat(:,:)* 2.d0*PI
      sr%d_yieqe%loc_mat(:,:) = sr%d_yieqe%loc_mat(:,:) * 2.d0*PI
      sr%e_yiwe%loc_mat(:,:) = sr%e_yiwe%loc_mat(:,:) * 2.d0*PI
      sr%g_yiwe%loc_mat(:,:) = sr%g_yiwe%loc_mat(:,:) * 2.d0*PI
    else
      sr%a_ye%loc_mat(:,:) = sr%a_ye%loc_mat(:,:) * 2.d0*PI
    endif

    if ( vacuum_debug .and. (my_id==0) ) write(*,*) 'Applied import normalization.'

    ! --- Compute ideal-wall and no-wall response matrices.
    if (.not. resistive_wall) then
      ! IF HALO+EDDY COUPLING
      if (halo_freeb .and. CARIDDI_mode) then
        call matrix_multiplication(my_id, sr%a_eyi, mat2=sr%a_yie, res_mat=sr%a_id)
        sr%a_id%loc_mat(:,:) = sr%a_ee%loc_mat(:,:) - sr%a_id%loc_mat(:,:)
      else
      ! IF ONLY EDDY COUPLING
        call matrix_multiplication(my_id,sr%a_ey,mat2=sr%a_ye, res_mat=sr%a_id )
        sr%a_id%loc_mat(:,:) = sr%a_ee%loc_mat(:,:) - sr%a_id%loc_mat(:,:)
      end if
    endif
    ! --- Transform STARWALL harmonics to account for periodicity
    if ( my_id==0 ) then
      j = 0
      do i = 1, sr%n_tor
        i_starw = sr%i_tor(i)
        n       = i_starw / 2
        is_sin  = i_starw - 2 * n
        i_starw = 2 * n/n_period + is_sin
        if ( (mod(n, n_period) /= 0) .or. (i_starw < 1) .or. (i_starw > n_tor) ) then
          write(*,*) 'WARNING: STARWALL harmonic has no JOREK equivalent!'
          write(*,*) 'i_starw    =', sr%i_tor(i)
          write(*,*) 'n_period   =', n_period
          write(*,*) 'n_tor      =', n_tor
        else
          j = j + 1
          sr%i_tor(j) = i_starw
        end if
      end do
     
      sr%n_tor = j
      if ( vacuum_debug) write(*,*) 'End of routine read_starwall_response.'
    end if

  end subroutine read_starwall_response
  
  
  
  !> Calculate c . M . D . alpha or c . D . M . alpha as extended matrix-vector product
  !! c is a constant
  !! alpha is a vector
  !! M is a response matrix (distributed and/or compressed)
  !! D is a diagonal matrix (represented as vector)
  !! if MD==.true., calculate c . M . D . alpha, otherwise c . D . M. alpha
  !! WARNING: res is assumed to be initialized
  subroutine extended_matrix_vector(my_id, c, M, D, alpha, MD, res, red)

    use mpi_mod
    implicit none
  
    ! --- Input parameters
    integer,              intent(in)    :: my_id
    real*8,               intent(in)    :: c
    type(t_distrib_mat),  intent(in)    :: M
    real*8, allocatable,  intent(in)    :: D(:)
    real*8,               intent(in)    :: alpha(:)
    logical,              intent(in)    :: MD
    real*8, allocatable,  intent(inout) :: res(:)
    logical,              intent(in)    :: red
  
    ! --- Local variables
    logical :: correct_dims
    integer :: dim(2)
    integer :: j, jloc, jglob
    integer :: ntasks
    integer :: ierr

    !> Find the number of ranks involved in the calculation
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
  
    !> WARNING: up to now there is no check regarding the matrix representation
    !!          but we would need to add that when we will implemet the compressed
    !!          matrices' treatment

    ! --- Check for correct dimensions
    if ( MD ) then
      correct_dims = ( M%dim(2)==size(D)       ) .and. &
                     ( M%dim(2)==size(alpha,1) ) .and. &
                     ( M%dim(1)==size(res)     )
    else
      correct_dims = ( M%dim(1)==size(D)       ) .and. &
                     ( M%dim(2)==size(alpha,1) ) .and. &
                     ( M%dim(1)==size(res)     )
    end if
    if ( .not. correct_dims ) then
      write(*,*) 'Error in extended_matrix_vector: wrong dimensions of input'
      write(*,*) 'M:      ', M%dim
      write(*,*) 'D:      ', size(D)
      write(*,*) 'alpha:  ', size(alpha)
      write(*,*) 'result: ', size(res)
      stop
    end if
  
    !> WARNING: Given that we are not checking what king of distribution is
    !!          for the matrix M, we here need to differentiate the tratment for
    !!          rowwise and columnwise JOREK's distributions
    !> NOTE: Recalling the indexing convention for JOREK's kind of distribution,
    !!       the relation between the local and global extrema of indices (for
    !!       rows in rowise or column in columnwise) is as follows:
    !!
    !!       -----------------------------------------------------------
    !!       | LOCAL |               GLOBAL                            |
    !!       -----------------------------------------------------------
    !!       |   1   | matrix%ind_start = my_id*matrix%step + 1        |
    !!       | last  | matrix%ind_end   = my_id*matrix%step + loc_size |
    !!       -----------------------------------------------------------
    !!
    !!       with loc_size = matrix%step or matrix%step + dim - ntasks * matrix%step
    !!       in the case of the last task.
    !!       Therefore, we could probably generalize with the following
    !!       definition of local indices from the global one:
    !!
    !!       jglob = my_id*matrix%step + jloc
    if (M%row_wise) then
      !> Here the matrix is distributed row-wisely, therefore the rows'
      !! index runs from matrix%ind_start to matrix%ind_end, and each MPI 
      !! rank and the result array has its part of each matrix column. The
      !! array can be populated in parallel if we apply the index conversion 
      !! reported above for the rows' indices
      do jglob = M%ind_start, M%ind_end
        jloc = jglob - my_id*M%step
        if ( MD ) then
          !> The general j-term of the result vector, may be written as
          !!
          !! result(j) = Sum_{k=1}^{M%dim(2)} [ M(j,k) * D(k,k) * alpha(k) ],
          !! with j = 1,...,M%dim(1)
          res(jglob) = res(jglob) + c * sum( M%loc_mat(jloc,:) * D(:) * &
                                             alpha(:) )
        else
          !> The general j-term of the result vector, may be written as
          !!
          !! result(j) = D(j,j) * { Sum_{k=1}^{M%dim(2)} [ M(j,k) * alpha(k) ] },
          !! with j = 1,...,M%dim(1)
          res(jglob) = res(jglob) + c * D(jglob) * sum( M%loc_mat(jloc,:) * &
                                                        alpha(:) )
        end if
      end do
    else ! M%row_wise
      !> Here the matrix is distributed column-wisely, therefore the rows'
      !! index runs from 1 to M%dim(1) on each MPI rank, but each rank has only
      !! loc_size columns, ranging from M%ind_start to M%ind_end. Therefore, we
      !! need a final MPI_SUM reduction of the array result, in order to make 
      !! it include the contributions from all the ranks
      do j = 1, M%dim(1)
        if ( MD ) then
          !> The general j-term of the result vector, may be written as
          !!
          !! result(j) = Sum_{k=1}^{M%dim(2)} [ M(j,k) * D(k,k) * alpha(k) ],
          !! with j = 1,...,M%dim(1)
          res(j) = res(j) + c * sum( M%loc_mat(j,:) * D(M%ind_start:M%ind_end) * &
                                     alpha(M%ind_start:M%ind_end) )
        else
          !> The general j-term of the result vector, may be written as
          !!
          !! result(j) = D(j,j) * { Sum_{k=1}^{M%dim(2)} [ M(j,k) * alpha(k) ]
          !},
          !! with j = 1,...,M%dim(1)
          res(j) = res(j) + c * D(j) * sum( M%loc_mat(j,:) * &
                                            alpha(M%ind_start:M%ind_end) )
        end if
      end do
    end if ! M%row_wise

    !> Only if red, we need to call MPI_ALLREDUCE here, in order to add the
    !! contributions from all the ranks
    if (red) call MPI_ALLReduce(MPI_IN_PLACE, res, size(res), MPI_DOUBLE_PRECISION, &
                                MPI_SUM, MPI_COMM_WORLD, ierr)


  end subroutine extended_matrix_vector 

  
  
  !> Routine for multiplying two (distributed) matrices.
  subroutine matrix_multiplication(my_id, mat1, mat2, mat2_not_distr, res_mat, res_mat_not_distr)
    
    use mpi_mod
  
    implicit none
  
    ! --- Routine parameters
    integer,                        intent(in)      :: my_id
    type(t_distrib_mat),            intent(in)      :: mat1
    type(t_distrib_mat), optional,  intent(in)      :: mat2
    real*8, allocatable, optional,  intent(inout)   :: mat2_not_distr(:,:)
    type(t_distrib_mat), optional,  intent(inout)   :: res_mat
    real*8, allocatable, optional,  intent(inout)   :: res_mat_not_distr(:,:)
  
    ! --- Local variables
    integer             :: ntasks, ierr, i, k, j, z, length, glob_index_i, glob_index_j 
    real*8, allocatable :: tmp(:,:)
    real*8              :: sum_element
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
  
    ! --- Multiplication mat1(distributed)*mat2(distributed) =res_mat(distributed)
    if ( present(res_mat) .and. present(mat2) ) then
  
      call  alloc_distr(my_id, res_mat, (/mat1%dim(1), mat2%dim(2)/), .true.)       
      res_mat%loc_mat=0.0
    
      if ( res_mat%row_wise .and. mat1%row_wise .and. (.not. mat2%row_wise) ) then
        
       if(mat1%dim(2) /= mat2%dim(1) ) then
          write(*,*) "ERROR: Dimensions of  multiplied matrices (1) #col =", mat1%dim(2), &
                     " #row =", mat2%dim(1) , "STOP (1)"
          stop
       endif   
 
        ! loop over all tasks
        do i = 1, ntasks
    
          ! Broadcasting size of the distributed column
          length = size(mat2%loc_mat,2)      
          call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr )
    
          ! temporal matrix for receiving part of the distributed matrix and
          ! multiply it
          if ( allocated(tmp) ) deallocate(tmp)
          allocate( tmp(size(mat2%loc_mat,1),length) )
    
          ! part of second matrix broadcasted to all tasks
          if(my_id == i-1)  tmp    = mat2%loc_mat  
          if(my_id == i-1 ) length = size(tmp)
          call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr)
          call MPI_BCAST( tmp, length, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, ierr)
       
          ! Matrix -matrix multiplications
          glob_index_j=(i-1)*mat2%step
          do z = 1, size(mat1%loc_mat,1)
            do k = 1, size(tmp,2)
    
              sum_element=0.0 
!$omp parallel private (j) shared (sum_element,mat1,z,k,tmp) 
!$omp do reduction (+:sum_element)
              do j = 1, size(mat1%loc_mat,2)
                sum_element = sum_element + mat1%loc_mat(z,j) * tmp(j,k)         
              end do
!$omp end do 
!$omp end parallel
              
              res_mat%loc_mat(z,k+glob_index_j) = sum_element  
        
            end do
          end do
    
        end do     
      
      else
        write(*,*) 'ERROR in matrix_multiplication: Combination of rowwise/columnwise not implemented (1).'
        stop
      end if 
    
    ! --- Multiplication mat1(distributed)*mat2(distributed) =res_mat(not_distributed)
    else if ( present(mat2) .and. present(res_mat_not_distr) ) then
      
      if (allocated(res_mat_not_distr)) deallocate(res_mat_not_distr)
      allocate( res_mat_not_distr(mat1%dim(1), mat2%dim(2)) )  

      if ( mat1%row_wise .and. (.not. mat2%row_wise) ) then
       
        if(mat1%dim(2) /= mat2%dim(1) ) then
          write(*,*) "ERROR: Dimensions of  multiplied matrices (2) #col =", mat1%dim(2), &
                     " #row =", mat2%dim(1) , "STOP (2)"
          stop
        endif
    
        res_mat_not_distr=0.0   
    
        ! loop over all tasks
        do i = 1, ntasks
    
          ! Broadcasting size of the distributed column
          length = size(mat2%loc_mat,2)
        
          call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr )
    
          ! temporal matrix for receiving part of the distributed matrix and
          ! multiply it
          if ( allocated(tmp) ) deallocate(tmp)
          allocate( tmp(size(mat2%loc_mat,1),length) )
    
          ! part of second matrix broadcasted to all tasks
          if(my_id == i-1)  tmp    = mat2%loc_mat
          if(my_id == i-1 ) length = size(tmp)
          call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr)
          call MPI_BCAST(tmp, length, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD,ierr)
    
          ! Matrix -matrix multiplications
    
          glob_index_j = (i-1)*mat2%step
          glob_index_i = my_id*mat1%step ! R.S> Correcting from mat2%step to mat1%step
          do z = 1, size(mat1%loc_mat,1)
            do k = 1, size(tmp,2)
    
              sum_element=0.0

!$omp parallel private (j) shared (sum_element,mat1,z,k,tmp) 
!$omp do reduction (+:sum_element)
              do j = 1, size(mat1%loc_mat,2)
                sum_element = sum_element + mat1%loc_mat(z,j) * tmp(j,k)
              end do
!$omp end do 
!$omp end parallel

              res_mat_not_distr(z + glob_index_i, k + glob_index_j) = sum_element
    
            end do
          end do
    
        end do
!       else if ( mat1%row_wise .and. mat2%row_wise ) then
       
!         if(mat1%dim(2) /= mat2%dim(1) ) then
!           write(*,*) "ERROR: Dimensions of  multiplied matrices (2) #col =", mat1%dim(2), &
!                      " #row =", mat2%dim(1) , "STOP (2)"
!           stop
!         endif
    
!         res_mat_not_distr=0.0   
    
!         ! loop over all tasks
!         do i = 1, ntasks
    
!           ! Broadcasting size of the distributed row
!           length = size(mat2%loc_mat,1)
        
!           call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr )
    
!           ! temporal matrix for receiving part of the distributed matrix and
!           ! multiply it
!           if ( allocated(tmp) ) deallocate(tmp)
!           allocate( tmp(length,size(mat2%loc_mat,2)) )
    
!           ! part of second matrix broadcasted to all tasks
!           if(my_id == i-1)  tmp    = mat2%loc_mat
!           if(my_id == i-1 ) length = size(tmp)
!           call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr)
!           call MPI_BCAST(tmp, length, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD,ierr)
    
!           ! Matrix -matrix multiplications
!           glob_index_j = (i-1)*mat2%step
!           glob_index_i = my_id*mat2%step

!           do z = 1, size(mat1%loc_mat,1) ! Dist
!             do k = 1, size(tmp,2) ! Full
    
!               sum_element=0.0

! !$omp parallel private (j) shared (sum_element,mat1,z,k,tmp) 
! !$omp do reduction (+:sum_element)
!               do j = 1, size(mat2%loc_mat,1) ! Dist
!                 sum_element = sum_element + mat1%loc_mat(z,j+glob_index_i) * tmp(j,k)
!               end do
! !$omp end do 
! !$omp end parallel

!               res_mat_not_distr(z + glob_index_i, k ) = sum_element
    
!             end do
!           end do
    
!         end do
      else
        write(*,*) 'ERROR in matrix_multiplication: Combination of rowwise/columnwise not implemented (2).'
        stop
      end if
     
    ! --- Multiplication mat1(distributed)*mat2(not_distributed) =res_mat(not_distributed)
    else if(present(mat2_not_distr) .AND. present(res_mat_not_distr)) then

       if (allocated(res_mat_not_distr)) deallocate(res_mat_not_distr)
       allocate( res_mat_not_distr(mat1%dim(1), size(mat2_not_distr,2)) )

       if(mat1%dim(2) /= size(mat2_not_distr,1) ) then
          write(*,*) "ERROR: Dimensions of  multiplied matrices (3) #col =", mat1%dim(2), &
                     " #row =", size(mat2_not_distr,1) , "STOP (3)"
          stop
       endif
       
      !Check distribution. Result matrix not distributed, first matrix row wise, second matrix not distributed
      if (mat1%row_wise) then
    
        res_mat_not_distr=0.0
        
        ! Matrix -matrix multiplications     
        glob_index_i = my_id*mat1%step
    
        do z = 1, size(mat1%loc_mat,1)
          do k = 1, size(mat2_not_distr,2)
    
            sum_element=0.0

!$omp parallel private (j) shared (sum_element,mat1,z,k,tmp,mat2_not_distr) 
!$omp do reduction (+:sum_element) 
            do j = 1, size(mat1%loc_mat,2)
              sum_element = sum_element + mat1%loc_mat(z,j) *mat2_not_distr(j,k)
            end do
!$omp end do 
!$omp end parallel

            res_mat_not_distr(z + glob_index_i, k) = sum_element
    
          end do
        end do
      else 
        write(*,*) 'ERROR in matrix_multiplication: Combination of rowwise/columnwise not implemented (3).'
        stop
      end if
        
!       else if (.not. mat1%row_wise) then   ! R.S>, need to check for correctness
!         ! mat1 is column-wise distributed
!         res_mat_not_distr=0.0
        
!         ! Matrix -matrix multiplications     
!         glob_index_i = my_id*mat1%step
!         do z = 1, size(mat1%loc_mat,1)
!           do k = 1, size(mat2_not_distr,2)
    
!             sum_element=0.0

! !$omp parallel private (j) shared (sum_element,mat1,z,k,tmp,mat2_not_distr) 
! !$omp do reduction (+:sum_element) 
!             do j = 1, size(mat1%loc_mat,2)
!               sum_element = sum_element + mat1%loc_mat(z,j) *mat2_not_distr(j+glob_index_i,k)
!             end do
! !$omp end do 
! !$omp end parallel

!             res_mat_not_distr(z , k) = sum_element
    
!           end do
!         end do
!     endif
    
    else
      write(*,*) 'ERROR in matrix_multiplication: Unsupported set of matrices provided.'
      stop

    end if
  
  end subroutine matrix_multiplication
  
  
  
  !> Broadcast the STARWALL response matrices to the other MPI procs.
  subroutine broadcast_starwall_response(my_id, sr)
    
    use mpi_mod
    
    implicit none
    
    
    ! --- Routine parameters
    integer,                   intent(in)    :: my_id
    type(t_starwall_response), intent(inout) :: sr
    
    ! --- Local parameters
    integer :: ierr
    real*8    :: test_sum, loc_sum    

    if ( vacuum_debug .and. (my_id == 0) ) write(*,*) my_id, 'Entering broadcast_starwall_response.'
    
    ! --- Broadcast parameters.
    call MPI_bcast(sr%file_version,            1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_bnd,                   1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%nd_bez,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    ! if (CARIDDI_mode .and. halo_freeb) then (We broadcast anyway?)
    call MPI_bcast(sr%n_pts_real_w,            1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_pts_real_J,            1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_electrode_w,            1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_electrode_J,            1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    !endif
    call MPI_bcast(sr%npot_w,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_w,                     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_nw,                    1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ntri_w,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_tor,                   1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_tor0,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)

    call MPI_bcast(sr%nv,                      1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_points,                1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%iwall,                   1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%nwu,                     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%nwv,                     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%mn_w,                    1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%max_mn_w,                1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)

    call MPI_bcast(sr%ntri_c,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_pol_coils,             1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_rmp_coils,             1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_voltage_coils,         1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_diag_coils,            1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ind_start_pol_coils,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ind_start_rmp_coils,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ind_start_voltage_coils, 1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ind_start_diag_coils,    1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%eta_thin_w,              1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    n_dof_starwall = sr%nd_bez
    n_wall_curr    = sr%n_w
    n_wall_curr_i  = sr%n_nw
    n_wall_curr_e  = sr%n_electrode_w-1
    n_shell_curr_e = sr%n_electrode_J-1

    ! --- Allocate matrices.
    if ( my_id /= 0 ) then
      if (allocated(sr%i_tor)   ) deallocate(sr%i_tor);    allocate(sr%i_tor(sr%n_tor))
      if (CARIDDI_mode .and. halo_freeb) then
        if (ALLOCATED(sr%dii_yy)) deallocate(sr%dii_yy);   allocate(sr%dii_yy(sr%n_nw))
      else
        if (allocated(sr%d_yy)    ) deallocate(sr%d_yy);     allocate(sr%d_yy(sr%n_w))
      endif
      if (CARIDDI_mode .and. halo_freeb) then
        if (allocated(sr%ielsphi_real_w)) deallocate(sr%ielsphi_real_w); allocate(sr%ielsphi_real_w(sr%n_pts_real_w,3))
        if (allocated(sr%ielsphi_real_J)) deallocate(sr%ielsphi_real_J); allocate(sr%ielsphi_real_J(sr%n_pts_real_J,3))
        if (allocated(sr%Jmat_real_w)) deallocate(sr%Jmat_real_w); allocate(sr%Jmat_real_w(sr%n_electrode_w,sr%n_pts_real_w))
        if (allocated(sr%Jmat_real_J)) deallocate(sr%Jmat_real_J); allocate(sr%Jmat_real_J(sr%n_electrode_J,sr%n_pts_real_J))
        if (allocated(sr%Jmat_w)) deallocate(sr%Jmat_w); allocate(sr%Jmat_w(sr%n_electrode_w,sr%nd_bez))
        if (allocated(sr%Jmat_J)) deallocate(sr%Jmat_J); allocate(sr%Jmat_J(sr%n_electrode_J,sr%nd_bez))
      endif
      if (.not. CARIDDI_mode) then
        if (allocated(sr%xyzpot_w)) deallocate(sr%xyzpot_w); allocate(sr%xyzpot_w(sr%npot_w,3))
        if (allocated(sr%jpot_w)  ) deallocate(sr%jpot_w);   allocate(sr%jpot_w(sr%ntri_w,3))
      end if
      if ( sr%file_version>=5 .and. .not. CARIDDI_mode) then 
        if (allocated(sr%m_w)        ) deallocate(sr%m_w);          allocate(        sr%m_w(sr%max_mn_w))
        if (allocated(sr%n_w_fourier)) deallocate(sr%n_w_fourier);  allocate(sr%n_w_fourier(sr%max_mn_w))

        if (allocated(sr%rc_w)       ) deallocate(sr%rc_w);         allocate(sr%rc_w(sr%max_mn_w))
        if (allocated(sr%rs_w)       ) deallocate(sr%rs_w);         allocate(sr%rs_w(sr%max_mn_w))
        if (allocated(sr%zc_w)       ) deallocate(sr%zc_w);         allocate(sr%zc_w(sr%max_mn_w))
        if (allocated(sr%zs_w)       ) deallocate(sr%zs_w);         allocate(sr%zs_w(sr%max_mn_w))
        if (allocated(sr%phi0_w)     ) deallocate(sr%phi0_w);       allocate(sr%phi0_w(sr%ntri_w,3))
      endif

      if ( sr%ncoil > 0) then
        if ( .not. CARIDDI_mode ) then
          if (allocated(sr%jtri_c)       ) deallocate(sr%jtri_c);        allocate(sr%jtri_c(sr%ncoil))
         if (allocated(sr%x_coil)       ) deallocate(sr%x_coil);        allocate(sr%x_coil(sr%ntri_c,3))
          if (allocated(sr%y_coil)       ) deallocate(sr%y_coil);        allocate(sr%y_coil(sr%ntri_c,3))
          if (allocated(sr%z_coil)       ) deallocate(sr%z_coil);        allocate(sr%z_coil(sr%ntri_c,3))
          if (allocated(sr%phi_coil)     ) deallocate(sr%phi_coil);      allocate(sr%phi_coil(sr%ntri_c,3))
          if (allocated(sr%eta_thin_coil)) deallocate(sr%eta_thin_coil); allocate(sr%eta_thin_coil(sr%ntri_c))
          if (allocated(sr%coil_resist)  ) deallocate(sr%coil_resist);   allocate(sr%coil_resist(sr%ncoil))
        end if
        if (allocated(sr%coil_name)  )   deallocate(sr%coil_name);     allocate(sr%coil_name(sr%ncoil))
      end if
    end if
    
    ! --- Broadcast matrices.
    call MPI_bcast(sr%i_tor,    sr%n_tor,                MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    if (CARIDDI_mode .and. halo_freeb) then
      call MPI_bcast(sr%dii_yy, sr%n_nw,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%ielsphi_real_w, sr%n_pts_real_w*3,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%ielsphi_real_J, sr%n_pts_real_J*3,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      !call MPI_bcast(sr%Jmat_real_w, sr%n_electrode_w*sr%n_pts_real_w,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      !call MPI_bcast(sr%Jmat_real_J, sr%n_electrode_J*sr%n_pts_real_J,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%Jmat_w, sr%n_electrode_w*sr%nd_bez,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%Jmat_J, sr%n_electrode_J*sr%nd_bez,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    else
      call MPI_bcast(sr%d_yy,     sr%n_w,                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    endif

    if (.not. CARIDDI_mode) then
      call MPI_bcast(sr%xyzpot_w, sr%npot_w*3,             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%jpot_w,   sr%ntri_w*3,             MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    end if
    
    if ( sr%file_version>=5 .and. .not. CARIDDI_mode) then 
      call MPI_bcast(sr%m_w,            sr%max_mn_w,     MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%n_w_fourier,    sr%max_mn_w,     MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%rc_w,           sr%max_mn_w,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%rs_w,           sr%max_mn_w,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%zc_w,           sr%max_mn_w,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%zs_w,           sr%max_mn_w,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%phi0_w,         sr%ntri_w*3,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    endif

    if ( sr%ncoil > 0 ) then
      if (.not. CARIDDI_mode) then
        call MPI_bcast(sr%jtri_c,   sr%ncoil,              MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
        call MPI_bcast(sr%x_coil,   sr%ntri_c*3,           MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_bcast(sr%y_coil,   sr%ntri_c*3,           MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_bcast(sr%z_coil,   sr%ntri_c*3,           MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_bcast(sr%phi_coil,      sr%ntri_c*3,      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_bcast(sr%eta_thin_coil, sr%ntri_c,        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_bcast(sr%coil_resist,   sr%ncoil,         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      end if
      call MPI_bcast(sr%coil_name,sr%ncoil*COIL_NAME_LEN,MPI_CHARACTER       , 0, MPI_COMM_WORLD, ierr)    
    end if
    
    if ( vacuum_debug ) then
      if (.not. (CARIDDI_mode .and. halo_freeb)) then
        loc_sum = sum(sr%a_ye%loc_mat) + sum(sr%a_ey%loc_mat) + sum(sr%a_ee%loc_mat) + &
                  sum(sr%s_ww_inv%loc_mat)   

        if (my_id==0) then

          loc_sum =   loc_sum + sum(sr%i_tor) + sum(sr%d_yy)     &
                    + sr%n_bnd + sr%nd_bez + sr%ncoil + sr%npot_w &
                    + sr%n_w + sr%ntri_w + sr%n_tor + sr%eta_thin_w                &
                    + sr%file_version + sr%ntri_c + sr%n_pol_coils                 &
                    + sr%n_rmp_coils + sr%n_voltage_coils + sr%n_diag_coils        &
                    + sr%ind_start_pol_coils + sr%ind_start_rmp_coils              &
                    + sr%ind_start_voltage_coils + sr%ind_start_diag_coils
        end if
    
        test_sum=0.0
        call MPI_Reduce(loc_sum, test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  
        if(my_id == 0)  write(*,'("Checksum",I4,ES24.16)') my_id, test_sum
      endif
      if(my_id == 0)  write(*,*) my_id, 'Exiting broadcast_starwall_response.'
   
    end if  
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  end subroutine broadcast_starwall_response
  
  

  !> Write out information about the STARWALL response matrices.
  subroutine log_starwall_response(my_id, sr)

    use mod_parameters, only: n_period
    use phys_module, only: resistive_wall
    use mpi_mod

    implicit none

    integer,                   intent(in) :: my_id
    type(t_starwall_response), intent(in) :: sr

    real*8  :: test_sum
    integer :: ierr

    32 format(3x,77('-'))
    33 format(3x,a,i8)
    34 format(3x,'sum(',a,')=',es24.16)
    35 format(3x,'sum(',a,')=',i24)
    36 format(3x,'sum(',a,')= ---not allocated---')
    37 format(3x,a,es25.15)
  
    if (my_id == 0) then
       write(*,*)
       write(*,32)
       write(*,33) 'STARWALL RESPONSE INFORMATION:'
       write(*,32)
       write(*,33) 'file_version            =', sr%file_version
       write(*,33) 'n_bnd                   =', sr%n_bnd
       write(*,33) 'nd_bez                  =', sr%nd_bez
       write(*,33) 'nv                      =', sr%nv    
       write(*,33) 'n_points                =', sr%n_points
       write(*,33) 'ncoil                   =', sr%ncoil
       write(*,33) 'npot_w                  =', sr%npot_w
       write(*,33) 'n_w                     =', sr%n_w
       if (CARIDDI_mode)        write(*,33) 'n_nw                     =', sr%n_nw

       write(*,33) 'ntri_w                  =', sr%ntri_w
       write(*,33) 'iwall                   =', sr%iwall  
       write(*,33) 'nwu                     =', sr%nwu    
       write(*,33) 'nwv                     =', sr%nwv    
       write(*,33) 'mn_w                    =', sr%mn_w   
       write(*,33) 'max_mn_w                =', sr%max_mn_w   
       write(*,33) 'n_tor                   =', sr%n_tor
       write(*,33) 'n_tor0                  =', sr%n_tor0
       write(*,33) 'n_pol_coils             =', sr%n_pol_coils
       write(*,33) 'n_rmp_coils             =', sr%n_rmp_coils
       write(*,33) 'n_voltage_coils         =', sr%n_voltage_coils
       write(*,33) 'n_diag_coils            =', sr%n_diag_coils
       write(*,33) 'ind_start_pol_coils     =', sr%ind_start_pol_coils
       write(*,33) 'ind_start_rmp_coils     =', sr%ind_start_rmp_coils
       write(*,33) 'ind_start_voltage_coils =', sr%ind_start_voltage_coils
       write(*,33) 'ind_start_diag_coils    =', sr%ind_start_diag_coils
      if ( sr%file_version >= 2) write(*,37) 'eta_thin_w        =', sr%eta_thin_w
      if (allocated(sr%i_tor)) write(*,33) 'i_tor               ='//trim(modes_to_str(sr%i_tor,sr%n_tor,n_period))
      if (vacuum_min) write(*,*) 'WARNING: vacuum_min = .true. this means you can only diagnose coil currents during the run'
    end if 


    if ( vacuum_debug ) then
      if (.not. (CARIDDI_mode .and. halo_freeb)) then
        if (my_id == 0) then
          write(*,32)
          if (allocated(sr%i_tor        )) then; write(*,35) 'i_tor        ',sum(sr%i_tor         ); else; write(*,36) 'i_tor        '; end if
          if (allocated(sr%d_yy         )) then; write(*,34) 'd_yy         ',sum(sr%d_yy          ); else; write(*,36) 'd_yy         '; end if
        end if 

        test_sum=0.0
        if (allocated(sr%a_ye%loc_mat)) then
          call MPI_Reduce(sum(sr%a_ye%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if (my_id == 0) write(*,34) 'a_ye         ', test_sum 
        else
          if (my_id == 0)  write(*,36) 'a_ye         '
        end if            
  
        test_sum=0.0
        if (allocated(sr%a_ey%loc_mat)) then 
          call MPI_Reduce(sum(sr%a_ey%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if (my_id == 0) write(*,34) 'a_ey         ', test_sum 
        else
          if (my_id == 0)  write(*,36) 'a_ey         ' 
        end if

        test_sum=0.0
        if (allocated(sr%a_ee%loc_mat)) then
          call MPI_Reduce(sum(sr%a_ee%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if (my_id == 0) write(*,34) 'a_ee         ', test_sum
        else
          if (my_id == 0) write(*,36) 'a_ee         '
        end if

        if (.not. resistive_wall) then
          test_sum=0.0
          if (allocated(sr%a_id%loc_mat)) then
            call MPI_Reduce(sum(sr%a_id%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if (my_id == 0) write(*,34) 'a_id         ', test_sum
          else
            if (my_id == 0) write(*,36) 'a_id         '
          end if
        end if

        test_sum=0.0
        if (allocated(sr%s_ww%loc_mat)) then
          call MPI_Reduce(sum(sr%s_ww%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if (my_id == 0) write(*,34) 's_ww         ', test_sum
        else
          if (my_id == 0) write(*,36) 's_nw         '
        end if

        test_sum=0.0
        if (allocated(sr%s_ww_inv%loc_mat)) then
          call MPI_Reduce(sum(sr%s_ww_inv%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if (my_id == 0) write(*,34) 's_ww_inv     ', test_sum
        else
          if (my_id == 0) write(*,36) 's_nw_inv         '
        end if
      endif
      if (my_id == 0 .and. .not. CARIDDI_mode) then
        if (allocated(sr%xyzpot_w     )) then; write(*,34) 'xyzpot_w     ',sum(sr%xyzpot_w      ); else; write(*,36) 'xyzpot_w     '; end if
        if (allocated(sr%jpot_w       )) then; write(*,35) 'jpot_w       ',sum(sr%jpot_w        ); else; write(*,36) 'jpot_w       '; end if
        if (allocated(sr%jtri_c       )) then; write(*,35) 'jtri_c       ',sum(sr%jtri_c        ); else; write(*,36) 'jtri_c       '; end if
        if (allocated(sr%x_coil       )) then; write(*,34) 'x_coil       ',sum(sr%x_coil        ); else; write(*,36) 'x_coil       '; end if
        if (allocated(sr%y_coil       )) then; write(*,34) 'y_coil       ',sum(sr%y_coil        ); else; write(*,36) 'y_coil       '; end if
        if (allocated(sr%z_coil       )) then; write(*,34) 'z_coil       ',sum(sr%z_coil        ); else; write(*,36) 'z_coil       '; end if
        if (allocated(sr%phi_coil     )) then; write(*,34) 'phi_coil     ',sum(sr%phi_coil      ); else; write(*,36) 'phi_coil     '; end if
        if (allocated(sr%eta_thin_coil)) then; write(*,34) 'eta_thin_coil',sum(sr%eta_thin_coil ); else; write(*,36) 'eta_thin_coil'; end if
        if (allocated(sr%coil_resist  )) then; write(*,34) 'coil_resist  ',sum(sr%coil_resist   ); else; write(*,36) 'coil_resist  '; end if
      end if
    end if

    if (my_id==0) then
      write(*,32)
      write(*,*)
    end if
   
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
  end subroutine log_starwall_response
  
  
  
  !> Write out resistive-wall data as a VTK-file.
  !!
  !! Scalar quantities (at wall triangle nodes):
  !! * pot_w: Wall current potentials
  !! * dpot_w: Change of wall current potentials in previous time-step
  !!
  !! Vector quantities (at triangles):
  !! * jsurf_w: Surface currents on the wall
  subroutine write_wall_vtk(index, resistive_wall, my_id)
    
    use phys_module, only: nout
    use constants
    use mpi_mod    

    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: index !< Time step index
    logical, intent(in) :: resistive_wall
    integer, intent(in) :: my_id
    
    ! --- Local variables
    real*8              :: phi1, phi2, phi3, r1(3), r2(3), r3(3), r21(3), r32(3), r13(3), r21_cross_r32(3)
    real*8              :: nm(3), n13(3), n32(3), n21(3), j_lin(3), j13, j32, j21
    real*8              :: ephi12(3),ephi13(3),ephi23(3), jsides(3)
    real*8              :: pol13(3), pol32(3), pol21(3)
    real*8              :: mid12(3), mid23(3), mid13(3), angle12, angle13, angle23
    real*8              :: Iw_net_tor 
    integer             :: filehandle = 60, i, maxcurr_pos
    logical             :: Iphi_max, Ipol_max, jphi_lin, jpol_lin
    character(len=18)   :: filename
    real*8, allocatable :: tripot_w(:), dtripot_w(:)
    integer :: ierr

    if ( mod(index,nout) /= 0 ) return

    if(my_id == 0) then
      if ( mod(index,nout) /= 0 ) return
      
      ! --- Preset values
      Iphi_max = .false.
      jphi_lin = .false.
      Ipol_max = .false.
      jpol_lin = .false.    
      
      ! --- VTK file header
      write(filename,'(A,I5.5,A)') 'wallcurr.',index,'.vtk'
      open(filehandle, file=filename, status='replace', action='write')
      140 format(a)
      141 format(a,i8,a)
      142 format(3es16.8)
      143 format(a,2i8)
      144 format(4i8)
      write(filehandle,140) '# vtk DataFile Version 2.0'
      write(filehandle,140) 'testdata'
      write(filehandle,140) 'ASCII'
      write(filehandle,140) 'DATASET POLYDATA'
      
      ! --- Triangle node positions
      write(filehandle,141) 'POINTS', sr%npot_w, ' float'
      do i = 1, sr%npot_w
        write(filehandle,142) sr%xyzpot_w(i,:)
      end do
      
      ! --- Node indices corresponding to triangles
      write(filehandle,143) 'POLYGONS', sr%ntri_w, sr%ntri_w * 4
      do i = 1, sr%ntri_w
        write(filehandle,144) 3, sr%jpot_w(i,:) - 1
      end do
      
      ! --- Wall current potentials
      write(filehandle,141) 'POINT_DATA', sr%npot_w
      write(filehandle,140) 'SCALARS pot_w float'
      write(filehandle,140) 'LOOKUP_TABLE default'
    end if
    
    call reconstruct_triangle_potentials(tripot_w, wall_curr, my_id, Iw_net_tor)

    if(my_id == 0) then
      do i = 1, sr%npot_w
        write(filehandle,142) tripot_w(i)
      end do
      
      ! --- Change of wall current potentials in previous time-step
      write(filehandle,140) 'SCALARS dpot_w float'
      write(filehandle,140) 'LOOKUP_TABLE default'
    end if

    call reconstruct_triangle_potentials(dtripot_w, dwall_curr, my_id)

    if(my_id == 0) then
  
      do i = 1, sr%npot_w
        write(filehandle,142) dtripot_w(i)
      end do

      deallocate(dtripot_w)
      
      ! --- Cell data variables
      write(filehandle,141) 'CELL_DATA', sr%ntri_w
     
      ! --- Maximum toroidal current flwoing on a triangle (kA)
      if (Iphi_max) then
          
        write(filehandle,140) 'SCALARS Iphi_max(kA) float'
        write(filehandle,140) 'LOOKUP_TABLE default'
         
        do i = 1, sr%ntri_w
         
          ! --- Wall potentials at triangle nodes
          phi1   = tripot_w(sr%jpot_w(i,1))
          phi2   = tripot_w(sr%jpot_w(i,2))
          phi3   = tripot_w(sr%jpot_w(i,3))
          
          ! --- Positions of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          
          ! --- Middle points on triangle sides
          mid12(:) =   (r1(:) + r2(:)) * 0.5d0
          mid13(:) =   (r1(:) + r3(:)) * 0.5d0
          mid23(:) =   (r2(:) + r3(:)) * 0.5d0      
               
          ! --- Toroidal angles on middle points, JOREK coordinate system  y --> 1, x --> 3
          angle12 = atan2(-mid12(1),mid12(3))
          angle13 = atan2(-mid13(1),mid13(3))
          angle23 = atan2(-mid23(1),mid23(3))
          
          ! --- Toroidal basis vectors on middle points, in clock-wise direction looking from above the torus
          ephi12(:) = (/- cos(angle12), 0.d0, - sin(angle12) /)
          ephi13(:) = (/- cos(angle13), 0.d0, - sin(angle13) /)
          ephi23(:) = (/- cos(angle23), 0.d0, - sin(angle23) /)
         
          ! --- Quantities needed to calculate the current density vector 
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r13(:) = r3(:)-r1(:)
         
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          !--- current density vector in kA/m, Merkel 2015
          j_lin = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2))
          j_lin = j_lin / mu_zero * 1.d-3  
            
          !--- Vector normal to triangle surface
          nm(:)  = r21_cross_r32(:) / sqrt(sum(r21_cross_r32**2))
          
          !--- Normal vectors to triangle sides
          n13(:) = (/ r13(2)*nm(3) - r13(3)*nm(2), r13(3)*nm(1) - r13(1)*nm(3),          &
            r13(1)*nm(2) - r13(2)*nm(1) /)
          n21(:) = (/ r21(2)*nm(3) - r21(3)*nm(2), r21(3)*nm(1) - r21(1)*nm(3),          &
            r21(1)*nm(2) - r21(2)*nm(1) /)
          n32(:) = (/ r32(2)*nm(3) - r32(3)*nm(2), r32(3)*nm(1) - r32(1)*nm(3),          &
            r32(1)*nm(2) - r32(2)*nm(1) /)
          
          !--- Normalized vectors normal to triangle sides
          n13(:) = n13(:)/ sqrt(sum(n13**2))
          n21(:) = n21(:)/ sqrt(sum(n21**2))
          n32(:) = n32(:)/ sqrt(sum(n32**2))
          
          !--- Toroidal current flowing per side of the triangle 
          j13    = abs(dot_product(n13,ephi13)) * sqrt(sum(r13**2)) * dot_product(ephi13,j_lin)
          j21    = abs(dot_product(n21,ephi12)) * sqrt(sum(r21**2)) * dot_product(ephi12,j_lin)
          j32    = abs(dot_product(n32,ephi23)) * sqrt(sum(r32**2)) * dot_product(ephi23,j_lin)
          
          jsides(:) = (/ j13, j21, j32 /)
              
          !--- Find the triangle side with maximum absolute current
          maxcurr_pos   = maxloc(abs(jsides),1)
          
          !--- Maximum toroidal current flowing on a triangle side, in kA
          write(filehandle,142) jsides(maxcurr_pos) 
        
        end do
      end if
      
      ! --- Maximum poloidal current flwoing on a triangle (kA)
      if (Ipol_max) then
     
        write(filehandle,140) 'SCALARS Ipol_max(kA) float'
        write(filehandle,140) 'LOOKUP_TABLE default'
         
        do i = 1, sr%ntri_w
         
          ! --- Wall potentials at triangle nodes
          phi1   = tripot_w(sr%jpot_w(i,1))
          phi2   = tripot_w(sr%jpot_w(i,2))
          phi3   = tripot_w(sr%jpot_w(i,3))
          
          ! --- Positions of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          
          ! --- Middle points on triangle sides
          mid12(:) =   (r1(:) + r2(:)) * 0.5d0
          mid13(:) =   (r1(:) + r3(:)) * 0.5d0
          mid23(:) =   (r2(:) + r3(:)) * 0.5d0      
               
          ! --- Toroidal angles on middle points, JOREK coordinate system  y --> 1, x --> 3
          angle12 = atan2(-mid12(1),mid12(3))
          angle13 = atan2(-mid13(1),mid13(3))
          angle23 = atan2(-mid23(1),mid23(3))
          
          ! --- Toroidal basis vectors on middle points, in clock-wise direction looking from above the torus
          ephi12(:) = (/- cos(angle12), 0.d0, - sin(angle12) /)
          ephi13(:) = (/- cos(angle13), 0.d0, - sin(angle13) /)
          ephi23(:) = (/- cos(angle23), 0.d0, - sin(angle23) /)
         
          ! --- Quantities needed to calculate the current density vector 
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r13(:) = r3(:)-r1(:)
         
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          !--- current density vector in kA/m, Merkel 2015
          j_lin = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2))
          j_lin = j_lin / mu_zero * 1.d-3  
            
          !--- Vector normal to triangle surface
          nm(:)  = r21_cross_r32(:) / sqrt(sum(r21_cross_r32**2))
          
          !--- Poloidal vectors by triangle side (e_phi x n), (at middle points)
          pol13(:) = (/ ephi13(2)*nm(3) - ephi13(3)*nm(2), ephi13(3)*nm(1) - ephi13(1)*nm(3),          &
            ephi13(1)*nm(2) - ephi13(2)*nm(1) /)
          pol21(:) = (/ ephi12(2)*nm(3) - ephi12(3)*nm(2), ephi12(3)*nm(1) - ephi12(1)*nm(3),          &
            ephi12(1)*nm(2) - ephi12(2)*nm(1) /)
          pol32(:) = (/ ephi23(2)*nm(3) - ephi23(3)*nm(2), ephi23(3)*nm(1) - ephi23(1)*nm(3),          &
            ephi23(1)*nm(2) - ephi23(2)*nm(1) /)
          
          !--- Normalized poloidal vectors by triangle side
          pol13(:) = pol13(:)/ sqrt(sum(pol13**2))
          pol21(:) = pol21(:)/ sqrt(sum(pol21**2))
          pol32(:) = pol32(:)/ sqrt(sum(pol32**2))
          
          !--- Normal vectors to triangle sides
          n13(:) = (/ r13(2)*nm(3) - r13(3)*nm(2), r13(3)*nm(1) - r13(1)*nm(3),          &
            r13(1)*nm(2) - r13(2)*nm(1) /)
          n21(:) = (/ r21(2)*nm(3) - r21(3)*nm(2), r21(3)*nm(1) - r21(1)*nm(3),          &
            r21(1)*nm(2) - r21(2)*nm(1) /)
          n32(:) = (/ r32(2)*nm(3) - r32(3)*nm(2), r32(3)*nm(1) - r32(1)*nm(3),          &
            r32(1)*nm(2) - r32(2)*nm(1) /)
          
          !--- Normalized vectors normal to triangle sides
          n13(:) = n13(:)/ sqrt(sum(n13**2))
          n21(:) = n21(:)/ sqrt(sum(n21**2))
          n32(:) = n32(:)/ sqrt(sum(n32**2))
          
          !--- Poloidal current flowing per side of the triangle 
          j13    = abs(dot_product(n13,pol13)) * sqrt(sum(r13**2)) * dot_product(pol13,j_lin)
          j21    = abs(dot_product(n21,pol21)) * sqrt(sum(r21**2)) * dot_product(pol21,j_lin)
          j32    = abs(dot_product(n32,pol32)) * sqrt(sum(r32**2)) * dot_product(pol32,j_lin)
          
          jsides(:) = (/ j13, j21, j32 /)
              
          !--- Find the triangle side with maximum absolute current
          maxcurr_pos   = maxloc(abs(jsides),1)
          
          !--- Maximum poloidal current flowing on a triangle side, in kA
          write(filehandle,142) jsides(maxcurr_pos) 
        
        end do
      end if
      
      ! --- Linear toroidal current density per triangle (kA/m)
      if (jphi_lin) then
     
        write(filehandle,140) 'SCALARS jphi_lin(kA/m) float'
        write(filehandle,140) 'LOOKUP_TABLE default'
         
        do i = 1, sr%ntri_w
         
          ! --- Wall potentials at triangle nodes
          phi1   = tripot_w(sr%jpot_w(i,1))
          phi2   = tripot_w(sr%jpot_w(i,2))
          phi3   = tripot_w(sr%jpot_w(i,3))
          
          ! --- Positions of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          
          ! --- Middle points on triangle sides
          mid12(:) =   (r1(:) + r2(:)) * 0.5d0
          mid13(:) =   (r1(:) + r3(:)) * 0.5d0
          mid23(:) =   (r2(:) + r3(:)) * 0.5d0      
               
          ! --- Toroidal angles on middle points, JOREK coordinate system  y --> 1, x --> 3
          angle12 = atan2(-mid12(1),mid12(3))
          angle13 = atan2(-mid13(1),mid13(3))
          angle23 = atan2(-mid23(1),mid23(3))
          
          ! --- Toroidal basis vectors on middle points, in clock-wise direction looking from above the torus
          ephi12(:) = (/- cos(angle12), 0.d0, - sin(angle12) /)
          ephi13(:) = (/- cos(angle13), 0.d0, - sin(angle13) /)
          ephi23(:) = (/- cos(angle23), 0.d0, - sin(angle23) /)
         
          ! --- Quantities needed to calculate the current density vector 
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r13(:) = r3(:)-r1(:)
         
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          !--- current density vector in kA/m, Merkel 2015
          j_lin = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2))
          j_lin = j_lin / mu_zero * 1.d-3  
          
          !--- Toroidal linear density current in each triangle side
          j13    =  dot_product(ephi13,j_lin)
          j21    =  dot_product(ephi12,j_lin)
          j32    =  dot_product(ephi23,j_lin)
          
          jsides(:) = (/ j13, j21, j32 /)
              
          !--- Find the triangle side with maximum current
          maxcurr_pos   = maxloc(abs(jsides),1)
          
          !--- Maximum toroidal linear density current flowing on a triangle side, in kA/m
          write(filehandle,142) jsides(maxcurr_pos) 
        
        end do
      end if
      
      ! --- Poloidal linear density current flowing on a triangle (kA/m)
      if (jpol_lin) then
     
        write(filehandle,140) 'SCALARS jpol_lin(kA/m) float'
        write(filehandle,140) 'LOOKUP_TABLE default'
         
        do i = 1, sr%ntri_w
         
          ! --- Wall potentials at triangle nodes
          phi1   = tripot_w(sr%jpot_w(i,1))
          phi2   = tripot_w(sr%jpot_w(i,2))
          phi3   = tripot_w(sr%jpot_w(i,3))
          
          ! --- Positions of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          
          ! --- Middle points on triangle sides
          mid12(:) =   (r1(:) + r2(:)) * 0.5d0
          mid13(:) =   (r1(:) + r3(:)) * 0.5d0
          mid23(:) =   (r2(:) + r3(:)) * 0.5d0      
               
          ! --- Toroidal angles on middle points, JOREK coordinate system  y --> 1, x --> 3
          angle12 = atan2(-mid12(1),mid12(3))
          angle13 = atan2(-mid13(1),mid13(3))
          angle23 = atan2(-mid23(1),mid23(3))
          
          ! --- Toroidal basis vectors on middle points, in clock-wise direction looking from above the torus
          ephi12(:) = (/- cos(angle12), 0.d0, - sin(angle12) /)
          ephi13(:) = (/- cos(angle13), 0.d0, - sin(angle13) /)
          ephi23(:) = (/- cos(angle23), 0.d0, - sin(angle23) /)
         
          ! --- Quantities needed to calculate the current density vector 
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r13(:) = r3(:)-r1(:)
         
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          !--- current density vector in kA/m, Merkel 2015
          j_lin = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2))
          j_lin = j_lin / mu_zero * 1.d-3  
            
          !--- Vector normal to triangle surface
          nm(:)  = r21_cross_r32(:) / sqrt(sum(r21_cross_r32**2))
          
          !--- Poloidal vectors by triangle side (e_phi x n), (at middle points)
          pol13(:) = (/ ephi13(2)*nm(3) - ephi13(3)*nm(2), ephi13(3)*nm(1) - ephi13(1)*nm(3),          &
            ephi13(1)*nm(2) - ephi13(2)*nm(1) /)
          pol21(:) = (/ ephi12(2)*nm(3) - ephi12(3)*nm(2), ephi12(3)*nm(1) - ephi12(1)*nm(3),          &
            ephi12(1)*nm(2) - ephi12(2)*nm(1) /)
          pol32(:) = (/ ephi23(2)*nm(3) - ephi23(3)*nm(2), ephi23(3)*nm(1) - ephi23(1)*nm(3),          &
            ephi23(1)*nm(2) - ephi23(2)*nm(1) /)
          
          !--- Normalized poloidal vectors by triangle side
          pol13(:) = pol13(:)/ sqrt(sum(pol13**2))
          pol21(:) = pol21(:)/ sqrt(sum(pol21**2))
          pol32(:) = pol32(:)/ sqrt(sum(pol32**2))
          
          !--- Poloidal current flowing per side of the triangle 
          j13    =  dot_product(pol13,j_lin)
          j21    =  dot_product(pol21,j_lin)
          j32    =  dot_product(pol32,j_lin)
          
          jsides(:) = (/ j13, j21, j32 /)
              
          !--- Find the triangle side with maximum absolute current
          maxcurr_pos   = maxloc(abs(jsides),1)
          
          !--- Maximum poloidal linear density current flowing on a triangle, in kA/m
          write(filehandle,142) jsides(maxcurr_pos) 
        
        end do
      end if
      
      ! --- Current density: Contribution of single valued potentials  
      write(filehandle,140) 'VECTORS jw_single_val(MA/m) float'

      do i = 1, sr%ntri_w
        ! --- Wall potential at triangle nodes
        phi1   = tripot_w(sr%jpot_w(i,1))
        phi2   = tripot_w(sr%jpot_w(i,2))
        phi3   = tripot_w(sr%jpot_w(i,3))
        ! --- Position of triangle nodes
        r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
        r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
        r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
        r21(:) = r1(:)-r2(:)
        r32(:) = r2(:)-r3(:)
        r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
          r21(1)*r32(2) - r21(2)*r32(1) /)
          
        ! Exports the linear wall density current in MA/m  
        write(filehandle,142) ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2)) &
                              / mu_zero * 1.d-6
      end do

      if (sr%file_version >= 5) then

        ! --- Current density: Contribution of net wall potential 
        write(filehandle,140) 'VECTORS jw_net(MA/m) float'
  
        do i = 1, sr%ntri_w
          ! --- Wall potential at triangle nodes
          phi1   = Iw_net_tor * sr%phi0_w(i,1) 
          phi2   = Iw_net_tor * sr%phi0_w(i,2) 
          phi3   = Iw_net_tor * sr%phi0_w(i,3)
          ! --- Position of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          ! Exports the linear wall density current in MA/m  
          write(filehandle,142) ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2)) &
                                / mu_zero * 1.d-6
        end do

        ! --- Total wall current density
        write(filehandle,140) 'VECTORS jw_tot(MA/m) float'
  
        do i = 1, sr%ntri_w
          ! --- Wall potential at triangle nodes
          phi1   = Iw_net_tor * sr%phi0_w(i,1) + tripot_w(sr%jpot_w(i,1)) 
          phi2   = Iw_net_tor * sr%phi0_w(i,2) + tripot_w(sr%jpot_w(i,2))
          phi3   = Iw_net_tor * sr%phi0_w(i,3) + tripot_w(sr%jpot_w(i,3))
          ! --- Position of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          ! Exports the linear wall density current in MA/m  
          write(filehandle,142) ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2)) &
                                / mu_zero * 1.d-6
        end do
      endif    ! file version 5

      
      ! --- Close file, clean up
      close(filehandle)
  
    end if  ! my_id==0

    if ( allocated(tripot_w) ) deallocate( tripot_w )

  end subroutine write_wall_vtk
  
  

  !> Implements the boundary integral in the current equation which vanishes for fixed boundary
  !! conditions.
  !!
  !! The magnetic field parallel to the interface (boundary of JOREK computational domain) is
  !! expressed by the STARWALL vacuum response in terms of the poloidal magnetic field at the
  !! interface (ideal and resistive wall) and the wall currents (resistive wall).
  subroutine vacuum_boundary_integral(my_id, bnd_node_list, node_list,bnd_elm_list, element_list,             &
    freeboundary_equil, resistive_wall, index_min, index_max, rhs_loc, tstep, index_now, a_mat)

    use data_structure, only: type_node_list, type_node, type_element_list, type_element, type_bnd_node_list,type_bnd_element_list, type_bnd_element, type_SP_MATRIX
    use mod_parameters, only: n_plane, n_var, n_tor, n_degrees
    use gauss,          only: n_gauss, xgauss, wgauss
    use global_distributed_matrix, only: det_row_col, det_sparse_pos
    use basis_at_gaussian, only: H1, H1_s, HZ
    use mod_basisfunctions
    use phys_module, only: t_now, t_start, PI, n_period, F0
    use mpi_mod
    use mod_settings
    use mod_integer_types
    use constants, only: MU_ZERO

    implicit none

    ! --- Routine parameters
    integer,                            intent(in)    :: my_id                !< MPI process ID
    type(type_node_list),               intent(in)    :: node_list            !< List of grid nodes
    type(type_bnd_node_list),           intent(in)    :: bnd_node_list        !< List of boundary grid nodes
    type(type_bnd_element_list),        intent(in)    :: bnd_elm_list         !< List of boundary elements
    type(type_element_list),            intent(in)    :: element_list
    logical,                            intent(in)    :: freeboundary_equil   !< Use free boundary equilibrium?
    logical,                            intent(in)    :: resistive_wall       !< Resistive or ideal wall?
    integer,                            intent(in)    :: index_min, index_max !< Responsibility of MPI proc
    real*8,                             intent(inout) :: rhs_loc(:)           !< Part of RHS of MPI proc 
    real*8,                             intent(in)    :: tstep                !< delta t, timestep
    integer,                            intent(in)    :: index_now            !< Current timestep index
    type(type_SP_MATRIX)                              :: a_mat

    ! --- Local variables
    real*8, allocatable :: psibnd_vec(:)    ! Vector of the values of Psi at the boundary
    real*8, allocatable :: dpsibnd_vec(:)   ! Vector of the values of deltaPsi at the boundary
    real*8, allocatable :: dzjbnd_vec(:)
    real*8, allocatable :: psibnd_coils(:)  ! Vector of the values of Psi_coil at the boundary
    real*8   :: amat_contrib, rhs_contrib   ! Vacuum response contribution to lhs and rhs
    real*8   :: amat_contrib_psi, amat_contrib_zj
    real*8   :: testfunc_l                  ! j^*_l in documentation
    real*8   :: basfunc_i                   ! b_i in documentation
    real*8   :: dA                          ! factor from definition of dA
    real*8   :: x_s(n_gauss), y_s(n_gauss)  ! values of dR/ds and dZ/ds at Gaussian points
    real*8   :: common_prefactor
    integer  :: m_bndelem, j_bndelem, s_bndelem , e_bndelem                 ! Boundary element index
    type(type_bnd_element) :: bndelem_m, bndelem_j, bndelem_e     ! Boundary element corresponding to index m_bndelem
    integer  :: ms                          ! Gauss point index
    integer  :: m_plane                     ! Toroidal plane index
    integer  :: sparsepos_jp, sparsepos_pp, sparsepos_jj  ! Position of lhs contribution in the sparse matrix
    integer  :: blockpos_jp, blockpos_pp    ! Position of respective block in sparse matrix
    !   --- Test function related quantities
    integer  :: l_vertex, l_dof, l_dir, l_node, l_node_bnd, l_index, l_tor, l_row_j, l_row_psi
    real*8   :: l_size
    !   --- Quantities related to the boundary dof at which response is calculated
    integer  :: i_vertex, i_dof, i_dir, i_node, i_node_bnd, i_index, i_starwall, i_tor, i_resp, i_resp_0
    real*8   :: i_size
    !   --- Quantities related to the boundary dof contributing to the response
    integer  :: j_dof, j_dir, j_node, s_node, j_node_bnd, j_index, j_starwall, j_tor, j_col_psi, j_col_zj, j_resp, j_vertex
    integer  :: s_vertex, s_dof, s_dir, s_tor
    real*8   :: s_size, j_size, fact_mu0, j_size1, j_size2
    !   --- Quantities for halo coupling implicit time stepping
    integer  :: e_node, e_node_bnd, e_dof, e_dir, e_index, e_vertex, e_harm, e_tor, real_index, e_resp
    real*8   :: e_size

    integer  :: i_resp_old

    integer  :: i_start_pf, i_end_pf     ! Indices for PF coils
    integer  :: npts_tot, jj , mp, npol,j
    real*8, allocatable :: identity_mat(:,:), real2bez(:,:)
    real*8   :: phi_vec(n_plane)
    real*8   :: sign_out, R_c, Z_c, R(n_gauss), Z(n_gauss), R_s(n_gauss), Z_s(n_gauss), vec_inside(2), grad_t(2), G(4, n_degrees)
    real*8   :: psi_s(n_gauss*n_plane), zj(n_gauss*n_plane), czj, cpsi
    real*8, allocatable :: response_R_zj(:,:), response_R_psis(:,:), J_normal(:), jn_bez(:), test_sig(:), bez_sig(:), phi_t
    real*8, allocatable :: response_m_t_h_zj(:,:) , response_m_t_h_psis(:,:)
    real*8, ALLOCATABLE :: contrib_iwe(:), contrib_ieleq(:), contrib_yi(:), contrib_dyi(:)
    type(type_bnd_element) :: bndelem_s
    type(type_element)     :: elem_s    
    type(type_node)        :: node_k, node_s
   
#ifdef __GFORTRAN__
    real*8 :: wgauss_copy(4)
#endif
    real*8 :: t_elaps_start, t_elaps_end !### timing ###
    logical, save  :: PF_perturbation = .true.
    logical  :: debug_field_contrib
    integer  :: ierr,i
    real*8, allocatable :: rhs_contrib_arr(:)
    real*8, allocatable :: diag_1(:), diag_e(:), diag_el_eq(:)

    if ( sr%n_tor == 0 ) then
      write(*,*) 'Skipping vacuum_boundary_integral since sr%n_tor==0.'
      return
    end if
    debug_field_contrib = .true.
    if ( vacuum_debug ) t_elaps_start = MPI_WTIME()  ! for timing

    if ( vacuum_debug ) write(*,*) my_id, 'Before:', sum(abs(rhs_loc)),sum(abs(a_mat%val))

    if (halo_freeb) then
    !!! R.S. > SHIFT ALL OF THIS TO A SUBROUTINE: MAKE A VERSION WHERE YOU ASSIGN DELTAZJ AND DELTAPSI AS
    ! OPTIONAL ARGUMENTS AND A VERSION FOR THE LHS ENTRIES WHERE YOU DON'T ASSIGN THEM

      ! if (.not. ALLOCATED(response_R_psis)) allocate(response_R_psis(sr%nd_bez, sr%nd_bez))
      ! if (.not. ALLOCATED(response_R_zj)) allocate(response_R_zj(sr%nd_bez, sr%nd_bez))
      ! if (.not. ALLOCATED(jn_bez)) allocate(jn_bez(sr%nd_bez))
      ! if (.not. ALLOCATED(response_m_t_h_zj)) allocate(response_m_t_h_zj(sr%nd_bez, sr%nd_bez))
      ! if (.not. ALLOCATED(response_m_t_h_psis)) allocate(response_m_t_h_psis(sr%nd_bez, sr%nd_bez))

      !call halo_time_stepping(my_id, node_list, element_list, bnd_elm_list, bnd_node_list, response_R_zj=response_R_zj, response_R_psis=response_R_psis, jn_bez=jn_bez, debug=.false.)
      
      ! Perform products with electrode projection matrices
      ! response_m_t_h_zj(:,:) = 0.d0 
      ! response_m_t_h_psis(:,:) = 0.d0
      ! response_m_t_h_zj = MATMUL(response_m_t_h, response_R_zj)
      ! response_m_t_h_psis = MATMUL(response_m_t_h, response_R_psis)

      ! if (.not. ALLOCATED(R_psis)) allocate(R_psis(npts_tot, sr%nd_bez))
      ! if (.not. ALLOCATED(R_zj)) allocate(R_zj(npts_tot, sr%nd_bez ))
      ! if (.not. ALLOCATED(fact_R_psis)) allocate(fact_R_psis(npts_tot))
      ! if (.not. ALLOCATED(fact_R_zj)) allocate(fact_R_zj(npts_tot))
      ! if (.not. ALLOCATED(J_normal)) allocate(J_normal(npts_tot))
      ! 
      ! if (.not. ALLOCATED(test_sig)) ALLOCATE(test_sig(npts_tot))
      ! if (.not. ALLOCATED(bez_sig)) ALLOCATE(bez_sig(n_dof_starwall))
      
      ! npol     = bnd_elm_list%n_bnd_elements * n_gauss
      ! if (my_id.eq. 0) write(*,*) 'NPOL', npol
      ! if (.not. allocated(identity_mat)) allocate(identity_mat(Npts_tot, Npts_tot))

      ! identity_mat = 0.d0
      ! do jj = 1 , Npts_tot
      !   identity_mat(jj,jj) = 1.d0
      ! end do
      ! do mp = 1, n_plane
      !   phi_vec   (mp)      = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
      ! enddo
      ! ! N.B. Real2Bez WILL BE ALLOCATED BASED ON THE NUMBER OF FREE BOUNDARY MODES IN STARWALL!    
      ! call fourier(identity_mat, Npts_tot, n_plane, n_gauss, Npts_tot, phi_vec, 'gauss', Real2Bez, my_id, bnd_node_list, bnd_elm_list)
      ! ! ! --- Write to file ---
      ! ! if (my_id .eq. 0) write(*,*) 'Writing real2bez.dat'
      ! ! open(unit=12, file="Real2Bez.dat", status="replace", action="write", form="formatted")
      ! ! do i = 1,n_dof_starwall
      ! !   write(12, '(100(1x,ES20.12))') (real2bez(i,j), j=1,npts_tot)
      ! ! end do
      ! ! close(12)
      ! ! if (my_id .eq. 0) write(*,*) 'Done writing real2bez.dat'

      ! ! test_sig = 0.d0
      ! ! bez_sig  = 0.d0
      ! ! ! Define test_sig
      ! ! do mp = 1,n_plane
      ! !   phi_t = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
      ! !   test_sig((mp-1)*npol+1:mp*npol) = sin(phi_t)
      ! ! enddo
      ! ! do i = 1,n_dof_starwall
      ! !   bez_sig(i) = sum(real2bez(i,:)*test_sig(:))
      ! ! enddo

      ! ! if (my_id .eq. 0) write(*,*) 'Writing test_sig.dat'
      ! ! open(unit=12, file="test_sig.dat", status="replace", action="write", form="formatted")
      ! ! do i = 1,npts_tot
      ! !   write(12,*) test_sig(i)
      ! ! end do
      ! ! close(12)
      ! ! if (my_id .eq. 0) write(*,*) 'Done writing test_sig.dat'  

      ! ! if (my_id .eq. 0) write(*,*) 'Writing bez_sig.dat'
      ! ! open(unit=12, file="bez_sig.dat", status="replace", action="write", form="formatted")
      ! ! do i = 1,n_dof_starwall
      ! !   write(12,*) bez_sig(i)
      ! ! end do
      ! ! close(12)
      ! ! if (my_id .eq. 0) write(*,*) 'Done writing bez_sig.dat'      

      ! ! Calculate real space coefficients for implicit time stepping
      ! fact_R_psis = 0.d0
      ! fact_R_zj   = 0.d0 
      ! J_normal    = 0.d0
      ! jn_bez      = 0.d0

      ! do s_bndelem = 1,bnd_elm_list%n_bnd_elements
      !   bndelem_s = bnd_elm_list%bnd_element(s_bndelem)
      !   elem_s    = element_list%element(bndelem_s%element)

      !   psi_s = 0.d0
      !   zj    = 0.d0
      !   R     = 0.d0
      !   Z     = 0.d0
      !   R_s   = 0.d0
      !   Z_s   = 0.d0

      !   do s_vertex = 1, 2 ! (loop over nodes in element s_bndelem)
      !     do s_dof = 1, 2 ! (loop over node dofs)
      !       s_node      = bndelem_s%vertex(s_vertex)
      !       node_s      = node_list%node(s_node)
      !       s_dir       = bndelem_s%direction(s_vertex,s_dof)
      !       s_size      = bndelem_s%size(s_vertex,s_dof)
            
      !       R  (:)   = R  (:)  + node_s%x(1,s_dir,1) * s_size * H1  (s_vertex,s_dof,:)
      !       Z  (:)   = Z  (:)  + node_s%x(1,s_dir,2) * s_size * H1  (s_vertex,s_dof,:)
      !       R_s(:)   = R_s(:)  + node_s%x(1,s_dir,1) * s_size * H1_s(s_vertex,s_dof,:)
      !       Z_s(:)   = Z_s(:)  + node_s%x(1,s_dir,2) * s_size * H1_s(s_vertex,s_dof,:)

      !       do mp=1,n_plane
      !         do s_tor = 1, n_tor ! (loop over toroidal harmonics)
      !           ! Chunks of size n_gauss * n_plane
      !           psi_s((mp-1)*n_gauss+1:mp*n_gauss) =    psi_s((mp-1)*n_gauss+1:mp*n_gauss) + &
      !           node_s%values(s_tor,s_dir,var_psi) * s_size*H1_s(s_vertex,s_dof,:) * HZ(s_tor,mp)

      !           zj((mp-1)*n_gauss+1:mp*n_gauss)    =    zj((mp-1)*n_gauss+1:mp*n_gauss) + &
      !           node_s%values(s_tor,s_dir,var_zj) * s_size*H1(s_vertex,s_dof,:) * HZ(s_tor,mp)
      !         enddo            
      !       enddo

      !     enddo
      !   enddo

      !   !---------------------- BOUNDARY ELEMENT ORIENTATION-------------------------- 
      !   ! --- Calculate an inside point on the element to then calculate the
      !   ! direction of bnd normals

      !   vec_inside(2) = 0.d0
      !   grad_t(2)     = 0.d0
      !   sign_out   = 0.d0

      !   call basisfunctions(xgauss(2),xgauss(2), G)  
      !   R_c = 0.d0 ;  Z_c = 0.d0 
      !   do i = 1, n_vertex_max
      !     do j = 1, n_degrees
      !       node_k = node_list%node(elem_s%vertex(i)) 
      !       R_c    = R_c + node_k%x(1,j,1) * elem_s%size(i,j) * G(i,j)
      !       Z_c    = Z_c + node_k%x(1,j,2) * elem_s%size(i,j) * G(i,j)
      !     enddo
      !   enddo  
      !   vec_inside = (/ R_c - R(2), Z_c - Z(2) /)       ! vector pointing towards the domain
      !   !--- Find out correct sign of the normal (it has to point outwards the domain)
      !   grad_t     = (/ -Z_s(2) , R_s(2) /)     ! gradient of the coordinate t (normal to the boundary here) 
      !   ! grad_t could point either outwards or inwards depending on the orientation of the bnd_element!
      !   sign_out   = -1.d0 * sign( 1.d0, ( vec_inside(1)*grad_t(1) + vec_inside(2)*grad_t(2) ) )  ! (+1 if normal outwards)

      !   do mp = 1, n_plane

      !     fact_R_zj ((s_bndelem-1)*n_gauss+1+(mp-1)*npol : s_bndelem*n_gauss+(mp-1)*npol)  = &
      !       psi_s((mp-1)*n_gauss+1:mp*n_gauss) &
      !       / sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out    

      !     fact_R_psis ((s_bndelem-1)*n_gauss+1+(mp-1)*npol : s_bndelem*n_gauss+(mp-1)*npol)  = &
      !         zj((mp-1)*n_gauss+1:mp*n_gauss) &
      !         / sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out  
                              
      !     J_normal ((s_bndelem-1)*n_gauss+1+(mp-1)*npol : s_bndelem*n_gauss+(mp-1)*npol) = &
      !       ( ( zj((mp-1)*n_gauss+1:mp*n_gauss) * psi_s((mp-1)*n_gauss+1:mp*n_gauss)) ) &
      !                                 /sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out
      !   enddo
      ! enddo
      

      ! ! Calculate the coefficients' matrix
      ! fact_mu0 = 1.d0/mu_zero
      ! R_zj = 0.d0 
      ! R_psis = 0.d0
      ! response_R_zj = 0.d0 
      ! response_R_psis = 0.d0

      ! do e_bndelem = 1,bnd_elm_list%n_bnd_elements 
      !   ! Select bnd element
      !   bndelem_e = bnd_elm_list%bnd_element(e_bndelem)
      !   do e_vertex = 1,2
      !     e_node      = bndelem_e%vertex(e_vertex)
      !     e_node_bnd  = bndelem_e%bnd_vertex(e_vertex)
      !     do e_dof = 1,2
      !       e_dir       = bndelem_e%direction(e_vertex,e_dof)
      !       e_index     = node_list%node(e_node)%index(e_dir)
      !       e_size      = bndelem_e%size(e_vertex,e_dof)
      !       do e_harm = 1, sr%n_tor
      !         e_tor  = sr%i_tor(e_harm)
      !         ! Loop over real pts
      !         do s_bndelem = 1,bnd_elm_list%n_bnd_elements
      !           do mp=1,n_plane
      !             do ms=1,n_gauss
      !               ! Build matrices
      !               real_index = n_gauss*(s_bndelem-1)+ ms + (mp-1)*npol
      !               ! First index to switch is the node dofs, then harmonics, the nodes
      !               e_resp   = (bnd_node_list%bnd_node(e_node_bnd)%index_starwall(1) - 1)*sr%n_tor0 &
      !                            + bnd_node_list%bnd_node(e_node_bnd)%n_dof*(e_harm-1)              &
      !                            + bnd_node_list%bnd_node(e_node_bnd)%index_starwall(e_dof)             &
      !                            - bnd_node_list%bnd_node(e_node_bnd)%index_starwall(1) + 1
                
      !               R_zj( real_index, e_resp ) = &
      !                     fact_mu0*fact_R_zj( real_index ) &
      !                   * e_size * H1(e_vertex,e_dof,ms) * HZ(e_tor,mp)                  

      !               R_psis( real_index, e_resp ) = &
      !                     fact_mu0*fact_R_psis( real_index ) &
      !                   * e_size * H1_s(e_vertex,e_dof,ms) * HZ(e_tor,mp)

      !             enddo 
      !           enddo 
      !         enddo 
      !       enddo 
      !     enddo 
      !   enddo 
      ! enddo

      ! ! Convert the matrices to Bezier space
      ! if (my_id .eq. 0) write(*,*) 'About to perform MATMUL from vacuum_boundary_integral...'
      ! response_R_zj = MATMUL(real2bez,R_zj) 
      ! response_R_psis = MATMUL(real2bez, R_psis )

      !   ! if (my_id.eq.0) then 
      !   !   open(unit=12, file="response_R_zj_0.dat", status="replace", action="write", form="formatted")
      !   !   do i=1,sr%nd_bez
      !   !     write(12,*) (response_R_zj(i,j), j=1,sr%nd_bez)
      !   !   enddo 
      !   !   close(12) 

      !   !   open(unit=12, file="response_R_psis_0.dat", status="replace", action="write", form="formatted")
      !   !   do i=1,sr%nd_bez
      !   !     write(12,*) (response_R_psis(i,j), j=1,sr%nd_bez)
      !   !   enddo 
      !   !   close(12)
      !   ! endif   
      !   ! if (my_id.eq.1) then 
      !   !   open(unit=12, file="response_R_zj_1.dat", status="replace", action="write", form="formatted")
      !   !   do i=1,sr%nd_bez
      !   !     write(12,*) (response_R_zj(i,j), j=1,sr%nd_bez)
      !   !   enddo 
      !   !   close(12) 

      !   !   open(unit=12, file="response_R_psis_1.dat", status="replace", action="write", form="formatted")
      !   !   do i=1,sr%nd_bez
      !   !     write(12,*) (response_R_psis(i,j), j=1,sr%nd_bez)
      !   !   enddo 
      !   !   close(12)
      !   ! endif        

      ! if (my_id .eq. 0) write(*,*) 'Done performing MATMUL from vacuum_boundary_integral.'

      ! RESCALE J_normal for SI units! 

      ! J_normal = J_normal * fact_mu0
      ! ! Calculate J_normal in Bezier
      ! do jj = 1,sr%nd_bez
      !   Jn_bez(jj) = sum(real2bez(jj,:) * J_normal(:))
      ! enddo

      ! if (my_id .eq. 0) write(*,*) 'Writing J_normal.dat'
      ! open(unit=12, file="J_normal.dat", status="replace", action="write", form="formatted")
      ! do i = 1,npts_tot
      !   write(12,*) J_normal(i)
      ! end do
      ! close(12)
      ! if (my_id .eq. 0) write(*,*) 'Done writing J_normal.dat'   

      !       if (my_id .eq. 0) write(*,*) 'Writing J_bez.dat'
      ! open(unit=12, file="J_bez.dat", status="replace", action="write", form="formatted")
      ! do i = 1,n_dof_starwall
      !   write(12,*) Jn_bez(i)
      ! end do
      ! close(12)
      ! if (my_id .eq. 0) write(*,*) 'Done writing Jn_bez.dat'   

    endif  ! endif halo

    ! --- Determine vectors of the psi and deltapsi boundary values.
    if (halo_freeb) then 
      call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec, dzjbnd_vec,  psibnd_coils)
    else 
      call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec,  psibnd_coils)
    endif

    ! --- Calculate current source term for PF coil currents imposition
    if ( (sr%ncoil /= 0) ) then    !avoid it when using COIL_FIELD
      call coil_current_source(my_id)
    else
      if (halo_freeb .and. CARIDDI_mode) then
        if (.not. allocated (Yi_coils0)) allocate(Yi_coils0(n_wall_curr_i))
        Yi_coils0(:) = 0.d0
      else
        if (.not. allocated (Y_coils0)) allocate(Y_coils0(n_wall_curr))
        Y_coils0(:) = 0.d0
      endif
    endif

    ! --- Update the derived response matrices
    call update_response(my_id, tstep, resistive_wall, node_list, element_list, bnd_elm_list, bnd_node_list )

    ! --- Perform the time-stepping for the wall currents.
    if ( resistive_wall .and.  (sr%n_tor>0) ) then 
      if (.not. halo_freeb) call evolve_wall_currents(my_id, psibnd_vec, dpsibnd_vec, A_MAT=a_mat)
      if (halo_freeb) then 
        call evolve_wall_currents(my_id, psibnd_vec, dpsibnd_vec, dzjbnd_vec, a_mat)
      endif
    endif

    ! --- Update old_dpsibnd_vec  (used for updating the wall currents in the next iteration)
    old_dpsibnd_vec(:) = dpsibnd_vec(:)
    if (halo_freeb) then 
      old_dIelec_J(:)   = dIelec_J(:)
      old_dIelec_w(:)   = dIelec_w(:)
    endif

    if ( vacuum_debug ) then
      write(*,*) my_id, 'psibnd_vec:  ', sum(abs(psibnd_vec)), sum(psibnd_vec)
      write(*,*) my_id, 'dpsibnd_vec: ', sum(abs(dpsibnd_vec)), sum(dpsibnd_vec)
    end if

    call boundary_check(my_id)
    if (.not. allocated(rhs_contrib_arr)) allocate(rhs_contrib_arr(n_dof_starwall))
    rhs_contrib_arr=0.0

    if (resistive_wall) then
    if (my_id .eq. 0) write(*,*) 'Performing halo RHS products'
    if (my_id .eq. 0) write(*,*) 'n_shell_curr_e=', n_shell_curr_e
    if (my_id .eq. 0) write(*,*) 'n_wall_curr_e=', n_wall_curr_e

      ! ===============================  Btan (eddy + halo) ===================================== 
      if (CARIDDI_mode .and. halo_freeb) then

        if (.not. allocated(diag_1)) allocate(diag_1(n_wall_curr_i))
        diag_1(:) = ( 1.d0 + response_d_b_i(:) )
        !> corresponds to matrix F_i in documentation
        call extended_matrix_vector(my_id,  1.0d0, sr%a_eyi, diag_1, &
            wall_curr_i(:), .true., rhs_contrib_arr, .false.)
        !> corresponds to matrix G_i in documentation
        call extended_matrix_vector(my_id,  1.0d0, sr%a_eyi, response_d_c_i, &
            dwall_curr_i(:), .true., rhs_contrib_arr, .false.)
        
        !!!!!!!!!!!!!!!!!!!!! EXPLICIT SCHEME TERMS !!!!!!!!!!!!!!!!!!!!!!
        if (.not. halo_imp) then
          !> corresponds to matrix V_i in documentation
          call extended_matrix_vector(my_id, -1.0d0, sr%a_eyi, response_d_b_i, &
              Yi_coils0(:), .true., rhs_contrib_arr, .false.)
          ! !> corresponds to matrix F_e in documentation
          ! if (.not. allocated(diag_e) ) allocate(diag_e(n_wall_curr_e))
          ! diag_e(:) = 1.d0
          ! call extended_matrix_vector(my_id, 1.d0, sr%a_eye, diag_e, &
          !       wall_curr_e(:), .true., rhs_contrib_arr, .false.)
          ! !> corresponds to matrix Qel_eq
          ! if (.not. allocated(diag_el_eq) ) allocate(diag_el_eq(n_shell_curr_e))
          ! diag_el_eq(:) = 1.d0
          ! call extended_matrix_vector(my_id, -1.d0, sr%aeq_el, diag_el_eq, &
          !       I_elec_J(:), .TRUE., rhs_contrib_arr, .true.)
        else
          !> corresponds to matrix V_i in documentation
          call extended_matrix_vector(my_id, -1.0d0, sr%a_eyi, response_d_b_i, &
              Yi_coils0(:), .true., rhs_contrib_arr, .true.)
        endif    

        if (debug_field_contrib) then   
          ! CALCULATE THE SEPARATE CONTRIBUTIONS
          if (.not. allocated(contrib_iwe)) allocate(contrib_iwe(n_dof_starwall))
          contrib_iwe = 0.d0 
          if (.not. allocated(diag_e) ) allocate(diag_e(n_wall_curr_e))
          diag_e(:) = 1.d0
          call extended_matrix_vector(my_id, 1.d0, sr%a_eye, diag_e, &
                wall_curr_e(:), .true., contrib_iwe, .true.)

          if (my_id .eq. 0) write(*,*) 'Writing contrib_iwe_Bez.dat'
          open(unit=12, file="contrib_iwe_bez.dat", status="replace", action="write", form="formatted")
          do i = 1,n_dof_starwall
            write(12,*) contrib_iwe(i)
          end do
          close(12)

          if (.not. allocated(contrib_ieleq)) allocate(contrib_ieleq(n_dof_starwall))
          contrib_ieleq = 0.d0 
          if (.not. allocated(diag_el_eq) ) allocate(diag_el_eq(n_shell_curr_e))
          diag_el_eq(:) = 1.d0
          call extended_matrix_vector(my_id, -1.d0, sr%aeq_el, diag_el_eq, &
                Ielec_J(:), .TRUE., contrib_ieleq, .true.)

          if (my_id .eq. 0) write(*,*) 'Writing contrib_ieleq_bez.dat'
          open(unit=12, file="contrib_ieleq_bez.dat", status="replace", action="write", form="formatted")
          do i = 1,n_dof_starwall
            write(12,*) contrib_ieleq(i)
          end do
          close(12)
          
          if (.not. allocated(contrib_yi) ) allocate(contrib_yi(n_dof_starwall))
          contrib_yi=0.d0
          if (.not. allocated(diag_1)) allocate(diag_1(n_wall_curr_i))
          diag_1(:) = ( 1.d0 + response_d_b_i(:) )
          call extended_matrix_vector(my_id,  1.0d0, sr%a_eyi, diag_1, &
            wall_curr_i(:), .true., contrib_yi, .true.)

          if (my_id .eq. 0) write(*,*) 'Writing contrib_yi_Bez.dat'
          open(unit=12, file="contrib_yi_bez.dat", status="replace", action="write", form="formatted")
          do i = 1,n_dof_starwall
            write(12,*) contrib_yi(i)
          end do
          close(12)

          if (.not. allocated(contrib_dyi) ) allocate(contrib_dyi(n_dof_starwall))
          contrib_dyi=0.d0
          call extended_matrix_vector(my_id,  1.0d0, sr%a_eyi, response_d_c_i, &
            dwall_curr_i(:), .true., contrib_dyi, .true.)
          if (my_id .eq. 0) write(*,*) 'Writing contrib_dyi_Bez.dat'
          open(unit=12, file="contrib_dyi_bez.dat", status="replace", action="write", form="formatted")
          do i = 1,n_dof_starwall
            write(12,*) contrib_dyi(i)
          end do
          close(12)

        endif
    
      else
      ! ===============================  Btan (eddy) ===================================== 

        if (.not. allocated(diag_1)) allocate(diag_1(n_wall_curr))
        diag_1(:) = ( 1.d0 + response_d_b(:) )
        !> corresponds to matrix F in documentation
        call extended_matrix_vector(my_id,  1.0d0, sr%a_ey, diag_1, &
            wall_curr(:), .true., rhs_contrib_arr, .false.)
        !> corresponds to matrix G in documentation
        call extended_matrix_vector(my_id,  1.0d0, sr%a_ey, response_d_c, &
            dwall_curr(:), .true., rhs_contrib_arr, .false.)
        !> corresponds to matrix V in documentation
        call extended_matrix_vector(my_id, -1.0d0, sr%a_ey, response_d_b, &
            Y_coils0(:), .true., rhs_contrib_arr, .true.)

      endif
    end if
  
   
    if (my_id.eq.0) write(*,*) 'About to enter bnd integral do loops'

#ifdef __GFORTRAN__
    wgauss_copy(1:4) = wgauss(1:4)
#endif
    ! --- Sum over boundary elements
    !$omp parallel do                                                                                         &
    !$omp default(none)                                                                                       &    
    !$omp shared(my_id, a_mat, rhs_loc, bnd_elm_list, element_list, bnd_node_list, node_list, index_min, index_max,          &
    !$omp   response_m_e, response_m_e_i, response_m_h, response_m_j, response_m_j_i, response_m_eh, response_m_zh, H1, H1_s, HZ, sr,              &
    !$omp   response_m_ph, response_m_qh, response_m_xh, response_m_vh,                                                                        &
    !$omp   bext_tan, I_coils, wall_curr, dwall_curr, psibnd_vec, dpsibnd_vec,psibnd_coils,                   &
    !$omp  wall_curr_i, dwall_curr_i, wall_curr_e, dwall_curr_e, halo_imp, F0, npol, npts_tot, J_normal, fact_mu0,                   &
    !$omp Ielec_w, Ielec_J, dIelec_w, dIelec_J, &
#ifdef __GFORTRAN__
    !$omp   wgauss_copy,                                                                                      &
#endif
    !$omp   starwall_equil_coils,  Y_coils0, Yi_coils0, rhs_contrib_arr, resistive_wall, CARIDDI_mode, halo_freeb) &
    !$omp private(m_bndelem, j_bndelem,  s_bndelem, bndelem_m, elem_s, bndelem_s, node_k, node_s, bndelem_j,  x_s, y_s, l_vertex, l_dof, l_node,l_dir, l_node_bnd,                  &
    !$omp   l_index, l_size, l_tor, l_row_j, l_row_psi, ms, dA, m_plane,common_prefactor, s_vertex, s_dof,                   &
    !$omp   testfunc_l, i_vertex, i_dof, i_node, i_dir, i_node_bnd, i_index, i_size, s_size, i_starwall,              &
    !$omp   mp,R_c, Z_c, R,Z,R_s,Z_s, G,   &
    !$omp   i_tor, i_resp, i_resp_old, i_resp_0, basfunc_i, j_node_bnd, j_dof, j_size, j_size1, j_size2, j_vertex, i,j,jj, j_node, s_node, j_dir,     &
    !$omp   j_index, j_starwall, j_tor, s_tor, j_resp, j_col_psi, j_col_zj, sparsepos_jp, sparsepos_pp, sparsepos_jj,                        &
    !$omp   amat_contrib, rhs_contrib, amat_contrib_psi, amat_contrib_zj, blockpos_jp, blockpos_pp, ierr )                                   &
    !$omp schedule(dynamic,1) collapse(4)
    L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements

      ! --- Select a test function (the weak form equation must hold for every test function)
      L_LS: do l_tor = a_mat%i_tor_min, a_mat%i_tor_max ! (loop over toroidal harmonics)
        L_LV: do l_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
          L_LD: do l_dof = 1, 2 ! (loop over node dofs)
     if (my_id .eq. 0 .and. m_bndelem.eq.1 .and. l_tor.eq.a_mat%i_tor_min .and. l_vertex.eq.1 .and. l_dof.eq.1) write(*,*) 'Inside loop'

            bndelem_m = bnd_elm_list%bnd_element(m_bndelem)

            ! --- Determine the values of R,s and Z,s at the Gaussian points.
            call det_coord_bnd(bndelem_m, node_list, R_S=x_s, Z_S=y_s)

            l_node      = bndelem_m%vertex(l_vertex)
            l_dir       = bndelem_m%direction(l_vertex,l_dof)
            l_node_bnd  = bndelem_m%bnd_vertex(l_vertex)
            l_index     = node_list%node(l_node)%index(l_dir)
            l_size      = bndelem_m%size(l_vertex,l_dof)

            if ( (l_index < index_min) .or. (l_index > index_max) ) cycle ! This MPI proc responsible?

            ! --- Determine the row in the main matrix.
            l_row_psi = det_row_col(l_index, var_psi, l_tor, a_mat%i_tor_min, a_mat%i_tor_max)
            l_row_j   = det_row_col(l_index, var_zj,  l_tor, a_mat%i_tor_min, a_mat%i_tor_max)

            ! --- Sum over boundary dofs at which response is calculated
            L_IV: do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)

              i_node      = bndelem_m%vertex(i_vertex)
              i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)

              L_ID: do i_dof = 1, 2 ! (loop over node dofs)

                i_dir       = bndelem_m%direction(i_vertex,i_dof)
                i_index     = node_list%node(i_node)%index(i_dir)
                i_size      = bndelem_m%size(i_vertex,i_dof)

                L_IS: do i_starwall = 1, sr%n_tor ! (loop over STARWALL harmonics)

                  i_tor    = sr%i_tor(i_starwall)

                  if ( (i_tor < a_mat%i_tor_min) .or. (i_tor > a_mat%i_tor_max) ) cycle    

                  i_resp_old   = response_index(i_node_bnd,i_starwall,i_dof)

                  i_resp   = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor0 &
                           + bnd_node_list%bnd_node(i_node_bnd)%n_dof*(i_starwall-1)              &
                           + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(i_dof)             &
                           - bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1

!                  if (i_resp_old .ne. i_resp) write(*,'(A,8i5)') 'PANIC! : ',i_node, i_starwall, &
!                  i_dof,bnd_node_list%bnd_node(i_node_bnd)%index_starwall,i_resp_old, i_resp

                  i_resp_0 = response_index_eq(i_node_bnd,i_dof)

                  ! --- Loop over Gaussian points -- integration in s-direction
                  common_prefactor = 0.
                  L_MS: do ms = 1, n_gauss

                    ! --- Integration factor from the definition of dA:
                    !     int dA = sum_{m_bndelem} int ds int dphi  sqrt{(R,s)^2 + (Z,s)^2}
                    dA = sqrt(x_s(ms)**2 + y_s(ms)**2)

                    ! --- Loop over toroidal planes -- integration in phi-direction
                    L_MP: do m_plane = 1, n_plane

                      ! --- Evaluate test function at current position
                      testfunc_l = H1(l_vertex,l_dof,ms) * l_size * HZ(l_tor,m_plane)

                      ! --- Determine basis function
                      basfunc_i = H1(i_vertex,i_dof,ms) * i_size * HZ(i_tor,m_plane)

#ifdef __GFORTRAN__
                      common_prefactor = common_prefactor + wgauss_copy(ms) * dA * testfunc_l * basfunc_i
#else
                      common_prefactor = common_prefactor + wgauss(ms) * dA * testfunc_l * basfunc_i
#endif

                    end do L_MP

                  end do L_MS

                  ! --- Sum over boundary dofs contributing to the response
                  L_JB: do j_node_bnd = 1, bnd_node_list%n_bnd_nodes ! (loop over boundary nodes)
                    
                    ! Index of the node in the node list
                    j_node      = bnd_node_list%bnd_node(j_node_bnd)%index_jorek
    
                    L_JD: do j_dof = 1, 2 ! (loop over node dofs)

                      j_dir       = bnd_node_list%bnd_node(j_node_bnd)%direction(j_dof)
                      j_index     = node_list%node(j_node)%index(j_dir)
                      L_JS: do j_starwall = 1, sr%n_tor ! (loop over STARWALL harmonics)
                      
                        j_tor  = sr%i_tor(j_starwall)
                    
                        if ( (j_tor < a_mat%i_tor_min) .or. (j_tor > a_mat%i_tor_max) ) cycle    

                        j_resp   = (bnd_node_list%bnd_node(j_node_bnd)%index_starwall(1) - 1)*sr%n_tor0 &
                                 + bnd_node_list%bnd_node(j_node_bnd)%n_dof*(j_starwall-1)              &
                                 + bnd_node_list%bnd_node(j_node_bnd)%index_starwall(j_dof)             &
                                 - bnd_node_list%bnd_node(j_node_bnd)%index_starwall(1) + 1                      

                        ! --- Option to switch off mode coupling due to a 3D  wall
                        if ( vacuum_decouple_modes .and. (j_tor /= i_tor) ) cycle

                        ! --- Determine the column in the main matrix
                        j_col_psi = det_row_col(j_index, var_psi, j_tor, a_mat%i_tor_min, a_mat%i_tor_max)
                        j_col_zj  = det_row_col(j_index, var_zj, j_tor, a_mat%i_tor_min, a_mat%i_tor_max)

                        ! --- Determine the position in the sparse matrix data structure
                        !     which corresponds to the matrix entry at  l_row_j, j_col_psi.
                        sparsepos_jp = det_sparse_pos(l_row_j,   j_col_psi, index_min, a_mat)
                        sparsepos_pp = det_sparse_pos(l_row_psi, j_col_psi, index_min, a_mat)
                        sparsepos_jj = det_sparse_pos(l_row_j,   j_col_zj, index_min, a_mat)

                        ! --- Vacuum response contribution to the lhs of the current equation
                        if (CARIDDI_mode .and. halo_freeb .and. .not. halo_imp ) then
                          amat_contrib = - common_prefactor * response_m_e_i(i_resp, j_resp)
                        elseif (CARIDDI_mode .and. halo_freeb .and. halo_imp) then

                        ! if (my_id.eq.0 .and. m_bndelem.eq.1 .and. l_tor.eq.1 .and. l_vertex.eq.1 &
                        !   .and. l_dof.eq.1 .and. i_vertex.eq.1 .and. i_dof.eq.1 .and. i_starwall.eq.1 &
                        !   .and. j_node_bnd.eq.1 .and. j_dof.eq.1 .and. j_starwall.eq.1) write(*,*) 'Performing products for halo implicit'
                
                          amat_contrib_psi = - common_prefactor * response_m_eh(i_resp,j_resp)              
                          amat_contrib_zj  = - common_prefactor * response_m_zh(i_resp,j_resp) 

                      ! if (my_id.eq.0 .and. m_bndelem.eq.1 .and. l_tor.eq.1 .and. l_vertex.eq.1 &
                      !   .and. l_dof.eq.1 .and. i_vertex.eq.1 .and. i_dof.eq.1 .and. i_starwall.eq.1 &
                      !   .and. j_node_bnd.eq.1 .and. j_dof.eq.1 .and. j_starwall.eq.1) write(*,*) 'Done Performing products for halo implicit'

                        else  ! Eddy coupling
                          amat_contrib = - common_prefactor * response_m_e(i_resp, j_resp)
                        endif

                        if (.not. halo_imp) then
                          !$omp atomic
                          a_mat%val(sparsepos_jp) = a_mat%val(sparsepos_jp) + amat_contrib
                        else
                          !$omp atomic
                          a_mat%val(sparsepos_jp) = a_mat%val(sparsepos_jp) + amat_contrib_psi
                          !$omp atomic
                          a_mat%val(sparsepos_jj) = a_mat%val(sparsepos_jj) + amat_contrib_zj
                        endif

                      end do L_JS
                    end do L_JD
                  end do L_JB

                  ! --- Contribution of vacuum response to the rhs of the current equation
                  if (CARIDDI_mode .and. halo_freeb .and. .not. halo_imp) then 
                    rhs_contrib = sum( response_m_h(i_resp,:) * psibnd_vec(:))                     &
                              + sum( response_m_j_i(i_resp,:) * dpsibnd_vec(:))
                  
                  elseif (CARIDDI_mode .and. halo_freeb .and. halo_imp) then

                    rhs_contrib = sum( response_m_h(i_resp,:) * psibnd_vec(:))                     &
                              + sum( response_m_j_i(i_resp,:) * dpsibnd_vec(:))                    &
                              + sum(response_m_ph(i_resp,:)*Ielec_w(:))                           &
                              + sum(response_m_qh(i_resp,:)*dIelec_w(:))                           &
                              + sum(response_m_xh(i_resp,:)*Ielec_J(:))                           &
                              + sum(response_m_vh(i_resp,:)*dIelec_J(:))         
                                                    
                  else  ! Eddy coupling
                  
                    rhs_contrib = sum( response_m_h(i_resp,:) * psibnd_vec(:))                     &
                              + sum( response_m_j(i_resp,:) * dpsibnd_vec(:))
                  endif

                  if ( (l_tor == 1) .and. (sr%i_tor(1) == 1) .and. (.not. starwall_equil_coils)) &
                    rhs_contrib = rhs_contrib - sum( bext_tan(i_resp_0, :) * I_coils(:) )        &
                                - sum( response_m_h(i_resp,:) *   psibnd_coils(:) )

                  if ( resistive_wall )rhs_contrib=rhs_contrib+rhs_contrib_arr(i_resp)

                  rhs_contrib = rhs_contrib * common_prefactor
                  !$omp atomic
                  rhs_loc(l_row_j) = rhs_loc(l_row_j) + rhs_contrib

                end do L_IS
              end do L_ID
            end do L_IV

          end do L_LD
        end do L_LV
      end do L_LS
    end do L_MB
    !$omp end parallel do
        if (my_id .eq. 0) write(*,*) 'End of vacuum_boundary_integral do loops'

    ! --- timing
    if ( vacuum_debug ) then
      t_elaps_end = MPI_WTIME()
      write(*,'(I4,A,F10.7,A)') my_id, '  Elapsed time vacuum_boundary_integral', t_elaps_end - t_elaps_start, ' s'
    end if
    
    if ( vacuum_debug ) write(*,*) my_id, 'After:', sum(abs(rhs_loc)),sum(abs(a_mat%val))

    if ( allocated(psibnd_vec ) ) deallocate( psibnd_vec  )
    if ( allocated(dpsibnd_vec) ) deallocate( dpsibnd_vec )
    if (allocated(dzjbnd_vec))    deallocate( dzjbnd_vec  )

  end subroutine vacuum_boundary_integral
  
  

  !> Determine the values of \f$R\f$, \f$Z\f$, \f$\partial R/\partial s\f$, and
  !! \f$ \partial Z/\partial s \f$ on the Gaussian points of a given boundary element.
  subroutine det_coord_bnd(bndelem, node_list, R, Z, R_s, Z_s)
    
    use gauss,             only: n_gauss, xgauss, wgauss
    use data_structure,    only: type_node, type_bnd_element, type_node_list
    use basis_at_gaussian, only: H1, H1_s, HZ
    
    implicit none
    
    ! --- Routine parameters
    type(type_bnd_element), intent(in)  :: bndelem       !< Boundary element to be considered.
    type(type_node_list),   intent(in)  :: node_list     !< List of grid nodes
    real*8, optional,       intent(out) :: R(n_gauss)    !< Values of R on Gaussian points
    real*8, optional,       intent(out) :: Z(n_gauss)    !< Values of Z on Gaussian points
    real*8, optional,       intent(out) :: R_s(n_gauss)  !< Values of R,s on Gaussian points
    real*8, optional,       intent(out) :: Z_s(n_gauss)  !< Values of Z,s on Gaussian points
    
    ! --- Local variables
    integer         :: k_vertex, k_dof, k_node, k_dir
    real*8          :: k_size
    type(type_node) :: node_k
    
    if ( present(R  ) ) R   = 0.d0
    if ( present(Z  ) ) Z   = 0.d0
    if ( present(R_s) ) R_s = 0.d0
    if ( present(Z_s) ) Z_s = 0.d0
    
    do k_vertex = 1, 2
      do k_dof = 1, 2
        k_node      = bndelem%vertex(k_vertex)
        k_dir       = bndelem%direction(k_vertex,k_dof)
        k_size      = bndelem%size(k_vertex,k_dof)
        node_k      = node_list%node(k_node)
        if ( present(R  ) ) R  (:)  = R  (:)  + node_k%x(1,k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
        if ( present(Z  ) ) Z  (:)  = Z  (:)  + node_k%x(1,k_dir,2) * k_size * H1  (k_vertex,k_dof,:)
        if ( present(R_s) ) R_s(:)  = R_s(:)  + node_k%x(1,k_dir,1) * k_size * H1_s(k_vertex,k_dof,:)
        if ( present(Z_s) ) Z_s(:)  = Z_s(:)  + node_k%x(1,k_dir,2) * k_size * H1_s(k_vertex,k_dof,:)
      end do
    end do
    
  end subroutine det_coord_bnd
  
  
  
  !> Determine vectors of the psi and deltapsi values at the boundary.
  subroutine det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec, dzjbnd_vec, psibnd_coils)

    use data_structure, only: type_node_list, type_bnd_node_list

    implicit none

    ! --- Routine parameters
    type(type_node_list),     intent(in)  :: node_list      !< List of grid nodes
    type(type_bnd_node_list), intent(in)  :: bnd_node_list  !< List of boundary grid nodes
    real*8, allocatable,      intent(out) :: psibnd_vec(:)  !< Vector of Psi boundary values
    real*8, ALLOCATABLE, optional, intent(out) :: dzjbnd_vec(:)
    real*8, allocatable,      intent(out) :: dpsibnd_vec(:) !< Vector of deltaPsi boundary values
    real*8, allocatable, optional, intent(out) :: psibnd_coils(:)!< Vector of Psi_coil boundary values

    ! --- Local variables
    integer :: jnode, jnode_glob, j_starwall, jtor, jbas, jdir, j_resp, j_resp_0, j_resp_old

    if ( allocated(psibnd_vec) ) deallocate(psibnd_vec)
    allocate( psibnd_vec(n_dof_starwall) )
    psibnd_vec(:) = 0.d0

    if ( allocated(dpsibnd_vec) ) deallocate(dpsibnd_vec)
    allocate( dpsibnd_vec(n_dof_starwall) )
    dpsibnd_vec(:) = 0.d0

    if ( present(psibnd_coils) ) then
      if ( allocated(psibnd_coils) ) deallocate(psibnd_coils)
      allocate( psibnd_coils(n_dof_starwall) )
      psibnd_coils(:) = 0.d0
    end if

    if (present(dzjbnd_vec)) then 
      if (ALLOCATED(dzjbnd_vec)) deallocate(dzjbnd_vec)
      ALLOCATE( dzjbnd_vec(n_dof_starwall))
      dzjbnd_vec(:) = 0.d0 
    endif

    ! --- Determine vector of (delta)psi boundary values.
    do jnode = 1, bnd_node_list%n_bnd_nodes       ! loop over nodes

      jnode_glob = bnd_node_list%bnd_node(jnode)%index_jorek

      do j_starwall = 1, sr%n_tor     ! loop over STARWALL harmonics

        jtor = sr%i_tor(j_starwall)   ! (mode corresponding to STARWALL harmonic)

        do jbas = 1, 2                ! loop over basis functions

          jdir         = bnd_node_list%bnd_node(jnode)%direction(jbas)
          j_resp_old   = response_index(jnode,j_starwall,jbas)
          j_resp_0     = response_index_eq(jnode,jbas)

          j_resp = (bnd_node_list%bnd_node(jnode)%index_starwall(1) - 1)*sr%n_tor0 &
                 +  bnd_node_list%bnd_node(jnode)%n_dof*(j_starwall-1) &
                 + (bnd_node_list%bnd_node(jnode)%index_starwall(jbas)-bnd_node_list%bnd_node(jnode)%index_starwall(1)) + 1

!          if (j_resp_old .ne. j_resp) write(*,'(A4i5)') 'PANIC jresp: ',j_resp_old,j_resp, &
!           bnd_node_list%bnd_node(jnode)%index_starwall(1), bnd_node_list%bnd_node(jnode)%index_starwall(jbas)

          psibnd_vec ( j_resp ) = node_list%node(jnode_glob)%values(jtor, jdir, var_psi)
          dpsibnd_vec( j_resp ) = node_list%node(jnode_glob)%deltas(jtor, jdir, var_psi)
          if (PRESENT(dzjbnd_vec)) then 
            dzjbnd_vec(j_resp) = node_list%node(jnode_glob)%deltas(jtor, jdir, var_zj)
          endif
          if ( (present(psibnd_coils)) .and. (allocated(I_coils)) .and. (jtor==1) .and. (.not. starwall_equil_coils) ) then
            j_resp_0 = 2*(jnode-1) + jbas
            psibnd_coils( j_resp ) = sum( bext_psi(j_resp_0,:) * I_coils(:) )
          end if

        end do
      end do
    end do
    
  end subroutine det_psibnd_vec
  
  
  
  !> Initialize the currents in the resistive wall and the external coil currents.
  subroutine init_wall_currents(my_id, resistive_wall)
    
    !use nodes_elements, only: bnd_node_list, node_list
    use phys_module, only: index_now, index_start, nstep
    use mpi_mod
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    real*8, allocatable :: psibnd_vec(:)    ! Vector of Psi values at the boundary
    real*8, allocatable :: dpsibnd_vec(:)   ! Vector of deltaPsi values at the boundary
    real*8, allocatable :: wall_and_coil_curr(:)
    real*8, allocatable :: tmp_coil_curr(:)
    integer             :: k, k2, ierr
    
    !call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec)
    
    if ( resistive_wall ) then
      

      if (halo_freeb .and. CARIDDI_mode)  then 
        ! ====================== EDDY+HALO COUPLING ================================
        if ( (allocated(wall_curr_i)) .and. (.not. wall_curr_initialized) ) &
          deallocate(wall_curr_i)
        if ( (allocated(dwall_curr_i)) .and. (.not. wall_curr_initialized) ) &
          deallocate(dwall_curr_i)
        if ( .not. allocated(wall_curr_i) ) then
          allocate( wall_curr_i(n_wall_curr_i) )
          do k = 1, n_wall_curr_i
            wall_curr_i(k) = 0.d0! 
          end do
        endif
        if ( .not. allocated(dwall_curr_i) ) then
          write(*,*) 'Allocating dwall_curr_i with', n_wall_curr_i
          allocate( dwall_curr_i(n_wall_curr_i) )
          do k = 1, n_wall_curr_i
            dwall_curr_i(k) = 0.d0!- sum( sr%a_yie(k,:) * dpsibnd_vec(:) )
          end do
        end if
        if ( (allocated(wall_curr_e)) .and. (.not. wall_curr_initialized) ) &
          deallocate(wall_curr_e)
        if ( (allocated(dwall_curr_e)) .and. (.not. wall_curr_initialized) ) &
          deallocate(dwall_curr_e)

        if ( .not. allocated(wall_curr_e) ) then
          allocate( wall_curr_e(n_wall_curr_e) )
          do k = 1, n_wall_curr_e
            wall_curr_e(k) = 0.d0! 
          end do
        endif

        if ( .not. allocated(dwall_curr_e) ) then
          write(*,*) 'Allocating dwall_curr_e with', n_wall_curr_e
          allocate( dwall_curr_e(n_wall_curr_e) )
          do k = 1, n_wall_curr_e
            dwall_curr_e(k) = 0.d0
          end do
        end if

        if ( .not. allocated(Ielec_w) ) then
          write(*,*) 'Allocating Ielec_w with', n_wall_curr_e
          allocate( Ielec_w(n_wall_curr_e) )
          do k = 1, n_wall_curr_e
            Ielec_w(k) = 0.d0
          end do
        end if       

        if ( .not. allocated(Ielec_J) ) then
          write(*,*) 'Allocating Ielec_J with', n_shell_curr_e
          allocate( Ielec_J(n_shell_curr_e) )
          do k = 1, n_shell_curr_e
            Ielec_J(k) = 0.d0
          end do
        end if 

        if ( .not. allocated(dIelec_w) ) then
          write(*,*) 'Allocating dIelec_w with', n_wall_curr_e
          allocate( dIelec_w(n_wall_curr_e) )
          do k = 1, n_wall_curr_e
            dIelec_w(k) = 0.d0
          end do
        end if       

        if ( .not. allocated(dIelec_J) ) then
          write(*,*) 'Allocating dIelec_J with', n_shell_curr_e
          allocate( dIelec_J(n_shell_curr_e) )
          do k = 1, n_shell_curr_e
            dIelec_J(k) = 0.d0
          end do
        end if        

      else
        ! ======================= EDDY COUPLING ===================================
        if ( (allocated(wall_curr)) .and. (.not. wall_curr_initialized) ) &
          deallocate(wall_curr)
        if ( (allocated(dwall_curr)) .and. (.not. wall_curr_initialized) ) &
          deallocate(dwall_curr)

        if ( .not. allocated(wall_curr) ) then
          allocate( wall_curr(n_wall_curr) )
          do k = 1, n_wall_curr
            wall_curr(k) = 0.d0! - sum( sr%a_ye(k,:) * psibnd_vec(:) )
          end do
        end if
      
        if ( .not. allocated(dwall_curr) ) then
          write(*,*) 'Allocating dwall_curr with', n_wall_curr
          allocate( dwall_curr(n_wall_curr) )
          do k = 1, n_wall_curr
            dwall_curr(:) = 0.d0!- sum( sr%a_ye(k,:) * dpsibnd_vec(:) )
          end do
        end if        
      

      endif
      
      ! ==========================================================================
       
    end if
    
    ! --- Also initialize the old_dpsibnd_vec
    if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
    allocate( old_dpsibnd_vec(n_dof_starwall) )
    old_dpsibnd_vec(:) = 0.d0

    if (halo_freeb) then 
      if (allocated(old_dIelec_w)) deallocate(old_dIelec_w)
      allocate(old_dIelec_w(n_wall_curr_e))
      old_dIelec_w(:) = 0.d0
      if (allocated(old_dIelec_J)) deallocate(old_dIelec_J)
      allocate(old_dIelec_J(n_shell_curr_e))
      old_dIelec_J(:) = 0.d0
    endif

    if (.not. CARIDDI_mode .and. .not. vacuum_min) call write_wall_vtk(0, resistive_wall, my_id)
    !deallocate( psibnd_vec, dpsibnd_vec )
    
    ! wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
    ! ====================================================================================================
    ! ============ HALO NOT IMPLEMENTED HERE SINCE THIS IS NOT FOR CARIDDI !!!!!! =====================  
    ! ====================================================================================================
    ! wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

    ! --- Initialize net toroidal wall current (for plot_live_data)
    if ( resistive_wall .and. (index_now>0) .and. (index_start+nstep>0)) then
      if ( .not. allocated(net_tor_wall_curr) ) then
        allocate( net_tor_wall_curr(index_start+nstep) )
        net_tor_wall_curr(:) = 0.d0
      end if
      if (.not. CARIDDI_mode .and. .not. vacuum_min) then
        k2 = sr%ncoil + 1 
        net_tor_wall_curr(index_now) = sum(sr%s_ww%loc_mat(k2,:) * wall_curr(sr%s_ww%ind_start:sr%s_ww%ind_end))
        call MPI_ALLReduce(MPI_IN_PLACE, net_tor_wall_curr(index_now),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if
    end if
    
    ! --- Initialize coil currents (for plot_live_data)
    if ( resistive_wall .and. (sr%ncoil > 0) .and. (index_now>0) .and. (index_start+nstep>0)) then
      ! Calculate coil currents for all coil types
      if  (.not. allocated(tmp_coil_curr) ) allocate(tmp_coil_curr(sr%ncoil))
      tmp_coil_curr = 0.d0
      if (.not. CARIDDI_mode) then
        do k = 1, sr%ncoil
          tmp_coil_curr(k) = sum(sr%s_ww%loc_mat(k,:) * wall_curr(sr%s_ww%ind_start:sr%s_ww%ind_end))
        end do
        call MPI_ALLReduce(MPI_IN_PLACE, tmp_coil_curr,sr%ncoil,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if

      ! Distribute the coil currents to the different coil types
      if ( sr%n_diag_coils > 0 ) then
        if ( (.not. allocated(diag_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( diag_coil_curr(index_start+nstep,sr%n_diag_coils) )
          diag_coil_curr(:,:) = 0.d0
        end if
        diag_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_diag_coils:sr%ind_start_diag_coils + sr%n_diag_coils -1)
      endif
      
      if ( sr%n_rmp_coils > 0 ) then
        if ( (.not. allocated(rmp_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( rmp_coil_curr(index_start+nstep,sr%n_rmp_coils) )
          rmp_coil_curr(:,:) = 0.d0
        end if
        rmp_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_rmp_coils : sr%ind_start_rmp_coils + sr%n_rmp_coils -1) 
      end if
      
      if ( sr%n_pol_coils > 0 ) then
        if ( (.not. allocated(pf_coil_curr)) .and. (index_start+nstep >0) ) then
           allocate( pf_coil_curr(index_start+nstep,sr%n_pol_coils) )
           pf_coil_curr(:,:) = 0.d0
        end if
        pf_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_pol_coils:sr%ind_start_pol_coils + sr%n_pol_coils -1) 
      end if

      ! Voltage coils not yet implemented
      ! if ( sr%n_voltage_coils > 0 ) then
      !  if ( (.not. allocated(voltage_coil_curr)) .and. (index_start+nstep >0) ) then
      !    allocate( voltage_coil_curr(index_start+nstep,sr%n_voltage_coils) )
      !    voltage_coil_curr(:,:) = 0.d0
      !  end if
      !  voltage_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_voltage_coils:sr%ind_start_voltage_coils + sr%n_voltage_coils -1) 
      !  endif
      !end if
    endif
          
 
    if ( vacuum_debug .and. (my_id == 0) ) write(*,*) 'Wall currents initialized.'
    wall_curr_initialized = .true.
    
  end subroutine init_wall_currents
  
  
  
  !> Imposing PF coil currents as a current source term
  subroutine coil_current_source(my_id)
    
    use constants
    use phys_module, only: t_now, index_now, tstep, nstep, index_start
    use equil_info, only: ES
    use profiles
    use mpi_mod    

    implicit none
   
    ! --- Routine parameters
    integer, intent(in)       :: my_id  
    integer                   :: i, i1, i2, ierr
    real*8, allocatable       :: potentials_real_0(:)
    real*8, save, allocatable :: delta_Icoils_0(:)
    real*8, save              :: t_last=-1e3, Z_p = 1.d10
    logical, save             :: initialized=.false.
    real*8                    :: z_ref_inter, Z_n
    if( my_id == 0 ) write(*,*) ' Imposing PF coil currents with a current source term '

    ! --- Calculate the difference between the input file coil currents and the ones coming from the restart file
    ! This is used below to avoid violent jumps in I_coils that can happen when restarting from a freeboudary_equil with feedback
    if(.not. initialized) then
      if (allocated(delta_Icoils_0)) deallocate(delta_Icoils_0)
      allocate(delta_Icoils_0(n_coils))
      delta_Icoils_0 = 0.d0
      do i=1, n_coils
          delta_Icoils_0(i) = I_coils(i) - interpolProf(coil_curr_time_trace(i)%time, coil_curr_time_trace(i)%curr, coil_curr_time_trace(i)%len, t_now)
      end do
      initialized = .true.
    end if
    if (t_last==-1e3) t_last=t_now

    ! --- Calculate the specified coil currents at present time 
    do i=1, n_coils
      I_coils(i) = interpolProf(coil_curr_time_trace(i)%time, coil_curr_time_trace(i)%curr, coil_curr_time_trace(i)%len, t_now) 
    end do

    ! -- During timestepping apply vertical feedback by the PF coils which were activated
    ! Add the feedback current to the difference from the restart file
    if (.not. allocated(vert_FB_response) .and. index_start+nstep>0) then
      allocate(vert_FB_response(index_start+nstep,4))
      vert_FB_response = 0.d0
    endif
    z_ref_inter = interpolProf(Z_axis_ref_ts%time, Z_axis_ref_ts%position ,  Z_axis_ref_ts%len, t_now)
    vert_FB_response(index_now,4) = z_ref_inter
    if ( t_now>start_VFB_ts .and. sum( abs(vert_FB_amp_ts(1:n_pf_coils)) )>1.e-6 .and. sr%i_tor(1) ==  1 ) then
      Z_n = ES%Z_axis
      if ( Z_p > 1.d9 ) Z_p = Z_n
      dZ_axis_integral = dZ_axis_integral + ( Z_n - z_ref_inter )*tstep
      if ( (t_now-t_last)>vert_FB_tact ) then 
        if (my_id==0) write(*,*) 'Vertical feedback active'
        t_last = t_now
        vert_FB_response(index_now,1:3) =  &
            (/  1.d3  * vert_FB_gain(1)* (Z_n - z_ref_inter), &
                1.d5  * vert_FB_gain(2)* (Z_n - Z_p )/ tstep, &
                1.d-3 * vert_FB_gain(3)* dZ_axis_integral /)
        delta_Icoils_0(1:n_pf_coils) = delta_Icoils_0(1:n_pf_coils) +  sum(vert_FB_response(index_now,1:3))*vert_FB_amp_ts(1:n_pf_coils)
      endif
      Z_p = Z_n
    endif

    

    ! -- Apply coil current limits
    do i=1, n_coils
      if ( abs( I_coils(i) + delta_Icoils_0(i)) .gt. I_coils_max(i) ) &
          delta_Icoils_0(i) = ( I_coils(i) + delta_Icoils_0(i) ) / abs( I_coils(i) + delta_Icoils_0(i))&
                                  *I_coils_max(i) - ( I_coils(i) + delta_Icoils_0(i) )
      I_coils(i) = I_coils(i)  + delta_Icoils_0(i)
    end do


    if (.not. allocated (potentials_real_0)) allocate(potentials_real_0(n_coils))
    ! --- Transform real currents into starwall currents
    potentials_real_0 = 0.d0
    potentials_real_0(:) = I_coils(:)  * mu_zero
    
    ! ================================ EDDY + HALO COUPLING ====================================================
    if (halo_freeb .and. CARIDDI_mode ) then
      if (.not. allocated (Yi_coils0)) allocate(Yi_coils0(n_wall_curr_i))
      Yi_coils0          = 0.d0
      i1 = sr%ind_start_coils;       i2 = sr%ind_start_coils + sr%ncoil -1
      do i = 1, n_wall_curr_i
        if ( (i>=sr%s_ww_ii_inv%ind_start) .and. (i<=sr%s_ww_ii_inv%ind_end) ) then
          Yi_coils0(i) = sum(sr%s_ww_ii_inv%loc_mat(i-my_id*sr%s_ww_ii_inv%step,i1:i2)*potentials_real_0(:))
        end if
      end do
      call MPI_ALLReduce(MPI_IN_PLACE, Yi_coils0, size(Yi_coils0),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    else
      ! ----------------------- EDDY COUPLING ---------------------------------------------------------------
      if (.not. allocated (Y_coils0)) allocate(Y_coils0(n_wall_curr))
      Y_coils0          = 0.d0
      i1 = sr%ind_start_coils;       i2 = sr%ind_start_coils + sr%ncoil -1
      do i = 1, n_wall_curr
        if ( (i>=sr%s_ww_inv%ind_start) .and. (i<=sr%s_ww_inv%ind_end) ) then
          Y_coils0(i) = sum(sr%s_ww_inv%loc_mat(i-my_id*sr%s_ww_inv%step,i1:i2)*potentials_real_0(:))
        end if
      end do
      call MPI_ALLReduce(MPI_IN_PLACE, Y_coils0, size(Y_coils0),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    endif
    ! ===========================================================================================================
  end subroutine coil_current_source
  
  
  
  !> Perform the time-evolution of the wall currents (resistive wall).
  subroutine evolve_wall_currents(my_id, psibnd_vec, dpsibnd_vec, dzjbnd_vec,  a_mat)
    
    use phys_module, only: index_now, index_start, nstep, resistive_wall, time_evol_theta, time_evol_zeta, tstep, nstep_n
    use mpi_mod
    use nodes_elements
    use data_structure, only: type_SP_MATRIX

    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    real*8,  intent(in) :: psibnd_vec (n_dof_starwall) !< Vector of Psi boundary values
    real*8,  intent(in) :: dpsibnd_vec(n_dof_starwall) !< Vector of deltaPsi boundary values
    real*8, optional, intent(in) :: dzjbnd_vec(n_dof_starwall)
    type(type_SP_MATRIX):: a_mat
    
    ! --- Local variables
    integer              :: k, k2, i, j
    integer              :: ierr, n_steps
    real *8, allocatable :: dwall_curr_2(:)
    real *8, allocatable :: dwall_curr_i2(:)   
    real *8, allocatable :: dwall_curr_e2(:)      
    real *8, allocatable :: tmp_coil_curr(:)
    real *8              :: zeta, theta
    real *8, allocatable :: tmp_d_s(:)
    real *8, allocatable :: inv_tmp_d_s(:)
    real *8, allocatable :: tmp_d_sii(:)
    real *8, allocatable :: inv_tmp_d_sii(:)
    real*8,  ALLOCATABLE :: diag_1(:)
    real*8,  ALLOCATABLE :: resp_R_zj, resp_R_psis
    real*8, ALLOCATABLE  :: response_R_psis(:,:), response_R_zj(:,:), wall_curr_i_dbg(:) 


    if (halo_freeb .and. CARIDDI_mode) then

      ! ===================================================================================================
      ! ================================ EDDY + HALO COUPLING =============================================
      ! ===================================================================================================

      if (my_id .eq. 0) then 

        OPEN(UNIT=12, FILE="dpsibnd_vec.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%nd_bez
          WRITE(12,*) dpsibnd_vec(k)
        END DO

        OPEN(UNIT=12, FILE="dzjbnd_vec.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%nd_bez
          WRITE(12,*) dzjbnd_vec(k)
        END DO
      endif

      ! Calculate deltaIel_w from deltazj and deltapsi and Iel_w at present
      if (.not. ALLOCATED(response_R_psis)) allocate(response_R_psis(sr%nd_bez, sr%nd_bez))
      if (.not. ALLOCATED(response_R_zj)) allocate(response_R_zj(sr%nd_bez, sr%nd_bez))
      
      ! --- EVOLVE ELECTRODE CURRENTS
      call halo_time_stepping(my_id, node_list, element_list, bnd_elm_list, bnd_node_list, &
                              response_R_zj=response_R_zj, response_R_psis=response_R_psis, debug=.false., &
                              dzjbnd_vec=dzjbnd_vec, dpsibnd_vec=dpsibnd_vec) 
      
      ! Starts from 1
      if (my_id .eq. 0) write(*,*) 'index_now=' , index_now

      ! --- Fill input matrices at the corresponding column

      ! Update the offline matrices
      n_steps = sum(nstep_n)
      if (.not. allocated (psibnd_mat)) then 
          write(*,*) 'Allocating psibnd_mat with', n_dof_starwall , n_steps
          allocate( psibnd_mat(n_dof_starwall , n_steps) )
          psibnd_mat(:,:) = 0.d0
      end if    
      if (.not. allocated (Ielecw_mat)) then 
          write(*,*) 'Allocating Ielecw_mat with', n_wall_curr_e , n_steps
          allocate( Ielecw_mat(n_wall_curr_e , n_steps) )
          Ielecw_mat(:,:) = 0.d0
      end if  
      if (.not. allocated (Ielecj_mat)) then 
          write(*,*) 'Allocating IelecJ_mat with', n_shell_curr_e , n_steps
          allocate( Ielecj_mat(n_shell_curr_e , n_steps) )
          Ielecj_mat(:,:) = 0.d0
      end if  
      if (.not. allocated (wall_curr_i_mat)) then 
          write(*,*) 'Allocating wall_curr_i_mat with', n_wall_curr_i , n_steps
          allocate( wall_curr_i_mat(n_wall_curr_i , n_steps) )
          wall_curr_i_mat(:,:) = 0.d0
      end if   
     
      if (.not. allocated (coil_curr_i_mat)) then
          write(*,*) 'Allocating coil_curr_i_mat with', n_wall_curr_i , n_steps
          allocate( coil_curr_i_mat(n_wall_curr_i , n_steps) )
          coil_curr_i_mat(:,:) = 0.d0
      end if

      psibnd_mat(:,index_now) = psibnd_vec(:)
      Ielecw_mat(:,index_now) = Ielec_w(:)
      Ielecj_mat(:,index_now) = Ielec_J(:)

      
      ! Print the offline matrices at the end of the sim
      !if (index_now .eq. n_steps ) then 

        if (my_id .eq. 0) then 

          open(unit=10, file="psibnd_mat.txt", status="replace")
          do i = 1, n_dof_starwall
            write(10,'(*(ES24.16,1X))') (psibnd_mat(i,j), j = 1, n_steps)
          end do
          close(10)

          open(unit=10, file="Ielecw_mat.txt", status="replace")
          do i = 1, n_wall_curr_e
            write(10,'(*(ES24.16,1X))') (Ielecw_mat(i,j), j = 1, n_steps)
          end do
          close(10)

          open(unit=10, file="IelecJ_mat.txt", status="replace")
          do i = 1, n_shell_curr_e
            write(10,'(*(ES24.16,1X))') (IelecJ_mat(i,j), j = 1, n_steps)
          end do
          close(10)

        endif
      !endif

      
      ! ---- Debug deltas
      if (my_id .eq. 0) then 

        if ( index_now .eq. 2 ) then

          OPEN(UNIT=12, FILE="I_elec_w_2.dat", ACTION="write", STATUS="replace")
          DO k=1,sr%n_electrode_w-1
            WRITE(12,*) Ielec_w(k)
          END DO

         OPEN(UNIT=12, FILE="delta_I_elec_w_2.dat", ACTION="write", STATUS="replace")
          DO k=1,sr%n_electrode_w-1
            WRITE(12,*) dIelec_w(k)
          END DO 

          OPEN(UNIT=12, FILE="I_elec_J_2.dat", ACTION="write", STATUS="replace")
          DO k=1,sr%n_electrode_J-1
            WRITE(12,*) Ielec_J(k)
          END DO

          OPEN(UNIT=12, FILE="delta_I_elec_J_2.dat", ACTION="write", STATUS="replace")
          DO k=1,sr%n_electrode_J-1
            WRITE(12,*) dIelec_J(k)
          END DO  

        endif 

        if ( index_now .eq. 3 ) then 

          OPEN(UNIT=12, FILE="I_elec_w_3.dat", ACTION="write", STATUS="replace")
          DO k=1,sr%n_electrode_w-1
            WRITE(12,*) Ielec_w(k)
          END DO

          OPEN(UNIT=12, FILE="I_elec_J_3.dat", ACTION="write", STATUS="replace")
          DO k=1,sr%n_electrode_J-1
            WRITE(12,*) Ielec_J(k)
          END DO      

        endif
      endif  

      if (my_id .eq. 0) then 
        OPEN(UNIT=12, FILE="I_elec_w_from_time_stepping.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) Ielec_w(k)
        END DO

        OPEN(UNIT=12, FILE="delta_I_elec_w_from_time_stepping.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) dIelec_w(k)
        END DO
      endif


      ! --------------------------------- INTERNAL DOFS ------------------------------------------------

      if  (.not. allocated(dwall_curr_i2)) allocate (dwall_curr_i2(n_wall_curr_i))
      if  (.not. allocated(dwall_curr_e2)) allocate (dwall_curr_e2(n_wall_curr_e))
      if ( vacuum_debug ) write(*,*) 'Halo: wall_curr_i and wall_curr_e(before)', sum(abs(wall_curr_i)), sum(wall_curr_i) ,sum(abs(wall_curr_e)), sum(wall_curr_e)  
      if ( (.not. allocated(old_dpsibnd_vec)) .or. (size(old_dpsibnd_vec,1)/= n_dof_starwall) ) then
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
      end if 
      if ( (.not. allocated(old_dIelec_J)) .or. (size(old_dIelec_J,1)/= n_shell_curr_e) ) then
        if ( allocated(old_dIelec_J) ) deallocate(old_dIelec_J)
        allocate( old_dIelec_J(n_shell_curr_e) )
      end if   
      if ( (.not. allocated(old_dIelec_w)) .or. (size(old_dIelec_w,1) /= n_wall_curr_e) ) then
        if ( allocated(old_dIelec_w) ) deallocate(old_dIelec_w)
        allocate( old_dIelec_w(n_wall_curr_e) )
      end if      

      !> WARNING: The following computation of dwall_curr_2 contains the term
      !!          dwall_curr which refers to the previous time-step
      if (my_id .eq. 0) write(*,*) 'n_wall_curr_i', n_wall_curr_i
      do k = 1, n_wall_curr_i
      
      ! --- contribution from not distributed matrices
      dwall_curr_i2(k) = response_d_b_i(k) * (wall_curr_i(k)  - Yi_coils0(k)) &  
          + response_d_c_i(k) * dwall_curr_i(k)

      end do

      ! Save coil current (transformed space) in  matrix

      coil_curr_i_mat(:,index_now) = Yi_coils0(:)

      if (my_id .eq. 0) then
        open(unit=10, file="coil_curr_i_mat.txt", status="replace")
        do i = 1, n_wall_curr_i
          write(10,'(*(ES24.16,1X))') (coil_curr_i_mat(i,j), j = 1, n_steps)
        end do
        close(10)
      endif

      theta = time_evol_theta
      zeta  = time_evol_zeta

      if  (.not. allocated(tmp_d_sii)     ) allocate (tmp_d_sii     (n_wall_curr_i))
      if  (.not. allocated(inv_tmp_d_sii) ) allocate (inv_tmp_d_sii (n_wall_curr_i))
      tmp_d_sii(:) = 1.d0 + zeta + tstep*theta*wall_resistivity*sr%dii_yy(:)

      inv_tmp_d_sii(:) = 1.0d0 / tmp_d_sii(:)

      ! !> WARNING: Now we can compute the dwall_curr term for the current time-step
      ! !!          to be added to the contribution from the previous time-step

      ! dwall_curr_i(:) = 0.d0

      ! ! N.B. Sign change is present and accounted for in read_starwall_response!
      ! ! Ci_star from CARIDDI already contains a minus sign! See documentation
      ! !> corresponds to matrix A_i in documentation
      ! if (my_id.eq.0) write(*,*) 'Calling for sr%a_yie*inv_tmp_d_sii*dpsibnd_vec'
      ! !if (my_id.eq.0) write(*,*) 'Size of dpsibnd_vec and old_dpsibnd_vec', size(dpsibnd_vec) , size(old_dpsibnd_vec)
      ! call extended_matrix_vector(my_id, -(1.d0+zeta), sr%a_yie, inv_tmp_d_sii, &
      !                           dpsibnd_vec(:), .false., dwall_curr_i, .false.) 
      ! if (my_id.eq.0) write(*,*) 'Calling for sr%a_yie*inv_tmp_d_sii*old_dpsibnd_vec'                                
      ! !> corresponds to matrix D_i in documentation
      ! call extended_matrix_vector(my_id, zeta, sr%a_yie, inv_tmp_d_sii, &
      !                           old_dpsibnd_vec(:), .false., dwall_curr_i, .true.) 
      ! !if (my_id.eq.0) write(*,*) 'DONE Calling for sr%a_yie*inv_tmp_d_sii*old_dpsibnd_vec' 

      ! if (my_id .eq. 0) then
      !   write(22,*) 'M_yie: infos', sr%a_yie%distrib, sr%a_yie%row_wise, sr%a_yie%ind_start, sr%a_yie%ind_end, sr%a_yie%step
      !   22 format(1x, 'distributed', L2, ',rowwise?', L2, ',ind_start of chunk:', I4, ',ind_end of chunk:', I4, ', step' , I4)
      ! endif
       
      ! dwall_curr_i(:) = dwall_curr_i(:) + dwall_curr_i2(:)
      ! if (index_now <= 1 ) dwall_curr_i = 0.d0
      ! wall_curr_i(:) = wall_curr_i(:) + dwall_curr_i(:) 
      ! if (my_id == 0) then
      !   OPEN(UNIT=12, FILE="wall_curr_i.dat", ACTION="write", STATUS="replace")
      !   DO k=1,n_wall_curr_i
      !       WRITE(12,*) wall_curr_i(k)
      !   END DO
      !   CLOSE(12)  
      ! endif

      dwall_curr_i = 0.d0
      !!! CALCULATE WALL_CURR_i WITH ALL CORRECTIONS
      ! N.B. Sign change is present and accounted for in read_starwall_response!
      ! Ci_star from CARIDDI already contains a minus sign! See documentation
      if (my_id.eq.0) write(*,*) 'CALCULATING COMPLETE WALL_CURR_I'

      !> corresponds to matrix A_i in documentation
      if (my_id.eq.0) write(*,*) 'Calling for sr%a_yie*inv_tmp_d_sii*dpsibnd_vec'
      call extended_matrix_vector(my_id, -(1.d0+zeta), sr%a_yie, inv_tmp_d_sii, &
                                dpsibnd_vec(:), .false., dwall_curr_i, .false.) 
      !> corresponds to matrix JE_i in documentation
      if (my_id.eq.0) write(*,*) "Calculating the new Di*dIelec_J contribution to wall_curr_i"
      call extended_matrix_vector(my_id, +(1.d0+zeta), sr%d_yieqe, inv_tmp_d_sii, &
                                dIelec_J(:), .false., dwall_curr_i, .false.) 
      !> corresponds to matrix WE_i in documentation
      if (my_id.eq.0) write(*,*) "Calculating the new Ei*dIelec_w contribution to wall_curr_i"
      call extended_matrix_vector(my_id, +(1.d0+zeta), sr%e_yiwe, inv_tmp_d_sii, &
                                dIelec_w(:), .false., dwall_curr_i, .false.) 
      !> corresponds to matrix WI_i in documentation
      if (my_id.eq.0) write(*,*) "Calculating the new Gi*dIelec_w contribution to wall_curr_i"
      call extended_matrix_vector(my_id, tstep*theta, sr%g_yiwe, inv_tmp_d_sii, &
                                dIelec_w(:), .false. , dwall_curr_i, .false. )
      !> corresponds to matrix WH_i in documentation
      if (my_id.eq.0) write(*,*) "Calculating the new Gi*Ielec_w contribution to wall_curr_i"
      call extended_matrix_vector(my_id, tstep, sr%g_yiwe, inv_tmp_d_sii, &
                                Ielec_w(:), .false. , dwall_curr_i, .false. )
      !> corresponds to matrix JP_i in documentation
      if (my_id.eq.0) write(*,*) 'Calling for -sr%d_yieqe*inv_tmp_d_sii*old_dIelec_J'                                
      call extended_matrix_vector(my_id, -zeta, sr%d_yieqe, inv_tmp_d_sii, &
                                old_dIelec_J(:), .false. , dwall_curr_i, .false. )
      !> corresponds to matrix WP_i in documentation
      if (my_id.eq.0) write(*,*) 'Calling for -sr%e_yiwe*inv_tmp_d_sii*old_dIelec_w'                                
      call extended_matrix_vector(my_id, -zeta, sr%e_yiwe, inv_tmp_d_sii, &
                                old_dIelec_w(:), .false. , dwall_curr_i, .false. )
      if (my_id.eq.0) write(*,*) 'Calling for sr%a_yie*inv_tmp_d_sii*old_dpsibnd_vec'                                
      !> corresponds to matrix D_i in documentation
      call extended_matrix_vector(my_id, zeta, sr%a_yie, inv_tmp_d_sii, &
                                old_dpsibnd_vec(:), .false., dwall_curr_i, .true.) 
      if (my_id.eq.0) write(*,*) 'DONE Calling for sr%a_yie*inv_tmp_d_sii*old_dpsibnd_vec' 
      
      dwall_curr_i(:) = dwall_curr_i(:) + dwall_curr_i2(:)
      if (index_now <= 1 ) dwall_curr_i = 0.d0
      wall_curr_i(:) = wall_curr_i(:) + dwall_curr_i(:) 


      if (my_id == 0) then
        OPEN(UNIT=12, FILE="wall_curr_i_full.dat", ACTION="write", STATUS="replace")
        DO k=1,n_wall_curr_i
            WRITE(12,*) wall_curr_i(k)
        END DO
        CLOSE(12)  
      endif

      ! Save wall current in matrix

      wall_curr_i_mat(:,index_now) = wall_curr_i(:)

      if (my_id .eq. 0) then 
        open(unit=10, file="wall_curr_i_mat.txt", status="replace")
        do i = 1, n_wall_curr_i
          write(10,'(*(ES24.16,1X))') (wall_curr_i_mat(i,j), j = 1, n_steps)
        end do
        close(10)
      endif
 
      ! ------------------------------------  CROSSING DOFS ------------------------------------------------

      ! iwe = (FwtFw)^-1 * iel_w

      if (vacuum_debug .and. my_id .eq. 0) then
        write(21,*) 'FwtFw: infos', sr%fwtfw_inv%distrib, sr%fwtfw_inv%row_wise, sr%fwtfw_inv%ind_start, sr%fwtfw_inv%ind_end, sr%fwtfw_inv%step
        21 format(1x, 'distributed', L2, ',rowwise?', L2, ',ind_start of chunk:', I4, ',ind_end of chunk:', I4, ', step' , I4)
      endif

      wall_curr_e(:) = 0.d0

      if (.not. allocated(diag_1)) ALLOCATE(diag_1(n_wall_curr_e))
      diag_1(:) = 1.d0
      call extended_matrix_vector(my_id, 1.d0, sr%fwtfw_inv, diag_1, &
                                  Ielec_w, .true., wall_curr_e, .true.)  

      ! Print out wall_curr_e
      if (my_id == 0) then
        OPEN(UNIT=12, FILE="wall_curr_e.dat", ACTION="write", STATUS="replace")
        DO k=1,sr%n_electrode_w-1
          WRITE(12,*) wall_curr_e(k)
        END DO
        CLOSE(12)      
        OPEN(UNIT=12, FILE="I_elec_w.dat", ACTION="write", STATUS="replace")
        DO k=1,n_wall_curr_e
          WRITE(12,*) Ielec_w(k)
        END DO
        CLOSE(12)
        OPEN(UNIT=12, FILE="I_elec_J.dat", ACTION="write", STATUS="replace")
        DO k=1,n_shell_curr_e
          WRITE(12,*) Ielec_J(k)
        END DO
        CLOSE(12)
      endif
        
    else

      ! ---------------------------------------------------------------------------------------------------
      ! ----------------------------           EDDY COUPLING          ------------------------------------- 
      ! ---------------------------------------------------------------------------------------------------
      if  (.not. allocated(dwall_curr_2)) allocate (dwall_curr_2(n_wall_curr))
      if ( vacuum_debug ) write(*,*) 'wall_curr(before)', sum(abs(wall_curr)), sum(wall_curr)    
        
      if ( (.not. allocated(old_dpsibnd_vec)) .or. (size(old_dpsibnd_vec,1)/= n_dof_starwall) ) then
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
      end if    

      !> WARNING: The following computation of dwall_curr_2 contains the term
      !!          dwall_curr which refers to the previous time-step
      write(*,*) 'n_wall_curr', n_wall_curr
      do k = 1, n_wall_curr
        ! --- contribution from not distributed matrices
        dwall_curr_2(k) = response_d_b(k) * (wall_curr(k) - Y_coils0(k)) &
            + response_d_c(k) * dwall_curr(k)

      end do

      theta = time_evol_theta
      zeta  = time_evol_zeta

      if  (.not. allocated(tmp_d_s)     ) allocate (tmp_d_s     (n_wall_curr))
      if  (.not. allocated(inv_tmp_d_s) ) allocate (inv_tmp_d_s (n_wall_curr))

      tmp_d_s(:) = 1.d0 + zeta + tstep * theta * wall_resistivity * sr%d_yy(:)
      inv_tmp_d_s(:) = 1.0d0 / tmp_d_s(:)

      dwall_curr(:) = 0.0d0

      !> corresponds to matrix A in documentation
      call extended_matrix_vector(my_id, -(1.d0+zeta), sr%a_ye, inv_tmp_d_s, &
                                  dpsibnd_vec(:), .false., dwall_curr, .false.)
      !> corresponds to matrix D in documentation
      call extended_matrix_vector(my_id, zeta, sr%a_ye, inv_tmp_d_s, &
                                  old_dpsibnd_vec(:), .false., dwall_curr, .true.) 
      
      !> Gathering all terms in dwall_curr
      dwall_curr(:) = dwall_curr(:)+dwall_curr_2(:)
      if (index_now <= 1 ) dwall_curr = 0.d0

      !> Evolve wall current
      wall_curr(:) = wall_curr(:) + dwall_curr(:)

    endif  ! ---- Type of current coupling if ------
    ! =================================================================================================


    if (.not. CARIDDI_mode  .and. .not. vacuum_min) then
      call write_wall_vtk(index_now, resistive_wall, my_id)

      if ( vacuum_debug .and. resistive_wall ) then
        call log_wall_curr(my_id)
!        call log_coil_curr()
      end if
    end if

    ! --- Extract net toroidal wall current such that it can be written to the macroscopic_vars.dat file (e.g., for ./util/plot_live_data.sh)
    if ( (.not. allocated(net_tor_wall_curr)) .and. (index_start+nstep >0) ) then
      allocate( net_tor_wall_curr(index_start+nstep) )
      net_tor_wall_curr(:) = 0.d0
    end if
    if (.not. CARIDDI_mode .and. .not. vacuum_min) then
      if (index_now>0) then
        k2 = sr%ncoil + 1
        net_tor_wall_curr(index_now) = sum(sr%s_ww%loc_mat(k2,:) * wall_curr(sr%s_ww%ind_start:sr%s_ww%ind_end))
      endif
      call MPI_ALLReduce(MPI_IN_PLACE, net_tor_wall_curr(index_now),1 ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    endif


    ! --- Extract coil currents such that they can be written to the macroscopic_vars.dat file (e.g., for ./util/plot_live_data.sh)
    if (sr%ncoil > 0) then
      ! Calculate coil currents for all coil types
      if  (.not. allocated(tmp_coil_curr) ) allocate(tmp_coil_curr(sr%ncoil))
      tmp_coil_curr = 0.d0
      if (.not. CARIDDI_mode) then
        do k = 1, sr%ncoil
          tmp_coil_curr(k) = sum(sr%s_ww%loc_mat(k,:) * wall_curr(sr%s_ww%ind_start:sr%s_ww%ind_end))
            
        end do
        call MPI_ALLReduce(MPI_IN_PLACE, tmp_coil_curr,sr%ncoil,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if

      ! Distribute the coil currents to the different coil types
      if ( sr%n_diag_coils > 0 ) then
        if ( (.not. allocated(diag_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( diag_coil_curr(index_start+nstep,sr%n_diag_coils) )
          diag_coil_curr(:,:) = 0.d0
        end if
        if (index_now>0) then
          diag_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_diag_coils:sr%ind_start_diag_coils + sr%n_diag_coils -1) 
        endif
      end if
      
      if ( sr%n_rmp_coils > 0 ) then
        if ( (.not. allocated(rmp_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( rmp_coil_curr(index_start+nstep,sr%n_rmp_coils) )
          rmp_coil_curr(:,:) = 0.d0
        end if
        if (index_now>0) then
          rmp_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_rmp_coils : sr%ind_start_rmp_coils + sr%n_rmp_coils -1) 
        endif
      end if
      
      if ( sr%n_pol_coils > 0 ) then
        if ( (.not. allocated(pf_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( pf_coil_curr(index_start+nstep,sr%n_pol_coils) )
          pf_coil_curr(:,:) = 0.d0
        end if
        if (index_now>0) then
          pf_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_pol_coils:sr%ind_start_pol_coils + sr%n_pol_coils -1) 
        endif
      end if

      ! Voltage coils not yet implemented
      !if ( sr%n_voltage_coils > 0 ) then
      !  if ( (.not. allocated(voltage_coil_curr)) .and. (index_start+nstep >0) ) then
      !    allocate( voltage_coil_curr(index_start+nstep,sr%n_voltage_coils) )
      !    voltage_coil_curr(:,:) = 0.d0
      !  end if
      !  if (index_now>0) then
      !    voltage_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_voltage_coils:sr%ind_start_voltage_coils + sr%n_voltage_coils -1) 
      !  endif
      !end if
    endif


    if ( vacuum_debug .and. (my_id ==0) ) write(*,*) 'wall_curr(after)', sum(abs(wall_curr)), sum(wall_curr)
    
  end subroutine evolve_wall_currents
  

  !> Reconstruct the potential values at the wall triangle nodes.
  !!
  !! To reconstruct the physical wall potentials from the ones we work 
  !! with we need to muliply them by the similarity transform matrix.
  !! The physical wall potentials include different types of potentials,
  !! and the way they are ordered in the array is
  !!
  !!   (I_coil_1, I_coil_2, ..., I_coil_ncoil, Iw_net_tor, Potw_1, Potw_2, ..., Potw_npotw-1)
  !!  
  !! where I_coil are the coil currents, Iw_net is the net wall current
  !! and Potw are the single valued wall potentials.  
  subroutine reconstruct_triangle_potentials(tripot_w, wall_curr, my_id, Iw_net_tor)
    
    use mpi_mod

    implicit none
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: tripot_w(:)
    real*8, allocatable, intent(in)    :: wall_curr(:)
    integer,             intent(in)    :: my_id
    real*8,  optional,   intent(inout) :: Iw_net_tor 

    ! --- Local variables
    integer              :: i, j, ierr
    integer              :: count=1
    real*8, allocatable  :: pot_tmp(:)

    if ( allocated(tripot_w) ) deallocate(tripot_w); allocate( tripot_w(sr%npot_w) )    
    if ( allocated(pot_tmp)  ) deallocate(pot_tmp);  allocate(  pot_tmp(sr%npot_w) )    
    tripot_w = 0.d0
    pot_tmp  = 0.d0

    ! --- multiply by the similarity transform matrix to get the physical wall potentials
    if ( allocated(wall_curr) ) then

      do i = 1, sr%npot_w
        j = i + sr%ncoil
        pot_tmp(i) = sum(sr%s_ww%loc_mat(j,:) * wall_curr(sr%s_ww%ind_start:sr%s_ww%ind_end))
      end do
     
    end if

    call MPI_AllREDUCE(MPI_IN_PLACE,pot_tmp,size(pot_tmp),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    if (present(Iw_net_tor)) then
      Iw_net_tor = pot_tmp(1)
    endif

    ! --- Correct indexing (1st potential is the net curret potential)
    ! --- Shift indexing so 1st wall node corresponds to 1st single valued potential (and so on)
    do i = 1, sr%npot_w-1
      tripot_w(i) = pot_tmp(i+1)
    enddo

    tripot_w(sr%npot_w) = 0.d0 ! STARWALL's BC for the single valued potentials (last potential is 0)
    !!! Net current potential contribution to be added soon!

    deallocate(pot_tmp) 

  end subroutine reconstruct_triangle_potentials
  





  !> Reconstruct the coil potentials
  !!
  !! To reconstruct the physical potentials from the ones we work 
  !! with we need to muliply them by the similarity transform matrix.
  !! The physicall potentials include different types of potentials,
  !! and the way they are ordered in the wall_curr array is
  !!
  !!   (I_coil_1, I_coil_2, ..., I_coil_ncoil, Iw_net_tor, Potw_1, Potw_2, ..., Potw_npotw-1)
  !!  
  !! where I_coil are the coil currents, Iw_net is the net wall current
  !! and Potw are the single valued wall potentials.  
  subroutine reconstruct_coil_potentials(pot_c, wall_curr, my_id)
    
    use mpi_mod

    implicit none
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: pot_c(:)
    real*8, allocatable, intent(in)    :: wall_curr(:)
    integer,             intent(in)    :: my_id

    ! --- Local variables
    integer              :: i, j, ierr
    integer              :: count=1

    if (sr%ncoil < 1 ) return

    if ( allocated(pot_c) ) deallocate(pot_c); allocate( pot_c(sr%ncoil) )    
    pot_c = 0.d0

    ! --- multiply by the similarity transform matrix to get the physical wall potentials
    if ( allocated(wall_curr) ) then

      do i = 1, sr%ncoil
        pot_c(i) = sum(sr%s_ww%loc_mat(i,:) * wall_curr(sr%s_ww%ind_start:sr%s_ww%ind_end))
      end do
     
    end if

    call MPI_AllREDUCE(MPI_IN_PLACE,pot_c,size(pot_c),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  end subroutine reconstruct_coil_potentials
  
  



  
  !> Write wall current potentials to logfile.
  subroutine log_wall_curr(my_id)
    
    implicit none
    
    integer,   intent(in) :: my_id

    ! --- Local variables
    real*8, allocatable :: tripot_w(:)
    
    call reconstruct_triangle_potentials(tripot_w, wall_curr, my_id)
    if (my_id==0) write(*,'(" Wall current potentials (min, max): ",ES16.8,"...",ES16.8)') minval(tripot_w), maxval(tripot_w)
    deallocate( tripot_w )
    
  end subroutine log_wall_curr

  
  !> Calculate matrix response_m_eq needed for freeboundary equilibrium
  subroutine init_vacuum_response(my_id,  freeboundary_equil)

    use mpi_mod
    use phys_module, only: restart
    implicit none

    ! --- Routine parameters
    integer,                     intent(in) :: my_id              !< MPI task ID
    logical,                     intent(in) :: freeboundary_equil !< Resistive or ideal wall?

    ! --- Local variables
    integer :: i, j, k, j2, k2
    integer :: ierr

    if ( sr%n_tor == 0 ) then
      write(*,*) 'Remark: Routine init_vacuum_response is not doing anything since sr%n_tor==0.'
      sr%initialized = .true.
      return
    end if

    if ( .not. freeboundary_equil ) then
      write(*,*) 'Remark: Routine init_vacuum_response is not doing anything since freeboundary_equil=.false..'
      sr%initialized = .true.
      return
    end if

    if ( restart ) then
      if (my_id .eq. 0) & 
          write(*,*) 'response_m_eq only required in equilibrium, return.'
      sr%initialized = .true.
      return
    end if
    
    if (sr%initialized) return

    if (my_id .eq. 0) write(*,*) 'INIT vacuum response, calculate response_m_eq'
    ! --- Allocate matrices if required
    if ( .not. allocated(response_m_eq)) allocate( response_m_eq(sr%nd_bez/sr%n_tor, sr%nd_bez/sr%n_tor) )
    ! --- Derived response matrix for equilibrium (extract n=0 part from STARWALL EE matrix)
    response_m_eq = 0.d0

    do j = 1, sr%nd_bez/sr%n_tor0, 2
      j2 = (j-1)*sr%n_tor0 + 1
      do k = 1, sr%nd_bez/sr%n_tor0, 2
        k2 = (k-1)*sr%n_tor0 + 1

        if ( sr%a_ee%row_wise ) then

          if ( (j2>=sr%a_ee%ind_start) .and. (j2<=sr%a_ee%ind_end) ) then
            response_m_eq(j,k:k+1) = sr%a_ee%loc_mat(j2-sr%a_ee%step*my_id,k2:k2+1)
          end if

          if( ((j2+1)>=sr%a_ee%ind_start) .and. ((j2+1)<=sr%a_ee%ind_end) ) then
            response_m_eq(j+1,k:k+1) = sr%a_ee%loc_mat((j2+1)-sr%a_ee%step*my_id,k2:k2+1)
          end if

        end if

      end do
    end do

    call MPI_AllREDUCE(MPI_IN_PLACE,response_m_eq,size(response_m_eq),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    sr%initialized=.true.

  end subroutine init_vacuum_response

 
  !> Update vacuum response
  !! This is necessary
  !! - right after the start or restart of the code
  !! - when wall resistivity, tstep, or some other parameters have changed.
  subroutine update_response(my_id, tstep, resistive_wall, node_list, element_list, bnd_elm_list, bnd_node_list)

    use phys_module, only: time_evol_theta, time_evol_zeta
    use mpi_mod
    use data_structure, only: type_node_list, type_node, type_element_list, type_element, type_bnd_node_list,type_bnd_element_list, type_bnd_element, type_SP_MATRIX

    implicit none

    ! --- Routine parameters
    integer,                     intent(in) :: my_id              !< MPI task ID
    real*8,                      intent(in) :: tstep              !< delta t,timestep
    logical,                     intent(in) :: resistive_wall     !< Resistive or ideal wall?
    type(type_node_list), OPTIONAL,       intent(in) :: node_list            !< List of grid nodes
    type(type_bnd_node_list), OPTIONAL,    intent(in) :: bnd_node_list        !< List of boundary grid nodes
    type(type_bnd_element_list),OPTIONAL, intent(in) :: bnd_elm_list         !< List of boundary elements
    type(type_element_list), OPTIONAL,    intent(in) :: element_list

    ! --- Local variables
    integer :: i, j, k, j2, k2, leq, lw,nnw
    real*8, allocatable :: tmp_d_s(:)
    real*8, ALLOCATABLE :: tmp_d_sii(:)
    integer :: dim_wh_i(2), dim_jp_i(2), dim_wp_i(2)
    real*8  :: theta, zeta
    logical :: update_required
    integer :: ierr,loc_size,ntasks
    real*8  :: test_sum,test_sum2
    real*8  :: Jj_vec_tot(n_shell_curr_e+1), Jj_vec(n_shell_curr_e)
    real*8,ALLOCATABLE :: diag_el_eq(:), u_h_vec(:), response_m_w_temp(:,:), Jmat_w_rank(:,:), Jmat_J_rank(:,:)
    real*8,allocatable :: response_m_w_h(:,:), response_m_u_h(:,:), response_m_t_h(:,:)
    real*8,allocatable :: response_m_jewj_i(:,:), response_m_wtotww_i(:,:), response_m_t_h_cpsi_i(:,:)
    real*8,ALLOCATABLE :: response_m_JEWJ_CPSI_I(:,:), response_m_wtotww_cpsi_i(:,:)
    real*8,ALLOCATABLE :: response_curr_contrib_dpsi(:,:), response_curr_contrib_dzj(:,:)
    real*8,ALLOCATABLE :: response_m_JEWJ_czj_i(:,:), response_m_wtotww_czj_i(:,:), response_m_t_h_czj_i(:,:)
    real*8,ALLOCATABLE :: response_R_psis(:,:), response_R_zj(:,:)
    real*8,ALLOCATABLE :: old_wh_i(:,:), new_wh_i(:,:)
    type(t_distrib_mat)  :: response_m_d                   !< \f$\hat{D}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_a                   !< \f$\hat{A}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_d_i                 !< \f$\hat{D_i}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_a_i                 !< \f$\hat{A_i}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_JE_i                !< \f$\hat{JE_i}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_WE_i                !< \f$\hat{WE_i}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_WI_i                !< \f$\hat{WI_i}\f$ in the documentation
    type(t_distrib_mat)  :: temp_m_WH_i
    type(t_distrib_mat)  :: temp_m_JP_i
    type(t_distrib_mat)  :: temp_m_WP_i
    type(t_distrib_mat)  :: response_m_WH_i                !< \f$\hat{WH_i}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_JP_i                !< \f$\hat{JP_i}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_WP_i                !< \f$\hat{WP_i}\f$ in the documentation
    type(t_distrib_mat)  :: response_m_Wtot_i   
    type(t_distrib_mat)  :: Jmat_J_par 
    type(t_distrib_mat)  :: Jmat_w_par
           

! --- Local variables to store the previous values of some parameters.
    real*8,  save :: old_thick   = 0.0
    real*8,  save :: old_res     = 0.0
    real*8,  save :: old_tstep   = 0.0
    real*8,  save :: old_theta   = 0.0
    real*8,  save :: old_zeta    = 0.0
    logical, save :: old_reswall = .false.

    if ( sr%n_tor == 0 ) then
      write(*,*) 'Remark: Routine update_response is not doing anything since sr%n_tor==0.'
      return
    end if

    theta = time_evol_theta
    zeta  = time_evol_zeta

    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

! --- Update response matrices only, if parameter values changed or matrices not allocated
    if (CARIDDI_mode .and. halo_freeb) then
      update_required = ( old_res     /= wall_resistivity       ) &
                 .or. ( old_tstep   /= tstep                  ) &
                 .or. ( old_theta   /= theta                  ) &
                 .or. ( old_zeta    /= zeta                   ) &
                 .or. ( old_reswall .neqv. resistive_wall     ) &
                 .or. ( .not. allocated(response_d_b_i)       ) &
                 .or. ( .not. allocated(response_d_c_i)       ) &
                 .or. ( .not. allocated(response_m_h)         ) &
                 .or. ( .not. allocated(response_m_j_i)         ) &
                 .or. ( .not. allocated(response_m_eh)         )   &
                 .or. ( .not. allocated(response_m_zh)           )   &
                 .or. ( .not. allocated(response_m_ph)           )   &
                 .or. ( .not. allocated(response_m_qh)           )   &
                 .or. ( .not. allocated(response_m_vh)           )   &
                 .or. ( .not. allocated(response_m_xh)           )   
      update_required = .true.
                       
    else
      update_required = ( old_res     /= wall_resistivity       ) &
                 .or. ( old_tstep   /= tstep                  ) &
                 .or. ( old_theta   /= theta                  ) &
                 .or. ( old_zeta    /= zeta                   ) &
                 .or. ( old_reswall .neqv. resistive_wall     ) &
                 .or. ( .not. allocated(response_d_b)         ) &
                 .or. ( .not. allocated(response_d_c)         ) &
                 .or. ( .not. allocated(response_m_e)         ) &
                 .or. ( .not. allocated(response_m_h)         ) &
                 .or. ( .not. allocated(response_m_j)         )

    endif


    if ( update_required ) then
      ! --- Remember parameter values.
      old_res     = wall_resistivity
      old_tstep   = tstep
      old_theta   = theta
      old_zeta    = zeta
      old_reswall = resistive_wall

      ! --- Allocate matrices if required
      if (CARIDDI_mode .and. halo_freeb) then
        if ( .not. allocated(response_d_b_i) ) allocate( response_d_b_i(n_wall_curr_i) )
        if ( .not. allocated(response_d_c_i) ) allocate( response_d_c_i(n_wall_curr_i) )
        if ( .not. allocated(response_m_h) ) allocate( response_m_h(n_dof_starwall, n_dof_starwall) )
        if( .not. allocated (response_m_a_i%loc_mat) ) &
            call  alloc_distr(my_id, response_m_a_i, (/n_wall_curr_i, n_dof_starwall/), .false.)
        if( .not. allocated (response_m_d_i%loc_mat) ) &
            call  alloc_distr(my_id, response_m_d_i, (/n_wall_curr_i, n_dof_starwall/), .false.)
        if( .not. allocated (response_m_JE_i%loc_mat) ) &
            call  alloc_distr(my_id, response_m_JE_i, (/n_wall_curr_i, n_shell_curr_e/), .true.)
        if( .not. allocated (response_m_WE_i%loc_mat) ) &
            call  alloc_distr(my_id, response_m_WE_i, (/n_wall_curr_i, n_wall_curr_e/), .true.)
        if( .not. allocated (response_m_WI_i%loc_mat) ) &
            call  alloc_distr(my_id, response_m_WI_i, (/n_wall_curr_i, n_wall_curr_e/), .true.)
        if( .not. allocated (response_m_Wtot_i%loc_mat) ) &
            call  alloc_distr(my_id, response_m_Wtot_i, (/n_wall_curr_i, n_wall_curr_e/), .true.)

        if( .not. allocated (temp_m_WH_i%loc_mat) ) &
            call  alloc_distr(my_id, temp_m_WH_i, (/n_wall_curr_i, n_wall_curr_e/), .true.)
        if (.not. allocated(old_wh_i)) allocate(old_wh_i(n_wall_curr_i,n_wall_curr_e))
        if (.not. allocated(new_wh_i)) allocate(new_wh_i(n_wall_curr_i,n_wall_curr_e))

        ! if( .not. allocated (response_m_WH_i%loc_mat) ) &
        !     call  alloc_distr(my_id, response_m_WH_i, (/n_wall_curr_i, n_wall_curr_e/), .false.)

        if( .not. allocated (temp_m_JP_i%loc_mat) ) &
            call  alloc_distr(my_id, temp_m_JP_i, (/n_wall_curr_i, n_shell_curr_e/), .true.)
        ! if( .not. allocated (response_m_JP_i%loc_mat) ) &
        !     call  alloc_distr(my_id, response_m_JP_i, (/n_wall_curr_i, n_shell_curr_e/), .false.)

        if( .not. allocated (temp_m_WP_i%loc_mat) ) &
            call  alloc_distr(my_id, temp_m_WP_i, (/n_wall_curr_i, n_wall_curr_e/), .true.)
        ! if( .not. allocated (response_m_WP_i%loc_mat) ) &
        !     call  alloc_distr(my_id, response_m_WP_i, (/n_wall_curr_i, n_wall_curr_e/), .false.)

        if( .not. allocated (Jmat_J_par%loc_mat) ) &
            call  alloc_distr(my_id, Jmat_J_par, (/ n_shell_curr_e , n_dof_starwall /), .true.)
        if( .not. allocated (Jmat_w_par%loc_mat) ) &
            call  alloc_distr(my_id, Jmat_w_par, (/ n_wall_curr_e , n_dof_starwall /), .true.)
      else
        if ( .not. allocated(response_d_b) ) allocate( response_d_b(n_wall_curr) )
        if ( .not. allocated(response_d_c) ) allocate( response_d_c(n_wall_curr) )
        if ( .not. allocated(response_m_h) ) allocate( response_m_h(n_dof_starwall, n_dof_starwall) )
        if( .not. allocated (response_m_a%loc_mat) ) &
            call  alloc_distr(my_id, response_m_a, (/n_wall_curr, n_dof_starwall/), .false.)
        if( .not. allocated (response_m_d%loc_mat) ) &
            call  alloc_distr(my_id, response_m_d, (/n_wall_curr, n_dof_starwall/), .false.)
      endif

      ! ======================================================================================
      ! ------------------ Derived response matrices for time-evolution ----------------------
      ! ======================================================================================

      if ( resistive_wall) then

        ! ======================================= EDDY+HALO COUPLING =====================================
        if (halo_freeb .and. CARIDDI_mode) then

          if (.not. ALLOCATED(tmp_d_sii)) ALLOCATE(tmp_d_sii(n_wall_curr_i))

          tmp_d_sii(:) = 1.d0 + zeta + tstep * theta * wall_resistivity * sr%dii_yy(:)

          response_d_b_i(:) = - tstep * wall_resistivity * sr%dii_yy / tmp_d_sii(:)
          response_d_c_i(:) = zeta / tmp_d_sii(:)

          !!!! CALL ROUTINE FOR HALO TIME STEPPING !!!
          if (present(node_list)) then 
            if (.not. ALLOCATED(response_R_psis)) allocate(response_R_psis(sr%nd_bez, sr%nd_bez))
            if (.not. ALLOCATED(response_R_zj)) allocate(response_R_zj(sr%nd_bez, sr%nd_bez))
            call halo_time_stepping(my_id, node_list, element_list, bnd_elm_list, bnd_node_list, response_R_zj=response_R_zj, response_R_psis=response_R_psis, debug=.false.)
          
            ! --- m_w_h = a_eye * fwtfw^-1 * J_w
            if (.not. allocated(Jmat_w_rank)) allocate(Jmat_w_rank(sr%n_electrode_w-1,sr%nd_bez) )
            Jmat_w_rank = sr%Jmat_w(1:sr%n_electrode_w-1,:)
            ! ! Define parallelized Jmat_w
            ! do lw=1,size(Jmat_w_par%loc_mat,1)
            !   Jmat_w_par%loc_mat(lw,:) = Jmat_w_rank(lw+my_id*Jmat_J_par%step,:)
            ! enddo 

            call matrix_multiplication(my_id, sr%w_h, mat2_not_distr=Jmat_w_rank, res_mat_not_distr=response_m_w_h)
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_w_h,size(response_m_w_h),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ! Minus sign as the wall-impinging current is considered with opposite normal
            response_m_w_h(:,:) = -response_m_w_h(:,:)

            ! --- m_u_h = - ael_eq * J_j
            if (.not. allocated(Jmat_J_rank)) allocate(Jmat_J_rank(n_shell_curr_e,sr%nd_bez) )
            Jmat_J_rank = sr%Jmat_J(1:n_shell_curr_e,:)
            
            !   ! Define parallelized Jmat_J
            !   do leq=1,size(Jmat_J_par%loc_mat,1)
            !     Jmat_J_par%loc_mat(leq,:) = Jmat_J_rank(leq+my_id*Jmat_J_par%step,:)
            !   enddo 

            call matrix_multiplication(my_id, sr%aeq_el, mat2_not_distr=Jmat_J_rank, res_mat_not_distr=response_m_u_h)
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_u_h,size(response_m_u_h),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            response_m_u_h(:,:) = -response_m_u_h(:,:)
            
            ! --- m_t_h = m_u_h + m_w_h
            if (.not. ALLOCATED(response_m_t_h)) allocate(response_m_t_h(sr%nd_bez, sr%nd_bez))
            response_m_t_h(:,:) = response_m_w_h(:,:) + response_m_u_h(:,:)

            ! N.B. Change of sign compared to eddy-only coupl. is performed after reading, the Ci_star in CARIDDI already contains "-"
            do j = 1, size(response_m_a_i%loc_mat,2) 
              response_m_a_i%loc_mat(:,j) = -(1.d0+zeta) * sr%a_yie%loc_mat(:,j) / tmp_d_sii(:)
            end do
            ! N.B. Change of sign compared to eddy-only coupl. is performed after reading, the Ci_star in CARIDDI already contains "-"
            do j = 1, size(response_m_d_i%loc_mat,2)
              response_m_d_i%loc_mat(:,j) = zeta * sr%a_yie%loc_mat(:,j) / tmp_d_sii(:)
            end do
            
            ! --- Calculate response_m_JE_i
            do nnw = 1,size(response_m_JE_i%loc_mat,1) 
              response_m_JE_i%loc_mat(nnw,:) = (1+zeta) * sr%d_yieqe%loc_mat(nnw,:) / tmp_d_sii(nnw+my_id*response_m_JE_i%step) 
            end do
            
            ! --- Calculate response_m_WE_i 
            do nnw = 1,size(response_m_WE_i%loc_mat,1)
              response_m_WE_i%loc_mat(nnw,:) = (1+zeta) * sr%e_yiwe%loc_mat(nnw,:) / tmp_d_sii(nnw+my_id*response_m_WE_i%step) 
            end do

            ! --- Calculate response_m_WI_i 
            do nnw = 1,size(response_m_WI_i%loc_mat,1)
              response_m_WI_i%loc_mat(nnw,:) = tstep*theta * sr%g_yiwe%loc_mat(nnw,:) / tmp_d_sii(nnw+my_id*response_m_WI_i%step) 
            end do

            ! --- Sum up WE and WI
            response_m_Wtot_i%loc_mat(:,:) = response_m_WE_i%loc_mat(:,:) + response_m_WI_i%loc_mat(:,:)

            ! ! --- Calculate response_m_WH_i 

            ! 1) Calculate temp. row-distributed matrix
            do nnw = 1,size(temp_m_WH_i%loc_mat,1)
              temp_m_WH_i%loc_mat(nnw,:) = tstep * sr%g_yiwe%loc_mat(nnw,:) / tmp_d_sii(nnw+my_id*temp_m_WH_i%step) 
            end do

            ! 2) Make matrix column-distributed
            dim_WH_i = (/  n_wall_curr_i , n_wall_curr_e  /)
            
            if (my_id .eq. 0) write(*,*) 'Calling row_col_distrib for WH...'
            write(*,*) 'temp_m_WH_i%step and my_id', temp_m_WH_i%step , my_id
            call row_col_distrib(my_id, ntasks, dim_WH_i, temp_m_WH_i, response_m_WH_i , row_wise = .false.)
            if (my_id .eq. 0) write(*,*) 'Done calling row_col_distrib for WH!'
            write(*,*) 'response_m_WH_i%step and my_id', response_m_WH_i%step , my_id
            write(*,*) 'response_m_WH_i%dims and my_id', response_m_WH_i%dim(1) , response_m_WH_i%dim(2),  my_id
          
            ! --- Calculate response_m_JP_i 
            do nnw = 1,size(temp_m_JP_i%loc_mat,1)
              temp_m_JP_i%loc_mat(nnw,:) = - zeta * sr%d_yieqe%loc_mat(nnw,:) / tmp_d_sii(nnw+my_id*temp_m_JP_i%step) 
            end do

            dim_JP_i = (/  n_wall_curr_i , n_shell_curr_e  /)
            
            if (my_id .eq. 0) write(*,*) 'Calling row_col_distrib for JP...'
            call row_col_distrib(my_id, ntasks, dim_JP_i, temp_m_JP_i, response_m_JP_i , row_wise = .false.)
            if (my_id .eq. 0) write(*,*) 'Done calling row_col_distrib for JP!'

            ! --- Calculate response_m_WP_i 
            do nnw = 1,size(temp_m_WP_i%loc_mat,1)
              temp_m_WP_i%loc_mat(nnw,:) = -zeta * sr%e_yiwe%loc_mat(nnw,:) / tmp_d_sii(nnw+my_id*temp_m_WP_i%step) 
            end do

            dim_WP_i = (/  n_wall_curr_i , n_wall_curr_e  /)
            
            if (my_id .eq. 0) write(*,*) 'Calling row_col_distrib for WP...'
            call row_col_distrib(my_id, ntasks, dim_WP_i, temp_m_WP_i, response_m_WP_i , row_wise = .false. )
            if (my_id .eq. 0) write(*,*) 'Done calling row_col_distrib for WP!'

            ! Calculate product JE*W^J
            if (my_id .eq. 0) write(*,*) 'sr%response_m_JE_i row_wise?', response_m_JE_i%row_wise
            call matrix_multiplication(my_id, response_m_JE_i, mat2_not_distr=Jmat_J_rank, res_mat_not_distr=response_m_JEWJ_i)
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_jewj_i,size(response_m_jewj_i),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ! Calculate product Wtot*W^w 
            if (my_id .eq. 0) write(*,*) 'sr%response_m_Wtot_i row_wise?', response_m_Wtot_i%row_wise
            call matrix_multiplication(my_id, response_m_Wtot_i, mat2_not_distr=Jmat_w_rank, res_mat_not_distr=response_m_wtotww_i)
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_wtotww_i,size(response_m_wtotww_i),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

            !! CALCULATING ALL THE PIECES MAKING UP M_EH !!!

            ! 1) response_m_JEWJ_i * response_R_psis
            if (.not. allocated(response_m_JEWJ_CPSI_I)) allocate(response_m_JEWJ_CPSI_I(n_wall_curr_i,sr%nd_bez))
            response_m_JEWJ_CPSI_I(:,:) = 0.d0
            response_m_JEWJ_CPSI_I = MATMUL(response_m_jewj_i, response_R_psis )
            ! 2) response_m_Wtotww_i * response_R_psis
            if (.not. allocated(response_m_wtotww_cpsi_i)) allocate(response_m_wtotww_cpsi_i(n_wall_curr_i,sr%nd_bez))
            response_m_wtotww_cpsi_i(:,:) = 0.d0
            response_m_wtotww_cpsi_i = MATMUL(response_m_wtotww_i, response_R_psis)
            ! 3) response_m_t_h * response_R_psis
            if (.not. allocated(response_m_t_h_cpsi_i)) allocate(response_m_t_h_cpsi_i(sr%nd_bez,sr%nd_bez))
            response_m_t_h_cpsi_i(:,:) = 0.d0
            response_m_t_h_cpsi_i = MATMUL(response_m_t_h , response_R_psis)

            ! Create response_curr_contrib_dpsi 
            if (.not. allocated(response_curr_contrib_dpsi)) allocate(response_curr_contrib_dpsi(n_wall_curr_i, n_dof_starwall))
            response_curr_contrib_dpsi = 0.d0
!            do j = 1, size(response_m_a_i%loc_mat,2)  ! m_a_i is column-wise
            response_curr_contrib_dpsi(:,sr%a_ee%ind_start:sr%a_ee%ind_end) = response_m_a_i%loc_mat(:,:) !+ response_m_JEWJ_CPSI_I(:,j+my_id*response_m_a_i%step)  &  
                                                                        !- response_m_wtotww_cpsi_i(:,j+my_id*response_m_a_i%step)
 !           enddo
            
            call MPI_AllREDUCE(MPI_IN_PLACE,response_curr_contrib_dpsi,size(response_curr_contrib_dpsi),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

            ! WARNING: The minus sign  " - response_m_wtotww_cpsi_i " exists because we assume a minus sign for the wall normal!           
            response_curr_contrib_dpsi(:,:) = response_curr_contrib_dpsi(:,:)  &
                                            + response_m_JEWJ_CPSI_I(:,:) - response_m_wtotww_cpsi_i(:,:)

            ! CALCULATE RESPONSE_M_EH
            ! --- m_eh_i = a_ee + matmul[ a_eyi, ( response_m_a_i + response_m_JEWJ_cpsi_i - response_m_wtotww_cpsi_i ) ] + matmul(response_m_t_h,Cpsi)

            if (my_id .eq. 0) write(*,*) 'Halo field: calculating response_m_eh...'
            call matrix_multiplication(my_id,sr%a_eyi,mat2_not_distr=response_curr_contrib_dpsi, res_mat_not_distr=response_m_eh)
            if (my_id .eq. 0) write(*,*) 'Matrix multiplication done!'
            do j = 1, size(sr%a_ee%loc_mat,1)
              response_m_eh(j+my_id*sr%a_ee%step,:) = response_m_eh(j+my_id*sr%a_ee%step,:) + sr%a_ee%loc_mat(j,:) + &
                                                      response_m_t_h_cpsi_i(j+my_id*sr%a_ee%step,:)
            end do
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_eh,size(response_m_eh),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            if (my_id .eq. 0) write(*,*) 'MPI_AllREDUCE done for response_m_eh!'

            if (my_id .eq. 0) write(*,*) 'Final sum done!'
            if (my_id .eq. 0) write(*,*) '...DONE calculating response_m_eh!'          


            !!! CALCULATING ALL THE PIECES MAKING UP M_Z!!!

            ! 1) response_m_JEWJ_i * response_R_zj
            if (.not. allocated(response_m_JEWJ_Czj_I)) allocate(response_m_JEWJ_Czj_I(n_wall_curr_i,sr%nd_bez))
            response_m_JEWJ_Czj_I(:,:) = 0.d0
            response_m_JEWJ_Czj_I = MATMUL(response_m_jewj_i, response_R_zj )
            ! 2) response_m_Wtotww_i * response_R_zj
            if (.not. allocated(response_m_wtotww_czj_i)) allocate(response_m_wtotww_czj_i(n_wall_curr_i,sr%nd_bez))
            response_m_wtotww_czj_i(:,:) = 0.d0
            response_m_wtotww_czj_i = MATMUL(response_m_wtotww_i, response_R_zj)
            ! 3) response_m_t_h * response_R_zj
            if (.not. allocated(response_m_t_h_czj_i)) allocate(response_m_t_h_czj_i(sr%nd_bez,sr%nd_bez))
            response_m_t_h_czj_i(:,:) = 0.d0
            response_m_t_h_czj_i = MATMUL(response_m_t_h , response_R_zj)
            
            ! Create response_curr_contrib_dzj 
            if (.not. allocated(response_curr_contrib_dzj)) allocate(response_curr_contrib_dzj(n_wall_curr_i, n_dof_starwall))
            response_curr_contrib_dzj = 0.d0
            ! WARNING: The minus sign  " - response_m_wtotww_czj_i " exists because we assume a minus sign for the wall normal!
            response_curr_contrib_dzj(:,:) = + response_m_JEWJ_Czj_I(:,:) - response_m_wtotww_czj_i(:,:)

            ! CALCULATE response_m_zh
            if (my_id .eq. 0) write(*,*) 'Halo field: calculating response_m_zh...'
            call matrix_multiplication(my_id,sr%a_eyi,mat2_not_distr=response_curr_contrib_dzj, res_mat_not_distr=response_m_zh)
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_zh,size(response_m_zh),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            response_m_zh(:,:) = response_m_zh(:,:) + response_m_t_h_czj_i(:,:)
            if (my_id .eq. 0) write(*,*) '...DONE calculating response_m_zh!'

            ! RESPONSE_M_P_H = matmul(M_eyi * WH) + M_e,elw
            if (my_id .eq. 0) write(*,*) 'Halo field: calculating response_m_ph...'
            if (my_id .eq. 0) write(*,*) 'response_m_wh_i dims', response_m_WH_i%dim(1) , response_m_wh_i%dim(2)
            call matrix_multiplication(my_id,sr%a_eyi,mat2=response_m_WH_i, res_mat_not_distr=response_m_PH)
            do j = 1, size(sr%a_ee%loc_mat,1)
              response_m_ph(j+my_id*sr%a_ee%step,:) = response_m_ph(j+my_id*sr%a_ee%step,:) + sr%w_h%loc_mat(j,:)
            end do
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_ph,size(response_m_ph),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            if (my_id .eq. 0) write(*,*) '...DONE calculating response_m_ph!'
    
            ! RESPONSE_M_Q_H = matmul(M_eyi * WP) 
            if (my_id .eq. 0) write(*,*) 'Halo field: calculating response_m_qh...'
            call matrix_multiplication(my_id,sr%a_eyi,mat2=response_m_WP_i, res_mat_not_distr=response_m_qh)
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_qh,size(response_m_qh),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            if (my_id .eq. 0) write(*,*) '...DONE calculating response_m_qh!'

            ! RESPONSE_M_X_H = - M_e,eleq 
            if (my_id .eq. 0) write(*,*) 'Halo field: calculating response_m_xh...'
            if (.not. allocated(response_m_xh)) allocate(response_m_xh(sr%nd_bez, n_shell_curr_e))
            response_m_xh(:,:) =0.0
            response_m_xh(sr%a_ee%ind_start:sr%a_ee%ind_end,:) = - sr%aeq_el%loc_mat(:,:)
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_xh,size(response_m_xh),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            if (my_id .eq. 0) write(*,*) '...DONE calculating response_m_xh!'

            ! RESPONSE_M_V_H = matmul (M_eyi * JP) 
            if (my_id .eq. 0) write(*,*) 'Halo field: calculating response_m_vh...'
            call matrix_multiplication(my_id,sr%a_eyi,mat2=response_m_JP_i, res_mat_not_distr=response_m_vh)        
            call MPI_AllREDUCE(MPI_IN_PLACE,response_m_vh,size(response_m_vh),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
            if (my_id .eq. 0) write(*,*) '...DONE calculating response_m_vh!'
          endif

          response_m_h =0.0
          response_m_h(sr%a_ee%ind_start:sr%a_ee%ind_end,:) = sr%a_ee%loc_mat(:,:)
          call MPI_AllREDUCE(MPI_IN_PLACE,response_m_h,size(response_m_h),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

          ! --- m_j_i =  matmul(a_eyi, m_d_i)
          call matrix_multiplication(my_id,sr%a_eyi,mat2=response_m_d_i,res_mat_not_distr=response_m_j_i)
          call MPI_AllREDUCE(MPI_IN_PLACE,response_m_j_i,size(response_m_j_i),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

          call  dealloc_distr(response_m_a_i)
          call  dealloc_distr(response_m_d_i)
          call  dealloc_distr(response_m_WE_i)
          call  dealloc_distr(response_m_WI_i)
          call  dealloc_distr(response_m_JE_i) 
          call  dealloc_distr(response_m_Wtot_i)  

          deallocate( tmp_d_sii)

        else
        ! ------------------------------- EDDY COUPLING ------------------------------------------------

          if (.not. ALLOCATED(tmp_d_s)) allocate( tmp_d_s(n_wall_curr) )

          tmp_d_s(:) = 1.d0 + zeta + tstep * theta * wall_resistivity * sr%d_yy(:)

          response_d_b(:) = - tstep * wall_resistivity * sr%d_yy(:) / tmp_d_s(:)
          response_d_c(:) = zeta / tmp_d_s(:)
          
          !response_d_b_i(:) = 0.d0
          !response_d_c_i(:) = 0.d0 

          do j = 1, size(response_m_a%loc_mat,2) 
            response_m_a%loc_mat(:,j) = -(1.d0+zeta) * sr%a_ye%loc_mat(:,j) / tmp_d_s(:)
          end do

          do j = 1, size(response_m_d%loc_mat,2)
            response_m_d%loc_mat(:,j) = zeta * sr%a_ye%loc_mat(:,j) / tmp_d_s(:)
          end do
          
          ! --- m_e = a_ee - matmul(a_ey, m_a)
          call matrix_multiplication(my_id,sr%a_ey,mat2=response_m_a, res_mat_not_distr=response_m_e)
      
          do j = 1, size(sr%a_ee%loc_mat,1)
            response_m_e(j+my_id*sr%a_ee%step,:) = response_m_e(j+my_id*sr%a_ee%step,:) + sr%a_ee%loc_mat(j,:)
          end do
          call MPI_AllREDUCE(MPI_IN_PLACE,response_m_e,size(response_m_e),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

          response_m_h =0.0
          response_m_h(sr%a_ee%ind_start:sr%a_ee%ind_end,:) = sr%a_ee%loc_mat(:,:)
          call MPI_AllREDUCE(MPI_IN_PLACE,response_m_h,size(response_m_h),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

          ! --- m_j =  matmul(a_ey, m_d)
          call matrix_multiplication(my_id,sr%a_ey,mat2=response_m_d,res_mat_not_distr=response_m_j)
          call MPI_AllREDUCE(MPI_IN_PLACE,response_m_j,size(response_m_j),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
          call  dealloc_distr(response_m_a)
          call  dealloc_distr(response_m_d)
          
          deallocate( tmp_d_s )

        endif
        ! ================================================================================================

        ! ===================================== IDEAL WALL ================================================

      else ! (Ideal wall) 

        response_m_a%loc_mat(:,:) = 0.d0
        response_d_b(:)   = 0.d0
        response_d_c(:)   = 0.d0
        response_d_b_i(:) = 0.d0
        response_d_c_i(:) = 0.d0
        
        ! response_m_e
        if (halo_freeb .and. CARIDDI_mode) then 
          response_m_e_i=0.0
          response_m_e_i(sr%a_id%ind_start:sr%a_id%ind_end,:) =sr%a_id%loc_mat(:,:)
          call MPI_AllREDUCE(MPI_IN_PLACE,response_m_e_i,size(response_m_e_i),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)        
        else
          response_m_e=0.0
          response_m_e(sr%a_id%ind_start:sr%a_id%ind_end,:) =sr%a_id%loc_mat(:,:)
          call MPI_AllREDUCE(MPI_IN_PLACE,response_m_e,size(response_m_e),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        endif

          
        response_m_h =0.0
        response_m_h(sr%a_id%ind_start:sr%a_id%ind_end,:) = sr%a_id%loc_mat(:,:)
        call MPI_AllREDUCE(MPI_IN_PLACE,response_m_h,size(response_m_h),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        
        if (halo_freeb .and. CARIDDI_mode) then
          response_m_j_i(:,:) = 0.d0
        else
          response_m_j(:,:) = 0.d0
        endif

      end if

      if ( vacuum_debug ) then
        if(my_id == 0)   write(*,*) 'DEBUG: Checksums'

        test_sum=0.0
        test_sum2=0.0
        if(my_id == 0)   write(*,*) 'd_b:', sum(abs(response_d_b)), sum(response_d_b)
        if(my_id == 0)   write(*,*) '1+d_b:',sum(abs(1.d0+response_d_b)), sum(1.d0+response_d_b)
        if(my_id == 0)   write(*,*) 'd_c:', sum(abs(response_d_c)), sum(response_d_c)
        if(my_id == 0)   write(*,*) 'd_b_i:', sum(abs(response_d_b_i)), sum(response_d_b_i)
        if(my_id == 0)   write(*,*) '1+d_b_i:',sum(abs(1.d0+response_d_b_i)), sum(1.d0+response_d_b_i)
        if(my_id == 0)   write(*,*) 'd_c_i:', sum(abs(response_d_c_i)), sum(response_d_c_i)

        test_sum=0.0
        test_sum2=0.0
        if(my_id == 0)   write(*,*) 'm_e:', sum(abs(response_m_e)), sum(response_m_e)

        if(my_id == 0)   write(*,*) 'm_h:', sum(abs(response_m_h)), sum(response_m_h)
        if(my_id == 0)   write(*,*) 'm_j:', sum(abs(response_m_j)), sum(response_m_j)
    
        if(my_id == 0)   write(*,*) 'END: Checksums'
      end if

    end if  

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine update_response
  
  
  !> Read a STARWALL response matrix from a file.
  subroutine read_response_matrix( matrix, dim_expected, filename, err )
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: matrix(:,:)     !< Matrix to be read
    integer,             intent(in)    :: dim_expected(2) !< Matrix dimension expected
    character(len=*),    intent(in)    :: filename        !< Filename to read from
    integer,             intent(out)   :: err             !< Error code
    
    ! --- Local variables
    integer :: dim(2), i, j, i2, j2
    
    err = 0
    
    open(42, FILE=trim(filename), status='old', action='read', iostat=err)
    if ( err /= 0 ) return

    read(42,*) dim

    if ( (dim(1) /= dim_expected(1)) .or. (dim(2) /= dim_expected(2)) ) then
      write(*,*) 'FATAL ERROR: Matrix dimension not as expected. Different resolutions?'
      write(*,'(1x,A,2I7)') 'dim_expected=',dim_expected
      stop
    end if
    
    if ( allocated(matrix) ) deallocate(matrix)
    allocate( matrix(dim(1),dim(2)) )
    matrix = 0.d0
    
    do i = 1, dim(1)
      do j = 1, dim(2)
        
        read (42,*) i2, j2, matrix(i,j)
        
        if ( ( i2 /= i ) .or. ( j2 /= j ) ) then
          write(*,*) 'FATAL ERROR: Matrix indices not as expected. Different resolutions?'
          stop
        end if
      
      end do
    end do
      
    close(42)

  end subroutine read_response_matrix
  
  
  
  !> Read a diagonal STARWALL response matrix from a file.
  subroutine read_response_diagonal( diagonal, dim_expected, filename, err )
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: diagonal(:)     !< Matrix to be read
    integer,             intent(in)    :: dim_expected    !< Matrix dimension expected
    character(len=*),    intent(in)    :: filename        !< Filename to read from
    integer,             intent(out)   :: err             !< Error code
    
    ! --- Local variables
    integer :: dim, i, i2
    
    open(42, FILE=trim(filename), status='old', action='read', iostat=err)
    if ( err /= 0 ) return

    read(42,*) dim

    if ( dim /= dim_expected ) then
      write(*,*) 'FATAL ERROR: Matrix dimension not as expected. Different resolutions?'
      stop
    end if
    
    if ( allocated(diagonal) ) deallocate(diagonal)
    allocate( diagonal(dim) )
    diagonal = 0.d0
    
    do i = 1, dim
        
      read (42,*) i2, diagonal(i)
      
      if ( i2 /= i ) then
        write(*,*) 'FATAL ERROR: Matrix indices not as expected. Different resolutions?'
        stop
      end if
        
    end do
      
    close(42)
    
  end subroutine read_response_diagonal
  
  
  
  !> Determine the index in the response matrix for a certain boundary degree of freedom.
  integer recursive function response_index(inode, i_starwall, ibas)
    
    ! --- Routine parameters
    integer, intent(in)    :: inode      !< Boundary index of the node
    integer, intent(in)    :: i_starwall !< STARWALL harmonic
    integer, intent(in)    :: ibas       !< Basis function (1 or 2)
    
    if ( (i_starwall < 0) .or. (i_starwall > sr%n_tor) ) then
      write(*,*) 'response_index: illegal value i_starwall=', i_starwall
      stop
    end if
    response_index = 2*sr%n_tor0*(inode-1) + 2*(i_starwall-1) + ibas
    
    if ( response_index < 1 ) then
      write(*,*) 'FATAL: RESPONSE_INDEX < 1 DETECTED'
      stop
    end if
    
  end function response_index
  
  
  
  !> Determine the index in the response matrix for a certain boundary degree of freedom.
  integer recursive function response_index_eq(inode, ibas)
    
    ! --- Routine parameters
    integer, intent(in)    :: inode      !< Boundary index of the node
    integer, intent(in)    :: ibas       !< Basis function (1 or 2)
    
    response_index_eq = 2*(inode-1) + ibas
    
    if ( response_index_eq < 1 ) then
      write(*,*) 'FATAL: RESPONSE_INDEX_EQ < 1 DETECTED'
      stop
    end if
    
  end function response_index_eq
  
  
  
  !> Return a description for a given toroidal mode index (mode number, cos/sin).
  character(len=12) function mode_to_str(i_tor, n_period)
    
    use phys_module, only: mode, mode_type
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_tor     !< Toroidal mode index
    integer, intent(in) :: n_period  !< Periodicity
    
    ! --- Local variables
    integer           :: n                ! Toroidal mode number
    character(len=3)  :: typ              ! sin or cos
    character(len=30) :: i_tor_str, n_str ! Character string representations for i_tor and n
    

    write(i_tor_str,'(I10)') i_tor   ! toroidal mode index
    n   = mode(i_tor)                ! toroidal mode number
    write(n_str,'(I10)') n           !  -"-
    typ = mode_type(i_tor)           ! sin or cos
    
    mode_to_str = trim(adjustl(i_tor_str))//' (n='//trim(adjustl(n_str))//' '//trim(typ)//')'
    
  end function mode_to_str
  
  
  
  !> Return a description of several toroidal modes.
  character(len=1400) function modes_to_str(i_tors, n_tor, n_period)
  
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: n_tor         !< Dimension of i_tors
    integer, intent(in) :: i_tors(n_tor) !< Toroidal mode numbers
    integer, intent(in) :: n_period      !< Periodicity
    
    ! --- Local variables
    integer :: i
    
    modes_to_str = ''
    do i = 1, n_tor
      modes_to_str = trim(modes_to_str)//' '//mode_to_str(i_tors(i), n_period)
      if ( i == n_tor ) exit
      modes_to_str = trim(modes_to_str)//', '
    end do
    
    modes_to_str = adjustl(modes_to_str)
    
  end function modes_to_str


  subroutine memory_prediction(n_w, nd_bez, ntasks)
  
    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer, intent(in) :: n_w           !< Number of wall 
    integer, intent(in) :: nd_bez        !< Number of bezier elements
    integer, intent(in) :: ntasks        !< Number of MPI tasks

    ! --- Local variables
    real*8    :: tot_mem, gbyte_conv, real_n_w, real_nd_bez
    real*8    :: s_ww_dim
    real_n_w = real(n_w,8); real_nd_bez = real(nd_bez,8)

    if (vacuum_min) then
      s_ww_dim = real(real_n_w * sr%ncoil,8)
    else
      s_ww_dim = real(real_n_w * real_n_w,8)
    end if

    100 format(80('-'))
    101 format(5x, a,'(',i6,',',i6,') = ',E11.3E2,' GB total / ',E11.3E2, ' GB per MPI task')
    102 format(10x,'TOTAL MEMORY CONSUMPTION PREDICTION = ', E12.3E2, ' GB total /', E11.3E2, 'GB per MPI task') 
    
    !1 DP converted to GB 
    gbyte_conv = 8.d0/1024.d0/1024.d0/1024.d0
    
    write(*,100)  
    write(*,*) 'Predicted memory consumption in the JOREK-STARWALL part:'
    write(*,*)
    write(*,101)'s_ww        ', n_w,      n_w,      s_ww_dim    *               gbyte_conv, s_ww_dim    *               gbyte_conv / ntasks
    write(*,101)'s_ww_inv    ', n_w,      n_w,      real_n_w    * real_n_w    * gbyte_conv, real_n_w    * real_n_w    * gbyte_conv / ntasks
    write(*,101)'a_ye        ', n_w,      nd_bez,   real_n_w    * real_nd_bez * gbyte_conv, real_n_w    * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'a_ey        ', nd_bez,   n_w,      real_nd_bez * real_n_w    * gbyte_conv, real_nd_bez * real_n_w    * gbyte_conv / ntasks
    write(*,101)'a_ee        ', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'a_id        ', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_a', n_w,      nd_bez,   real_n_w    * real_nd_bez * gbyte_conv, real_n_w    * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_h', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_j', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_e', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks

    tot_mem =  (2.d0*real_n_w**2 + 7.d0*real_n_w*real_nd_bez + 6.d0*real_nd_bez**2) * gbyte_conv

    write(*,*)
    write(*,102) tot_mem, tot_mem / ntasks
    write(*,100)

  end subroutine memory_prediction 

  subroutine distribute_coil_names
    use vacuum

    if ((sr%n_diag_coils .gt. 0) .and. (.not. allocated(diag_coil_name)) )  then
      allocate( diag_coil_name(sr%n_diag_coils) )
      diag_coil_name(:) = sr%coil_name(sr%ind_start_diag_coils:sr%ind_start_diag_coils +  sr%n_diag_coils -1)
    end if

    if (( sr%n_rmp_coils .gt. 0) .and. (.not. allocated( rmp_coil_name)) ) then
      allocate( rmp_coil_name(sr%n_rmp_coils) )
      rmp_coil_name(:)  = sr%coil_name(  sr%ind_start_rmp_coils:sr%ind_start_rmp_coils +  sr%n_rmp_coils  -1)
    end if
    
    if (( sr%n_pol_coils .gt. 0) .and. (.not. allocated(  pf_coil_name)) )  then
      allocate( pf_coil_name(sr%n_pol_coils) )
      pf_coil_name(:)   = sr%coil_name(  sr%ind_start_pol_coils:sr%ind_start_pol_coils +  sr%n_pol_coils  -1)
    end if
    
  end subroutine distribute_coil_names

  ! subroutine row_col_distrib(my_id, ntasks, dim , old_mat, new_mat, row_wise )

  !   ! -- Converts a row-distributed matrix to a column distributed one (viceversa to be implemented)
  !   use mpi_mod
  !   implicit none
  !   ! --- Routine parameters
  !   integer, intent(in) :: my_id         !< Process ID
  !   integer, intent(in) :: ntasks        !< Number of MPI tasks
  !   integer, intent(in) :: dim(2)        !< Dimensions of matrix
  !   logical, intent(in) :: row_wise      !< How to distribute the new matrix?
  !   type(t_distrib_mat) , intent(in)    :: old_mat
  !   type(t_distrib_mat) , intent(inout) :: new_mat
  !   ! --- Local variables
  !   integer             :: step_row , step_col , loc_size
  !   integer             :: indmin_row(ntasks) , indmax_row(ntasks) , indmin_col(ntasks) , indmax_col(ntasks)
  !   integer             :: status(MPI_STATUS_SIZE) , ierr, i , j
  !   real*8, ALLOCATABLE :: buf_send(:,:), buf_recv(:,:)

  !   ! Calculate indmin_row, indmax_row, indmin_col, indmax_col (all size ntasks)

  !   !!! N.B. THIS WILL DIVIDE THE DIMENSIONS UNEVENLY FOR RECTANGULAR MATRICES

  !   step_row = dim(1) / ntasks
  !   step_col = dim(2) / ntasks

  !   ! step_row = old_mat%step
  !   ! step_col = old_mat%step

  !   do i =1,ntasks
  !     indmin_row(i) = 1 + (i-1) * step_row
  !     indmax_row(i) = i * step_row 
  !     ! The last task gets the larger chunk
  !     if (i==ntasks) indmax_row(i) = dim(1)
  !     indmin_col(i) = 1 + (i-1) * step_col 
  !     indmax_col(i) = i * step_col
  !     ! The last task gets the larger chunk
  !     if (i==ntasks) indmax_col(i) = dim(2)
  !   enddo
    
  !   ! ALLOCATE NEW_MAT MANUALLY

  !   call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
  !   if (allocated(new_mat%loc_mat)) call dealloc_distr(new_mat)

  !   new_mat%step      = step_col

  !   loc_size = indmax_col(my_id+1) - indmin_col(my_id+1) + 1
  !   new_mat%row_wise  = .false.
  !   new_mat%ind_start = indmin_col(my_id+1)
  !   new_mat%ind_end   = indmax_col(my_id+1)
  !   allocate( new_mat%loc_mat(dim(1), loc_size) )
    
  !   new_mat%loc_mat(:,:) = 0.d0
  !   write(*,*) 'size new_mat%loc_mat , my_id', size(new_mat%loc_mat,1) , size(new_mat%loc_mat,2) , my_id
  !   write(*,*) 'size old_mat%loc_mat , my_id', size(old_mat%loc_mat,1) , size(old_mat%loc_mat,2) , my_id

  !   if (my_id .eq. 0) write(*,*) 'row_col_distrib: global dimensions=', dim(1) , dim(2)
    
  !   do i = 1, ntasks
      
  !     do j = 1, ntasks

  !       if ( (i == j) .and. (i == my_id + 1) ) then
          
  !         ! We can just locally copy
  !         new_mat%loc_mat(indmin_row(i):indmax_row(i),:) = &
  !           old_mat%loc_mat(:,indmin_col(j):indmax_col(j))

  !       else if ( my_id + 1 == i) then
  !         ! We need to send

  !         ! Allocate buf_send block (buffer for row-distributed matrix)
  !         if (allocated(buf_send)) deallocate(buf_send)

  !         ! The largest dimensions at last process will be on the row for old_mat and the column for new_mat
  !         write(*,*) 'buf_send dimensions, my_id, indmax_row(i), indmin_row(i), indmax_col(j), indmin_col(j)', indmax_row(i)-indmin_row(i)+1 , indmax_col(j)-indmin_col(j)+1, my_id, indmax_row(i), indmin_row(i), indmax_col(j), indmin_col(j)
  !         ! For i=1,j=4:
  !         ! We are allocating as 7201 x 7202
  !         allocate(buf_send(indmax_row(i)-indmin_row(i)+1, indmax_col(j)-indmin_col(j)+1))
          
  !         buf_send = old_mat%loc_mat(:,indmin_col(j):indmax_col(j))
  !         call MPI_Send (buf_send, size(buf_send) , MPI_DOUBLE_PRECISION, j-1, 0 , MPI_COMM_WORLD, ierr ) 

  !       else if ( my_id + 1 == j) then

  !         ! We need to receive; note that i and j are used exchanged now!

  !         if (allocated(buf_recv)) deallocate(buf_recv)
  !         write(*,*) 'i,j,buf_recv dimensions, my_id, indmin_row(j), indmax_row(j), indmin_row(i), indmax_row(i), indmin_col(i), indmax_col(i), indmin_col(j), indmax_col(j)', i,j,indmax_row(j)-indmin_row(j)+1 , indmax_col(i)-indmin_col(i)+1, my_id, indmin_row(j), indmax_row(j), indmin_row(i), indmax_row(i), indmin_col(i), indmax_col(i), indmin_col(j), indmax_col(j)
  !         allocate(buf_recv(indmax_row(j)-indmin_row(j)+1, indmax_col(i)-indmin_col(i)+1))
  !         ! ex. i=1,j=4, 4 MPI procs
  !         ! At proc. 4 buf_recv gets allocated as 7204 x 7199
  !         ! But we sent something that is 7201 x 7202
  !         buf_recv = 0.d0

  !         call MPI_Recv( buf_recv, size(buf_recv), MPI_DOUBLE_PRECISION, i-1, 0, MPI_COMM_WORLD, status, ierr)

  !         new_mat%loc_mat(indmin_row(j):indmax_row(j),indmin_col(i):indmax_col(i)) = buf_recv

  !       end if

  !     end do
      
  !   end do

  ! end subroutine row_col_distrib

  !> Converts a row-distributed matrix to a column distributed one (or the other way around)
  subroutine row_col_distrib(my_id, ntasks, dim, old_mat, new_mat, row_wise)

    use mpi_mod
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in)                :: my_id    !< MPI rank
    integer, intent(in)                :: ntasks   !< Number of MPI ranks
    integer, intent(in)                :: dim(2)   !< Global dimensions of the matrix ! ############## not needed, should be taken from old_mat
    logical, intent(in)                :: row_wise !< How to distribute the new matrix? ! ############### not needed, we just convert into the other format, respectively
    type(t_distrib_mat), intent(in)    :: old_mat
    type(t_distrib_mat), intent(inout) :: new_mat
    
    ! --- Local variables
    integer             :: step_row, step_col, loc_size
    integer             :: indmin_row(ntasks), indmax_row(ntasks), indmin_col(ntasks), indmax_col(ntasks)
    integer             :: status(MPI_STATUS_SIZE), ierr, i ,j
    real*8, allocatable :: buf_send(:,:), buf_recv(:,:)
    
    ! --- Convert row-wise distributed to colum-wise distributed representation
    if ( old_mat%row_wise ) then
      
      ! --- Calculate distribution across matrices indmin_row, indmax_row, indmin_col, indmax_col (all size ntasks)
      step_row = dim(1) / ntasks
      step_col = dim(2) / ntasks
      do i =1, ntasks
        ! --- Matrix distribution in row-wise format
        indmin_row(i) = 1 + (i-1) * step_row
        indmax_row(i) = i * step_row
        if (i==ntasks) indmax_row(i) = dim(1) ! The last task gets the larger chunk

        ! --- Matrix distribution in colum-wise format
        indmin_col(i) = 1 + (i-1) * step_col 
        indmax_col(i) = i * step_col
        if (i==ntasks) indmax_col(i) = dim(2) ! The last task gets the larger chunk
      enddo
      
      ! --- Allocate and initialize the new matrix
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
      if (allocated(new_mat%loc_mat)) call dealloc_distr(new_mat)
      new_mat%step         = step_col
      new_mat%distrib      = .true.
      loc_size             = indmax_col(my_id+1) - indmin_col(my_id+1) + 1
      new_mat%row_wise     = .false.
      new_mat%ind_start    = indmin_col(my_id+1)
      new_mat%ind_end      = indmax_col(my_id+1)
      new_mat%dim          = (/ dim(1) , dim(2)   /)
      allocate( new_mat%loc_mat(dim(1), loc_size) )
      new_mat%loc_mat(:,:) = 0.d0

      ! write(*,*) 'size new_mat%loc_mat , my_id', size(new_mat%loc_mat,1) , size(new_mat%loc_mat,2) , my_id
      ! write(*,*) 'size old_mat%loc_mat , my_id', size(old_mat%loc_mat,1) , size(old_mat%loc_mat,2) , my_id


      ! --- Go through the matrix block-wise to convert
      do i = 1, ntasks ! i=block row index
        do j = 1, ntasks ! j=block column index
          
          ! --- This diagonal block remains on the same MPI rank, so we can just copy it locally
          if ( (i == j) .and. (i == my_id + 1) ) then
            
            new_mat%loc_mat(indmin_row(i):indmax_row(i),:) = &
              old_mat%loc_mat(:,indmin_col(j):indmax_col(j))
  
          ! --- This off-diagonal block is local on the MPI rank i-1, we send to MPI rank j-1
          else if ( my_id + 1 == i ) then

            if (allocated(buf_send)) deallocate(buf_send)
            !write(*,*) 'send: i, j, indmax_row(i)-indmin_row(i)+1 , indmax_col(j)-indmin_col(j)+1, my_id', i , j, indmax_row(i)-indmin_row(i)+1, indmax_col(j)-indmin_col(j)+1, my_id 
            allocate(buf_send(indmax_row(i)-indmin_row(i)+1, indmax_col(j)-indmin_col(j)+1))
            buf_send = old_mat%loc_mat(:,indmin_col(j):indmax_col(j))
            call MPI_Send (buf_send, size(buf_send) , MPI_DOUBLE_PRECISION, j-1, 0 , MPI_COMM_WORLD, ierr ) 
            deallocate(buf_send)

          ! --- This off-diagonal block is needed on MPI rank j-1, we receive from MPI rank i-1  
          else if ( my_id + 1 == j ) then

            if (allocated(buf_recv)) deallocate(buf_recv)
            !write(*,*) 'recv: i, j, indmax_row(i)-indmin_row(i)+1 , indmax_col(j)-indmin_col(j)+1, my_id', i , j, indmax_row(i)-indmin_row(i)+1, indmax_col(j)-indmin_col(j)+1, my_id 
            allocate(buf_recv(indmax_row(i)-indmin_row(i)+1, indmax_col(j)-indmin_col(j)+1))
            buf_recv = 0.d0
            call MPI_Recv( buf_recv, size(buf_recv), MPI_DOUBLE_PRECISION, i-1, 0, MPI_COMM_WORLD, status, ierr)
            ! if (my_id .eq. ntasks - 1) write(*,*) 'last rank: indmax_row(i)-indmin_row(i)+1, indmax_row(j)-indmin_row(j)+1', indmax_row(i)-indmin_row(i)+1, indmax_row(j)-indmin_row(j)+1
            ! write(*,*) 'i,j,buf_recv dimensions, my_id, indmin_row(j), indmax_row(j), indmin_row(i), indmax_row(i), indmin_col(i), indmax_col(i), indmin_col(j), indmax_col(j)', i,j,indmax_row(j)-indmin_row(j)+1 , indmax_col(i)-indmin_col(i)+1, my_id, indmin_row(j), indmax_row(j), indmin_row(i), indmax_row(i), indmin_col(i), indmax_col(i), indmin_col(j), indmax_col(j)          
            new_mat%loc_mat(indmin_row(i):indmax_row(i),:) = buf_recv
            deallocate(buf_recv)
  
          end if
  
        end do ! j loop
      end do ! i loop
    
    ! --- Convert column-wise to row-wise representation
    else ! (i.e., old matrix was column wise)
      
      write(*,*) 'Conversion of column wise to row wise matrix not yet implemented.'
      stop
      
    end if
  
  end subroutine row_col_distrib


end module vacuum_response

  
