! Calculate and output screening factors for the reactions in a network

program screening_factors
    
    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use probin_module, only: run_prefix, small_temp, small_dens
    use eos_type_module, only : eos_get_small_temp, eos_get_small_dens
    use burn_type_module
    use reaclib_rates, only: init_reaclib, net_screening_init
    use table_rates, only: init_tabular
    use screening_module
    use microphysics_module
    use runtime_init_module
    use network
    
    implicit none
    
    type(burn_t) :: bstate
    type(plasma_state) :: pstate
    real(rt) :: density, temperature, massfractions(nspec)
    real(rt) :: Y(nspec)
    
    character(len=256) :: params_file
    integer :: params_file_unit
    
    real(rt) :: scor, dscor_dt, dscor_dd
    integer :: i
    
    namelist /params/ density, temperature, massfractions
    
    ! runtime
    call runtime_init(.true.)
    
    ! microphysics
    call microphysics_init(small_temp=small_temp, small_dens=small_dens)
    call eos_get_small_temp(small_temp)
    call eos_get_small_dens(small_dens)
    
    ! Defaults
    density = 3.43e6
    temperature = 1e9
    ! Zero mass fractions
    massfractions = 0.0
    
    ! Get initial conditions for the burn
    call get_command_argument(1, value = params_file)

    open(newunit=params_file_unit, file=params_file, status="old", action="read")
    read(unit=params_file_unit, nml=params)
    close(unit=params_file_unit)
    
    ! Make sure user set all the mass fractions to values in the interval [0, 1]
    do i = 1, nspec
       if (massfractions(i) .lt. 0 .or. massfractions(i) .gt. 1) then
          call amrex_error('mass fraction for ' // short_spec_names(i) // ' not initialized in the interval [0,1]!')
       end if
    end do
    
    ! Echo state
    write(*,*) 'State Density (g/cm^3): ', density
    write(*,*) 'State Temperature (K): ', temperature
    do i = 1, nspec
       write(*,*) 'Mass Fraction (', short_spec_names(i), '): ', massfractions(i)
    end do
    
    ! Set conditions
    bstate % T = temperature
    bstate % rho = density
    bstate % xn(:) = massfractions(:)
    
    ! Fill plasma state
    Y(:) = bstate % xn(:) * aion_inv(:)
    call fill_plasma_state(pstate, bstate % T, bstate % rho, Y)
    
    ! Compute and print screening factors
    do i = 1, nscreen
        call screen5(pstate, i, scor, dscor_dt, dscor_dd)
        write(*,*) 'fscr(', i, '): ', scor
    end do
    
end program screening_factors
