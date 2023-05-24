! Module to calculate the radiation emitted by particles during motion.
! Following the algorithm given by M. Pardal, et al. in Comp. Phys. Comm., Vol.
! 285, 2023, 108634, ISSN 0010-4655, (https://doi.org/10.1016/j.cpc.2022.108634).

! Written by B. Unzicker

MODULE deck_calc_radiation_block
    
    USE strings_advanced
    USE utilities
    USE calc_radiation

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE calc_radiation_deck_initialise
#ifdef CALC_RADIATION
            IF (deck_state /= c_ds_first) RETURN

            use_calc_radiation = .FALSE.
            radiation_species = ''
            calc_radiation_start_time = 0.0_num
        
            detector_type = 'point'
            x_det_min = 1.0_num
            x_det_max = 0.0_num
            y_det_min = 0.0_num
            y_det_max = 0.0_num
            z_det_min = 0.0_num
            z_det_max = 0.0_num
            nx_det = 1
            ny_det = 1
            nz_det = 1
            detector_pos = 1.0_num
        
            t_det_min = 0.0_num
            t_det_max = 1.0_num
            dt_det = 0.001_num

            calc_rad_E_min = 0.0_num
            calc_rad_gamma_min = 1.0_num
        
#endif
    
    END SUBROUTINE calc_radiation_deck_initialise

    ! Write function to check for errors in initializing the calc_radiation block
    SUBROUTINE calc_radiation_deck_finalise

        REAL(num) :: spec_mc_sq

        IF (deck_state == c_ds_first) RETURN

        IF (use_calc_radiation) THEN
            CALL get_radiation_species_int

            spec_mc_sq= species_list(rad_species_int)%mass * c**2
            calc_rad_gamma_min = (calc_rad_E_min*ev*1.0e-6_num)/spec_mc_sq + 1

            ! Generate arrays
            det_times = arange(t_det_min, t_det_max, dt_det)
            nt_det = SIZE(det_times)
            ALLOCATE(field_at_detector(nt_det, 3))
            field_at_detector = 0.0_num

            ! Detector location
            IF (detector_type == 'plane_x') THEN
                ALLOCATE(x_det_array(1), y_det_array(ny_det), z_det_array(nz_det))
                x_det_array = (/detector_pos/)
                y_det_array = linspace(y_det_min, y_det_max, ny_det)
                z_det_array = linspace(z_det_min, z_det_max, nz_det)
            ELSE IF (detector_type == 'plane_y') THEN
                ALLOCATE(x_det_array(nx_det), y_det_array(1), z_det_array(nz_det))
                x_det_array = linspace(x_det_min, x_det_max, nx_det)
                y_det_array = (/detector_pos/)
                z_det_array = linspace(z_det_min, z_det_max, nz_det)
            ELSE IF (detector_type == 'plane_z') THEN
                ALLOCATE(x_det_array(nx_det), y_det_array(ny_det), z_det_array(1))
                x_det_array = linspace(x_det_min, x_det_max, nx_det)
                y_det_array = linspace(y_det_min, y_det_max, ny_det)
                z_det_array = (/detector_pos/)
            ELSE
                ALLOCATE(x_det_array(1), y_det_array(1), z_det_array(1))
                x_det_array = (/x_det_min/)
                y_det_array = (/y_det_min/)
                z_det_array = (/z_det_min/)
            END IF

            nx_det = SIZE(x_det_array)
            ny_det = SIZE(y_det_array)
            nz_det = SIZE(z_det_array)
        END IF

    END SUBROUTINE calc_radiation_deck_finalise


    SUBROUTINE calc_radiation_block_start
        use_calc_radiation = .TRUE.
    END SUBROUTINE calc_radiation_block_start


    SUBROUTINE calc_radiation_block_end

    END SUBROUTINE calc_radiation_block_end


    FUNCTION calc_radiation_block_handle_element(element, value) RESULT(errcode)

        CHARACTER(*), INTENT(IN) :: element, value
        INTEGER :: errcode
    
        errcode = c_err_none
        IF (deck_state == c_ds_first) RETURN
        IF (element == blank .OR. value == blank) RETURN
        
#ifdef CALC_RADIATION

        IF (str_cmp(element, 'use_calc_radiation')) THEN
            use_calc_radiation = as_logical_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'include_species')) THEN
            radiation_species = TRIM(value)
            RETURN
        END IF

        IF (str_cmp(element, 'calc_radiation_start_time') & 
                .OR. str_cmp(element, 'start_time')) THEN
            calc_radiation_start_time = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'detector_type') & 
                .OR. str_cmp(element, 'type')) THEN
            detector_type = TRIM(ADJUSTL(value))
            RETURN
        END IF 

        IF(str_cmp(element, 'x_det_pos') & 
                .OR. str_cmp(element, 'y_det_pos') &
                .OR. str_cmp(element, 'z_det_pos')) THEN
            detector_pos = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'x_det_min')) THEN
            x_det_min = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'x_det_max')) THEN
            x_det_max = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'y_det_min')) THEN
            y_det_min = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'y_det_max')) THEN
            y_det_max = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'z_det_min')) THEN
            z_det_min = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'z_det_max')) THEN
            z_det_max = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'nx_det')) THEN
            nx_det = as_integer_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'ny_det')) THEN
            ny_det = as_integer_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'nz_det')) THEN
            nz_det = as_integer_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 't_det_min')) THEN
            t_det_min = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 't_det_max')) THEN
            t_det_max = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'dt_det')) THEN
            dt_det = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'calc_radiation_energy_min') & 
                .OR. str_cmp(element, 'E_min')) THEN
            calc_rad_E_min = as_real_print(value, element, errcode)
            RETURN
        END IF

        errcode = c_err_unknown_element
#endif

    END FUNCTION calc_radiation_block_handle_element

    FUNCTION calc_radiation_block_check() RESULT(errcode)

        INTEGER :: errcode
#ifdef CALC_RADIATION
        INTEGER :: io, iu
#endif
        errcode = c_err_none

#ifdef CALC_RADIATION
    IF (radiation_species == '' .OR. rad_species_int == -1) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You must set a calc_radiation species.', &
              'Please set calc_radiation species to a valid species.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_missing_elements + c_err_terminate
    END IF
#endif

    END FUNCTION calc_radiation_block_check


END MODULE deck_calc_radiation_block