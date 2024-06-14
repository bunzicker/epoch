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
            n_det_hor = 1
            n_det_ver = 1
            det_center = (/0.0_num, 0.0_num, 0.0_num/)
            det_angle = 0.0_num
        
            t_det_min = 0.0_num
            t_det_max = 1.0_num
            dt_det = 0.001_num

            calc_rad_E_min = 0.0_num
            calc_rad_gamma_min = 1.0_num
        
#endif
    
    END SUBROUTINE calc_radiation_deck_initialise

    SUBROUTINE calc_radiation_deck_finalise

        REAL(num) :: spec_mc_sq
        CHARACTER(LEN=8) :: string
        INTEGER :: arr_con

        IF (deck_state == c_ds_first) RETURN

        IF (use_calc_radiation) THEN
            CALL get_radiation_species_int
            IF (rad_species_int >= 1 .AND. rad_species_int <= n_species) THEN
                IF (rank == 0) THEN
                    CALL integer_as_string(rad_species_int, string)
                    PRINT *
                    PRINT *, 'Using Species ', TRIM(ADJUSTL(string)), &
                             ' (', TRIM(species_list(rad_species_int)%name) , &
                             ') to calculate radiation.' 
                END IF
            END IF

            spec_mc_sq= species_list(rad_species_int)%mass * c**2
            calc_rad_gamma_min = (calc_rad_E_min*ev*1.0e6_num)/spec_mc_sq + 1

            IF (rank == 0) THEN
                PRINT *, 'gamma_min for calc_radiation = ', calc_rad_gamma_min
                PRINT *
            END IF

            ! Generate arrays
            det_times = arange(t_det_min, t_det_max, dt_det)
            nt_det = SIZE(det_times)
            ALLOCATE(field_at_detector(nt_det, n_det_hor, n_det_ver, 3))
            field_at_detector = 0.0_num

            ! Detector location
            IF (detector_type == 'plane') THEN
                ALLOCATE(det_dir1(n_det_hor), det_dir2(n_det_hor), &
                            det_dir3(n_det_ver))
                ! det_dir1 and det_dir2 lie in the xy plane. 
                ! det_dir3 is in the z_hat direction.
                det_dir1 = linspace(-det_width/2.0_num, & 
                                det_width/2.0_num, n_det_hor) + det_center(1)
                det_dir2 = (/(det_center(2), arr_con = 1, n_det_hor)/)
                det_dir3 = linspace(-det_height/2.0_num, & 
                                det_height/2.0_num, n_det_ver) + det_center(3)
            ELSE ! use point by default
                ALLOCATE(det_dir1(1), det_dir2(1), det_dir3(1))
                det_dir1 = (/det_center(1)/)
                det_dir2 = (/det_center(2)/)
                det_dir3 = (/det_center(3)/)
            END IF

            n_det_hor = SIZE(det_dir1)
            n_det_ver = SIZE(det_dir3)

            ! Convert det_angle to radians
            det_angle = det_angle*pi/180.0_num

            ! Copy data to field_at_detector_output for dumping
            IF (rank == 0) THEN
                ALLOCATE(field_at_detector_output(nt_det, n_det_hor, &
                            n_det_ver, 3))
                field_at_detector_output = 0.0_num
            END IF
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

        IF (str_cmp(element, 'det_center')) THEN
            CALL get_vector(value, det_center, errcode)
            RETURN
          END IF
      
        IF (str_cmp(element, 'det_width')) THEN
            det_width = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'det_height')) THEN
            det_height = as_real_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'n_pixels_hor')) THEN
            n_det_hor = as_integer_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'n_pixels_ver')) THEN
            n_det_ver = as_integer_print(value, element, errcode)
            RETURN
        END IF

        IF (str_cmp(element, 'det_angle')) THEN
            det_angle = as_real_print(value, element, errcode)
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

        IF (use_calc_radiation) THEN 
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
                errcode = c_err_bad_value + c_err_terminate
            END IF
        END IF
    
    END FUNCTION calc_radiation_block_check


END MODULE deck_calc_radiation_block