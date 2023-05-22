MODULE calc_radiation

#ifdef CALC_RADIATION
    
    USE mpi
    USE sdf
    USE shared_data

    CONTAINS

    FUNCTION arange(start, fin, delta)
    ! Return evenly spaced values within a given interval, [start, fin) 
    ! with spacing delta.
    
        REAL(num), INTENT(IN) :: start, fin, delta
        INTEGER :: n, i
        REAL(num), DIMENSION(:), ALLOCATABLE :: arange
        
        n = FLOOR((fin - start)/delta)
        ALLOCATE(arange(n))
        
        arange = delta*(/(i, i=0, n-1, 1)/) + start
        
    END FUNCTION arange
    
    FUNCTION linspace(start, fin, n)
    ! Create an array of n evenly spaced points over the range [start, fin]
    ! This is equivalent to numpy.linspace
    
        REAL(num), INTENT(IN) :: start, fin
        INTEGER, INTENT(IN) :: n
        REAL(num), DIMENSION(n) :: linspace
        REAL(num) :: delta_x
        INTEGER :: i
        
        delta_x = (fin - start)/(n - 1)
        
        DO i = 0, n - 1
            linspace(i + 1) = start + i*delta_x
        END DO
        
        linspace(n) = fin
    END FUNCTION linspace
                        
    FUNCTION cross(x,y)
    ! Compute the cross product of two vectors, x cross y
    
        REAL(num), DIMENSION(3) :: cross, x, y
    
        cross(1) = x(2)*y(3) - x(3)*y(2)
        cross(2) = x(3)*y(1) - x(1)*y(3)
        cross(3) = x(1)*y(2) - x(2)*y(1)
    
    END FUNCTION cross
    
        
    SUBROUTINE get_radiation_species_int
        INTEGER :: ispecies

        DO ispecies = 1, n_species
            IF (species_list(ispecies)%name == radiation_species) THEN
                rad_species_int = ispecies
                RETURN
            ELSE 
                rad_species_int = -1
                RETURN
            END IF
        END DO
    END SUBROUTINE get_radiation_species_int
            
    FUNCTION field(r_part, r_det, beta, beta_dot)
    ! Compute the electric field at r_det due to a particle at r_part.
    
        REAL(num) :: field_coeff
        REAL(num), DIMENSION(3), INTENT(IN) :: r_part, r_det, beta, beta_dot
        REAL(num), DIMENSION(3) :: R_vec, n_hat, numerator, field
        REAL(num) :: R_mag, denominator
    
        field_coeff = qe/(4*pi*epsilon0*c)
    
        ! Separation vector between particle and detector
        R_vec = r_det - r_part
        R_mag = SQRT(R_vec(1)**2 + R_vec(2)**2 + R_vec(3)**2)
        n_hat = R_vec/R_mag
        
        numerator = cross(n_hat, cross(n_hat - beta, beta_dot))
        denominator = R_mag*(1 - dot_product(n_hat, beta))**3
        field = field_coeff*numerator/denominator  
        
    END FUNCTION field
        
    SUBROUTINE interp_field(t, t_prev, field)
    ! Interpolate field onto t_det_array 
        REAL(num), INTENT(IN) :: t, t_prev
        REAL(num), DIMENSION(3), INTENT(IN) :: field
        INTEGER :: n_slot, n_slot_prev, n_iter
        REAL(num) :: scale_fac, t_temp
                    
        ! Define the temporary variables
        n_slot = FLOOR((t - t_det_min)/dt_det)
        n_slot_prev = FLOOR((t_prev - t_det_min)/dt_det)
        n_iter = n_slot_prev
        t_temp = t_prev
        
        ! Interpolate field onto t_det_array
        DO WHILE (n_iter < n_slot)
            scale_fac = (det_times(n_iter + 2) - t_temp)/dt_det
            field_at_detector(n_iter + 1, :) = + &
                            field_at_detector(n_iter + 1, :) + scale_fac*field
            n_iter = n_iter + 1
            t_temp = det_times(n_iter + 1)
        END DO
        
        ! Use a different scale to prevent double counting
        scale_fac = (t - det_times(n_slot + 1))/dt_det
        field_at_detector(n_slot + 1, :) = &
                            field_at_detector(n_slot + 1, :) + scale_fac*field
    END SUBROUTINE interp_field        
        

    SUBROUTINE deallocate_calc_radiation
        ! Deallocate all variables used by calc_radiation
        IF (use_calc_radiation) THEN
            DEALLOCATE(x_det_array, y_det_array, z_det_array)
            DEALLOCATE(det_times)
            DEALLOCATE(field_at_detector)
        END IF
    END SUBROUTINE deallocate_calc_radiation

#endif
END MODULE calc_radiation