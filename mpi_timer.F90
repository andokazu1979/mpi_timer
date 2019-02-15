module mpi_timer
  use mpi
  implicit none

  integer, private :: i, myrank, nprocs, ierr
  integer, private, parameter :: max_rec_size = 5000
  integer, private, parameter :: iunit = 97
  real(8), private, save :: start_time(1:max_rec_size) = 0.0_8
  real(8), private, save :: end_time(1:max_rec_size) = 0.0_8
  real(8), private, save :: elapse(1:max_rec_size) = 0.0_8
  real(8), private, save :: min_elapse = 0.0_8
  real(8), private, save :: max_elapse = 0.0_8
  real(8), private, save :: sum_elapse = 0.0_8
  character(80), private, save :: comment(1:max_rec_size)
  integer(8), private, save :: num_calls(1:max_rec_size) = 0
  integer(8), private, save :: all_num_calls = 0

contains

  subroutine timer_sta(num, str)
    integer, intent(in) :: num
    character(*), intent(in) :: str
    !$omp master
    start_time(num) = MPI_WTIME()
    num_calls(num) = num_calls(num) + 1
    if(myrank == 0) then
      write(comment(num), *) str
    endif
    !$omp end master
  end subroutine timer_sta

  subroutine timer_end(num)
    integer, intent(in) :: num
    !$omp master
    end_time(num) = MPI_WTIME()
    elapse(num) = elapse(num) + (end_time(num) - start_time(num))
    !$omp end master
  end subroutine timer_end

  subroutine timer_summary
    !$omp master
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    open(iunit, file='timer.out.summary')

    do i = 1, max_rec_size
       call MPI_REDUCE(elapse(i), max_elapse, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, 0, MPI_COMM_WORLD, ierr)
       call MPI_REDUCE(elapse(i), min_elapse, 1, MPI_DOUBLE_PRECISION, &
            MPI_MIN, 0, MPI_COMM_WORLD, ierr)
       call MPI_REDUCE(elapse(i), sum_elapse, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call MPI_REDUCE(num_calls(i), all_num_calls, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)

       if(myrank == 0) then
          if(max_elapse /= 0.0_8 .and. min_elapse /= 0.0_8) then
             write(iunit, '(i4,1x,a,1x,a,f13.3,1x,a,1x,f13.3,1x,a,1x,f13.3,1x,a,1x,i8)') &
              i, comment(i), 'min = ', min_elapse, &
              'max =', max_elapse, 'ave =', sum_elapse / nprocs, &
              'calls(all procs) =', &
              all_num_calls
          end if
       end if
    end do

    close(iunit)
    !$omp end master
  end subroutine timer_summary

  subroutine timer_result
    character(100) :: fname

    !$omp master
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    write(fname,'(''timer.out.'',I6.6)') myrank
    open(iunit, file=fname)

    do i = 1, max_rec_size
      if(elapse(i) /= 0.0_8) then
        write(iunit,'(" ",i4,1x,a,1x,f13.3,1x,i8)') &
          i, comment(i), elapse(i), num_calls(i)
      endif
    enddo

    close(iunit)
    !$omp end master
  end subroutine timer_result

end module mpi_timer
