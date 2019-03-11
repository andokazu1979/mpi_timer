module mpi_timer
  use mpi
  use omp_lib
  implicit none

  integer, private :: i, j, myrank, nprocs, ierr
  integer, private, parameter :: max_rec_size = 5000
  integer, private, parameter :: max_num_threads = 8
  integer, private, parameter :: iunit = 97
  real(8), private, save :: start_time(1:max_rec_size,1:max_num_threads) = 0.0_8
  real(8), private, save :: end_time(1:max_rec_size,1:max_num_threads) = 0.0_8
  real(8), private, save :: elapse(1:max_rec_size,1:max_num_threads) = 0.0_8
  real(8), private, save :: min_elapse(1:max_num_threads) = 0.0_8
  real(8), private, save :: max_elapse(1:max_num_threads) = 0.0_8
  real(8), private, save :: sum_elapse(1:max_num_threads) = 0.0_8
  character(80), private, save :: comment(1:max_rec_size) = ''
  integer(8), private, save :: num_calls(1:max_rec_size) = 0
  integer(8), private, save :: all_num_calls = 0

contains

  subroutine timer_sta(num, str)
    integer, intent(in) :: num
    character(*), intent(in) :: str
    integer id
    id = 0
    !$ id = omp_get_thread_num()
    start_time(num,id+1) = MPI_WTIME()
    !$omp master
    num_calls(num) = num_calls(num) + 1
    if(myrank == 0) then
      write(comment(num), *) str
    endif
    !$omp end master
  end subroutine timer_sta

  subroutine timer_end(num)
    integer, intent(in) :: num
    integer id
    id = 0
    !$ id = omp_get_thread_num()
    end_time(num,id+1) = MPI_WTIME()
    elapse(num,id+1) = elapse(num,id+1) + (end_time(num,id+1) - start_time(num,id+1))
  end subroutine timer_end

  subroutine timer_summary
    !$omp master
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    open(iunit, file='timer.out.summary')

    do i = 1, max_rec_size
       call MPI_REDUCE(elapse(i,1), max_elapse, max_num_threads, MPI_DOUBLE_PRECISION, &
            MPI_MAX, 0, MPI_COMM_WORLD, ierr)
       call MPI_REDUCE(elapse(i,1), min_elapse, max_num_threads, MPI_DOUBLE_PRECISION, &
            MPI_MIN, 0, MPI_COMM_WORLD, ierr)
       call MPI_REDUCE(elapse(i,1), sum_elapse, max_num_threads, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call MPI_REDUCE(num_calls(i), all_num_calls, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)

       if(myrank == 0) then
          if(comment(i) /= '') then
             write(iunit, '(i4,1x,a,1x,i8,1x)') i, comment(i), all_num_calls
             write(iunit, '(a,1x)', advance='no') 'min = '
             do j = 1, max_num_threads
               write(iunit, '(f13.3,1x)', advance='no') min_elapse(j)
             enddo
             write(iunit, '()')
             write(iunit, '(a,1x)', advance='no') 'max = '
             do j = 1, max_num_threads
               write(iunit, '(f13.3,1x)', advance='no') max_elapse(j)
             enddo
             write(iunit, '()')
             write(iunit, '(a,1x)', advance='no') 'ave = '
             do j = 1, max_num_threads
               write(iunit, '(f13.3,1x)', advance='no') sum_elapse(j) / nprocs
             enddo
             write(iunit, '()')
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
      if(comment(i) /= '') then
        write(iunit,'(i4,",",a,",",i8,",")', advance='no') i, comment(i), num_calls(i)
        do j = 1, max_num_threads
          write(iunit,'(f13.3)', advance='no') elapse(i,j)
          if(j == max_num_threads) then
            write(iunit, '()')
          else
            write(iunit,'(",")', advance='no')
          endif
        enddo
      endif
    enddo

    close(iunit)
    !$omp end master
  end subroutine timer_result

end module mpi_timer
