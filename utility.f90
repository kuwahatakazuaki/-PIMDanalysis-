module utility
  implicit none
  real(8) :: pi = atan(1.0d0)*4.0d0
contains

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_cumulative(out_cumulative)
    use input_parameter, only: TNstep, graph_step
    use calc_parameter,  only: data_step
    integer :: i, cumu_step = 100
    real(8) :: data_dev, data_err
    character(len=128), intent(in) :: out_cumulative

    open(20, file=out_cumulative, status='replace')
    do i = 1, TNstep
      if (mod(i,graph_step*cumu_step) == 0) then
        call calc_deviation(data_dev, data_err, end_step=i)
        write(20,'(I6,3F13.6)') i, sum(data_step(1:i))/dble(i), data_dev, data_err
      end if
    end do
    close(20)
  end subroutine
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_deviation(data_dev, data_err, end_step)
    use input_parameter, only: Nbeads, TNstep
    use calc_parameter,  only: data_beads
    implicit none
    integer i, j, Nstep
    real(8) :: data_ave
    integer, intent(in), optional :: end_step
    real(8), intent(out) :: data_dev
    real(8), intent(out), optional :: data_err

    if (present(end_step)) then
      Nstep = end_step
    else
      Nstep = ubound(data_beads,2)
    end if

    data_ave = sum(data_beads)/size(data_beads)
    data_dev = 0.0d0
    do i = lbound(data_beads,2), Nstep
      do j = 1, Nbeads
        data_dev = data_dev + (data_beads(j,i) - data_ave)**2
      end do
    end do
    data_dev = data_dev / dble(Nstep*Nbeads)
    data_dev = dsqrt(data_dev)
    if (present(data_err)) data_err = data_dev / sqrt(dble(Nstep))
  end subroutine calc_deviation
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(8) function norm(x)
    implicit none
    integer i
    real(8) :: x(3)

    norm=0.0d0
    do i = 1, 3
      norm = norm + x(i)**2
    end do
    norm = dsqrt(norm)
  end function
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function outer_product(a,b) result(vec)
    real(8), intent(in) :: a(3), b(3)
    real(8) :: vec(3)

    vec(1) = a(2) * b(3) - a(3) * b(2)
    vec(2) = a(3) * b(1) - a(1) * b(3)
    vec(3) = a(1) * b(2) - a(2) * b(1)
  end function outer_product
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module utility

