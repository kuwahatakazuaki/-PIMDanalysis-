module calc_histogram1D
  use calc_parameter, only: r, data_beads !, data_dev
  use input_parameter, &
      only: Natom, Nbeads, TNstep, Nhist, Nfile, Nbond, &
            hist_min, hist_max, hist_margin, Usave, out_hist
  implicit none
  integer, private :: j, k, l
  real(8), private :: Dhist
  real(8), public, save, allocatable :: histogram(:,:) ! 1:Xaxis, 2:Yaxis

contains
! NEED: data_beads(Nbeads,TNstep)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_1Dhist(hist_min_ex, hist_max_ex, out_hist_ex)
    use utility
    real(8), intent(in), optional :: hist_min_ex, hist_max_ex
    character(len=128), intent(in), optional :: out_hist_ex
    real(8) :: data_max, data_min, data_ave, data_dev, data_err

    if (present(hist_min_ex)) hist_min(1) = hist_min_ex
    if (present(hist_max_ex)) hist_max(1) = hist_max_ex
    if (present(out_hist_ex)) out_hist = out_hist_ex
    allocate(histogram(Nhist,2))

    print '(" *****START calculation 1D histgram******")'

    data_max = maxval(data_beads)
    data_min = minval(data_beads)
    data_ave = sum(data_beads)/size(data_beads)

    if ( hist_min(1) == 0.0d0 .and. hist_max(1) == 0.0d0 ) then
      print *, "   Using the margin parameter"
      hist_min(1) = data_min - hist_margin
      hist_max(1) = data_max + hist_margin
    else
      print *, "   Using the hist_X_min and hist_X_max"
    end if
    Dhist = (hist_max(1) - hist_min(1)) / dble(Nhist)

    print '("    X range max =", F13.6)', hist_max(1)
    print '("    X range min =", F13.6)', hist_min(1)
    print '("    Number hist =", I8)',    Nhist
    print '("    Delta hist  =", F13.6)', Dhist

    call calc_1Dhist_sub
    call calc_deviation(data_dev, data_err)

    open(Usave, file=trim(out_hist), status='replace')
      write(Usave,'(" # ", a)') trim(out_hist)
      write(Usave,'(" # Maximum hist =", F13.6)')  data_max
      write(Usave,'(" # Minimum hist =", F13.6)')  data_min
      write(Usave,'(" # Average hist =", F13.6)')  data_ave
      write(Usave,'(" # St. deviation=", F13.6)')  data_dev
      write(Usave,'(" # St. error    =", F13.6)')  data_err
      write(Usave,'(" # X range max  =", F13.6)')  hist_max(1)
      write(Usave,'(" # X range min  =", F13.6)')  hist_min(1)
      write(Usave,'(" # Max of hist  =", F13.6, I3)') maxval(histogram(:,2)), int(maxval(histogram(:,2)))+1 !maxloc(histogram(:,2))
      write(Usave,'(" # Number hist  =", I8)'   )  Nhist
      write(Usave,'(" # Delta hist   =", F13.6)')  Dhist

      do l = 1, Nhist
        write(Usave,'(F13.6,E13.4)') histogram(l,:)
      end do
    close(Usave)
    deallocate(histogram)
    print '(a,a,/)', "    Hist data is saved in ", '"'//trim(out_hist)//'"'
    print '(a,/)', " *****END calculation 1D histgram*****"
  end subroutine calc_1Dhist
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_1Dhist_sub
    histogram(:,:) = 0.0d0
    do l = 1, Nhist
      histogram(l,1) = hist_min(1) + Dhist*dble(l)
    end do
  
    do k = 1, TNstep
      do j = 1, Nbeads
        do l = 1, Nhist
          if (data_beads(j,k) <= histogram(l,1)) then
            histogram(l,2) = histogram(l,2) + 1.0d0
            goto 100
          end if
        end do
        100 continue
      end do
    end do
    histogram(:,1) = histogram(:,1) - 0.5d0 * Dhist
    histogram(:,2) = histogram(:,2)/(dble(TNstep*Nbeads)*Dhist)
  end subroutine calc_1Dhist_sub
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module calc_histogram1D

