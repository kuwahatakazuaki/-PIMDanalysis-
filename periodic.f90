module mod_periodic
  use input_parameter
  use calc_parameter, only: data_beads, data_step
  use calc_histogram1D
  use utility
  implicit none
  private
  integer :: i, j, k, l, Ihist
  real(8) :: Dhist
  public periodic

contains

  subroutine periodic

    select case(jobtype)
      case(81)
        call RDF1
      case(82)
        call RDF2
    end select
  end subroutine periodic

! ++++++++++++++++++
! +++ Start RDF1 +++
! ++++++++++++++++++
  subroutine RDF1
    integer :: Nelement
    real(8) :: r12(3), minbox, d12, rho
    allocate(histogram(Nhist,2))
    histogram(:,:) = 0.0d0

    if ( trim(out_hist) == "0" ) then
      write(out_hist, '(a,I0,a,I0,a)') "rdf1_", Ielement1, "-", Felement1, ".out"
    end if

    minbox = minval(Lbox(:))
    Nelement = Felement1 - Ielement1 + 1
    rho = dble(Nelement) / (Lbox(1)*Lbox(2)*Lbox(3))
    Dhist = minbox / dble(Nhist)
    histogram(:,:) = 0.0d0
    do Ihist = 1, Nhist
      histogram(Ihist,1) = Dhist * dble(Ihist)  ! not dble(Ihist-1)
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        do k = Ielement1, Felement1
          do l = k+1, Felement1
            r12(:) = r(:,k,j,i) - r(:,l,j,i)
            r12(:) = r12(:) - Lbox(:) * nint(r12(:)/Lbox(:))
            d12 = dsqrt( sum( r12(:)*r12(:) ) )
            do Ihist = 1, Nhist
              if ( d12 <= histogram(Ihist,1) ) then
                histogram(Ihist,2) = histogram(Ihist,2) + 1.0d0
                goto 100
              end if
            end do
            100 continue
          end do
        end do
      end do
    end do
    histogram(:,1) = histogram(:,1) - 0.5d0 * Dhist
    histogram(:,2) = histogram(:,2) / (4*pi*rho*Dhist*TNstep*Nbeads)
    do Ihist = 1, Nhist
      histogram(Ihist,2) = histogram(Ihist,2) / (histogram(Ihist,1)*histogram(Ihist,1))
    end do

    open(Usave, file=trim(out_hist), status='replace')
      do Ihist = 1, Nhist
        write(Usave,'(F13.6, E13.4)') histogram(Ihist,:)
      end do
    close(Usave)

  end subroutine RDF1
! ++++++++++++++++
! +++ End RDF1 +++
! ++++++++++++++++

! ++++++++++++++++++
! +++ Start RDF2 +++
! ++++++++++++++++++
  subroutine RDF2
    integer :: Nelement
    real(8) :: r12(3), minbox, d12, rho
    allocate(histogram(Nhist,2))

    if ( trim(out_hist) == "0" ) then
      write(out_hist, '(a,I0,a,I0,a,I0,a,I0,a)') & 
         "rdf2_", Ielement1, "-", Felement1, "_",Ielement2, "-",Felement2, ".out"
    end if

    minbox = minval(Lbox(:))
    Nelement = (Felement1 - Ielement1 + 1) + (Felement2 - Ielement2 + 1)
    rho = dble(Nelement) / (Lbox(1)*Lbox(2)*Lbox(3))
    Dhist = minbox / dble(Nhist)
    histogram(:,:) = 0.0d0
    do Ihist = 1, Nhist
      histogram(Ihist,1) = Dhist * dble(Ihist)  ! not dble(Ihist-1)
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        do k = Ielement1, Felement1
          do l = Ielement2, Felement2
            r12(:) = r(:,k,j,i) - r(:,l,j,i)
            r12(:) = r12(:) - Lbox(:) * nint(r12(:)/Lbox(:))
            d12 = dsqrt( sum( r12(:)*r12(:) ) )
            do Ihist = 1, Nhist
              if ( d12 <= histogram(Ihist,1) ) then
                histogram(Ihist,2) = histogram(Ihist,2) + 1.0d0
                goto 100
              end if
            end do
            100 continue
          end do
        end do
      end do
    end do
    histogram(:,1) = histogram(:,1) - 0.5d0 * Dhist
    histogram(:,2) = histogram(:,2) / (4*pi*rho*Dhist*TNstep*Nbeads)
    do Ihist = 1, Nhist
      histogram(Ihist,2) = histogram(Ihist,2) / (histogram(Ihist,1)*histogram(Ihist,1))
    end do

    open(Usave, file=trim(out_hist), status='replace')
      do Ihist = 1, Nhist
        write(Usave,'(F13.6, E13.4)') histogram(Ihist,:)
      end do
    close(Usave)

  end subroutine RDF2
! ++++++++++++++++
! +++ End RDF2 +++
! ++++++++++++++++
end module mod_periodic




