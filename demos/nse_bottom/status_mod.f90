module status_mod
implicit none
logical, private, save :: initialized = .FALSE.
contains

subroutine get_status(my_status)
implicit none
logical, intent(out) :: my_status
my_status = initialized
end subroutine get_status

subroutine set_status(my_status)
implicit none
logical, intent(in) :: my_status
initialized = my_status
end subroutine set_status

end module status_mod
