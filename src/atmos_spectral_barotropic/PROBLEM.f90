! ---------------------------------- ADVECTION TERMS ---------------------------------------
Dyn%Grid%vq(:,:,future) = Dyn%grid%pv*Dyn%Grid%v(:,:,current) ! Full advection terms
Dyn%Grid%uq(:,:,future) = Dyn%grid%pv*Dyn%Grid%u(:,:,current)! Note the positive sign
!print *, "VQ", Dyn%Grid%vq(10,10,future)
!print *, "-UQ", -Dyn%Grid%uq(10,10,future)

! WOULD LIKE TO BREAK INTO RELATIVE AND ABSOLUTE COMPONENTS
! First compute just RELATIVE vorticity advection
! FIRST WAY (check consistency)
 call vor_div_from_uv_grid (Dyn%grid%vor(:,:,current)*Dyn%Grid%v(:,:,current), -Dyn%grid%vor(:,:,current)*Dyn%Grid%u(:,:,current), dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, rel_advec)
print *, vors_rub(10,10)
do j = js, je
 A(:,j) = zerog(:,j)+coriolis(j)
end do
 call vor_div_from_uv_grid (A*Dyn%Grid%v(:,:,current), -A*Dyn%Grid%u(:,:,current), dt_vors_rub, dt_divs_rub) 
 call trans_spherical_to_grid(dt_vors_rub, plan_advec)
! SECOND WAY
print *, vors_rub(10,10)
vors_rub = zeros
print *, vors_rub(10,10)
 call trans_grid_to_spherical  (Dyn%grid%vor(:,:,current), vors_rub)
 call horizontal_advection     (vors_rub, Dyn%Grid%u(:,:,current), Dyn%Grid%v(:,:,current), rel_advecII)
print *, vors_rub(10,10)
!print *, "RELATIVE VORTICITY", Dyn%grid%vor(10,10,current)
 call trans_grid_to_spherical  (A, vors_rub)
 call horizontal_advection     (vors_rub, Dyn%Grid%u(:,:,current), Dyn%Grid%v(:,:,current), plan_advecII)

print *, "REL_ADVEC", rel_advec(10,10)
print *, "REL_ADVECII", rel_advecII(10,10)

print *, "PLAN_ADVEC", plan_advec(10,10)
print *, "PLAN_ADVECII", plan_advecII(10,10)
