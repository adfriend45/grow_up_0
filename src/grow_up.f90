!======================================================================!
subroutine grow_up
!----------------------------------------------------------------------!
! Advance individual tree dimensions by dt.
! Inputs:
! Cup (kg[C] ind=1 s-1)
!----------------------------------------------------------------------!
use params
use vars
!----------------------------------------------------------------------!
implicit none
real :: dMfol_max, drb_max, dht_max, dMfro_max
real :: dMa_dht, dMa_drb
real :: dMa_max, ddp_max, dMb_ddp, dMb_drb, dMb_max
real :: fr, GMfol_max, M0, GMfro_max
!----------------------------------------------------------------------!
M0 = Ma
Ahw = pi * rb_hw ** 2
Astem = pi * rb ** 2
Asw = Astem - Ahw
!----------------------------------------------------------------------!
Chia = one + alpha_ag + alpha_ag ** 2
Chib = one + alpha_bg + alpha_bg ** 2
!----------------------------------------------------------------------!
! Litter fluxes from each compartment.
!----------------------------------------------------------------------!
LMfol = Mfol / tau_fol
LMa = Ma / tau_w
LMb = Mb / tau_w
LMfro = Mfro / tau_fro
!----------------------------------------------------------------------!
! Remove litter from each compartment
!----------------------------------------------------------------------!
!Mfol = Mfol - dt * LMfol
!Ma   = Ma   - dt * LMa
!Mb   = Mb   - dt * Lmb
!Mfro = Mfro - dt * LMfro
!----------------------------------------------------------------------!
! Rates of growth for each compartment at current size ()
! Ignore fIAA for now. maths of dht and drb together OK?
!----------------------------------------------------------------------!
GMfol_max = (lasa * Asw / sla - Mfol) / dt + LMfol
GMfol = fW_gr * fT_gr * fa_Su * dMfol_base * Asw
GMfol = min (GMfol, GMfol_max)
GMfol = max (zero, GMfol)
drb   = fW_gr * fT_gr * fa_Su * drb_base
dht   = fW_gr * fT_gr * fa_Su * dht_base
dMa_dht = rho_wood * (pi / 3.0) * (rb ** 2) * Chia
dMa_drb = 2.0 * rho_wood * (pi / 3.0) * ht * rb * Chia
GMa = dMa_dht * dht + dMa_drb * drb
ddp = fW_gr * fT_gr * fa_Su * ddp_base
!if (dp > 1.5) ddp = 0.0 !****adf
dMb_ddp = rho_wood * (pi / 3.0) * (rb ** 2) * Chib
dMb_drb = 2.0 * rho_wood * (pi / 3.0) * dp * rb * Chib
GMb = dMb_ddp * ddp + dMb_drb * drb
GMfro_max = (frasa * Asw / sfra - Mfro) / dt + LMfro
GMfro = fW_gr * fT_gr * fa_Su * dMfro_base * Asw
GMfro = min (GMfro, GMfro_max)
GMfro = max (zero, GMfro)
!----------------------------------------------------------------------!
! Rate of change of each compartment ()
!----------------------------------------------------------------------!
dMfol = GMfol - LMfol
dMa = GMa - LMa
dMb = GMb - LMb
dMfro = GMfro - LMfro
!----------------------------------------------------------------------!
! Increment compartment masses ()
!----------------------------------------------------------------------!
Mfol = Mfol + dt * dMfol
Ma = Ma + dt * dMa
!write (*,*) Ma-M0,dt*dMa
Mb = Mb + dt * dMb
Mfro = Mfro + dt * dMfro
!----------------------------------------------------------------------!
if (dMa > zero) then
  ! Basal radius and height increases due to non-infill growth ()
  fr = (dMa_drb * drb) / GMa
  drb = (fr * dMa) / dMa_drb
  rb = rb + dt * drb
  dht = ((1.0 - fr) * dMa) / dMa_dht
  ht = ht + dt * dht
else
  ! Above-ground woody mass shrinking. So, keep rb but reduce ht
  ht = Ma / (rho_wood * (pi / 3.0) * (rb ** 2) * Chia)
endif
! At this rb determine dp from Mb.
dp = Mb / (rho_wood * (pi / 3.0) * (rb ** 2) * Chib)
!write (*,*) Ma-M0,dt*dMa,fr
!write (*,*) Ma,Mb, dp
!if (it == 20) stop
!----------------------------------------------------------------------!
end subroutine grow_up
!======================================================================!
