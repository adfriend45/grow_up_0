!======================================================================!
subroutine grow_up
!----------------------------------------------------------------------!
! Advance individual tree dimensions by dt.
! Inputs:
! Cup (kg[C] ind-1 s-1)
!----------------------------------------------------------------------!
use params
use vars
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
! Basal heartwood area                                              (m2)
!----------------------------------------------------------------------!
Ahw = pi * rb_hw ** 2
!----------------------------------------------------------------------!
! Basal stem area                                                   (m2)
!----------------------------------------------------------------------!
Astem = pi * rb ** 2
!----------------------------------------------------------------------!
! Basal sapwood area                                                (m2)
!----------------------------------------------------------------------!
Asw = Astem - Ahw
!----------------------------------------------------------------------!
! Crown height effect on basal radial and height growth       (fraction)
!----------------------------------------------------------------------!
fIAA = one - ht / 100.0
fIAA = max (zero, fIAA)
fIAA = min (one , fIAA)
!----------------------------------------------------------------------!
! Residual respiration. Approach from Thornley and Cannell (2000),
! Eqn. (17), with temperature effect from .                                           (kg[C] ind-1 s-1)
!----------------------------------------------------------------------!
Rresidual = fT_R * kmmax * (fSu / (Km_R + fSu)) * &
           (fNfol * Mfol + fNstm * Ma_sw + fNcro * Mb_sw + fNfro * Mfro)
!----------------------------------------------------------------------!
! Litter fluxes from each compartment                 (kg[DM] ind-1 s-1)
!----------------------------------------------------------------------!
LMfol = Mfol / tau_fol
LMa   = Ma   / tau_w
LMb   = Mb   / tau_w
LMfro = Mfro / tau_fro
!----------------------------------------------------------------------!
! Total litter flux                                   (kg[DM] ind-1 s-1)
!----------------------------------------------------------------------!
LMtot = LMfol + LMa + LMb + LMfro
!----------------------------------------------------------------------!
! Total biomass                                           (kg[DM] ind-1)
!----------------------------------------------------------------------!
Mtot = Mfol + Ma + Mb + Mfro
!----------------------------------------------------------------------!
! Sucrose litter flux                                 (kg[Su] ind-1 s-1)
!----------------------------------------------------------------------!
LSu = Su * LMtot / (Mtot + eps)
!----------------------------------------------------------------------!
! Starch litter flux                                  (kg[St] ind-1 s-1)
!----------------------------------------------------------------------!
LSt = St * LMtot / (Mtot + eps)
!----------------------------------------------------------------------!
! Rate of growth of starch pool                       (kg[St] ind-1 s-1)
!----------------------------------------------------------------------!
GSt = alpha_St * Mtot_nhw * (fa_Su * fi_St - fi_Su * fa_St)
!----------------------------------------------------------------------!
! Starch mass derivative                             (kg[St] ind-1 s-1))
!----------------------------------------------------------------------!
dSt = GSt - LSt
!----------------------------------------------------------------------!
! Integrate starch pool                                   (kg[St] ind-1)
!----------------------------------------------------------------------!
St = St + dt * dSt
!----------------------------------------------------------------------!
! Rates of growth for each compartment at current size ()
! Ignore fIAA for now. maths of dht and drb together OK?
!----------------------------------------------------------------------!
! Sapwood-area-limited possible foliage growth rate over timestep
!                                                     (kg[DM] ind-1 s-1)
!----------------------------------------------------------------------!
GMfol_max = (lasa * Asw / sla - Mfol) / dt + LMfol
!----------------------------------------------------------------------!
GMfol = fW_gr * fT_gr * fa_Su * dMfol_base * Asw
GMfol = min (GMfol, GMfol_max)
GMfol = max (zero, GMfol)
drb   = fW_gr * fT_gr * fIAA * fa_Su * drb_base
dht   = fW_gr * fT_gr * fIAA * fa_Su * dht_base
!----------------------------------------------------------------------!
! Rate of change in above-ground woody biomass with height growth
!                                                           (kg[DM] m-1)
!----------------------------------------------------------------------!
dMa_dht = rho_wood * (pi / 3.0) * (rb ** 2) * Chia
!----------------------------------------------------------------------!
! Rate of change in above-ground woody biomass with radial growth
!                                                           (kg[DM] m-1)
!----------------------------------------------------------------------!
dMa_drb = 2.0 * rho_wood * (pi / 3.0) * ht * rb * Chia
!----------------------------------------------------------------------!
GMa = dMa_dht * dht + dMa_drb * drb
ddp = fW_gr * fT_gr * fa_Su * ddp_base
!----------------------------------------------------------------------!
! Rate of change in below-ground woody biomass with depth growth
!                                                           (kg[DM] m-1)
!----------------------------------------------------------------------!
dMb_ddp = rho_wood * (pi / 3.0) * (rb ** 2) * Chib
!----------------------------------------------------------------------!
! Rate of change in below-ground woody biomass with radial growth
!                                                           (kg[DM] m-1)
!----------------------------------------------------------------------!
dMb_drb = 2.0 * rho_wood * (pi / 3.0) * dp * rb * Chib
!----------------------------------------------------------------------!
GMb = dMb_ddp * ddp + dMb_drb * drb
!----------------------------------------------------------------------!
! Sapwood-area-limited possible fine root growth rate over timestep
!                                                     (kg[DM] ind-1 s-1)
!----------------------------------------------------------------------!
GMfro_max = (frasa * Asw / sfra - Mfro) / dt + LMfro
!----------------------------------------------------------------------!
GMfro = fW_gr * fT_gr * fa_Su * dMfro_base * Asw
GMfro = min (GMfro, GMfro_max)
GMfro = max (zero, GMfro)
!----------------------------------------------------------------------!
! Rate of change of each compartment                  (kg[DM] ind-1 s-1)
!----------------------------------------------------------------------!
dMfol = GMfol - LMfol
dMa   = GMa   - LMa
dMb   = GMb   - LMb
dMfro = GMfro - LMfro
!----------------------------------------------------------------------!
! Increment compartment masses                           (kg[DM] ind-1))
!----------------------------------------------------------------------!
Mfol = Mfol + dt * dMfol
Ma   = Ma   + dt * dMa
Mb   = Mb   + dt * dMb
Mfro = Mfro + dt * dMfro
!----------------------------------------------------------------------!
! Basal radius can only increase or remain unchanged. If unchanged,
! height compensates any change in biomass.
!----------------------------------------------------------------------!
if (dMa > zero) then
  ! Basal radius and height increases due to non-infill growth ()
  !--------------------------------------------------------------------!
  ! Fraction of above-ground woody growth due to radial growth
  !                                                           (fraction)
  !--------------------------------------------------------------------!
  fr = (dMa_drb * drb) / GMa
  !--------------------------------------------------------------------!
  drb = (fr * dMa) / dMa_drb
  rb = rb + dt * drb
  dht = ((1.0 - fr) * dMa) / dMa_dht
  ht = ht + dt * dht
else
  !--------------------------------------------------------------------!
  ! Above-ground woody mass shrinking or unchanged. So, keep rb but
  ! derive ht                                                        (m)
  !--------------------------------------------------------------------!
  ht = Ma / (rho_wood * (pi / 3.0) * (rb ** 2) * Chia)
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
! Compute depth given rb and Mb                                      (m)
!----------------------------------------------------------------------!
dp = Mb / (rho_wood * (pi / 3.0) * (rb ** 2) * Chib)
!----------------------------------------------------------------------!
! Basal C balance decays towards 0                                    ()
!----------------------------------------------------------------------!
An_base = An_base - dt * An_base / tau_Anbase
!----------------------------------------------------------------------!
! Percentage loss of conductivity decays towards 0                   (%)
!----------------------------------------------------------------------!
PLC = PLC - dt * PLC / tau_PLC
!----------------------------------------------------------------------!
! Total structural growth rate                        (kg[DM] ind-1 s-1)
!----------------------------------------------------------------------!
Gtot = GMfol + GMa + GMb + GMfro
!----------------------------------------------------------------------!
! Growth respiration. From Cannell and Thornley (2000) (kg[C] ind-1 s-1)
!----------------------------------------------------------------------!
Rgrowth = CDM * Gtot * (one - YG) / YG
!----------------------------------------------------------------------!
! Respiration associated with phloem loading           (kg[C] ind-1 s-1)
! Taken form Cannell and Thornley (2000), Eqn. (11)
! N.B. need to make foliage resp. consistent.
!----------------------------------------------------------------------!
Rphloem = cphloem * (Cup - Rfol)
!----------------------------------------------------------------------!
! Total respiration. Approach from Thornley and Cannell (2000),
! Eqn. (17).                                           (kg[C] ind-1 s-1)
!----------------------------------------------------------------------!
Rtot = Rphloem + Rresidual + Rgrowth
!----------------------------------------------------------------------!
! Sucrose mass derivatve                              (kg[Su] ind-1 s-1)
!----------------------------------------------------------------------!
dSu = (Cup - Rtot - CDM * Gtot - fCStarch * GSt) / fCsucrose - LSu
!----------------------------------------------------------------------!
! Integrate sucrose pool                                  (kg[Su] ind-1)
!----------------------------------------------------------------------!
Su = Su + dt * dSu
!----------------------------------------------------------------------!
end subroutine grow_up
!======================================================================!
