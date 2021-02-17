begin

print("Reading the input file...")

f1 = addfile("data/u.col.anom.erainterim.nc", "r")
f2 = addfile("data/duepy.col.anom.erainterim.nc", "r")

u = f1->u(time|:, {lat|-80:-20}, lon|0)
duepy = f2->duepy(time|:, {latitude|-80:-20}, lon|0)

u@_FillValue = -999
duepy@_FillValue = -999

print("Weighting latitude...")

rad = 4.0 * atan(1.0) / 180.0
wgt_u = sqrt(cos(u&lat * rad))
wgt_duepy = sqrt(cos(duepy&latitude * rad))

u = (/u * conform(u, wgt_u, 1)/)
duepy = (/duepy * conform(duepy, wgt_duepy, 1)/)

print("Calculating EOFs...")

neof = 2

eof = eofunc_n_Wrap(u, neof, False, 0)
pc = eofunc_ts_n_Wrap(u, eof, False, 0)

eof(0, :) = (/-1.0 * eof(0, :)/)
pc(0, :) = (/-1.0 * pc(0, :)/)

print("Projecting...")

m_total = eofunc_ts_n_Wrap(duepy, eof, False, 0)

delete(duepy)

m_total_1 = m_total(0, :)

z1 = pc(0, :)

print("Calculating spectra...")

nSampTot = dimsizes(m_total_1)
nSampWin = 256
nSampOverlap = 128

nWindow = (nSampTot - nSampWin) / (nSampWin - nSampOverlap) + 1

cospc_total = new(nSampWin / 2, float, -999)

quspc_total = new(nSampWin / 2, float, -999)

coher_total = new(nSampWin / 2, float, -999)

phase_total = new(nSampWin / 2, float, -999)

cospc_total = 0.0

quspc_total = 0.0

coher_total = 0.0

phase_total = 0.0

spectra_z1 = new(nSampWin / 2, float, -999)

spectra_z1 = 0.0

ntStrt = 0
ntLast = nSampWin - 1

do nw = 0, (nWindow - 1)
	work_total = dtrend(m_total_1(ntStrt:ntLast), False)
	work_z1 = dtrend(z1(ntStrt:ntLast), False)
	
	work_total = wgt_runave(work_total, (/0.25, 0.5, 0.25/), 1)
	work_z1 = wgt_runave(work_z1, (/0.25, 0.5, 0.25/), 1)
	
	specxy_total = specxy_anal(work_z1, work_total, 0, 0, 0.0)
	specxy_total_smooth = specxy_anal(work_z1, work_total, 0, 3, 0.0)
	
	cospc_total = cospc_total + specxy_total@cospc
	
	quspc_total = quspc_total + specxy_total@quspc
	
	coher_total = coher_total + specxy_total_smooth@coher
	
	phase_total = phase_total + specxy_total@phase
	
	spectra_z1 = spectra_z1 + specxy_total@spcx
	
	ntStrt = ntLast - nSampOverlap + 1
	ntLast = ntStrt + nSampWin - 1
end do

cospc_total = cospc_total / nWindow

quspc_total = quspc_total / nWindow

coher_total = coher_total / nWindow

phase_total = abs(phase_total / nWindow)

spectra_z1 = spectra_z1 / nWindow

cospc_total = abs(cospc_total / spectra_z1)

quspc_total = abs(quspc_total / spectra_z1)

freq = specxy_total@frq

pi = acos(-1.0)

omega = 2.0 * pi * freq

freq_max = 0.05

omega!0 = "freq"
omega&freq = freq

omega_lowfreq = omega({freq|0:freq_max})

cospc_total!0 = "freq"
quspc_total!0 = "freq"

cospc_total&freq = freq
quspc_total&freq = freq

cospc_total_lowfreq = cospc_total({freq|0:freq_max})
quspc_total_lowfreq = quspc_total({freq|0:freq_max})

rc_real = regCoef(omega_lowfreq, cospc_total_lowfreq)
rc_imag = regCoef(omega_lowfreq, quspc_total_lowfreq)

alpha = rc_real@yintercept
beta = rc_imag

tau = beta / alpha

spectra_total = new((/3, nSampWin / 2/), float, -999)
spectra_total(0, :) = cospc_total
spectra_total(1, :) = quspc_total
spectra_total(2, :) = omega

rad = 180.0 / pi

phase_diff = new((/2, nSampWin / 2/), float, -999)
phase_diff(0, :) = phase_total
phase_diff(1, :) = atan(tau * omega) * rad

print("Plotting...")

wks = gsn_open_wks("pdf", "cross_spectrum_LH01")
plot = new(3, graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@vpWidthF = 1.0
res@vpHeightF = 0.6
res@trXMinF = 0.0
res@trXMaxF = 0.25
res@tmYLMajorOutwardLengthF = 0.0
res@tmYLMinorOutwardLengthF = 0.0
res@tmXBMajorOutwardLengthF = 0.0
res@tmXBMinorOutwardLengthF = 0.0
res@tmYLLabelFontHeightF = 0.04
res@tmXBLabelFontHeightF = 0.04
res@tiXAxisFontHeightF = 0.04
res@gsnStringFontHeightF = 0.04

res@pmLegendDisplayMode = "Always"
res@pmLegendSide = "top"
res@pmLegendWidthF = 0.2
res@pmLegendHeightF = 0.2
res@pmLegendOrthogonalPosF = -0.5
res@pmLegendParallelPosF = 0.25
res@lgPerimOn = False
res@lgItemOrder = (/2, 0, 1/)
res@lgLabelFontHeightF = 0.03
res@xyExplicitLegendLabels = (/"   Real", "   Imag", "   2*pi*freq"/)

res@trYMinF = 0.0
res@trYMaxF = 2.0
res@xyDashPatterns = (/0, 1, 16/)
res@xyLineColors = (/"black", "black", "black"/)
res@xyLineThicknesses = (/5.0, 3.0, 3.0/)

res@gsnCenterString = "(a) (Z*M)/(Z*Z)"
plot(0) = gsn_csm_xy(wks, freq, spectra_total, res)

delete(res@xyDashPatterns)
delete(res@xyLineColors)
delete(res@xyLineThicknesses)

res@trYMinF = 0.0
res@trYMaxF = 1.0
res@xyDashPattern = 0
res@xyLineColor = "black"
res@xyLineThicknessF = 5.0
res@pmLegendDisplayMode = "NoCreate"

res@gsnCenterString = "(b) Coherence Squared"
plot(1) = gsn_csm_xy(wks, freq, coher_total, res)

res@trYMinF = 0.0
res@trYMaxF = 100.0
res@xyDashPatterns = (/0, 16/)
res@xyLineColors = (/"black", "black"/)
res@xyLineThicknesses = (/5.0, 3.0/)
res@tiXAxisString = "Frequency [cpd]"

delete(res@lgItemOrder)
delete(res@xyExplicitLegendLabels)

res@pmLegendDisplayMode = "Always"
res@pmLegendWidthF = 0.2
res@pmLegendHeightF = 0.15
res@pmLegendOrthogonalPosF = -0.9
res@pmLegendParallelPosF = 0.75
res@lgItemOrder = (/1, 0/)
res@lgLabelFontHeightF = 0.03
res@xyExplicitLegendLabels = (/"   Observed", "   arctan(~F33~t*w~F21~)"/)

res@gsnCenterString = "(c) Phase Difference"
plot(2) = gsn_csm_xy(wks, freq, phase_diff, res)

resP = True
resP@gsnMaximize = True
resP@gsnPanelBottom = 0.05
resP@gsnPanelYWhiteSpacePercent = 5.0

gsn_panel(wks, plot, (/3, 1/), resP)

end
