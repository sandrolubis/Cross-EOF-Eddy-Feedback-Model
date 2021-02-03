; Author: Dr. Sandro Lubis (Jan 2020)
;         Calculate one-point lag correlation of the vertically averaged 
;         zonal-mean zonal wind anomalies reconstructed from projections onto two leading EOFs
;         Lubis and Hassanzadeh, 2020 JAS (Fig. 1 and Fig. 6)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

print("Reading the input file...")

f = addfile("data/u.col.anom.erainterim.nc", "r")

u = f->u(time|:, {lat|-80:-25}, lon|0)
u@_FillValue = -999


print("Weighting latitude...")

rad = 4.0 * atan(1.0) / 180.0
wgt_u = sqrt(cos(u&lat * rad))
u = (/u * conform(u, wgt_u, 1)/)

print("Calculating EOFs...")

neof = 2
eof = eofunc_n_Wrap(u, neof, False, 0)
pc = eofunc_ts_n_Wrap(u, eof, False, 0)
u = (/eof2data_n(eof, pc, 0)/)

nrun = 21
u = runave_n_Wrap(u, nrun, 0, 0)


print("Calculating lag correlation...")

lat_ref = -30.0
maxlag = 60
lag = ispan(-maxlag, maxlag, 1)

r_u_neg = esccr(u(lat|:, time|:), u({lat|lat_ref}, time|:), maxlag)
r_u_pos = esccr(u({lat|lat_ref}, time|:), u(lat|:, time|:), maxlag)

r_u = new((/dimsizes(u&lat), 2 * maxlag + 1/), float, -999)
r_u(:, 0:maxlag - 1) = r_u_neg(:, 1:maxlag:-1)
r_u(:, maxlag:) = r_u_pos

r_u!0 = "lat"
r_u!1 = "lag"

r_u&lat = u&lat
r_u&lag = lag

t_u = new((/dimsizes(u&lat), 2 * maxlag + 1/), float, -999)
p_u = t_u

do i = 0, dimsizes(lag) - 1
	do y = 0, dimsizes(u&lat) - 1
		if (r_u(y, i).eq.1) then
			p_u(y, i) = 0.0
		else
			t_u(y, i) = r_u(y, i) * sqrt((dimsizes(u&time) - 2) / (1 - r_u(y, i) ^ 2))
			p_u(y, i) = student_t(t_u(y, i), dimsizes(u&time) - 2)
		end if
	end do
end do

r_u_all = r_u
r_u = (/where(p_u.le.(0.05), r_u, r_u@_FillValue)/)


print("Plotting...")

wks = gsn_open_wks("pdf", "one_point_lag_correlation_map")

cmap1 = read_colormap_file("BlWhRe")
cmap2 = read_colormap_file("BlAqGrWh2YeOrReVi22")
cmap3 = read_colormap_file("BlueWhiteOrangeRed")

res = True
res@gsnDraw = False
res@gsnFrame = False
res@vpWidthF = 0.6
res@vpHeightF = 0.4
res@cnFillOn = True
res@cnFillPalette = cmap3
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1.0
res@cnMaxLevelValF = 1.0
res@cnLevelSpacingF = 0.1
res@lbLabelBarOn = False
res@lbBoxLinesOn = False
res@lbLabelFontHeightF = 0.02
res@trYReverse = True
res@tmYLMode = "Explicit"
res@tmYLValues = (/max(u&lat), -40.0, -60.0, -79.5/)
res@tmYLLabels = (/"25S", "40S", "60S", "80S"/)
res@tmYLMajorOutwardLengthF = 0.0
res@tmYLMinorOutwardLengthF = 0.0
res@tmYRMajorOutwardLengthF = 0.0
res@tmYRMinorOutwardLengthF = 0.0
res@tmXBMajorOutwardLengthF = 0.0
res@tmXBMinorOutwardLengthF = 0.0
res@tmXTMajorOutwardLengthF = 0.0
res@tmXTMinorOutwardLengthF = 0.0
res@tmXMajorGrid = True
res@tmXMajorGridThicknessF = 0.25
res@tmXMajorGridLineDashPattern = 1
res@tmYMajorGrid = True
res@tmYMajorGridThicknessF = 0.25
res@tmYMajorGridLineDashPattern = 1
res@tmYLLabelFontHeightF = 0.02
res@tmXBLabelFontHeightF = 0.02
res@tiYAxisFontHeightF = 0.02
res@tiXAxisFontHeightF = 0.02
res@gsnLeftString = ""
res@gsnRightString = ""
res@gsnStringFontHeightF = 0.02

resAll = True
resAll@gsnDraw = False
resAll@gsnFrame = False
resAll@cnFillOn = False
resAll@cnLinesOn = True
resAll@cnLineLabelsOn = False
resAll@cnInfoLabelOn = False
resAll@cnLineColor = "black"
resAll@cnLineThicknessF = 1.0
resAll@gsnContourNegLineDashPattern = 1
resAll@gsnContourZeroLineThicknessF = 0.0
resAll@cnLevelSelectionMode = "ManualLevels"
resAll@cnMinLevelValF = res@cnMinLevelValF
resAll@cnMaxLevelValF = res@cnMaxLevelValF
resAll@cnLevelSpacingF = res@cnLevelSpacingF * 1.0

res@tiYAxisString = "Lat [deg]"
res@tiXAxisString = "Lag [days]"
res@gsnLeftString = "ERA-I"
res@gsnRightString = "Base Lat 30~S~0~N~S"
plot = gsn_csm_contour(wks, r_u({lat|-80:-20}, lag|:), res)
plot_all = gsn_csm_contour(wks, r_u_all({lat|-80:-20}, lag|:), resAll)
overlay(plot, plot_all)

resP = True
resP@gsnMaximize = True
resP@gsnPanelLabelBar = True
resP@lbOrientation = "vertical"
resP@lbLabelStride = 2
resP@lbBoxLinesOn = False
resP@lbLabelFontHeightF = 0.015
;resP@pmLabelBarWidthF = 0.6
;resP@pmLabelBarHeightF = 0.075
;resP@pmLabelBarOrthogonalPosF = -0.01
resP@pmLabelBarWidthF = 0.05
resP@pmLabelBarHeightF = 0.275
resP@pmLabelBarOrthogonalPosF = 0.025

gsn_panel(wks, plot, (/1, 1/), resP)

end
