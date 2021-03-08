begin

print("Reading the input file...")

f1 = addfile("data/column/u.anom.1979-2013.ori.nc", "r")
f2 = addfile("data/column/duepy.anom.1979-2013.ori.nc", "r")

u = f1->u(time|:, {lat|-80:-20}, lon|0)
m = f2->duepy(time|:, {latitude|-80:-20}, lon|0)

print("Weighting...")

rad = 4.0 * atan(1.0) / 180.0
wgt = sqrt(cos(u&lat * rad))

u = (/u * conform(u, wgt, 1)/)
m = (/m * conform(m, wgt, 1)/)

print("Calculating EOFs...")

neof = 2

eof = eofunc_n_Wrap(u, neof, False, 0)

printVarSummary(eof)

print("Projecting on EOFs...")

Z = eofunc_ts_n_Wrap(u, eof, False, 0)
M = eofunc_ts_n_Wrap(m, eof, False, 0)

printVarSummary(Z)
printVarSummary(M)

print("Calculating cross correlation...")

u1 = Z(0, :)
u2 = Z(1, :)
m1 = M(0, :)
m2 = M(1, :)

maxlag = 30
lag = ispan(-maxlag, maxlag, 1)

ccr_u1_u2_pos = esccr(u1 * -1.0, u2, maxlag)
ccr_u1_u2_neg = esccr(u2, u1 * -1.0, maxlag)
ccr_u1_u2 = new(2 * maxlag + 1, float)
ccr_u1_u2(0:(maxlag - 1)) = ccr_u1_u2_neg(1:maxlag:-1)
ccr_u1_u2(maxlag:) = ccr_u1_u2_pos(0:maxlag)
ccr_u1_u2!0 = "lag"
ccr_u1_u2&lag = lag

ccr_u1_m1_pos = esccr(u1, m1, maxlag)
ccr_u1_m1_neg = esccr(m1, u1, maxlag)
ccr_u1_m1 = new(2 * maxlag + 1, float)
ccr_u1_m1(0:(maxlag - 1)) = ccr_u1_m1_neg(1:maxlag:-1)
ccr_u1_m1(maxlag:) = ccr_u1_m1_pos(0:maxlag)
ccr_u1_m1!0 = "lag"
ccr_u1_m1&lag = lag

ccr_u2_m2_pos = esccr(u2, m2, maxlag)
ccr_u2_m2_neg = esccr(m2, u2, maxlag)
ccr_u2_m2 = new(2 * maxlag + 1, float)
ccr_u2_m2(0:(maxlag - 1)) = ccr_u2_m2_neg(1:maxlag:-1)
ccr_u2_m2(maxlag:) = ccr_u2_m2_pos(0:maxlag)
ccr_u2_m2!0 = "lag"
ccr_u2_m2&lag = lag

ccr_u2_m1_pos = esccr(u2, m1, maxlag)
ccr_u2_m1_neg = esccr(m1, u2, maxlag)
ccr_u2_m1 = new(2 * maxlag + 1, float)
ccr_u2_m1(0:(maxlag - 1)) = ccr_u2_m1_neg(1:maxlag:-1)
ccr_u2_m1(maxlag:) = ccr_u2_m1_pos(0:maxlag)
ccr_u2_m1!0 = "lag"
ccr_u2_m1&lag = lag

ccr_u1_m2_pos = esccr(u1, m2, maxlag)
ccr_u1_m2_neg = esccr(m2, u1, maxlag)
ccr_u1_m2 = new(2 * maxlag + 1, float)
ccr_u1_m2(0:(maxlag - 1)) = ccr_u1_m2_neg(1:maxlag:-1)
ccr_u1_m2(maxlag:) = ccr_u1_m2_pos(0:maxlag)
ccr_u1_m2!0 = "lag"
ccr_u1_m2&lag = lag

acr_u1_pos = esacr(u1, maxlag)
acr_u1_neg = esacr(u1, maxlag)
acr_u1 = new(2 * maxlag + 1, float)
acr_u1(0:(maxlag - 1)) = acr_u1_neg(1:maxlag:-1)
acr_u1(maxlag:) = acr_u1_pos(0:maxlag)

acr_u2_pos = esacr(u2, maxlag)
acr_u2_neg = esacr(u2, maxlag)
acr_u2 = new(2 * maxlag + 1, float)
acr_u2(0:(maxlag - 1)) = acr_u2_neg(1:maxlag:-1)
acr_u2(maxlag:) = acr_u2_pos(0:maxlag)

acr_m1_pos = esacr(m1, maxlag)
acr_m1_neg = esacr(m1, maxlag)
acr_m1 = new(2 * maxlag + 1, float)
acr_m1(0:(maxlag - 1)) = acr_m1_neg(1:maxlag:-1)
acr_m1(maxlag:) = acr_m1_pos(0:maxlag)

acr_m2_pos = esacr(m2, maxlag)
acr_m2_neg = esacr(m2, maxlag)
acr_m2 = new(2 * maxlag + 1, float)
acr_m2(0:(maxlag - 1)) = acr_m2_neg(1:maxlag:-1)
acr_m2(maxlag:) = acr_m2_pos(0:maxlag)

var_u1_u2 = new(2 * maxlag + 1, float)

do i = 0, 2 * maxlag
	var_u1_u2(i) = (1.0 / (dimsizes(u1) - abs(lag(i)))) * sum(acr_u1 * acr_u2)
end do

var_u1_m1 = new(2 * maxlag + 1, float)

do i = 0, 2 * maxlag
	var_u1_m1(i) = (1.0 / (dimsizes(u1) - abs(lag(i)))) * sum(acr_u1 * acr_m1)
end do

var_u2_m2 = new(2 * maxlag + 1, float)

do i = 0, 2 * maxlag
	var_u2_m2(i) = (1.0 / (dimsizes(u1) - abs(lag(i)))) * sum(acr_u2 * acr_m2)
end do

var_u1_m2 = new(2 * maxlag + 1, float)

do i = 0, 2 * maxlag
	var_u1_m2(i) = (1.0 / (dimsizes(u1) - abs(lag(i)))) * sum(acr_u1 * acr_m2)
end do

var_u2_m1 = new(2 * maxlag + 1, float)

do i = 0, 2 * maxlag
	var_u2_m1(i) = (1.0 / (dimsizes(u1) - abs(lag(i)))) * sum(acr_u2 * acr_m1)
end do

acr = new((/2, 2 * maxlag + 1/), float)
acr(0, :) = acr_u1
acr(1, :) = acr_u2

ccr1 = ccr_u1_m1

ccr2 = ccr_u2_m2

ccr3 = ccr_u2_m1

ccr4 = ccr_u1_m2

print("Plotting...")

wks = gsn_open_wks("x11", "ccr_erainterim_full_grid")
plot = new(6, graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@vpWidthF = 0.8
res@vpHeightF = 0.5
res@trXMinF = -30.0
res@trXMaxF = 30.0
res@gsnYRefLine = 0.0
res@gsnXRefLine = 0.0
res@tmYLMode = "Manual"
res@tmYLMajorOutwardLengthF = 0.0
res@tmYLMinorOutwardLengthF = 0.0
res@tmYRMajorOutwardLengthF = 0.0
res@tmYRMinorOutwardLengthF = 0.0
res@tmXBMajorOutwardLengthF = 0.0
res@tmXBMinorOutwardLengthF = 0.0
res@tmXTMajorOutwardLengthF = 0.0
res@tmXTMinorOutwardLengthF = 0.0
res@tmYLLabelFontHeightF = 0.0325
res@tmXBLabelFontHeightF = 0.0325
res@tmXMajorGrid = True
res@tmXMajorGridThicknessF = 0.5
res@tmXMajorGridLineDashPattern = 1
res@tmYMajorGrid = True
res@tmYMajorGridThicknessF = 0.5
res@tmYMajorGridLineDashPattern = 1
res@tiXAxisString = ""
res@tiXAxisFontHeightF = 0.04
res@gsnRightString = ""
res@gsnStringFontHeightF = 0.06

res@xyLineColors = (/"red", "blue4"/)
res@xyLineThicknesses = (/5.0, 5.0/)
res@xyDashPatterns = (/0, 0/)

res@trYMinF = -0.2
res@trYMaxF = 1.0
res@tmYLTickStartF = 0.0
res@tmYLTickEndF = 1.0
res@tmYLTickSpacingF = 0.25
res@gsnCenterString = "z~B~1~N~z~B~1~N~, z~B~2~N~z~B~2~N~"
plot(0) = gsn_csm_xy(wks, lag, acr, res)

delete(res@xyLineColors)
delete(res@xyLineThicknesses)
delete(res@xyDashPatterns)

res@xyLineColor = "black"
res@xyLineThicknessF = 5.0
res@xyDashPattern = 0

res@trYMinF = -0.2
res@trYMaxF = 0.15
res@tmYLTickStartF = -0.2
res@tmYLTickEndF = 0.15
res@tmYLTickSpacingF = 0.1
res@gsnCenterString = "z~B~1~N~z~B~2"
plot(1) = gsn_csm_xy(wks, lag, ccr_u1_u2, res)

res@trYMinF = -0.3
res@trYMaxF = 0.6
res@trXMinF = -30.0
res@trXMaxF = 30.0
res@tmYLTickStartF = -0.2
res@tmYLTickEndF = 0.6
res@tmYLTickSpacingF = 0.2
res@gsnCenterString = "m~B~1~N~z~B~1"
plot(2) = gsn_csm_xy(wks, lag, ccr1, res)

res@trYMinF = -0.3
res@trYMaxF = 0.6
res@trXMinF = -30.0
res@trXMaxF = 30.0
res@tmYLTickStartF = -0.2
res@tmYLTickEndF = 0.6
res@tmYLTickSpacingF = 0.2
res@gsnCenterString = "m~B~2~N~z~B~2"
plot(3) = gsn_csm_xy(wks, lag, ccr2, res)

res@trYMinF = -0.15
res@trYMaxF = 0.15
res@trXMinF = -30.0
res@trXMaxF = 30.0
res@tmYLTickStartF = -0.2
res@tmYLTickEndF = 0.2
res@tmYLTickSpacingF = 0.1
res@tiXAxisString = "Lag [days]"
res@gsnCenterString = "m~B~1~N~z~B~2"
plot(4) = gsn_csm_xy(wks, lag, ccr3, res)

res@trYMinF = -0.15
res@trYMaxF = 0.15
res@trXMinF = -30.0
res@trXMaxF = 30.0
res@tmYLTickStartF = -0.2
res@tmYLTickEndF = 0.2
res@tmYLTickSpacingF = 0.1
res@tiXAxisString = "Lag [days]"
res@gsnCenterString = "m~B~2~N~z~B~1"
plot(5) = gsn_csm_xy(wks, lag, ccr4, res)

xp_1 = new(2 * dimsizes(lag), float)
yp_1 = new(2 * dimsizes(lag), float)

do i = 0, (dimsizes(lag) - 1)
	dy = sqrt(var_u1_u2(i))
	xp_1(i) = lag(i)
	yp_1(i) = dy
	xp_1(2 * dimsizes(lag) - 1 - i) = lag(i)
	yp_1(2 * dimsizes(lag) - 1 - i) = -dy
end do

xp_2 = new(2 * dimsizes(lag), float)
yp_2 = new(2 * dimsizes(lag), float)

do i = 0, (dimsizes(lag) - 1)
	dy = sqrt(var_u1_m1(i))
	xp_2(i) = lag(i)
	yp_2(i) = dy
	xp_2(2 * dimsizes(lag) - 1 - i) = lag(i)
	yp_2(2 * dimsizes(lag) - 1 - i) = -dy
end do

xp_3 = new(2 * dimsizes(lag), float)
yp_3 = new(2 * dimsizes(lag), float)

do i = 0, (dimsizes(lag) - 1)
	dy = sqrt(var_u2_m2(i))
	xp_3(i) = lag(i)
	yp_3(i) = dy
	xp_3(2 * dimsizes(lag) - 1 - i) = lag(i)
	yp_3(2 * dimsizes(lag) - 1 - i) = -dy
end do

xp_4 = new(2 * dimsizes(lag), float)
yp_4 = new(2 * dimsizes(lag), float)

do i = 0, (dimsizes(lag) - 1)
	dy = sqrt(var_u2_m1(i))
	xp_4(i) = lag(i)
	yp_4(i) = dy
	xp_4(2 * dimsizes(lag) - 1 - i) = lag(i)
	yp_4(2 * dimsizes(lag) - 1 - i) = -dy
end do

xp_5 = new(2 * dimsizes(lag), float)
yp_5 = new(2 * dimsizes(lag), float)

do i = 0, (dimsizes(lag) - 1)
	dy = sqrt(var_u1_m2(i))
	xp_5(i) = lag(i)
	yp_5(i) = dy
	xp_5(2 * dimsizes(lag) - 1 - i) = lag(i)
	yp_5(2 * dimsizes(lag) - 1 - i) = -dy
end do

resPolygon = True
resPolygon@tfPolyDrawOrder = "PreDraw"
resPolygon@gsFillColor = "gray80"

polygon_1 = gsn_add_polygon(wks, plot(1), xp_1, yp_1, resPolygon)
polygon_2 = gsn_add_polygon(wks, plot(2), xp_2, yp_2, resPolygon)
polygon_3 = gsn_add_polygon(wks, plot(3), xp_3, yp_3, resPolygon)
polygon_4 = gsn_add_polygon(wks, plot(4), xp_4, yp_4, resPolygon)
polygon_5 = gsn_add_polygon(wks, plot(5), xp_5, yp_5, resPolygon)

resText = True
resText@txFontColor = "black"
resText@txFontHeightF = 0.045

text_1 = gsn_add_text(wks, plot(1), "z~B~1~N~ leads", 15.0, -0.06, resText)
text_2 = gsn_add_text(wks, plot(2), "z~B~1~N~ leads", 15.0, -0.13, resText)
text_3 = gsn_add_text(wks, plot(3), "z~B~2~N~ leads", 15.0, -0.13, resText)
text_4 = gsn_add_text(wks, plot(4), "z~B~2~N~ leads", 15.0, -0.075, resText)
text_5 = gsn_add_text(wks, plot(5), "z~B~1~N~ leads", 15.0, -0.075, resText)

resP = True
resP@gsnMaximize = True
resP@gsnPanelYWhiteSpacePercent = 5.0

gsn_panel(wks, plot, (/3, 2/), resP)

end
