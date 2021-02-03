; Author: Dr. Sandro Lubis (Jan 2020)
;         An eddy-zonal flow feedback model that accounts for 
;         the coupling of the leading two EOFs and their feedbacks
;         Lubis and Hassanzadeh, 2020 JAS (Section 3; Eqs. 4-10)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


begin

print("Reading the input file...")

f1 = addfile("data/u.col.anom.erainterim.nc", "r")
f2 = addfile("data/duepy.col.anom.erainterim.nc", "r")

u = f1->u(time|:, {lat|-80:-20}, lon|0)
u@_FillValue = -999
duepy = f2->duepy(time|:, {latitude|-80:-20}, lon|0)
duepy@_FillValue = -999


print("Weighting latitude...")


rad = 4.0 * atan(1.0) / 180.0
wgt_u = sqrt(cos(u&lat * rad))
u = (/u * conform(u, wgt_u, 1)/)
wgt_duepy = sqrt(cos(duepy&latitude * rad))
duepy = (/duepy * conform(duepy, wgt_duepy, 1)/)


print("Calculating EOFs...")

neof = 2
eof = eofunc_n_Wrap(u, neof, False, 0)
pc = eofunc_ts_n_Wrap(u, eof, False, 0)


print("Projection...")

m = eofunc_ts_n_Wrap(duepy, eof, False, 0)

z1 = pc(0, :)
z2 = pc(1, :)

m1 = m(0, :)
m2 = m(1, :)


;;;;;;;;;;;;;;;;;;;;;; Calculating cross-EOF Eddy Feedbacks ;;;;;;;;;;;;;;;;;;;;;; 

maxlag = 20

b1 = new((/maxlag + 1, 2/), float, -999)
b2 = new((/maxlag + 1, 2/), float, -999)

do i = 0, maxlag
	reg_z1_z1 = regCoef_n(z1(0:(dimsizes(z1) - 1 - i)), z1(i:(dimsizes(z1) - 1)), 0, 0)
	reg_z2_z1 = regCoef_n(z1(0:(dimsizes(z1) - 1 - i)), z2(i:(dimsizes(z2) - 1)), 0, 0)
	reg_z1_z2 = regCoef_n(z2(0:(dimsizes(z2) - 1 - i)), z1(i:(dimsizes(z1) - 1)), 0, 0)
	reg_z2_z2 = regCoef_n(z2(0:(dimsizes(z2) - 1 - i)), z2(i:(dimsizes(z2) - 1)), 0, 0)
	
	reg_m1_z1 = regCoef_n(z1(0:(dimsizes(z1) - 1 - i)), m1(i:(dimsizes(m1) - 1)), 0, 0)
	reg_m2_z1 = regCoef_n(z1(0:(dimsizes(z1) - 1 - i)), m2(i:(dimsizes(m2) - 1)), 0, 0)
	reg_m1_z2 = regCoef_n(z2(0:(dimsizes(z2) - 1 - i)), m1(i:(dimsizes(m1) - 1)), 0, 0)
	reg_m2_z2 = regCoef_n(z2(0:(dimsizes(z2) - 1 - i)), m2(i:(dimsizes(m2) - 1)), 0, 0)
	
	A_b1 = (/	(/reg_z1_z1, reg_z2_z1/),	\
				(/reg_z1_z2, reg_z2_z2/)	/)
	B_b1 = 		(/reg_m1_z1, reg_m1_z2/)
	
	b1(i, :) = solve_linsys(A_b1, B_b1)
	
	A_b2 = (/	(/reg_z2_z1, reg_z1_z1/),	\
				(/reg_z2_z2, reg_z1_z2/)	/)
	B_b2 = 		(/reg_m2_z1, reg_m2_z2/)
	
	b2(i, :) = solve_linsys(A_b2, B_b2)
end do

b11 = b1(:, 0)
b12 = b1(:, 1)
b21 = b2(:, 0) * -1.0
b22 = b2(:, 1) * -1.0

b11!0 = "lag"
b11&lag = ispan(0, maxlag, 1)

copy_VarCoords(b11, b12)
copy_VarCoords(b11, b21)
copy_VarCoords(b11, b22)

lagi = 7
lagf = 15

b11_avg = avg(b11({lag|lagi:lagf}))
b12_avg = avg(b12({lag|lagi:lagf}))
b21_avg = avg(b21({lag|lagi:lagf}))
b22_avg = avg(b22({lag|lagi:lagf}))


print("Feedback strengths averaging over the Lag " + lagi + "-" + lagf+": ")

print("b11: " + b11_avg)
print("b12: " + b12_avg)
print("b21: " + b21_avg)
print("b22: " + b22_avg)


; Periodicity and Decay Time Scale from Theory (See Section 4, Lubis and Hassanzadeh, 2020)

tau1 = 8.3
tau2 = 8.4

sigma = 0.5 * (1.0 / tau1 + 1.0 / tau2 - b11_avg - b22_avg)
sigma = 1.0 / sigma

omega = 0.5 * sqrt( abs( ( (1.0 / tau1 - 1.0 / tau2) - (b11_avg - b22_avg) ) ^ 2.0 + 4.0 * b12_avg * b21_avg ) )
pi = acos(-1.0)
omega = 2.0 * pi / omega


print("Periodicity and decay time scale from the theory: ")

print("sigma: " + sigma)
print("omega: " + omega)

end
