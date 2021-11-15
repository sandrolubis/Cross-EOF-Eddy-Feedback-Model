begin

f1 = addfile("data/u.col.anom.erainterim.nc", "r")
f2 = addfile("data/duepy.col.anom.erainterim.nc", "r")

time = f1->time
lat = f1->lat({-80:-20})

u = f1->u(time|:, {lat|-80:-20}, lon|0)
duepy = f2->duepy(time|:, {latitude|-80:-20}, lon|0)

u@_FillValue = -999
duepy@_FillValue = -999

rad = get_d2r("float")

wgt = sqrt(cos(lat * rad))

u = (/u * conform(u, wgt, 1)/)
duepy = (/duepy * conform(duepy, wgt, 1)/)

neof = 2

eof = eofunc_n_Wrap(u, neof, False, 0)

pc = eofunc_ts_n_Wrap(u, eof, False, 0)
m = eofunc_ts_n_Wrap(duepy, eof, False, 0)

z1 = pc(0, :)
z2 = pc(1, :)
m1 = m(0, :)
m2 = m(1, :)

maxlag = 20

lag = ispan(0, maxlag, 1)

b1_samp = new((/2, dimsizes(lag)/), float, -999)
b2_samp = new((/2, dimsizes(lag)/), float, -999)

do l = 0, (dimsizes(lag) - 1)
	reg_z1_z1 = regCoef(z1(0:(dimsizes(z1) - 1 - l)), z1(l:(dimsizes(z1) - 1)))
	reg_z2_z1 = regCoef(z1(0:(dimsizes(z1) - 1 - l)), z2(l:(dimsizes(z2) - 1)))
	reg_z1_z2 = regCoef(z2(0:(dimsizes(z2) - 1 - l)), z1(l:(dimsizes(z1) - 1)))
	reg_z2_z2 = regCoef(z2(0:(dimsizes(z2) - 1 - l)), z2(l:(dimsizes(z2) - 1)))
	
	reg_m1_z1 = regCoef(z1(0:(dimsizes(z1) - 1 - l)), m1(l:(dimsizes(m1) - 1)))
	reg_m1_z2 = regCoef(z2(0:(dimsizes(z2) - 1 - l)), m1(l:(dimsizes(m1) - 1)))
	reg_m2_z1 = regCoef(z1(0:(dimsizes(z1) - 1 - l)), m2(l:(dimsizes(m2) - 1)))
	reg_m2_z2 = regCoef(z2(0:(dimsizes(z2) - 1 - l)), m2(l:(dimsizes(m2) - 1)))
	
	A1 = (/ (/reg_z1_z1, reg_z2_z1/), \
	        (/reg_z1_z2, reg_z2_z2/) /)
	B1 = (/reg_m1_z1, reg_m1_z2/)
	b1_samp(:, l) = solve_linsys(A1, B1)
	
	A2 = (/ (/reg_z1_z1, reg_z2_z1/), \
	        (/reg_z1_z2, reg_z2_z2/) /)
	B2 = (/reg_m2_z1, reg_m2_z2/)
	b2_samp(:, l) = solve_linsys(A2, B2)
end do

b11_samp = b1_samp(0, :)
b12_samp = b1_samp(1, :)
b21_samp = b2_samp(0, :)
b22_samp = b2_samp(1, :)

lagi = 7
lagf = 20

b11 = avg(b11_samp(lagi:lagf))
b12 = avg(b12_samp(lagi:lagf))
b21 = avg(b21_samp(lagi:lagf))
b22 = avg(b22_samp(lagi:lagf))

print("ERA-Interim: Lag " + lagi + "-" + lagf)

print("b11    = " + sprintf("%8.4f", b11))
print("b12    = " + sprintf("%8.4f", b12))
print("b21    = " + sprintf("%8.4f", b21))
print("b22    = " + sprintf("%8.4f", b22))

tau1 = 8.3
tau2 = 8.4

sigma = -0.5 * (1.0 / tau1 + 1.0 / tau2 - b11 - b22)

D = ((1.0 / tau1 - 1.0 / tau2) - (b11 - b22)) ^ 2.0 + 4.0 * b12 * b21

pi = acos(-1.0)

if (D.ge.0) then
	omega = 0.5 * sqrt(D)
	
	sigma1 = sigma + omega
	sigma2 = sigma - omega
	
	sigma1 = -1.0 / sigma1
	sigma2 = -1.0 / sigma2
	
	print("sigma1 = " + sprintf("%8.4f", sigma1))
	print("sigma2 = " + sprintf("%8.4f", sigma1))
else
	omega = 0.5 * sqrt(-1.0 * D)
	
	sigma = -1.0 / sigma
	omega = 2.0 * pi / omega
	
	print("sigma  = " + sprintf("%8.4f", sigma))
	print("omega  = " + sprintf("%8.4f", omega))
end if

end