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

cov_z1_z1 = esccv(z1, z1, maxlag)
cov_z2_z1 = esccv(z1, z2, maxlag)
cov_z1_z2 = esccv(z2, z1, maxlag)
cov_z2_z2 = esccv(z2, z2, maxlag)

cov_m1_z1 = esccv(z1, m1, maxlag)
cov_m1_z2 = esccv(z2, m1, maxlag)
cov_m2_z1 = esccv(z1, m2, maxlag)
cov_m2_z2 = esccv(z2, m2, maxlag)

b1_samp = new((/2, dimsizes(lag)/), float, -999)
b2_samp = new((/2, dimsizes(lag)/), float, -999)

do l = 0, (dimsizes(lag) - 1)
	A1 = (/ (/cov_z1_z1(l), cov_z2_z1(l)/), \
	        (/cov_z1_z2(l), cov_z2_z2(l)/) /)
	B1 = (/cov_m1_z1(l), cov_m1_z2(l)/)
	b1_samp(:, l) = solve_linsys(A1, B1)
	
	A2 = (/ (/cov_z1_z1(l), cov_z2_z1(l)/), \
	        (/cov_z1_z2(l), cov_z2_z2(l)/) /)
	B2 = (/cov_m2_z1(l), cov_m2_z2(l)/)
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


end
