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

cospc_total = 0.0

quspc_total = 0.0

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
	
	cospc_total = cospc_total + specxy_total@cospc
	
	quspc_total = quspc_total + specxy_total@quspc
	
	spectra_z1 = spectra_z1 + specxy_total@spcx
	
	ntStrt = ntLast - nSampOverlap + 1
	ntLast = ntStrt + nSampWin - 1
end do

cospc_total = cospc_total / nWindow

quspc_total = quspc_total / nWindow

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

print("alpha = " + alpha)
print("beta = " + beta)

tau = beta / alpha

print("tau = " + tau)

end
