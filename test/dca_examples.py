from petbox import dca

#
ftime = dca.get_time()


#
thm = dca.THM(qi=750, Di=.8, bi=2, bf=.5, telf=28)
q_trans = thm.transient_rate(ftime)
N_trans = thm.transient_cum(ftime)
D_trans = thm.transient_D(ftime)
b_trans = thm.transient_b(ftime)
beta_trans = thm.transient_beta(ftime)

#
q_thm = thm.rate(ftime)
N_thm = thm.cum(ftime)
D_thm = thm.D(ftime)
b_thm = thm.b(ftime)
beta_thm = thm.beta(ftime)

#
mh = dca.MH(qi=725, Di=0.85, bi=0.6, Dterm=0.2)
q_mh = mh.rate(ftime)
N_mh = mh.cum(ftime)
N_mh = mh.monthly_vol(ftime)
N_mh = mh.interval_vol(ftime)
D_mh = mh.D(ftime)
b_mh = mh.b(ftime)
beta_mh = mh.beta(ftime)


#
ple = dca.PLE(qi=750, Di=.1, Dinf=.00001, n=.5)
q_ple = ple.rate(ftime)
N_ple = ple.cum(ftime)
N_ple = ple.monthly_vol(ftime)
N_ple = ple.interval_vol(ftime)
D_ple = ple.D(ftime)
b_ple = ple.b(ftime)
beta_ple = ple.beta(ftime)


#
se = dca.SE(qi=715, tau=90.0, n=.5)
q_se = se.rate(ftime)
N_se = se.cum(ftime)
N_se = se.monthly_vol(ftime)
N_se = se.interval_vol(ftime)
D_se = se.D(ftime)
b_se = se.b(ftime)
beta_se = se.beta(ftime)


#
dg = dca.Duong(qi=715, a=2.8, m=1.4)
q_dg = dg.rate(ftime)
N_dg = dg.cum(ftime)
N_dg = dg.monthly_vol(ftime)
N_dg = dg.interval_vol(ftime)
D_dg = dg.D(ftime)
b_dg = dg.b(ftime)
beta_dg = dg.beta(ftime)


#
thm.add_secondary(dca.Yield(c=1000, m0=-.1, m=.5, t0=2 * 365.25 / 12))
g = thm.secondary.rate(ftime)
y = thm.secondary.gor(ftime)
g = thm.secondary.cum(ftime)
g = thm.secondary.monthly_vol(ftime)
g = thm.secondary.interval_vol(ftime, t0=0.)
g = thm.secondary.D(ftime)
g = thm.secondary.beta(ftime)
g = thm.secondary.b(ftime)
