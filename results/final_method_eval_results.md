# Final method, all samples -- condor cluster 2951035 (2026-08-30)
# 8 cpus, 39.0/48 GB peak, 2d07h user CPU

[job] host workers-ads-69cf697649-lrg6q, start Sun Aug 30 11:59:40 UTC 2026

===== vbf =====
  e2e fold 0: val-best 94.40%
  e2e fold 1: val-best 94.57%
  e2e fold 2: val-best 94.53%
  e2e fold 3: val-best 94.31%

vbf: TRKPTZ 90.45  ctime 93.77  base 94.99  t0f 93.55  t0e2e 94.45  union 95.91

===== zjets =====
  e2e fold 0: val-best 69.06%
  e2e fold 1: val-best 68.68%
  e2e fold 2: val-best 68.97%
  e2e fold 3: val-best 68.79%

zjets: TRKPTZ 61.87  ctime 67.78  base 69.43  t0f 67.56  t0e2e 68.72  union 72.79

===== dijet =====
  e2e fold 0: val-best 90.35%
  e2e fold 1: val-best 90.05%
  e2e fold 2: val-best 90.27%
  e2e fold 3: val-best 89.80%

dijet: TRKPTZ 86.88  ctime 89.87  base 90.78  t0f 89.52  t0e2e 90.12  union 92.25

===== ttbar =====
  e2e fold 0: val-best 90.95%
  e2e fold 1: val-best 90.92%
  e2e fold 2: val-best 90.88%
  e2e fold 3: val-best 90.94%

ttbar: TRKPTZ 87.08  ctime 90.42  base 91.47  t0f 89.97  t0e2e 90.92  union 92.77

===== vbf_mu0 =====
vbf_mu0: TRKPTZ 99.89  ctime 99.92  base 99.90  t0f 99.88  t0e2e nan  union 99.91

===== zeejets_mu0 =====
zeejets_mu0: TRKPTZ 99.91  ctime 99.93  base 99.93  t0f 99.91  t0e2e nan  union 99.98

===== ttbar_mu0 =====
ttbar_mu0: TRKPTZ 99.91  ctime 99.94  base 99.92  t0f 99.91  t0e2e nan  union 99.93

===== pooled meta (mu200), per-sample evaluation =====
  vbf          TRKPTZ  90.45  [0.35]  95.19  [0.4]  95.22  [0.45]  95.26  [0.5]  95.26  [0.55]  95.25  [0.6]  95.22
  zjets        TRKPTZ  61.87  [0.35]  69.51  [0.4]  69.69  [0.45]  69.76  [0.5]  69.79  [0.55]  69.73  [0.6]  69.64
  dijet        TRKPTZ  86.88  [0.35]  90.90  [0.4]  90.96  [0.45]  91.00  [0.5]  91.00  [0.55]  91.02  [0.6]  90.98
  ttbar        TRKPTZ  87.08  [0.35]  91.64  [0.4]  91.70  [0.45]  91.74  [0.5]  91.75  [0.55]  91.74  [0.6]  91.71
  vbf_mu0      TRKPTZ  99.89  [0.35]  99.91  [0.4]  99.91  [0.45]  99.91  [0.5]  99.91  [0.55]  99.90  [0.6]  99.90
  zeejets_mu0  TRKPTZ  99.91  [0.35]  99.96  [0.4]  99.93  [0.45]  99.93  [0.5]  99.93  [0.55]  99.93  [0.6]  99.93
  ttbar_mu0    TRKPTZ  99.91  [0.35]  99.93  [0.4]  99.93  [0.45]  99.93  [0.5]  99.93  [0.55]  99.93  [0.6]  99.92
wrote fe_results.json
[job] done Sun Aug 30 20:17:19 UTC 2026
