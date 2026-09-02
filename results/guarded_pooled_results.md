# pooled_switch run log (guarded metric: per-sample stage-1+agg, own vs pooled switch)
# generated 2026-08-28 23:02 by /tmp/pooled_switch.py (committed as python/guarded_pooled.py)
vbf: pairs 1,018,253
zjets: pairs 71,418
dijet: pairs 158,146
ttbar: pairs 1,104,657

===== switch training: own =====
  fold 0 done
  fold 1 done
  fold 2 done
  fold 3 done
  vbf    TRKPTZ  90.45%  [t=0.6]  93.93 (-2851/+20883)  [t=0.9]  92.95 (-338/+13302)  [t=0.95]  92.45 (-130/+10527)  [t=0.98]  91.84 (-35/+7234)
  zjets  TRKPTZ  61.87%  [t=0.6]  66.58 (-1057/+2779)  [t=0.9]  63.36 (-54/+598)  [t=0.95]  62.56 (-10/+263)  [t=0.98]  61.99 (-0/+45)
  dijet  TRKPTZ  86.88%  [t=0.6]  89.80 (-809/+3154)  [t=0.9]  88.54 (-63/+1394)  [t=0.95]  87.96 (-20/+881)  [t=0.98]  87.38 (-4/+404)
  ttbar  TRKPTZ  87.08%  [t=0.6]  90.48 (-4930/+24028)  [t=0.9]  88.86 (-384/+10397)  [t=0.95]  88.22 (-111/+6559)  [t=0.98]  87.64 (-21/+3177)

===== switch training: pooled =====
  fold 0 done
  fold 1 done
  fold 2 done
  fold 3 done
  vbf    TRKPTZ  90.45%  [t=0.6]  93.95 (-2793/+20976)  [t=0.9]  92.93 (-303/+13143)  [t=0.95]  92.43 (-119/+10395)  [t=0.98]  91.81 (-27/+7055)
  zjets  TRKPTZ  61.87%  [t=0.6]  66.98 (-898/+2767)  [t=0.9]  63.27 (-25/+538)  [t=0.95]  62.54 (-4/+250)  [t=0.98]  62.07 (-0/+74)
  dijet  TRKPTZ  86.88%  [t=0.6]  89.99 (-656/+3154)  [t=0.9]  88.50 (-40/+1340)  [t=0.95]  87.99 (-11/+901)  [t=0.98]  87.51 (-2/+509)
  ttbar  TRKPTZ  87.08%  [t=0.6]  90.49 (-4785/+23958)  [t=0.9]  88.83 (-355/+10204)  [t=0.95]  88.25 (-109/+6685)  [t=0.98]  87.70 (-23/+3497)

## v2: + per-candidate tagger aggregates + cross-candidate time consistency (2026-08-28 23:59)

vbf: pack built, 1,018,253 pairs
zjets: pack built, 71,418 pairs
dijet: pack built, 158,146 pairs
ttbar: pack built, 1,104,657 pairs

===== pooled switch, full feature set =====
  fold 0 done
  fold 1 done
  fold 2 done
  fold 3 done
  vbf    TRKPTZ  90.45%  [t=0.6]  94.02 (-2817/+21335)  [t=0.9]  93.01 (-326/+13597)  [t=0.95]  92.52 (-133/+10859)  [t=0.98]  91.88 (-34/+7427)
  zjets  TRKPTZ  61.87%  [t=0.6]  67.49 (-987/+3040)  [t=0.9]  63.50 (-34/+629)  [t=0.95]  62.73 (-8/+323)  [t=0.98]  62.17 (-0/+111)
  dijet  TRKPTZ  86.88%  [t=0.6]  90.12 (-695/+3295)  [t=0.9]  88.56 (-48/+1396)  [t=0.95]  88.06 (-13/+959)  [t=0.98]  87.54 (-4/+528)
  ttbar  TRKPTZ  87.08%  [t=0.6]  90.61 (-4877/+24723)  [t=0.9]  88.96 (-385/+10938)  [t=0.95]  88.35 (-107/+7240)  [t=0.98]  87.76 (-29/+3856)
