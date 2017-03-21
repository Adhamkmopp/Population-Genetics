# ADMIXTURE

```
mkdir admixture
cd admixture
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do admixture --cv ../locality_grants.bed $i; done > locality_cvoutput
grep -i 'CV error' locality_cvoutput
```

```
CV error (K=2): 0.41512
CV error (K=3): 0.37669
CV error (K=4): 0.35243
CV error (K=5): 0.35868
CV error (K=6): 0.34950
CV error (K=7): 0.37009
CV error (K=8): 0.33952
CV error (K=9): 0.44841
CV error (K=10): 0.42623
CV error (K=11): 0.41128
CV error (K=12): 0.42997
CV error (K=13): 0.47111
CV error (K=14): 0.46702
CV error (K=15): 0.54337
```
