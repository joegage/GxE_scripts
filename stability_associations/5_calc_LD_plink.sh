#!/usr/bin/env bash

plink --file ./hybrid_geno_g2f \
  --make-founders first \
  --r2 square  \
  --ld-window-kb 5 \
  --ld-window 9999999 \
  --ld-window-r2 0 \
  --out hybrid_geno_g2f
