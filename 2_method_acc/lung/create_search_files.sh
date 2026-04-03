#!/usr/bin/bash

for f in {scissor,scab,scpas,scipac,scpp,degas,lp_sgl,pipet}_search.R; do
  cp template.R "$f"
done
