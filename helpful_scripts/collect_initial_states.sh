#!/usr/bin/env bash

  for ii in {0..254}
  do
    mv gaussians_$ii.dat ../../Initial_states/0_5/
    mv remnants_$ii.dat ../../Initial_states/0_5/
  done

  for ii in {255..509}
  do
    mv gaussians_$ii.dat ../../Initial_states/5_10/
    mv remnants_$ii.dat ../../Initial_states/5_10/
  done

  for ii in {510..764}
  do
    mv gaussians_$ii.dat ../../Initial_states/10_15/
    mv remnants_$ii.dat ../../Initial_states/10_15/
  done

  for ii in {765..1019}
  do
    mv gaussians_$ii.dat ../../Initial_states/15_20/
    mv remnants_$ii.dat ../../Initial_states/15_20/
  done

  for ii in {1020..1274}
  do
    mv gaussians_$ii.dat ../../Initial_states/20_25/
    mv remnants_$ii.dat ../../Initial_states/20_25/
  done

  for ii in {1275..1529}
  do
    mv gaussians_$ii.dat ../../Initial_states/25_30/
    mv remnants_$ii.dat ../../Initial_states/25_30/
  done

  for ii in {1530..1784}
  do
    mv gaussians_$ii.dat ../../Initial_states/30_35/
    mv remnants_$ii.dat ../../Initial_states/30_35/
  done

  for ii in {1785..2039}
  do
    mv gaussians_$ii.dat ../../Initial_states/35_40/
    mv remnants_$ii.dat ../../Initial_states/35_40/
  done

  for ii in {2040..2294}
  do
    mv gaussians_$ii.dat ../../Initial_states/40_45/
    mv remnants_$ii.dat ../../Initial_states/40_45/
  done

  for ii in {2295..2549}
  do
    mv gaussians_$ii.dat ../../Initial_states/45_50/
    mv remnants_$ii.dat ../../Initial_states/45_50/
  done
