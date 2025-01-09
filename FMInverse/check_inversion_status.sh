#!/bin/bash
julia step7a_inverse_status.jl
grep --color=yes -B 2 -A 2 '%'  inverse_status.md

