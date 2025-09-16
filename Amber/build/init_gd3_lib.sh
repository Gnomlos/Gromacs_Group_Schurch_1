cat > build/gd3.lib <<'EOF'
!entry.GD3.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elm  dbl chg
!entry.GD3.unit
!entry.GD3.unit.atoms
 "GD" "Gd3+" 0 1 131073 1 64 3.0000
!entry.GD3.unit.connect array int
!entry.GD3.unit.residues table  str name  int seq  int chain  int res
!entry.GD3.unit.residues
 "GD3" 1 1 1
!entry.GD3.unit.positions table dbl x dbl y dbl z
!entry.GD3.unit.positions
  0.0  0.0  0.0
!entry.GD3.unit.connectivity table int atom1x int atom2x int flags
EOF
