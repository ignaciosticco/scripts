# Lammps For

variable imax equal 14

variable i loop ${imax}

label a

variable b equal $i-$i/2

print "Iteration i = $b"

next i

jump SELF a