import math
import subprocess
import numpy as np
import time


start_time = time.time()



Deme_dimension = 4**2 
if Deme_dimension < 20:
    Augmentation_factor = 10
elif Deme_dimension < 50:
    Augmentation_factor = 3
else:
    Augmentation_factor = 1



print(f"Augmentation_factor:{Augmentation_factor}")
SampleCount = Augmentation_factor*Deme_dimension 
pop_deme = int(100*2000/Deme_dimension)  


Trial = 5
mutarate = 1e-6
lociCount = 5000
theta = 4*pop_deme*mutarate*lociCount
Mirgation_rate = 4*pop_deme*0.15
Recombination_rate = 4*pop_deme*0.05
Homo_select = 2*pop_deme*0.02
Hetero_select = 2*pop_deme*0.01

## Migration matrix
z = math.sqrt(Deme_dimension)
M = np.zeros((Deme_dimension, Deme_dimension)).astype(int)
# Fill in the migration matrix.
for i in range(Deme_dimension):
    for j in range(Deme_dimension):
        # If deme j is immediately to the right of deme i and i is not on the right edge of the grid,
        # then set the migration rate from i to j to be the defined migration rate.
        if j == i+1 and (i+1) % z != 0:
            M[i,j] = int(Mirgation_rate)
        # If deme j is directly below deme i, then set the migration rate from i to j to be the defined migration rate.
        elif j == i+z:
            M[i,j] = int(Mirgation_rate)
        
        # Since migration is bidirectional, set the migration rate from j to i to be the same as from i to j.
        M[j,i] = M[i,j]

Ma = M.tolist()
for i in range(Deme_dimension):
    for j in range(Deme_dimension):
        if i == j:
            Ma[i][j] = 'x'
## end of Migration matrix

newMa = [" ".join(str(x) for x in row) for row in Ma]
newMa = " ".join(newMa)

command = f"msms {SampleCount} {Trial} -N {pop_deme} -t {theta} -r {Recombination_rate} {lociCount} -Sp 0.8 -I {Deme_dimension} "

for i in range(Deme_dimension):
    command += f"{Augmentation_factor} "
command += "-ma " + newMa + f" -SAA {Homo_select} -SaA {Hetero_select} -SF 0"


# Run the command and save the output
print(len(command))
print(command)
result = subprocess.run(command, stdout=subprocess.PIPE)



# Save the output to a file
with open('4x4spatial_test2.txt', 'w') as f:
    f.write(result.stdout.decode())

end_time = time.time()
run_time = end_time - start_time
print(f"Run time: {run_time:.2f} seconds")

with open('msmsrun.txt', 'w') as f:
    f.write(f"Run time: {run_time:.2f} seconds")