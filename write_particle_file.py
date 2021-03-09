import argparse
import numpy as np
from math import sqrt, log, prod

# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-o", "--output", required = True,
                help = "path to output file")
#ap.add_argument("-n", "--number_of_particles", type = int, required = True,
#                help = "number of particles")
ap.add_argument("-n", "--nb_of_particles", type = int, required = True,
                nargs = 3,
                help = "number of particles in x, y, z direction")
ap.add_argument("-s", "--start_point", type = float, default = [0.0, 0.0, 0.0],
                nargs = 3,
                help = "center of particle mass in x, y, z direction")
ap.add_argument("-e", "--epsilon", type = float, default = 1.0,
                help = "depth of Lennard-Jones potential")
ap.add_argument("-r", "--radius", type = float, default = 1.0,
                help = "radius of particle")
ap.add_argument("-m", "--mass", type = float, default = 1.0,
                help = "mass of particle")
ap.add_argument("-v", "--velocity", type = float, default = [0.0, 0.0, 0.0],
                nargs = 3, help = "General velocity of all particles")
ap.add_argument("-vt", "--thermal_noise", type = float, default = 0.0,
                help = "Amount of thermal noise of particles, i.e. temperature")
args = vars(ap.parse_args())

# compute r_min, the position of the minimum of the Lennard-Jones potential
r_min = args["radius"] * 2.0**(1.0/6.0)
start_point = np.array(args["start_point"])
velocity = np.array(args["velocity"])
total_particles = prod(args["nb_of_particles"])
vt = args["thermal_noise"]

def thermal_velocity_noise():
    while True:
        thermal_noise = np.random.uniform(size = 3)
        r = np.sum(thermal_noise**2)
        if r > 0.0 and r < 1.0:
            break

    r = sqrt(-2.0 * log(r)/r)
    return thermal_noise * r

# open output file and write to it
with open(args["output"], "w") as file:

    file.write(str(total_particles))
    file.write("\n")

    id = 0
    for z in range(args["nb_of_particles"][2]):
        for y in range(args["nb_of_particles"][1]):
            for x in range(args["nb_of_particles"][0]):

                id += 1
                position = start_point + (r_min * np.array([x, y, z]))
                velocity = velocity + vt * thermal_velocity_noise()

                file.write(str(id) + " ")
                np.savetxt(file, position, newline=" ", fmt ='%.4f')
                np.savetxt(file, velocity, newline=" ", fmt = "%.4f")
                file.write("\n")
