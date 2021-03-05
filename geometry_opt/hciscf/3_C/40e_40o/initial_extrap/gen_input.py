"""
Generates input.dat files for Dice calculations.
"""

input_header = """
#system
nocc 40
0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 
end
orbitals ../FCIDUMP_INITIAL
nroots 1

#variational
schedule 
"""

input_footer = """

davidsonTol 5e-05
dE 1e-08
maxiter 21

#pt
nPTiter 1000
epsilon2 1e-07
epsilon2Large 1e-5
targetError 2e-5
sampleN 1000

#misc
noio
"""


def write_schedule(eps1: float) -> str:
    schedule = ""
    for i in range(10):
        schedule += f"{i}   {(9-i)*1e-5 + eps1:.2e}\n"
    return schedule + "end\n\n"


def write_input(input_name: str, eps1: float):

    with open(input_name, "w") as f:
        f.write(input_header)
        f.write(write_schedule(eps1))
        f.write(input_footer)


if __name__ == "__main__":
    import os
    import numpy as np

    eps1s = np.linspace(5e-5, 1e-5, num=10)
    print(eps1s)

    for i, eps1 in enumerate(eps1s):
        write_input(f"input_{i}.dat", eps1)

