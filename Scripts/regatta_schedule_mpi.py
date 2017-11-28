import os

threads = [1, 2, 4, 8, 16]
sizes = [512, 1024, 1536, 2048, 2560, 3072, 4096, 8192, 16384]

path = "MPI"

os.system("mkdir -p " + path)
for thread in threads:
    for size in sizes:
        for i in range(3):
             os.system("mpisubmit -n " + str(thread) + " -w 15:00"
	                   " -stdout " + path + "/" + str(thread) + "x" + str(size) + "_" + str(i) + ".out"
		               " -stderr " + path + "/" + str(thread) + "x" + str(size) + "_" + str(i) + ".err"
                       " sor_2d_mpi " + " " + str(size + 2))