import os

threads = [1, 2, 4]
sizes = [512, 1024, 1536, 2048, 2560, 3072, 4096, 8192, 16384]

path = "OpenMP"

os.system("mkdir -p " + path)
for thread in threads:
    for size in sizes:
        for i in range(3):
             os.system("mpisubmit.bg -n 1 -m smp -w 15:00"
	                   " -stdout " + path + "/" + str(thread) + "x" + str(size) + "_" + str(i) + ".out"
		               " -stderr " + path + "/" + str(thread) + "x" + str(size) + "_" + str(i) + ".err"
                       " sor_2d_omp " + str(thread) + " " + str(size + 2) + " " + str(size / 64))