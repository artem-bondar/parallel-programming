import os

threads = [1, 2, 4, 8]
sizes = [512, 1024, 1536, 2048, 2560, 3072, 4096, 8192, 16384, 32768]

path = "OpenMP"

os.system("mkdir -p " + path)
for thread in threads:
    for size in sizes:
        for i in range(3):
             os.system("./sor_2d_omp " + str(thread) + " " + str(size + 2) + " " + str(size / 64) +
	                   " 1> " + path + "/" + str(thread) + "x" + str(size) + "_" + str(i) + ".out"
		               " 2> " + path + "/" + str(thread) + "x" + str(size) + "_" + str(i) + ".err")