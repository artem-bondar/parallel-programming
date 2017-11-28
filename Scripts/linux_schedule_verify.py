import os

sizes = [512, 1024, 1536, 2048, 2560, 3072, 4096, 8192, 16384, 32768]
path = "Verify"

os.system("mkdir -p " + path)
for size in sizes:
     os.system("./sor_2d_verify " + str(size + 2) +
	           " 1> " + path + "/" + str(size) + ".out"
		       " 2> " + path + "/" + str(size) + ".err")