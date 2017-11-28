result = open("results.txt","w")

threads = [1, 2, 4, 8]
sizes = [512, 1024, 1536, 2048, 2560, 3072, 4096, 8192]

path = "MPI"

for thread in threads:
    for size in sizes:
        for i in range(3):
			file = open(path + "/" + str(thread) + "x" + str(size) + "_" + str(i) + ".out", "r")
			out = file.readlines()
			result.write(str(thread) + " " + str(size) + " " + out[1].split(" ")[1][:-1] + " " + out[2].split(" ")[1][:-1] + "\n")
			file.close()
result.close()