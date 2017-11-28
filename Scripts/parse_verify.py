result = open("verify.txt","w")

sizes = [512, 1024, 1536, 2048, 2560, 3072, 4096, 8192, 16384, 32768]

path = "Verify"

for size in sizes:
	file = open(path + "/" + str(size) + ".out", "r")
	out = file.readlines()
	result.write(str(size) + " " + out[1].split(" ")[1][:-1] + " " + out[2].split(" ")[1][:-1] + "\n")
	file.close()
result.close()