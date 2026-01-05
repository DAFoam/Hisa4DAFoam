f = open("log.txt", "r")
lines = f.readlines()
f.close()

CD = []
for line in lines:
    if "Cd:" in line:
        cols = line.split()
        CD.append(float(cols[1]))
print(CD)
CD_Final = CD[-1]

if abs(CD_Final - 0.018972611) / 0.018972611 > 1e-7:
    print("HiSA test failed!")
    exit(1)
else:
    print("HiSA test passed!")
