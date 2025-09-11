f = open("log_jst.txt", "r")
lines = f.readlines()
f.close()

CD = []
for line in lines:
    if "Cd       :" in line:
        cols = line.split()
        CD.append(float(cols[2]))
print(CD)
CD_Final = CD[-1]

if abs(CD_Final - 0.27387) / 0.27387 > 1e-7:
    print("HiSA test failed!")
    exit(1)
else:
    print("HiSA test passed!")
