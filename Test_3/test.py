mismatch = 0
for line in open('QueryOutput.txt').readlines():
    counts = line.split(' ')
    if counts[2] != counts[4] or counts[2] != counts[6] or counts[6] != counts[8]:
        print(line)
        mismatch+=1
print(mismatch)
    