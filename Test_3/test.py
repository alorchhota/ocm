mismatch = 0
for line in open('QueryOutput.txt').readlines():
    counts = line.split(' ')
    if counts[4] != counts[6] or counts[4] != counts[8] or counts[8] != counts[10]:
        print(line)
        mismatch+=1
print(mismatch)
    