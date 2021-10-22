def AllPossibleKmars(k):
    '''Provides a list of all possible kmars'''
    lst = ['A', 'G','T','C']
    for j in range(k-1):
        print('Working on: ',j)
        lst = lst* 4
        for i in range(0,int(len(lst)/4)):
            lst[i]+='A'
        for i in range(int(len(lst)/4),int(len(lst)/2)):
            lst[i]+='G'
        for i in range(int(len(lst)/2),int(3*len(lst)/4)):
            lst[i]+='T'
        for i in range(int(3*len(lst)/4),int(len(lst))):
            lst[i]+='C'
    return lst

lst = AllPossibleKmars(10)
print('10 mers are generated')
print(len(lst))

fhand = open('all13mers.txt','w')
for kmer in lst:
    fhand.write(kmer+'\n')

fhand.close()