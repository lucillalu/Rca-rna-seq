import sys

def main():
    fp_in = sys.argv[1]
    fp_out = fp_in[:-5] + 'fasta'

    data = {}
    with open(fp_in,'r') as f:
        i = 0
        j = 0
        for line in f:
            line = line.strip()
            i += 1
            if line[0] == '@':
                data[j] = [0,0]
                data[j][0] = '>' + line[1:]
            elif (i % 4 == 2):
                data[j][1] = line
                j += 1

    with open(fp_out,'w') as o:
        for key,value in data.items():
            print(value[0],file=o)
            print(value[1],file=o)

if __name__ == "__main__":
    main()