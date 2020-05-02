import sys

def main():
    fp_in = sys.argv[1]
    fp_out = fp_in[:-3] + '_3.fa'

    data = {}
    with open(fp_in,'r') as f:
        i = 0
        for line in f:
            line = line.strip()
            if line[0] == '>':
                data[i] = [0,0]
                data[i][0] = line
            else:
                new_line = line + line + line
                data[i][1] = new_line
                i += 1

    with open(fp_out,'w') as o:
        for key,value in data.items():
            print(value[0],file=o)
            print(value[1],file=o)

if __name__ == "__main__":
    main()