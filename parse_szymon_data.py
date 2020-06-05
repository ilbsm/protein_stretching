import sys
path = sys.argv[1]
f_prefix = sys.argv[2]
start = sys.argv[3]
end = sys.argv[4]
f_suffix = sys.argv[5]
if len(sys.argv) > 6:
    cg = bool(int(sys.argv[6]))
else:
    cg = True
if len(sys.argv) > 7:
    output_name = sys.argv[7]
else:
    output_name = f_prefix.strip('_') + '.csv'

if cg:
    cols = [6, 7]
else:
    cols = [1, 3]


max_lines = 0
for k in range(int(start), int(end) + 1):
    if cg and k < 10:
        number = '0' + str(k)
    else:
        number = str(k)
    fname = path + '/' + f_prefix + number + f_suffix
    count = 0
    with open(fname, 'r') as myfile:
        for line in myfile.readlines():
            count += 1
    if count > max_lines:
        max_lines = count

data = [[] for k in range(max_lines)]

for k in range(int(start), int(end) + 1):
    if cg and k < 10:
        number = '0' + str(k)
    else:
        number = str(k)
    fname = path + '/' + f_prefix + number + f_suffix
    count = 0
    with open(fname, 'r') as myfile:
        for line in myfile.readlines():
            data[count].append(line.split()[cols[0]].strip())
            data[count].append(line.split()[cols[1]].strip())
            count += 1
        for l in range(count, max_lines):
            data[l].append('')
            data[l].append('')

result = '\n'.join([';'.join(line) for line in data])
with open(output_name, 'w') as ofile:
    ofile.write(result)
