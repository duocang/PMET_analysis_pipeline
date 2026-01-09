import sys
import os

def find_isoforms(input_file, longest_output):
    longest_isoforms = {}
    shorter_isoforms = {}

    # 创建较短的 isoforms 文件的路径
    dir_name, file_name = os.path.split(longest_output)
    shorter_output = os.path.join(dir_name, "shorter_" + file_name)

    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            isoform_key = parts[3].split('.')[0]
            length = int(parts[2]) - int(parts[1])
            if isoform_key not in longest_isoforms or length > longest_isoforms[isoform_key][1]:
                if isoform_key in longest_isoforms:
                    shorter_isoforms.setdefault(isoform_key, []).append(longest_isoforms[isoform_key][0])
                longest_isoforms[isoform_key] = (line.strip(), length)
            else:
                shorter_isoforms.setdefault(isoform_key, []).append(line.strip())

    with open(longest_output, 'w') as f:
        for line, _ in longest_isoforms.values():
            f.write(line + '\n')

    with open(shorter_output, 'w') as f:
        for lines in shorter_isoforms.values():
            for line in lines:
                f.write(line + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <longest_output_file>")
        sys.exit(1)

    input_file, longest_output = sys.argv[1], sys.argv[2]
    find_isoforms(input_file, longest_output)
