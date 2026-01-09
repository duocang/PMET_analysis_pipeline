import argparse
import numpy as np
import os

def count_motifs_in_file(file_content):
    """Count the number of occurrences of 'MOTIF' in the file content. 计算文件内容中 'MOTIF' 出现的次数。
    Args:
        file_content (list of str): Content of the file as a list of lines. 文件内容，作为行列表。

    Returns:
        int: The number of motifs found.
                找到的motifs数量。
    """
    return sum('MOTIF' in line for line in file_content)


def distribute_motifs_evenly(bigfile, threads, outdir):
    """Distribute motifs across multiple files as evenly as possible. 尽可能均匀地将motifs分布到多个文件中。
    Args:
        bigfile (np.ndarray): The content of the original file loaded into a NumPy array. 原始文件内容加载到NumPy数组中。
        threads (int): The number of separate files to distribute the motifs across. 要将motifs分布到的独立文件数量。
        outdir (str): The directory where the output files will be saved. 输出文件将被保存的目录。

    Returns:
        None
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Find the location of the first MOTIF.
    headind = next(i for i, line in enumerate(bigfile) if 'MOTIF' in line)
    meme_header_lines = bigfile[0:headind]

    # Initially create all the output files and write the header lines.
    files = [open(os.path.join(outdir, f"{i}.txt"), 'w') for i in range(threads)]
    for f in files:
        f.writelines(meme_header_lines)
        f.close()

    file_counter = 0
    ind1 = headind
    for i in range(ind1 + 1, len(bigfile)):
        if 'MOTIF' in bigfile[i] or i == len(bigfile) - 1:
            ind2 = i if 'MOTIF' in bigfile[i] else i + 1  # Include the last line.

            # Convert the found MOTIF interval to uppercase.
            bigfile[ind1] = bigfile[ind1].upper()

            lines_to_write = bigfile[ind1:ind2]
            with open(os.path.join(outdir, f"{file_counter}.txt"), 'a') as writer:
                writer.writelines(lines_to_write)

            # Rotate the file_counter within the range [0, threads-1].
            file_counter = (file_counter + 1) % threads

            ind1 = ind2


def main():
    """The main function to parse arguments and call the distribution function."""
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=argparse.FileType('r'), help="Input file path. 输入文件路径")
    parser.add_argument('outdir', type=str, help="Output directory. 输出目录")
    parser.add_argument('batch', type=int, default=8, help="The number of motifs in each batch. 每批次中的motifs数量")
    args = parser.parse_args()

    bigfile = np.asarray(args.file.readlines())
    distribute_motifs_evenly(bigfile, args.batch, args.outdir)

if __name__ == "__main__":
    main()
