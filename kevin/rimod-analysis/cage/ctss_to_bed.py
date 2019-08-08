import os
import glob

out_dir = "/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/CAGEseq/bed_files/"

ctss_path = "/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/CAGEseq/all_ctss/*fro*.ctss"
ctss_files = glob.glob(ctss_path)



# Loop over all files
count = 0
for file in ctss_files:
    count += 1
    f = open(file)
    bn = os.path.basename(file)
    line = f.readline()

    # open out file
    outfile = out_dir + bn + ".bed"
    out = open(outfile, "w")

    print(f"Processing file {count}")
    while line:
        # write new line into out file
        ls = line.split()
        chr = ls[0]
        start = ls[1]
        end = str(int(ls[1]) + 1)
        strand = ls[2]
        score = ls[3]
        nl = "\t".join([chr, start, end, bn, score, strand]) + "\n"

        out.write(nl)
        line = f.readline()

    f.close()
    out.close()
