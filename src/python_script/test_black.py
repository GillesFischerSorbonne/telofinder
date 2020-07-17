with open("telom_length.csv", "a") as filout:
    if file_exists:
        filout.write(
            "{0}\t{1}\tL\t{2}\t{4}\n{0}\t{1}\tR\t{3}\t{5}\n".format(
                strain, chrom, left_tel, right_tel, offset, revoffset
            )
        )
    else:
        filout.write(
            "{0}\t{1}\t{2}\t{3}\t{4}\n{5}\t{6}\tL\t{7}\t{9}\n{5}\t{6}\tR\t{8}\t{10}\n".format(
                "Strain",
                "Chromosome",
                "Contig_side",
                "Telom_length",
                "Offset",
                strain,
                chrom,
                left_tel,
                right_tel,
                offset,
                revoffset,
            )
        )

        print(
            "fdqf f qzd fqdf qf qzf zf zf zf zfz rf z zrf rzf zrf zr ",
            "frzfzrfzarfrafarfregfregerag\
            eageargaegaeg a gar gar g areg ae g ge gear g",
        )
