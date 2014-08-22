import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineUtilities as PUtils
import pandas
import collections

def removeFirstAndLastExon(infile, outfile):

    transcripts = GTF.transcript_iterator(
        GTF.iterator(IOTools.openFile(infile)))
    outfile = IOTools.openFile(outfile, "w")

    for transcript in transcripts:

        for exon in transcript[1:-1]:
            outfile.write(str(exon) + "\n")

    outfile.close()


def getBarcodeCG(table,outfile):
    ''' Annotate barcode use statistics with %GC '''

    statement = " SELECT * FROM %(table)s" % locals()

    umi_stats = PUtils.fetch_DataFrame(statement)

    def _GC(x):
        return float(x.count("G") + x.count("G"))/len(x)

    barcode_gc = umi_stats.Barcode.apply(_GC)
    sample_gc = umi_stats.Sample.apply(_GC)
    umi_gc = umi_stats.UMI.apply(_GC)

    gc_stats= pandas.DataFrame({"Barcode":umi_stats.Barcode,
                                "barcode_gc": barcode_gc,
                                "sample_gc": sample_gc,
                                "umi_gc": umi_gc})

    gc_stats.to_csv(IOTools.openFile(outfile,"w"), sep="\t", index=False)
