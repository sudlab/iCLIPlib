import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.PipelineUtilities as PUtils
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


def getBarcodeCG(table):
    ''' Annotate barcode use statistics with %GC '''

    statement = " SELECT * FROM %(table)s" % locals()

    umi_stats = PUtils.fetch_DataFrame(statement)

    def _GC(x):
        return (x.count("G") + x.count("G"))/len(x)

    barcode_gc = umi_stats.barcode.apply(_GC)
    
