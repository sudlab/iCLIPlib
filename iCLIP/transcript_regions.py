'''Functions in this module all take a list of CGAT.GTF.Entry()
and return a list of the same. Return lists always contain Entries
that are always marked as "exon". They are used for extracting
particular parts of transcripts. For example, CDS extracts
the parts of the transcript that correspond to the coding 
sequence of the transcript.

flank5 and flank3 also take a length setting for how much 
should be returned. 

CDS and the UTRs require there to be at least one CDS
entry in the transcripts.
'''

import cgat.GTF as GTF
import cgat.Intervals as Intervals


def flank5(transcript, length=500):
    
    exons = [e for e in transcript if e.feature == "exon"]
    
    if exons[0].strand == "-":
        start = max(x.end for x in exons)
        end = start + length
    else: 
        end = min(x.start for x in exons)
        start = end - length

    returned_exon = GTF.Entry().fromGTF(exons[0])
    returned_exon.start = start
    returned_exon.end = end
    
    return [returned_exon]


def flank3(transcript, length=500):
    
    exons = [e for e in transcript if e.feature == "exon"]
    
    if exons[0].strand == "+":
        start = max(x.end for x in exons)
        end = start + length
    else: 
        end = min(x.start for x in exons)
        start = end - length
        
    returned_exon = GTF.Entry().fromGTF(exons[0])
    returned_exon.start = start
    returned_exon.end = end
    
    return [returned_exon]


def CDS(transcript):
    
    CDS = [e for e in transcript if e.feature == "CDS"]

    if len(CDS) == 0:
        return list()
    
    returned_exons = [GTF.Entry().fromGTF(e) for e in CDS]
    for e in returned_exons:
        e.feature = "exon"
        
    return returned_exons


def UTR3(transcript):
    
    exons = GTF.asRanges(transcript, "exon")
    cds = GTF.asRanges(transcript, "CDS")

    if len(cds) == 0:
        return list()
    
    utrs = Intervals.truncate(exons, cds)

    if transcript[0].strand == "+":
        utr3 = [exon for exon in utrs
                if exon[0] >= cds[-1][1]]
    else:
        utr3 = [exon for exon in utrs
                if exon[-1] <= cds[0][0]]

    for e in transcript:
        if e.feature == "exon":
            template_exon = e
            break
            
    returned_exons = []     
    for e in utr3:
        gtf = GTF.Entry().fromGTF(template_exon)
        gtf.start = e[0]
        gtf.end = e[1]
        returned_exons.append(gtf)
        
    return returned_exons


def UTR5(transcript):
    
    exons = GTF.asRanges(transcript, "exon")
    cds = GTF.asRanges(transcript, "CDS")

    utrs = Intervals.truncate(exons, cds)

    if len(cds) == 0:
        return list()
    
    if transcript[0].strand == "-":
        utr3 = [exon for exon in utrs
                if exon[0] >= cds[-1][1]]
    else:
        utr3 = [exon for exon in utrs
                if exon[-1] <= cds[0][0]]

    for e in transcript:
        if e.feature == "exon":
            template_exon = e
            break
            
    returned_exons = []     
    for e in utr3:
        gtf = GTF.Entry().fromGTF(template_exon)
        gtf.start = e[0]
        gtf.end = e[1]
        returned_exons.append(gtf)
        
    return returned_exons


def introns(transcript):
    
    introns = GTF.toIntronIntervals(transcript)
    
    for e in transcript:
        if e.feature == "exon":
            template_exon = e
            break
            
    returned_exons = []     
    for e in introns:
        gtf = GTF.Entry().fromGTF(template_exon)
        gtf.start = e[0]
        gtf.end = e[1]
        returned_exons.append(gtf)
        
    return returned_exons


def first_exon(transcript):

    exons = [e for e in transcript if e.feature == "exon"]

    if len(exons) < 3:
        return list()
    
    if exons[0].strand == "-":
        exons.sort(key=lambda x: x.end, reverse=True)
    else:
        exons.sort(key=lambda x: x.start)
    return [exons[0]]


def last_exon(transcript):

    exons = [e for e in transcript if e.feature == "exon"]

    if len(exons) < 3:
        return list()

    if exons[0].strand == "-":
        exons.sort(key=lambda x: x.end, reverse=True)
    else:
        exons.sort(key=lambda x: x.start)
    return [exons[-1]]


def middle_exons(transcript):

    exons = [e for e in transcript if e.feature == "exon"]

    if len(exons) < 3:
        return list()

    exons.sort(key=lambda x: x.start)
    return exons[1:-1]


def exons(transcript):

    e = [e for e in transcript if e.feature == "exon"]
    return (e)

def primary_transcript(transcript):

    ex = exons(transcript)

    transcript = ex[0]
    transcript.start = min([e.start for e in ex])
    transcript.end = max(e.end for e in ex)

    return [transcript]

    
def tss(transcript, upstream=500, downstream=500):

    exons = [e for e in transcript if e.feature == "exon"]
    
    if exons[0].strand == "-":
        start = max(x.end for x in exons) - downstream
        end = start + upstream + downstream
    else: 
        end = min(x.start for x in exons) + downstream
        start = end - upstream - downstream

    returned_exon = GTF.Entry().fromGTF(exons[0])
    returned_exon.start = start
    returned_exon.end = end
    
    return [returned_exon]


def tts(transcript, upstream=500, downstream=500):

    exons = [e for e in transcript if e.feature == "exon"]
     
    if exons[0].strand == "+":
        start = max(x.end for x in exons) - upstream
        end = start + upstream + downstream
    else: 
        end = min(x.start for x in exons) + upstream
        start = end - upstream - downstream

    returned_exon = GTF.Entry().fromGTF(exons[0])
    returned_exon.start = start
    returned_exon.end = end
    
    return [returned_exon]
