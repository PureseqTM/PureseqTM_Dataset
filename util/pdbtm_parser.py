#!/usr/bin/env python
from xml.dom import minidom
import os, sys

#------- usage -------#
def Usage():
    print 'python pdbtm_parser.py <pdbtm_xml> <output_dir>'


#------- main --------#
def main(argv):
    if len(argv) != 2:
        Usage()
        sys.exit(-1)

    #-> input arguments
    xmlFile = argv[0]
    outFolder = argv[1]

    #-> mkdir
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)

    #-> parse input XML file
    xmldoc = minidom.parse(xmlFile)
    pdbtm_all = xmldoc.documentElement.getElementsByTagName("pdbtm")

    pSummary = []
    for pdbtm in pdbtm_all:
        #-> get ID, CHAIN, modification
        pName = pdbtm.getAttribute("ID")
        chains = pdbtm.getElementsByTagName("CHAIN")
        modification = pdbtm.getElementsByTagName("MODIFICATION")
        if len(modification) > 0:
            lastModification = modification[-1]
            lastModificationTime = lastModification.getElementsByTagName("DATE")[0].childNodes[0].data
        else:
            lastModificationTime = pdbtm.getElementsByTagName("CREATE_DATE")[0].childNodes[0].data

        #-> process each chain
        for chain in chains:
            chainID = chain.getAttribute("CHAINID")
            numOfTM = int( chain.getAttribute("NUM_TM") )
            chainType = chain.getAttribute("TYPE")
            chainSeq = chain.getElementsByTagName("SEQ")[0].childNodes[0].data.strip()
            chainSeq = ''.join(chainSeq.split())
            regions = chain.getElementsByTagName("REGION")
            chainLabel = ['X'] * len(chainSeq) 
            for region in regions:
                regionBegin = int( region.getAttribute("seq_beg") )
                regionEnd = int( region.getAttribute("seq_end") )
                regionType = region.getAttribute("type")
                chainLabel[regionBegin-1: regionEnd] = [regionType] * (regionEnd - regionBegin + 1)
            chainLabel = ''.join(chainLabel)
            if 'X' in chainLabel:
                missingPDB = True
            else:
                missingPDB = False

            #-> record to pSummary
            pSummary.append( (pName + chainID, len(chainSeq), numOfTM, chainType, missingPDB, lastModificationTime) )
            with open(outFolder + '/' + pName + chainID + '.pdbtm', 'w') as f:
                print >> f, '>' + pName + chainID
                print >> f, chainSeq
                print >> f, chainLabel

    #-> output to file
    pSummary.sort(key= lambda x: x[0])
    with open('pdbtm.summary', 'w') as f:
        print >> f, "#Target        Length             TM          type           MISS      Date"
        for query in pSummary:
            print >> f, "%s%15d%15d%15s%15s%15s"%query


#-------------- python main ---------------#
if __name__ == "__main__":
    main(sys.argv[1:])


