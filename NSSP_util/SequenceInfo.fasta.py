
def run_clc_sequence_info(file):
    outfile = "{0}.SeqInfo".format(file)
    if not os.path.isfile(outfile):
        cmd = "/BiO/BioProjects/TBD180101-SCHU-Fungi-smallRNA-20180322/NSSP_util/"\
              "clc_sequence_info -r -n {0} > {1}".format(file, outfile)
        #print (cmd)
        os.system(cmd)
    else:
        pass
        #print ('file is exist')
    return outfile

def make_infoDict(seqinfo_file, infoDict):
    for line in open(seqinfo_file):
        line = line.strip()
        if line == '':
            continue
        items = line.split() 
        line_tag = items[0] 
        line_value = items[-1]
        if line_tag == "Number" :
            if not "%" in line_value:
                infoDict[seqinfo_file].setdefault('No_Reads', line_value)
            base_tag = items[2] 
            if "A" in base_tag:
                a = items[3]
            elif "C" in base_tag:
                c = items[3]
            elif "G" in base_tag:
                g = items[3]
            elif "T" in base_tag:
                t = items[3]
            elif "N" in base_tag:
                n = items[3]

        elif line_tag == "Total" :
            infoDict[seqinfo_file].setdefault("Residues", line_value)
            total = line_value
        elif line_tag == "Minimum" : 
            infoDict[seqinfo_file].setdefault("Minimum", line_value)
        elif line_tag == 'Maximum' : 
            infoDict[seqinfo_file].setdefault("Maximum", line_value)
        elif line_tag == 'Average' : 
            infoDict[seqinfo_file].setdefault("Average", line_value)
        elif line_tag == 'N50' : 
            infoDict[seqinfo_file].setdefault("N50", line_value)
        else :
            pass 

    gc = (int(g)+int(c))/float(total)*100
    npct = int(n)/float(total)*100
    infoDict[seqinfo_file].setdefault("Npct", npct)
    infoDict[seqinfo_file].setdefault("GC", gc)

    return infoDict
    
def makeFileList(path, extensions):
    files = []
    for extension in extensions:
        newfiles = glob.glob('{0}/*.{1}'.format(path, extension))
        files.extend(newfiles)
    return files

def number_format(n):
    s = str(n)
    if s.find('.') < 0:
        e = re.compile(r"(\d)(\d\d\d)$")
        s, cnt = re.subn(e, r"\1,\2", s)
    e = re.compile(r"(\d)(\d\d\d([.,]))")
    while 1:
        s, cnt = re.subn(e, r"\1,\2", s)
        if not cnt: break
    return s

def main(outprefix, files):
    infoDict = dict()
    for file in files:
        seqinfo_file = run_clc_sequence_info(file)
        infoDict.setdefault(seqinfo_file, {})
        infoDict = make_infoDict(seqinfo_file, infoDict)

    outfile = open("{0}.SeqInfo.Report".format(outprefix), "w")
    outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(\
            'SampleName', 'No.Seq', 'Residues', 'Average',\
            'Minimum', 'Maximum', 'N50', 'Npct', 'GC'
            ))
    for info_file, ainfoDict in sorted(infoDict.items()):
        outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(\
                info_file.split('/')[-1].split('.')[0],
                number_format(ainfoDict['No_Reads']),
                number_format(ainfoDict['Residues']),
                number_format(ainfoDict['Average']),
                number_format(ainfoDict['Minimum']),
                number_format(ainfoDict['Maximum']),
                number_format(ainfoDict['N50']),
                round(ainfoDict['Npct'], 1), round(ainfoDict['GC'], 1)
                ))

if __name__=='__main__':
    import glob
    import os
    import sys
    import re
    import argparse
    parser = argparse.ArgumentParser(description='Run clc_sequence_info and Merge')
    parser.add_argument('-p', '--outprefix', type=str, help='outprefix as full PATH')
    parser.add_argument('-fs', '--files', nargs='+')
    args = parser.parse_args()
    main(args.outprefix, args.files)
