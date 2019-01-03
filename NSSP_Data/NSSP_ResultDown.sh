

datapath="siyoo@tiger:/BiO/BioProjects/TBD180101-SCHU-Fungi-smallRNA-20180322/do_summary"

echo "mkdir Raw" | /bin/bash
echo "mkdir Clean" | /bin/bash
echo "mkdir Mapping" | /bin/bash
echo "mkdir ExtConSeq" | /bin/bash

cmd="sshpass -p 'a;sldkfj' scp ${datapath}/Raw/Raw.SeqInfo.Report Raw/Raw.SeqInfo.Report"
echo ${cmd} | /bin/bash
cmd="sshpass -p 'a;sldkfj' scp ${datapath}/Clean/Clean.SeqInfo.Report Clean/Clean.SeqInfo.Report"
echo ${cmd} | /bin/bash
cmd="sshpass -p 'a;sldkfj' scp ${datapath}/Mapping/Mapping.MapInfo.Report Mapping/Mapping.MapInfo.Report"
echo ${cmd} | /bin/bash
cmd="sshpass -p 'a;sldkfj' scp ${datapath}/ExtConSeq/extconseq.SeqInfo.Report ExtConSeq/extconseq.SeqInfo.Report"
echo ${cmd} | /bin/bash
cmd="sshpass -p 'a;sldkfj' scp ${datapath}/smRNAs.xls ./smRNAs.xls"
echo ${cmd} | /bin/bash
