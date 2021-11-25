import os,sys,re
import subprocess
import shutil
from tqdm import tqdm as tq
from operator import itemgetter
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser as SFP
from collections import defaultdict

def grep_pattern(line,P):
    if not P:
        yield None
    else:
        line=line[:20]+'_'+line[20]+'_'+line[21:]
        p=re.compile(P,re.IGNORECASE)
        yield p.search(line)

def find_p(name):
    # check the program path
    loc = shutil.which(name)
    if loc != None:
        try:
            help_m='-h'
            if 'samtools' in name:
                help_m='depth'
            o=subprocess.check_output([loc,help_m],stderr=subprocess.STDOUT,shell=False,universal_newlines=True,encoding='utf-8')
            return loc
        except:
            print("Please check the {0} !\n".format(name))
            sys.exit()
    else:
        print("Please set the {0} path !\n".format(name))
        sys.exit()

def Rename(Motif):
    return Motif.replace('[A,T,G,C]','N').replace('[A,T]','W').replace('[A,G]','R').replace('[C,T]','Y').replace('[G,C]','S').replace('[G,T]','K').replace('[A,C]','M').replace('[A,C,T]','H').replace('[A,C,G]','V').replace('[A,G,T]','D').replace('[C,G,T]','B')

def Rename2(Motif):
    return Motif.replace('N','[A,T,G,C]').replace('W','[A,T]').replace('R','[A,G]').replace('Y','[C,T]').replace('S','[G,C]').replace('K','[G,T]').replace('M','[A,C]').replace('H','[A,C,T]').replace('V','[A,C,G]').replace('D','[A,G,T]').replace('B','[C,G,T]')

def rep_name(p,m):
    if not p:
        motif_rep=('None/'+Rename(m))
    elif not m:
        motif_rep=(Rename(p)+'/None')
    else:
        motif_rep=(Rename(p)+'/'+Rename(m))
    yield motif_rep

def make_range(ls,p):
    if ls:
        last=ls[-1][1]
        if p-last<=1:
            ls[-1]=(ls[-1][0],p)
        else:
            ls.append((p,p+1))
        yield (ls)
    else:
        yield ([(p,p+1)])

def Same_Range(List1,List2):
    List1=sorted(List1, key=itemgetter(0))
    List2=sorted(List2, key=itemgetter(0))
    for S,E in List1:
        for S2,E2 in List2:
            if int(S2) <= int(S) <= int(E2):
                Start=int(S)
                if int(E) <= int(E2):
                    End=int(E)
                    yield Start,End
                else:
                    End=int(E2)
                    yield Start,End
            elif int(S2) <= int(E) <= int(E2):
                End=int(E)
                if int(S) <= int(S2):
                    Start=int(S2)
                    yield Start,End
                else:
                    Start=int(S)
                    yield Start,End

def Motif_count1(seq,Motifs,strand=0,motif_search=False):
    for i in Motifs:
        p=i[0];m=i[1] # p: plus strand motif / m: minus strand motif
        motif_rep=next(rep_name(p,m))
        if i[strand]: #filter None
            if motif_search:
                motif=i[strand].replace('_','',2)
            else:
                motif=i[strand]
            score=len(re.findall(
                      '(?={0})'.format(motif),seq,re.IGNORECASE))
        else:
            score=0
        yield motif_rep,score

def Motif_count2(seq,Motifs,strand=0,motif_search=False):
    tmp={}
    for i in Motifs:
        p=i[0];m=i[1] # p: plus strand motif / m: minus strand motif
        motif_rep=next(rep_name(p,m))
        if i[strand]: #filter None
            offset=i[strand].find('_')
            if motif_search:
                motif=i[strand].replace('_','',2)
            else:
                motif=i[strand]
            if strand==0:
                pos=set([j.start()+1+offset for j in re.finditer(
                   '(?={0})'.format(motif), seq, re.IGNORECASE)])
            elif strand==1:
                pos=set([len(seq)-j.start()-offset for j in re.finditer(
                   '(?={0})'.format(motif), seq, re.IGNORECASE)])
        tmp[next(rep_name(p,m))]=pos
    return tmp

def MAGs_contig1(fileP,MAGs,contigs,methyl,M_Motifs):
    genomeS=0
    fracS=[];fracL=0
    motifD=defaultdict(int)
    methylD=defaultdict(int)
    with open(os.path.join(fileP,MAGs),'r') as F:
        for title,seq in SFP(F):
            title=title.split()[0]
            if title in contigs:
                for s,e in contigs[title]:
                    # mapped fraction of MAG
                    if e-s <4:
                        continue
                    fracS.append(seq[s:e])
                    fracL+=len(seq[s:e])
            if title in methyl:
                for motif_rep in methyl[title]:
                    methylD[motif_rep]+=methyl[title][motif_rep]
            # MAG genome size
            genomeS+=len(seq)

    # Count motifs in frac
    fracS_rev=[str(Seq(i).reverse_complement()) for i in fracS]
    fracS_rev='-'.join(fracS_rev)
    fracS='-'.join(fracS)
    # Count motif
    for name,count in Motif_count1(fracS,M_Motifs,strand=0,motif_search=True):
        motifD[name]+=count
    # Count motif in reverse
    for name,count in Motif_count1(fracS_rev,M_Motifs,strand=1,motif_search=True):
        motifD[name]+=count

    cover=100.0*fracL/float(genomeS)
    valueD=defaultdict(str)
    for motif_rep in methylD:
        methylScore=methylD[motif_rep]
        motifScore=motifD[motif_rep]
        if motifScore>0:
            valueD[motif_rep]=str(100*methylScore/float(motifScore))
    return methylD,motifD,valueD,cover

def MAGs_contig2(fileP,M_Motifs):
    motifD=defaultdict(dict)
    with open(fileP,'r') as F:
        for title,seq in SFP(F):
            title=title.split()[0]
            seq_rev=Seq(seq).reverse_complement()
            # Count motif
            motifD[title]['+']=Motif_count2(seq,M_Motifs,strand=0,motif_search=True)
            # Count motif in reverse
            motifD[title]['-']=Motif_count2(str(seq_rev),M_Motifs,strand=1,motif_search=True)
        return motifD

def Get_Met(fileP,M_Motifs,dp):
    metD={}
    total_seq={}
    with open(fileP,'r') as M:
        while True:
           line=M.readline()
           if not line: break
           line=line.rstrip('\n').split()
           name=line[0]
           if not name.startswith('#'):
               if not name in metD:
                   metD[name]=defaultdict(int)
                   total_seq[name]=defaultdict(list)
               info=line[8]
               #filter read coverage !
               if int(info.split('coverage=')[1].split(';')[0])>=dp and int(line[5])>=30:
                   seq=info.split('context=')[1].split(';')[0]
                   if line[6]=='+':
                       strand=0
                   else:
                       strand=1
                   total_seq[name][strand].append(seq[:20]+'_'+seq[20]+'_'+seq[21:])

    for name in tq(total_seq,total=len(total_seq),desc='Methylation motif count', position=0):
        for strand in total_seq[name]:
            seq='-'.join(total_seq[name][strand])
            for met,count in Motif_count1(seq,M_Motifs,strand=strand):
                metD[name][met]+=count
    return metD

def Get_range(fileP,dp,Samtools,Bedtools):
    rangeD=defaultdict(list)
    cmd1=Samtools+" depth "+fileP
    cmd2=r"""awk -F '\t' '{{if(int($3)>={}) printf("%s\t%d\t%d\n",$1,int($2)-1,$2);}}'""".format(str(dp))
    cmd3=Bedtools+" merge -d 0"
    P1=subprocess.Popen(cmd1,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    P2=subprocess.Popen(cmd2,stdin=P1.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    P3=subprocess.Popen(cmd3,stdin=P2.stdout,stdout=subprocess.PIPE,shell=True)#.decode('utf-8')
    while True:
        line = P3.stdout.readline().decode('utf-8')
        if not line: break
        line=line.rstrip('\n').split()
        rangeD[line[0]].append((int(line[1]),int(line[2])))
    return rangeD

def nMotif_filter(Sample,MAG,Motif,motif_n):
    try:
        if Motif in motif_n[Sample][MAG]:
            return True
        else:
            return False
    except:
        return False

def get_Cov(COV_info_d,a,b):
    return COV_info_d[a][b]

def same_df(df,t_motif,t_MAGs,tree_order=False):
    for col in t_motif:
        if not col in df:
            df[col]=-50
    for row in t_MAGs:
        if not row in df.index:
            df.loc[row]=-50
    if tree_order:
        df=df.reindex(tree_order)
    df=df[t_motif]
    return df
