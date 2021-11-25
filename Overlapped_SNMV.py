import os
import optparse
import pickle
import copy
import gzip as gz
import pandas as pd
import utils.MEpi_fun as MEpi
from tqdm import tqdm as tq
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as SFP

def main(Path, Genome_path, Dp, Score, Search_motif, Output):
    # get genomes files 
    fMAGs=os.listdir(Genome_path)
    fMAGs=list(filter(lambda x:x.endswith('.fa'),fMAGs))
    ext='.fa'
    if len(fMAGs)==0:
        fMAGs=list(filter(lambda x:x.endswith('.fna'),fMAGs))
        ext='.fna'
    assert len(fMAGs)!=0, 'there is no genome file (*.fa or *.fna)'
    fMAGs.sort()
    # get contig info
    contigs=set()
    MAG_ctg={}
    for i in fMAGs:
        with open(Genome_path+'/'+i,'r') as F:
            for t,s in SFP(F):
                header=t.split()[0]
                assert not header in contigs, 'there are duplicate contigs name. please make an unique contig name beteween MAGs'
                contigs.add(header)
                MAG_ctg[header]=i.split(ext)[0]

    # get IPDsummary files 
    tmp=os.listdir(Path)
    tmp1=list(filter(lambda x:x.endswith('.csv'),tmp))
    tmp2=list(filter(lambda x:x.endswith('.gff'),tmp))
    tmp3=set()
    fIPDs=[]
    for i in tmp1:
        with open(Path+'/'+i,'r') as I:
            if next(I)=='refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage,frac,fracLow,fracUp\n':
                tmp3.add(i.split('.csv')[0])
    for i in tmp2:
        with open(Path+'/'+i,'r') as I:
            next(I)
            if 'ipdSummary' in next(I):
                 sample=i.split('.gff')[0]
                 if sample in tmp3:
                     fIPDs.append(sample)
    assert len(fIPDs)!=0, 'there is no ipdsummary results (csv/gff)'
    fIPDs.sort()

    Depth_total={}
    for i in fIPDs:
        n_line=sum(1 for i in open(Path+'/'+i+'.csv'))
        Depth=defaultdict(dict)
        print ('loading a sample ({0}) mapped region'.format(i))
        flag=0
        with open(Path+i+'.csv','r') as D:
            n_line-=1
            next(D)
            for line in tq(D,total=n_line):
            #for line in D:
                n_line-=1
                line=line.rstrip('\n').split(',')
                pos=int(line[1])
                depth=int(line[9])
                if depth < Dp:
                    continue
                ctg=line[0].replace('"','')
                if not ctg in contigs: #filter of indifference contigs
                    continue
                if line[2]=='0':
                    strand='+'
                else:
                    strand='-'
                if not ctg in Depth:
                    Depth[ctg]['+']=[]
                    Depth[ctg]['-']=[]
                Depth[ctg][strand]=next(MEpi.make_range(Depth[ctg][strand],pos))
        Depth_total[i]=Depth
    print ('calculated overlapped region between {0} samples'.format(len(Depth_total)))
    remove_ctg=set()
    for n,i in tq(enumerate(Depth_total),total=len(Depth_total)):
        # Overlapped region
        if n==0:
            Depth_ov=copy.deepcopy(Depth_total[i]) # make overlapped region
        else:
            ov_contig=set(Depth_total[i].keys())&set(Depth_ov.keys())-remove_ctg
            for c in ov_contig:
                ov1=MEpi.Same_Range(list(Depth_ov[c]['+']),list(Depth_total[i][c]['+']))
                ov2=MEpi.Same_Range(list(Depth_total[i][c]['+']),list(Depth_ov[c]['+']))
                Depth_ov[c]['+']=set()
                for r in ov1:
                    Depth_ov[c]['+'].add(r)
                for r in ov2:
                    Depth_ov[c]['+'].add(r)
                ov1=MEpi.Same_Range(list(Depth_ov[c]['-']),list(Depth_total[i][c]['-']))
                ov2=MEpi.Same_Range(list(Depth_total[i][c]['-']),list(Depth_ov[c]['-']))
                Depth_ov[c]['-']=set()
                for r in ov1:
                    Depth_ov[c]['-'].add(r)
                for r in ov2:
                    Depth_ov[c]['-'].add(r)
        # remove non union contigs
        for c in Depth_ov:
            if not c in Depth_total[i]:
                remove_ctg.add(c)
    for c in remove_ctg:
        del(Depth_ov[c])
    #with gz.open(Ov_pickle,'wb') as D:
    #    pickle.dump(Depth_ov,D)
    assert len(Depth_ov) != 0, 'there is no overlapped region between samples'

    # get methylated fraction dataframe
    metD=defaultdict(dict)
    Fwd_re=MEpi.Rename2(Search_motif.split('/')[0])
    Rev_re=MEpi.Rename2(Search_motif.split('/')[1])
    for i in fIPDs:
        with open(Path+i+'.gff','r') as M:
            while True:
                line=M.readline()
                if not line: break
                line=line.rstrip('\n').split()
                name=line[0]
                if not name in Depth_ov:
                    continue
#                if not name.startswith('#'):
                info=line[8]
                if int(info.split('coverage=')[1].split(';')[0])>=Dp and line[2]!='modified_base' and int(line[5])>=Score: #filter option
                    Pos=int(line[3])
                    if line[6]=='+':
                        Search_re=Fwd_re
                    else:
                        Search_re=Rev_re
                    for s,e in Depth_ov[name][line[6]]:
                        if s<=Pos<=e:
                            try: 
                                if next(MEpi.grep_pattern(info.split('context=')[1].split(';')[0],Search_re)):
                                    metD[i][name+':'+str(Pos)+':'+line[6]]=float(info.split(';frac=')[1].split(';')[0]) # get fraction
                                    #metD[i][Pos]=float(info.split('IPDRatio=')[1].split(';')[0])*flag*Order[name][1] # get IPDratio
                            except:
                                continue
    df=pd.DataFrame.from_dict(metD)
    df=df.fillna(0).T
    assert len(df.columns) !=0 , 'there is modification signal of the motif'
    # Add non methylated position
    for i in fMAGs:
        motifs_genome=MEpi.MAGs_contig2(Genome_path+'/'+i,[(Fwd_re,Rev_re)])
    
    for i in motifs_genome:
        for j in motifs_genome[i]['+'][Search_motif]:
            for s,e in Depth_ov[i]['+']:
                if j>=s and j<=e:
                    if i+':'+str(j)+':+' not in df.columns:
                        df[i+':'+str(j)+':+']=0
                    break
        for j in motifs_genome[i]['-'][Search_motif]:
            for s,e in Depth_ov[i]['-']:
                if j>=s and j<=e:
                    if i+':'+str(j)+':-' not in df.columns:
                        df[i+':'+str(j)+':-']=0
                    break
    # add genome name into column
    multi_col=[]
    for i in df.columns:
        multi_col.append((MAG_ctg[i.split(':')[0]],i))
    df.columns = pd.MultiIndex.from_tuples(multi_col)
    # save target genome methylome
    df.to_csv(Output, index=True, header=True) #header=True, index=True, index_label='MAGs')

if __name__=="__main__":
    usage = """Overlapped_SNMV.py -i [folder of ipdSummary files; Extension: csv and gff] -g [folder of MAGs; Extension: fa or fna] -m motif [fwd/rev] -o [output]"""
    parser = optparse.OptionParser(usage)
    parser.add_option("-g","--genome",dest="genome",
        help="metagenome-assembled genome file path",type="string")
    parser.add_option("-o","--output",dest="output",
        help="output file name",type="string")
    parser.add_option("-i","--ipdSummary",dest="ipdSummary",
        help="ipdSummary file (gff and csv files) folder",type="string")
    parser.add_option("-m","--motif",dest="motif",
        help="Motifs for search by each strand | e.g. 'G_A_NTC/G_A_NTC'(fwd/rev); Modificated base should be indicated by an underscore.",type="string")
    parser.add_option("-d","--depth",dest="depth",
        help="mapped read depth| default 20",type="string")
    parser.add_option("-s","--score",dest="score",
        help="modified base score by ipdSummary| default 30",type="string")
    (opts,args)=parser.parse_args()

    if opts.genome is None or opts.output is None or opts.motif is None or opts.ipdSummary is None:
        parser.print_help()
    else:
        if opts.depth is None:
            opts.depth=20
        if opts.score is None:
            opts.score=30
        main(opts.ipdSummary, opts.genome, int(opts.depth), int(opts.score), opts.motif, opts.output)
