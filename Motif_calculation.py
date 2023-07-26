import os
import itertools
import optparse
import pathlib
import gzip
import pickle
import time
import pandas as pd
import utils.MEpi_fun as MEpi
from tqdm import tqdm as tq
from collections import defaultdict

#ipdSummary gff format nucleotid position is 1-index based.
def main(Genome, Bam, Ipd, Out, Motif, dp):
    # check DB and program path
    info_P=str(pathlib.Path(__file__).parent.absolute())+'/programs.txt'
    assert os.path.isfile(info_P), "Please check program.txt file"
    with open(info_P,'r') as I:
        for i in I:
            i=i.rstrip('\n').split(':')
            if i[0]=='bedtools':
                Bedtools=MEpi.find_p(i[1])
            if i[0]=='samtools':
                Samtools=i[1] #find_p(i[1])
    
    # make folder 
    if not os.path.exists(Out+'/RSII_pickle/'):
        os.makedirs(Out+'/RSII_pickle/')
    if not os.path.exists(Out+'/met_result/'):
        os.makedirs(Out+'/met_result/')

    # only use samples with co-exists in bam and ipdsummary folder
    Bams=set([x.split('.bam')[0] for x in filter(lambda x:x.endswith('.bam'),os.listdir(Bam))])
    Ipds=set([x.split('.gff')[0] for x in filter(lambda x:x.endswith('.gff'),os.listdir(Ipd))])
    Sample=list(Bams&Ipds)
    Sample.sort()
    assert len(Sample)!=0, 'there is no sample with ipdSummary and bam file'

    # get genomes files 
    fMAGs=os.listdir(Genome)
    fMAGs=list(filter(lambda x:x.endswith('.fa'),fMAGs))
    ext='.fa'
    if len(fMAGs)==0:
        fMAGs=list(filter(lambda x:x.endswith('.fna'),fMAGs))
        ext='.fna'
    assert len(fMAGs)!=0, 'there is no genome file (*.fa or *.fna)'
    fMAGs.sort()
    
    # get Motif with 5bp flanking 4mC or 6mA
    if Motif=='5bp':
        M_Motifs={}
        for i in map(''.join, itertools.product('GATC', repeat=4)):
            #j=str(Seq(i).reverse_complement())
            fA=''.join(i[:2])+'_A_'+''.join(i[2:])
            #rA=''.join(j[:2])+'_A_'+''.join(j[2:])
            fC=''.join(i[:2])+'_C_'+''.join(i[2:])
            #rC=''.join(j[:2])+'_C_'+''.join(j[2:])
            M_Motifs[(fA,fA)]=''
            M_Motifs[(fC,fC)]=''
    else:
        with gzip.open(Motif, 'rb') as f:
            M_Motifs = pickle.load(f)

    # calculate frequency of motifs on genomes with specific depth
    start_t=time.time()
    totalC=defaultdict(dict)
    print ('Sample lists: '+ ', '.join(Sample))
    for s in Sample:
        print (s)
        if os.path.isfile(Out+'/RSII_pickle/'+s+'.d'+str(dp)+'.pickle'):
            with open(Out+'/RSII_pickle/'+s+'.d'+str(dp)+'.pickle', 'rb') as f:
                ContigR = pickle.load(f)
        else:
            ContigR=MEpi.Get_range(Bam+'/'+s+'.bam',dp,Samtools,Bedtools)
            with open(Out+'/RSII_pickle/'+s+'.d'+str(dp)+'.pickle', 'wb') as f:
                pickle.dump(ContigR, f, pickle.HIGHEST_PROTOCOL)

        if os.path.isfile(Out+'/met_result/'+s+'.d'+str(dp)+'.pickle'):
            with open(Out+'/met_result/'+s+'.d'+str(dp)+'.pickle','rb') as f:
                metD = pickle.load(f)
        else:
            metD=MEpi.Get_Met(Ipd+'/'+s+'.gff',M_Motifs,dp)
            with open(Out+'/met_result/'+s+'.d'+str(dp)+'.pickle', 'wb') as f:
                pickle.dump(metD, f, pickle.HIGHEST_PROTOCOL)

        totalD=defaultdict(dict)
        totalM1=defaultdict(dict)
        totalM2=defaultdict(dict)
        num=len(fMAGs)
        for b in tq(fMAGs,total=num,desc='finding Motifs in MAGs'):
            totalM1[b.rstrip(ext)],totalM2[b.rstrip(ext)],totalD[b.rstrip(ext)],totalC[s][b.rstrip(ext)]=MEpi.MAGs_contig1(Genome,b,ContigR,metD,M_Motifs)

        df=pd.DataFrame(totalD)
        dfM1=pd.DataFrame(totalM1)
        dfM2=pd.DataFrame(totalM2)
        dfC=pd.DataFrame(totalC)
        df.T.to_csv(Out+'/met_result/'+s+'.d'+str(dp)+'.frac.csv',index_label='methylome')
        dfM1.T.to_csv(Out+'/met_result/'+s+'.d'+str(dp)+'.methyl.csv',index_label='methylome')
        dfM2.T.to_csv(Out+'/met_result/'+s+'.d'+str(dp)+'.motif.csv',index_label='methylome')
        dfC.T.to_csv(Out+'/met_result/'+s+'.d'+str(dp)+'.cov.csv',index_label='MAGs')

    df=pd.DataFrame(totalC)
    COV_info_d=df.to_dict()
    df.T.to_csv(Out+'/met_result/coverage.total.csv',index_label='MAGs')

    # get ID of motifs and MAGs 
    mm=os.listdir(Out+'/met_result/')
    mm=list(filter(lambda x:x.endswith(str(dp)+'.frac.csv'),mm))
    mm.sort()
    t_MAGs=set()
    t_motif=set()
    for i in mm:
        df=pd.read_csv(Out+'/met_result/'+i,header=0,index_col=0)
        # Motif filtering by methylated sum > 10
        df=df[df.columns[df.sum()>10]]
        # MAG filtering by methylated sum > 10
        df=df[df.sum(axis=1)>10]
        # add motifs and MAGs
        t_MAGs=t_MAGs|set(list(df.index))
        t_motif=t_motif|set(list(df))

    # filter with motif counts <- more than 10 #motifs in mapped region
    mm=os.listdir(Out+'/met_result/')
    mm=list(filter(lambda x:x.endswith(str(dp)+'.motif.csv'),mm))
    mm.sort()
    t_MAGs=list(t_MAGs)
    t_motif=list(t_motif)
    t_motif.sort()
    motif_n={}
    for i in mm:
        ID=i.split('.')[0]
        motif_n[ID]=defaultdict(set)
        df=pd.read_csv(Out+'/met_result/'+i,header=0,index_col=0)
        tmp=(df
         .reset_index()
         .melt(id_vars='methylome')
         .loc[lambda df:df.value>=10]
        )
        for index, row in tmp.iterrows():
            motif_n[ID][row['methylome']].add(row['variable'])

    # total methylome data frame from MAGs
    mm=os.listdir(Out+'/met_result/')
    mm=list(filter(lambda x:x.endswith(str(dp)+'.frac.csv'),mm))
    mm.sort()
    df_T=pd.DataFrame()
    for i in mm:
        df_t=pd.read_csv(Out+'/met_result/'+i,header=0,index_col=0)
        df_t=MEpi.same_df(df_t,t_motif,t_MAGs) # filtering/ ordering MAGs and Motifs #t_motif_strict or t_motif
        df_t['MAGs']=df_t.index
        df_t['Sample']=str(i.split('.')[0])
        df_t=df_t.set_index(df_t.index.astype(str) + '|'+str(i.split('.')[0]))
        df_T=pd.concat([df_T,df_t])
    df_T=df_T.fillna(-50)
    df_T['Coverage'] = df_T.apply(lambda x: MEpi.get_Cov(COV_info_d= COV_info_d, a = x['Sample'], b = x['MAGs']), axis=1)

    tmp=pd.melt(df_T, id_vars=['MAGs','Sample','Coverage'])
    tmp=tmp[tmp.apply(lambda x: MEpi.nMotif_filter(Sample = x['Sample'], MAG = x['MAGs'], Motif = x['variable'], motif_n=motif_n), axis=1)] # nMotif filter (n>=10)
    # mapped coverage filter
    df_F=(tmp[(tmp.Coverage>=20)]
     .pivot_table(columns='variable',values='value',index=['MAGs','Sample','Coverage'], fill_value=0)
     #.assign(Domain=lambda df:df.reset_index()['MAGs'].apply(get_domain).values)
    )
    df_F.to_csv(Out+'/met_result/methylome.dp'+str(dp)+'.total.csv',header=True, index=True) # mapped coverage and number of motifs

if __name__=="__main__":
    usage = """Motif_calculation.py -b [bam; folder of bamfiles; Extension: bam] -i [folder of ipdSummary files; Extension: gff] -g [folder of MAGs; Extension: fa or fna] -o [output]"""
    parser = optparse.OptionParser(usage)
    parser.add_option("-g","--genome",dest="genome",
        help="metagenome-assembled genome file path",type="string")
    parser.add_option("-o","--output",dest="output",
        help="output folder",type="string")
    parser.add_option("-b","--bam",dest="bam",
        help="bam file folder| extension: bam",type="string")
    parser.add_option("-i","--ipdSummary",dest="ipdSummary",
        help="ipdSummary file folder| extension: gff",type="string")
    parser.add_option("-m","--motif",dest="motif",
        help="MTase_motif.pickle from REBASE| default: use denovo 5bp NNANN, NNCNN",type="string")
    parser.add_option("-d","--depth",dest="depth",
        help="mapped read depth| default 10",type="string")
    (opts,args)=parser.parse_args()

    if opts.genome is None or opts.output is None or opts.bam is None or opts.ipdSummary is None:
        parser.print_help()
    else:
        if opts.motif is None:
            opts.motif='5bp'
        elif not os.path.exists(opts.motif):
            print ('Motif.pickle file not found... use the denovo motifs with 5bp (NNANN, NNCNN)')
            opts.motif='5bp'
        if opts.depth is None:
            opts.depth='10'
        main(opts.genome, opts.bam, opts.ipdSummary, opts.output, opts.motif, int(opts.depth))
