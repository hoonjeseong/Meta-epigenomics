#### __load library__


```python
import urllib.request
import re
import gzip
import pickle
import pandas as pd
from collections import defaultdict
from Bio.Seq import Seq
```

#### __output variable__


```python
mslist='./msublist.html'
mprot='./msubprolist.html'
mprot_seq='./msub_prot.faa'
mtase_pickle='./MTase.pickle'
mtase_df='./MTase_motifs.tsv'
```

#### __get rebase data__


```python
urllib.request.urlretrieve("http://rebase.neb.com/cgi-bin/msubprolist", mprot)
urllib.request.urlretrieve("http://rebase.neb.com/cgi-bin/msublist", mslist)
```




    ('./msublist.html', <http.client.HTTPMessage at 0x7f077d0b1110>)



#### __load functions__


```python
def cleanhtml(raw_html):
    cleanr = re.compile('<.*?>')
    cleantext = re.sub(cleanr, '', raw_html)
    return cleantext

def Change_Form(Seq):
    Seq=re.sub('R','[A,G]',Seq)
    Seq=re.sub('Y','[C,T]',Seq)
    Seq=re.sub('S','[G,C]',Seq)
    Seq=re.sub('W','[A,T]',Seq)
    Seq=re.sub('K','[G,T]',Seq)
    Seq=re.sub('M','[A,C]',Seq)
    Seq=re.sub('B','[C,G,T]',Seq)
    Seq=re.sub('D','[A,G,T]',Seq)
    Seq=re.sub('H','[A,C,T]',Seq)
    Seq=re.sub('V','[A,C,G]',Seq)
    Seq=re.sub('N','.',Seq)
    yield Seq

def check_palindrom(seq):
    my_dna = Seq(seq)
    rv_seq=str(my_dna.reverse_complement())
    if seq==rv_seq:
        return True

def Rename(Motif):
    yield Motif.replace('.','N').replace('[A,T]','W').replace('[A,G]','R')\
    .replace('[C,T]','Y').replace('[G,C]','S').replace('[G,T]','K')\
    .replace('[A,C]','M').replace('[A,C,T]','H').replace('[A,C,G]','V')\
    .replace('[A,G,T]','D').replace('[C,G,T]','B')

def make_motif(seq,num1=None,num2=None):
    if num1 and num2:
        my_dna = Seq(seq)
        rv_seq=str(my_dna.reverse_complement())
        num2=len(seq)-num2
        seq1=seq[:num1-1]+'_'+seq[num1-1]+'_'+seq[num1:]
        seq2=rv_seq[:num2]+'_'+rv_seq[num2]+'_'+rv_seq[num2+1:]
        #print (seq,num1,num2,next(Change_Form(seq1)),next(Change_Form(seq2)))
        yield (next(Change_Form(seq1)),next(Change_Form(seq2)))
    elif num1:
        seq1=seq[:num1-1]+'_'+seq[num1-1]+'_'+seq[num1:]
        seq2=None
        yield (next(Change_Form(seq1)),None)
    elif num2:
        num2=len(seq)-num2
        my_dna = Seq(seq)
        rv_seq=str(my_dna.reverse_complement())
        seq1=None
        seq2=rv_seq[:num2]+'_'+rv_seq[num2]+'_'+rv_seq[num2+1:]
        yield (None,next(Change_Form(seq2)))
```

#### __make faa file and save Mtase informantion__


```python
Mtase=defaultdict(dict)
flag=0
temp=[]
with open(mslist,'r') as M:
    for line in M:
        i=cleanhtml(line).rstrip('\n')
        if line.startswith("</tr>"):
            flag+=1
            if flag>1:# and len(temp[1:])==4:
                if '6-methyladenosine' in temp[-1]:
                    Mtase[temp[0]]=temp[1:]
                elif '4-methylcytosine' in temp[-1]:
                    Mtase[temp[0]]=temp[1:]
                elif '5-methylcytosine' in temp[-1]:
                    Mtase[temp[0]]=temp[1:]
            temp=[]
        elif not i:
            continue
        else:
            temp.append(i.strip())

IDs=set()
with open(mprot,'r') as M, open(mprot_seq,'w') as F:
    for line in M:
        i=cleanhtml(line).rstrip('\n')
        if line.startswith("</tr>"):
            flag+=1
            if flag>1:
                if temp[0] in Mtase:
                    seq=''.join(temp[4:]).replace(' ','')
                    seq='\n'.join([seq[start:start+80] for start in range(0, len(seq), 80)])
                    if not seq=='-':
                        F.write('>'+temp[0]+'\n'+seq+'\n')
                        IDs.add(temp[0])
                    else:
                        del Mtase[temp[0]]
            temp=[]
        elif not i:
            continue
        else:
            temp.append(i.strip())

for i in (set(Mtase.keys())-IDs):
    del Mtase[i]

print (Mtase['M.Sdy9764Dam'])
    
# save total Mtase data and compress.
with gzip.open(mtase_pickle, 'wb') as f:
    pickle.dump(Mtase, f)
```

    ['Orphan  methyltransferase', 'a', 'GATC', 'Base (Type of modification):2 (6-methyladenosine)']


#### __Save Mtase information by methylated motifs__


```python
Mtase_Type=defaultdict(dict)
p=re.compile(r"\(?\d+[-]methyl[^\d]*sine\)")
for i in Mtase:
    # filter Mtase info which have modification motif seq and with only 4 columns
    if not '?' in Mtase[i][-1] and not '-' in Mtase[i][1] and len(Mtase[i])==4:
        line=Mtase[i][-1].split(':')[1]
        if '(Complementary strand)' in line:
            if len(p.findall(line)) ==1:
                pos1=None
                pos2=int(line.split('(Complementary strand)')[1].split()[0])
            else:
                pos1=int(line.split()[0])
                pos2=int(line.split('(Complementary strand)')[1].split()[0])
            seq=next(make_motif(Mtase[i][2],pos1,pos2))
            #print (Mtase[i])
        else:
            pos=int(line.split()[0])
            if check_palindrom(Mtase[i][2]):
                seq=next(make_motif(Mtase[i][2],pos,len(Mtase[i][2])-pos+1))
            else:
                seq=next(make_motif(Mtase[i][2],pos))
        if 'Type' not in Mtase_Type[seq]:
            Mtase_Type[seq]['Type']=set()
        if 'geneID' not in Mtase_Type[seq]:
            Mtase_Type[seq]['geneID']=set()
        # type name change
        Type=Mtase[i][0].replace('Putative','').replace('methyltransferase','').strip()
        if 'restriction' in Type:
            Type=Type.split('restriction')[0].strip()
        Mtase_Type[seq]['Type'].add(Type.replace(' ','_'))
        Mtase_Type[seq]['geneID'].add(i)

df=pd.DataFrame(Mtase_Type).T
new_id=[]
for a,b in zip(df.index.get_level_values(0),df.index.get_level_values(1)):
    if pd.notnull(a):
        A=next(Rename(a))
    else:
        A='None'
    if pd.notnull(b):
        B=next(Rename(b))
    else:
        B='None'
    new_id.append((A,B))
df.index=new_id
df.Type=df.Type.str.join(', ')
df.geneID=df.geneID.str.join(', ')
df.index.names = ['Methylated motifs (Forward,Reverse)']
df.to_csv(mtase_df,sep='\t')
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Type</th>
      <th>geneID</th>
    </tr>
    <tr>
      <th>Methylated motifs (Forward,Reverse)</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>(G_A_TC, G_A_TC)</th>
      <td>Type_II, Orphan</td>
      <td>M.Sdy9764Dam, M.Eco12650Dam, M.RaqMEM40Dam, M....</td>
    </tr>
    <tr>
      <th>(_A_GCNNNNRTCA, TG_A_YNNNNGCT)</th>
      <td>Type_I</td>
      <td>M.AalSMS7II</td>
    </tr>
    <tr>
      <th>(G_A_GNNNNNNNGTG, C_A_CNNNNNNNCTC)</th>
      <td>Type_I</td>
      <td>M.AalSMS7III</td>
    </tr>
    <tr>
      <th>(ACCG_A_G, None)</th>
      <td>Type_III, Type_IIG</td>
      <td>Aam10684II, M.SspLM7II</td>
    </tr>
    <tr>
      <th>(G_A_NTC, G_A_NTC)</th>
      <td>Type_II, Orphan</td>
      <td>M.PinP74I, M.HpyAs005IX, M.PgaP129I, M.OpsK8I,...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>(GC_A_NNNNNNCTGG, CC_A_GNNNNNNTGC)</th>
      <td>Type_I</td>
      <td>M.Yre11967I, M.Yre11968I, M.Yre11966I</td>
    </tr>
    <tr>
      <th>(TG_A_CNNNNNNGTC, G_A_CNNNNNNGTCA)</th>
      <td>Type_I</td>
      <td>M.Yre11967II, M.Yre11966II</td>
    </tr>
    <tr>
      <th>(G_A_CNNNNNNGTCA, TG_A_CNNNNNNGTC)</th>
      <td>Type_I</td>
      <td>M.Yre11968II</td>
    </tr>
    <tr>
      <th>(GGC_A_GG, None)</th>
      <td>Type_III</td>
      <td>M.YroYRAII</td>
    </tr>
    <tr>
      <th>(G_A_GNNNNNNNTTCC, GGA_A_NNNNNNNCTC)</th>
      <td>Type_I</td>
      <td>M.ZalSM2I</td>
    </tr>
  </tbody>
</table>
<p>2697 rows × 2 columns</p>
</div>


