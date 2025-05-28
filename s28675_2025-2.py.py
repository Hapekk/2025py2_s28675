#!/usr/bin/env python3
import sys
from Bio import Entrez,SeqIO
import pandas as pd,matplotlib.pyplot as plt

Entrez.email,Entrez.api_key = "s28675@pjwstk.edu.pl","0c603b5f507827ac323147e96a5485f91f08"

taxid,lo,hi,mx = (sys.argv[1:5] if len(sys.argv)==5 else input("Podaj taxid min max max_rec: ").split())
lo,hi,mx = map(int,(lo,hi,mx))

r = Entrez.read(Entrez.esearch(db="nucleotide",term=f"txid{taxid}[Organism]",usehistory="y"))
cnt,we,qn = int(r["Count"]),r["WebEnv"],r["QueryKey"]

rows=[]
for i in range(0,min(cnt,mx),400):
    h=Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retstart=i,retmax=min(400,mx-i),webenv=we,query_key=qn)
    for x in SeqIO.parse(h,"genbank"):
        L=len(x.seq)
        if lo<=L<=hi: rows.append((x.id,L,x.description))

df = pd.DataFrame(rows,columns=["Acc","Len","Desc"]).sort_values(by="Len",ascending=False)
csvf = f"gb_{taxid}.csv"; df.to_csv(csvf,index=False)

plt.plot(df["Acc"],df["Len"],marker='o')
plt.xticks(rotation=45,fontsize=6);plt.tight_layout()
plt.savefig(f"gb_{taxid}.png")

print(f"Zapisano: {csvf}, wierszy: {len(df)}")
