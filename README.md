# Adv_Python_Portfolio
This is William Payne Smith's portfolio of python code that I learned during BISC 450C (ADV Python)


## SEQUENCE OBJECTS 1-4
```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence

print(len(my_seq))
```

    5



```python
#We can print specific nucleotides of a sequence

print (my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
Seq("AAAA").count('AA')
```




    2




```python
my_seq = Seq("GATCAGATGGGCCTATAAGGATCGAAAATCGC")
```


```python
#We can determine the number of nucleotides in sequence

len(my_seq)
```




    32




```python
my_seq.count("G")
```




    9




```python
#We can calculate the G-C content of a sequence

100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    46.875




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("CATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
#We can use an imported equation to calculate G-C content

gc_fraction(my_seq)
```




    0.46875




```python
my_seq[4:12]
```




    Seq('GATGGGCC')




```python
#We can isolate specific intervals of nucleotides, necessary to choose where we start to read codons.

my_seq[0::3]
```




    Seq('CCTGTAGTAAG')




```python
my_seq[2:3]
```




    Seq('T')




```python
#We can print the sequence backwards

my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAC')




```python
str(my_seq)
```




    'CATCGATGGGCCTATATAGGATCGAAAATCGC'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
#We can insert a verb placeholder to indicate what the fasta string is

print(fasta_format_string)
```

    >Name
    CATCGATGGGCCTATATAGGATCGAAAATCGC
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
#We can use this spacer object made above and join it with "contigs".
#This places the ten Ns between each contig sequence.

spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTTGCA')




```python
dna_seq = Seq("acgtACGT")


```


```python
dna_seq
```




    Seq('acgtACGT')




```python
#We can use this function to make sure all of your sequence is uppercase (or lowercase: "dna_seq.lower()")

dna_seq.upper()
```




    Seq('ACGTACGT')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
#We can print the complement sequence of our sequence original sequence

my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
#We can print the reverse complement of our sequence.

my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
#We can identify a protein sequence, then determine its complement and print it.
#The letters are ambiguity codes
#V means A C or G, N mean any, etc

protein_seq = Seq("EVRANAK")
protein_seq.complement()
```




    Seq('EBYTNTM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAGGGTGCCCGATAG')




```python
# We use the reverse complement function to generate our template strand based on coding_dna

template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAGGGTGCCCGATAG')




```python
# We use the transcribe function to replace "T" with "U"; transcribe coding strand to mRNA
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAGGGUGCCCGAUAG')




```python
# We use this function to see the true transcribtion that takes place in biology, template strand to mRNA.

template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAGGGUGCCCGAUAG')




```python
# We can reverse the mRNA back to the original sequence

messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAGGGTGCCCGATAG')




```python
# Transcription DNA to RNA (T to U)
# RNA to ribosome - ribosome takes codons to identify amino acids; Translation
# This function translates our RNA sequence into a protein
# "*" means stop codon, in this case premature
messenger_rna.translate()

```




    Seq('MAIVMGR*RVPD')




```python
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRW*VPD')




```python
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRW*VPD')




```python
# We can make the protein output stop at the prescence of the stop codon '*'

coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
# We can make the protein output stop at the prescence of the stop codon '*'

coding_dna.translate(table= 2, to_stop=True)
```




    Seq('MAIVMGRW')




```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRW!VPD')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGTCCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLSPWQHRLRKLR*SRQ*NYR*AIVIIVAITGMEVTGATTA...SPL')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLSPWQHRLRKLR')




```python
gene.translate(table = "Bacterial", cds = True)
```


    ---------------------------------------------------------------------------

    TranslationError                          Traceback (most recent call last)

    <ipython-input-153-9b7bfe3f3060> in <module>
    ----> 1 gene.translate(table = "Bacterial", cds = True)
    

    ~/anaconda3/lib/python3.7/site-packages/Bio/Seq.py in translate(self, table, stop_symbol, to_stop, cds, gap)
       1446 
       1447         return self.__class__(
    -> 1448             _translate_str(str(self), table, stop_symbol, to_stop, cds, gap=gap)
       1449         )
       1450 


    ~/anaconda3/lib/python3.7/site-packages/Bio/Seq.py in _translate_str(sequence, table, stop_symbol, to_stop, cds, pos_stop, gap)
       2791         if n % 3 != 0:
       2792             raise CodonTable.TranslationError(
    -> 2793                 f"Sequence length {n} is not a multiple of three"
       2794             )
       2795         if str(sequence[-3:]).upper() not in stop_codons:


    TranslationError: Sequence length 292 is not a multiple of three



```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
# Vertebrate Mitochondrial is table 2

mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAACCTGAATGTGAGATCAGTCAAGGATAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CTGAATGTGA')




```python
seq[117512670:]
```




    Seq({13: 'TTGAAACCTGAATGTGAGATCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length = 10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAGGTGCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAGGTGCCGA')




```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAGGTGCCGA')




```python
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAGGTGCCGA')




```python
mutable_seq.reverse
```




    <bound method MutableSeq.reverse of MutableSeq('GCCACGTAATGGGCCGCTGAAGGTGCCGA')>




```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAGGTGCCGA')




```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```


```python
'CTAACCAGCAGCACGACCACCTTCCAACGACCCATAACAGC'
```




    'CTAACCAGCAGCACGACCACCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```

```python
'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUAG'
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUAG'




```python
back_transcribe(my_string)
```



```python
'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```




```python
'AVMGRWKGGRAAG*'
```


## SEQUENCE ANNOTATIONS 1-4

```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequnece for the VDB Computational Biology Class"
```


```python
print(simple_seq_r.description)
```

    Made up sequnece for the VDB Computational Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequnece for the VDB Computational Biology Class', dbxrefs=[])




```python
# Above gives us a way to store seqences with annotation: "name", "id", "description"
```


```python
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequnece for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# description prints the organism name
```


```python
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
```


```python
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python
# ^ End position of the protein feature is somewhere between position 5 and 9
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
from Bio.SeqRecord import SeqRecord
```


```python
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGSETHLDSMVGQALFGD"
...          "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLDFVPGLISK"
...          "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
...          "SSAC"),
                                id="gi|14150838|gb|AAK54648.1|AF376133_1" ,
...     description="chalcone synthase [Cucumis sativus]",
... )
```


```python
print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGSETHLDSMVGQ...SAC')



```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# We can bring up specific details or "features about a locations that we input: "condon_start", "type", etc
```


```python
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
len(sub_record.features)
```




    2




```python
# We can identify how many features within a subunit of a record
```


```python
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
sub_record.annotations["topology"] = "linear"
```


```python
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_record.id
```




    'NC_005816.1'




```python
sub_record.name
```




    'NC_005816'




```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
sub_record.description =  'Yersinia pestus buivar Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python
sub_record.description
```




    'Yersinia pestus buivar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestus buivar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...



```python
# We can bring up genbank records about our sequence: "bp", "accesion", and more
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
len(shifted)
```




    9609




```python
len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
# We can print the first value as a string, then present the next 4 as integers.
```


```python
rc = record.reverse_complement(id = "Testing")
```


```python
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



```python
# We can establish the id name to avoid sequence confusion.
```


## SEQUENCE I/O

```python
from Bio import SeqIO
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```


```python
for seq_record in Seq_IO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```



```python
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
# We can pull out all identifiers for every seq record, which could be describing species name in 
# practical example.
```


```python
record_iterator = SeqIO.parse("ls_orchid.gbk.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


    ---------------------------------------------------------------------------

    StopIteration                             Traceback (most recent call last)

    <ipython-input-8-898c3713621c> in <module>
    ----> 1 first_record = next(record_iterator)
    

    ~/anaconda3/lib/python3.7/site-packages/Bio/SeqIO/Interfaces.py in __next__(self)
         70         """Return the next entry."""
         71         try:
    ---> 72             return next(self.records)
         73         except Exception:
         74             if self.should_close_stream:


    StopIteration: 



```python
print(first_record.id)
```


```python
print(first_record.description)
```



```python
second_record = next(record_iterator)
```


    ---------------------------------------------------------------------------

    StopIteration                             Traceback (most recent call last)

    <ipython-input-11-40ef276aaa04> in <module>
    ----> 1 second_record = next(record_iterator)
    

    ~/anaconda3/lib/python3.7/site-packages/Bio/SeqIO/Interfaces.py in __next__(self)
         70         """Return the next entry."""
         71         try:
    ---> 72             return next(self.records)
         73         except Exception:
         74             if self.should_close_stream:


    StopIteration: 



```python
print(second_record.id)
```



```python
print(second_record.description)
```



```python
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
print("The last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records[6]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730


from Bio import SeqIO


```python
from Bio import SeqIO
```


```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
all_species= []
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank" )
]
```


```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
print(all_species)
```

    []



```python
record_iterator = SeqIO.parse("ls_orchid.gbk.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


    ---------------------------------------------------------------------------

    StopIteration                             Traceback (most recent call last)

    <ipython-input-35-898c3713621c> in <module>
    ----> 1 first_record = next(record_iterator)
    

    ~/anaconda3/lib/python3.7/site-packages/Bio/SeqIO/Interfaces.py in __next__(self)
         70         """Return the next entry."""
         71         try:
    ---> 72             return next(self.records)
         73         except Exception:
         74             if self.should_close_stream:


    StopIteration: 



```python
first_record.id
```




    'Z78533.1'




```python
first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='Z78533', description='C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
first_record.description = first_record.id + " " + "Mutations induced randomly"
```


```python
print(first_record.format("fasta"[:200]))
```

    >new_id Mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
# We can make changes in the ID to help us identify strands
```


```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord(
Seq("LDKFJGNASDLKJGEASDLKGJDNSGLKDJSNLDKLASKJDGLKAALDFKJHIUHADKJFU"
    "AKJHDSJHOUHASDLKVJBEIUHGFBVLAIOUABSDFOIUBADFUVNHAPIUHWEOIDSFF"
   "SLAISDUBNLIUEWAIOUBVPIUBWENAERGOIAIUDVBNIUBASDIUBNPIWEUNPIUSED"
   "SSAC",
),
    id = "gi | 14150838 |gb| AAK54658.1|",
    description = "chalcone synthase [Cumumis sativus]",
)
```


```python
rec2 = SeqRecord(
Seq("JUWEHYBFIWBNOWIUNFOUNEUIUCNWOEIUNCOWNOUBOUNOUIIUWSBNCIUNBF"
    "IWJDUENVIKJWEBNFIKJENCOUNWIUBOXIQWEUN",
   ),
    id = "gi | 13919613 |gb| AAK33142.1|",
    description = "chalcone synthase [Fragaria vesca subsp. bracteata]",
)
```


```python
rec3 = SeqRecord(
Seq("ALKJNIUBIUABNIDJUVNIAUBEROUABJUBIJABWEIKFJBEIKJBNJBA"
   "IKAWNEOVINOINQOINOPQISMCOIJPEWIKNLOKNLOIKJCLOIWJAWKLOKSKWJES"
   "KWJUESNOLIKNAIUWEBNOINIVUJWOEFUJOIUHJNAOWIEHFOPIKNHGOANSROLK"
    "AKJWENOINPIJUWBSEFIJUOINWQOIBCOINOINWOSBAONOINBOIBNAOWPSIKN"
    "AOENOLIKNPOIAUBIOUBNOIUNAOUEBIOUBOOUHOUOISCOIUNEVOBOUBIOUBWE"
   "AUJWENFFOUNAOUWBIKJUBAOUAWEIUJFNWKJESFBNIKABJNOJUBNWOESFINWES"
   "OAEIKRNGLOIKNOLIKEOGINW",
   ),
    id = "gi|13925890|gb|AAK49457.1|",
    description = "chalcone synthase [Nicotiana tabacum]"
)
```


```python
my_records = [rec1, rec2, rec3]
```


```python
from Bio import SeqIO
SeqIO.write(my_records, "my_example.faa", "fasta")
```




    3




```python
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
count = SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
# We can save our sequences from one file as a fasta file
```


```python
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement())
```

    Z78533.1
    GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGTCAACCGCCCTCCATTGGGGTTCGGAAGGGTTCTCCTGCCTGTCCGACAAGCACGACAACATGGGGGATTCGCACGGCAGCTGCTGCCAGCATCCGTCCACCTCTTGCCGGGTTCCGGGCCATCAAAACACCAGCTCTTGGACCCGCCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACCACGCCGGCCTGGCTGTATGCCGGGCAAGCATTGGCAGGAGAGAGACGAAGCGCGACGCCCAAGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCAAGATTCCCGTTTGCGAGAGTCATCAAAATTCATTGGGCACGCGACAGCACGCCGCCGCTCCGGGTTTTGGGGAAGACAATGCCATTCGCCGGTGATGCTTTCATATGGCTTGGCGCCCAAACTGCGCCGGGCTAGAGGTTCAAACCCGCCATGGACGCTCCCGAGGCGGCCCAACAACAAATCAGGGTCACCACGGGAGCAATGCCCCCGGTGAGCTGAGTACACCGGTCCTCCGGATTCACTCGATCGTTTATTCCACGGTCTCATCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78532.1
    GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAATCAACTGCCCAATGGTTATTTTAGCTCCATTGGGGTTCAATTAGGTTCTTGTGTAGGTTCGAAAAAATACAACAACATGGGGGATTCAAATAGCAGCCTTATGACTGTTAGCATTCTCCACCTCGTGCCACATTCCTACCCATCAAAGCAACAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCATTCACATCCGTATAATGCCAGCTTAGCGATATGCCAAGCAAGCATTGGTAGGAGAGACGCAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTAGCTGCGTTCTTCATCGATGTTAGAGCCAAGATATCCGTTGCGAGAGTCATAAAAATTCACTGGCATGCTCAGTAGCATACTGCCCCTCTGATTTTTTCTGACAATAATGTCATTCATCAGTGATGCTTTATATGACTTGGCGCAAACTGCACCGTACTAAAGTTCAAACCTGCCATGAAAGCTCTTGAGGAGGCCCAACAACAAAGCAGGGTCACGACAAAAGCAGTGCCACGACGAGCTGAGTTACCACAGGTCCTCCAGATTCACTCGATCATATATTCTGTTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78531.1
    TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCCCGAGGGCGGTTTTGCGTCCACTGGGGGTTCAAACAGGTTCTTCTGTAGGCCTGACGAGTACGACAACACGGGGGATTCGAACGGCAGCCTTGCGGCTGCCAGCATCCTCCACCTACTGCCAAGTTCCGGGCCATCAGAACACCGATGCTTAGACCCGCCGCACCTAGGCGCAAGGGGCCAATCTTTCACGTCCTCACAATGACAGACTAGCCGCATGCCAATCAAGCATTATCAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAGATATCCCGTTGCCGAGAGTCGTCATAAATTCATTGGCACGCGCAACAGACGCCGCCCCTCCGAGTTTTTCTTGACAAAAATGCCATCCATCGGTGACGCTCCATATGACTTGGCGCAAACTGCGCCGCGCTAGACGTTCAAACCCGCCATGAAAGTTCCCGAGGCGGCCCGGTCGCAAACCGGGTTCACCACGAGAGCAAAGCCACGGTGAGCCGTGTAACCACGGGTCCTCCGGATTCACTCGATCGTATGTTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78530.1
    ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAACATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGCGTTCTTCATCGATGCAAGAGTCAAGATATCCGTCTGCCGAGAGTCATCAAAATTCATTGGAATGAACATTAACACACCACACCTCCGACGTTGTGTGGGGCAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGGACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGATATCCCTAGCGAGCCAAATTACCACAAGTCCTCCAGATTCACTCAATCGTTTATTATGTTGTTTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78529.1
    TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAATATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGGGATTCTTCATCGATGCAAGAGCCAAGATATCCCTCGGCGAGAGTCATCAACAATTCATTGGCATGACCAATAACACACCACACCTCCGACTTTTTTGACAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGCACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGAAATCCCTAGCGAGCCAAATAACCACAAGTCCTCCAGATTCACTCAATCGTATATTCTGCTGTCTCAACAATGTCCTTCGGCAGCTCGCCGT
    Z78527.1
    GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACAATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACCCAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAAGAGCCAGATATCCGTTGGCCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAACACACTGCCACTCCAACTTTTGGTGACAGTAGTGCCATTTGCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGCCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGCGAGCCGAGTAACCACAGGTCATCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78526.1
    ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGGCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAGCACACTGCCACTCCAACTTTTGGTGACAATAATGCCATTGTTCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGGCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGGGAGCCGAGTAACCACAGGTCCTCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78525.1
    TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATCCAAGGGTCATTTTGCCTCCACTGGGGTTCAAACAAGTTATTATATAGGTCCGATAAATACGATAACATGGGGGATTTAAATGACAACCTTGTGGCAGCCAATGTCCTCCACCTCCTGCCAAGTTTCGGGCCATCAAAACATTCCTTAGACCCAACGCGCCAAGGCACAAGGGGCCAATCTTTCGCATCCGCACAATGCAAGCCTAGATATATGGACAGGCATTGGCAAGAGAGACACAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCCGATGCAAGAGCCAAGATATCCGCTTGTCGAGAGTCATCAAAATTATATTCGGCATGCACAGGAGGACACCGCCCCTCCAACTTTTTTGACAATAATGTCATTTGTCGGTGATGCTTCATATGACTTGGCACAAACTGTGCCGTGCTAGAGTTTGAAACCTGCCATGAAAGCTCCCGAGGAGGCCCAACGACAAATTTGGGTCACCACAAAAGCAAAGCCCTCGGCAAGCCGAATAACCACAGGTCCTCCGGATTCACTCGATGTATATTCTGCTATCTCAACA
    Z78524.1
    GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAACCGCCCAAGGGTCATTTTTCCTTCATTGGGGTTCAAACAAGTTCTTTTATAGGTTCAACAAATACGACAACATGGGGAATTCAAATGGCAGCCTTGCAGCTACCAACATTCTCCACCTCCTGCCGAGTTCCGAGCTATCAAAATACTAATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATATTTCATCCGCACAATGCTAACCTAGCTATATGCTGGGCAAGCATTGGNAGGAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGGATGTAAAGAGCCAAGTTCCGTTGCGAGAGTCATCAAAAAATTCATTAGCATGCACAACAGCATGCCGCCCCTCCGACTTTTTTAACAACAATGCCATTCATCAGGGATGCTTATATGACTTGGCACAAACTGCGCCGTGCTAAAGGTTTGAACCCGCCATGAAAGCTCCTAAGGAGGCCCAATTAAGGTCACCACAAAAGCAAAGCACTGACGAGCCAAGTAACCACATGTCCTCCATATTCACTCAATCATATATTCTACTATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78523.1
    CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGGTTCAAACGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCGCATCCACGCAATGCAAACCAATCTATATGATAGGAAAGCATTGATAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACATATGTCGAGAGTCAATAAAAAATTCATTATGCATGCATGCAAGAGCACGCTTCCACTCTAGCTTTTGGGACAATAATGTCATTCCTCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCACACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAGTTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTGCCGCAGGTTCACCTACGGAAACCTGGTTACG
    Z78522.1
    CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCCCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACGCAATGCAAACCTGTCTATATGACGTATAACCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTTCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCGATCACACCACGATATGTCGAGAGTCATAAAAAATTCATTTGCATGTATGCAATAGCACGCTTCCACTCCAACTTTTTGGACAATAATGTCATTCATCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGGGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78521.1
    GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAGTTGTTTGCCTCCATTGGGATTCAAACATGTTTTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGTACCTAGGCACAAGGGGCCATTCTTTCACATCCACGCAATGCAAACCTATCTATATGACATGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTTGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAATTCATTTGCATGCATGCAATAGCACGCTTCCACTCCAACTTTTTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGGCAAATTAGGGTCACCGCCAAAGCCAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78520.1
    AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCCAGGTCCGACAAACACGACAATATGGGGGATTCACATGACAACCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCGATATGACAGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCACATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTTCGCTGCGTTCTTCATCCGATGCAAGAGCCAAGATATCCCTTGTCGAGAGTCATAAAATATTCTATTGCATGCATGCAATAGCACGCTTCTACTCCAACTTTTTGGACAATAATGTCATTCATCGGTGATGATTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78519.1
    TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCAACAAACGCGACAATATGGGGGATTCAAATGATAGCCTTATATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACATTAAGTTCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCTATATGATGGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAAATTCATTTGCATGCATGCAGTAGCACGCTTCCACTCCAACTTTTGGACAATAATGTCATTCATCGGCAATGCTTCATATGACTTGGTGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCAGATTCACTCGATCATAT
    Z78518.1
    GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCAAGGGTTGTTTTGCCTCCATAGGGGTTCAAATAGGTTCTTCTGTGGGTCCGGCAAATACGACAACATGTGGGATTTGAACGACAGCCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTTGGGCCATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGAGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGTATGCGTAACACACACCGTCCATCCGAGTTTGGGGGACAATAATGTAATTCGTCGGCGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGCTCAAACCCGCCATAAAAGCTCTCGAGGATGCCCAACTACAAATCAGGGTCACCACAAAAGTTAAGCCTTCGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCGAATATTCTACTATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78517.1
    GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAATGACAGTCTTGCGACTATCAACATGATCCACCTCCTGCCGGGGTTCGGGCCACCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGTCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGGAAGAGCAAGATATCCGTTGGCCGAGAGTCATCAAAAATTTATTGGCATGCACAACAGCACACCGCCCATCCGACTTTTGTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCAAACTGCACCGTGATAGAGGTTCAAACCCGCCATAAAATCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCTGGATTCACTCGATCGTATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78516.1
    TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCAATCACTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAACGATAGTCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTCGGGACATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTAAGATATCCGTTGACGAGAGTCATCAAAAAATTATTGGTATGCACAATAGCACACCGCCCATCCGACTTTGGTGACAATAATGTCATTCGTCGGTGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGTTCAAACCCGCCATAAAAGCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCATATATTATACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78515.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAATCAACCGCCCATGGGTTGGGTTGATGTGGCACCAATGGGGTTCCGAAGGGTCCTACGAATTATGATGGCCCATCAAATCGGACAACATACGGCATTCGCACAACAGCTTTGGCTCAGCCATTGCCCATCCACCTCTTGTCGAATTTTGGACCATCAAAAGCCCGAGGTCTTAGACCCATTGAACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGACCGAACGGTATTCTTGGACAACCATTGGTGGGAGAGACACAGTACACAACGCCCAGGCAGGCGTGCCCTTAGCCTGACGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTAAAAATGTAACTTGGGGGGCACGAATCAAGTGCCACCCCCACCGGCTACAAGGGAACACCCCTATGCCGATTCTTGTCTCAAAAGACTTGGCGCGAAGCTGCGCCGGGCTAGAAGTTTGATCGCCTTCAAGAACACTCCCAAGGAGGCCCGATGGCAGATCAAAGGCCACCACAGTAGCGAAAGCCCCTGGGAGATCGAAGTACACCGGTCCTCCAGATTAACTCGATCGTTTATTCTACGGTCTCAGCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78514.1
    TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCCATGGGTTGGGTCGACGTGCCTCCGACGGGGTTCGAAAGGGATCATTCCGGCCCATCAAATACGACAACATAGGGCATTCGCACAACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGACGAGAGTTGTCAGACGAATCTTAAGGGGGCACAGCAGCACGCCGGCCCTCCCCGGTCCATCCATGGAGGACGAAGACCCTTGTCGGTTCTATCTCATTTGACTTGGAGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATAGAGAACTCTCCCAAGGAGGTTCGATGGGAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTATCTCGATCGTATATTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78513.1
    CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGGTTGGGTCGACGTGCCTCCGAAGGGGTTCGAAAGGGATCGTTCTGGCACATCAAATACGACAACATAGGGCATTCGCACAAAAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACAGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCGAAAATTCTTTCGGGCACAGCAGCATGCCGGCCTCCCCGGCTCCATGGAGGACGATGCCCTTGCCGGTTCTATCTCATTTGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATGGAGAAATCTCCCAATGAGGCTCGATGGCAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTAACTCGATCGTGTATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78512.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAGCTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCGAGAGTTGTCAAAAAATAATTTCATTGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGGGGACGATGCCCTCCGCCGGTTCCATCTCATTTGATTTGGGGCGAAACTGCGCCGGGCCAAGGGTTTCATTTGCCATCAAGAATTCCTCCCAGAGGGTTCGACGGAAATTCACGGTCACCACAGTAGCCAAGGCCCCTGGGAGACCAAACTACACCGGCCCTCCGGATTAATTCGATCGTTTTTTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78511.1
    TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGCTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTTCATCGGATGCAAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAAAAATTCATGGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAATTCTCCCAAGGAGGTTCGACGGAAATTCACGGCCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78510.1
    TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAAAACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAATCTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCCATCCGATGCAAAGAGGCGAGATATCCGTTGCGAGAGTTGTCAAAAAATCCGAGGGGGGAACGGAAGAATGCCGGCCCTCCCCGGTTCCATGGGGGGGGACGCCCCCTGCCGGTCCTATCTCATTTGATTGGGGGGGAAATTGCGCCGGGCCAAGGGTCCATTGGCCACCAAGAATTCTCCCAGGGGGGTTCGATGGAAATTCACGGCCACCAAGGGGGGAAAGCCCCTGGGGAGACCAAATTAAACCGGTCCTCCGGTTTAATTCGATCGTTTTTTCGGGGGCTTAAAAAGGAATCCTCCCGAAGGTCACCTCGGAACCCTGGTTAG
    Z78509.1
    TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACCCCCCATTGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAAAACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAACCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCTGCTTATCACATTTAAACCTTGTTACGTCCGTTGCCGAGAGTTGTCAGAAAAATTCATGGGGGGGCACGGAGCATGCCGGCCCTCCCCGCTCCATGGAGGACGACGCCCTCCGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78508.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAATTCATTGGGGGCACGGCAGCATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78507.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGTGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAAATTCATTGGGGGACGGAAGCACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGAAATTCACGGTCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78506.1
    TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTATGAGCTGAGGGCTCCTTGAACGAGAGTTGTCAGCCCGATTCAAAGGGGGTACGGCGGCATGCCGGTCCTCCCCGGTTCCATGGAGGACGAAGCCCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78505.1
    AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGACCGAACAGTATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTGACGAGAGTTGTCACAAAGATTCATTGGGGGTATGGTGGTATGCCGGGCCTCCCCGGTTCCATGGAGGACGAAGACCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCACGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAATATCATGGTCACCGCAGTAGCCACCGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78504.1
    TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGTTTTGTCGCAGCCACCGTCCGTCCACCTCTTGTCGGATTTGGTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGTCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAGATTCATTGGGGGCACGGCGGGCATGCCGGCCCTCCCCGGCTCCATGGAGGGACGATGCCCTCTGACGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78503.1
    TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAAAGAGCTGAGACTCCGTTGCCGAGAGTTGTCAAAATAATTCATTGGGGGTACGGTAGCATGCCGGCCCTCCCCGGTTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGTAAATCACGGACACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAAGGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78502.1
    GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTACGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTTCGCTGCGTCCTTCCATCGATGCAAGAGCCGAGATTCCGGCCGAGAGTTGTCAATATAAAATTCATTGGGGGGACGGTAGTATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGACCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGTCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGGGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78501.1
    TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCTCGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACAACACTTATCACAATTTCGCTGCGTCCTTCCATCCGATGCAAGGAGCCGAGATTTCCCGTTTGGCCGAGAGTTGCCCAAAAAAAAATTCATTGGGGGCACGGCAACACGCCGGCCCTCCCCGGCTCCATGGAGGACGATTCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTTCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78500.1
    CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCGAGATTCCCGTTGGCCGAGAAGTTTTCCCAAAGAAAAATTCATTGGGGGCACGGCAGGACGCCGGCCCTCCCCGGCTCCATGGAGGGCGATGCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCAATTATTTTGCGGTCTCAACAATGAGCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78499.1
    GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACACTTATCGCATTTCGCTGCGTCCTCCATCCGATGCAAAGAGCCGAGATATCCGTTGGCTGAGAGTTGTCAAAAAAATAATTGGGAGAGGTCAAGGCCACATGCCGTCCCCTCCTCGAATTCAACGAAGGGGGTATGTCCTTCACCATTTATGTGTCATATGACTTGGTGCAAAACTGCGACGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGAAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGCGATCTCAACAATGATCCCTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78498.1
    GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAACTGATCACGGGTTGATGTGCCTCCAGTGGGGTTCACAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAGTTTTGTTGTAATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGTACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCCAAAAATAATTTGGAGAGGGGGAGTGATTGTAGGGTCAAGGGCAAAATGCCGCCCCCTCCTTCGCCATTATTGTGTCATATGACTTGGGGCAAAACTGCCCCGGGCAAGGGTAAGATCGCCGGCAAGAAAACTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCACCCATGGGTGACCGGAGTAAACTGATCCTCTGGAGTAACTCGATCAATTATTATGTGATCTCAACAATGACCTTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78497.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTGTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCGCGTTGCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACACCCGCACCGGCTGAACAGTATGCCTGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAGAGCTGAGATCTCCGTTGCTGAGAGTTGTTTCACAAAAAATAATTTGGAGAGGGGGAGGGAGTGTAGGCTCAAGGCCACATGCCCTTCACCTCCTTGAATTCACAAAGGGCCGTGCCCTTCACCATTTATGTGTCATGTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCGCCCATGGGCGACCAAAGTAAACTGATCCTCTGGATTCACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78496.1
    GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGGTGCGTCCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTCAAAAAATTAATTGGAGAGGTCAAGGCCACATGCCGGCCCCTCCTCGAATTCACGAAGGGATATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGCAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78495.1
    CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTGGATTATGCCGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAACAGTGTTGTTGTAGCCTGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCTGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTGAGGGAGAGAGATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTCATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGGATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTCAGAAAAATACTTTGGAGAGGGGGAGAGCGAGTGCAGTGTCAAGGCCACATGCCGCCCCCCCTCCTTGAATACACGAAGGGCTATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACCGCGCCGGACAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCTCCATGGCAACTCTAGGTCACCGCAATAGCAAGCGCCCATGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78494.1
    CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTACGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGACAGGCGTTCCCTTGGCCTGTTGGCCTCGGGCCCAACTTGCGTTCAAAGACTCGATGTTAACGGGATTCTGCAATTCACACCTGAGATATCCGAAGTTGAGAGTGTTACCTCAATAATGGGGAGTTGGTCAAGGCAGTATGCCGTCCCCCTTCATTAATTATGGGTCATTTGATTTGGCGCATTACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGTTCCCAAGGAGGCCCGATGGTAAATCTAGGTCACTGCAACAGTAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGACCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78493.1
    GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCACACCGGATGACCAGTACGTTTGGACAGCCTCATTAAGGGAGAGATATGAGTCGCAATGTCCAGGCAGGCGTCCCCTTGGGCTGATGGCCTCGCGCGCAACTTGCGTTCAAAGATTCGATGGTTCACGGGATTCTGCAAGGCACACCATTTATCGAATCTCGCTGCGTTCCTTCATCGATGCATGAGCCGAGATATCCGTTGCTGAGAGTTGTTACCTAAATACTTGGGGGCTGGTCAAGGAAGAATGCCGCCCCCCTTCACCACTTATGTAAAATTTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGCTCCCAGGGAGGCCCGATGGAAATTCTAGGTCACTGCAACAGAAAATGCCCATGGGTGACCAAAGTAAACCGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78492.1
    TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCTGATTAAGCTGGCCCATGTAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAATTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCGTGCAGCTCTTAGACCCCCCGTCCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTCCACGGGATTCTGCAAATCACACCAGTATCGCATTTCGCTGGGTTCTACATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGAGTGGGTCAAGGCAGTACGCTCCCCCCCTTCACCAATTATGTGACATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78491.1
    GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTCCAAAGGGTCTTCTGATTATGCTGGTCCATGTAACACAACAACCCGGGGCATTCGCACAACATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAATATAATTTGGAGCGGGTCAAGGCAGCACGCCTCCCCCCAACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCTAAGGAGGCCCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78490.1
    TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGACGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAAATTATGCTGGCCCATCTAATACGACAACGCGGGGCATTCGCACAACACTTGTTGCAGAACAGTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCGTTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCGAGGCAGCATGCCGCCCCCTTCACCAATTATGTGCCATATGACTTGGCGCAAAACTGCGCCGGACTAGGGTTCAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78489.1
    GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCATGGGGTTCAAAAGGGTCTTCAGATATGCTGGGCCCATCTAATACGACAACTCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAATCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACCCACATCCCCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATCTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTCGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCTCGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78488.1
    AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGGGTTGATGGGCCTCCGATGGGGTTCAAAAGGGTCTTCAGATTATGATGGCCCATCTAATACGACAACCCGGGGGATCCGCACAACAGTTTGGTTGTGGCGTTGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGGCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTCCATCCCGTTGATGAGAGTTGTTAGAAAATAATTAGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCACCAATTATGTGTCGTATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTAAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAACAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTGCGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACAG
    Z78487.1
    TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTACGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGCATTTAGGACCATCGATAGCCCATGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCGAACTCACATCCGCACCGGCTGAACAATATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTAAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGAGTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78486.1
    TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGGCACTGCAATAGCAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATAATCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGATTCACCTACGGAAACCTCGTGACG
    Z78485.1
    TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAGGGGTCTTTAGATTGTGCTGGCCCGTCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGGCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGCCCAAACTCACATCCGCACCGGCTGCACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCCCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGATATCTCAACAATGATCCATCCGCAGATTCACCTTCGGACACCAGGTTCAG
    Z78484.1
    AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTGTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGGGAGGGTAGGAAACATGCCGCCCCTTCACAATTATGTGTCATATGACTTGGGCCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGGAGGCTCGATGGTAAACTCTCGGTCACTGCAATAGCCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCCCAGGTTCACCTACGGAAACCTTGTTACG
    Z78483.1
    TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGATATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATTGTTCACGGGATTCTGCAATTCACACCATTTATGATTTCGGTTTGGTCTTCATCCATGCAAGAGGCCAGATTTCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78482.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACGCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGNTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGAAGAGGCGAGATATCCGTTGCNGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTACGGTTTAGATCGCCAGCAAGAAAGCTCCCAGGAGGCTCGATGGCAAATCTCGGTCACTGCAGTAGA
    Z78481.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTCCGTTGCTGAGAGTTGTTAAAAAAATGATTTGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78480.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCGCACAACATTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78479.1
    ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTGATCTTGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAGACTCGATGGTTCACGGGATCTCTGTATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGGCAAGGCAAGGCAGCATCCCTCCCCTTCCAGTAATGTCTCATATAAATTGGCGCAAATCTGCGCCGGGCAAGGGTTTAGTTCGGCTCGATGGCAAATCTAGGTCACTGAAACAGCAAATGCCCATGGGTGACCAAAGTAACCTGATCCTCCTGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78478.1
    GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAACTGATAATGGGTTGATGGGCCTCCAATGGGGGTCAAAAGGGTCTTTAGATTATGCTGACCCATCTAAACGACAACCGGGGCATTCGCACAACAGTTCTGTTGTCACATAGCTCGTCCACCTCTTGCCGTATTTAGGACCATCCAACCCATACAGCTCTTAGACCCCCGTACCAAAGGACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGCAGCCTCATTAAGGGAGAGATATGTCTCACAATTTCCAGGCAGGCGTCCCCTTGACCTGATGGCCTCGGACGCAACATGCGTTCAAAACTCGATGTTCAGGGATCTGCTATCACGCCATTATCAATATTTGGGGGAGCCATGGCAGCATGCCGCCCCTTCCAGTTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAACTCTAGGCCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCCCCAGATTAACTCGATCAATTATTATGTGATCTCAACACTGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78477.1
    GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCCGCAACTTGCGTTCCAAAGACTCGATGGTTCAGGGATTCTTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAGAGCCGAGATACTCGTTGCTGAGAGTTGGTTAACAAATATTTGGGGAGGGCAAGGCAGCATGCCGCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGCAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78476.1
    GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGATTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCGCAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTGGATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCGGATGGCCTAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGCCAAGGCAGCATGCCGCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78475.1
    ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCATTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGACGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78474.1
    AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAGAAGGGTCTGTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGTGTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTCAAAAATATAATTTGGGGAGGGTCAAGGCAGGATGCCGCCCTTCCAATTATGTGTCATATGATTGGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGAAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTACGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78473.1
    CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGGCAAGAACAAGGGGCCAAACTCACATCCGCACCGACTGAACAGTATGTATGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAACAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCACCAAGAAAGCTCCCATGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78472.1
    GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATTGGCCTCTAAAGGGGTTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGAACTACAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAAGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATAGACAGCCTCATTAAAGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATTCCGTTGCTGAGAGTTATCAAAAAATAATTTGGGGAGGGTCTTAGACTGCATGCCGCCCCCTTCCAATTATGTGTGAAATGACTTGGCGCAAGACTGCGCCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78471.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTACTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGATGTTGAGAGTTGTCAAAAAATTAATTTGGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCATTATGTGTCGGGGTGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGACAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACCGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78470.1
    AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTAGTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTAAAAAAATAATTTGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78469.1
    AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGCACTATCATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCTAGACCGCATGCCGCCCCCTTCCAATTATGTGTGATATGACTTGGCGCAAGACTGCGCCGGGCAACGAATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78468.1
    AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCACACAATAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTCCATCAAAAGCCCAGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGAGTTCCGGAGCTGAGAGTCGTTAAAAAAATGATATGGGGAGGGTCAAGGCAGCATGCTGCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTCAGATCGCCATCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATTCCCATGAGTGACCAAAGTAAACTGATCCTCCAGATTAACTCAATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78467.1
    TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAGAGGGTCTTTAGATTATGCTGGCCCAGCTAATACGACAACCCGGGTCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCACCAAAGGCACATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTACCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78466.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACGCAAGAACAAGGGGCCAAACTCACATCCGCACAGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCAAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTACAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGAAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78465.1
    GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACTCGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTTCGTTCAAAGACTCGATGGTTTACGGGATTCTTCAATTCACACCAATTATCGCATTTCGCTGCGTTCCTCATCGATTCAAGAACCGAGAAATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78464.1
    GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATTACGAAAGCCCGGGGAATCCGAACAGAAATTTGGTTGAAGAATAGTCCGCCCACCTCTTGCCGTTTTTAGGACCACCAAAAGCCCATGAAGTTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGAACCGGTTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATTTGATTCGAAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGAACTTGGCGTTCAAAGACTCGATGGTTCACGGGATTCTGAAATTCACACCATTTATCGGATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCTCCGGGCAAGGGTTTAGATCTCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGATCAATTATTATGTGATCTCAACAATGACCCTTCCGCTCACCTACGGAAACCTTGTTACG
    Z78463.1
    GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGTATTCGCACAGCAATTTTGTTGCAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAACATGCCGCCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTCCCCCGGGCGAGGGTTTAGATCCCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGAACAATGATTATGTGATCTCAACAATGAACCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78462.1
    ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTGCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGGTCACATCCGGAGACCTCGTGACG
    Z78461.1
    TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCATGGGTTGATGGGCTCCAATGGGGTTCAAAGGGTCTTCAGATTACGGTGGCCCATCTTATACGACAACCCGGGGGCTTTGGACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGGATGTATGGACAGGCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGGAATTCACACCATTTATCGGATTTCGGTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTGAAAAAATTATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGTCATTTGACTGGGGGAAAAATTGCGCCGGGTAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTCCAGCAGCAACTGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78460.1
    TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTGGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78459.1
    AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAATTTATCGCAATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78458.1
    CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCGTACGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78457.1
    CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACGATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGGTCCTTCATCCGATGAAAGAGCCGAGATATCCGTTGTTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGCCATTTGACTGGGGGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTACAACAGAAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78456.1
    GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGACCGGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTCCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78455.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGAAGCTCTTAGACCCCCCGTGCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCAACCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGGATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGAAAGGGTTTAGTTCTCGCCAAAAAGAAAGTTCCCAAGGGGGTTCGATGGAAATTCTAGGTCACTCCAAAAGAAATTGTTCATGGGTGACCAAAGTAAAAAGTTCCTCCAGATTAACTCGATCAATTTTTATGTGATCTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTGGTTACG
    Z78454.1
    GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGGTTGATGGGCCTCCAATGGGTTCAAAAGGGGCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAGAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGTCTGATGGCCTCGGCCGCAACTTGCGTTCACAGACTCGATGGTTCACGGGATTCTGCAAATCACACCATTTATCGCATGAGTTCCGTTGATGAGAGTTGTTAACAATAGTGGGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCGATTATGTGTCAAATGACTTGGAGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGTAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78453.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAGTTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78452.1
    TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACAGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCATCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATTCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAAAAAGAAAGCTCCCAAGGAGGCTTGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAATAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78451.1
    GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCTATAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGGCTTCGCACAACAATTTTATTGTAGGGATAGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAATCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGGTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATTTGACTCGCGATGCCCAGGGAGGGGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGATTTCGCTGGTTCTTCATCGATGAAAGAGCGAGATTTCCGTTGTTGAGAGTTGCAAAAAATACTTTGGGGAGGGTCAGGGCAGGATGCCCGCCCCCTACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCACGGGTTTAGATCTCGCCAGCAAGAAAGCTCCCAGGGGGGCTCCATGGAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTACACCTACGGAAACCTTGTTACG
    Z78450.1
    CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGTCCTCCAATGGGGTTCAACAGGGTCTACAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGCACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCAGTTTTGAGGACCATCATATGCCCATGCAGTTCTTAGACCCCCCGGACCAAAGAACAAGGGGCCACACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAGGGGAGAGATATGACTCGGAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGGCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGTAAATGCTCATGGATGACCAAAGTAAACAGATCCTCCAGCTTAACTCGATCAATTATTATGTGATATCAGCAATGATCCTTCC
    Z78449.1
    GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTAGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGGAGGGTCAGGGCAGGATGCCACCCTTCACCAATTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78448.1
    CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTGATCATGGGCTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTNTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGCGTTCTTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAAGGGAAGGATGCCACCCCCCTTCACCAATTATGTGTCATACGACTTGGCGCAAAACTGCGCCGGGAAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78447.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATAGGATGCCACCCCCTTCACCACTTATGTGTCATTTGACTGGGGGCAAAATTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAAAGCTCCCAGGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGGAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78446.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGAGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGGCCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCTGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTATCGCAATTTCGCTGGTCCTCCATCGATGAAGAGCCGAGTATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGGAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78445.1
    ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTGCACGGGATTCTGCAATTCACACCATTTATCGCATTTGCCGAAGACTGAGTCTCCGTTGCTGAGAGTTGTTAAAAAAATAGTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCAAACGACTTGGAGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78444.1
    AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCTTCAGATATCGGTGGCCCATCTTATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCATCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTAAGGGAAGAATGCCACCCCCTCCACCAATTTTGTGTAATTTGATTGGGGGAAAAATTGAGGGTTTAGATATCGCCAACAAGGATGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGTTCACCCTACGGAAACCTTGTTACG
    Z78443.1
    CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTATGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGGTGTGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCCCAAAACTCTGCCGGGCAAGGGTTTAGATATCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78442.1
    ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTTGTTGTGGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCGGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTTTGGATAGCCTCGTTAGGGGAGAGATATGATTCCCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCAAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78441.1
    CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCACAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAACTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATATCGGCAATGACCTTCC
    Z78440.1
    TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGCTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCAACGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGATATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAATGCTCCCAAGGAGGGTCGATGGAAATTCTAGGTCACTACAACAGAAATTGCCCATGGGTGACCAAAATATGCAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGGAGGTCCACCTACGGAAACCTTGTTACG
    Z78439.1
    GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTTCCGTATTTAGGACCATCCAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATTACCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCCAACTCCGCCGGGCAGGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATG



```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.gbk.fasta.txt", "fasta")
]
```


```python
# We can take the reverse complement, add rc_, add a description for each complement; this is done for each record.
```


```python
print(records)
```

    []



```python
len(records)
```




    0




```python
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.gbk.fasta.txt", "fasta")
    if len(rec) < 700
    
]
```


```python
len(records)
```




    0




```python
records = (
rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.gbk.fasta.txt", "fasta")
    if len(rec) < 700
)
SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    0

## MULTIPLE SEQUENCE ALIGNMENT 1-5

    ```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
aligner = Align.PairwiseAligner(match_score = 1.0)
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    3.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
aligner.mode = "local"
```


```python
target = "AGAACTC"
```


```python
query = "GAACT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    5.0




```python
aligments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
aligner.algorithm
```




    'Smith-Waterman'




```python
aligner.epsilon
```




    1e-06




```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
alignments = aligner.align(target, query)
```


```python
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7f63c2b77850>




```python
alignment.score
```




    3.0




```python
alignment.target
```




    'GAACT'




```python
alignment.query
```




    'GAT'




```python
print(alignment)
```

    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
alignment.coordinates
```




    array([[0, 1, 2, 3, 4, 5],
           [0, 1, 1, 2, 2, 3]])




```python
len(alignment)
```




    2




```python
alignment.shape
```




    (2, 5)




```python
alignment.mode = "local"
```


```python
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
local_alignment = local_alignments[0]
```


```python
print(local_alignment)
```

    target            0 TGAACT 6
                      0 -||-|- 6
    query             0 -GA-C- 3
    



```python
local_alignment.shape
```




    (2, 6)




```python
aligner.mode = "global"
```


```python
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
```


```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
print(local_alignment)
```

    target            0 TGAACT 6
                      0 -||-|- 6
    query             0 -GA-C- 3
    



```python
local_alignment.sort()
```


```python
print(local_alignment)
```

    target            0 -GA-C- 3
                      0 -||-|- 6
    query             0 TGAACT 6
    



```python
from Bio import Align
```


```python
from Bio.Seq import reverse_complement
```


```python
target = "AAACCC"
```


```python
query = "AACC"
```


```python
aligner = Align.PairwiseAligner(mismatch_score = -1, internal_gap_score = -1)
```


```python
aligner.score(target, query)
```




    4.0




```python
aligner.score(target, reverse_complement(query))
```




    0.0




```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
aligner.score(target, query, strand = "-")
```




    0.0




```python
# We can manipulate the scores based on different alignments with reverse complement and target
```


```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, query, strand = "-")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
print(alignments[1])
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
aligner.left_gap_score = -0.5
```


```python
aligner.right_gap_score = -0.2
```


```python
aligner.score(target, query)
```




    3.3




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
print(alignments)
```

    <Bio.Align.PairwiseAlignments object at 0x7f63b0729910>



```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
aligner.score(target, reverse_complement(query), strand = "+")
```




    -2.0




```python
# WE can check the number of alignments, print the best alignment, and calculate alignment scores for both the reverse "-"
# and forward "+" strands. The higher score (3.3) suggests a better match than (-2.0).
```
## BLAST

```python
from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = ""
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("m_cold.fasta2.txt", format = "fasta")
```


```python
print(record)
```

    ID: gi|8332116|gb|BE037100.1|BE037100
    Name: gi|8332116|gb|BE037100.1|BE037100
    Description: gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
    Number of features: 0
    Seq('CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...TTC')



```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
with open("m_cold.fasta.txt", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520B), microRNA
    length: 61
    e value: 4.91307e-23
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b (MIR520B), microRNA
    length: 60
    e value: 1.71483e-22
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133242|ref|NR_032573.1| Macaca mulatta microRNA mir-519a (MIR519A), microRNA
    length: 85
    e value: 2.54503e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| |||||||||||||||||||||||||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTGGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171267|ref|NR_035851.1| Pan troglodytes microRNA mir-519b (MIR519B), microRNA
    length: 80
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205302|ref|NR_030191.1| Homo sapiens microRNA 519b (MIR519B), microRNA
    length: 81
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171259|ref|NR_035850.1| Pan troglodytes microRNA mir-519a (MIR519A), microRNA
    length: 86
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205451|ref|NR_030222.1| Homo sapiens microRNA 519a-2 (MIR519A2), microRNA
    length: 87
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171276|ref|NR_035852.1| Pan troglodytes microRNA mir-519c (MIR519C), microRNA
    length: 86
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205290|ref|NR_030188.1| Homo sapiens microRNA 519c (MIR519C), microRNA
    length: 87
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171354|ref|NR_035860.1| Pan troglodytes microRNA mir-520f (MIR520F), microRNA
    length: 86
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205281|ref|NR_030186.1| Homo sapiens microRNA 520f (MIR520F), microRNA
    length: 87
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 4.60152e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||| ||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAGAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTCTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171394|ref|NR_035865.1| Pan troglodytes microRNA mir-522 (MIR522), microRNA
    length: 86
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205429|ref|NR_030218.1| Homo sapiens microRNA 519a-1 (MIR519A1), microRNA
    length: 85
    e value: 1.60609e-16
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| |||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205423|ref|NR_030217.1| Homo sapiens microRNA 522 (MIR522), microRNA
    length: 87
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171401|ref|NR_035866.1| Pan troglodytes microRNA mir-523 (MIR523), microRNA
    length: 79
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133247|ref|NR_032574.1| Macaca mulatta microRNA mir-519b (MIR519B), microRNA
    length: 81
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  |||||||||||||||||||||||||||||||||| |||| ||| |||||||||...
    CCCTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGTGCATCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205309|ref|NR_030193.1| Homo sapiens microRNA 523 (MIR523), microRNA
    length: 87
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 1.95662e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| ||||||||||||||||  |||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAATAAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| ||   | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTTATTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171437|ref|NR_035870.1| Pan troglodytes microRNA mir-526a-2 (MIR526A-2), microRNA
    length: 68
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||||||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||| |||||||||| | |||||| |||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGGAAGAAAAGAATGCGCTTCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 9.49283e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||| |||| || |   | | || ||||||||||||| | |||||||...
    CCCTCTAAAGGGAAGCGCATTCTTTTCTTCCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171211|ref|NR_035845.1| Pan troglodytes microRNA mir-518b (MIR518B), microRNA
    length: 82
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| ||||||||||| |||||||||||||||||| || |||||||||...
    CCCTCTAGAAGGAAGCACTTTCTGTTGTTTGAAAGAAAAGAAAGTGCATCATTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 3.31332e-06
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |  || || |||||| || | | ||| || ||||||||||||||| |||||||...
    CCTCTAAAATGATGCACTTTCTTTTCTTTCAAACAACAGAAAGTGCTTCCTTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205321|ref|NR_030196.1| Homo sapiens microRNA 518b (MIR518B), microRNA
    length: 83
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171342|ref|NR_035859.1| Pan troglodytes microRNA mir-520e (MIR520E), microRNA
    length: 86
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205265|ref|NR_030183.1| Homo sapiens microRNA 520e (MIR520E), microRNA
    length: 87
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133223|ref|NR_032569.1| Macaca mulatta microRNA mir-518c (MIR518C), microRNA
    length: 101
    e value: 1.01355e-12
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  ||||||||||||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133232|ref|NR_032571.1| Macaca mulatta microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||| |||||||  |||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGATTTCTGTGATCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171300|ref|NR_035855.1| Pan troglodytes microRNA mir-520a (MIR520A), microRNA
    length: 84
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171161|ref|NR_035839.1| Pan troglodytes microRNA mir-516b-1 (MIR516B-1), microRNA
    length: 89
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||| |||||||||    |||||| |||||| |||||||||||||||||||...
    CCCTCCAAAGGGAAGCACTTTCTGTTTGTTGTCTGAGAGAAAACAAAGTGCTTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | |||||| |||| | || ||| |   || ||| ||||||||||||| ||| |||||...
    CTCTAAAAGGAAGCACTTTGTTTTCTCTCAGACAACAAACAGAAAGTGCTTCCCTTTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | |||||||| |||||||||||||||||||||  |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAGCAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||| ||||||||   || | | | | || ||||||||||||| | | |||||...
    CCTCTAAAGGGGAGCGCTTTGCTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132725|ref|NR_032718.1| Macaca mulatta microRNA mir-516 (MIR516), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205396|ref|NR_030212.1| Homo sapiens microRNA 516b-1 (MIR516B1), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205294|ref|NR_030189.1| Homo sapiens microRNA 520a (MIR520A), microRNA
    length: 85
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    CTCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 0.00599063
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGAG...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| |||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171454|ref|NR_035872.1| Pan troglodytes microRNA mir-527 (MIR527), microRNA
    length: 86
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171333|ref|NR_035858.1| Pan troglodytes microRNA mir-520d (MIR520D), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    CTCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171291|ref|NR_035854.1| Pan troglodytes microRNA mir-519e (MIR519E), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||||||||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTTGTCTGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|269847477|ref|NR_031696.1| Homo sapiens microRNA 1283-2 (MIR1283-2), microRNA
    length: 87
    e value: 4.30974e-11
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||  || |||||||||||||||||||||||||||||||  |||||| ||| |||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAATCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205435|ref|NR_030219.1| Homo sapiens microRNA 527 (MIR527), microRNA
    length: 85
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|270133264|ref|NR_032577.1| Macaca mulatta microRNA mir-520a (MIR520A), microRNA
    length: 85
    e value: 1.50425e-10
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||   ||||||||||||||| ||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTATTTTCTGTTGTCTGAAGGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGA...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| ||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGA...
    ****ALIGNMENT****
    sequence: gi|262205358|ref|NR_030204.1| Homo sapiens microRNA 520d (MIR520D), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    TCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| ||||||||||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTTTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171169|ref|NR_035840.1| Pan troglodytes microRNA mir-516b-2 (MIR516B-2), microRNA
    length: 89
    e value: 5.25034e-10
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| |||||||||||||||||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTGAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTCAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171381|ref|NR_035863.1| Pan troglodytes microRNA mir-521-1 (MIR521-1), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133201|ref|NR_032565.1| Macaca mulatta microRNA mir-517 (MIR517), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||| |||| ||||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTGGTCTAAAAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205417|ref|NR_030216.1| Homo sapiens microRNA 521-1 (MIR521-1), microRNA
    length: 87
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133214|ref|NR_032567.1| Macaca mulatta microRNA mir-518a (MIR518A), microRNA
    length: 87
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTC-CTTT...
    ||||  | |||||||| ||||||||||||| || |||||||||||||||| ||||...
    CCCTACAAAGGGAAGCCCTTTCTGTTGTCTAAACGAAAAGAAAGTGCTTCTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205407|ref|NR_030214.1| Homo sapiens microRNA 517c (MIR517C), microRNA
    length: 95
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||||||||  |||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTTGTCT--AAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205377|ref|NR_030208.1| Homo sapiens microRNA 526a-2 (MIR526A2), microRNA
    length: 65
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| ||||||||||    |||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTG----AAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133311|ref|NR_032588.1| Macaca mulatta microRNA mir-523b (MIR523B), microRNA
    length: 89
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCT--GAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| || |||||||||||||||| ||   |||||| || ||| |||||| |||||||...
    CCCTCTAGAGCGAAGCGCTTTCTGTTGGCTAGAAAAGAATAGGAAGCGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133297|ref|NR_032585.1| Macaca mulatta microRNA mir-521 (MIR521), microRNA
    length: 87
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| ||| |||...
    CCCTCCAAAGGGAAGTACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 2.2325e-08
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132721|ref|NR_032717.1| Macaca mulatta microRNA mir-524 (MIR524), microRNA
    length: 85
    e value: 2.71974e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||  | ||| |||| |||||| |||||| ||||||||||| |||||||| |||...
    CCCTACAAAGGCAAGCACTTTCTCTTGTCTAAAAGAAAAGAAGGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205276|ref|NR_030185.1| Homo sapiens microRNA 519e (MIR519E), microRNA
    length: 84
    e value: 9.49283e-07
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||   |||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTT---TGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171363|ref|NR_035861.1| Pan troglodytes microRNA mir-520g (MIR520G), microRNA
    length: 89
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| | || |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-ACGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133196|ref|NR_032564.1| Macaca mulatta microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133191|ref|NR_032563.1| Macaca mulatta microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171387|ref|NR_035864.1| Pan troglodytes microRNA mir-521-2 (MIR521-2), microRNA
    length: 86
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205353|ref|NR_030203.1| Homo sapiens microRNA 521-2 (MIR521-2), microRNA
    length: 87
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133237|ref|NR_032572.1| Macaca mulatta microRNA mir-518f (MIR518F), microRNA
    length: 87
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||  ||||||||||||||| |  | | | | || | ||||||||||| | |||||||...
    CCTCTGAAGGGAAGCGCTTTCTTTCCTTTCACACAAGATAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171283|ref|NR_035853.1| Pan troglodytes microRNA mir-519d (MIR519D), microRNA
    length: 87
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTCAAACAAAGTGCCTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133288|ref|NR_032583.1| Macaca mulatta microRNA mir-520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||  |||| |||||||||    ||||||   |||| ||||||||||| || ||||...
    CCCTCTAGAGA-AAGCACTTTCTGTTTGTTGTCTGAGGAAAAACAAAGTGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205371|ref|NR_030207.1| Homo sapiens microRNA 516b-2 (MIR516B2), microRNA
    length: 85
    e value: 0.00171634
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||     |||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAG-----AAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205349|ref|NR_030202.1| Homo sapiens microRNA 519d (MIR519D), microRNA
    length: 88
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTTAAACAAAGTGCCTCCCTTTAGAG...



```python
# We can perform a BLAST using NCBI's online service, importing NCBIWWW to submita blastn query agaist
# the nucleotide data (nt) for a given sequence ID.
#The retrieved sequence is read, then used to do another BLAST search.
# NCBIXML is used to aprse the BLAST results from a file.
```

## CHALLENGE

```python
from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = ""
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("BRAF(2).fasta", format = "fasta")
```


```python
print(record)
```

    ID: BRAF
    Name: BRAF
    Description: BRAF Homo sapiens
    Number of features: 0
    Seq('CTTCCCCCAATCCCCTCAGGCTCGGCTGCGCCCGGGGCCGCGGGCCGGTACCTG...ml>')



```python
with open("BRAF.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520B), microRNA
    length: 61
    e value: 4.91307e-23
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b (MIR520B), microRNA
    length: 60
    e value: 1.71483e-22
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133242|ref|NR_032573.1| Macaca mulatta microRNA mir-519a (MIR519A), microRNA
    length: 85
    e value: 2.54503e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| |||||||||||||||||||||||||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTGGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171267|ref|NR_035851.1| Pan troglodytes microRNA mir-519b (MIR519B), microRNA
    length: 80
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205302|ref|NR_030191.1| Homo sapiens microRNA 519b (MIR519B), microRNA
    length: 81
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171259|ref|NR_035850.1| Pan troglodytes microRNA mir-519a (MIR519A), microRNA
    length: 86
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205451|ref|NR_030222.1| Homo sapiens microRNA 519a-2 (MIR519A2), microRNA
    length: 87
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171276|ref|NR_035852.1| Pan troglodytes microRNA mir-519c (MIR519C), microRNA
    length: 86
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205290|ref|NR_030188.1| Homo sapiens microRNA 519c (MIR519C), microRNA
    length: 87
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171354|ref|NR_035860.1| Pan troglodytes microRNA mir-520f (MIR520F), microRNA
    length: 86
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205281|ref|NR_030186.1| Homo sapiens microRNA 520f (MIR520F), microRNA
    length: 87
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 4.60152e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||| ||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAGAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTCTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171394|ref|NR_035865.1| Pan troglodytes microRNA mir-522 (MIR522), microRNA
    length: 86
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205429|ref|NR_030218.1| Homo sapiens microRNA 519a-1 (MIR519A1), microRNA
    length: 85
    e value: 1.60609e-16
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| |||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205423|ref|NR_030217.1| Homo sapiens microRNA 522 (MIR522), microRNA
    length: 87
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171401|ref|NR_035866.1| Pan troglodytes microRNA mir-523 (MIR523), microRNA
    length: 79
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133247|ref|NR_032574.1| Macaca mulatta microRNA mir-519b (MIR519B), microRNA
    length: 81
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  |||||||||||||||||||||||||||||||||| |||| ||| |||||||||...
    CCCTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGTGCATCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205309|ref|NR_030193.1| Homo sapiens microRNA 523 (MIR523), microRNA
    length: 87
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 1.95662e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| ||||||||||||||||  |||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAATAAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| ||   | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTTATTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171437|ref|NR_035870.1| Pan troglodytes microRNA mir-526a-2 (MIR526A-2), microRNA
    length: 68
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||||||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||| |||||||||| | |||||| |||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGGAAGAAAAGAATGCGCTTCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 9.49283e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||| |||| || |   | | || ||||||||||||| | |||||||...
    CCCTCTAAAGGGAAGCGCATTCTTTTCTTCCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171211|ref|NR_035845.1| Pan troglodytes microRNA mir-518b (MIR518B), microRNA
    length: 82
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| ||||||||||| |||||||||||||||||| || |||||||||...
    CCCTCTAGAAGGAAGCACTTTCTGTTGTTTGAAAGAAAAGAAAGTGCATCATTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 3.31332e-06
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |  || || |||||| || | | ||| || ||||||||||||||| |||||||...
    CCTCTAAAATGATGCACTTTCTTTTCTTTCAAACAACAGAAAGTGCTTCCTTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205321|ref|NR_030196.1| Homo sapiens microRNA 518b (MIR518B), microRNA
    length: 83
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171342|ref|NR_035859.1| Pan troglodytes microRNA mir-520e (MIR520E), microRNA
    length: 86
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205265|ref|NR_030183.1| Homo sapiens microRNA 520e (MIR520E), microRNA
    length: 87
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133223|ref|NR_032569.1| Macaca mulatta microRNA mir-518c (MIR518C), microRNA
    length: 101
    e value: 1.01355e-12
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  ||||||||||||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133232|ref|NR_032571.1| Macaca mulatta microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||| |||||||  |||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGATTTCTGTGATCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171300|ref|NR_035855.1| Pan troglodytes microRNA mir-520a (MIR520A), microRNA
    length: 84
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171161|ref|NR_035839.1| Pan troglodytes microRNA mir-516b-1 (MIR516B-1), microRNA
    length: 89
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||| |||||||||    |||||| |||||| |||||||||||||||||||...
    CCCTCCAAAGGGAAGCACTTTCTGTTTGTTGTCTGAGAGAAAACAAAGTGCTTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | |||||| |||| | || ||| |   || ||| ||||||||||||| ||| |||||...
    CTCTAAAAGGAAGCACTTTGTTTTCTCTCAGACAACAAACAGAAAGTGCTTCCCTTTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | |||||||| |||||||||||||||||||||  |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAGCAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||| ||||||||   || | | | | || ||||||||||||| | | |||||...
    CCTCTAAAGGGGAGCGCTTTGCTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132725|ref|NR_032718.1| Macaca mulatta microRNA mir-516 (MIR516), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205396|ref|NR_030212.1| Homo sapiens microRNA 516b-1 (MIR516B1), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205294|ref|NR_030189.1| Homo sapiens microRNA 520a (MIR520A), microRNA
    length: 85
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    CTCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 0.00599063
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGAG...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| |||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171454|ref|NR_035872.1| Pan troglodytes microRNA mir-527 (MIR527), microRNA
    length: 86
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171333|ref|NR_035858.1| Pan troglodytes microRNA mir-520d (MIR520D), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    CTCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171291|ref|NR_035854.1| Pan troglodytes microRNA mir-519e (MIR519E), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||||||||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTTGTCTGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|269847477|ref|NR_031696.1| Homo sapiens microRNA 1283-2 (MIR1283-2), microRNA
    length: 87
    e value: 4.30974e-11
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||  || |||||||||||||||||||||||||||||||  |||||| ||| |||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAATCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205435|ref|NR_030219.1| Homo sapiens microRNA 527 (MIR527), microRNA
    length: 85
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|270133264|ref|NR_032577.1| Macaca mulatta microRNA mir-520a (MIR520A), microRNA
    length: 85
    e value: 1.50425e-10
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||   ||||||||||||||| ||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTATTTTCTGTTGTCTGAAGGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGA...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| ||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGA...
    ****ALIGNMENT****
    sequence: gi|262205358|ref|NR_030204.1| Homo sapiens microRNA 520d (MIR520D), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    TCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| ||||||||||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTTTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171169|ref|NR_035840.1| Pan troglodytes microRNA mir-516b-2 (MIR516B-2), microRNA
    length: 89
    e value: 5.25034e-10
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| |||||||||||||||||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTGAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTCAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171381|ref|NR_035863.1| Pan troglodytes microRNA mir-521-1 (MIR521-1), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133201|ref|NR_032565.1| Macaca mulatta microRNA mir-517 (MIR517), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||| |||| ||||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTGGTCTAAAAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205417|ref|NR_030216.1| Homo sapiens microRNA 521-1 (MIR521-1), microRNA
    length: 87
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133214|ref|NR_032567.1| Macaca mulatta microRNA mir-518a (MIR518A), microRNA
    length: 87
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTC-CTTT...
    ||||  | |||||||| ||||||||||||| || |||||||||||||||| ||||...
    CCCTACAAAGGGAAGCCCTTTCTGTTGTCTAAACGAAAAGAAAGTGCTTCTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205407|ref|NR_030214.1| Homo sapiens microRNA 517c (MIR517C), microRNA
    length: 95
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||||||||  |||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTTGTCT--AAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205377|ref|NR_030208.1| Homo sapiens microRNA 526a-2 (MIR526A2), microRNA
    length: 65
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| ||||||||||    |||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTG----AAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133311|ref|NR_032588.1| Macaca mulatta microRNA mir-523b (MIR523B), microRNA
    length: 89
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCT--GAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| || |||||||||||||||| ||   |||||| || ||| |||||| |||||||...
    CCCTCTAGAGCGAAGCGCTTTCTGTTGGCTAGAAAAGAATAGGAAGCGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133297|ref|NR_032585.1| Macaca mulatta microRNA mir-521 (MIR521), microRNA
    length: 87
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| ||| |||...
    CCCTCCAAAGGGAAGTACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 2.2325e-08
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132721|ref|NR_032717.1| Macaca mulatta microRNA mir-524 (MIR524), microRNA
    length: 85
    e value: 2.71974e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||  | ||| |||| |||||| |||||| ||||||||||| |||||||| |||...
    CCCTACAAAGGCAAGCACTTTCTCTTGTCTAAAAGAAAAGAAGGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205276|ref|NR_030185.1| Homo sapiens microRNA 519e (MIR519E), microRNA
    length: 84
    e value: 9.49283e-07
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||   |||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTT---TGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171363|ref|NR_035861.1| Pan troglodytes microRNA mir-520g (MIR520G), microRNA
    length: 89
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| | || |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-ACGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133196|ref|NR_032564.1| Macaca mulatta microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133191|ref|NR_032563.1| Macaca mulatta microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171387|ref|NR_035864.1| Pan troglodytes microRNA mir-521-2 (MIR521-2), microRNA
    length: 86
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205353|ref|NR_030203.1| Homo sapiens microRNA 521-2 (MIR521-2), microRNA
    length: 87
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133237|ref|NR_032572.1| Macaca mulatta microRNA mir-518f (MIR518F), microRNA
    length: 87
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||  ||||||||||||||| |  | | | | || | ||||||||||| | |||||||...
    CCTCTGAAGGGAAGCGCTTTCTTTCCTTTCACACAAGATAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171283|ref|NR_035853.1| Pan troglodytes microRNA mir-519d (MIR519D), microRNA
    length: 87
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTCAAACAAAGTGCCTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133288|ref|NR_032583.1| Macaca mulatta microRNA mir-520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||  |||| |||||||||    ||||||   |||| ||||||||||| || ||||...
    CCCTCTAGAGA-AAGCACTTTCTGTTTGTTGTCTGAGGAAAAACAAAGTGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205371|ref|NR_030207.1| Homo sapiens microRNA 516b-2 (MIR516B2), microRNA
    length: 85
    e value: 0.00171634
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||     |||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAG-----AAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205349|ref|NR_030202.1| Homo sapiens microRNA 519d (MIR519D), microRNA
    length: 88
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTTAAACAAAGTGCCTCCCTTTAGAG...



```python
# the e-value is 1.71483e-22 when compared to chimpanzee (pan troglodytes)
```


```python

```


## OPEN CV 1-3

```python
import numpy as np 
import matplotlib.pyplot as plt
%matplotlib inline
import cv2
```


```python
img = cv2.imread("Dawg.jpg")
```


```python
type(img)
```




    numpy.ndarray




```python
img_wrong = cv2.imread('wrong/path/doesnot/abcdegh.jpg')
```


```python
type(img_wrong)
```




    NoneType




```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f789fca0150>




![png](output_5_1.png)



```python
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7f789fc96890>




![png](output_7_1.png)



```python
img_gray = cv2 = cv2.imread("Dawg.jpg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (1440, 1440)




```python
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7f789fbe1d10>




![png](output_9_1.png)



```python
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f789fb73610>




![png](output_10_1.png)



```python
fix_img.shape
```




    (1440, 1440, 3)




```python
new_img = cv2.resize(fix_img, (1000,670))
plt.imshow(new_img)
```




```python
new_img.shape
```


```python
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img (0,0), fix_img, w_ratio, h_ratio)
```


```python
plt.imshow(new_img)
```



```python
new_img.shape
```



```python
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```


```python
type(fix_img)
```




    numpy.ndarray




```python
cv2.imwrite("Dawg.jpg", flip_img)
```


    


```python
True
```

```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
img = cv2.imread("Dawg.jpg")
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f887f448b50>




![png](output_3_1.png)



```python
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f887f38a450>




![png](output_5_1.png)



```python
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f887f2fe6d0>




![png](output_7_1.png)



```python
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f887f26e610>




![png](output_9_1.png)



```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f887f1e0210>




![png](output_10_1.png)



```python
img2 = cv2.imread("Walt.jpg")
img1 = cv2.imread("Dawg.jpg")
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f887f1ce290>




![png](output_12_1.png)



```python
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f887f13e410>




![png](output_14_1.png)



```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f887c295a90>




![png](output_15_1.png)



```python
img1 = cv2.resize(img1, (1200, 1200))
img2 = cv2.resize(img2, (1200, 1200))
```


```python
alpha = 0.5
beta = 0.5
```


```python
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma = 0)
```


```python
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7f887c2804d0>




![png](output_19_1.png)



```python
alpha = 0.6
beta = 0.4
blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7f887c1eca10>




![png](output_20_1.png)



```python
img1 = cv2.imread("Walt.jpg")
img2 = cv2.imread("Dawg.jpg")

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (100,100))
```


```python
large_img = img1
small_img = img2 

x_offset = 0
y_offset = 0

x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]

large_img[y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```


```python
True
```

```python
import cv2 
import matplotlib.pyplot as plt 
%matplotlib inline
```


```python
img = cv2.imread('rainbow.jpg')
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f2c1a3d4290>




![png](output_3_1.png)



```python
img = cv2.imread('rainbow.jpg', 0)
```


```python
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f2c1a2f86d0>




![png](output_5_1.png)



```python
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
```


```python
ret1
```




    127.0




```python
plt.imshow(thresh1, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f2c1825f890>




![png](output_8_1.png)



```python
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1,cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f2c18245a90>




![png](output_9_1.png)



```python
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f2c181b1650>




![png](output_10_1.png)



```python
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f2c18190e90>




![png](output_11_1.png)



```python
def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = "gray")
```


```python
show_pic(img_r)
```


![png](output_13_0.png)



```python
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_14_0.png)



```python
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_15_0.png)



```python
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
```


```python
show_pic(th1)
```


![png](output_17_0.png)



```python
blended = cv2.addWeighted(src1 = th1, alpha =0.6,
                          src2 = th2, beta = 0.4, gamma = 0)
show_pic(blended)
```


![png](output_18_0.png)



```python
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6, 
                         src2 = th3, beta = 0.4, gamma = 0)
```


```python
show_pic(blended)
```


![png](output_20_0.png)

```python
# This section uses OpenCV for fundamental image processing operations. It includes loading images in different color spaces (BGR, RGB, HSV, HLS, and grayscale), displaying them using #Matplotlib, and performing transformations such as resizing, flipping, and blending multiple images. Additionally, it covers image overlay techniques, where a smaller image is placed onto a #larger one using NumPy indexing. The section also explores thresholding methods, including binary, truncation, and adaptive thresholding, to enhance contrast and segment specific parts of an #image, making it useful for object detection and image preprocessing.
```

## CORNER DETECTION

```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
green_chess = cv2.imread("green_chess.jpg")
green_chess = cv2.cvtColor(green_chess, cv2.COLOR_BGR2RGB)
plt.imshow(green_chess)
```




    <matplotlib.image.AxesImage at 0x7fc61819d090>




![png](output_1_1.png)



```python
green_chess = cv2.cvtColor(green_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(green_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc6180c2a10>




![png](output_2_1.png)



```python
real_chess = cv2.imread("real_chess.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fc618030190>




![png](output_4_1.png)



```python
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc609d06310>




![png](output_5_1.png)



```python
gray = np.float32(green_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)
dst = cv2.dilate(dst, None)

green_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(green_chess)
```


   

```python
corners = cv2.goodFeaturesToTrack(green_chess, 64, 0.01, 10)
```


```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(green_chess, (x,y),3,(255,0,0), -1)

    plt.imshow(green_chess)
```


![png](output_8_0.png)



```python
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fc609af2390>




![png](output_9_1.png)

```python
# This section detects "corners" in two different chessboard images using OpenCV.

# Load and convert images** to grayscale for better corner detection. Apply the Harris Corner Detector** on the
# green chessboard image to find corners, but an error occurs because the image format isn't compatible with color marking
# Use the Good Features to Track method** to detect and mark corners in both chessboard images. Draw small circles on detected corners to highlight them.
# Display the images** with detected corners marked in red for the green chessboard and green for the real chessboard.
# This technique is useful in **computer vision tasks** where recognizing shapes, edges, and objects is needed.
```

## EDGE DETECTION


```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread("Dawg.jpg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7ff733953ed0>




![png](output_3_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7ff73393f9d0>




![png](output_4_1.png)



```python
med_value = np.median(img)
med_value
```




    80.0




```python
lower = int(max(0, 0.7*med_value))
upper = int(max(255, 1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7ff7338a7d90>




![png](output_6_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper + 100)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7ff733812690>




![png](output_7_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                  threshold1 = lower,
                  threshold2 = upper)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7ff7337f98d0>




![png](output_8_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7ff73376d350>




![png](output_9_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 50)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7ff7336dd1d0>




![png](output_10_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 100)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7ff7336bee90>




![png](output_11_1.png)

## FEATURE MATCHES

```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
def display(img, cmap = 'gray'):
    fig = plt.figure(figsize = (12,10))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
froot_loops = cv2.imread("froot_loops.jpg", 0)
display(froot_loops)
```


![png](output_2_0.png)



```python
cereals = cv2.imread('All_Cereals.jpg', 0)
display (cereals)
```


![png](output_3_0.png)



```python
orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(froot_loops, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches= bf.match(des1,des2)
```


```python
matches = sorted(matches, key = lambda x:x.distance)
```


```python
froot_loops_matches = cv2.drawMatches(froot_loops, kp1, cereals, kp2, matches[:25], None)
```


```python
display(froot_loops_matches)
```


![png](output_8_0.png)



```python
sift = cv2.SIFT_create()
```


```python
kp1, des1 = sift.detectAndCompute(froot_loops, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1,des2, k=2)
```


```python
good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
print('length of total mathces:', len(matches))
print('length of good matches:', len(good))
```

    length of total mathces: 12315
    length of good matches: 201



```python
sift_matches = cv2.drawMatchesKnn(froot_loops, kp1, cereals, kp2, good, None, flags = 2)
display(sift_matches)
```


![png](output_14_0.png)



```python
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(froot_loops, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm=flann_index_KDtree, trees = 5)
search_params =dict(checks = 50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)

good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
flann_matches = cv2.drawMatchesKnn(froot_loops, kp1, cereals, kp2, good, None, flags = 0)
display(flann_matches)
```


![png](output_18_0.png)



```python
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(froot_loops, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm = flann_index_KDtree, trees = 5)
search_params =dict(checks = 50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k =2)
```


```python
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
for i, (match1, match2) in enumerate(matches):
    if match1.distance <0.75*match2.distance:
        matchesMask[i] = [1,0]
        
draw_params = dict(matchColor = (0,255,0),
                  singlePointColor = (255,0,0),
                  matchesMask = matchesMask,
                  flags = 0)
```


```python
# flann based matching; fast librarying for nearst neighbors; faster than sift
flann_matches = cv2.drawMatchesKnn(froot_loops, kp1,  cereals, kp2, matches, None, **draw_params)
display(flann_matches)
```


![png](output_24_0.png)



```python

```

```python
# This section compares two images using OpenCV to find where the Froot Loops cereal appears in a larger cereal image. 

# Steps:
# Loading and display images** in grayscale for better feature detection. Find key features in both images using ORB and SIFT, which detect important points in an image. Match
# features between the two images using different matching methods (Brute Force Matcher and FLANN).
# Filter out weak matches** to keep only the strongest matches. Draw lines to show matching points between the two images.
```






## OBJECT DETECTION

```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
full = cv2.imread("sunflower_training.jpg")
```


```python
full =cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7eff8d9a5c50>




![png](output_3_1.png)



```python
test = cv2.imread('sunflower_testing.jpg')
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7eff8d8df710>




![png](output_4_1.png)



```python
print('Full image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Full image shape: (610, 612, 3)
    Training image shape: (177, 285, 3)



```python
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
```


```python
for m in methods:
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
        
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(top_left, top_left, bottom_right, (255,0,0),10)
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title('Detection of template')
    
    plt.suptitle(m)
    
    plt.show()
    print('\n')
    print('\n')
```


![png](output_7_0.png)


    
    
    
    



![png](output_7_2.png)


    
    
    
    



![png](output_7_4.png)


    
    
    
    



![png](output_7_6.png)


    
    
    
    



![png](output_7_8.png)


    
    
    
    



![png](output_7_10.png)


    
    
    
    



```python

```

```python
# This section is trying to find a smaller image (the "template") inside a bigger image by using a technique called **template matching** in OpenCV.  
# Here’s what’s happening step by step:
# - First, the images are loaded and converted so they display properly.
# - The size (width, height, and color channels) of both images is printed to make sure they match up correctly.
# - Different methods are tested to compare the smaller image with different parts of the bigger image.
# - The program looks at the results and decides where the best match is by checking for the highest or lowest value (depending on the method used).
# - A rectangle is drawn around the area where the best match was found.
# - Finally, two images are shown: 
#   1. A **heatmap** that shows how well different parts of the big image match the small one.
#   2. The **big image with the detected template** outlined in a box.
# - This process is repeated for multiple matching methods to see which one works best.
```
