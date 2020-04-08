# How Corona Virus Pandemic Got Me Into Bioinformatics

I knew nothing about bioinformatics till now apart from its a field in biology that deals with genetic codes in organisms. As all the activities are at stop now due to the corona virus pandemic. I had some time to waste. Well if you waste your time learning then its not wasted time, right? So, I went on to learn about corona virus. Internet is the oracle of information. Google is the key. Well this new virus originated from Wuhan, China. It transmitted from bats to human in the wild animal markets. This is not the first time some virus from other animal mutated to infect humans. Eating wild animals is very popular in China as well as other Southeast Asian countries. There is a global movement to close down these wild animal markets, but they work through mob. So, Closing them is not so easy even if China banned them officially. [You can watch the documentary here.](https://www.youtube.com/watch?v=Y7nZ4mw4mXw) At first, The Wuhan virus (first they didn't know of the type of virus it is, so they gave a location based naming) only transmitted from bats to humans, but soon it mutated to transmit from human to human. Wuhan went into lockdown. China sequenced the virus' complete genome and send it to rest of the world to work on the vaccine. They recognized it's form of 'Corona Virus' which is common in humans that causes seasonal flu. But its more leathal than seasonal flu virus. Where only 0.1% people die of seasonal flu, this virus kills at least 1-2% people. 


At this point I had some questions that needed solving. 
- **Whats in the genome; How did they sequence the genome of the virus.**
- **How did they know its similar to other Corona virus'.**
- **How to analyze it or what it tell us about the virus**

## Genome; How to sequence it

From intermediate biology, I knew every living oranisms have their genetic material in the cell called DNA or sometimes RNA. Segments of DNA or RNA is called gene that regulates certain traits in the organisms. And in reproduction it passes down to the offspring. Thats why children gets the some traits from their parents. Its called heredity. and all those readings about Mendel and his discoveries with those pea plants. I knew modern science we can sequence the genes. There was a billion dollar project to sequence human genome. We have around 3 billion of code in our genome. Some months earlier a bangladeshi scientist sequenced the genome of jute. And some information here and there. But nothing concrete. Now its time to dig deep. 

You can get the information from some google searches and wikipedia. Here I will only talk about the things thats needed. Genome of any living organisms (theres debate whether virus is living or not you can google about it) consists of four [nucleotides](https://en.wikipedia.org/wiki/Nucleotide) molecule. These are Adenine (A), Cytosine (C), Guanine(G), and Thymine(T). Sometimes Thymine is replaced with Uracil(U). These nucleotides joins together one after one in a long chain, as in human its 3 billions long. This chain is called the genome. [This video will help to understand the structure of DNA](https://www.youtube.com/watch?v=o_-6JXLYS-k)

So, How do they sequence this tiny series of molecules? These two videos describes the process flow where I got the gist of it. [How to sequence human genome](https://www.youtube.com/watch?v=MvuYATh7Y74&t=15s) and [DNA Sequencing](https://www.youtube.com/watch?v=ONGdehkB8jU). I wanted to know what happens in a laboratory setting. [This video helped me](https://www.youtube.com/watch?v=KfyAwAtyUQE). After all the process a sequenced genome look like this. 

'CTGTGGCCCTGATGGCTACCCTCTTGAGTGCATTAAAGACCTTCTAGCACGTGCTGGTAAAGCTTCATGC'

Lets see the genome of Wuhan virus now. You can go to this [link](https://www.ncbi.nlm.nih.gov/nuccore/MN908947). Copy the genome from the gene bank site and after some processing we will get the sequence that we can work on. First we copied the code in this [ file](SARS-Cov-2-genome-genebank-MN908947.3.txt). I will use python process and analyze the genome. Lets import and get on with it.


```python
# opening the genome file
genome = open("SARS-Cov-2-genome-genebank-MN908947.3.txt", "r").read()
# showing only the first 1000 letter
print(genome[:1000])
```

            1 attaaaggtt tataccttcc caggtaacaa accaaccaac tttcgatctc ttgtagatct
           61 gttctctaaa cgaactttaa aatctgtgtg gctgtcactc ggctgcatgc ttagtgcact
          121 cacgcagtat aattaataac taattactgt cgttgacagg acacgagtaa ctcgtctatc
          181 ttctgcaggc tgcttacggt ttcgtccgtg ttgcagccga tcatcagcac atctaggttt
          241 cgtccgggtg tgaccgaaag gtaagatgga gagccttgtc cctggtttca acgagaaaac
          301 acacgtccaa ctcagtttgc ctgttttaca ggttcgcgac gtgctcgtac gtggctttgg
          361 agactccgtg gaggaggtct tatcagaggc acgtcaacat cttaaagatg gcacttgtgg
          421 cttagtagaa gttgaaaaag gcgttttgcc tcaacttgaa cagccctatg tgttcatcaa
          481 acgttcggat gctcgaactg cacctcatgg tcatgttatg gttgagctgg tagcagaact
          541 cgaaggcatt cagtacggtc gtagtggtga gacacttggt gtccttgtcc ctcatgtggg
          601 cgaaatacca gtggcttacc gcaaggttct tcttcgtaag aacggtaata aaggagctgg
          661 tggccatagt tacggcgccg atctaaagtc atttgactta ggcgacgagc ttggcactga
          721 tccttatgaa gattttcaag aaaactggaa cactaaacat agcagtggtg ttacccgtga
          781 ac



```python
# we need to remove the numbers and white spaces between the letters 
for x in " 0123456789\n":
    genome = genome.replace(x, "")
print(genome[:1000]) 
# the length of the genome
print("\n length of the genome: "+ str(len(genome)))
```

    attaaaggtttataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttctctaaacgaactttaaaatctgtgtggctgtcactcggctgcatgcttagtgcactcacgcagtataattaataactaattactgtcgttgacaggacacgagtaactcgtctatcttctgcaggctgcttacggtttcgtccgtgttgcagccgatcatcagcacatctaggtttcgtccgggtgtgaccgaaaggtaagatggagagccttgtccctggtttcaacgagaaaacacacgtccaactcagtttgcctgttttacaggttcgcgacgtgctcgtacgtggctttggagactccgtggaggaggtcttatcagaggcacgtcaacatcttaaagatggcacttgtggcttagtagaagttgaaaaaggcgttttgcctcaacttgaacagccctatgtgttcatcaaacgttcggatgctcgaactgcacctcatggtcatgttatggttgagctggtagcagaactcgaaggcattcagtacggtcgtagtggtgagacacttggtgtccttgtccctcatgtgggcgaaataccagtggcttaccgcaaggttcttcttcgtaagaacggtaataaaggagctggtggccatagttacggcgccgatctaaagtcatttgacttaggcgacgagcttggcactgatccttatgaagattttcaagaaaactggaacactaaacatagcagtggtgttacccgtgaactcatgcgtgagcttaacggaggggcatacactcgctatgtcgataacaacttctgtggccctgatggctaccctcttgagtgcattaaagaccttctagcacgtgctggtaaagcttcatgcactttgtccgaacaactggactttattgacactaagaggggtgtatactgctgccgtgaacatgagcatgaaattgcttggtacacggaacgttct
    
     length of the genome: 29903


Well the genome has **29903** number of nucleotides in the sequence. You can google and match it with the other sources. 


```python
# checking if the string consist only the "atcg" no other letters and assign as the nucleotides
nucleotides = set(genome)
print(nucleotides)
```

    {'a', 'c', 'g', 't'}


At this stage I want to know some basic statistics of the sequence. Like how many "A", "T", "C" and "G" are there. Their percentage. 


```python
# Counting each nucleotide in the sequence
# creating a counter dictionary from the nucleotide set
counter = dict.fromkeys(nucleotides, 0)

# counting how many of each nucleotides are there
for i in genome:
    counter[i] = counter[i] +1 

print("Number of each Nucleotides in the sequence: " + str(counter))

# sum of all the nucleotides number is equal to the total number
length = 0
percentage = {}
for i in counter:
    length = length + counter[i]
    percentage[i] = str(round(counter[i]/len(genome)*100, 2))+"%"

    
print("Percentage of each Nucleotides in the sequence: " + str(percentage))
print("sum of all: " + str(length))
```

    Number of each Nucleotides in the sequence: {'a': 8954, 'c': 5492, 'g': 5863, 't': 9594}
    Percentage of each Nucleotides in the sequence: {'a': '29.94%', 'c': '18.37%', 'g': '19.61%', 't': '32.08%'}
    sum of all: 29903


Now I don't know what else to do with the genome. But there can be some improvements made to the process of getting the genomic data. Here I did it manually. While looking through the data sources, I saw only few has the Genome data. Usually researchers of various laboratory sequences the genome and then send it to some genome data storing organizations. They store it and make it available for other researcher to access it. This [article](https://www.yourgenome.org/facts/how-are-sequenced-genomes-stored-and-shared) describes the process and oranizations involve in genomic data storage. Each data store has their tools to sort and search those data. 



### Getting the genome data automatically

Now we will try to automate the process of getting the genomic data from those data sources. There is a good python library for bioinformatics called **Biopython**. We will use it to load data and do some basic analysis to understand the capabilities of the library. 


```python
# importing some utility module
from pprint import pprint as print

# importing the Entrez module from the Biopython library
from Bio import Entrez
# you need to setup email address to the Entreze module
Entrez.email = "abulkalamfaruk@gmail.com"
```

There are two helper function in python for helping to use any library. **dir()** function is used to know the modules and function inside any library and **help()** function is there to know the arguments and usage of the function. Also in jupyter notebook there is a hack for help(). you just have to put a "?" mark after the function to view the help. 

For the getting the genome data from the gene banks we need to use ***efetch()*** from the **Entrez** module. As argument to the function we need, which database to access which is nucleotide database **(db="nucleotide")**, **'id'** is the most important argument to put. If you go to any data sources online for gene information. for each genome you will get an ***accession number***. This accession number is the **id** in the *efecth* function. For our case its *MN908947*

There are two important **rettype** as in file type of genomic data. One is FASTA format and the other is GENE BANK (gb) format. You can learn about those formats from [here](https://www.genomatix.de/online_help/help/sequence_formats.html). For our case we would use FASTA format. 


```python
# getting the data from data bank
covid_fasta = Entrez.efetch(db="nucleotide", id="MN908947", rettype="fasta", retmode="text")
covid_genome = covid_fasta.read().split("\n")[1:]
covid_genome = "".join(covid_genome).lower()
print(covid_genome[:1000])
```

    'attaaaggtttataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttctctaaacgaactttaaaatctgtgtggctgtcactcggctgcatgcttagtgcactcacgcagtataattaataactaattactgtcgttgacaggacacgagtaactcgtctatcttctgcaggctgcttacggtttcgtccgtgttgcagccgatcatcagcacatctaggtttcgtccgggtgtgaccgaaaggtaagatggagagccttgtccctggtttcaacgagaaaacacacgtccaactcagtttgcctgttttacaggttcgcgacgtgctcgtacgtggctttggagactccgtggaggaggtcttatcagaggcacgtcaacatcttaaagatggcacttgtggcttagtagaagttgaaaaaggcgttttgcctcaacttgaacagccctatgtgttcatcaaacgttcggatgctcgaactgcacctcatggtcatgttatggttgagctggtagcagaactcgaaggcattcagtacggtcgtagtggtgagacacttggtgtccttgtccctcatgtgggcgaaataccagtggcttaccgcaaggttcttcttcgtaagaacggtaataaaggagctggtggccatagttacggcgccgatctaaagtcatttgacttaggcgacgagcttggcactgatccttatgaagattttcaagaaaactggaacactaaacatagcagtggtgttacccgtgaactcatgcgtgagcttaacggaggggcatacactcgctatgtcgataacaacttctgtggccctgatggctaccctcttgagtgcattaaagaccttctagcacgtgctggtaaagcttcatgcactttgtccgaacaactggactttattgacactaagaggggtgtatactgctgccgtgaacatgagcatgaaattgcttggtacacggaacgttct'



```python
# Checking the length to match the previous result
print("Length of the genome: " +str(len(covid_genome)))
```

    'Length of the genome: 29903'


## How did they know its similar to  Coronavirus family

So, the genome is sequenced now. or you can get any genome data from the NCBI database. But one of the important task when finding new virus is that what type of virus it is. The classification of virus doesn't follow the other life classification. It is done by simmilarities in its genome. For this there is an algorithm/tool. [BLAST](https://en.wikipedia.org/wiki/BLAST_(biotechnology)) or Basic Local Alignment Search Tool is used to match genomic sequence to the previous databases to see if the similarities found. There is a [online version of the tool](https://blast.ncbi.nlm.nih.gov/Blast.cgi) in the NCBI site where you can blast various genome. You can download it in the desktop too to use it in the command line. you can download the executables from [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/). The instructions can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/)

I am using Biopython package to access the NCBI Blast tool. I have matched our virus to SARS virus from the genebank accession number. It is sequenced in 2017 in USA. I created a function to do the Blast search between two sequence from the accession number. 
- query genome is "MN908947" wuhan virus
- subject to search is "MK062179" sars virus

There can be whole database search which takes some time to finish. Or you can search between two sequence. Here I searched only two sequence. In the result the its important to see the "e value". The lower the e value the better the match. Here our e value is zero. So, Wuhan virus is quite similar to the SARS virus. 


```python
import Bio
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from io import StringIO
```


```python
def blast_two_accession(a, b):
    for each in [a,b]:
        fasta = Entrez.efetch(db="nucleotide", id=each, rettype="fasta", retmode="text")
        query = SeqRecord(SeqIO.read(fasta, "fasta").seq, id= each)
        SeqIO.write(query, (each +".fasta"), "fasta")
        fasta.close()
    output = NcbiblastnCommandline(query= (a+".fasta"), subject= (b+ ".fasta"), outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))
    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:
            print('****Alignment****')
            print('sequence:' + str(alignment.title))
            print ('length:'+ str(alignment.length))
            print( 'e value:'+str(hsp.expect))  
```


```python
blast_two_accession("MN908947","MK062179" )
```

    '****Alignment****'
    'sequence:MK062179 MK062179 <unknown description>'
    'length:29727'
    'e value:0.0'
    '****Alignment****'
    'sequence:MK062179 MK062179 <unknown description>'
    'length:29727'
    'e value:0.0'
    '****Alignment****'
    'sequence:MK062179 MK062179 <unknown description>'
    'length:29727'
    'e value:0.0'
    '****Alignment****'
    'sequence:MK062179 MK062179 <unknown description>'
    'length:29727'
    'e value:0.0'



```python
blast_two_accession("MN908947","KF600628" )
```

## How to analyze it or what it tell us about the virus

Now its time to decode the genome of Corona Virus(SARS-CoV2). But first some basic understanding of molecular biology is required. Mainly what does the genome/code actually do in the cell/organism. Well they are the blueprint of protein. A sequence of nucleotides sometimes knowns as gene corresponds to certain amino acid and a chain of these amino acid forms into peptide chain also knows as protein thats necessary for certain function in organism. In short genome is the blueprint to produce protein in organism. To understand more about it you can watch [this video.](https://www.youtube.com/watch?v=wvTv8TqWC48). There are 20 different amino acids that can be found in living organisms. Various chained structure of these amino acids forms various protein structure. How the genomic code transfers into protein, [this video](https://www.youtube.com/watch?v=gG7uCskUOrA&t) is good to understand.

Transcription and Translation are two process in the pipeline of building protein in organisms. This [article from Khan academy ](https://www.khanacademy.org/science/high-school-biology/hs-molecular-genetics/hs-rna-and-protein-synthesis/a/intro-to-gene-expression-central-dogma) is good resource to understand the expression in molecular biology and the process of translation and transcription. Now we will try to translate our genome into decoded protein chain. For translation, we need a decoder to figure out the cipher in the genomic sequence. The decoder is called **codon**. The translation process can be understood by [this article](https://www.khanacademy.org/science/biology/gene-expression-central-dogma/translation-polypeptides/a/translation-overview). The codon table is follwoing. ![Codon](https://cdn.kastatic.org/ka-perseus-images/f5de6355003ee322782b26404ef0733a1d1a61b0.png)

Three letters in genome corresponds to a single amino acid out of 20. From the table we can get a codon collection that can be used to decode the genome of Corona Virus. Following table is taken from wikipedia. After some cleanup. A collection of amino acid and codon is generated to use in translation process. 


```python
codon = """Ala / A	GCU, GCC, GCA, GCG
Ile / I	AUU, AUC, AUA
Arg / R	CGU, CGC, CGA, CGG; AGA, AGG
Leu / L	CUU, CUC, CUA, CUG; UUA, UUG
Asn / N	AAU, AAC
Lys / K	AAA, AAG
Asp / D	GAU, GAC
Met / M	AUG
Phe / F	UUU, UUC
Cys / C	UGU, UGC
Pro / P	CCU, CCC, CCA, CCG
Gln / Q	CAA, CAG
Ser / S	UCU, UCC, UCA, UCG; AGU, AGC
Glu / E	GAA, GAG
Thr / T	ACU, ACC, ACA, ACG
Trp / W	UGG
Gly / G	GGU, GGC, GGA, GGG
Tyr / Y	UAU, UAC
His / H	CAU, CAC
Val / V	GUU, GUC, GUA, GUG
STOP	UAA, UGA, UAG""".strip()
```


```python
code = {}
for each in codon.split("\n"):
    a,b = each.split("\t")
    if "/" in a:
        a = a.split("/")[1]
    b = b.replace(" ", "").replace(";", ",").split(",")
    for bb in b:
        code[bb] = a
```

Now lets split the genome in tiplets. But where to start. We can start from first, second or third letter. As the first 6 letters in the genome are "attaaa". We can divide into three from "a", "t" or second "t". All of them can be a possible message. So we will take all three. these are called reading frame. And the first one is Open Reading Frame. You can learn more about reading frame [here](https://en.wikipedia.org/wiki/Reading_frame)


```python
amino_acid= []
for r in range(3):
    for i in range(r, len(covid_genome)-3, 3):
        gen = covid_genome.upper().replace("T", "U")
        amino_acid.append(code[gen[i:i+3]])
amino_acid = "".join(amino_acid).replace(" ", "" )
```


```python
print("First 1000 Amino Acid Sequence: "+ amino_acid[:1000])
```

    ('First 1000 Amino Acid Sequence: '
     'IKGLYLPRSTOPQTNQLSISCRSVLSTOPTNFKICVAVTRLHASTOPCTHAVSTOPLITNYCRSTOPQDTSNSSIFCRLLTVSSVLQPIISTSRFRPGVTERSTOPDGEPCPWFQRENTRPTQFACFTGSRRARTWLWRLRGGGLIRGTSTSSTOPRWHLWLSRSSTOPKRRFASTSTOPTALCVHQTFGCSNCTSWSCYGSTOPAGSRTRRHSVRSSTOPWSTOPDTWCPCPSCGRNTSGLPQGSSSSTOPERSTOPSTOPRSWWPSTOPLRRRSKVISTOPLRRRAWHSTOPSLSTOPRFSRKLEHSTOPTSTOPQWCYPSTOPTHASTOPASTOPRRGIHSLCRSTOPQLLWPSTOPWLPSSTOPVHSTOPRPSSTCWSTOPSFMHFVRTTGLYSTOPHSTOPEGCILLPSTOPTSTOPASTOPNCLVHGTFSTOPKELSTOPIADTFSTOPNSTOPIGKEISTOPHLQWGMSKFCISLKFHNQDYSTKGSTOPKEKASTOPWLYGSTOPNSICLSSCVTKSTOPMQPNVPFNSHEVSTOPSLWSTOPNFMADGRFCSTOPSHLRILWHSTOPEFDSTOPRRCHYLWLLTPKCCCSTOPNLLSSMSQFRSRTSTOPASTOPSCRIPSTOPSTOPIWLENHSSSTOPGWSHYCLWRLCVLLCWLPSTOPQVCLLGSTCSTOPRSTOPHRLSTOPPYRCCWRRFRRSSTOPSTOPQPSSTOPNTPKRESQHQYCWSTOPLSTOPTSTOPSTOPRDRHYFGIFFCFHKCFCGNCERFGLSTOPSIQTNCSTOPILWSTOPFSTOPSYKRKSSTOPKRCLEYWSTOPTEINTESSLCICIRGCSCCTINFLPHSSTOPNCSKFCACFTEGRYNNTRWNFTVFTETHSTOPCYDVHISTOPFGYSTOPQSSCNGLHYRWCCSVDFAVANSTOPHLWHCLSTOPKTQTRPSTOPLASTOPREVSTOPGRCRVSSTOPRRLGNCSTOPIYLNLCLSTOPNCRWT')


Now before spliting the peptide chains I got the gene bank data for the genes related to each type of protein in SARS-Cov-2 virus. I used efetch function to get the data from the ncbi nucleotide database. as rettype I gave gb as genebank format. Then I did some parsing to get only the peptide chain from the file. There are 10 peptide chain / gene in the genome. These will be used to compare our own spliting data. 


```python
genebank = Entrez.efetch(db="nucleotide", id="MN908947", rettype="gb", retmode="text").read()
```


```python
protein_gb = []
for each in genebank.split("3'UTR")[0].split("/translation=")[1:]:
    pep = each.split("gene", 1)[0].replace("\n", "").replace(" ", "").replace("\"", "")
    protein_gb.append(pep)
print("Number of gene: " + str(len(protein_gb)))
```

    'Number of gene: 10'


In the amino acid sequence there are 20 acid letters and "STOP" is there. But for splitting some information is important regarding the chain. 
- "Methionine/ Met" or the letter "M" is the start code for any peptide chain. Almost all the protein chain starts with "M"
- "STOP" is used to stop any peptides sequence. 
- Sometimes multiple sequence is added and combined to produce a new peptide sequence. 
- "M" Can be found in the middle too. 

I still have questions which "M" I should start. Thats also a learning scope. lets now just get on with it. For initial test we would take for each stop segments, start would be the first M and the second M. Then take all in our possible peptide chain.


```python
pep_ch = []
for each in amino_acid.split("STOP"):
    if "M" in each:
        pp = each.split("M",1)[1]
        pep_ch.append("M" + pp)
        if "M" in pp:
            pep_ch.append("M"+ pp.split("M", 1)[1])
        
# first filtered all those chain with less than 35 acids, then removed the duplicates with set() then sorted the list        
pep_ch = sorted(set((list(filter(lambda x: len(x) > 35, pep_ch)))), key=len)
```

From the data from the genebank we will match if those matches with our spliting. I have found that 9 out of 10 of the gene chain is matched with the gene bank data. Now lets see which one didn't matched. Only one chain didn't matched. probly because that is a combination of two peptide chain combined and mixed. 


```python
print("Matched genes: " + str(len(list(set(protein_gb).intersection(set(pep_ch))))))
print("Unmatched genes: " + str(list(set(protein_gb).difference(set(pep_ch)))))
```

    'Matched genes: 9'
    ('Unmatched genes: '
     "['MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGGAYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGGAPTKVTFGDDTVIEVQGYKSVNITFELDERIDKVLNEKCSAYTVELGTEVNEFACVVADAVIKTLQPVSELLTPLGIDLDEWSMATYYLFDESGEFKLASHMYCSFYPPDEDEEEGDCEEEEFEPSTQYEYGTEDDYQGKPLEFGATSAALQPEEEQEEDWLDDDSQQTVGQQDGSEDNQTTTIQTIVEVQPQLEMELTPVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLRVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQVEQKIAEIPKEEVKPFITESKPSVEQRKQDDKKIKACVEEVTTTLEETKFLTENLLLYIDINGNLHPDSATLVSDIDITFLKKDAPYIVGDVVQEGVLTAVVIPTKKAGGTTEMLAKALRKVPTDNYITTYPGQGLNGYTVEEAKTVLKKCKSAFYILPSIISNEKQEILGTVSWNLREMLAHAEETRKLMPVCVETKAIVSTIQRKYKGIKIQEGVVDYGARFYFYTSKTTVASLINTLNDLNETLVTMPLGYVTHGLNLEEAARYMRSLKVPATVSVSSPDAVTAYNGYLTSSSKTPEEHFIETISLAGSYKDWSYSGQSTQLGIEFLKRGDKSVYYTSNPTTFHLDGEVITFDNLKTLLSLREVRTIKVFTTVDNINLHTQVVDMSMTYGQQFGPTYLDGADVTKIKPHNSHEGKTFYVLPNDDTLRVEAFEYYHTTDPSFLGRYMSALNHTKKWKYPQVNGLTSIKWADNNCYLATALLTLQQIELKFNPPALQDAYYRARAGEAANFCALILAYCNKTVGELGDVRETMSYLFQHANLDSCKRVLNVVCKTCGQQQTTLKGVEAVMYMGTLSYEQFKKGVQIPCTCGKQATKYLVQQESPFVMMSAPPAQYELKHGTFTCASEYTGNYQCGHYKHITSKETLYCIDGALLTKSSEYKGPITDVFYKENSYTTTIKPVTYKLDGVVCTEIDPKLDNYYKKDNSYFTEQPIDLVPNQPYPNASFDNFKFVCDNIKFADDLNQLTGYKKPASRELKVTFFPDLNGDVVAIDYKHYTPSFKKGAKLLHKPIVWHVNNATNKATYKPNTWCIRCLWSTKPVETSNSFDVLKSEDAQGMDNLACEDLKPVSEEVVENPTIQKDVLECNVKTTEVVGDIILKPANNSLKITEEVGHTDLMAAYVDNSSLTIKKPNELSRVLGLKTLATHGLAAVNSVPWDTIANYAKPFLNKVVSTTTNIVTRCLNRVCTNYMPYFFTLLLQLCTFTRSTNSRIKASMPTTIAKNTVKSVGKFCLEASFNYLKSPNFSKLINIIIWFLLLSVCLGSLIYSTAALGVLMSNLGMPSYCTGYREGYLNSTNVTIATYCTGSIPCSVCLSGLDSLDTYPSLETIQITISSFKWDLTAFGLVAEWFLAYILFTRFFYVLGLAAIMQLFFSYFAVHFISNSWLMWLIINLVQMAPISAMVRMYIFFASFYYVWKSYVHVVDGCNSSTCMMCYKRNRATRVECTTIVNGVRRSFYVYANGGKGFCKLHNWNCVNCDTFCAGSTFISDEVARDLSLQFKRPINPTDQSSYIVDSVTVKNGSIHLYFDKAGQKTYERHSLSHFVNLDNLRANNTKGSLPINVIVFDGKSKCEESSAKSASVYYSQLMCQPILLLDQALVSDVGDSAEVAVKMFDAYVNTFSSTFNVPMEKLKTLVATAEAELAKNVSLDNVLSTFISAARQGFVDSDVETKDVVECLKLSHQSDIEVTGDSCNNYMLTYNKVENMTPRDLGACIDCSARHINAQVAKSHNIALIWNVKDFMSLSEQLRKQIRSAAKKNNLPFKLTCATTRQVVNVVTTKIALKGGKIVNNWLKQLIKVTLVFLFVAAIFYLITPVHVMSKHTDFSSEIIGYKAIDGGVTRDIASTDTCFANKHADFDTWFSQRGGSYTNDKACPLIAAVITREVGFVVPGLPGTILRTTNGDFLHFLPRVFSAVGNICYTPSKLIEYTDFATSACVLAAECTIFKDASGKPVPYCYDTNVLEGSVAYESLRPDTRYVLMDGSIIQFPNTYLEGSVRVVTTFDSEYCRHGTCERSEAGVCVSTSGRWVLNNDYYRSLPGVFCGVDAVNLLTNMFTPLIQPIGALDISASIVAGGIVAIVVTCLAYYFMRFRRAFGEYSHVVAFNTLLFLMSFTVLCLTPVYSFLPGVYSVIYLYLTFYLTNDVSFLAHIQWMVMFTPLVPFWITIAYIICISTKHFYWFFSNYLKRRVVFNGVSFSTFEEAALCTFLLNKEMYLKLRSDVLLPLTQYNRYLALYNKYKYFSGAMDTTSYREAACCHLAKALNDFSNSGSDVLYQPPQTSITSAVLQSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQSAVKRTIKGTHHWLLLTILTSLLVLVQSTQWSLFFFLYENAFLPFAMGIIAMSAFAMMFVKHKHAFLCLFLLPSLATVAYFNMVYMPASWVMRIMTWLDMVDTSLSGFKLKDCVMYASAVVLLILMTARTVYDDGARRVWTLMNVLTLVYKVYYGNALDQAISMWALIISVTSNYSGVVTTVMFLARGIVFMCVEYCPIFFITGNTLQCIMLVYCFLGYFCTCYFGLFCLLNRYFRLTLGVYDYLVSTQEFRYMNSQGLLPPKNSIDAFKLNIKLLGVGGKPCIKVATVQSKMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQGAVDINKLCEEMLDNRATLQAIASEFSSLPSYAAFATAQEAYEQAVANGDSEVVLKKLKKSLNVAKSEFDRDAAMQRKLEKMADQAMTQMYKQARSEDKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIVQLSEISMDNSPNLAWPLIVTALRANSAVKLQNNELSPVALRQMSCAAGTTQTACTDDNALAYYNTTKGGRFVLALLSDLQDLKWARFPKSDGTGTIYTELEPPCRFVTDTPKGPKVKYLYFIKGLNNLNRGMVLGSLAATVRLQAGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPEANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVCGMWKGYGCSCDQLREPMLQSADAQSFLNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLDNQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVNCLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEGSSVELKHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGKARLYYDSMSYEDQDALFAYTKRNVIPTITQMNLKYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRGATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLANECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYECLYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNSTYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTETDLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQEYADVFHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQAVGACVLCNSQTSLRCGACIRRPFLCCKCCYDHVISTSHKLVLSVNPYVCNAPGCDVTDVTQLYLGGMSYYCKSHKPPISFPLCANGQVFGLYKNTCVGSDNVTDFNAIATCDWTNAGDYILANTCTERLKLFAAETLKATEETFKLSYGIATVREVLSDRELHLSWEVGKPRPPLNRNYVFTGYRVTKNSKVQIGEYTFEKGDYGDAVVYRGTTTYKLNVGDYFVLTSHTVMPLSAPTLVPQEHYVRITGLYPTLNISDEFSSNVANYQKVGMQKYSTLQGPPGTGKSHFAIGLALYYPSARIVYTACSHAAVDALCEKALKYLPIDKCSRIIPARARVECFDKFKVNSTLEQYVFCTVNALPETTADIVVFDEISMATNYDLSVVNARLRAKHYVYIGDPAQLPAPRTLLTKGTLEPEYFNSVCRLMKTIGPDMFLGTCRRCPAEIVDTVSALVYDNKLKAHKDKSAQCFKMFYKGVITHDVSSAINRPQIGVVREFLTRNPAWRKAVFISPYNSQNAVASKILGLPTQTVDSSQGSEYDYVIFTQTTETAHSCNVNRFNVAITRAKVGILCIMSDRDLYDKLQFTSLEIPRRNVATLQAENVTGLFKDCSKVITGLHPTQAPTHLSVDTKFKTEGLCVDIPGIPKDMTYRRLISMMGFKMNYQVNGYPNMFITREEAIRHVRAWIGFDVEGCHATREAVGTNLPLQLGFSTGVNLVAVPTGYVDTPNNTDFSRVSAKPPPGDQFKHLIPLMYKGLPWNVVRIKIVQMLSDTLKNLSDRVVFVLWAHGFELTSMKYFVKIGPERTCCLCDRRATCFSTASDTYACWHHSIGFDYVYNPFMIDVQQWGFTGNLQSNHDLYCQVHGNAHVASCDAIMTRCLAVHECFVKRVDWTIEYPIIGDELKINAACRKVQHMVVKAALLADKFPVLHDIGNPKAIKCVPQADVEWKFYDAQPCSDKAYKIEELFYSYATHSDKFTDGVCLFWNCNVDRYPANSIVCRFDTRVLSNLNLPGCDGGSLYVNKHAFHTPAFDKSAFVNLKQLPFFYYSDSPCESHGKQVVSDIDYVPLKSATCITRCNLGGAVCRHHANEYRLYLDAYNMMISAGFSLWVYKQFDTYNLWNTFTRLQSLENVAFNVVNKGHFDGQQGEVPVSIINNTVYTKVDGVDVELFENKTTLPVNVAFELWAKRNIKPVPEVKILNNLGVDIAANTVIWDYKRDAPAHISTIGVCSMTDIAKKPTETICAPLTVFFDGRVDGQVDLFRNARNGVLITEGSVKGLQPSVGPKQASLNGVTLIGEAVKTQFNYYKKVDGVVQQLPETYFTQSRNLQEFKPRSQMEIDFLELAMDEFIERYKLEGYAFEHIVYGDFSHSQLGGLHLLIGLAKRFKESPFELEDFIPMDSTVKNYFITDAQTGSSKCVCSVIDLLLDDFVEIIKSQDLSVVSKVVKVTIDYTEISFMLWCKDGHVETFYPKLQSSQAWQPGVAMPNLYKMQRMLLEKCDLQNYGDSATLPKGIMMNVAKYTQLCQYLNTLTLAVPYNMRVIHFGAGSDKGVAPGTAVLRQWLPTGTLLVDSDLNDFVSDADSTLIGDCATVHTANKWDLIISDMYDPKTKNVTKENDSKEGFFTYICGFIQQKLALGGSVAIKITEHSWNADLYKLMGHFAWWTAFVTNVNASSSEAFLIGCNYLGKPREQIDGYVMHANYIFWRNTNPIQLSSYSLFDMSKFPLKLRGTAVMSLKEGQINDMILSLLSKGRLIIRENNRVVISSDVLVNN']")



```python
matched_gene = list(set(protein_gb).intersection(set(pep_ch)))
```

### Now Lets see 3D view of each of the protein structure

One of the cool thing I discovered that you can easily show the 3d structure of each protein sequence that I just translated and filtered. For this theres a web application from Swiss Institute of Bioinformatics: [Swiss Model](https://swissmodel.expasy.org/). All you need to do is:
- copy the protein structure you want to build 3d model. 
- Go to [Swiss Model](https://swissmodel.expasy.org/) site
- click **start modeling**
- paste the sequence you copied and **validate** for errors
- Click **build model** and **VOILA!**

Now download the pdb model file from the left side of the model result dropdown. pdb is protein database format. used to store protein structure data. You can use any protein only limit is the length should be less than 5000. 


```python
# first protein structure to copy
print(matched_gene[0])
```

    'MKFLVFLGIITTVAAFHQECSLQSCTQHQPYVVDDPCPIHFYSKWYIRVGARKSAPLIELCVDEAGSKSPIQYIDIGNYTVSCLPFTINCQEPKLGSLVVRCSFYEDFLEYHDVRVVLDFI'


To view the pdb file. There a package called NGL. It has a package that can view pdb file directly in jupyter notebook. I installed the required package from [here](https://github.com/arose/nglview). You need to jupyterlab installed from to use this. you can install it by ***"pip install jupyterlab"***


```python
import nglview
# view the downloaded model file from the swiss model application
show = nglview.show_file("model_01.pdb")
show
```


    NGLWidget()



```python
show.render_image()
show._display_image()
```




![png](how-corona-virus-got-me-into-bioinformatics_files/how-corona-virus-got-me-into-bioinformatics_52_0.png)




```python
show.add_surface(selection="protein", color="sstruc")
show.render_image()
show._display_image()
```




![png](how-corona-virus-got-me-into-bioinformatics_files/how-corona-virus-got-me-into-bioinformatics_53_0.png)



#### Its soooo coooool !!! 


```python

```
