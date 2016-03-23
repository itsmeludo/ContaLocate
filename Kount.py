#!/usr/bin/env python3

### Author: Ludovic V. Mallet, PhD
### 2016.03.22
### licence: GPLv3
### Version: Alpha.1
### Garanteed with misatkes. <- Including this one.


# This program computes oligonucleotide profiles (microcomposition) for sequences given as arguments.
# Depending on the set of parameters used, the program guess what is possible to do and will perform accordingly. Use the proper (and minimal) parameter set you need to achieve what you want.
# Luckily, the program will tell you what it does, and how you can change its behaviour by providing him with supplementary parameters.
# See the help by calling the program without any argument.


#dependencies: 
  #Biopython
  #numpy
  #cython

#as root/admin if wanted installed globally
#aptitude/yum install python3-dev python3-setuptools


#easy_install -U setuptools
#easy_install3 -U setuptools

#pip install biopython
#pip3 install biopython 

#pip install cython
#pip3 install cython

#pip install numpy
#pip3 install numpy



from optparse import OptionParser
import os
import re
import math
import numpy

from Bio import SeqIO
from Bio.Seq import Seq

from collections import Counter
import multiprocessing

numpy.seterr(divide='ignore', invalid='ignore')
minimum_number_of_windows_per_fasta_entry_to_use_multiple_cpu_for_this_entry=20
word_order =( "CCCC", "GCCC", "CGCC", "GGCC", "CCGC", "GCGC", "CGGC", "GGGC", "CCCG", "GCCG", "CGCG", "GGCG", "CCGG", "GCGG", "CGGG", "GGGG", "ACCC", "TCCC", "AGCC", "TGCC", "ACGC", "TCGC", "AGGC", "TGGC", "ACCG", "TCCG", "AGCG", "TGCG", "ACGG", "TCGG", "AGGG", "TGGG", "CACC", "GACC", "CTCC", "GTCC", "CAGC", "GAGC", "CTGC", "GTGC", "CACG", "GACG", "CTCG", "GTCG", "CAGG", "GAGG", "CTGG", "GTGG", "AACC", "TACC", "ATCC", "TTCC", "AAGC", "TAGC", "ATGC", "TTGC", "AACG", "TACG", "ATCG", "TTCG", "AAGG", "TAGG", "ATGG", "TTGG", "CCAC", "GCAC", "CGAC", "GGAC", "CCTC", "GCTC", "CGTC", "GGTC", "CCAG", "GCAG", "CGAG", "GGAG", "CCTG", "GCTG", "CGTG", "GGTG", "ACAC", "TCAC", "AGAC", "TGAC", "ACTC", "TCTC", "AGTC", "TGTC", "ACAG", "TCAG", "AGAG", "TGAG", "ACTG", "TCTG", "AGTG", "TGTG", "CAAC", "GAAC", "CTAC", "GTAC", "CATC", "GATC", "CTTC", "GTTC", "CAAG", "GAAG", "CTAG", "GTAG", "CATG", "GATG", "CTTG", "GTTG", "AAAC", "TAAC", "ATAC", "TTAC", "AATC", "TATC", "ATTC", "TTTC", "AAAG", "TAAG", "ATAG", "TTAG", "AATG", "TATG", "ATTG", "TTTG", "CCCA", "GCCA", "CGCA", "GGCA", "CCGA", "GCGA", "CGGA", "GGGA", "CCCT", "GCCT", "CGCT", "GGCT", "CCGT", "GCGT", "CGGT", "GGGT", "ACCA", "TCCA", "AGCA", "TGCA", "ACGA", "TCGA", "AGGA", "TGGA", "ACCT", "TCCT", "AGCT", "TGCT", "ACGT", "TCGT", "AGGT", "TGGT", "CACA", "GACA", "CTCA", "GTCA", "CAGA", "GAGA", "CTGA", "GTGA", "CACT", "GACT", "CTCT", "GTCT", "CAGT", "GAGT", "CTGT", "GTGT", "AACA", "TACA", "ATCA", "TTCA", "AAGA", "TAGA", "ATGA", "TTGA", "AACT", "TACT", "ATCT", "TTCT", "AAGT", "TAGT", "ATGT", "TTGT", "CCAA", "GCAA", "CGAA", "GGAA", "CCTA", "GCTA", "CGTA", "GGTA", "CCAT", "GCAT", "CGAT", "GGAT", "CCTT", "GCTT", "CGTT", "GGTT", "ACAA", "TCAA", "AGAA", "TGAA", "ACTA", "TCTA", "AGTA", "TGTA", "ACAT", "TCAT", "AGAT", "TGAT", "ACTT", "TCTT", "AGTT", "TGTT", "CAAA", "GAAA", "CTAA", "GTAA", "CATA", "GATA", "CTTA", "GTTA", "CAAT", "GAAT", "CTAT", "GTAT", "CATT", "GATT", "CTTT", "GTTT", "AAAA", "TAAA", "ATAA", "TTAA", "AATA", "TATA", "ATTA", "TTTA", "AAAT", "TAAT", "ATAT", "TTAT", "AATT", "TATT", "ATTT", "TTTT")


def KL(a,b):
  #with numpy.errstate(invalid='ignore'):
  d = a * numpy.log(a/b)
  d[numpy.isnan(d)]=0 
  d[numpy.isinf(d)]=0
  return (numpy.sum(d))*10000

def Eucl(a,b):
  #with numpy.errstate(invalid='ignore'):
  d = pow(a-b,2)
  d[numpy.isnan(d)]=0
  d[numpy.isinf(d)]=0
  return numpy.sqrt(numpy.sum(d))*10000

def JSD(a,b):
  #with numpy.errstate(invalid='ignore'):
  h = (a + b)/2
  d = (KL(a,h)/2)+(KL(b,h)/2)
  return d
    
   
def vector_to_matrix(profile):
  return list((zip(*(iter(profile),)*int(math.sqrt(len(profile))))))
  
  
  
def chunkitize(liste, chunks):
  out=list()
  chunk_size= int((len(liste) / float(chunks)))
  for i in range(0,chunks):
    out.append(liste[chunk_size*i:(chunk_size*i+chunk_size)if(i!=chunks-1)else len(liste)])
  return out

#chunkitize([1, 2, 3, 4, 5], 2)




def select_strand (seq,strand="both"):
  Bioseq_record=Seq(seq)
  if(strand =="both"):
    return str(str(seq)+str(Bioseq_record.reverse_complement())).upper()
  elif(strand == "minus"):
    return str(Bioseq_record.reverse_complement()).upper()
  elif(strand == "plus"):
    return str(seq).upper()
  else:
    return 1
   
   
def frequency (seq,ksize=4,strand="both"):
  seq=select_strand(seq,strand)
  seq_words=list()
  d=dict()
  for s in re.split('[^ACGTacgt]+',seq): #excludes what is not a known characterised nucleotide 
    seq_letters=list(s)
    ##print(len(s))
    if (len(seq_letters) >=ksize):
      # launch k-1 times the word generation with 1 nucleotide shift every iteration to have overlapping words.
      for i in range(ksize-1):
        #generate words from seq_letters
        words=list((zip(*(iter(seq_letters),)*ksize)))
        seq_words.extend(list(map(''.join,words))) # adds the words for this subsequence and frame to the total list of words
        seq_letters.pop(0) # shift one to compute overlapping words at the next iteration
  c = Counter(seq_words)
  word_count=sum(c.values())
  if(word_count > 0):
    ret=list()
    for w in word_order:
      if(c.get(w) != None):
        ret.append(c.get(w)/word_count)
      else:
        ret.append(0)
    return ret
  else:
    return 1


def parallel_subwin_dist(args):
  res=list()
  for window,start,stop,seq_id,mcp_comparison,windows_size,windows_step,ksize,dist,position,contig_size in args:
    #print(window,start,stop,mcp_comparison,windows_size,windows_step,ksize,dist,position)
    #to avoid border effects being a problem in the code, we use only the simple formula start+windows_size/2 (+/-) windows_step/2 to find the significant center part of the windows. when several windows overlap, this centered part, long as the windows step is the most representative of the windows, not representing as much other part of this window that are overlapped by other windows. BUT: this simple formula has border effects, so we manually correct the start of the first window and the stop of the last window to match the contig borders.
    if(start == (windows_size/2-windows_step/2)):
      displayed_start=1
    else:
       displayed_start=start
    if(stop-windows_step/2+windows_size/2>=contig_size-windows_step and stop-windows_step/2+windows_size/2<=contig_size):
      displayed_stop=contig_size
    else:
      displayed_stop=stop
    if((window.count('N')/windows_size) <= float(options.n_max_freq_in_windows)):
      if(position==False):
        res.append(globals()[dist](mcp_comparison,numpy.array(frequency(seq=window,ksize=ksize))))
      else:
        res.append([seq_id,displayed_start,displayed_stop,globals()[dist](mcp_comparison,numpy.array(frequency(seq=window,ksize=ksize)))])
    else:
      if(position==False):
        res.append(numpy.nan)
      else:
        res.append([seq_id,displayed_start,displayed_stop,numpy.nan])
  return(res)

#def subwin_dist(args):
  #window,start,stop,mcp_comparison,windows_size,windows_step,ksize,dist,position = args
  ##print(threading.currentThread().getName(), 'Starting')
  ##print("Process "+str(start)+"\n")
  #if((window.count('N')/windows_size) <= float(options.n_max_freq_in_windows)):
    #if(position==False):
      #res=globals()[dist](mcp_comparison,numpy.array(frequency(seq=window,ksize=ksize)))
    #else:
      #res=[start,stop,globals()[dist](mcp_comparison,numpy.array(frequency(seq=window,ksize=ksize)))]
  #else:
    #if(position==False):
      #res=numpy.nan
    #else:
      #res=[start,stop,numpy.nan]
  #return(res)

def sliding_windows_distances(seq,mcp_comparison,seq_id,dist="KL",windows_size=5000,windows_step=500,ksize=4,position=False):
  #seq=str(Bioseq_record.seq)
  ret=list()
  if(len(seq)<windows_size): #only enough to compute one window, no sliding,
    if((seq.count('N')/windows_size) <= float(options.n_max_freq_in_windows)):
      if(position==False):
        ret.append( globals()[dist](mcp_comparison,numpy.array(frequency(seq=seq,ksize=ksize))))
      else:
        ret.append([seq_id,0,int(len(seq)),globals()[dist](mcp_comparison,frequency(seq=seq,ksize=ksize))])
    else:
      if(position==False):
        ret.append(numpy.nan)
      else:
        ret.append([seq_id,0,int(len(seq)),numpy.nan])
  elif(len(seq)<minimum_number_of_windows_per_fasta_entry_to_use_multiple_cpu_for_this_entry*windows_step): #not many windows in this contig, so better launching it in serial rather than in parallel
    tmp_out=list()
    if(position==False):
      for s in range(0,len(seq)-windows_size,windows_step):
        window=seq[s:s+windows_size]
        if((window.count('N')/windows_size) <= float(options.n_max_freq_in_windows)):
          tmp_out.append( globals()[dist](mcp_comparison,numpy.array(frequency(seq=window,ksize=ksize))))
        else:
          tmp_out.append(numpy.nan)
    else:
      for s in range(0,len(seq)-windows_size,windows_step):
        if(s==0):# to avoid border effects being a problem in the code, we use only the simple formula start+windows_size/2 (+/-) windows_step/2 to find the significant center part of the windows. when several windows overlap, this centered part, long as the windows step is the most representative of the windows, not representing as much other part of this window that are overlapped by other windows. BUT: this simple formula has border effects, so we manually correct the start of the first window and the stop of the last window to match the contig borders.
          displayed_start=1
        else:
          displayed_start=int(s+windows_size/2-windows_step/2)

        if(s==len(seq)-windows_size):
          displayed_stop=len(seq)
        else:
          displayed_stop=int(s+windows_size/2+windows_step/2)

        window=seq[s:s+windows_size]
        if((window.count('N')/windows_size) <= float(options.n_max_freq_in_windows)):
          tmp_out.append([seq_id,displayed_start,displayed_stop,globals()[dist](mcp_comparison,frequency(seq=window,ksize=ksize))])
        else:
          tmp_out.append([seq_id,displayed_start,displayed_stop,numpy.nan])
    for i in tmp_out:
      ret.append(i)
  else:
    args = [(seq[s:s+windows_size],int(s+windows_size/2-windows_step/2),int(s+windows_size/2+windows_step/2),seq_id,mcp_comparison,windows_size, windows_step, ksize,dist, position,len(seq)) for s in range(0,len(seq)-windows_size,windows_step)]
    parallel_args_set=chunkitize(args,options.threads_max) 
    pool= multiprocessing.Pool(processes=options.threads_max)
    res=pool.map(parallel_subwin_dist, parallel_args_set)
    pool.close()
    pool.join()
    for i in res:
      for g in i:
        ret.append(g)
  return ret

#def subwin_feq(window,windows_size=5000,windows_step=500,ksize=4):
  #if((window.count('N')/windows_size)>= options.n_max_freq_in_windows):
    #q.put(numpy.array(numpy.nan))
  #else:
    #q.put(numpy.array(frequency(seq=window,ksize=ksize)))

#def sliding_windows_frequencies(seq,windows_size=5000,windows_step=500,ksize=4,position=False):
  ##seq=str(Bioseq_record.seq)
  #if(len(seq)<windows_size): #only enough to compute one window, no sliding,
    #return numpy.array(frequency(seq=seq,ksize=ksize))
  #else:
    #ret = list()
    #q = Queue()
    #with concurrent.futures.ThreadPoolExecutor(max_workers=options.threads_max) as executor:
      #for s in range(0,len(seq)-windows_size,windows_step):
        ##print(s,s+windows_size)
        #window=seq[s:s+windows_size]
        #executor.submit(subwin_feq, window,windows_size,windows_step,ksize)
    #while not q.empty():
      #ret.append(q.get())
    #return ret

#def sliding_windows_target_regions(seq, windows_size=5000, windows_step=500):
  ##seq=str(Bioseq_record.seq)
  #if(len(seq)<windows_size): #only enough to compute one window, no sliding,
    #return list(0,len(seq))
  #else:
    #ret=list()
    #for s in range(0,len(seq)-windows_size,windows_step):
      ##print(s,s+windows_size)
      #window=[int(s+windows_size/2-windows_step/2),int(s+windows_size/2+windows_step/2)]
      #ret.append(window)
  #return ret


if __name__ == "__main__":

  Utilisation = "%prog [-i FILE] [options]"
  parser = OptionParser(usage = Utilisation)
  parser.add_option("-i","--assembly", dest = "genome", help = "multifasta of the genome assembly")
  parser.add_option("-c","--conta", dest = "conta", help = "multifasta of the contaminant species training set")
  parser.add_option("-r","--host", dest = "host", help = "optional multifasta of the host species training set")
  parser.add_option("-n","--n_max_freq_in_windows", type = "float", dest = "n_max_freq_in_windows", default = "0.4", help = "maximum proportion of N tolerated in a window to compute the microcomposition anyway [0~1]. Too much 'N's will reduce the number of kmer counts and will artificially coarse the resolution of the frequencies. If your assembly contains many stretches of 'N's, consider rising this parameter and shortening the windows step in order to allow computation and signal in the output, it might cost you computational time though. Windows which do not meet this criteria will be affected a distance of 'nan'")
  parser.add_option("-k","--lgMot", dest = "k", type = "int", default = 4, help = "word wise / kmer lenght / k [default:%default]")
  parser.add_option("-w","--windows_size", dest ="windows_size", type = "int", help = "Sliding windows size (bp)[default:%default]")
  parser.add_option("-t","--windows_step", dest ="windows_step", type = "int", help = "Sliding windows step size(bp)[default:%default]")
  parser.add_option("-s","--strand", default = "double", help = "strand used to compute microcomposition. leading, lagging ou double [default:%default]")
  parser.add_option("-d","--distance", dest = "dist", default = "JSD", help = "méthode de distance entre 2 signatures : KL: Kullback-Leibler, Eucl : Euclidienne[default:%default], JSD : Jensen-Shannon divergence")
  parser.add_option("-u","--cpu", dest = "threads_max", type = "int", default = 4, help = "how maany threads to use for windows microcomposition computation[default:%default]")

  #parser.add_option("-q","--qqconque", action = "store_true", dest="afaire")

  (options,argument) = parser.parse_args()




  if(not options.genome):
    parser.error("An input fasta file (-i ) is mandatory")
    exit()

  if(not options.conta):
    if(not options.windows_size and not options.windows_step):
      #parser.error("An input fasta file (-i ) is mandatory")
      print("A genome was provided but no training set to learn the contamination profile (-c). You didn't provide sliding window parameters (-w -t) so I will just compute the signature of the whole genome ")
      whole_seq=Seq("")
      for record in SeqIO.parse(options.genome, "fasta"):
        whole_seq.seq=str(whole_seq)+"N"+str(record.seq)
      genome=numpy.array(frequency(seq=str(whole_seq.seq),ksize=option.k if options.k else 4))
      f = open(str(os.path.basename(options.genome))+".microcomposition.mat", 'w')
      f.write(str(vector_to_matrix(genome)))
      f.close()
      
      
      
    elif(options.windows_size or options.windows_step):
      print("A genome was provided but no training set to learn the contamination profile (-c). You provided sliding window parameters (either -w -t) so I will compute the genome microcomposition signature and I will compute the distance of every window microcomposition to that of the whole genome ")
      whole_seq=Seq("")
      for record in SeqIO.parse(options.genome, "fasta"):
        whole_seq.seq=str(whole_seq)+"N"+str(record.seq)
      genome=numpy.array(frequency(seq=str(whole_seq.seq),ksize=options.k if options.k else 4))
      
      #frq_wins=list()
      f = open(str(os.path.basename(options.genome))+".mcp_windows_vs_whole_"+options.dist+".dist", 'w')
      for record in SeqIO.parse(options.genome, "fasta"):
        windows_distances=sliding_windows_distances(str(record.seq),seq_id=record.id,mcp_comparison=genome,position=True,dist=options.dist, windows_size=options.windows_size if options.windows_size else 5000, windows_step=options.windows_step if options.windows_step else 500, ksize=options.k if options.k else 4)
        for t in windows_distances:
          f.write(str("\t".join(map(str,t)))+"\n")
      f.close()

  if(options.genome and options.conta):
    print("A genome and a training set to learn the contamination profile was provided so I will compute the microcomposition signature of the genome (-i) or host training set if provided (-r), that of the contamination training set (-c) and those of the sliding windows along the assembly (-i). with all of that I will output the distances of the microcomposition of every windows compared to both microcompositions of the contaminant and the host genome.")
    whole_seq=Seq("")
    if(options.host):
      print("using the provided host training set to learn the microcomposition")
      target=options.host
    else:
      print("No host training set provided, I will use the whole genome to learn the microcomposition")
      target=options.genome
    whole_seq=Seq("")
    for record in SeqIO.parse(target, "fasta"):
      whole_seq.seq=str(whole_seq)+"N"+str(record.seq)
    genome=numpy.array(frequency(seq=str(whole_seq.seq),ksize=options.k if options.k else 4))

    whole_conta=Seq("")
    for record in SeqIO.parse(options.conta, "fasta"):
      whole_conta.seq=str(whole_conta)+"N"+str(record.seq)
    conta=numpy.array(frequency(seq=str(whole_conta.seq),ksize=options.k if options.k else 4))

    fname=str(os.path.basename(options.genome))+".mcp_hostwindows_vs_"
    if options.host:
      fname=fname+ "host_"+str(os.path.basename(options.host))+"_"+options.dist+".dist"
    else:
      fname=fname+"wholegenome_"+options.dist+".dist"
    #print(fname+" totot")
    f = open(fname, 'w')
    for record in SeqIO.parse(options.genome, "fasta"):
      #f.write('>'+str(record.id)+"\n")
      #frq_wins = sliding_windows_frequencies(str(record.seq), windows_size=options.windows_size if options.windows_size else 5000, windows_step=options.windows_step if options.windows_step else 500, ksize=options.k if options.k else 4)
      windows_distances=sliding_windows_distances(str(record.seq),seq_id=record.id,mcp_comparison=genome,position=True,dist=options.dist, windows_size=options.windows_size if options.windows_size else 5000, windows_step=options.windows_step if options.windows_step else 500, ksize=options.k if options.k else 4)
      #windows_distances=list()
      for t in windows_distances:
        f.write(str("\t".join(map(str,t)))+"\n")
    f.close()
    
    fname=str(os.path.basename(options.genome))+".mcp_hostwindows_vs_"+"conta_"+str(os.path.basename(options.conta))+"_"+options.dist+".dist"
    #print(fname+" totowdawdt")
    f = open(fname, 'w')
    for record in SeqIO.parse(options.genome, "fasta"):
      #f.write('>'+str(record.id)+"\n")
      #frq_wins = sliding_windows_frequencies(str(record.seq), windows_size=options.windows_size if options.windows_size else 5000, windows_step=options.windows_step if options.windows_step else 500, ksize=options.k if options.k else 4)
      windows_distances=sliding_windows_distances(str(record.seq),seq_id=record.id,mcp_comparison=conta,position=True,dist=options.dist, windows_size=options.windows_size if options.windows_size else 5000, windows_step=options.windows_step if options.windows_step else 500, ksize=options.k if options.k else 4)
      #windows_distances=list()
      for t in windows_distances:
        f.write(str("\t".join(map(str,t)))+"\n")
    f.close()


































#for record in SeqIO.parse("host_test.fa", "fasta"):
  #host=numpy.array(frequency(seq=str(record.seq)))
  #frq_wins=sliding_windows_frequencies(str(record.seq))
 
  
  
#for record in SeqIO.parse("conta_test.fa", "fasta"):
  #conta=numpy.array(frequency(seq=str(record.seq)))
  ##vector_to_matrix(frequency(seq))
  
  
#windows_distances=list()
#for profile in frq_wins:
  #windows_distances.append(JSD(conta,profile))
  
#print(windows_distances)
