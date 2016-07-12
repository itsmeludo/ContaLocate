#!/usr/bin/env python3
"""
This program computes oligonucleotide profiles (microcomposition) for sequences given as arguments.
Depending on the set of parameters used, the program guess what is possible to do and will perform accordingly. Use the proper (and minimal) parameter set you need to achieve what you want.
Luckily, the program will tell you what it does, and how you can change its behaviour by providing him with supplementary parameters.
See the help by calling the program without any argument.


dependencies: 
  Biopython
  numpy
  cython

as root/admin if wanted installed globally
aptitude/yum install python3-dev python3-setuptools


easy_install -U setuptools
easy_install3 -U setuptools

pip install biopython
pip3 install biopython 

pip install cython
pip3 install cython

pip install numpy
pip3 install numpy
"""

__author__ = "Ludovic V. Mallet, PhD"
__copyright__ = ""
__date__ = "2016.03.22"
__licence__ = "GPLv3"
__version__ = "0.1"
__status__ = "alpha"
__email__ = ""

import os, re, math, sys, argparse
from collections import Counter
import multiprocessing
from itertools import product

import numpy
from Bio import SeqIO
from Bio.Seq import Seq


numpy.seterr(divide='ignore', invalid='ignore')
#TBF: wtf way to long ! (-:
#minimum_number_of_windows_per_fasta_entry_to_use_multiple_cpu_for_this_entry=20
#LM Fine by me, this very variable makes the implementation very much complicated and triggers the multicore or serial mode depending on the number of windows in a contig. therefore I wanted to be explicit on what it does. But I guess this mere sentence explqins it now ^^.
min_nb_w_per_fasta_for_mul_cpu=20

# TBF too ugly on one line 
#LM Actually this only needed to represent the matrix according to the spatialisation of the chaos game representation. It only works with 4-mers and I'm too lazy to port the algo to make it generic because it is not used in the program. In PhylOligo, I generalised with "word_universe" at the loss of the CGR representation (order is wrong to represent as a matrix with kmer spacial coherence).
#so, I comment it, it is only  necessary for vector_to_matrix() to be consistent with the GOHTAM database which is not connected here anyway so that the matrix can be sought against the DB instead of having to recompute the profile first. It would have been great, but recomputing profile within GOHTAM is cheaper in brain time, and out of scope for this Contalocate in this version.
#I changed word_order to word_universe in frequency() as I implemented in PhylOligo afterwards.
#word_order =( "CCCC", "GCCC", "CGCC", "GGCC", "CCGC", "GCGC", "CGGC", "GGGC", 
             #"CCCG", "GCCG", "CGCG", "GGCG", "CCGG", "GCGG", "CGGG", "GGGG", 
             #"ACCC", "TCCC", "AGCC", "TGCC", "ACGC", "TCGC", "AGGC", "TGGC", 
             #"ACCG", "TCCG", "AGCG", "TGCG", "ACGG", "TCGG", "AGGG", "TGGG", 
             #"CACC", "GACC", "CTCC", "GTCC", "CAGC", "GAGC", "CTGC", "GTGC", 
             #"CACG", "GACG", "CTCG", "GTCG", "CAGG", "GAGG", "CTGG", "GTGG", 
             #"AACC", "TACC", "ATCC", "TTCC", "AAGC", "TAGC", "ATGC", "TTGC", 
             #"AACG", "TACG", "ATCG", "TTCG", "AAGG", "TAGG", "ATGG", "TTGG", 
             #"CCAC", "GCAC", "CGAC", "GGAC", "CCTC", "GCTC", "CGTC", "GGTC", 
             #"CCAG", "GCAG", "CGAG", "GGAG", "CCTG", "GCTG", "CGTG", "GGTG", 
             #"ACAC", "TCAC", "AGAC", "TGAC", "ACTC", "TCTC", "AGTC", "TGTC", 
             #"ACAG", "TCAG", "AGAG", "TGAG", "ACTG", "TCTG", "AGTG", "TGTG", 
             #"CAAC", "GAAC", "CTAC", "GTAC", "CATC", "GATC", "CTTC", "GTTC", 
             #"CAAG", "GAAG", "CTAG", "GTAG", "CATG", "GATG", "CTTG", "GTTG", 
             #"AAAC", "TAAC", "ATAC", "TTAC", "AATC", "TATC", "ATTC", "TTTC", 
             #"AAAG", "TAAG", "ATAG", "TTAG", "AATG", "TATG", "ATTG", "TTTG", 
             #"CCCA", "GCCA", "CGCA", "GGCA", "CCGA", "GCGA", "CGGA", "GGGA", 
             #"CCCT", "GCCT", "CGCT", "GGCT", "CCGT", "GCGT", "CGGT", "GGGT", 
             #"ACCA", "TCCA", "AGCA", "TGCA", "ACGA", "TCGA", "AGGA", "TGGA", 
             #"ACCT", "TCCT", "AGCT", "TGCT", "ACGT", "TCGT", "AGGT", "TGGT", 
             #"CACA", "GACA", "CTCA", "GTCA", "CAGA", "GAGA", "CTGA", "GTGA", 
             #"CACT", "GACT", "CTCT", "GTCT", "CAGT", "GAGT", "CTGT", "GTGT", 
             #"AACA", "TACA", "ATCA", "TTCA", "AAGA", "TAGA", "ATGA", "TTGA", 
             #"AACT", "TACT", "ATCT", "TTCT", "AAGT", "TAGT", "ATGT", "TTGT", 
             #"CCAA", "GCAA", "CGAA", "GGAA", "CCTA", "GCTA", "CGTA", "GGTA", 
             #"CCAT", "GCAT", "CGAT", "GGAT", "CCTT", "GCTT", "CGTT", "GGTT", 
             #"ACAA", "TCAA", "AGAA", "TGAA", "ACTA", "TCTA", "AGTA", "TGTA", 
             #"ACAT", "TCAT", "AGAT", "TGAT", "ACTT", "TCTT", "AGTT", "TGTT", 
             #"CAAA", "GAAA", "CTAA", "GTAA", "CATA", "GATA", "CTTA", "GTTA", 
             #"CAAT", "GAAT", "CTAT", "GTAT", "CATT", "GATT", "CTTT", "GTTT", 
             #"AAAA", "TAAA", "ATAA", "TTAA", "AATA", "TATA", "ATTA", "TTTA", 
             #"AAAT", "TAAT", "ATAT", "TTAT", "AATT", "TATT", "ATTT", "TTTT")


def KL(a,b):
    """ KL measure 
    """
    #with numpy.errstate(invalid='ignore'):
    d = a * numpy.log(a/b)
    d[numpy.isnan(d)]=0 
    d[numpy.isinf(d)]=0
    return (numpy.sum(d))*10000

def Eucl(a,b):
    """ Euclidean distance between two vectors 
    """
    #with numpy.errstate(invalid='ignore'):
    d = pow(a-b,2)
    d[numpy.isnan(d)]=0
    d[numpy.isinf(d)]=0
    return numpy.sqrt(numpy.sum(d))*10000

def JSD(a,b):
    """ JSD distance
    """
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
        # TBF: should an exit or an Exception be raised here? I guess the program should never go here right?
        #LM yep
        raise ValueError('the strand parameter you set up on the command line can not be interpreted, please use "both","minus" or "plus" , case sensitive, no extra space.')
        exit()
        return 1
   
   
def frequency(seq, ksize=4, strand="both"):
    """ Computes the frequencies of the k-long overlapping oligonucleotides in a sequence. only word composed exclusively of ACGT are taken into account.
    """
    seq = select_strand(seq, strand)
    seq_words = list()
    d = dict()
    #excludes what is not a known characterised nucleotide 
    for s in re.split('[^ACGTacgt]+', seq): 
        seq_letters = list(s)
        ##print(len(s))
        if (len(seq_letters) >= ksize):
            # launch k-1 times the word generation with 1 nucleotide shift
            # every iteration to have overlapping words.
            for i in range(ksize-1):
                #generate words from seq_letters
                words = list((zip(*(iter(seq_letters),)*ksize)))
                # adds the words for this subsequence and 
                # frame to the total list of words
                seq_words.extend(list(map(''.join, words))) 
                # shift one to compute overlapping words at the next iteration
                seq_letters.pop(0) 
    c = Counter(seq_words)
    word_count = sum(c.values())
    if (word_count > 0):
        ret=list()
        word_universe=list(map("".join,(list(product(("C","G","A","T"),("C","G","A","T"),("C","G","A","T"),("C","G","A","T"))))))
        for w in word_universe:
        #for w in word_order:
            if(c.get(w) != None):
                ret.append(c.get(w)/word_count)
            else:
                ret.append(0)
        return ret
    else:
        # TBF should an exit or an Exception be raised here? I guess the program should never go here right?
        # LM Right, it should not go there unless the sequence is strictly shorter than k. there is no reason to exit the program for this if only a few sequences are in this case, although their frequencies is trash. By default, the contig sequence would have to be 3bp long to arrive here. Albeit it can happen, I decided to set the frequency to 1 instead of a vector of values. I haven't tested what happen in that case, TODO. I'm also curious about how long it will take for someone to report this because they tried to compute an oligonucleotide profile of k length with a sequence of k-1 oo'
        return 1


def parallel_subwin_dist(args):
    res=list()
    for window,start,stop,seq_id,mcp_comparison,windows_size,windows_step,ksize,dist,position,contig_size in args:
        #print(window,start,stop,mcp_comparison,windows_size,windows_step,ksize,dist,position)
        #to avoid border effects being a problem in the code, 
        # we use only the simple formula start+windows_size/2 (+/-) windows_step/2 
        # to find the significant center part of the windows. 
        # when several windows overlap, this centered part, 
        # long as the windows step is the most representative of the windows, 
        # not representing as much other part of this window that are 
        # overlapped by other windows. 
        # BUT: this simple formula has border effects, 
        # so we manually correct the start of the first window and
        # the stop of the last window to match the contig borders.
        if(start == (windows_size/2-windows_step/2)):
            displayed_start=1
        else:
             displayed_start=start
        if (stop-windows_step/2+windows_size/2>=contig_size-windows_step and 
                stop-windows_step/2+windows_size/2<=contig_size):
            displayed_stop=contig_size
        else:
            displayed_stop=stop
        if ((window.count('N')/windows_size) <= float(options.n_max_freq_in_windows)):
            if (position==False):
                res.append(globals()[dist](mcp_comparison,
                        numpy.array(frequency(seq=window,ksize=ksize))))
            else:
                res.append([seq_id,displayed_start,displayed_stop,
                        globals()[dist](mcp_comparison, 
                            numpy.array(frequency(seq=window,ksize=ksize)))])
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


def sliding_windows_distances(seq, mcp_comparison, seq_id, dist="KL", windows_size=5000, windows_step=500, ksize=4, position=False):
    #seq=str(Bioseq_record.seq)
    ret=list()
    if(len(seq)<windows_size): #only enough to compute one window, no sliding,
        if((seq.count('N')/windows_size) <= float(options.n_max_freq_in_windows)):
            if(position==False):
                ret.append( globals()[dist](mcp_comparison,numpy.array(frequency(seq=seq,ksize=ksize,strand=strand))))
            else:
                ret.append([seq_id,0,int(len(seq)),globals()[dist](mcp_comparison,frequency(seq=seq,ksize=ksize,strand=strand))])
        else:
            if(position==False):
                ret.append(numpy.nan)
            else:
                ret.append([seq_id,0,int(len(seq)),numpy.nan])
    elif(len(seq)<min_nb_w_per_fasta_for_mul_cpu*windows_step): #not many windows in this contig, so better launching it in serial rather than in parallel
        tmp_out=list()
        if(position==False):
            for s in range(0,len(seq)-windows_size,windows_step):
                window=seq[s:s+windows_size]
                if((window.count('N')/windows_size) <= float(options.n_max_freq_in_windows)):
                    tmp_out.append( globals()[dist](mcp_comparison,numpy.array(frequency(seq=window,ksize=ksize,strand=strand))))
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
                    tmp_out.append([seq_id,displayed_start,displayed_stop,globals()[dist](mcp_comparison,frequency(seq=window,ksize=ksize,strand=strand))])
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

def get_cmd():
    """ read command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--assembly", action="store", required=True, 
                    dest="genome", help="multifasta of the genome assembly")
    parser.add_argument("-c", "--conta", action="store", dest="conta", 
                    help="multifasta of the contaminant species training set")
    parser.add_argument("-r", "--host", action="store", dest="host", 
                    help="optional host species training set in multifasta")
    parser.add_argument("-n", "--n_max_freq_in_windows", action="store", 
                    type=float, dest="n_max_freq_in_windows", default = 0.4,
                    help = "maximum proportion of N tolerated in a window to "
                    " compute the microcomposition anyway [0~1]. "
                    "Too much 'N's will reduce the number of kmer counts and "
                    "will artificially coarse the frequency resolutions")
                    #"If your assembly contains many stretches of 'N's, "
                    #"consider rising this parameter and shortening the windows"
                    #" step in order to allow computation and signal in the "
                    #"output, it might cost you computational time though. "
                    #"Windows which do not meet this criteria will be affected "
                    #"a distance of 'nan'")
    parser.add_argument("-k", "--lgMot", action="store", dest="k", type=int, 
                default=4, help="word wise/ kmer lenght/ k [default:%(default)d]")
    parser.add_argument("-w", "--windows_size", action="store", 
                    dest="windows_size", type=int, 
                    help="Sliding windows size (bp)")
    parser.add_argument("-t", "--windows_step", action="store", 
                    dest="windows_step", type=int, 
                    help="Sliding windows step size(bp)")
    parser.add_argument("-s", "--strand", action="store", default="double", 
                    choices=["double", "leading", "lagging"],
                    help="strand used to compute microcomposition. leading, "
                    "lagging ou double [default:%(default)s]")
    parser.add_argument("-d", "--distance", action="store", dest="dist", 
                  default="JSD", help="distance method between two signatures:"
                    "KL: Kullback-Leibler, "
                    "Eucl : Euclidienne[default:%(default)s], "
                    "JSD : Jensen-Shannon divergence")
    parser.add_argument("-u", "--cpu", action="store", dest="threads_max", 
                    type=int, default=4, 
                    help="how maany threads to use for windows "
                    "microcomposition computation[default:%(default)d]")

    options = parser.parse_args()

    return options

def main():

    # get parameters
    options = get_cmd()

    if (not options.conta):
        whole_seq = Seq("")
        for record in SeqIO.parse(options.genome, "fasta"):
            # TBF: why add a N?
            #LM N is not interpreted to compute frequencies, so by putting one between two sequences, I avoid to account for the chimeric words that would be created by juxtaposing the 2 sequences.
            whole_seq.seq = str(whole_seq)+"N"+str(record.seq)
        genome = numpy.array(frequency(seq=str(whole_seq.seq),ksize=option.k))
        
        if (not options.windows_size and not options.windows_step):
            #parser.error("An input fasta file (-i ) is mandatory")
            print("Warning: A genome was provided but no training set to learn the contamination profile (-c). "
                  "You didn't provide sliding window parameters (-w -t), "
                  "the signature will be computed from the whole genome ", file=sys.stderr)
            
            with open(str(os.path.basename(options.genome))+".microcomposition.mat", 'w') as outf:
                outf.write(str(vector_to_matrix(genome)))
            
        elif (options.windows_size or options.windows_step):
            print("Warning: A genome was provided but no training set to learn the contamination profile (-c). "
                  "You provided sliding window parameters (either -w -t), "
                  "I will compute the genome microcomposition signature and "
                  "I will compute the distance of every window microcomposition to that of the whole genome ", file=sys.stderr)           
            #frq_wins=list()
            
            path = str(os.path.basename(options.genome)) + ".mcp_windows_vs_whole_" + options.dist+".dist"
            with open(path, 'w') as outf:
                for record in SeqIO.parse(options.genome, "fasta"):
                    windows_distances = sliding_windows_distances(str(record.seq), seq_id=record.id, mcp_comparison=genome, 
                                                                  position=True, dist=options.dist, 
                                                                  windows_size=options.windows_size if options.windows_size else 5000, 
                                                                  windows_step=options.windows_step if options.windows_step else 500, 
                                                                  ksize=options.k)
                    for t in windows_distances:
                        outf.write(str("\t".join(map(str,t)))+"\n")
        # TBF: else?
        #LM no else here, if the exact set of parameters given by the user: no conta, no window parameter, nothing can be done... if the sliding window parameters are given, we can scan the genome windows against the whole genome profile, aka HGT detection. Otherwise, we skip to the next set of command line parameters that could make sense.

    else: # TBF: equivalent no? genome is always here (required argument) and conta is present
        #LM yep ^^ you confirm that the parser will enforce the presence of the genome? the lines testing that were removed
    #if(options.genome and options.conta):
        print("A genome and a training set to learn the contamination profile was provided so "
              "I will compute the microcomposition signature of the genome (-i) or host training set if provided (-r), "
              "that of the contamination training set (-c) and those of the sliding windows along the assembly (-i). "
              "with all of that I will output the distances of the microcomposition of every windows compared "
              "to both microcompositions of the contaminant and the host genome.")
        
        fname = str(os.path.basename(options.genome))+".mcp_hostwindows_vs_"
        if(options.host):
            print("using the provided host training set to learn the microcomposition")
            target = options.host
            fname = fname+ "host_"+str(os.path.basename(options.host))+"_"+options.dist+".dist"
        else:
            print("No host training set provided, I will use the whole genome to learn the microcomposition")
            target = options.genome
            fname = fname+"wholegenome_"+options.dist+".dist"
        #print(fname+" totot")
            
        whole_seq=Seq("")
        for record in SeqIO.parse(target, "fasta"):
            whole_seq.seq = str(whole_seq)+"N"+str(record.seq)
            
        genome = numpy.array(frequency(seq=str(whole_seq.seq),ksize=options.k))

        whole_conta = Seq("")
        for record in SeqIO.parse(options.conta, "fasta"):
            whole_conta.seq=str(whole_conta)+"N"+str(record.seq)
            
        conta = numpy.array(frequency(seq=str(whole_conta.seq),ksize=options.k))
  
        with open(fname, 'w') as outf:
            for record in SeqIO.parse(options.genome, "fasta"):
                #f.write('>'+str(record.id)+"\n")
                #frq_wins = sliding_windows_frequencies(str(record.seq), windows_size=options.windows_size if options.windows_size else 5000, windows_step=options.windows_step if options.windows_step else 500, ksize=options.k if options.k else 4)
                windows_distances = sliding_windows_distances(str(record.seq), seq_id=record.id, mcp_comparison=genome, 
                                                              position=True, dist=options.dist, 
                                                              windows_size=options.windows_size if options.windows_size else 5000, 
                                                              windows_step=options.windows_step if options.windows_step else 500, 
                                                              ksize=options.k)
                #windows_distances=list()
                for t in windows_distances:
                    outf.write(str("\t".join(map(str,t)))+"\n")
        
        fname = str(os.path.basename(options.genome))+".mcp_hostwindows_vs_"+"conta_"+str(os.path.basename(options.conta))+"_"+options.dist+".dist"
        #print(fname+" totowdawdt")
        with open(fname, 'w') as outf:
            for record in SeqIO.parse(options.genome, "fasta"):
                #f.write('>'+str(record.id)+"\n")
                #frq_wins = sliding_windows_frequencies(str(record.seq), windows_size=options.windows_size if options.windows_size else 5000, windows_step=options.windows_step if options.windows_step else 500, ksize=options.k if options.k else 4)
                windows_distances = sliding_windows_distances(str(record.seq), seq_id=record.id, mcp_comparison=conta,
                                                              position=True, dist=options.dist, 
                                                              windows_size=options.windows_size if options.windows_size else 5000, 
                                                              windows_step=options.windows_step if options.windows_step else 500, 
                                                              ksize=options.k if options.k else 4)
                #windows_distances=list()
                for t in windows_distances:
                    outf.write(str("\t".join(map(str,t)))+"\n")

    sys.exit(0)

if __name__ == "__main__":
    main()

