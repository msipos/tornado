# -*- coding: iso-8859-1 -*-
#VOLTRON
#A script to merge alignments from NAST and RDPv10/Infernal
# Written by Patricio Jeraldo
# Edited by Maksim Sipos

#Import stuff
import sys, re, string
from Bio import SeqIO, Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bitarray import bitarray

#define a useful string
inf_consensus_seq= '#=GC_RF'
global_log_file = ''

class MergeError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

def pad_with_gaps(d):
  #get the maximum column number
  maximum= -2
  for l in d.itervalues():
    maximum= (maximum, l[2])[l[2]> maximum]
  #now add the right-hand-side gaps
  for k in d.iterkeys():
    if len(d[k][0]) > 0:
      d[k][0] += ''.center(maximum-d[k][2], '-')
  #now look for the minimum
  minimum= 100000000
  for l in d.itervalues():
    minimum= (minimum, l[1])[(l[1]< minimum) and (l[1] != -1)]
  #add the LHS gaps, and deal with the empties
  for k in d.iterkeys():
    if len(d[k][0]) > 0:
      d[k][0] = ''.center(d[k][1]-minimum, '-') + d[k][0]
    else:
      d[k][0] = ''.center(maximum-minimum, '-')
  #and I think we're done, fone now
  return d



#put the consensus sequence in the record list
def add_consensus(recs):
  #get the first Seq from the list
  length= len(recs.pop().seq)
  #create new SeqRecord with gaps and return it
  return SeqRecord(seq=Seq(''.center(length,'-'), Alphabet.generic_nucleotide), id='#=GC_RF', name='', description='')





#create an alignment object from a ser_record list
from Bio.Align.Generic import Alignment
def aln_from_record(record_list):
  aln= Alignment(Alphabet.generic_nucleotide)
  #the following is a hack, documented in
  #bugs 2553/2554 of biopython
  #A workaround should be possible
  aln_records= aln._records
  for rec in record_list:
    aln_records.append(rec)
  #return the alignment
  return aln

def remove_gapped_cols(seq_records):
  #create alignment
  length= len(seq_records[0].seq.tostring())
  #using bitarrays.... generate translation tables
  trans_table= string.maketrans(string.uppercase + string.lowercase + '.-', ''.join(['0' for i in xrange(52)]) + '11')
  gap_bitarray= length*bitarray('1')

  #translate gaps to '1', characters to '0'.. bitwise AND wit
  #previous colums
  for record in seq_records:
      gap_bitarray &= bitarray(record.seq.tostring().translate(trans_table))

  #list of non-gap columns, the ones we want to keep
  nongaps= gap_bitarray.search('0')
#  new_rec_list= []
  #now that we have the gap list, remove the guilty columns
  for rec in seq_records:
    seq_s= rec.seq.tostring()
    new_s= ''.join([seq_s[i] for i in nongaps])
    yield SeqRecord(seq=Seq(new_s, Alphabet.generic_nucleotide), id=rec.id, name='', description='')
    #think of this as a generator... let's see what happens
#    yield seqr
#    new_rec_list.append(seqr)

#  return new_rec_list

def find_and_replace(nast_d, inf_d, common, columns):
  '''Replace the unalignable columns in the infernal alignment
   with the corresponding alignment coming from NAST.

  '''
  #compile the anti-gap regexp
  antigap_re= re.compile('[.-]*')
  #do this for every suspect column
  # a list of the dictionaries for each column
  nast_master= []
  for col_start,col_end in columns:
    #create the dict for NAST columns
    nast_seq= dict()
    #look in the infernal dictionary
    for key in common:
      sequence= inf_d[key].seq.tostring()[col_start:col_end]
      #lookahead the next 10
      lookahead= inf_d[key].seq.tostring()[col_end:]
      lookahead=antigap_re.sub('', lookahead.upper())
      lookahead=lookahead[:5]
      #format the sequence
      #get rid of gaps
      #also, uppercase the sequence
      sequence= antigap_re.sub('', sequence.upper())
      #requirement= antigap_re.sub('', requirement.upper())
      #if seq string is empty, then add a special one in the cols
      if len(sequence)==0 :
	nast_seq[key]= ['',-1,-1]
	continue
      #now, create the gap string
      gap_str='(' + '[.-]*'.join(sequence) + ')'
      if len(lookahead) >0:
	gap_str += '(?=[.-]*' + '[.-]*'.join(lookahead) + '[.-]*)'
      #find the corresponding sequence in NAST
      mo= re.search(gap_str, nast_d[key].seq.tostring())
      if mo is None:
	print 'something happened, nothing matched'
	print gap_str
	print key
	sys.exit(1)
      #good, now store it
      nast_seq[key]= [mo.group(0), mo.start(), mo.end()]
    #done
    #now insert the padding in the gaps
    nast_seq= pad_with_gaps(nast_seq)
    #the padded sequence is now stored for later use
    nast_master.append(nast_seq)
  ###
  #now that the sequences were identified and processed,
  #let's create the new Seq and SeqRecord objects
  #record list... this will be eventually returned by
  #this function
  record_list= []
  #now, for every common sequence present in infernal
  for key in common:
    #if key not in common: continue
    #get the raw sequence string
    orig_seq= inf_d[key].seq.tostring()
    new_seq=''
    #splice the sequence into the different sections
    spliced= []
    aux= 0
    for col_start,col_end in columns:
      spliced.append(orig_seq[aux:col_start])
      aux= col_end
    spliced.append(orig_seq[aux:])
    #now add the different pieces
    for i in range(0,len(nast_master)):
      new_seq+= spliced[i] + nast_master[i][key][0]
    new_seq+= spliced[-1]
    #finally, create a SeqRecord object
    rec= SeqRecord(seq=Seq(new_seq, Alphabet.generic_nucleotide), id=key, name='', description='')
    #and add the SeqRecord to the list
    record_list.append(rec)
  #are we done?
  return record_list


def find_infernal_cols(d, t):
  '''Find the "unalignable" columns in the infernal alignment, for future
  replacement with the NAST alignment.

  This function takes the infernal alignment dictionary, and looks for a
  continuous subset of lowercase characters and/or periods (a signature
  that infernal leaves to indicate a section which it couldn't not align to
  the structure. The function returns a list of tuples containing the
  initial and final column numbers delimiting the unaligned region.

  This is, by far, the most important function of the merging script.
  '''
  #look for consecutive lowercase letters or dots
  #but only a minimum of 't' are required to get a result
  temp= '[a-z.]{%d,}' % t
  rexp= re.compile(temp)
  #get a sequence string. Use the first one in the dict,
  #as any sequence should work (FIXME: baaaaaaad assumption)
  #get the list of sequence ids
  seq_ids= d.keys()
  #see how big it is, pick the one in the middle :)
  seq_at_middle= len(seq_ids)/2
  #get the sequence at the middle
  seq_str= d[seq_ids[seq_at_middle]].seq.tostring()
  #match the RegExp in the seq string
  mo= rexp.finditer(seq_str)
  #iterate over the matched strings
  #add the start/end columns to the list of columns
  #and return them
  return [(m.start(), m.end()) for m in mo]

def trim_alns(nast_d, inf_d, common):
    global global_log_file

    '''Trim alignments to the shortest common set for each sequence. '''
    antigap_re= re.compile('[.-]*')
    #look through the commons
    for key in common:
	try:
	    #get a NAST sequence, the first 10 non-gaps
	    nast= nast_d[key].seq.tostring()
	    #remove the gaps charachters
	    nast= antigap_re.sub('', nast)
	    inf= inf_d[key].seq.tostring()
	    inf= antigap_re.sub('', inf)
	    if len(nast) < 16:
		raise MergeError('sequence %s too short: %s!' % (key, nast))
	    #get the first and last 15 letters of the sequence
	    m_nast= re.search('^[A-Z]{15}',nast).group(0)
	    m_end=  re.search('[A-Z]{15}$',nast).group(0)
	    #create the regexps
	    gap_start= '[.-]*'.join(m_nast)
	    gap_end= '[.-]*'.join(m_end)
	    #match it in the corresponding infernal
	    mo_s= re.search(gap_start, inf_d[key].seq.tostring(), re.IGNORECASE)
	    mo_e= re.search(gap_end, inf_d[key].seq.tostring(), re.IGNORECASE)
	    #make the sequence mutable
	    inf_d[key].seq= inf_d[key].seq.tomutable()
	    nast_d[key].seq= nast_d[key].seq.tomutable()
	    #check if there's actually a match in the infernal sequence
	    #and delete the leftovers
	    if mo_s is not None:
		for i in range(0,mo_s.start()):
		    inf_d[key].seq[i]='-'
	    #try doing the substitution the other way around
	    else:
		#get the 15 first letters of the infernal sequence
		m_nast= re.search('^[A-Z]{15}',inf,re.IGNORECASE).group(0)
		#create the regexp
		gap_start= '[.-]*'.join(m_nast)
		#search for it in the NAST alignment
		mo_s= re.search(gap_start, nast_d[key].seq.tostring(), re.IGNORECASE)
		#and delete the leftovers
		for i in range(0,mo_s.start()):
		    nast_d[key].seq[i]='-'
	    #now, the same thing but with the end of the sequence
	    if mo_e is not None:
		for i in range(mo_e.end(), len(inf_d[key].seq)):
		    inf_d[key].seq[i]='-'
	    #reverse trim
	    else:
		#get the 15 last letters
		m_end=  re.search('[A-Z]{15}$',inf, re.IGNORECASE).group(0)
		#create the regexp
		gap_end= '[.-]*'.join(m_end)
		#search for the pattern
		mo_e= re.search(gap_end, nast_d[key].seq.tostring(), re.IGNORECASE)
		#and delete the leftovers
		for i in range(mo_e.end(), len(nast_d[key].seq)):
		    nast_d[key].seq[i]='-'
	    #we're done now with this sequence
	    inf_d[key].seq= inf_d[key].seq.toseq()
	    nast_d[key].seq= nast_d[key].seq.toseq()
	except:
	    global_log_file.write("Problem with trimming alignments of '%s'\n" % key)

  #and we return the modified dictionaries
  #shouldn't be necessary
  #return nast_d, inf_d

def get_common(nast_d, inf_d):
  '''Return a set of the sequence IDs common to both aliignment objects.

  '''
  nast_key_set= set(nast_d.keys())
  inf_key_set= set(inf_d.keys())
  #the common set is the intersection
  common_set= nast_key_set & inf_key_set
  return common_set

def fasta_dict(nast_file, inf_file):
  '''Make dictionaries of sequence objects.

  I'll assume they are untouched by jalview
  but sanity code should be here, before loading the files as dicts.
  '''
  nast_d= SeqIO.to_dict(SeqIO.parse(open(nast_file, 'rU'), 'fasta'))
  inf_d= SeqIO.to_dict(SeqIO.parse(open(inf_file, 'rU'), 'fasta'))

  return nast_d, inf_d

def merge(nast_file, inf_file, result_file, log_file, threshold):
  global global_log_file
  
  global_log_file = open(log_file, 'w')

  #load the files, obtain dictionaries
  nast_d, inf_d= fasta_dict(nast_file, inf_file)
  #get common sequence IDs
  common_set= get_common(nast_d, inf_d)
  #trim the infernal seqs with the common set
  trim_alns(nast_d, inf_d, common_set)
  #get the list of special columns from the infernal alignment
  cols= find_infernal_cols(inf_d, threshold)
  #do the find and replace
  new_records= find_and_replace(nast_d, inf_d, common_set, cols)
  #save some memory
  del nast_d, inf_d, common_set, cols
  #remove gapped columns
  new_records= remove_gapped_cols(new_records)

  global_log_file.close()

  #save the records as a fasta file
  with open(result_file, "w") as ofile:
      SeqIO.write(new_records, ofile, 'fasta')

#For testing purposes
if __name__ == '__main__':
    nast_file = sys.argv[1]
    inf_file = sys.argv[2]
    out_file = sys.argv[3]
    log_file = sys.argv[4]
    merge(nast_file, inf_file, out_file, log_file, 10)
