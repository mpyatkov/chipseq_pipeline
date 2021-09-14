#Amy Li 08/16/2013
#revised 08/26/2013
#this script is used to filter straight peaks from BED files
#Straight peaks are defined as peaks with identical chromosome number, start, end, and strandness that appears at least five times and are greater than [READ_LENGTH] apart from the closest read with the same strandness
#|==========INSTRUCTIONS FOR RUNNING THIS SCRIPT=================|
#running with paired bash script:
#must have original .bed files in current directory
#./remove_straight_peaks.sh

#running this file directly: 
#must have sorted bed file in current directory with file format: *_sort.bed
#python remove_straight_peaks.py

#Sample output file will appear the same as the input file except with straight peaks removed
#output files: 
#[INPUT_FILE_NAME]_StrgtPks.bed: contains the straight peaks
#[INPUT_FILE_NAME]_StrgtPksRm.bed: contains the non-straight peaks

import os, fnmatch, re, sys

read_length = sys.argv[1]

#printing a list to output stream in tab seperated format
#content: list to be printed
#out: the output file stream
def write_output(out, content):
  for print_line in content:
    for item in print_line:
      out.write(item)
      out.write('\t')
    out.write('\n')
    
#checks if speculated straight peak is within read_length away from neighboring read
#returns true if is within read_length away, false if not within read_length away
def check_distance(chr1, chr2, left, right):
  if (chr1 != chr2):
    return False
  elif int(right) - int(left) > int(read_length):
    return False
  else: 
    return True

#read from input sorted bed file
#identify straight peaks
#write straight peaks into straight peak output file
#write non-straight peaks into straight peaks removed output file


for curr_file in os.listdir('.'):
  if fnmatch.fnmatch(curr_file, '*_sort.bed'):
    f_in = open(curr_file)
    body = f_in.readlines()
    base_name = os.path.splitext(curr_file)[0]
    base_name_trimmed = base_name[:-5]
    output_name = base_name_trimmed+'_StrgtPksRm.bed'
    output_straightpeaks_name = base_name_trimmed + '_StrgtPks.bed'
    f_out = open(output_name, 'w')
    f_straight_peak = open(output_straightpeaks_name, 'w')

    #fields on the current stack
    temp_chrom = ' '
    temp_start = ' '
    temp_end = ' '
    temp_strand = ' '
    temp_index = 0
    stack = []

    index = 0
    num_lines = len(body)
    
    #process each line in input bed file
    while index<num_lines:

      #read current line
      line = body[index]
      array = line.split()
      chrom = array[0]
      start = array[1]
      end = array[2]
      strand = array[5]     
 
      #if current line is the same as lines in stack, add current line to stack  
      if chrom == temp_chrom and start == temp_start and end == temp_end and strand == temp_strand:
        stack += [array]

      else:
      #if not the same, process current stack depending on the size of the stack and reset the stack information as current line
            
        #print stack to non-straight peak file if size of stack is less than five (not straight peak)
        if len(stack) < 5 and len(stack) != 0: 
          write_output(f_out, stack)
        
        else:
        #possible straight-peak, check distance to neighbors
          index_left = temp_index
          left_chr = temp_chrom
          left_start = temp_start

          while (index_left>0 and left_start==temp_start):
            index_left-=1
            left_line = body[index_left].split()
            left_chr = left_line[0]
            left_start = left_line[1]
           
          index_right = temp_index
          right_chr = temp_chrom
          right_start = temp_start

          while (index_right<num_lines-1 and right_start == temp_start):
            index_right+=1
            right_line = body[index_right].split()
            right_chr = right_line[0]
            right_start = right_line[1]

          if len(stack)!=0 and (check_distance(left_chr, temp_chrom, left_start, temp_start) or check_distance(temp_chrom, right_chr, temp_start, right_start)):
            #has neighboring read within read_length away
            #not classified as straight peak
            write_output(f_out, stack)
            
          else:
            #is a straight_peak
            write_output(f_straight_peak, stack)

        #reset stack as current line 
        stack = [array]
        temp_chrom = chrom
        temp_start = start
        temp_end = end
        temp_strand = strand
        temp_index = index

      index+=1    

    #at the end of the file, print stack if not straight peak, check left neighbor only
    if len(stack)<5 and len(stack) != 0:
      write_output(f_out, stack)
    else:     
      #check previous
      index_left = temp_index
      left_chr = temp_chrom
      left_start = temp_start
     
      while (index_left>0 and left_start==temp_start):
        index_left-=1
        left_line = body[index_left].split()
        left_chr = left_line[0]
        left_start = left_line[1]
           
        if len(stack)!=0:
          if check_distance(left_chr, temp_chrom, left_start, temp_start):
            #has neighboring left read within read_length away
            #not classified as straight peak
            write_output(f_out, stack)
            
          else:
            #is a straight_peak
            write_output(f_straight_peak, stack)

    #close input and output files
    f_in.close()
    f_out.close()
    f_straight_peak.close()
