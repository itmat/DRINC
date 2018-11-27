import sys

#System arguments: 1. stranded experiment or not  2. epsilon 3. input sam file 4. Equivalence class (all dups) log file 
#5. Equivalence class (optical dups) log file 6. Stats log file 7. filtering mode (optional 0: just output new sam with tag, without filtering 
#1: filter all dups 2: filter optical dups) 8. new sam file with XD label (optional)

if(len(sys.argv) < 7):
    print("\nusage python DRINK5Cluster.py <stranded info> <radius> <input sam file> <output all dups log file> <output optical dups log file> <output stats file> "
        "[filtering mode] [output sam file w/ tags]"
        "\n\nFor the stranded info, \n - enter 0 for a nonstranded experiment.\n - enter 1 for a stranded experiment.\n"
        "\n\nThe radius is the max possible distance between two optical duplicates.\n - If you set this too high you might call legitimate dups as optical dups.\n - If you set it too low you might miss optical dups.\n"
        "\n\nThe filtering mode is optional.\n - Enter 0 for just outputing new sam with tags, without filtering.\n - Enter 1 for filtering all dups.\n - Enter 2 for filtering optical dups only.\n"
        "\n\nThe output sam file w/ tags field is to be specified if and only if the filtering mode is specified.\n")
    sys.exit()
    
if int(sys.argv[1]) not in {0, 1}:
    print("\nPlease enter either 0 (for nonstranded experiment) or 1 (for stranded experiment) for the first argument\n")
    sys.exit()
    
if len(sys.argv) >= 8 and int(sys.argv[7]) not in {0, 1, 2}:
    print("\nPlease enter either 0 (just output a new sam file with tags, do not filter), 1 (output a sam file filtering all types of dups) or 2 (output a sam file filtering only optical dups)\n")
    sys.exit()
    
with open(sys.argv[3]) as f:     #Opens input sam file to read
#with open('/Users/songpeng/Desktop/Projects/DRINC/Test3.txt') as f: #Opens input sam file to read

    epsilon = int( sys.argv[2] )   #specifies an epsilon. Allow users to input this epsilon themselves
    #epsilon = 1   #specifies an epsilon. Later allow users to input this epsilon themselves
    
    stranded = int( sys.argv[1] )   #specifies whether the sam file is stranded or unstranded. If =0, unstranded. If =1, stranded. 
    #stranded = 1
    
    #####Starts the first phase: constructing a dictionary that contains all readIDs as keys and their unsorted signatures as values
    dic1 = {}    #creates a new empty dictionary, its key is readID, and its value is the signature for that readID
    
    for line in f:                                      #read one line at a time from the input sam file.                                    
        if line.startswith('@'):                        #skip header lines                                            
            continue
        
        Elements = line.split("\t")     #For a non-header line, split its tab-delimited elements into a list  
        readID = Elements[0]    #The 0th entry in this list is the readID
        chromosome = Elements[2]    #the 2nd entry is the chromosome
        start_position = Elements[3]    #3rd entry is the start position
        Cigar = Elements[5] #5th entry is the Cigar string
        ATCGsequence = Elements[9] #9th entry is the ATCG sequence of the read
        
        #construct the signature
        if chromosome == '*':   #if the read did not map, so the chromosome is *,
            sig = (ATCGsequence[:15],)   #then the signature is the first 15 bases in its ATCG sequence
        elif chromosome != '*' and stranded == 0:   #if the read mapped and we have a non-stranded experiment,
            sig = (chromosome, start_position, Cigar)   #the the signature consists of the tuple (chromosome, start_position, Cigar)
        elif chromosome != '*' and stranded == 1:   # if the read mapped and we have a stranded experiment, 
            flag_binary = format(int(Elements[1]), '012b')            #Convert sam flag from decimal to binary 
            strand_info = str(flag_binary[-5])  #The 5th field counting backwards from the right tells us if the read is from the + or - strand. Want it to be a string so that it can be sorted with unmapped reads which have a string of bases as its signature
            sig = (strand_info, chromosome, start_position, Cigar)  ##the the signature consists of the tuple (strand_info, chromosome, start_position, Cigar) 
        else:
            print('Please enter either 0 = unstranded or 1 = stranded for the first argument of the command')
            sys.exit()
            
        if readID in dic1:
            dic1[readID].append(sig)
        else:
            dic1[readID] = [sig]
        #end of constructing the signature
    
#####End of the first phase: constructing a dictionary that contains all readIDs as keys and their unsorted signatures as values stored as a list

#####Now start the second phase: creates another dictionary to store sorted signatures as keys and all readIDs with this signature as values stored as a set
    
dic2 = {}    #creates another dictionary which is going to store sorted signatures as keys and all readIDs with this signature as values stored as a set
    
for read_ID in dic1:    #for every key in dic1,
    signature_list_form = sorted(dic1[read_ID]) #sort its value, which is an unsorted signature. Now the signature is sorted, and its type is a list
    signature = str( signature_list_form )  #turn the signature into a string, since we cannot use lists as keys for dictionaries
    
    if signature in dic2:   #if the signature in dic1 is already in dic2, add the corresponding readID to dic2 as a value
        dic2[signature].add(read_ID)
    else:   #if signature in dic1 not already in dic2, add it to dic2 as a key with its correcponding readID as a value
        dic2[signature] = {read_ID}
        
Tot_number_of_readIDs = len(dic1)   #####this is for the stats part. The number of keys in dic1 is the total number of readIDs in the sam file
        
dic1.clear()    #clear dic1 to free memory
import gc 
gc.collect()    #This frees memory according to python documentation but some people said it depends on the operating system you are using
        
#####End of second phase: creates another dictionary to store sorted signatures as keys and all readIDs with this signature as values stored as a set
#####Dic2 is a partition of the sam file into equivalence classes of duplicates (be it real, PCR or optical duplicates)

#####Start phase 3: generate output file which lists all equivalence classes of duplicates (all kinds of duplicates) which contain more than 1 representative

with open(sys.argv[4], "w") as All_Duplicates:
#with open('/Users/songpeng/Desktop/Projects/DRINC/AllDuplicates.txt', "w") as All_Duplicates:
    
    Dup_sets_decomposing_into_1_subset_of_op_dups = 0   #####This is for the stats part. Set the counter of Dup sets that does not further partition into subsets of optical dups to be 0.
    
    for signtr in dic2: 
        if len( dic2[signtr] ) > 1: #if a particular equivalence class in dic2 contains more than 1 readID, 
            All_Duplicates.write( str( dic2[signtr] ) + '\n' )  #then print it to the log file showing all dups 
        else:   #this is for the stats part. if a equivalence class just contains 1 readID,
            Dup_sets_decomposing_into_1_subset_of_op_dups = Dup_sets_decomposing_into_1_subset_of_op_dups + 1   #then increment the counter of Dup sets that does not further partition into subsets of optical dups by 1
            
dic2.clear()    #clear dic2 to free memory
gc.collect()    #This frees memory according to python documentation but some people said it depends on the operating system you are using
#####End of phase 3: generate output file which lists all equivalence classes of duplicates (all kinds of duplicates) which contain more than 1 representative

#####Start phase 4: generate output file which lists all equivalence classes of optical duplicates only which contain more than 1 representative
with open(sys.argv[4]) as All_Duplicates_read, open(sys.argv[5], "w") as Optical_Duplicates:
#with open('/Users/songpeng/Desktop/Projects/DRINC/AllDuplicates.txt') as All_Duplicates_read, open('/Users/songpeng/Desktop/Projects/DRINC/OpticalDuplicates.txt', "w") as Optical_Duplicates:
    
    dic3 = {}   #####This is for the stats part. Creates an empty dictionary that stores the breakdown of dup sets further into optical dup sets
    Total_number_of_all_dups = 0    ###this is new in DRINC6. For each all dup equivalence class of size n, n-1 is the number of all dups. Then sum over all equiv class.
    
    for equivalence_class in All_Duplicates_read:
        equivalence_class_no_braces = equivalence_class.replace("\n","").replace("{","").replace("}","").replace("'","") #for each equivalence class, remove the braces etc
        representatives = equivalence_class_no_braces.split(", ") #split each equivalence class into its representatives by comma followed by a space
        
        number_of_all_dups_for_this_line = len(representatives) - 1  #this is new in DRINC6. counts the number of all dups for this line
        Total_number_of_all_dups = Total_number_of_all_dups + number_of_all_dups_for_this_line  #this is new in DRINC6. updates counter
        
        list_of_equiv_sets = []     #Now we are within one set of dups(all types). This creates an empty list to store the equivalent sets of optical dups for this signature
        
        
        for i in range( 0, len(representatives)-1 ):   #Now we are within one set of dups(all types). This for loop is for every readID in this set
            base_readID = representatives[i]           #Take the first value(readID) associated with this signature 
            equiv_set = {base_readID}   #start to create a set that contains all readIDs within an epsilon distance to the base readID
            base_readID_split = base_readID.split(":")  #split this base readID by colon
            base_readID_tile = int(base_readID_split[-3])   #this is the tile number of the base readID
            base_readID_x_coord = int(base_readID_split[-2])    #find the x and y coords associated with this base readID
            base_readID_y_coord = int(base_readID_split[-1])
            
            for j in range( i+1, len(representatives) ):  #for subsequent readIDs to the base readID
                next_readID = representatives[j]
                next_readID_split = next_readID.split(":")  #split this subsequent readID by colon
                next_readID_tile = int(next_readID_split[-3])   #this is the tile number of the next readID
                next_readID_x_coord = int(next_readID_split[-2])    #find the x and y coords associated with this subsequent readID
                next_readID_y_coord = int(next_readID_split[-1])
                dif_in_x = abs( base_readID_x_coord - next_readID_x_coord ) #find the difference bwteen the x-coordinates of the baseID and this subsequentID in question
                dif_in_y = abs( base_readID_y_coord - next_readID_y_coord )
                max_dif_in_dist = max(dif_in_x, dif_in_y)
                
                if max_dif_in_dist <= epsilon and base_readID_tile == next_readID_tile:      #If the next_readID is within an epsilon radius to the base readID and they are on the same tile
                    equiv_set.add(next_readID)      #add this next_readID to the equivalent set
                else:
                    continue
            #by this step, we are done making one equivalent set      
            list_of_equiv_sets.append( frozenset(equiv_set) )  #append the equivalent set we just made to the list of equivalent sets
        #by this step, we are done making all the equivalent sets for this signature, except for the last equivalent set which is the singleton set containing the last readID in this signature, because it has no subsequent readIDs to compare radius with
        list_of_equiv_sets.append( frozenset({ representatives[-1] }) )   #append the last singleton equivalent set to the list of equivalent sets
        #by this step, the list of equivalent sets is complete
        #Next, equivalent sets with nonempty intersections should be unioned, to create equivalence classes.
        for k in representatives:  #for each readID in the set of all dups
            equivalence_class = set()   #creates an empty set which is to eventually become a equivalent class of readIDs that are optical dups of each other
            List_to_be_minused = []     #creates an empty list to contain equivalent sets already unioned into an equavalent class, these equivalent sets are to be removed from the list of equivalent sets
            
            for l in list_of_equiv_sets:    #for each equivalent set in the list of equivalent sets
                if k in l:      #if the current readID in the outer-by-1 for loop is in this equivalent set
                    equivalence_class = equivalence_class.union(l)  #union this equivalent set with the current equivalent class
                    List_to_be_minused.append(l)    #also throw this equivalent set into the list of equivalent sets to be discarded from the master list of equivalent sets
                else:
                    continue
                
            interm = list(set(list_of_equiv_sets) - set(List_to_be_minused))    #from the master list of equivalent sets, remove those sets which are already merged into the equivalence class
            interm.append( frozenset(equivalence_class) )   #add the equivalence class to the master list of equivalent sets
            list_of_equiv_sets = interm
        
        #Now, list_of_equiv_sets is finally the list of equivalence classes for this signature!
        
        for EC in list_of_equiv_sets:   #Now for every equivalence class in this master list of equivalent sets
            if len(EC) > 1:     #if the equivalence class contains more than one readID, then write it to output
                Optical_Duplicates.write( str( set(EC) ) + '\n' )
            
        number_of_subsets_of_op_dups = len(list_of_equiv_sets)  #####This is for the stats part. The number of subsets of optical dups within this set of all dups, is the nuber of equivalence classes
            
        if number_of_subsets_of_op_dups in dic3:   #####This is for the stats part. Set the number of subsets of optical dups to be keys in the dictionary. If the key is already present, 
            dic3[number_of_subsets_of_op_dups] = dic3[number_of_subsets_of_op_dups] + 1    #####This is for the stats part. update the number of subsets of optical dups by 1
        else:    #####This is for the stats part. If key is not already in the dictionary,
            dic3[number_of_subsets_of_op_dups] = 1  #####This is for the stats part. update the number of subsets of optical dups to be 1
                
    if 1 in dic3:   #####This is for the stats part. Check if dic3 has any set of all dups that does not decompose. If it has, 
        dic3[1] = dic3[1] + Dup_sets_decomposing_into_1_subset_of_op_dups #####This is for the stats part. then add the number of sets that did not decompose from sets of all dups of length 1 to it
    else:   #####This is for the stats part. If it does not, 
        dic3[1] = Dup_sets_decomposing_into_1_subset_of_op_dups     #####This is for the stats part. creates the key of sets that did not decompose, and add the number of sets that did not decompose from sets of all dups of length 1 to it
#####End of phase 4: generate output file which lists all equivalence classes of optical duplicates only which contain more than 1 representative

#####Start phase 5: generate stats file
with open(sys.argv[6], "w") as Stats, open(sys.argv[5]) as Optical_Duplicates_read, open(sys.argv[4]) as All_Duplicates_read_2:  #Creates a stats log file to write to, also opens the optical duplicate log file to read
#with open('/Users/songpeng/Desktop/Projects/DRINC/DRINCStats.txt', "w") as Stats, open('/Users/songpeng/Desktop/Projects/DRINC/OpticalDuplicates.txt') as Optical_Duplicates_read:   #Creates a stats log file to write to, also opens the optical duplicate log file to read
    
    Stats.write('Total number of readIDs = ' + str(Tot_number_of_readIDs) + '\n' + '\n' )  #writes the total number of readIDs in this sam file
    Stats.write('Number m of sets of all duplicates that further partition into n subsets of optical duplicates, (m,n):' + '\n' )
    
    import collections #use this to sort dic3 by keys
    dic3_sorted_by_key = collections.OrderedDict(sorted(dic3.items()))  #use this to sort dic3 by keys
    
    for partitions in dic3_sorted_by_key:
        Stats.write('(' + str( dic3_sorted_by_key[partitions] ) + ',' + str(partitions) + ')' + '\n')
        
    dic3.clear()    #clear dic3 to free memory
    gc.collect()    #This frees memory according to python documentation but some people said it depends on the operating system you are using
    
    Stats.write('\n')   #writes an emtpy line that separates the next chunk of data
    
    dic4 = {}   #creates a new dictionary that stores the sizes of sets of optical dups. key = optical equivalence class size, value = number of optical equivalence classes with this size
    
    for equivalence_class_2 in Optical_Duplicates_read:
        equivalence_class_2_no_braces = equivalence_class_2.replace("\n","").replace("{","").replace("}","").replace("'","") #for each equivalence class, remove the braces etc
        representatives_2 = equivalence_class_2_no_braces.split(", ") #split each equivalence class into its representatives by comma followed by a space
        size_of_op_dup_Equivalence_Class = len(representatives_2)
        
        if size_of_op_dup_Equivalence_Class in dic4:   #Set the size of optical dup class to be keys in the dictionary. If the key is already present, 
            dic4[size_of_op_dup_Equivalence_Class] = dic4[size_of_op_dup_Equivalence_Class] + 1    #update the number of optical equivalence classes of this size by 1
        else:    #If key is not already in the dictionary,
            dic4[size_of_op_dup_Equivalence_Class] = 1  #updates the number of optical equivalence classes of this size to be 1
            
    dic4_sorted_by_key = collections.OrderedDict(sorted(dic4.items()))  #use this to sort dic4 by keys
    
    Stats.write('Number i of sets of optical duplicates of size j, (i,j):' + '\n' )
    
    To_be_subtracted_from_tot_number_of_readIDs = 0 #set a counter, which is later going to be (number of optical equivalence classes)x(its size) to be subtracted from total number of readIDs to get number of optical equivalence classes of size 1
    
    Total_number_of_opt_dups = 0    ###this is new in DRINC6. For each opt dup equivalence class of size n, n-1 is the number of opt dups. Use n-1 to multiply to the number of opt equi classes of this size to get total number of opt dups for opt equiv classes of this size. Then sum over all equiv class sizes.
    
    for sizes in dic4_sorted_by_key:
        number_of_opt_dups_for_this_size = (sizes - 1) * dic4_sorted_by_key[sizes]  #This new in DRINC6. Find the total number of opt dups from all equiv classes of this size
        Total_number_of_opt_dups = Total_number_of_opt_dups + number_of_opt_dups_for_this_size #This is new in DRINC6. Updates counter
        Stats.write('(' + str( dic4_sorted_by_key[sizes] ) + ',' + str(sizes) + ')' + '\n')
        To_be_subtracted_from_tot_number_of_readIDs = To_be_subtracted_from_tot_number_of_readIDs + sizes * dic4_sorted_by_key[sizes]
    
    Optical_EC_of_size_1 = Tot_number_of_readIDs - To_be_subtracted_from_tot_number_of_readIDs
    Stats.write('(' + str( Optical_EC_of_size_1 ) + ',' + str(1) + ')' + '\n' + '\n')
    
    dic4.clear()    #clear dic4 to free memory
    gc.collect()    #This frees memory according to python documentation but some people said it depends on the operating system you are using

    #####This part is new to DRINC6. Finds the % all dups for this sam file
    percent_all_dup = Total_number_of_all_dups * 100 / Tot_number_of_readIDs #this new in DRINC6. Write total % of all dups in this sam file.
    Stats.write('%_all_dups = ' + str(percent_all_dup) + '%' + '\n' + '\n')   #this new in DRINC6. Write total % of opt dups in this sam file.
        
    percent_opt_dup = Total_number_of_opt_dups * 100 / Tot_number_of_readIDs #this new in DRINC6. Write total % of opt dups in this sam file.
    Stats.write('%_optical_dups = ' + str(percent_opt_dup) + '%' + '\n' + '\n')   #this new in DRINC6. Write total % of opt dups in this sam file.
    #####End of new part to DRINC6
#####End of phase 5: generate stats file
#####Start of phase 6: filter sam file
if(len(sys.argv) < 8):
    sys.exit()
with open(sys.argv[4]) as All_Duplicates_read2, open(sys.argv[5]) as Optical_Duplicates_read2, open(sys.argv[3]) as sam_in, open(sys.argv[8], 'w') as sam_out:
#with open('/Users/songpeng/Desktop/Projects/DRINC/AllDuplicates.txt') as All_Duplicates_read2, open('/Users/songpeng/Desktop/Projects/DRINC/OpticalDuplicates.txt') as Optical_Duplicates_read2, open('/Users/songpeng/Desktop/Projects/DRINC/Test3.txt') as sam_in, open('/Users/songpeng/Desktop/Projects/DRINC/sam_out.txt', 'w') as sam_out:
    
    filtering_mode = int( sys.argv[7] )  #set filtering mode. 0: just output new sam with tag, without filtering #1: filter all types of dups, only letting one representative remain 2: filter optical dups, only letting one representative remain
    #filtering_mode = 2  #set filtering mode. 0: just output new sam with tag, without filtering #1: filter all types of dups, only letting one representative remain 2: filter optical dups, only letting one representative remain
    
    dic5 = {}   #creates an empty dictionary to store the XD tags
    
    lines_all_dup = All_Duplicates_read2.readlines()    #break the log file with all dups into lines
    
    for n in range( len(lines_all_dup) ):   #for each line in the log file that contains all dups, make all readIDs in the line into a list
        output_line = lines_all_dup[n].replace("\n","").replace("{","").replace("}","").replace("'","") 
        output_line_split = output_line.split(", ") 
        
        for m in range( len(output_line_split) ):   #for each representative in this equiv class
                dic5[ output_line_split[m] ] = (n+1, m+1)    #record its row(line) number and column number into dic5. Plus 1 because python counts from 0 but on sam file we want to count from 1.
                
    lines_opt_dup = Optical_Duplicates_read2.readlines()    #break the log file with opt dups into lines
    
    for x in range( len(lines_opt_dup) ):   #for each line in the log file that contains all dups, make all readIDs in the line into a list
        output_line2 = lines_opt_dup[x].replace("\n","").replace("{","").replace("}","").replace("'","") 
        output_line_split2 = output_line2.split(", ") 
        
        for y in range( len(output_line_split2) ):   #for each representative in this equiv class
                dic5[ output_line_split2[y] ] = dic5[ output_line_split2[y] ] + (x+1, y+1)    #append its row(line) number and column number into dic5. Plus 1 because python counts from 0 but on sam file we want to count from 1.
    
    for sam_line in sam_in: #for every line in the original sam file
        
        if sam_line.startswith('@'):    # for header lines, just write them to the new sam file
            sam_out.write( sam_line )
            continue
        
        sam_Elements = sam_line.split("\t") #again split the line by \t
        sam_readID = sam_Elements[0]    #the first field is the readID
        
        if sam_readID in dic5 and filtering_mode == 0: #if a readID is in dic5 i.e. tagged but the filtering mode is 0, i.e. do not filter, then just output it with tag
            sam_out.write( sam_line.rstrip() )
            if len( dic5[sam_readID] ) < 4: #if this readID is tagged as a dup but not an opt dup, output its row number in the all dup file
                sam_out.write( '\t' + 'XD:s:' + str( dic5[sam_readID][0] ) + '\n' )
            else:    #if this readID is tagged as a opt dup, output its row number in both the all dup log file and the opt log file
                sam_out.write( '\t' + 'XD:s:' + str( dic5[sam_readID][0] ) + '.' + str( dic5[sam_readID][2] ) + '\n' )
        elif sam_readID in dic5 and filtering_mode == 1:    #if a readID is in dic5 i.e. tagged and the filtering mode is 1, i.e. filter all dups, only output it when it is the first column in each row of the all dup log file.
            if dic5[sam_readID][1] == 1:
                sam_out.write( sam_line.rstrip() )
                sam_out.write( '\t' + 'XD:s:' + str( dic5[sam_readID][0] ) + '\n' )
            else: 
                continue
        elif sam_readID in dic5 and filtering_mode == 2:    #if a readID is in dic5 i.e. tagged and the filtering mode is 2, i.e. filter opt dups only,
            if len( dic5[sam_readID] ) < 4: #then if this readID is not in the opt dup file, just output it with all dup tag
                sam_out.write( sam_line.rstrip() )
                sam_out.write( '\t' + 'XD:s:' + str( dic5[sam_readID][0] ) + '\n' )
            elif len( dic5[sam_readID] ) == 4 and dic5[sam_readID][3] == 1: #then if this readID is in the opt dup file, only output it if it is in the first column
                sam_out.write( sam_line.rstrip() )
                sam_out.write( '\t' + 'XD:s:' + str( dic5[sam_readID][0] ) + '.' + str( dic5[sam_readID][2] ) + '\n' )
            else:
                continue
        else:   #if a readID is not in dic5 i.e. not tagged, write it as it is
            sam_out.write( sam_line )

#dic5.clear()    #clear dic5 to free memory
#gc.collect()    #This frees memory according to python documentation but some people said it depends on the operating system you are using
#####End of phase 6: filter sam file