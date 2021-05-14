#!/usr/bin/python3

# imports
import re
import math
import requests
import sys
from subprocess import call
import os

######################## 1) LOAD FILE #####################################

def load_pdb(pdb_id):
    """ Function outputs amino acid sequence from the PDB ID
    that is requested by the user.
    A comment will be made if the program removes any unknown amino acids
    """
    # if removes any spaces written by the user
    pdb_id = pdb_id.replace(' ', '')
    # test if the length of the input is correct.
    if len(pdb_id) != 4:
        print('The input ', pdb_id, ' is not a PDB ID', '\n')
        sys.exit(1)

    # imports pdb file from the PDB website
    try:
        pdb_web = requests.get('http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=PDB&compression=NO&structureId=' + pdb_id)
        outfile = open('pdb_file.pdb', 'w')
    except requests.ConnectionError as error:
        sys.stdout.write('Error when loading file from the PDB server: '
                            + str(error) + '\n')
        sys.exit(1)
    except IOError as error:
        sys.stdout.write('Cannot write file, reason:' + str(error) + '\n')
        sys.exit(1)
    except:
        sys.stdout.write('Unknown error has occurred when downloading the PDB file from the web site')

    # changes the format to txt file and downloads it on the computer
    for char in pdb_web.text:
        if char == '\n':
            outfile.write('\n')
        else:
            outfile.write(char)
    outfile.close()

    # attempt to open the pdb_file
    try:
        infile = open('pdb_file.pdb', 'r')
        # if the PDB file downloaded is empty and error will be raised
        if not os.stat('pdb_file.pdb').st_size:
            sys.stdout.write('The pdb_file downloaded is empty. This may be because no PDB file exist for this PDB ID\n')
            sys.exit(1)
    except IOError as error:
        sys.stdout.write('Cannot open file, reason:' + error.strerror + '\n')
        sys.exit(1)

    # used to insure that the sequence only contains the 20 amino acids
    aa_verify = (
    'ALA', 'ARG', 'ASN','ASP', 'CYS', 'GLU', 'GLN',
    'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
    'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
                )
    re_prot_term = None
    re_aa_pdb = None
    aa_num = 0
    aa_seq = []
    stop_flag = False
    # find all Asn og Gln in sequnece
    for line in infile:
        # find the amino acids in pdb file
        re_aa_pdb = re.search('^ATOM\s+\d+\s+\w.+?\s+(\w+\s\w\s+\d+)\s+.+$', line)
        re_prot_term = re.search('^TER\s+\d*\s*(ASN)\s*\w\s*\d*', line)

        # if there is a terminal Asn on one of the chains in the
        # PDB file, then that is added to the correct position
        if re_prot_term is not None:
            aa_seq.append('ASN')

        if 'ENDMDL   ' in line:
            stop_flag = True

        # splits the amino acids into [amino_acid, chain, amino acid number]
        if re_aa_pdb is not None and not stop_flag:
            aa_pdb = re_aa_pdb.group(1)
            aa_pos = aa_pdb.split()

            # if amino acid number is different from earlier shown and
            # the earlier number is greater than the new number.
            if aa_pos[2] != aa_num:
                aa_num = aa_pos[2]
                # checks that all amino acids are the 20 canonical amino acids
                if aa_pos[0] in aa_verify:
                    aa_seq.append(aa_pos[0])
                else:
                    print('amino acid ', aa_pos[0], 'at position', aa_num,
                            'is non canonical amino acid')
                    aa_seq.append(aa_pos[0])
    infile.close()
    return aa_seq

#####################  2) IMPORT DATA FROM DSSP  ############################

def dssp_import():
    """ Function runes the DSSP program and saves the result in
    a file named dssp.txt ready to be used by function named HH_bind()
    """
    # sends the PDB file to the DSSP server and retrieves DSSP file
    call('dssp-2.0.4 pdb_file.pdb > dssp.txt', shell = True)

    # opens DSSP file if it can be found
    try:
        dssp_input = open('dssp.txt', 'r')
    except IOError as error:
        sys.stdout.write("Cannot open file, reason:" + error.strerror + "\n")
        sys.exit(1)
    return dssp_input

#######################  3) PRIMARY STRUCTURE  ##############################

def prim_calc(next_aa, steric, CONS):
    """ Function calculates the half times of a asparagine or glutamine
    from the primary structure
    Use in primary_half_time()
    """
    # make calculation of each amines half time
    half_time = math.log(2) / 86400 * math.exp(steric / 100 + CONS)
    # correcting for hydrolysis
    half_time_hydro = 1 / ((1 / 8000) + (1 / half_time ))
    return half_time_hydro

###########################################################################

def primary_half_time(aa_seq):
    """ Function finds steric hindrance values for the amino acids that fit the
    criteria for deamination. A function is called to calculates the
    deamination half time from the primary structure.
    The asparagine and glutamine data is collected by
    this function and append to lists.
    """
    ## dictionaries for all defined steric hindrance values for the Asn-Xxx
    ##  and Gln-Xxx.
    # P,Y and W are not included for Glutamine as there is no experimental
    # data available
    steric_asn_info = {
    		        'GLY':0.0, 'HIS':219.5, 'SER':262.0, 'ALA':306.8,
                    'ASP':333.7, 'THR':370.6, 'CYS':378.8, 'LYS':393.1,
                    'MET':396.5, 'GLU':401.1, 'ARG':400.7, 'PHE':411.9,
                    'TYR':425.1, 'TRP':444.0, 'LEU':466.3, 'VAL':538.5,
                    'ILE':561.1,
                    'PRO':500, 'ASN':40, 'GLN':60,
                    }
    steric_gln_info = {
                    'GLY':0.0, 'HIS':350.4, 'SER':334.4, 'ALA':347.1,
                    'ASP':562.3, 'THR':305.3, 'CYS':127.6, 'LYS':353.7,
                    'MET':273.5, 'GLU':482.0, 'ARG':459.5, 'PHE':602.0,
                    'LEU':367.7, 'VAL':399.7,'ILE':379.0
    		        }
    sum_asn_data, sum_gln_data, asn_pos_lst = [], [], []
    count = 0
    AA_LENGTH = len(aa_seq)
    # finds all Asn and Gln in sequence and appends the amino acids on the
    # C-terminal of the Asn, Gln
    for num in range(AA_LENGTH):
        count += 1
        # every amino acids in the sequnece is defined as aa
        aa = aa_seq[num]
        # the previous amino acid is stored as previous_aa
        previous_aa = aa_seq[num-1]
        # finds asparagine residues
        if aa == 'ASN' and num < AA_LENGTH - 1:
            # finds amino acid on C-term of Asn and test conditions
            next_aa = aa_seq[num+1]
            # inputs the positions of the asparagine residues
            asn_pos_lst.append(count)
            # find steric hindrance value from dict of
            steric_asn = steric_asn_info[next_aa]
            # specific constant for asparagine deamination
            ASN_CONS = 11.863
            # calls function to make primary calculations
            half_time_hydro = prim_calc(next_aa, steric_asn, ASN_CONS)
            # appends half time to list
            sum_asn_data.append(half_time_hydro)

        # if the terminal amino acid is an asparagine
        # the value will be calculated
        elif num == AA_LENGTH - 1 and aa == 'ASN':
            # find steric hindrance value from dict of
            steric_asn = steric_asn_info[aa]
            # specific constant for asparagine deamination
            ASN_CONS = 11.863
            # calls function to make primary calculations
            half_time_hydro = prim_calc(aa, steric_asn, ASN_CONS)
            # appends half time to list
            sum_asn_data.append(half_time_hydro)

        # finds Glutamine residues
        elif (aa == 'GLN' and previous_aa != 'ASN' and previous_aa != 'GLN' and
                num < AA_LENGTH - 1):
            # finds amino acid on C-term of Gln and test conditions
            next_aa = aa_seq[num+1]
            # different conditions are present due to insufficient data
            exclusion_aa = ['PRO', 'TRP', 'TYR', 'GLN', 'ASN']
            if next_aa not in exclusion_aa:
                # find steric hindrance value from dict
                steric_gln = steric_gln_info[next_aa]
                # specific constant for glutamine deamination
                GLN_CONS = 18.311
                # calls function to make primary calculations
                half_time_hydro = prim_calc(next_aa, steric_gln, GLN_CONS)
                # appends half time to list
                sum_gln_data.append(half_time_hydro)
    return asn_pos_lst, sum_asn_data, sum_gln_data

##########################  3) S-VALUE FIND   ##############################

def structure_finder(structure):
    """ Function simplifies the DSSP output for the secondary structure.
    The function output is a string consisting of alpha helix = H,
    beta-sheet = B and coil = T.
    Function used in HH_bind()
    """
    if structure == 'H' or structure == 'I' or structure == 'G':
        structure_point = 'H'
    # finds beta-sheet structure
    elif structure == 'B' or structure == 'E':
        structure_point = 'B'
    # finds coil structure
    elif structure == 'T' or structure == 'S' or structure == ' ':
        structure_point = 'T'
    return structure_point

##################################################################################

def S8_finder(line):
    """ Function finds the S8 value.
    S8 value is the number of interactions with the Asn side chain.
    Function used by HH_bind()
    """
    NH_bind1 = line[38:45].split( maxsplit = 0)
    NH_bind2 = line[63:67].split( maxsplit = 0)

    # find S8. number of H-H to the backbone of the Neighboring  N
    if NH_bind1[0] != '0' or NH_bind2[0] != '0':
        S8 = 1
    else:
        S8 = 0
    return S8

##################################################################################

def S7_finder(line):
    """ Function finds the S7 value from the DSSP_input.
    The S7 is the number of interactions to the NH2 on the COOH side of Asn
    Function is called in the HH_bond() function.
    """
    # find S7. number of H-H bonds on Asn NH2 side chain
    N_chain_H_bind = line[25:33].split()
    # if 0 in list there is 1 or 0 H-H
    if '0' in N_chain_H_bind:
        # if S7 = 0 then there is no H-H to amine
        if N_chain_H_bind[0] == '0' and N_chain_H_bind[1] == '0':
            S7 = 0
        else:
            S7 = 1
    # if there is no 0 in list there is two H-H on NH2
    else:
        S7 = 2
    return S7

##################################################################################

def HH_bond(dssp_input):
    """ Function reads the DSSP input file and extracts information
    to calculate S7 and S8. The function also outputs a simplification of
    the structure_data
    """
    S7_collect, S8_collect = [], []
    previous_aa = None
    aa = None
    flag_start = False
    structure_data = ''
    amin_search = False
    for line in dssp_input:
        # finds header line
        if '  #  RESIDUE AA STRUCTURE ' in line:
            flag_start = True
        # finds all lines with structural information
        elif flag_start:
            # find all the secondary structures from DSSP
            CHAIN = 11
            # amino acid defined as aa
            aa = line[CHAIN+2]
            # structure information position
            structure = line[CHAIN+5]
            # function that simplifies the DSSP output. The function output is
            # appended to a string
            structure_point = structure_finder(structure)
            structure_data += (structure_point)

            # finds all Asn residues in the sequence
            if aa == 'N':
                # calls function to find the S7 value
                S7 = S7_finder(line)
                S7_collect.append(S7)
                amin_search = True
                # special case for S8 where to asn follow each other.
                # finds all the residues that are on the COOH side
                #  of a asparagine
                if previous_aa == 'N':
                        # Function finds first and second N-H binding
                        # positon for backbone
                        S8 = S8_finder(line)
                        S8_collect.append(S8)

            # all the C-terminal proteins that fit the criteria are
            # found downstream of Asn
            elif amin_search:
                # reset flag to false
                amin_search = False
                # Function finds first and second N-H binding
                # positon for backbone
                S8 = S8_finder(line)
                S8_collect.append(S8)

            # saves the previous amino acid
            previous_aa = aa

    if len(structure_data) != len(aa_seq):
        print('WORNING: The length of the DSSP import file is not',
         'the same as the length of the amino acid list found in the PDB file\n',
         'This may lead to incorrect results', 'Control the PDB and DSSP files\n',
         'length data from PDB file = ', len(aa_seq),
         ' length data from DSSP = ', len(structure_data), sep='')
    return S7_collect, S8_collect, structure_data

#############################################################################

def helix_S(asn_pos_lst, structure_data):
    """ Function finds the asparagine positioned in helix structures.
    All asparagine that are in a helix are given a one zero value
    on ether S1,S2 or S3.
    """
    S1_collect, S2_collect, S3_collect = [], [], []
    # loops through the positions of the asparagine hits, to find the S1-3 values
    for element in asn_pos_lst:
        # finds the asparagines that are in a helix
        if structure_data[element] == 'H':

            # find the structural information of the amino acids adjacent to
            # the asparagine
            sub_string = structure_data[element-2:element+3]
            # if there are only helixes in the sub string then the asparagine is
            # placed deep inside of a helix. therefor S3 = 1 , S1 = 0, S2 = 0
            if 'T' in sub_string or 'B' in sub_string:
                # count the 'H' in the string before and after the hit.
                # if there can only be 2 scenarios for each of the two sub_strings
                before_element = sub_string[:2].count('H')
                after_element = sub_string[-2:].count('H')
                # special case where S1 = 1
                if sub_string[:2] == 'HT' or sub_string[:2] == 'HB':
                    S1 = 1
                    S2 = 0
                # special case where S2 = 1
                elif sub_string[-2:] == 'TH' or sub_string[-2:] == 'BH':
                    S2 = 1
                    S1 = 0
                # if a helix is flanked by other structures
                elif before_element == after_element:
                    S1 = before_element
                    S2 = 0
                # if there are more helix structures upstream than downstream
                elif before_element <= after_element:
                    S1 = before_element + 1
                    S2 = 0
                # if there are more helix structures downstream than upstream
                elif after_element < before_element:
                    S2 = after_element + 1
                    S1 = 0
                S1_collect.append(S1)
                S2_collect.append(S2)
                S3_collect.append(0)
            # then S3 = 1, S1=0, S2=0
            else:
                S1_collect.append(0)
                S2_collect.append(0)
                S3_collect.append(1)
        # if the asparagine is not placed in a helix then S1,S2,S3 = 0
        else:
            S1_collect.append(0)
            S2_collect.append(0)
            S3_collect.append(0)
    return S1_collect, S2_collect, S3_collect

###############################################################################

def S11S12_find(strukture_string):
    """ Function finds a value for ether S11 or S12.
    A op to 5 long string is inputted for upstream or downstream of a Asn.
    Function is called in function end_chain.
    """
    S = 0
    for char in strukture_string:
        if char == 'T':
            S += 1
        else:
            break
    return S

###############################################################################

def end_chain(asn_pos_lst, structure_data):
    """ Function calculates S10-S12 from the structural data
    and from a list of Asn positions and a string of the structure_data
    """
    S10_collect, S11_collect, S12_collect = [], [], []
    STRUCTURE_LENGTH = len(structure_data)
    # loops through all the Asn residues found
    for asn_pos in asn_pos_lst:
        # looks at the last 20 amino acids in the peptide chain
        if asn_pos >= STRUCTURE_LENGTH - 20 and structure_data[asn_pos] == 'T':
            # the structure data from the residues +3 and -3 from the Asn
            downstream3 = structure_data[asn_pos:asn_pos+3]
            upstream3 = structure_data[asn_pos:asn_pos+3]

            # if the structures is 3 from the end of the chain and
            # contains only coil structures (T) then S10 = 1
            if ('TTT' in downstream3 or 'TTT' in upstream3):
                S10 = 1
                S10_collect.append(S10)
            # if the chain is further than 3 from helix or beta sheet
            # then S10 = 0 structures
            else:
                S10_collect.append(0)
            # the structure data from +5 and -5 of the Asn
            downstream = structure_data[asn_pos:asn_pos+5]
            upstream_rev = structure_data[asn_pos-1:asn_pos-6:-1]
            # find S11
            # looks at the secondary structures of the Asn and exceeding -5 on N-term
            # by calling function S11S12_find
            S11 = S11S12_find(upstream_rev)
            S11_collect.append(S11)
            # find S12
            # looks at the secondary structures of the Asn and following +5 on C-term
            # by calling function S11S12_find
            S12 = S11S12_find(downstream)
            S12_collect.append(S12)
        # if the Asn is not placed at the end of the chain and not in a coil
        # then S10,S11 and S12 are 0
        else:
            S10_collect.append(0)
            S11_collect.append(0)
            S12_collect.append(0)
    return S10_collect, S11_collect, S12_collect

#######################  4) DEAMIDATION CALCULATOR  ##########################

def half_time_calculator(S1_collect, S2_collect, S3_collect, S7_collect, S8_collect, S10_collect, S11_collect, S12_collect):
    """Function outputs a list of deamination half times for each glutamine and asparagine and
    a list of with information that can be converted to degradation
    """
    # calculation deamidation for each Gln residues from the primary structure
    ID_lst = [1 / CD for CD in sum_gln_data]

    # calculates the deamination half time for all Asn in the 3D structure
    CD_lst = []
    ASN_LIST_LENGTH = len(asn_pos_lst)
    for element in range(ASN_LIST_LENGTH):
        # if S5 = 0 a division by zero error is avoided by inserting a 1 in the denominator

        f = 0.48 * ( (1.0) * (S1_collect[element]) + (3.1) * (S2_collect[element])
                + 10 * (S3_collect[element])
                + 0.5 * (S7_collect[element]) + 3.2 * (S8_collect[element])
                + 2.0 * (1 - S10_collect[element]) + 0.26 * (5 - S11_collect[element])
                + 0.62 * (5 - S12_collect[element])
                )
        # calculate the individual halftime(CD) of each Asn
        prim_half_time = sum_asn_data[element]
        CD = (0.01) * prim_half_time * math.exp(f)
        CD_lst.append(CD)

        # prepares the individual
        ID_val = 1 / CD
        ID_lst.append(ID_val)
    return ID_lst

##############################  MAIN SCRIPT #################################################

# user defined PDB ID
pdb_id = input('Pleas write a PDB ID: ')
# loads PDB file from the web creates a list of the amino acids in the protein
aa_seq = load_pdb(pdb_id)
# caluclates the primary half time
(asn_pos_lst, sum_asn_data, sum_gln_data) = primary_half_time(aa_seq)
# sends a request to the DSSP program on the computer.
# DSSP predicts the secondary structure this is given as output
dssp_input = dssp_import()
# finds the S7-S8 values and finds structal_data
S7_collect, S8_collect, structure_data = HH_bond(dssp_input)
# finds the S values related to helix structe
S1_collect, S2_collect, S3_collect = helix_S(asn_pos_lst, structure_data)
# finds the S values related to the end of the chain
S10_collect, S11_collect, S12_collect = end_chain(asn_pos_lst, structure_data)

LENGTH = len(asn_pos_lst)
if ( len(S1_collect) != LENGTH and len(S2_collect) != LENGTH and
            len(S3_collect) != LENGTH and len(S7_collect) != LENGTH and
            len(S8_collect) != LENGTH and len(S10_collect) != LENGTH and
            len(S11_collect) != LENGTH and len(S12_collect) != LENGTH ):
    sys.stdout.write('The S-value lists are not the same length\n')
    sys.exit(1)

# calculates the prerequisites for the 3D half time of the protein
ID_lst = half_time_calculator(S1_collect, S2_collect, S3_collect, S7_collect, S8_collect, S10_collect, S11_collect, S12_collect)

# the total deamidation half time of total protein:
try:
    # calculates the half time of the protein if ID list is not zero
    half_time_of_total_protein = 100 / sum(ID_lst)

except ZeroDivisionError as error:
    sys.stdout.write('The peptide has no asparagine residues, therefor the half time cannot be calculated\n')
    sys.exit(1)

print('Protein half time is estimated to ', '%.1f' % half_time_of_total_protein, 'days')
