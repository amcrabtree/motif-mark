import argparse
import cairo
import seaborn
import re

######################### ARGPARSE #########################
parser = argparse.ArgumentParser(description='Visualize sequence motifs on DNA sequences.')
parser.add_argument('-f', '--fasta', help='fasta file (required)')
parser.add_argument('-m', '--motifs', help='motif file')
args = parser.parse_args()

###################### GLOBAL FUNCTIONS ######################

def store_motif_set(txt_file: str) -> set:
    '''Store motif sequences from a file as a set.'''
    motifs = set()
    with open(txt_file, "r") as motif_h:
        for line in motif_h:
            line = line.strip("\n")
            motifs.add(line)
    return motifs

def store_gene_dict(fasta_file: str) -> dict:
    gene_dict = {}
    with open(fasta_file, "r") as gene_h:
        gene_seq = ''
        gene_header = ''
        i=0
        for line in gene_h:
            i+=1
            line = line.strip()
            if line.startswith('>'): # if header line
                if i!=1:
                    gene_dict[gene_header] = gene_seq
                gene_header = line.strip('>')
                gene_seq = ''
            else:
                gene_seq = ''.join([gene_seq, line])
                if line != "\n": # if last line in file
                    gene_dict[gene_header] = gene_seq
    return gene_dict

def gen_exon_dict(seq: str) -> dict:
    '''Generates a dictionary of all exon positions and lengths, using uppercase to denote exons.'''
    # exon dictionary (key=exon_id, value=[exon_position, exon_length])
    exon_dict = {}
    exon_matches = [m for m in re.finditer('[A-Z]+', seq)]
    for i in range(len(exon_matches)):
        ex_id = 'exon' + str(i+1)
        ex_pos = exon_matches[i].start()+1
        ex_len = len(exon_matches[i].group())
        exon_dict[ex_id] = [ex_pos, ex_len]
    return exon_dict

def degen_pattern(seq: str) -> str:
    '''Generates a string used for regex recognition of all possible dna sequences''' 
    '''from one degenerate sequence.'''
    # degenerate dictionary
    degen_dict = {
        'U': 'T', # included in case motif is in RNA form
        'M': '[AC]',
        'R': '[AG]',
        'W': '[AT]',
        'S': '[CG]',
        'Y': '[CT]',
        'K': '[GT]',
        'V': '[ACG]',
        'H': '[ACT]',
        'D': '[AGT]',
        'B': '[CGT]',
        'N': '[ACGT]'
    }
    # create output list of all possible sequences
    seq = seq.upper() # convert sequence to uppercase
    pattern = seq
    # store set of degenerate nucleotides
    degen_matches = {x.group() for x in re.finditer('[^ATGC]', seq)}
    if degen_matches:  # if there ARE degenerate nts,
        for degen_nt in degen_matches:
            degen_pattern = degen_dict[degen_nt]  # store degenerate pattern
            pattern = pattern.replace(degen_nt, degen_pattern)
    return pattern

def draw_gene(context, gene_obj, dist_from_top=60, dist_from_left=50):
    context.select_font_face("Ariel", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    context.set_font_size(30)
    context.move_to(20, dist_from_top-10)
    context.set_source_rgba(0,0,0, 0.8)
    gene_title = re.sub(r'chr\S+', "", gene_obj.description) #remove chr position from title
    context.show_text(gene_title)

    # draw horizontal line (lay genome track)
    x_final = len(gene_obj.seq)+dist_from_left
    if x_final > 1000:
        print("CAUTION: You have a sequence > 1000 bp, so it will run off the page.")
    x, neg_y, neg_y_final = [dist_from_left, dist_from_top+40, dist_from_top+40]
    context.set_line_width(5)
    context.set_line_cap(cairo.LINE_CAP_ROUND)
    context.move_to(x, neg_y)
    context.line_to(x_final, neg_y_final)
    context.set_source_rgb(0,0,0)
    context.stroke()

    # draw rectangle (gene element)
    ex_pos, ex_len = [x for x in gene_obj.exons.values()][0]
    x, neg_y, wid, hig = [ex_pos+dist_from_left, dist_from_top, ex_len, 80]
    context.set_line_width(0.02)
    context.rectangle(x, neg_y, wid, hig)
    context.set_source_rgba(0,0,0, 0.8)
    context.fill()
    context.stroke()

def draw_motifs(context, motif_dict, dist_from_top=60, dist_from_left=50):
    # create list of colors contingent with number of motifs present
    palette = [list(x) for x in seaborn.color_palette("husl", len(motif_dict))]
    p=0 # initialize palette counter
    # draw rectangle of unique color for each motif type
    for motif in motif_dict:
        r,g,b = palette[p] # set motif color
        for motif_pos in motif_dict[motif]:
            motif_len = len(motif)
            x, neg_y, wid, hig = [motif_pos+dist_from_left, dist_from_top, motif_len, 80]
            context.set_line_width(0.02)
            context.rectangle(x, neg_y, wid, hig)
            context.set_source_rgba(r, g, b, 0.8)
            context.fill()
            context.stroke()
        p+=1

def draw_legend(context, motif_dict, dist_from_top=60, dist_from_left=1070):
    # create list of colors corresponding to motifs
    palette = [list(x) for x in seaborn.color_palette("husl", len(motif_dict))]
    p=0 # initialize palette counter
    for motif in motif_dict:
        # draw square of unique color for each motif type
        r,g,b = palette[p] # motif color
        context.set_source_rgba(r, g, b, 0.8)
        p+=1
        x, neg_y, wid, hig = [dist_from_left, dist_from_top-10, 40, 40]
        context.set_line_width(0.02)
        context.rectangle(x, neg_y, wid, hig)
        context.fill()
        context.stroke()
        # draw text
        context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        context.set_font_size(30)
        context.move_to(dist_from_left+50, dist_from_top+20)
        context.set_source_rgba(0,0,0, 0.8)
        context.show_text(motif)
        dist_from_top+=50
        
######################### CLASSES #########################

class Motif:
    def __init__(self, id, seq):
        '''This stores motifs.'''
    ## Data ##
        self.id = id
        self.seq = seq
        self._owner = None
    ## Methods ##
    
class Gene:
    def __init__(self, description, seq):
        '''This stores gene info.'''
    ## Data ##
        self.description = description
        self.seq = seq
        self.exons = gen_exon_dict(seq)
        self._owner = None
    ## Methods ##
    def motif_dict (self, motif_set: set) -> dict:
        '''Store dictionary with all motifs found within gene.'''
        mdict = {}
        for motif in motif_set:
            pattern = degen_pattern(motif)
            match_list = [x for x in re.finditer(pattern, self.seq.upper())]
            positions_list = [p.start()+1 for p in match_list]
            mdict[motif] = positions_list # key=motif_seq, value=nt_positions
        #print(f'Motifs found in {self.description}:')
        #for k in mdict:
        #    print(k,":",mdict[k])
        return mdict

######################### MAIN #########################

## open motif file and store all motifs seqs as a set
my_motif_set = store_motif_set(args.motifs)

## open fragment file and store genes as dictionary elements
my_gene_dict = store_gene_dict(args.fasta)

## convert gene dictionary to list of gene objects
gene_objects = []
for key, val in my_gene_dict.items():
    gene_objects.append(Gene(key, val))

## write cairo img file    
with cairo.ImageSurface(cairo.FORMAT_ARGB32, 1400, 1600) as surface:
    context = cairo.Context(surface)
    d = 60 # distance from top of page to start drawing genes
    for g in gene_objects:
        draw_gene(context, g, d)
        my_motif_dict = g.motif_dict(my_motif_set)
        draw_motifs(context, my_motif_dict, d)
        d+=140 # move next gene down the page
        draw_legend(context, my_motif_dict)
    # save as png
    png_name = args.fasta.split('.')[0]+'.png'
    print('wrote',png_name,'to file.')
    surface.write_to_png(png_name)



