import argparse
import sys
import pickle
import fasta, fastq
import os
import re
from collections import defaultdict

class BWTMatcher:
    def __init__(self, f, rank_table, firstIndexList, alphadic):
        self.rank_table = rank_table
        self.f = f
        self.firstIndexList = firstIndexList
        self.alphadic = alphadic

def main():
    argparser = argparse.ArgumentParser(
        description="Readmapper",
        usage="\n\treadmap -p genome\n\treadmap -d dist genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "-d", type=int, metavar="integer",
        default=1, help="max edit distance."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    if not(args.p):
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        
    if args.genome is None:
        argparser.print_help()
        sys.exit(1)

    if args.p:
        genomes = fasta.fasta_parse(args.genome)
        genomes_to_file(args.genome.name, genomes)
    else:
        k = args.d
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        
        genomes = fasta.fasta_parse(args.genome)
        bwtList = None
        # Check if we have a .dat file
        datFile = args.genome.name+".dat"
        if os.path.isfile(datFile):
            datFileStream = open(datFile, "rb")
            bwtList = pickle.load(datFileStream)
        else:
            bwtList = genomes_to_file(args.genome.name, genomes)
        
        reads = fastq.fastq_parser(args.reads)

        out = []

        for i, g in enumerate(genomes):
            if len(g[1]) == 0:
                continue
            for r in reads:
                length = len(r[1])
                if length == 0:
                    continue

                matches = searchPattern(r[1], bwtList[i], k)
                for m in matches:
                    out.append((getTrailingNumber(r[0]), getTrailingNumber(g[0]), m[0], m[1], r[1]))

        for t in out:
            print(f"{t[0][0]}{t[0][1]}\t{t[1][0]}{t[1][1]}\t{t[2]}\t{t[3]}\t{t[4]}")

def getTrailingNumber(s):
    m = re.search(r'\d+$', s)
    return (s[:m.start()], int(s[m.start():]))

def compactCigar(u_cigar):
    c_cigar = ""
    prevChar = ""
    counter = 0
    for c in u_cigar:
        if c==prevChar:
            counter += 1
        else:
            if counter > 0:
                c_cigar += str(counter)+prevChar
            counter = 1
        
        prevChar = c
    
    if counter > 0:
        c_cigar += str(counter)+prevChar

    return c_cigar

def genomes_to_file(filename, genomes):
    bwtList = preprocess_genomes(genomes)
    outputFile = open(str(filename)+".dat", "wb")
    pickle.dump(bwtList, outputFile)
    return bwtList

def preprocess_genomes(genomes):
    bwtList = []
    for gen in genomes:
        string = gen[1]+"$"
        alphadic = {a: i for i, a in enumerate(set(string))}

        x = memoryview(string.encode())
        suf = getSuffixes(x)
        f = radix_sort(suf)
        bwt = [(i-1)%len(f) for i in f]
        rank_table = build_rank_table(x, alphadic, bwt)
        firstIndexList = getFirstIndexList(x, f, alphadic)
        bwtList.append(BWTMatcher(f, rank_table, firstIndexList, alphadic))
    return bwtList

def getSuffixes(x):
    """
    Gets all suffixes from a string
    """
    suffixes = [x[i:] for i in range(0, len(x))] 
    # print("list of unordered suffixes: ",suffixes) 
    return suffixes

def radix_sort(lst: list[memoryview]):
    # print("Radix sort", len(lst))
    maximum = max(lst, key = len)

    for place in range(len(maximum),0, -1):
        lst = counting_sort(lst, place - 1)
    
    lst = [len(maximum) - len(suffix) for suffix in lst]
    return lst

def counting_sort(lst: list[memoryview], place):
    maximum = max(lst, key = len)
    
    counts = defaultdict(list)
    for key in ["$",*sorted(maximum)]:
        counts[key] = []

    for string_index in range(len(lst)):
        if place >= len(lst[string_index]):
            counts["$"].append(lst[string_index])
        else:
            counts[lst[string_index][place]].append(lst[string_index])

    ordered = []
    for key in counts:
        ordered += counts[key]
    return ordered

def build_rank_table(x, alphadic, bwt):
    table = [[0 for _ in alphadic] for _ in range(0, len(bwt)+1)]

    for i in range(1, len(bwt)+1):        
        for j in range(0, len(alphadic)):
            table[i][j] = table[i-1][j]
        
        bwtValue = bwt[i-1]
        c = chr(x[bwtValue])

        index = alphadic.get(c)
        table[i][index] += 1
    return table

def getFirstIndexList(x, f, alphadic):
    firstIndexList = {a: -1 for a in alphadic}
    indexesFound = 0

    for i, xIndex in enumerate(f):   
        c = chr(x[xIndex])
        if firstIndexList.get(c) < 0:
            firstIndexList[c] = i
            indexesFound += 1

            if indexesFound >= len(alphadic):
                return firstIndexList
    
    return firstIndexList

def getrank(alphadic, index, c, rank_table):
    return rank_table[index][alphadic.get(c)]

def get_d_table(p, bwtMatcher):
    d_table = [0 for _ in range(len(p))]
    edits = 0

    left, right = 0, len(bwtMatcher.f)
    for i, c in enumerate(reversed(p)):
        left, right = fm_match(bwtMatcher, left, right, c)
        if left >= right:
            edits += 1
            d_table[i] = edits
            left, right = 0, len(bwtMatcher.f)
        else:
            d_table[i] = edits
    
    return d_table

class EditNode:
    def __init__(self, index, edits, cigar, left, right):
        self.index = index
        self.edits = edits
        self.cigar = cigar
        self.left = left
        self.right = right

    def __str__(self):
        return str(self.index) + ' ' +  str(self.edits) + ' ' + self.cigar + ' ' + str(self.left) + ' ' + str(self.right)

def fm_match(bwtMatcher, left, right, c):
    if c == "$" or c not in bwtMatcher.alphadic:
        return left, left
    left = bwtMatcher.firstIndexList[c] + bwtMatcher.rank_table[left][bwtMatcher.alphadic.get(c)]
    right = bwtMatcher.firstIndexList[c] + bwtMatcher.rank_table[right][bwtMatcher.alphadic.get(c)]
    return left, right

def searchPattern(p, bwtMatcher, k):
    if p == "":
        return
    
    d_table = get_d_table(p, bwtMatcher)
   

    if k < d_table[len(p)-1]: return

    # We need index, number edits, cigar, left, right
    stack = [EditNode(len(p)-1, k, "", 0, len(bwtMatcher.f))]

    sols = []
    while stack:
        node = stack.pop()

        if node.edits < 0:
            continue
        if node.index < 0:
            for i in range(node.left, node.right):
                compactedCigar = compactCigar(node.cigar)
                yield [bwtMatcher.f[i]+1, compactedCigar]
            continue

        if node.edits < d_table[node.index]:
            continue

        # Insertion
        newNode = EditNode(node.index-1, node.edits-1, 'I'+node.cigar, node.left, node.right)
        stack.append(newNode)

        # Deletion
        if node.cigar != "":
            for c in bwtMatcher.alphadic:
                left, right = fm_match(bwtMatcher, node.left, node.right, c)

                if left < right:
                    newNode = EditNode(node.index, node.edits-1, 'D'+node.cigar, left, right)
                    stack.append(newNode)

        # Match/Substitution for each char
        for c in bwtMatcher.alphadic:
            left, right = fm_match(bwtMatcher, node.left, node.right, c)
            
            if left < right:
                if c == p[node.index]:
                    # Matching
                    newNode = EditNode(node.index-1, node.edits, 'M'+node.cigar, left, right)
                    stack.append(newNode)
                else:
                    # Substitution
                    newNode = EditNode(node.index-1, node.edits-1, 'M'+node.cigar, left, right)
                    stack.append(newNode)
        
    return sols

if __name__ == '__main__':
    main()