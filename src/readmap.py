import argparse
import sys
import pickle
import fasta, fastq
import os
import re

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

    if args.p:
        print(f"Preprocess {args.genome}")
    else:
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        print(
            f"Search {args.genome} for {args.reads} within distance {args.d}")
        
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
            bwtList = genomes_to_file(args.genome, genomes)
        
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
                    out.append((getTrailingNumber(r[0]), getTrailingNumber(g[0]), m+1, length, r[1]))

        for t in sorted(out, key=lambda x: (x[0], x[1], x[2])):
            print(f"{t[0][0]}{t[0][1]}\t{t[1][0]}{t[1][1]}\t{t[2]}\t{t[3]}M\t{t[4]}")

def getTrailingNumber(s):
    m = re.search(r'\d+$', s)
    return (s[:m.start()], int(s[m.start():]))

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
    d_table = [0 * len(p)]
    edits = 0

    left, right = 0, len(bwtMatcher.f)
    for i, c in enumerate(reversed(p)):
        left = bwtMatcher.firstIndexList[c] + bwtMatcher.rank_table[left][bwtMatcher.alphadic.get(c)]
        right = bwtMatcher.firstIndexList[c] + bwtMatcher.rank_table[right][bwtMatcher.alphadic.get(c)]
        if left >= right:
            edits += 1
            d_table[i] = edits
            left, right = 0, len(bwtMatcher.f)
        else:
            d_table[i] = edits
    
    return d_table

class EditNode:
    def __init__(self, index, edits, p_prime, left, right):
        self.index = index
        self.edits = edits
        self.p_prime = p_prime
        self.left = left
        self.right = right

def searchPattern(p, bwtMatcher, k):
    if p == "":
        return
    
    d_table = get_d_table(p, bwtMatcher)

    # if k < d_table[len(p)-1]: return

    # We need index, number edits, p', left, right
    stack = [EditNode(len(p)-1, k, "", 0, len(bwtMatcher.f))]

    #TODO: Check when do we need to build the CIGAR
    while stack:
        node = stack.pop()
        if node.edits < 0 or node.edits < d_table[len(d_table)-node.index]: # TODO: Check conditions
            continue
        if node.index < 0:
            report

        # Deletion
        stack.append(EditNode(node.index-1, k-1, node.p_prime, node.left, node.right))

        # Addition
        for c in bwtMatcher.alphadic:
            if c == p[node.index]:
                stack.append(EditNode(node.index-1, k-1, node.p_prime+p[node.index], node.left, node.right))
            else:
                stack.append(EditNode(node.index, k-1, node.p_prime+p[node.index], node.left, node.right))
        # Match/Substitution for each char
        for c in bwtMatcher.alphadic:
            left = node.left, right = node.right
            left = bwtMatcher.firstIndexList[c] + bwtMatcher.rank_table[left][bwtMatcher.alphadic.get(c)]
            right = bwtMatcher.firstIndexList[c] + bwtMatcher.rank_table[right][bwtMatcher.alphadic.get(c)]
            
            if left >= right:
                # Substitution
                left = 0
                right = len(bwtMatcher.f)
                stack.append(EditNode(node.index-1, k-1, node.p_prime+=c, left, right))
            else:
                # Matching
                stack.append(EditNode(node.index-1, k, node.p_prime+c, left, right))

if __name__ == '__main__':
    main()
