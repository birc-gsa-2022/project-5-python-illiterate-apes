def fastq_parser(file):
    out = []
    for l in (l for l in file):
        if l.startswith('@'):
            out.append([l[1:].strip()])
        else:
            out[-1].append(l.strip())
    return out