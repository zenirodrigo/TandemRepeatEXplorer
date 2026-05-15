awk '
/^>/ {
    if (seq) {
        print seq > outfile
        close(outfile)
    }
    header = substr($0,2)
    outfile = header ".fasta"
    print ">" header > outfile
    seq=""
    next
}
{
    seq = seq $0
}
END {
    if (seq) {
        print seq > outfile
        close(outfile)
    }
}
' satdnas.fasta
