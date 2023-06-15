table bigLollyExample3
"bigLolly Example 3"
    (
    string   chrom;            "Reference sequence chromosome or scaffold"
    uint     chromStart;       "Start position in chrom"
    uint     chromEnd;         "End position in chrom"
    string   name;             "dbSNP Reference SNP (rs) identifier or <chromNum>:<pos>"
    uint     score;            "Score from 0-1000, derived from p-value"
    char[1]  strand;           "Unused.  Always '.'"
    uint     thickStart;       "Start position in chrom"
    uint     thickEnd;         "End position in chrom"
    uint     color;            "Red (positive effect) or blue (negative). Brightness reflects pvalue"
    double   lollySize;        "Size of lollipop"
    )
