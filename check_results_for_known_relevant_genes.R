# 2018-09-07
# Katie Saund

# Check the results of the three GWAS methods for genes known to be relevant to 
# the phenotypes. This serves as a sanity check for the method and will tell us
# if we can trust our results or not. 

# Load in and sort genes into groups by phenotype ------------------------------
# these lists below are cobbled together from: 
# /nfs/esnitkin/Project_Cdiff/Analysis/Project_crystallography/2018-07-12_curated_genes/data/
# and
# 
expected_sporulation_genes <- c("CD1352", "CD1492", "CD1579", "CD1949", 
                                "CD2492", "Spo0A", "SpoIIAA", "SpoIIAB", 
                                "SpoIIE", "Sigma factor F", "Sigma factor G", 
                                "Sigma factor E", "Sigma factor K", "SpoIVA", 
                                "SpoVM", "SipL", "CdeC", "BclA1", "bclA2", 
                                "bclA3", "CD3563", "CspBA", "SleC", "CspC", 
                                "CotA", "CotB", "CotCB", "CotD", "CotE", "CotF", 
                                "CotJB2", "CotG", "SodA", "cheB", "sigE", 
                                "sigF", "sigG", "spoIIAA", "spoIIAB", "spoIVB", 
                                "gpr", "spoIIE", "spoIIR", "spoIVA", "spoIVB2")

expected_germination_genes <- c("gerP", "cwlJ", "gerKA", "gerKB", "gerKC", 
                                "gerLA", "gerLB", "gerLC", "spoVA", "gerS", 
                                "cspA", "cspC", "cspB", "cspBA", "sipL", "ssb", 
                                "CD3571", "asd", "CD2894", "spoIIIAA", "CD0972")

expected_toxin_genes       <- c("rho", "tpi", "CD2808", "tcdA", "tcdC", "codY", 
                                "ldh", "polA", "CD1263", "tcdE", "tcdB", "cdtR", 
                                "ermB", "gluD", "gyrA", "CD2894A")
expected_fqR_genes         <- c("gyrA", "gyrB")

#expected_growth_genes      <- c()


# Load in GWAS results for 3, 2+, or 1+ methods by phenotype -------------------



# Convert LJPAEC### locus_tags into CD630_##### locus tags? --------------------

# Search GWAS results for expected genes ---------------------------------------


# Print table of results and visualize -----------------------------------------


# End of script ----------------------------------------------------------------
