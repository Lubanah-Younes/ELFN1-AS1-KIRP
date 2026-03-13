import gzip

print("Scanning for ELFN1-AS1 (ENSG00000238117)...")

found = False
with gzip.open("tcga_RSEM_gene_tpm.gz", 'rt') as f:
    header = f.readline()  # Read header
    print("Header read successfully")
    
    for i, line in enumerate(f):
        gene = line.split('\t')[0]
        if i % 10000 == 0:
            print(f"Checked {i} genes...")
        
        if gene.startswith("ENSG00000238117"):
            print(f"✅ FOUND GENE at line {i}: {gene}")
            found = True
            break

if not found:
    print("❌ Gene NOT found in file")