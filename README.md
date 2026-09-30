## PPI screen: Human-1433Z-dimer x Lpn-effectors

### Steps 
- Get **amino acid sequences** in fasta format 
    - 14-3-3-zeta from human from UniProt  
    - 368 Legionella pneumophial effector IDs (UniProt) from pathogen3d db
        - Download the pdf containing all 368 effector IDs 
        - Convert pdf into csv with AI tools 
            - Check random IDs to verify 
        - Get fasta sequences for all IDs with `fetch_fasta_csv.py`
            - Fasta for each protein and a multifasta is made 
- Make bait:protein **paired-fasta** 
    - Make a bait:bait:protein pair up with `pair_bait_protein.py`
- **MSA** w `colabfold_search` on BIOMIX 
    - Make job script to run MSA (`run_msa_cf_search.sh`)
        - Make local dir on working node to avoid I/O bottleneck 
        - Use fewer threads to avoid needing more memory/diskspace to r/w
    - Reserve hardware and time, and submit the job `submit_biomix_cf_search_msa_no_cpu.slurm`


**Commands**
```bash
# get aa seqs - on local machine 
python3 fetch_fasta_csv.py -i legionella_effectors_368.csv

# make bait:protein pairs - on local machine 
python3 pair_bait_protein.py -i ./proteins/multi_fasta.fasta -bf 1433Z_HUMAN.fasta --bait-copy 2

# submit job on BIOMIX/HPC cluster 
sbatch submit_biomix_cf_search_msa_no_cpu.slurm

# check run status 
squeue -u $USER
tail -f colabfold_search_biomix.log

```

