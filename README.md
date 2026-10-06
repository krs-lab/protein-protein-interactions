## Protein-protein interaction

Screening interactions between Human 1433Z-dimer & Legionella pneumophila effectors computationally. 

**Steps** 
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
        - Use fewer threads to avoid needing more memory/diskspace to r/w
    - Reserve hardware and time, and submit the job `submit_biomix_cf_search_msa_no_cpu.slurm`
- Generate AlphaFold **multimer models** w `colabfold_batch`
    - Make job script to feed the MSA files to af_multimer `run_multimer_cf_batch.sh`
    - Submit job to teh HPC cluster with `submit_biomix_cf_batch_multimer.slurm`

**Commands**
```bash
# get aa seqs - on local machine 
python3 fetch_fasta_csv.py -i legionella_effectors_368.csv

---

# make bait:protein pairs - on local machine 
python3 pair_bait_protein.py -i ./proteins/multi_fasta.fasta -bf 1433Z_HUMAN.fasta --bait-copy 2

# submit job on BIOMIX/HPC cluster 
sbatch submit_biomix_cf_search_msa_no_cpu.slurm

# check run status 
squeue -u $USER
tail -f colabfold_search_biomix.log

# check output file (.a3m) integrity 
for f in *.a3m; do
  # Check if empty
  if [ ! -s "$f" ]; then
      echo "EMPTY: $f"
      continue
  fi

  # Check if a sequence header exists anywhere in the file
  if ! grep -q "^>" "$f"; then
      echo "INVALID FORMAT (No '>' headers): $f"
      continue
  fi

  # Count MSA depth (number of sequences)
  seq_count=$(grep -c "^>" "$f")
  if [ "$seq_count" -lt 2 ]; then
      echo "LOW DEPTH ($seq_count sequences): $f"
  fi
done

# check MSA depth for bait and prey 
echo "filename,file_size_bytes,sequence_count_bait,sequence_count_protein" > msa_summary.csv

for f in *.a3m; do
  awk -v f="$f" -v size="$(stat -c%s "$f")" '
    function tally(   a, b) {
      gsub(/[a-z]/, "", seq)            # drop lowercase insertions
      a = substr(seq, 1, L); b = substr(seq, L + 1)
      gsub(/-/, "", a); gsub(/-/, "", b)
      nb += (a != "" && b == "")        # residues only in bait part
      np += (a == "" && b != "")        # residues only in protein part
    }
    NR == 1 { split($0, h, /[#,\t]/); L = h[2]; next }   # bait length from header
    /^>/    { tally(); seq = ""; next }
    { seq = seq $0 }
    END { tally(); printf "%s,%s,%d,%d\n", f, size, nb, np }
  ' "$f" >> msa_summary.csv
done

cat msa_summary.csv
```

