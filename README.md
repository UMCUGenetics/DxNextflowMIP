# DxNextflowMIP
Genome Diagnostics Nextflow MIP Fingerprint workflow

## Get Nextflow Modules
```bash
git submodule update --init --recursive
```

## Install Nextflow
```bash
mkdir tools && cd tools
curl -s https://get.nextflow.io | bash
```

## Running MIP fingerprint workflow
```bash
nextflow run MipFingerprint.nf -c MipFingerprint.config --fastq_path <fastq_dir_path> --outdir <output_dir_path> --email <email> [-profile slurm|mac]
```
