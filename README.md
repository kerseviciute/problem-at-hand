# Adaptive neurotechnologies

## How to run this analysis

1. Install conda (and mamba) and create the _snakemake_ conda environment. Activate the environment.

```shell
mamba env create -f env/snakemake.yml
conda activate snakemake
```

2. Run ``snakemake``:

```shell
snakemake --cores 1 --conda-frontend mamba --use-conda
```

## Sample sheet information

- **SID**: File identifier, used to locate the data files.
- **SampleID**: Sample identifier, used to determine the person.
- **URL**: URL used to download the data.
- **Age**: Age of the person.
- **Sex**: Sex of the person.
- **DominantHand**: Dominant hand of the person.

## Data aquisition

- **Product**: g.Pangolin (g.tec medical engineering GmbH, Austria)
- **Number of channels**: 256 channels
- **Sampling rate**: 600Hz

## To Do

- [ ] Data filtering reports
- [ ] Feature extraction reports
- [ ] Finalize the classifier analysis
- [ ] Fit models for all samples
- [ ] Which frequencies are the best for finger classification? -> create topographic maps using these frequencies?
- [ ] Create figures for the paper
  - [ ] Data preprocessing
  - [ ] Feature extraction
  - [ ] Model results
  - [ ] Topographic maps
