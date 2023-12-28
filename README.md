# Problem at Hand

## How to run this analysis

1. Install conda and mamba. Create the _snakemake_ (`env/snakemake.yml`) conda environment
(required to run the analysis) and the _problem-at-hand-mne_ (`env/mne.yml`) conda environment
(required to generate the reports).

```shell
mamba env create -f env/snakemake.yml
mamba env create -f env/mne.yml
```

2. Activate the snakemake environment.

```shell
conda activate snakemake
```

3. Run ``snakemake``:

```shell
snakemake --conda-frontend mamba --use-conda --cores 1 -p all
```

## Information about the data

The original dataset (available [here](https://osf.io/4dwjt/?view_only=d23acfd50655427fbaae381a17cbfbcc))
had data for 5 persons, one of which was left-handed. In this analysis, only the data from
4 subjects, all right-handed, were used.

### Sample sheet information

- **SID**: File identifier, used to locate the data files.
- **SampleID**: Sample identifier, used to determine the person.
- **URL**: URL used to download the data.
- **Age**: Age of the person.
- **Sex**: Sex of the person.
- **DominantHand**: Dominant hand of the person.

### Data acquisition

- **Product**: g.Pangolin (g.tec medical engineering GmbH, Austria)
- **Number of channels**: 256 channels
- **Sampling rate**: 600Hz

## To Do

- [x] Data filtering reports
- [ ] Feature extraction reports
- [ ] Finalize the classifier analysis
- [ ] Fit models for all samples
- [ ] Which frequencies are the best for finger classification? -> create topographic maps using these frequencies?
- [ ] Create figures for the paper
  - [ ] Data preprocessing
  - [ ] Feature extraction
  - [ ] Model results
  - [ ] Topographic maps
- [ ] Try to deny the results.
