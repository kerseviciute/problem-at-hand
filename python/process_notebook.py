import json
import os

snakemake_info = {
    "input": snakemake.input,
    "output": snakemake.output,
    "params": snakemake.params,
    "wildcards": snakemake.wildcards
}

filename = f'.{os.path.basename(snakemake.input.notebook)}_snakemake.json'

print(f'Saving the snakemake object to {filename}')

with open(filename, "w") as json_file:
    json.dump(snakemake_info, json_file)

os.system(f'jupyter nbconvert --execute --to html --TagRemovePreprocessor.remove_input_tags=\'{{"hide_code"}}\' {snakemake.input.notebook} --output {snakemake.output.report}')
