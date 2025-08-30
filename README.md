# NW4BIBMS

**Needleman-Wunsch Implementation for BIBMS**

NW4BIBMS is a Shiny web application designed to align two protein sequences using the **Needleman-Wunsch algorithm**. The app uses the **BLOSUM62** substitution matrix (downloaded from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/matrices/)) with gap opening and extension penalties of **-11** and **-1**, respectively, emulating NCBI's defaults for the *Needleman-Wunsch Global Align Protein Sequences* tool.

## Features

- **Global Protein Alignment:** Computes the Needleman-Wunsch alignment between two user-provided protein sequences.
- **Alignment Score:** Reports the alignment score for the submitted sequences.
- **Randomized Control Scores:** Generates 200 alignments of randomly shuffled versions of the input sequences to provide a background distribution of NW scores.
- **Downloadable Results:** Allows users to download the list of random alignment scores as a CSV file for further analysis.

## Intended Use

This app is intended for **undergraduate students in the BIBMS program** to help them understand the principles of sequence alignment. It can be used as part of course exercises or practical demonstrations in bioinformatics lectures.

## Usage

1. Enter two protein sequences in the provided text boxes.
2. Click **Run Alignment** to compute the alignment and view the Needleman-Wunsch score.
3. Review the NW scores of random shufflings to understand the significance of your alignment.
4. (Optional) Click **Download Random Scores** to export the 200 randomized scores as a CSV file.

## Technical Details

- Built using **R Shiny** for interactive web visualization.
- Uses the **Biostrings** and **pwalign** packages for sequence alignment computations.
- Substitution matrix: **BLOSUM62** (standard amino acids only), downloaded from NCBI.

## License

This project is licensed under the MIT License.
