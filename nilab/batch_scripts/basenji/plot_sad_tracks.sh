conda activate kidney_finemapping

cd /home/rshuai/research/ni-lab/kidney_finemapping/nilab || exit

python3 basenji/plot/plot_sad_tracks.py /home/rshuai/research/ni-lab/analysis/kidney_data/220513_variants/sad_pos_shifts/all_chrs/sad.h5 -t /home/rshuai/research/ni-lab/kidney_finemapping/nilab/data/targets/kidney_sc_samples.txt -o /home/rshuai/research/ni-lab/kidney_finemapping/nilab/out_dir/220513_variants/plots/sad_tracks
