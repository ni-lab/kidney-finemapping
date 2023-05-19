conda activate kidney_finemapping

cd /home/rshuai/research/ni-lab/kidney_finemapping/kidney_finemapping || exit

python3 basenji/plot/plot_sad_tracks.py \
/home/rshuai/research/ni-lab/analysis/kidney_data/220513_variants/sad_pos_shifts/all_chrs/sad.h5 \
-t /home/rshuai/research/ni-lab/kidney_finemapping/kidney_finemapping/data/targets/kidney_sc_samples.txt \
--overlay \
--overlay_lines_only \
-o /home/rshuai/research/ni-lab/kidney_finemapping/kidney_finemapping/out_dir/220513_variants/plots/sad_tracks_overlay_lines_only
