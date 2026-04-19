#!/usr/bin/env bash

for tile_dim in {32..1024..32}
do
	echo "Running tile ${tile_dim}"
	echo "${tile_dim} ${tile_dim} ${tile_dim}" > tiles/dense_1024
	./scache dense_1024 dense_1024 config/sweep_config.json
	mv output_sweep/CGustBase_2.000000MB_68.000000GBs_64PEs_64sbanks__dense_1024_dense_1024_RR_0_data_4_coord_4.txt "output_sweep/${tile_dim}.txt"
done
