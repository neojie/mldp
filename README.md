# mldp
machine learning data preparation

recal protocal
1- convergence test, 1 single point with CONTCAR is enough
2- the previous step also test for single submission for the 
3-test for environment works for both python 3 and vasp
4-check for array submission script
5- recal/post_recal/merge_out/extract_deepmd/pre_train
6-run training
7-freeze model
8-analyze with stat.py 


Rescale training
1-2-3-4
5- recal/post_recal/merge_out/extract_deepmd
6- stat.py , no model mode to get the scale function model
7- scale deepmd
8- pre_train.py scale deepmd
9- stat.py
10-post_stat.py

pbs-lammps training
job_array.py => like dsq


order	code		function		arguments								ATTENTION	
1	recal	recal_relax	submit jobs		ip 		if						check if poscar generated make sense	
2	post_recal	post_recal_relax	check how many submission finished, collect the failed runs, re-submit		ip 			f					check the flag, should not overlap with the existing ones	
3	merge_out		 extract the OUTCAR		ip 	o				r			MUST wait until #2 output flags, START to remove files this step!	
4	extract_deepmd		 extract the deepmd		ip 	o			bs				May experience force assertion error, NOT solved yet!	
5	pre_train		output dp train folder		ip 						d		d could be deepmd, deepmd_relax2, whatever exists, or even deepmd-deepmd_relax2	
					if not, then current folder		vasp inputs	finish flag	batch size	remove file?				
6	pre_train													
7	 													
														
														
														
	mis	remove_recal.py			we randomly generate files, repeatedly generate files, run this to remove the files which have not been submitted yet									
		adjust_outcar.py			when P> 10000 GPa, VASP output misfunctions									
		model_test.py			calculate errors with model			output all errors in as a file						
		std_vs_t2.py			calculate std			output vol, sigma, e,f,v, stds						
		dp_test.py			function used by std_vs_t2									
