# PyRosetta_TACC_MPI module


### class PyRosetta_TACC_MPI.AbstractExperimentRunner(start_pose_pdbs, rosetta_init_options, comm, restrict_to_chain, max_pack_rounds, min_cst_sd, min_restrict_radius, PDB_res, dump_ref_pdb, dump_mut_pdb, pdb_base, verbose, constraint_file)
Bases: `object`


#### dump_ddg_results(outfile_name, header=True)
Specification Required!


#### gather_results()

#### run_all_jobs_serial(outfile='ddg_out.txt')

#### scatter_job()

#### setup_MPI(comm)

### class PyRosetta_TACC_MPI.AbstractPackerJob(convergence_fn=None, conv_threshold=0.1, repack_radius=10, scorefn='mm_std', mintype='dfpmin_armijo_nonmonotone', max_rounds=100, restrict_to_chain=None, verbose=1, constraint_file=None)
Bases: `object`


#### add_constraints_to_pose(pose, constraint_file)

#### add_constraints_to_scorefxn(constraint_types=None, weights=None, default_weight=0.1)

#### build_movemap(pose, restrict_radius_center=None, restrict_radius=None, restrict_residues=None)

#### make_CompoundMover(movers, repeats_per_round)

#### make_minmover(mintype, movemap, pose)

#### make_packer_task_with_residues(pose, residues=None)
Builds a packer task with the specified residues activated for repacking.

Did this to avoid PackerTask.temporarily_\* methods which apparently we’re not supposed to use


#### make_packmin_mover(pose, packertask, minmover, n_packing_steps, n_minimize_steps, kT=1.0)
Build a TrialMover with n_packing_steps RotamerTrialMover moves and

    n_minimization_steps MinMover moves executed sequentially


#### mutate_aa(pose, residue, aa_name, orient_bb=True, repack_sidechain=True, clone_pose=True)
Swap w/t AA at residue number ‘residue’ in ‘pose’ with ‘ncaa_name’ (3-letter code)

Return a new Pose object

Note that this assumes the ncaa .params and .rotlib files have been permanently added to the database


#### packer_task_repack_in_radius(pose, residue, radius)
Build a packer task that repacks all residues with CAlphas within distance radius

Might want to remake to take centroid side-chains?  Looks like SwitchResidueTypeMover has trouble w/ nsAAs…


#### residue_CAs_in_radius(pose, centerAA, radius)
Get a list of residues with C-alpha atoms within ‘radius’ distance from centerAA’a C-alpha


#### std_dev_threshold_fn_builder(threshold, last_n=5)

### class PyRosetta_TACC_MPI.FastRelaxPoseJob(in_pdb, MPI_rank, scorefn, restrict_to_chain=None, constraint_file=None)
Bases: `PyRosetta_TACC_MPI.AbstractPackerJob`


#### dump_pose(outfile)

#### pack_pose()

### class PyRosetta_TACC_MPI.MultiMutantPackerJob(start_pose_pdb, mutation_list, replicate, convergence_fn=None, conv_threshold=0.1, repack_radius=10, scorefn='mm_std', mintype='dfpmin_armijo_nonmonotone', n_pack_steps=3, n_min_steps=1, max_rounds=100, min_cst_sd=None, min_restrict_radius=False, restrict_to_chain=False, PDB_res=False, verbose=1, constraint_file=None)
Bases: `PyRosetta_TACC_MPI.AbstractPackerJob`


#### dump_mut_pdb(outfile)

#### dump_ref_pdb(outfile)

#### get_result()

#### run()

#### setup_poses()

### class PyRosetta_TACC_MPI.MutagenesisExperimentRunner(start_pose_pdbs, rosetta_init_options, comm, residue_list=[], AA_list=[], nreps=50, restrict_to_chain=False, max_pack_rounds=25, min_cst_sd=None, min_restrict_radius=False, PDB_res=False, dump_ref_pdb=False, dump_mut_pdb=False, pdb_base='', verbose=1, constraint_file=None)
Bases: `PyRosetta_TACC_MPI.AbstractExperimentRunner`


#### build_job_list(residues, AAs, replicates, shuffle_jobs=True)

#### dump_ddg_results(outfile_name, header=True)
Dump results to a .tsv file


#### setup_jobs(residue_list, AA_list, nreps)

### class PyRosetta_TACC_MPI.MutantCombinationsExperimentRunner(start_pose_pdbs, rosetta_init_options, comm, mutant_list=[], nreps=50, restrict_to_chain=False, max_pack_rounds=25, min_cst_sd=None, min_restrict_radius=False, PDB_res=False, dump_ref_pdb=False, dump_mut_pdb=False, pdb_base='', verbose=1, constraint_file=None)
Bases: `PyRosetta_TACC_MPI.AbstractExperimentRunner`


#### build_job_list(mutant_list, replicates, shuffle_jobs=True)

#### dump_ddg_results(outfile_name, header=True)
Dump results to a .tsv file


#### setup_jobs(mutant_list, nreps)

### class PyRosetta_TACC_MPI.MutantddGPackerJob(start_pose_pdb, residue, chain, AA, replicate, convergence_fn=None, conv_threshold=0.1, repack_radius=10, scorefn='mm_std', mintype='dfpmin_armijo_nonmonotone', n_pack_steps=3, n_min_steps=1, max_rounds=100, min_cst_sd=None, min_restrict_radius=False, restrict_to_chain=False, PDB_res=False, verbose=1, constraint_file=None)
Bases: `PyRosetta_TACC_MPI.AbstractPackerJob`


#### dump_mut_pdb(outfile)

#### dump_ref_pdb(outfile)

#### get_result()

#### run()

### PyRosetta_TACC_MPI.NSAAS_PATCH( = {})
‘ACK’:{‘cognateAA’:’LYS’,
‘type’:chemical.VariantType.ACETYLATION},
‘PHS’:{‘cognateAA’:’SER’,
‘type’:chemical.VariantType.PHOSPHORYLATION},
‘PHT’:{‘cognateAA’:’THR’,
‘type’:chemical.VariantType.PHOSPHORYLATION},
‘PHY’:{‘cognateAA’:’TYR’,
‘type’:chemical.VariantType.PHOSPHORYLATION}}


### class PyRosetta_TACC_MPI.PackSinglePoseJob(in_pdb, MPI_rank, convergence_fn=None, conv_threshold=0.1, repack_radius=10, scorefn='mm_std', mintype='dfpmin_armijo_nonmonotone', n_pack_steps=1, n_min_steps=1, max_rounds=100, restrict_to_chain=None, kT=1.0, MCmover=True, verbose=1, constraint_file=None)
Bases: `PyRosetta_TACC_MPI.AbstractPackerJob`


#### dump_pose(outfile)

#### pack_pose()

### exception PyRosetta_TACC_MPI.PyRosettaError()
Bases: `Exception`


### class PyRosetta_TACC_MPI.SecondaryMutantScanExperimentRunner(start_pose_pdbs, rosetta_init_options, comm, center_residue, center_residue_AA, center_residue_chain, mutate_radius, center_res_ref_pose, mutate_secondary_to_AAs, nreps=50, restrict_to_chain=False, max_pack_rounds=25, min_cst_sd=None, min_restrict_radius=False, PDB_res=False, dump_ref_pdb=False, dump_mut_pdb=False, pdb_base='', verbose=1)
Bases: `PyRosetta_TACC_MPI.MutantCombinationsExperimentRunner`


#### residue_CAs_in_radius(center_pose, center_chain, centerAA, radius)
Get a list of residues with C-alpha atoms within ‘radius’ distance from centerAA’a C-alpha
