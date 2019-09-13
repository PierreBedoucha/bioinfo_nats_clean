import os
import numpy as np
import Bio.PDB

def read_pdb_starts():
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_starts_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_starts_dict[line[0:4]] = int(line_array[1])
    return pdb_starts_dict

def read_pdb_ends():
    file_path = os.path.join("../data/input/etc", "pdb_ends.txt")
    pdb_ends_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_ends_dict[line[0:4]] = int(line_array[1])
    return pdb_ends_dict

def test_missing_res(filename, file):
    parser = Bio.PDB.PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(id=filename[0:4], file=file)
    old_resid = 1
    for model in structure:
        for ch in model:
            for res in ch:
                resid = res.id[1]
                if not resid == old_resid + 1 and not old_resid == 1:
                    print("Missing res for {}".format(filename[0:4]))
                    stop_resid = resid
                    output_name_1 = file.replace(".pdb", "-a.pdb")  # Create the name of the output file
                    output_1 = open(output_name_1, "w")  # Open the output file
                    output_name_2 = file.replace(".pdb", "-a2.pdb")  # Create the name of the output file
                    output_2 = open(output_name_2, "w")  # Open the output file
                    is_pdb_restarted = False
                    atom_cpt = 0
                    with open(file, "r") as input_file:  # open the pdb file
                        for line in input_file:
                            if line[0:4] == "ATOM":
                                if not int(line[22:26].strip()) >= stop_resid:
                                    output_1.write("%s" % line)
                                else:
                                    atom_cpt += 1
                                    if not is_pdb_restarted:
                                        start_stop_resid_dict[filename[0:4]] = (old_resid, int(line[22:26].strip()))
                                        is_pdb_restarted = True
                                    if atom_cpt >= 100:
                                        line = line.replace(line[0:11], "ATOM    {0}{1}{2}".format(str(atom_cpt)[0],
                                                                                                    str(atom_cpt)[1],
                                                                                                    str(atom_cpt)[2]))
                                    elif atom_cpt >= 10:
                                        line = line.replace(line[0:11], "ATOM     {0}{1}".format(str(atom_cpt)[0],
                                                                                                    str(atom_cpt)[1]))
                                    else:
                                        line = line.replace(line[0:11], "ATOM      {0}".format(str(atom_cpt)[0]))
                                    output_2.write("%s" % line)
                    output_1.write("%s" % 'TER\n')
                    output_1.write("%s" % 'END')
                    output_1.close()
                    output_2.write("%s" % 'TER\n')
                    output_2.write("%s" % 'END')
                    output_2.close()
                    return [output_name_1, output_name_2]
                old_resid = resid


# GLOBAL VARS
start_stop_resid_dict = {}

if __name__ == '__main__':

    import os

    directory_list = list()
    files_list = list()
    paths_list = list()
    scripts_list = list()
    scripts_build_list = list()
    scripts_mini_list = list()
    for root, dirs, files in os.walk("../../compnma_api/data/output/python_output", topdown=False):
        for name in dirs:
            if name == "dcd_SC":
                directory_list.append(os.path.join(root, name))

    pdb_align_list = ['3tfy', '5isv', '4pv6', '2z0z', '1s7l', '2x7b', '3igr', '5k18',
                      '2cns', '5hh0', '5wjd', '5icv', '4kvm', '4u9v']

    for root in directory_list:
        # files_list.extend([x for x in os.listdir(dir) if os.path.isfile(x)])
        for x in os.listdir(root):
            if os.path.isfile(os.path.join(root, x)) and not x.startswith(".")\
                    and x.split("_")[0] in pdb_align_list:
                # Test missing Residues
                split_file_list = test_missing_res(x, os.path.join(root, x))
                if split_file_list:
                    files_list.append(split_file_list[0].split("/")[-1])
                    paths_list.append(split_file_list[0])
                    files_list.append(split_file_list[1].split("/")[-1])
                    paths_list.append(split_file_list[1])
                else:
                    files_list.append(x)
                    paths_list.append(os.path.join(root, x))


    # files_list = [x for x in os.listdir(".") if x.endswith(".pdb") and 'minimized' not in x]
    for file in files_list:
        line_list = list()
        if not file.endswith("-a.pdb") and not file.endswith("-a2.pdb"):
            with open("../data/input/etc/build_mini_nats.inp") as f1:
                for line in f1.readlines():
                    if line.startswith("set protein"):
                        line = "set protein {}\n".format(file.replace(".pdb", ""))
                    if line.lower().startswith("read coor pdb") and read_pdb_starts()[file[0:4]] != 1:
                        # line = "read coor pdb offset -{} unit 1\n".format(str(read_pdb_starts()[file[0:4]] - 1))
                        line = "read coor pdb unit 3 offset -{}\n".format(str(read_pdb_starts()[file[0:4]] - 1))
                    elif line.lower().startswith("read coor pdb"):
                        line = "read coor pdb unit 3\n"
                    if line.lower().startswith("cons fix sele resid 1:"):
                        line = "cons fix sele resid 1:{} .and. -\n".format(str(read_pdb_ends()[file[0:4]] -
                            (read_pdb_starts()[file[0:4]] - 1)))
                    line_list.append(line)
            with open("build_mini_{}.inp".format(file.replace(".pdb", "")), "w") as f2:
                f2.writelines(line_list)
            scripts_list.append("build_mini_{}.inp".format(file.replace(".pdb", "")))
            scripts_build_list.append("build_mini_{}.inp".format(file.replace(".pdb", "")))
        else:
            name = ''
            with open("../data/input/etc/build_mini_nats_residue-gap.inp") as f1:
                for line in f1.readlines():
                    if line.startswith("set protein"):
                        if file.endswith("-a.pdb"):
                            line = "set protein {}\n".format(file.replace("-a.pdb", ""))
                        if file.endswith("-a2.pdb"):
                            line = "set protein {}\n".format(file.replace("-a2.pdb", ""))
                    if line.lower().startswith("read coor pdb unit 10") and read_pdb_starts()[file[0:4]] != 1:
                        line = "read coor pdb unit 10 offset -{}\n".format(str(read_pdb_starts()[file[0:4]] - 1))
                    if line.lower().startswith("read coor pdb unit 10") and read_pdb_starts()[file[0:4]] == 1:
                        line = "read coor pdb unit 10"
                    if line.lower().startswith("read coor pdb unit 11"):
                        line = "read coor pdb unit 11 offset -{}\n".format(str(start_stop_resid_dict[file[0:4]][1]
                                                                               - start_stop_resid_dict[file[0:4]][0]))
                    if line.lower().startswith("cons fix sele resid 1:"):
                        line = "cons fix sele resid 1:{} .and. -\n".format(str(read_pdb_ends()[file[0:4]] -
                            (read_pdb_starts()[file[0:4]] - 1)))
                    line_list.append(line)
            if file.endswith("-a.pdb"):
                name = file.replace("-a.pdb", "")
                with open("build_mini_{}.inp".format(name), "w") as f2:
                    f2.writelines(line_list)
            if file.endswith("-a2.pdb"):
                name = file.replace("-a2.pdb", "")
                with open("build_mini_{}.inp".format(name), "w") as f2:
                    f2.writelines(line_list)
            if "build_mini_{}.inp".format(name) not in scripts_list:
                scripts_list.append("build_mini_{}.inp".format(name))
                scripts_build_list.append("build_mini_{}.inp".format(name))


        # line_list2 = list()
        # if not file.endswith("-a.pdb") and not file.endswith("-a2.pdb"):
        #     with open("../data/input/etc/minimize_nats.inp") as f3:
        #         for line in f3.readlines():
        #             if line.startswith("set protein"):
        #                 line = "set protein {}\n".format(file.replace(".pdb", ""))
        #             if line.lower().startswith("read coor pdb") and read_pdb_starts()[file[0:4]] != 1:
        #                 line = "read coor pdb offset -{} unit 1\n".format(str(read_pdb_starts()[file[0:4]] - 1))
        #             line_list2.append(line)
        #     with open("minimize_{}.inp".format(file.replace(".pdb", "")), "w") as f4:
        #         f4.writelines(line_list2)
        #     scripts_list.append("minimize_{}.inp".format(file.replace(".pdb", "")))
        #     scripts_mini_list.append("minimize_{}.inp".format(file.replace(".pdb", "")))
        # else:
        #     name = ''
        #     with open("../data/input/etc/minimize_nats.inp") as f3:
        #         for line in f3.readlines():
        #             if line.startswith("set protein"):
        #                 if file.endswith("-a.pdb"):
        #                     line = "set protein {}\n".format(file.replace("-a.pdb", ""))
        #                 if file.endswith("-a2.pdb"):
        #                     line = "set protein {}\n".format(file.replace("-a2.pdb", ""))
        #             if line.lower().startswith("read coor pdb") and read_pdb_starts()[file[0:4]] != 1:
        #                 line = "read coor pdb offset -{} unit 1\n".format(str(read_pdb_starts()[file[0:4]] - 1))
        #             line_list2.append(line)
        #     if file.endswith("-a.pdb"):
        #         name = file.replace("-a.pdb", "")
        #         with open("minimize_{}.inp".format(name), "w") as f2:
        #             f2.writelines(line_list2)
        #     if file.endswith("-a2.pdb"):
        #         name = file.replace("-a2.pdb", "")
        #         with open("minimize_{}.inp".format(name), "w") as f2:
        #             f2.writelines(line_list2)
        #     if "minimize_{}.inp".format(name) not in scripts_list:
        #         scripts_list.append("minimize_{}.inp".format(name))
        #         scripts_mini_list.append("minimize_{}.inp".format(name))
        # # elif file.endswith("-a.pdb"):
        # #     with open("../data/input/etc/minimize_nats.inp") as f3:
        # #         for line in f3.readlines():
        # #             if line.startswith("set protein"):
        # #                 line = "set protein {}\n".format(file.replace("-a.pdb", ""))
        # #             if line.lower().startswith("read coor pdb") and read_pdb_starts()[file[0:4]] != 1:
        # #                 line = "read coor pdb offset -{} unit 1\n".format(str(read_pdb_starts()[file[0:4]] - 1))
        # #             line_list2.append(line)
        # #     with open("minimize_{}.inp".format(file.replace("-a.pdb", "")), "w") as f4:
        # #         f4.writelines(line_list2)
        # #     scripts_list.append("minimize_{}.inp".format(file.replace("-a.pdb", "")))
        # #     scripts_mini_list.append("minimize_{}.inp".format(file.replace("-a.pdb", "")))
        # # elif file.endswith("-a2.pdb"):
        # #     with open("../data/input/etc/minimize_nats.inp") as f3:
        # #         for line in f3.readlines():
        # #             if line.startswith("set protein"):
        # #                 line = "set protein {}\n".format(file.replace("-a2.pdb", ""))
        # #             if line.lower().startswith("read coor pdb") and read_pdb_starts()[file[0:4]] != 1:
        # #                 line = "read coor pdb offset -{} unit 1\n".format(str(read_pdb_starts()[file[0:4]] - 1))
        # #             line_list2.append(line)
        # #     with open("minimize_{}.inp".format(file.replace("-a2.pdb", "")), "w") as f4:
        # #         f4.writelines(line_list2)
        # #     scripts_list.append("minimize_{}.inp".format(file.replace("-a2.pdb", "")))
        # #     scripts_mini_list.append("minimize_{}.inp".format(file.replace("-a2.pdb", "")))

    # files_list = [x.replace(".pdb", "") for x in files_list]
    with open("transfer_script_login.sh", "w") as f5:
        str_line = "scp {} pierreb@login.ii.uib.no:~".format(' '.join(paths_list))
        f5.write(str_line)
    with open("transfer_script_doggpil.sh", "w") as f6:
        str_line = "scp {} pierreb@doggpil.cbu.uib.no:/net/orinoco/pierreb/charmm/nathalie_charmm/struct/".format(' '.join(files_list))
        f6.write(str_line)

    # transfer mini pdbs back
    files_mini_list = [x.replace(".pdb", "_mini.pdb").lower() for x in files_list]
    files_mini_list = [x.replace("-a2", "") for x in files_mini_list]
    files_mini_list = [x.replace("-a", "") for x in files_mini_list]
    with open("transfer_struct_mini_doggpil.sh", "w") as f7:
        str_line = "scp {} pierreb@login.ii.uib.no:./\n".format(' '.join(files_mini_list))
        f7.write(str_line)
    with open("transfer_struct_mini_login.sh", "w") as f8:
        str_line = "scp pierreb@login.ii.uib.no:\{{{0}\}} ./".format(','.join(files_mini_list))
        f8.write(str_line)

    with open("transfer_build-mini_login.sh", "w") as f5:
        str_line = "scp {} pierreb@login.ii.uib.no:~".format(' '.join(scripts_list))
        f5.write(str_line)
    with open("transfer_build-mini_doggpil.sh", "w") as f6:
        str_line = "scp {} pierreb@doggpil.cbu.uib.no:/net/orinoco/pierreb/charmm/nathalie_charmm/calc/".format(' '.join(scripts_list))
        f6.write(str_line)

    # run build script
    with open("run_build.sh", "w") as f5:
        f5.write("cd ../struct/\n")
        f5.write("for f in `find`; do mv -v \"$f\" \"`echo $f | tr '[A-Z]' '[a-z]'`\"; done\n")
        f5.write("cd ../calc/\n")
        for script in scripts_build_list:
            str_line = "/net/orinoco/apps/charmm/c38b2/exec/gnu_xxlarge/charmm_xxlarge < {0} > {1} | ".format(script,
                                                                                                           script.replace(".inp", ".out"))
            f5.write(str_line)
    with open("run_build.sh", 'rb+') as filehandle:
        filehandle.seek(-1, os.SEEK_END)
        filehandle.truncate()

    ## run minimize script
    # with open("run_minimize.sh", "w") as f5:
    #     for script in scripts_mini_list:
    #         str_line = "/net/orinoco/apps/charmm/c38b2/exec/gnu_xxlarge/charmm_xxlarge < {0} > {1} | ".format(script,
    #                                                                                                        script.replace(".inp", ".out"))
    #         f5.write(str_line)
    # with open("run_minimize.sh", 'rb+') as filehandle:
    #     filehandle.seek(-1, os.SEEK_END)
    #     filehandle.truncate()
