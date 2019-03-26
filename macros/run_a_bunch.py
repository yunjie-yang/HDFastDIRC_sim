import subprocess, shutil, os, csv

def read_control_in(infile):
    dict_to_return = {}
    with open(infile,'r') as f:
        reader = csv.reader(f, delimiter=' ')
        for line in reader:
            if len(line) >= 2 and line[0][0]!='c' and line[0][0]!='C':
                dict_to_return.update({line[0].strip():line[1].strip()})
    return dict_to_return


def main():
    output_dir_top = "/media/sf_SharedFolderVM/FastDIRC_geometry/outputs"
    exec_dir = "/home/yunjiey/Documents/HDFastDIRC_sim"

    for bar_i in range(12):
        print "Processing Bar #%d" % (bar_i)
        # check/create output dir
        output_dir_bar = "%s/bar_%d"%(output_dir_top,bar_i)
        if not os.path.exists(output_dir_bar):
            os.mkdir(output_dir_bar)
        # write control.in and copy over 
        control_in_file = "%s/control.in"%output_dir_bar
        outfile = open(control_in_file,'w')
        lines = []

	##############################################################
	#		Write desired control.in file here	     #	
	##############################################################

        lines.append("SIM_ONLY  1\n")
        lines.append("OUTFILE  '%s/hist_bar_%d.root'\n"%(output_dir_bar,bar_i))
        lines.append("GEOMETRY_INFILE  '/media/sf_SharedFolderVM/FastDIRC_geometry/FastDIRC_HDDS_Nominal.csv'\n")
        lines.append("GEOMETRY_OUTFILE  '%s/dirc_model_geometry_bar_%d.csv'\n"%(output_dir_bar,bar_i))
        lines.append("PARTICLE_BAR  %.01f\n"%bar_i)
	lines.append("END \n")

	##############################################################

        outfile.writelines(lines)
        outfile.close()
        shutil.copyfile(control_in_file,"%s/control.in"%exec_dir)
        # move to directory and execute
        os.chdir(exec_dir)
        command = "./hdfastdirc"
        subprocess.call(command,shell=True)

if __name__=="__main__":
    main()
