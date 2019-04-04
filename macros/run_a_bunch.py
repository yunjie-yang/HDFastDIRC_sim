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

    #bar_theta_phi = [[1,4,32],[3,2,45],[3,3,165],[5,1,-75],[6,7,-110],[13,5,141]]
    #bar_theta_phi_x_y = [[14,4,40,-250,-68.8493],[2,2.5,68.4,120.8,-32.9],[4,1.2,113.5,180.2,-25.1]]
    #bar_theta_phi_x_y = [[1,4,32,12.47,-13.8358]]
    #bar_theta_phi_x_y = [[6,2.5,68.4,120.8,-32.9]]
    bar_theta_phi_x_y = [[3,1.2,113.5,180.2,-25.1]]

    #for bar_i in range(12):
    #for item in bar_theta_phi:
    for item in bar_theta_phi_x_y:
	bar_i = int(item[0])
	if len(item)>1:
	    theta_i = float(item[1])
	    phi_i = float(item[2])
	if len(item)>3:
	    x_i = float(item[3])
	    y_i = float(item[4])

        print "Processing Bar #%d" % (bar_i)
        # check/create output dir

	if len(item)==3:
	    output_filelabel = "bar_%d_%d_%d"%(bar_i,theta_i,phi_i)
	if len(item)>3:
            output_filelabel = "bar_%d_%d_%d_%d_%d"%(bar_i,int(round(theta_i)),int(round(phi_i)),int(round(x_i)),int(round(y_i)))

        output_dir_bar = "%s/%s"%(output_dir_top,output_filelabel)

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
        lines.append("OUTFILE  '%s/hist.root'\n"%(output_dir_bar))
        lines.append("GEOMETRY_INFILE  '/media/sf_SharedFolderVM/FastDIRC_geometry/FastDIRC_HDDS_Nominal.csv'\n")
        lines.append("GEOMETRY_OUTFILE  '%s/dirc_model_geometry_bar_%d.csv'\n"%(output_dir_bar,bar_i))
        lines.append("PARTICLE_BAR  %.01f\n"%bar_i)
 	if len(item)>1:
	    lines.append("PARTICLE_THETA  %.01f\n"%theta_i)
	    lines.append("PARTICLE_PHI  %.01f\n"%(phi_i))
	if len(item)>3:
            lines.append("PARTICLE_X  %.05f\n"%(x_i))
            lines.append("PARTICLE_Y  %.05f\n"%(y_i))
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
