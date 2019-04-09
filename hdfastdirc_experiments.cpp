#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "src/dirc_point.h"
#include "src/dirc_threesegbox_sim.h"
#include "src/dirc_rect_digitizer.h"
#include "src/dirc_spread_gaussian.h"
#include "src/GlueXUserOptions.h"

#include <TFile.h>
#include <TTree.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMinuit.h>

#include <TChain.h>
#include <TClonesArray.h>
#include "src/DrcHit.h"
#include "src/DrcEvent.h"

#define rad2deg 57.2958

int main(int nargs, char* argv[])
{
	/**************************************************************************/
	/***********               INITIALIZATION              ********************/
	/**************************************************************************/

	bool SIM_ONLY = 1;

	const char* config_str;


	double energy = 5.;
	double kmass = .4937;
	double pimass = .1396;

	double particle_x = 0;
	double particle_y = 0;
	double particle_theta = 0;
	double particle_phi = 0;
	double particle_bar = 2;


	bool kaleidoscope_plot = false;
	bool use_moliere_scattering = false;

	double wedge_uncertainty = 0/57.3;
	double mirror_angle_change = 0;
	double mirror_angle_change_unc = 0;
	double mirror_angle_change_yunc = 0;
	double box_rot = 0;
	double bar_box_box_angle = 0/57.3;
	double mirror_r_difference = 400;//1200 - 400 = 800.  Changed 05/09/2016.  Does not affect threeseg mirror reconstruction as far as I can tell - this was known.
	double wedge_non_uniformity = 0;
	double pmt_offset = 0;
	double pmt_angle_offset = 0;
	double main_mirror_nonuniformity = 0;
	double foc_mirror_size = 288;

	double main_mirror_angle_off = 0;
	double main_mirror_yangle_off = 0;
	double main_mirror_zangle_off = 0;
	double main_mirror_yoff = 0;
	double main_mirror_zoff = 0;

	double bar_box_xoff = 0;
	double bar_box_yoff = 0;
	double bar_box_zoff = 0;

	int rseed = 1337;

	//double tracking_unc = .0000*57.3; //mrad
	double ckov_unc = .003*57.3; //transport = 3mrad

	double resx = 6;
	double resy = 6;
	double rest = 1;
	double minx = -1500;
	double maxx = 1500;
	double miny = -500;
	double maxy = 500;
	double mint = 0;
	double maxt = 1000;
	double t_unc = .27;
	double t_bin_size = 1;

	double digit_miny = -50;
	double digit_maxy = 300;

	digit_miny = miny;
	digit_maxy = maxy;


	int n_sim_phots = 40;

	//int n_phi_phots = 150000;
	//int n_phi_phots = 15000;
	int n_phi_phots = 5000;
	int n_z_phots = 4;

	bool use_quartz_for_liquid = false;

	double foc_mirror_yoff = 0;
	double foc_mirror_zoff = 0;


	double liquid_absorbtion = 0*-log(.7)/1000;
	double liquid_index = 1.33;

        double s_func_x = 6;
        double s_func_y = s_func_x;
        double s_func_t = 1.0;
        double sfunc_sig = 1;

	char* geometry_infilename = new char[256];
	sprintf(geometry_infilename, "FastDIRC_geometry_input.csv");

	char* geometry_outfilename = new char[256];
	sprintf(geometry_outfilename,"dirc_model_geometry.csv");

	char* root_outfilename = new char[256];
	sprintf(root_outfilename,"fitdirc.root");

	char* dirctree_filename = new char[256];
	sprintf(dirctree_filename,"dirc_tree.root");

	int Nrandom_points = 7;

	double res_enhance = 1;
	TH1F *pion_dist_x = new TH1F("pion_dist_x","x val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *pion_dist_y = new TH1F("pion_dist_y","y val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *pion_dist_t = new TH1F("pion_dist_t","t val of intercepted points - pion",(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *pion_dist_xt = new TH2F("pion_dist_xt","xt val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *pion_dist_yt = new TH2F("pion_dist_yt","yt val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *pion_dist_rowcol = new TH2F("pion_dist_rowcol","hit pattern - pion; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);
	TH3F *pion_dist_3D = new TH3F("pion_dist_3D","(x,y,t) - pion; Pixel Row ; Pixel Column ; Hit Time (ns)",144,-0.5,143.5,48,-0.5,47.5,300,0,300);


	TH1F *kaon_dist_x = new TH1F("kaon_dist_x","x val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *kaon_dist_y = new TH1F("kaon_dist_y","y val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *kaon_dist_t = new TH1F("kaon_dist_t","t val of intercepted points - kaon",(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_xy = new TH2F("kaon_dist_xy","xy val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *kaon_dist_xt = new TH2F("kaon_dist_xt","xt val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_yt = new TH2F("kaon_dist_yt","yt val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_rowcol = new TH2F("kaon_dist_rowcol","hit pattern - kaon; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);
	TH3F *kaon_dist_3D = new TH3F("kaon_dist_3D","(x,y,t) - kaon; Pixel Row ; Pixel Column ; Hit Time (ns)",144,-0.5,143.5,48,-0.5,47.5,300,0,300);


        TH1F *ll_diff_pion = new TH1F("ll_diff_pion","Difference of log likelihood real = pion",200000,-200,200);
        TH1F *ll_diff_kaon = new TH1F("ll_diff_kaon","Difference of log likelihood real = kaon",200000,-200,200);
        TH1F *phot_found_pion = new TH1F("phot_found_pion","number of photons found on pion angle", 1001,-.5,1000.5);
        TH1F *phot_found_kaon = new TH1F("phot_found_kaon","number of photons found on kaon angle", 1001,-.5,1000.5);


	TH2F *hit_dist_rowcol_tree = new TH2F("hit_dist_rowcol_tree","hit pattern; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);
	TH1F *hit_dist_t_tree = new TH1F("hit_dist_t_tree","hit time ; hit time (ns);",(maxt-mint)/(res_enhance*rest),mint,maxt);

        TH1F *ll_diff_tree = new TH1F("ll_diff_tree","Difference of log likelihood",200000,-200,200);

        TH1F *ll_diff_pion_tree = new TH1F("ll_diff_pion_tree","Difference of log likelihood (pion selection)",200000,-200,200);
        TH1F *ll_diff_kaon_tree = new TH1F("ll_diff_kaon_tree","Difference of log likelihood (kaon selection)",200000,-200,200);

        TH1F *ll_diff_pion_inexact = new TH1F("ll_diff_pion_inexact","Difference of log likelihood (inexact) real = pion",200000,-200,200);
        TH1F *ll_diff_kaon_inexact = new TH1F("ll_diff_kaon_inexact","Difference of log likelihood (inexact) real = kaon",200000,-200,200);

        TH1F *Nph_tree = new TH1F("Nph_tree","number of photons", 1001,-.5,1000.5);
        TH1F *Nph_pion = new TH1F("Nph_pion","number of photons (real = pion)", 1001,-.5,1000.5);
        TH1F *Nph_kaon = new TH1F("Nph_kaon","number of photons (real = kaon)", 1001,-.5,1000.5);


	std::vector<std::vector<double>> PBins =  {
								{0.,1.5},
								{1.5, 2.5},
								{2.5, 3.5},
								{3.5, 4.5},
								{4.5, 12.}
							 };

	int NumBinsLL = 200000;
	double MinLL  = -200.;
	double MaxLL  = 200.;

	std::map<int,TH1F*> histMap_ll_diff_pion_tree;
	std::map<int,TH1F*> histMap_ll_diff_kaon_tree;
	std::map<int,TH1F*> histMap_ll_diff_pion_sim_exact;
	std::map<int,TH1F*> histMap_ll_diff_kaon_sim_exact;
	std::map<int,TH1F*> histMap_ll_diff_pion_sim_inexact;
	std::map<int,TH1F*> histMap_ll_diff_kaon_sim_inexact;

	int NumBinsNhits = 1000;
	double MinNhits  = -0.5;
	double MaxNhits  = 999.5;
	std::map<int,TH1F*> histMap_Nhits_pion_tree;
	std::map<int,TH1F*> histMap_Nhits_kaon_tree;
	std::map<int,TH1F*> histMap_Nhits_pion_sim_inexact;
	std::map<int,TH1F*> histMap_Nhits_kaon_sim_inexact;
	for (size_t locPBin = 0 ; locPBin < PBins.size() ; locPBin++)
	{
		histMap_ll_diff_pion_tree[int(locPBin)] = new TH1F(Form("ll_diff_pion_tree_PBin%d",int(locPBin)),\
								   Form("DLL (#pi selection), P: [%.1f,%.1f];DLL;",\
								        PBins[locPBin][0],PBins[locPBin][1]),\
								   NumBinsLL,MinLL,MaxLL);
		histMap_ll_diff_kaon_tree[int(locPBin)] = new TH1F(Form("ll_diff_kaon_tree_PBin%d",int(locPBin)),\
								   Form("DLL (K selection), P: [%.1f,%.1f];DLL;",\
								        PBins[locPBin][0],PBins[locPBin][1]),\
								   NumBinsLL,MinLL,MaxLL);
		histMap_ll_diff_pion_sim_exact[int(locPBin)] = new TH1F(Form("ll_diff_pion_sim_exact_PBin%d",int(locPBin)),\
								        Form("DLL (#pi sim), P: [%.1f,%.1f];DLL;",\
								             PBins[locPBin][0],PBins[locPBin][1]),\
								        NumBinsLL,MinLL,MaxLL);
		histMap_ll_diff_kaon_sim_exact[int(locPBin)] = new TH1F(Form("ll_diff_kaon_sim_exact_PBin%d",int(locPBin)),\
								        Form("DLL (K sim), P: [%.1f,%.1f];DLL;",\
								             PBins[locPBin][0],PBins[locPBin][1]),\
								        NumBinsLL,MinLL,MaxLL);
		histMap_ll_diff_pion_sim_inexact[int(locPBin)] = new TH1F(Form("ll_diff_pion_sim_inexact_PBin%d",int(locPBin)),\
									  Form("DLL (#pi sim), P: [%.1f,%.1f];DLL;",\
									       PBins[locPBin][0],PBins[locPBin][1]),\
									  NumBinsLL,MinLL,MaxLL);
		histMap_ll_diff_kaon_sim_inexact[int(locPBin)] = new TH1F(Form("ll_diff_kaon_sim_inexact_PBin%d",int(locPBin)),\
									  Form("DLL (K sim), P: [%.1f,%.1f];DLL;",\
									       PBins[locPBin][0],PBins[locPBin][1]),\
									  NumBinsLL,MinLL,MaxLL);

		histMap_Nhits_pion_tree[int(locPBin)] = new TH1F(Form("Nhits_pion_tree_PBin%d",int(locPBin)),\
								   Form("Nhits (#pi selection), P: [%.1f,%.1f];DLL;",\
								        PBins[locPBin][0],PBins[locPBin][1]),\
								   NumBinsNhits,MinNhits,MaxNhits);
		histMap_Nhits_kaon_tree[int(locPBin)] = new TH1F(Form("Nhits_kaon_tree_PBin%d",int(locPBin)),\
								   Form("Nhits (K selection), P: [%.1f,%.1f];DLL;",\
								        PBins[locPBin][0],PBins[locPBin][1]),\
								   NumBinsNhits,MinNhits,MaxNhits);
		histMap_Nhits_pion_sim_inexact[int(locPBin)] = new TH1F(Form("Nhits_pion_sim_inexact_PBin%d",int(locPBin)),\
									Form("Nhits (#pi sim), P: [%.1f,%.1f];DLL;",\
									     PBins[locPBin][0],PBins[locPBin][1]),\
									NumBinsNhits,MinNhits,MaxNhits);
		histMap_Nhits_kaon_sim_inexact[int(locPBin)] = new TH1F(Form("Nhits_kaon_sim_inexact_PBin%d",int(locPBin)),\
									Form("Nhits (K sim), P: [%.1f,%.1f];DLL;",\
									     PBins[locPBin][0],PBins[locPBin][1]),\
									NumBinsNhits,MinNhits,MaxNhits);
	}






	double pion_beta, kaon_beta;
	pion_beta=kaon_beta=-1;
	//double pion_angle, kaon_angle;
	//pion_angle=kaon_angle = -1;

	/**************************************************************************/
	/***********               READ CONFIG                 ********************/
	/**************************************************************************/

	if (nargs==2) config_str = argv[1];
	else config_str = "control.in";

	GlueXUserOptions user_opts;
	if (user_opts.ReadControl_in(config_str) == 0)
	{
		std::cerr << "Reading control.in failed" << std::endl;
		exit(-1);
	}

	std::map<int, int> opt_int;
	std::map<int, std::string> opt_str;
	std::map<int, double> opt_val;

	if (user_opts.Find("SIM_ONLY", opt_int)) SIM_ONLY = opt_int[1];
	if (user_opts.Find("OUTFILE", opt_str)) sprintf(root_outfilename,"%s",opt_str[1].c_str());
	if (user_opts.Find("DIRCTREE_INFILE", opt_str)) sprintf(dirctree_filename,"%s",opt_str[1].c_str());
	if (user_opts.Find("GEOMETRY_INFILE", opt_str)) sprintf(geometry_infilename,"%s",opt_str[1].c_str());
	if (user_opts.Find("GEOMETRY_OUTFILE", opt_str)) sprintf(geometry_outfilename,"%s",opt_str[1].c_str());

	if (user_opts.Find("N_PHI_PHOTS", opt_int)) n_phi_phots = opt_int[1];

	printf("SIM_ONLY         = %d \n",SIM_ONLY);
	printf("OUTFILE          = %s \n", root_outfilename);
	printf("DIRCTREE_INFILE  = %s \n", dirctree_filename);
	printf("GEOMETRY_INFILE  = %s \n", geometry_infilename);


	/**************************************************************************/
	/***********              APPLY CONFIG                 ********************/
	/**************************************************************************/

	TFile* tfile = new TFile(root_outfilename,"RECREATE");

	/**************************************************************************/
	/***********           INITIALIZE DIRC MODEL           ********************/
	/**************************************************************************/

	DircRectDigitizer digitizer(\
			minx,\
			maxx,\
			resx,\
			digit_miny,\
			digit_maxy,\
			resy,\
			t_unc,\
			t_bin_size);


	double main_mirror_angle_nominal = 74.0197;
	double mirror_r_nominal  = -1200 + mirror_r_difference;
	double pmt_angle_nominal = 47.87 + box_rot + mirror_angle_change;

	// suppose this is the truth/reality
	DircThreeSegBoxSim *dirc_model_truth = new DircThreeSegBoxSim(\
			rseed,\
			mirror_r_nominal,\
			foc_mirror_size,\
			main_mirror_angle_nominal,\
			600,\
			pmt_angle_nominal,\
			geometry_infilename);


	dirc_model_truth->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	dirc_model_truth->set_wedge_mirror_rand(wedge_non_uniformity);


	//various running condition controls
	dirc_model_truth->set_store_traveled(false);// uses LOTS of memory if set to true.
	dirc_model_truth->set_liquid_index(liquid_index);
	dirc_model_truth->set_kaleidoscope_plot(kaleidoscope_plot);
	dirc_model_truth->set_use_quartz_n_for_liquid(use_quartz_for_liquid);
	dirc_model_truth->set_use_moliere(use_moliere_scattering);
	dirc_model_truth->set_moliere_p(energy*1000);//assume momentum is the same for both for now - high energy;
	dirc_model_truth->set_liquid_absorbtion(liquid_absorbtion);
	dirc_model_truth->set_liquid_index(liquid_index);


	// suppose this is what we thought/used as the truth
	DircThreeSegBoxSim *dirc_model_used = new DircThreeSegBoxSim(\
			rseed,\
			mirror_r_nominal,\
			foc_mirror_size,\
			main_mirror_angle_nominal,\
			600,\
			pmt_angle_nominal,\
			geometry_infilename);


	dirc_model_used->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	dirc_model_used->set_wedge_mirror_rand(wedge_non_uniformity);


	//various running condition controls
	dirc_model_used->set_store_traveled(false);// uses LOTS of memory if set to true.
	dirc_model_used->set_liquid_index(liquid_index);
	dirc_model_used->set_kaleidoscope_plot(kaleidoscope_plot);
	dirc_model_used->set_use_quartz_n_for_liquid(use_quartz_for_liquid);
	dirc_model_used->set_use_moliere(use_moliere_scattering);
	dirc_model_used->set_moliere_p(energy*1000);//assume momentum is the same for both for now - high energy;
	dirc_model_used->set_liquid_absorbtion(liquid_absorbtion);
	dirc_model_used->set_liquid_index(liquid_index);

	/**************************************************************************/
	/***********       APPY OFFSETS TO NOMINAL GEOMETRY       *****************/
	/**************************************************************************/

	TRandom3 spread_ang(rseed+3);

/*
	dirc_model_truth->set_focus_mirror_angle(			\
					   spread_ang.Gaus(main_mirror_angle_nominal,mirror_angle_change_unc)+main_mirror_angle_off, \
					   3., \
					   main_mirror_zangle_off);
*/

/*
	dirc_model_truth->set_pmt_offset(pmt_offset);
	dirc_model_truth->set_upper_wedge_angle_diff(wedge_uncertainty);
	dirc_model_truth->set_bar_box_angle(bar_box_box_angle);

	dirc_model_truth->set_mirror_plane_offsets(foc_mirror_yoff,foc_mirror_zoff);


	dirc_model_truth->set_pmt_offset(pmt_offset);
	dirc_model_truth->set_pmt_angle(47.87+pmt_angle_offset);

	dirc_model_truth->set_focus_mirror_angle(			\
					   spread_ang.Gaus(main_mirror_angle_nominal,mirror_angle_change_unc)+main_mirror_angle_off, \
					   spread_ang.Gaus(0,mirror_angle_change_yunc)+main_mirror_yangle_off, \
					   main_mirror_zangle_off);

	dirc_model_truth->set_mirror_plane_offsets(			\
						 main_mirror_yoff,	\
						 main_mirror_zoff);

	dirc_model_truth->set_bar_box_offsets(\
						 bar_box_xoff,\
						 bar_box_yoff,\
						 bar_box_zoff);

	//dirc_model_truth->set_geometry_outfile(geometry_outfilename);
	//dirc_model_truth->print_model();
*/

	// Construct the "used" model by applying some offsets
	dirc_model_used->set_pmt_offset(pmt_offset);
	dirc_model_used->set_upper_wedge_angle_diff(wedge_uncertainty);
	dirc_model_used->set_bar_box_angle(bar_box_box_angle);

	dirc_model_used->set_mirror_plane_offsets(foc_mirror_yoff,foc_mirror_zoff);


	dirc_model_used->set_pmt_offset(pmt_offset);
	dirc_model_used->set_pmt_angle(47.87+pmt_angle_offset);

	dirc_model_used->set_focus_mirror_angle(			\
					   spread_ang.Gaus(main_mirror_angle_nominal,mirror_angle_change_unc)+main_mirror_angle_off, \
					   spread_ang.Gaus(0,mirror_angle_change_yunc)+main_mirror_yangle_off, \
					   main_mirror_zangle_off);

	dirc_model_used->set_mirror_plane_offsets(			\
						 main_mirror_yoff,	\
						 main_mirror_zoff);

	dirc_model_used->set_bar_box_offsets(\
						 bar_box_xoff,\
						 bar_box_yoff,\
						 bar_box_zoff);

	dirc_model_used->set_geometry_outfile(geometry_outfilename);
	dirc_model_used->print_model();


	printf("\n DIRC model all set. \n");

	/**************************************************************************/
	/***********           INITIALIZE INPUT TREE           ********************/
	/**************************************************************************/

	TChain* event_chain = new TChain("dirc");

	event_chain -> Add(dirctree_filename);

	DrcEvent* evt = new DrcEvent();

	TClonesArray* particle_array = new TClonesArray("DrcEvent");

	event_chain->SetBranchAddress("DrcEvent", &particle_array);

	double particle_mom_tree   = 0.;
	double particle_theta_tree = 0.;
	double particle_phi_tree   = 0.;
	double particle_x_tree     = 0.;
	double particle_y_tree     = 0.;
	double particle_t_tree     = 0.;
	int    particle_pid_tree   = 0.;
	
	int hit_ChannelId   = -1;
	int hit_pixelrow = -1;
	int hit_pixelcol = -1;
	double hit_t   = -1;

	TVector3 particle_p3_tree, particle_x3_tree;

	double loc_energy_pion(0.);
	double loc_energy_kaon(0.);

	std::vector<dirc_point> sim_points_pion;
	std::vector<dirc_point> sim_points_kaon;

	std::vector<dirc_point> sim_points_pion_inexact;
	std::vector<dirc_point> sim_points_kaon_inexact;

	std::vector<dirc_point> hit_points_pion;
	std::vector<dirc_point> hit_points_kaon;

	std::vector<dirc_point> tree_points;

	TRandom3 rndm(rseed+3);

	double ll_pion(0.),ll_kaon(0.);
	double ll_pion_tree(0.),ll_kaon_tree(0.);

	int counter = 0;

	int locBin_PID = -1;
	int locBin_P   = -1;

	bool locGoodParticle = false;
	
	printf("\n\n Total Entries = %d\n",int(event_chain->GetEntries()));

	//for (int event_i = 0 ; event_i < 15001 ; event_i++)
	for (int event_i = 0 ; event_i < event_chain->GetEntries() ; event_i++)
	{
	event_chain->GetEntry(event_i);
	if(event_i%1000==0)
		printf("Event #%d Counter = %d\n",event_i,counter);

	for (int particle_i = 0; particle_i < particle_array->GetEntriesFast(); particle_i++)
	{
		
		evt = (DrcEvent*) particle_array->At(particle_i);
		particle_p3_tree    = evt->GetMomentum();
		particle_x3_tree    = evt->GetPosition();
		particle_t_tree     = evt->GetTime();
		particle_pid_tree   = evt->GetPdg();

		particle_mom_tree   = particle_p3_tree.Mag();
		particle_theta_tree = particle_p3_tree.Theta()*rad2deg; 
		particle_phi_tree   = particle_p3_tree.Phi()*rad2deg; 
		particle_x_tree     = particle_x3_tree.X();
		particle_y_tree     = particle_x3_tree.Y(); 

		for (size_t locPBin = 0; locPBin < PBins.size(); locPBin++)
		{
			if (particle_mom_tree > PBins[locPBin][0] && particle_mom_tree < PBins[locPBin][1])
				locBin_P = locPBin; 
		}
		if (locBin_P == -1)
			continue;


		//if (particle_mom_tree < 2.8 || particle_mom_tree > 3.2) continue;
		//if (particle_mom_tree < 3.5 || particle_mom_tree > 4.5) continue;
		//if (particle_mom_tree < 4) continue;
		//if (particle_mom_tree > 3) continue;

		//if (particle_x_tree > 5. || particle_x_tree < 0.) continue;
		//if (particle_y_tree > -19.1158 || particle_y_tree < -22.6158) continue;//bar 3

		if (abs(particle_pid_tree)==211)
			locBin_PID = 0;
		else if (abs(particle_pid_tree)==321)
			locBin_PID = 1;
		else
			continue;

		locGoodParticle = dirc_model_truth->convert_particle_kinematics(particle_x,particle_y,particle_theta,particle_phi,particle_bar,particle_x_tree,particle_y_tree,particle_theta_tree,particle_phi_tree);

		if (!locGoodParticle)
			continue;

		counter++;

		for (auto hit : evt->GetHits())
		{
			dirc_point tree_point;
			hit_ChannelId   = hit.GetChannel();
			hit_t           = hit.GetLeadTime() - particle_t_tree;	
			hit_pixelrow    = digitizer.GetPixelRow(hit_ChannelId);
			hit_pixelcol    = digitizer.GetPixelColumn(hit_ChannelId);

			hit_dist_rowcol_tree -> Fill(hit_pixelrow,hit_pixelcol); 
			hit_dist_t_tree -> Fill(hit_t); 
			tree_point.pixel_row = hit_pixelrow;
			tree_point.pixel_col = hit_pixelcol;
			tree_point.t         = hit_t;

			digitizer.undigitize_point(tree_point);

			tree_points.push_back(tree_point);

		}

		Nph_tree -> Fill(int(evt->GetHitSize()));
	
		loc_energy_pion = sqrt(particle_mom_tree*particle_mom_tree + pimass*pimass);
		loc_energy_kaon = sqrt(particle_mom_tree*particle_mom_tree + kmass*kmass);

		pion_beta = dirc_model_truth->get_beta(loc_energy_pion,pimass);
		kaon_beta = dirc_model_truth->get_beta(loc_energy_kaon,kmass);


		/*********************************************************************************************/
		/************   Generate FastDIRC-predicted Hits Given Particle Kinematics  ******************/
		/*********************************************************************************************/

		if (locBin_PID == 0)
		{
			dirc_model_truth->sim_rand_n_photons(\
					sim_points_pion_inexact,\
					n_sim_phots,\
					-1,\
					particle_bar,\
					particle_x,\
					particle_y,\
					0.,\
					particle_theta,\
					particle_phi,\
					0.,\
					ckov_unc/1.,\
					pion_beta);
			digitizer.digitize_points(sim_points_pion_inexact);

			dirc_model_used->sim_rand_n_photons(\
					sim_points_pion,\
					n_sim_phots,\
					-1,\
					particle_bar,\
					particle_x,\
					particle_y,\
					0.,\
					particle_theta,\
					particle_phi,\
					0.,\
					ckov_unc/1.,\
					pion_beta);
			digitizer.digitize_points(sim_points_pion);



		}

		if (locBin_PID == 1)
		{
			dirc_model_truth->sim_rand_n_photons(\
					sim_points_kaon_inexact,\
					n_sim_phots,\
					-1,\
					particle_bar,\
					particle_x,\
					particle_y,\
					0.,\
					particle_theta,\
					particle_phi,\
					0.,\
					ckov_unc/1.,\
					kaon_beta);
			digitizer.digitize_points(sim_points_kaon_inexact);

			dirc_model_used->sim_rand_n_photons(\
					sim_points_kaon,\
					n_sim_phots,\
					-1,\
					particle_bar,\
					particle_x,\
					particle_y,\
					0.,\
					particle_theta,\
					particle_phi,\
					0.,\
					ckov_unc/1.,\
					kaon_beta);
			digitizer.digitize_points(sim_points_kaon);

		}

		// Add random noise
		for (int loc_i = 0 ; loc_i < Nrandom_points ; loc_i++)
		{
			dirc_point random_point;
			digitizer.get_random_point(random_point);
			if (locBin_PID == 0)
				sim_points_pion_inexact.push_back(random_point);
			if (locBin_PID == 1)
				sim_points_kaon_inexact.push_back(random_point);

		}

        	double x,y,t_ns;
        	int pixel_row;
        	for (unsigned int i = 0; i < sim_points_pion_inexact.size(); i++)
        	{
                	x = sim_points_pion_inexact[i].x;
               		y = sim_points_pion_inexact[i].y;
                	t_ns = sim_points_pion_inexact[i].t;
                	pixel_row = sim_points_pion_inexact[i].pixel_row;
                	if (pixel_row<0) continue;
                	pion_dist_x->Fill(x);
                	pion_dist_y->Fill(y);
                	pion_dist_t->Fill(t_ns);
                	pion_dist_xy->Fill(x,y);
                	pion_dist_xt->Fill(x,t_ns);
                	pion_dist_yt->Fill(y,t_ns);
                	pion_dist_t->Fill(t_ns);

                	pion_dist_rowcol->Fill(sim_points_pion_inexact[i].pixel_row,sim_points_pion_inexact[i].pixel_col);
                	pion_dist_3D->Fill(sim_points_pion_inexact[i].pixel_row,sim_points_pion_inexact[i].pixel_col,t_ns);
        	}
		for (unsigned int i = 0; i < sim_points_kaon_inexact.size(); i++)
		{
			x = sim_points_kaon_inexact[i].x;
			y = sim_points_kaon_inexact[i].y;
			t_ns = sim_points_kaon_inexact[i].t;
			pixel_row = sim_points_kaon_inexact[i].pixel_row;
			if (pixel_row<0) continue;
			kaon_dist_x->Fill(x);
			kaon_dist_y->Fill(y);
			kaon_dist_t->Fill(t_ns);
			kaon_dist_xy->Fill(x,y);
			kaon_dist_xt->Fill(x,t_ns);
			kaon_dist_yt->Fill(y,t_ns);
			kaon_dist_t->Fill(t_ns);

			kaon_dist_rowcol->Fill(sim_points_kaon_inexact[i].pixel_row,sim_points_kaon_inexact[i].pixel_col);
                	kaon_dist_3D->Fill(sim_points_kaon_inexact[i].pixel_row,sim_points_kaon_inexact[i].pixel_col,t_ns);
		}

		/*********************************************************************************/
		/************      Generate PDFs Given Particle Kinematics      ******************/
		/*********************************************************************************/

	        dirc_model_used->sim_reg_n_photons(\
        	                hit_points_pion,\
                	        n_phi_phots,\
                        	n_z_phots,\
                        	-1,\
                 		particle_bar,\
                        	particle_x,\
                        	particle_y,\
                        	0.,\
                        	particle_theta,\
                        	particle_phi,\
                        	0,\
                        	ckov_unc/1.,\
                        	pion_beta,\
                        	-1);
 
		dirc_model_used->sim_reg_n_photons(\
				hit_points_kaon,\
				n_phi_phots,\
				n_z_phots,\
				-1,\
				particle_bar,\
				particle_x,\
				particle_y,\
				0.,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/1.,\
				kaon_beta,\
				-1); 

        	DircSpreadGaussian* pdf_pion = new DircSpreadGaussian(\
                	        sfunc_sig,\
                        	hit_points_pion,\
                        	s_func_x,\
                        	s_func_y,\
                        	s_func_t);

        	DircSpreadGaussian* pdf_kaon = new DircSpreadGaussian(\
                	        sfunc_sig,\
                        	hit_points_kaon,\
                        	s_func_x,\
                        	s_func_y,\
                        	s_func_t);

		/***************************************************************************************************/
		/************      Evaluate Likelihoods for both observed and generated hits      ******************/
		/***************************************************************************************************/

		//observed hits
		ll_pion_tree = pdf_pion->get_log_likelihood(tree_points);
		ll_kaon_tree = pdf_kaon->get_log_likelihood(tree_points);

		ll_diff_tree->Fill(1*(ll_pion_tree-ll_kaon_tree));

		if (locBin_PID == 0)
		{
			ll_diff_pion_tree->Fill(1*(ll_pion_tree-ll_kaon_tree));
			histMap_ll_diff_pion_tree[locBin_P]->Fill(1*(ll_pion_tree-ll_kaon_tree));
			histMap_Nhits_pion_tree[locBin_P]->Fill(int(tree_points.size()));
		}
		if (locBin_PID == 1)
		{
			ll_diff_kaon_tree->Fill(1*(ll_pion_tree-ll_kaon_tree));
			histMap_ll_diff_kaon_tree[locBin_P]->Fill(1*(ll_pion_tree-ll_kaon_tree));
			histMap_Nhits_kaon_tree[locBin_P]->Fill(int(tree_points.size()));
		}
		//generated hits
		// if we know the PDF exactly
		ll_pion = pdf_pion->get_log_likelihood(sim_points_pion);
		ll_kaon = pdf_kaon->get_log_likelihood(sim_points_pion);
		if (locBin_PID == 0)
		{
			ll_diff_pion->Fill(1*(ll_pion-ll_kaon));
			histMap_ll_diff_pion_sim_exact[locBin_P]->Fill(1*(ll_pion-ll_kaon));

		}
		ll_pion = pdf_pion->get_log_likelihood(sim_points_kaon);
		ll_kaon = pdf_kaon->get_log_likelihood(sim_points_kaon);
		if (locBin_PID == 1)
		{
			ll_diff_kaon->Fill(1*(ll_pion-ll_kaon));
			histMap_ll_diff_kaon_sim_exact[locBin_P]->Fill(1*(ll_pion-ll_kaon));
		}
		// if 1) we don't know the geometry exactly, 2)additional hits etc.
		ll_pion = pdf_pion->get_log_likelihood(sim_points_pion_inexact);
		ll_kaon = pdf_kaon->get_log_likelihood(sim_points_pion_inexact);
		if (locBin_PID == 0)
		{
			ll_diff_pion_inexact->Fill(1*(ll_pion-ll_kaon));
			Nph_pion -> Fill(int(sim_points_pion_inexact.size()));
			histMap_ll_diff_pion_sim_inexact[locBin_P]->Fill(1*(ll_pion-ll_kaon));
			histMap_Nhits_pion_sim_inexact[locBin_P]->Fill(int(sim_points_pion_inexact.size()));
		}

		ll_pion = pdf_pion->get_log_likelihood(sim_points_kaon_inexact);
		ll_kaon = pdf_kaon->get_log_likelihood(sim_points_kaon_inexact);
		if (locBin_PID == 1)
		{
			ll_diff_kaon_inexact->Fill(1*(ll_pion-ll_kaon));
			Nph_kaon -> Fill(int(sim_points_kaon_inexact.size()));
			histMap_ll_diff_kaon_sim_inexact[locBin_P]->Fill(1*(ll_pion-ll_kaon));
			histMap_Nhits_kaon_sim_inexact[locBin_P]->Fill(int(sim_points_kaon_inexact.size()));
		}


		sim_points_pion.clear();
		sim_points_kaon.clear();
		sim_points_pion_inexact.clear();
		sim_points_kaon_inexact.clear();
		tree_points.clear();
		delete pdf_pion;
		delete pdf_kaon;
	
	}// END OF PARTICLE LOOP

	}// END OF EVENT_CHAIN LOOP



	printf(" Run Completed \n");

	/**************************************************************************/
	/***********               Write histograms             *******************/
	/**************************************************************************/

	tfile->cd();
	pion_dist_x->Write();
	pion_dist_y->Write();
	pion_dist_xy->Write();
	pion_dist_xt->Write();
	pion_dist_yt->Write();
	pion_dist_t->Write();
	pion_dist_rowcol->Write();
	pion_dist_3D->Write();

	kaon_dist_x->Write();
	kaon_dist_y->Write();
	kaon_dist_xy->Write();
	kaon_dist_xt->Write();
	kaon_dist_yt->Write();
	kaon_dist_t->Write();
	kaon_dist_rowcol->Write();

        phot_found_pion->Write();
        phot_found_kaon->Write();

	hit_dist_rowcol_tree->Write();
	hit_dist_t_tree->Write();

	Nph_tree -> Write();
	Nph_pion -> Write();
	Nph_kaon -> Write();


        ll_diff_pion->Write();
        ll_diff_kaon->Write();

        ll_diff_pion_inexact->Write();
        ll_diff_kaon_inexact->Write();

        ll_diff_tree->Write();
        ll_diff_pion_tree->Write();
        ll_diff_kaon_tree->Write();

	for (size_t locPBin = 0 ; locPBin < PBins.size() ; locPBin++)
	{
		histMap_ll_diff_pion_tree[locPBin]->Write();
		histMap_ll_diff_kaon_tree[locPBin]->Write();
		histMap_ll_diff_pion_sim_exact[locPBin]->Write();
		histMap_ll_diff_kaon_sim_exact[locPBin]->Write();
		histMap_ll_diff_pion_sim_inexact[locPBin]->Write();
		histMap_ll_diff_kaon_sim_inexact[locPBin]->Write();

		histMap_Nhits_pion_tree[locPBin]->Write();
		histMap_Nhits_kaon_tree[locPBin]->Write();
		histMap_Nhits_pion_sim_inexact[locPBin]->Write();
		histMap_Nhits_kaon_sim_inexact[locPBin]->Write();
	}	

	tfile->Close();

	int status = 0;
	return status;

}
