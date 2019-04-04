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

	double particle_x_hall     = 0;
	double particle_y_hall     = 0;
	double particle_theta_hall = 0;
	double particle_phi_hall   = 0;


	double particle_flight_distance = 0;

	bool kaleidoscope_plot = false;
	bool use_moliere_scattering = false;

	int num_runs = 1000;

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

	double tracking_unc = .0000*57.3; //mrad
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

	int n_phi_phots = 150000;
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

	double res_enhance = 1;
	TH1F *pion_dist_x = new TH1F("pion_dist_x","x val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *pion_dist_y = new TH1F("pion_dist_y","y val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *pion_dist_t = new TH1F("pion_dist_t","t val of intercepted points - pion",(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *pion_dist_xt = new TH2F("pion_dist_xt","xt val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *pion_dist_yt = new TH2F("pion_dist_yt","yt val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *pion_dist_rowcol = new TH2F("pion_dist_rowcol","hit pattern - pion; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);
	//TH2F *pion_dist_rowcol = new TH2F("pion_dist_rowcol","hit pattern - pion; Pixel Row ; Pixel Column",190,-10.5,180.5,70,-9.5,60.5);
	TH3F *pion_dist_3D = new TH3F("pion_dist_3D","(x,y,t) - pion; Pixel Row ; Pixel Column ; Hit Time (ns)",144,-0.5,143.5,48,-0.5,47.5,300,0,300);


	TH1F *kaon_dist_x = new TH1F("kaon_dist_x","x val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *kaon_dist_y = new TH1F("kaon_dist_y","y val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *kaon_dist_t = new TH1F("kaon_dist_t","t val of intercepted points - kaon",(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_xy = new TH2F("kaon_dist_xy","xy val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *kaon_dist_xt = new TH2F("kaon_dist_xt","xt val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_yt = new TH2F("kaon_dist_yt","yt val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_rowcol = new TH2F("kaon_dist_rowcol","hit pattern - kaon; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);

        TH1F *ll_diff_pion = new TH1F("ll_diff_pion","Difference of log likelihood real = pion",200000,-200,200);
        TH1F *ll_diff_kaon = new TH1F("ll_diff_kaon","Difference of log likelihood real = kaon",200000,-200,200);
        TH1F *phot_found_pion = new TH1F("phot_found_pion","number of photons found on pion angle", 1001,-.5,1000.5);
        TH1F *phot_found_kaon = new TH1F("phot_found_kaon","number of photons found on kaon angle", 1001,-.5,1000.5);

	TH2F *hit_dist_rowcol_tree = new TH2F("hit_dist_rowcol_tree","hit pattern; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);


	double pion_beta, kaon_beta;
	pion_beta=kaon_beta=-1;
	double pion_angle, kaon_angle;
	pion_angle=kaon_angle = -1;

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

	if (user_opts.Find("PARTICLE_BAR", opt_val)) particle_bar = double(opt_val[1]);
	if (user_opts.Find("E", opt_val)) energy = opt_val[1];
	if (user_opts.Find("PARTICLE_THETA", opt_val)) particle_theta_hall = opt_val[1];
	if (user_opts.Find("PARTICLE_PHI", opt_val)) particle_phi_hall = opt_val[1];
	if (user_opts.Find("PARTICLE_X", opt_val)) particle_x_hall = opt_val[1];
	if (user_opts.Find("PARTICLE_Y", opt_val)) particle_y_hall = opt_val[1];

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

	DircThreeSegBoxSim *dirc_model = new DircThreeSegBoxSim(\
			rseed,\
			mirror_r_nominal,\
			foc_mirror_size,\
			main_mirror_angle_nominal,\
			600,\
			pmt_angle_nominal,\
			geometry_infilename);


	dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	dirc_model->set_wedge_mirror_rand(wedge_non_uniformity);


	//various running condition controls
	dirc_model->set_store_traveled(false);// uses LOTS of memory if set to true.
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_kaleidoscope_plot(kaleidoscope_plot);
	dirc_model->set_use_quartz_n_for_liquid(use_quartz_for_liquid);
	dirc_model->set_use_moliere(use_moliere_scattering);
	dirc_model->set_moliere_p(energy*1000);//assume momentum is the same for both for now - high energy;
	dirc_model->set_liquid_absorbtion(liquid_absorbtion);
	dirc_model->set_liquid_index(liquid_index);


	/**************************************************************************/
	/***********       APPY OFFSETS TO NOMINAL GEOMETRY       *****************/
	/**************************************************************************/

	TRandom3 spread_ang(rseed+3);

	dirc_model->set_pmt_offset(pmt_offset);
	dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
	dirc_model->set_bar_box_angle(bar_box_box_angle);

	dirc_model->set_mirror_plane_offsets(foc_mirror_yoff,foc_mirror_zoff);


	dirc_model->set_pmt_offset(pmt_offset);
	dirc_model->set_pmt_angle(47.87+pmt_angle_offset);

	dirc_model->set_focus_mirror_angle(			\
					   spread_ang.Gaus(main_mirror_angle_nominal,mirror_angle_change_unc)+main_mirror_angle_off, \
					   spread_ang.Gaus(0,mirror_angle_change_yunc)+main_mirror_yangle_off, \
					   main_mirror_zangle_off);

	dirc_model->set_mirror_plane_offsets(			\
						 main_mirror_yoff,	\
						 main_mirror_zoff);

	dirc_model->set_bar_box_offsets(\
						 bar_box_xoff,\
						 bar_box_yoff,\
						 bar_box_zoff);

	dirc_model->set_geometry_outfile(geometry_outfilename);
	//dirc_model->print_model();

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
	double particle_z_tree     = 0.;
	
	int hit_ChannelId   = -1;
	int hit_pixelrow = -1;
	int hit_pixelcol = -1;

	TVector3 particle_p3_tree, particle_x3_tree;

	double loc_energy_pion(0.);
	double loc_energy_kaon(0.);

	std::vector<dirc_point> sim_points_pion;
	std::vector<dirc_point> sim_points_kaon;

	std::vector<dirc_point> hit_points_pion;
	std::vector<dirc_point> hit_points_kaon;

	std::vector<dirc_point> tree_points;

	double ll_pion(0.),ll_kaon(0.);

	int counter = 0;

	//for (int event_i = 0 ; event_i < 2000 ; event_i++)
	for (int event_i = 0 ; event_i < event_chain->GetEntries() ; event_i++)
	{
	event_chain->GetEntry(event_i);
	if(event_i%1000==0)
		printf("Event #%d Counter = %d\n",event_i,counter);

	for (int particle_i = 0; particle_i < particle_array->GetEntriesFast(); particle_i++)
	{
		
		evt = (DrcEvent*) particle_array->At(particle_i);
		particle_p3_tree = evt->GetMomentum();
		particle_x3_tree = evt->GetPosition();
		particle_mom_tree   = particle_p3_tree.Mag();
		particle_theta_tree = particle_p3_tree.Theta()*rad2deg; 
		particle_phi_tree   = particle_p3_tree.Phi()*rad2deg; 
		particle_x_tree     = particle_x3_tree.X();
		particle_y_tree     = particle_x3_tree.Y(); 
		particle_z_tree     = particle_x3_tree.Z(); 

		//if (particle_mom_tree < 3.8 || particle_mom_tree > 4.2) continue;
		if (particle_mom_tree < 4) continue;
		if (particle_x_tree > 5. || particle_x_tree < 0.) continue;
		//if (particle_y_tree > -12.0858 || particle_y_tree < -15.5858) continue;//bar 1 
		if (particle_y_tree > -19.1158 || particle_y_tree < -22.6158) continue;//bar 3

		counter++;
/*
		if (counter > 10)
		{
			printf("Counter = %d ; event_i = %d ; skipping...\n",counter,event_i);
			continue;
		}
*/
		for (auto hit : evt->GetHits())
		{
			dirc_point tree_point;
			hit_ChannelId   = hit.GetChannel();
			hit_pixelrow    = digitizer.GetPixelRow(hit_ChannelId);
			hit_pixelcol    = digitizer.GetPixelColumn(hit_ChannelId);
			hit_dist_rowcol_tree -> Fill(hit_pixelrow,hit_pixelcol); 
			tree_point.pixel_row = hit_pixelrow;
			tree_point.pixel_col = hit_pixelcol;
			tree_points.push_back(tree_point);
		}
	
		dirc_model->convert_particle_kinematics(particle_x,particle_y,particle_theta,particle_phi,particle_bar,\
							particle_x_tree,particle_y_tree,particle_theta_tree,particle_phi_tree);

	
		loc_energy_pion = sqrt(particle_mom_tree*particle_mom_tree + pimass*pimass);
		loc_energy_kaon = sqrt(particle_mom_tree*particle_mom_tree + kmass*kmass);

		pion_beta = dirc_model->get_beta(loc_energy_pion,pimass);
		kaon_beta = dirc_model->get_beta(loc_energy_kaon,pimass);

                dirc_model->sim_rand_n_photons(\
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

        	double x,y,t_ns;
        	int pixel_row;
        	for (unsigned int i = 0; i < sim_points_pion.size(); i++)
        	{
                	x = sim_points_pion[i].x;
               		y = sim_points_pion[i].y;
                	t_ns = sim_points_pion[i].t;
                	pixel_row = sim_points_pion[i].pixel_row;
                	if (pixel_row<0) continue;
                	pion_dist_x->Fill(x);
                	pion_dist_y->Fill(y);
                	pion_dist_t->Fill(t_ns);
                	pion_dist_xy->Fill(x,y);
                	pion_dist_xt->Fill(x,t_ns);
                	pion_dist_yt->Fill(y,t_ns);
                	pion_dist_t->Fill(t_ns);

                	pion_dist_rowcol->Fill(sim_points_pion[i].pixel_row,sim_points_pion[i].pixel_col);
                	pion_dist_3D->Fill(sim_points_pion[i].pixel_row,sim_points_pion[i].pixel_col,t_ns);
        	}
		sim_points_pion.clear();

/*
	        dirc_model->sim_reg_n_photons(\
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
 
		dirc_model->sim_reg_n_photons(\
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

		digitizer.undigitize_points(tree_points);

		ll_pion = pdf_pion->get_log_likelihood(tree_points);
		ll_kaon = pdf_kaon->get_log_likelihood(tree_points);

		ll_diff_pion->Fill(1*(ll_pion-ll_kaon));

		phot_found_pion->Fill(tree_points.size());
*/		


	}

	}



	/**************************************************************************/
	/***********        Define Particle Kinematics         ********************/
	/**************************************************************************/

/*
	dirc_model->convert_particle_kinematics(particle_x,particle_y,particle_theta,particle_phi,particle_bar,\
						particle_x_hall,particle_y_hall,particle_theta_hall,particle_phi_hall);
	

	printf("particle_bar     = %8.02f \n",particle_bar);
	printf("energy           = %8.02f \n",energy);
	printf("theta_hall       = %8.02f \n",particle_theta_hall);
	printf("phi_hall         = %8.02f \n",particle_phi_hall);
	printf("particle_x_hall  = %8.02f \n",particle_x_hall);
	printf("particle_y_hall  = %8.02f \n",particle_y_hall);
	printf("theta            = %8.02f \n",particle_theta);
	printf("phi              = %8.02f \n",particle_phi);
	printf("particle_x       = %8.02f \n",particle_x);
	printf("particle_y       = %8.02f \n",particle_y);
*/



	/**************************************************************************/
	/***********                      STEP 1               ********************/
	/**************************************************************************/
	/*
		Given particle track information (hit position, angle and momentum),
			generate expected hits under different hypotheses.
	*/

/*
	pion_beta = dirc_model->get_beta(energy,pimass);
	kaon_beta = dirc_model->get_beta(energy,kmass);

	double pion_time = particle_flight_distance/(pion_beta*.3);
	double kaon_time = particle_flight_distance/(kaon_beta*.3);

	std::vector<dirc_point> sim_points;
	std::vector<dirc_point> hit_points_pion;
	std::vector<dirc_point> hit_points_kaon;

	dirc_model->sim_reg_n_photons(\
			hit_points_pion,\
			n_phi_phots,\
			n_z_phots,\
			-1,\
			particle_bar,\
			particle_x,\
			particle_y,\
			pion_time,\
			particle_theta,\
			particle_phi,\
			0,\
			ckov_unc/1.,\
			pion_beta,\
			-1); 

	dirc_model->sim_reg_n_photons(\
			hit_points_kaon,\
			n_phi_phots,\
			n_z_phots,\
			-1,\
			particle_bar,\
			particle_x,\
			particle_y,\
			kaon_time,\
			particle_theta,\
			particle_phi,\
			0,\
			ckov_unc/1.,\
			kaon_beta,\
			-1); 

	if (SIM_ONLY)
		num_runs = 0;
*/


	/**************************************************************************/
	/***********                      STEP 2               ********************/
	/**************************************************************************/
	/*
		Given generated expected hits, construct corresponding 
			PDFs using KDE.
	*/

/*
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
*/

	/**************************************************************************/
	/***********                      STEP 3               ********************/
	/**************************************************************************/
	/*
		Given PDFs under different particle hypothesis, and the observed hits, 
			calculate the log-likelihoods(LLs) given each PDF.
		Then, one can form the difference between different hypotheses and 
			do particle identification.


		Here, we don't have real observed hits, we therefore simulate observed hits.
	*/

/*
	double ll_pion(0.),ll_kaon(0.);


	printf("\n\n");
	for (int i = 0; i < num_runs; i++)
	{

		printf("\r                                                                      ");
		printf("\r Running particle %8d/%d  ",i+1,num_runs);
		fflush(stdout);

		//suppose this is a pion
		dirc_model->sim_rand_n_photons(\
				sim_points,\
				n_sim_phots,\
				pion_angle,\
				particle_bar,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				tracking_unc,\
				ckov_unc,\
				pion_beta);
		digitizer.digitize_points(sim_points);

		ll_pion = pdf_pion->get_log_likelihood(sim_points);
		ll_kaon = pdf_kaon->get_log_likelihood(sim_points);

		ll_diff_pion->Fill(1*(ll_pion-ll_kaon));
		phot_found_pion->Fill(sim_points.size());

		//suppose this is a kaon
		dirc_model->sim_rand_n_photons(\
				sim_points,\
				n_sim_phots,\
				kaon_angle,\
				particle_bar,\
				particle_x,\
				particle_y,\
				kaon_time,\
				particle_theta,\
				particle_phi,\
				tracking_unc,\
				ckov_unc,\
				kaon_beta);
		digitizer.digitize_points(sim_points);

		ll_pion = pdf_pion->get_log_likelihood(sim_points);
		ll_kaon = pdf_kaon->get_log_likelihood(sim_points);

		ll_diff_kaon->Fill(1*(ll_pion-ll_kaon));
		phot_found_kaon->Fill(sim_points.size());



	}
*/
	printf(" Run Completed \n");

	/**************************************************************************/
	/***********              Fill hit histograms           *******************/
	/**************************************************************************/
/*
	if (SIM_ONLY)
	{
		digitizer.digitize_points(hit_points_pion);
		digitizer.digitize_points(hit_points_kaon);
	}

	double x,y,t_ns;
	int pixel_row;
	for (unsigned int i = 0; i < hit_points_pion.size(); i++)
	{
		x = hit_points_pion[i].x;
		y = hit_points_pion[i].y;
		t_ns = hit_points_pion[i].t;
		pixel_row = hit_points_pion[i].pixel_row;
		if (SIM_ONLY && pixel_row<0) continue;
		pion_dist_x->Fill(x);
		pion_dist_y->Fill(y);
		pion_dist_t->Fill(t_ns);
		pion_dist_xy->Fill(x,y);
		pion_dist_xt->Fill(x,t_ns);
		pion_dist_yt->Fill(y,t_ns);
		pion_dist_t->Fill(t_ns);

		pion_dist_rowcol->Fill(hit_points_pion[i].pixel_row,hit_points_pion[i].pixel_col);
		pion_dist_3D->Fill(hit_points_pion[i].pixel_row,hit_points_pion[i].pixel_col,t_ns);
	}

	for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
	{
		x = hit_points_kaon[i].x;
		y = hit_points_kaon[i].y;
		t_ns = hit_points_kaon[i].t;
		pixel_row = hit_points_kaon[i].pixel_row;
		if (SIM_ONLY && pixel_row<0) continue;
		kaon_dist_x->Fill(x);
		kaon_dist_y->Fill(y);
		kaon_dist_t->Fill(t_ns);
		kaon_dist_xy->Fill(x,y);
		kaon_dist_xt->Fill(x,t_ns);
		kaon_dist_yt->Fill(y,t_ns);
		kaon_dist_t->Fill(t_ns);

		kaon_dist_rowcol->Fill(hit_points_kaon[i].pixel_row,hit_points_kaon[i].pixel_col);
	}
*/



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

        ll_diff_pion->Write();
        ll_diff_kaon->Write();
        phot_found_pion->Write();
        phot_found_kaon->Write();

	hit_dist_rowcol_tree->Write();

	tfile->Close();

	int status = 0;
	return status;

}
