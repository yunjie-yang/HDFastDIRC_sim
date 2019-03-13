#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "src/dirc_optical_sim.h"
#include "src/dirc_threesegbox_sim.h"
#include "src/dirc_point.h"
#include "src/dirc_probability_spread.h"
#include "src/dirc_probability_separation.h"
#include "src/dirc_spread_radius.h"
#include "src/dirc_spread_relative.h"
#include "src/dirc_spread_linear_soft.h"
#include "src/dirc_spread_gaussian.h"
#include "src/dirc_digitizer.h"
#include "src/dirc_rect_digitizer.h"
#include "src/dirc_babar_digitizer.h"
#include "src/dirc_babar_sim.h"
#include "src/dirc_progressive_separation.h"
#include "src/dirc_gluex_lut_enum.h"
#include "src/dirc_lut_enum.h"
#include "src/dirc_lut.h"
#include <TFile.h>
#include <TTree.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMinuit.h>


int main(int nargs, char* argv[])
{

	const char* in_str;

	double energy = 5.0;
	double kmass = .4937;
	double pimass = .1396;

	double particle_x = 0;
	double particle_y = 0;
	double particle_theta = 0;
	double particle_phi = 0;

	double particle_flight_distance = 0;

	bool kaleidoscope_plot = false;
	bool use_moliere_scattering = false;

	int num_runs = 1000;

	double mean_n_phot = 40;
	double spread_n_phot = 0;

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


	double pmt_min_z = -1000;
	double pmt_max_z = 1000;
	double large_mirror_min_z = -1000;
	double large_mirror_max_z = 1000;

	pmt_min_z = -559;
	pmt_max_z = -329;
	large_mirror_min_z = -559;
	large_mirror_max_z = -130;

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

	//Sets the side boundarys of the distributions
	double sm_xl = -10000000;
	double sm_xr = -sm_xl;


	int n_sim_phots = 40;

	int n_phi_phots = 150000;
	int n_z_phots = 4;

/*
	sm_xl = -50;
	sm_xr = sm_xl + 440;
*/

	sm_xl = -45;
	sm_xr = sm_xl + 1000;

	bool use_quartz_for_liquid = false;
	bool three_seg_mirror = true;

	double foc_mirror_yoff = 0;
	double foc_mirror_zoff = 0;


	double liquid_absorbtion = 0*-log(.7)/1000;
	double liquid_index = 1.33;

	char* rootfilename = new char[256];
	sprintf(rootfilename,"fitdirc.root");

	printf("Arguments Passed=%d\n",nargs);

	if(nargs==2){
		in_str = argv[1];
		printf("%s\n",in_str);
		printf("nargs=2\n");
	}
	else{
		for (int i = 1; i < nargs; i++)
		{
			if (strcmp(argv[i], "-of") == 0)
			{
				i++;
				sprintf(rootfilename,"%s",argv[i]);
			}
			else if (strcmp(argv[i], "-particle_phi") == 0)
			{
				i++;
				particle_phi = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_flight_distance") == 0)
			{
				//meters
				i++;
				particle_flight_distance = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-tracking_unc") == 0)
			{
				i++;
				tracking_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-ckov_unc") == 0)
			{
				i++;
				ckov_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-side_mirror") == 0)
			{
				i++;
				sm_xl = atof(argv[i]);
				i++;
				sm_xr = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_theta") == 0)
			{
				i++;
				particle_theta = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-use_moliere_scattering") == 0)
			{
				use_moliere_scattering = true;
				printf("Enabling Moliere Scattering - only implemented in loop mode\n");
			}
			else if (strcmp(argv[i], "-kaleidoscope_plot") == 0)
			{
				kaleidoscope_plot = true;
			}
			else if (strcmp(argv[i], "-open_image_plane") == 0)
			{
				pmt_min_z = -1000;
				pmt_max_z = 1000;
				large_mirror_min_z = -1000;
				large_mirror_max_z = 1000;
			}
			else if (strcmp(argv[i], "-mean_n_phot") == 0)
			{
				i++;
				mean_n_phot = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-t_unc") == 0)
			{
				i++;
				t_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-t_bin_size") == 0)
			{
				i++;
				t_bin_size = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-n_phi_phots") == 0)
			{
				i++;
				n_phi_phots = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-n_z_phots") == 0)
			{
				i++;
				n_z_phots = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_y") == 0)
			{
				i++;
				particle_y = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_x") == 0)
			{
				i++;
				particle_x = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-E") == 0)
			{
				i++;
				energy = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-n") == 0)
			{
				i++;
				num_runs = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-pmt_res") == 0)
			{
				i++;
				resx = atof(argv[i]);
				resy = resx;
			}
			else if (strcmp(argv[i], "-liquid_index") == 0)
			{
				i++;
				liquid_index = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-wedge_uncertainty") == 0)
			{
				i++;
				wedge_uncertainty = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change") == 0)
			{
				i++;
				mirror_angle_change = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change_unc") == 0)
			{
				i++;
				mirror_angle_change_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change_yunc") == 0)
			{
				i++;
				mirror_angle_change_yunc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_rot") == 0)
			{
				i++;
				box_rot = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_box_angle") == 0)
			{
				i++;
				bar_box_box_angle = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_r_difference") == 0)
			{
				i++;
				mirror_r_difference = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-wedge_non_uniformity") == 0)
			{
				i++;
				wedge_non_uniformity = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-pmt_offset") == 0)
			{
				i++;
				pmt_offset = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-pmt_angle_offset") == 0)
			{
				i++;
				pmt_angle_offset = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_nonuniformity") == 0)
			{
				i++;
				main_mirror_nonuniformity = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_angle_off") == 0)
			{
				i++;
				main_mirror_angle_off = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_yangle_off") == 0)
			{
				i++;
				main_mirror_yangle_off = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_zangle_off") == 0)
			{
				i++;
				main_mirror_zangle_off = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_yoff") == 0)
			{
				i++;
				main_mirror_yoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_zoff") == 0)
			{
				i++;
				main_mirror_zoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-rseed") == 0)
			{
				i++;
				rseed = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_xoff") == 0)
			{
				i++;
				bar_box_xoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_yoff") == 0)
			{
				i++;
				bar_box_yoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_zoff") == 0)
			{
				i++;
				bar_box_zoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-foc_mirror_yoff") == 0)
			{
				i++;
				foc_mirror_yoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-foc_mirror_zoff") == 0)
			{
				i++;
				foc_mirror_zoff = atof(argv[i]);
			}
			else
			{
				printf("Unrecognized argument: %s\n",argv[i]);
			}
		}
	}


	double main_mirror_angle = 74.11+mirror_angle_change;


	double rad_to_deg = 57.2958;

	double res_enhance = 1;


	TRandom3 spread_ang(rseed+3);

	DircThreeSegBoxSim *dirc_model = new DircThreeSegBoxSim(\
			rseed,\
			-1200 + mirror_r_difference,\
			foc_mirror_size,\
			main_mirror_angle,\
			600,\
			47.87 + box_rot + mirror_angle_change);

	dirc_model->set_store_traveled(false);// uses LOTS of memory if set to true.
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_wedge_mirror_rand(wedge_non_uniformity);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_kaleidoscope_plot(kaleidoscope_plot);
	dirc_model->set_pmt_plane_zs(pmt_min_z,pmt_max_z);
	dirc_model->set_large_mirror_zs(large_mirror_min_z,large_mirror_max_z);
	dirc_model->set_use_quartz_n_for_liquid(use_quartz_for_liquid);


	double pion_beta, kaon_beta;
	pion_beta=kaon_beta=-1;
	double pion_angle, kaon_angle;
	pion_angle=kaon_angle = -1;

	TFile* tfile = new TFile(rootfilename,"RECREATE");

	TH1F *pion_dist_x = new TH1F("pion_dist_x","x val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *pion_dist_y = new TH1F("pion_dist_y","y val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *pion_dist_t = new TH1F("pion_dist_t","t val of intercepted points - pion",(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *pion_dist_xt = new TH2F("pion_dist_xt","xt val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *pion_dist_yt = new TH2F("pion_dist_yt","yt val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *pion_dist_rowcol = new TH2F("pion_dist_rowcol","hit pattern - pion; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);

	DircRectDigitizer digitizer(\
			minx,\
			maxx,\
			resx,\
			digit_miny,\
			digit_maxy,\
			resy,\
			t_unc,\
			t_bin_size);
	printf("Beginning Run\n");
	std::vector<dirc_point> sim_points;
	dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	if (num_runs > 0)
	{

		printf("No input file specified.  Running in loop mode\n");


		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		std::vector<dirc_point> hit_points_pion;

		dirc_model->set_use_moliere(use_moliere_scattering);
		dirc_model->set_moliere_p(energy*1000);//assume momentum is the same for both for now - high energy;

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);


		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		dirc_model->set_mirror_plane_offsets(foc_mirror_yoff,foc_mirror_zoff);


		//ns
		double pion_time = particle_flight_distance/(pion_beta*.3);


		//--------------- Applying offsets ----------------//

		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_pmt_angle(47.87+pmt_angle_offset);

		dirc_model->set_focus_mirror_angle(			\
						   spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc)+main_mirror_angle_off, \
						   spread_ang.Gaus(0,mirror_angle_change_yunc)+main_mirror_yangle_off, \
						   main_mirror_zangle_off);

		dirc_model->set_mirror_plane_offsets(			\
					 		 main_mirror_yoff,	\
					 		 main_mirror_zoff);

		dirc_model->set_bar_box_offsets(\
					 		 bar_box_xoff,\
					 		 bar_box_yoff,\
					 		 bar_box_zoff);



		dirc_model->sim_reg_n_photons(\
				hit_points_pion,\
				n_phi_phots,\
				n_z_phots,\
				-1,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/1.,\
				pion_beta,\
				-1); //stores values for next generated model
		/*
		dirc_model->sim_rand_n_photons(\
				//sim_points,
				hit_points_pion,\
				n_sim_phots,\
				pion_angle,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				tracking_unc,\
				ckov_unc,\
				pion_beta);
		*/

		digitizer.digitize_points(hit_points_pion);


		printf("\nRun Completed\n");

		double x,y,t_ns;
		for (unsigned int i = 0; i < hit_points_pion.size(); i++)
		{
			x = hit_points_pion[i].x;
			y = hit_points_pion[i].y;
			t_ns = hit_points_pion[i].t;
			pion_dist_x->Fill(x);
			pion_dist_y->Fill(y);
			pion_dist_t->Fill(t_ns);
			pion_dist_xy->Fill(x,y);
			pion_dist_xt->Fill(x,t_ns);
			pion_dist_yt->Fill(y,t_ns);
			pion_dist_t->Fill(t_ns);

			pion_dist_rowcol->Fill(hit_points_pion[i].pixel_row,hit_points_pion[i].pixel_col);
		}

	}

	tfile->cd();
	pion_dist_x->Write();
	pion_dist_y->Write();
	pion_dist_xy->Write();
	pion_dist_xt->Write();
	pion_dist_yt->Write();
	pion_dist_t->Write();
	pion_dist_rowcol->Write();

	tfile->Close();

	int status = 0;
	return status;

}
