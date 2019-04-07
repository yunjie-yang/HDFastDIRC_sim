#include <vector>
#include <utility>
#include <TRandom3.h>

#include "dirc_base_sim.h"
#include "dirc_point.h"

#ifndef DIRC_THREESEGBOX_SIM
#define DIRC_THREESEGBOX_SIM 
class DircThreeSegBoxSim : public DircBaseSim
{
protected:
	// three-seg mirrors / focMirror

	bool three_seg_mirror;

	double foc_r;
	double foc_mirror_size;
	double foc_rot;

        double focMirrorBottom;
        double focMirrorTop;
        double focMirrorZDim;

        double threeSeg1Nx,threeSeg1Ny,threeSeg1Nz,threeSeg1D;
        double threeSeg2Nx,threeSeg2Ny,threeSeg2Nz,threeSeg2D;
        double threeSeg3Nx,threeSeg3Ny,threeSeg3Nz,threeSeg3D;

        double threeSeg1_2dny,threeSeg1_2dnz,threeSeg1_2dd;
        double threeSeg2_2dny,threeSeg2_2dnz,threeSeg2_2dd;
        double threeSeg3_2dny,threeSeg3_2dnz,threeSeg3_2dd;

        double threeSeg1Y,threeSeg1Z;
        double threeSeg2Y,threeSeg2Z;
        double threeSeg3Y,threeSeg3Z;
        double threeSeg1Y_end,threeSeg1Z_end;
        double threeSeg2Y_end,threeSeg2Z_end;
	double threeSeg3Y_end,threeSeg3Z_end;

	double threeSeg_theta_1;
	double threeSeg_theta_2;
	double threeSeg_theta_3;
	double seg_h;
	
	double focYoff;
	double focZoff;

	double focYoff_threeSeg1;
	double focZoff_threeSeg1;

	double focYoff_threeSeg2;
	double focZoff_threeSeg2;

	double focYoff_threeSeg3;
	double focZoff_threeSeg3;

	double foc_yrot;
	double foc_zrot;

	double foc_yrot_threeSeg1;
	double foc_zrot_threeSeg1;

	double foc_yrot_threeSeg2;
	double foc_zrot_threeSeg2;

	double foc_yrot_threeSeg3;
	double foc_zrot_threeSeg3;

	bool nonUniformFocMirror;
	double foc_mirror_nonuni;

	// PMT Plane / sensPlane
	double sens_size;
	double sens_rot;

	double pmtPlaneMinZ;
	double pmtPlaneMaxZ;

	double pmtPlane_xl;
	double pmtPlane_xr;

	double sensPlaneYdistConversion;
	double sensPlaneZdistConversion;

	double sensPlaneNx;
	double sensPlaneNy;
	double sensPlaneNz;
	double sensPlaneD;
	double sensPlaneY;
	double sensPlaneZ;
	
	double unReflSensPlaneNx;
	double unReflSensPlaneNy;
	double unReflSensPlaneNz;
	double unReflSensPlaneD;
	double unReflSensPlaneY;
	double unReflSensPlaneZ;


	// large Planar mirror / FTMS(N)
        double largePlanarMirrorNx;
        double largePlanarMirrorNy;
        double largePlanarMirrorNz;
        double largePlanarMirrorD;
        double largePlanarMirrorMinZ;
        double largePlanarMirrorMaxZ;


	// side mirrors
	double sidemirror_xr;
	double sidemirror_xl;
	double sidemirror_reflectivity;


	//additions to adapt to GlueX geometry:
	double upperWedgeMirrorTop;
	double largePlanarMirrorY;

	double distDCBR10DCBR11;

	// I/O
	char* geometry_outfile  = new char[256];
	char* geometry_infile   = new char[256];





	/* ----------------------------------------------------- */	

        double focPlaneNx;
        double focPlaneNy;
        double focPlaneNz;
        double focPlaneD;
        double focPlaneMinZ;

	double boxCloseZ;
	double reflOff;
	double baseReflOff;
	
	double focMirrorY;
	double focMirrorZ;
	
	double quartzLiquidY;
	
	double box_angle_off_cval;
	double box_angle_off_sval;

	double liquidAbsorbtion;
	std::vector<double> dist_traveled;
	bool store_traveled;
	bool kaleidoscope_plot;
	
	bool store_refraction;
	std::vector<double> refraction_before;
	std::vector<double> refraction_after;

	bool storeOpticalAngles;
	std::vector<double> focus_photon_angles;
	std::vector<double> side_photon_angles;
	std::vector<double> large_flat_photon_angles;
	
	double min_QE,max_QE,sep_QE;
	int num_QE;
	std::vector<double> vals_QE;

	double min_transmittance,max_transmittance,sep_transmittance;
	int num_transmittance;

	bool absorbtion_mc(double dx, double dy);
	void build_readout_box();
	void fill_sens_plane_vecs();
	void fill_threeseg_plane_vecs();
	void fill_other_mirrors_plane_vecs();
	void fill_foc_mirror_vecs();
	void sidemirror_reflect_points(std::vector<dirc_point> &points);
	void spread_wedge_mirror();

	void warp_readout_box(\
		dirc_point &out_val,\
		int particle_bar,\
		double &mm_index,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	double cylindrical_reflect(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	double three_seg_reflect(\
                double &x,\
                double &y,\
                double &z,\
                double &dx,\
                double &dy,\
                double &dz); 


	//Technically one of the "warp" functions, but can and should be optimized somehow
	double warp_sens_plane(\
		dirc_point &fill_val,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	double warp_box(\
                double &x,\
                double &y,\
                double &z,\
                double &dx,\
                double &dy,\
                double &dz);

public:
	double get_cerenkov_angle_rand(double beta, double additional_spread, double &wavelength);
	
	void set_focmirror_nonuniformity(double nonuni_deg);
	void set_foc_mirror_r(double ifoc_r);
	void set_sidemirror(double ixr, double ixl);
	void set_sidemirror_reflectivity(double isr);
	void sidemirror_reflect_point(dirc_point &ipt);
	void set_three_seg_mirror(bool itsm);
	void set_pmt_offset(double r);
	void set_liquid_absorbtion(double iabs);
	std::vector<double> get_dist_traveled();
	void set_store_traveled(bool sst = true);
	void set_focus_mirror_angle(double ang,double yang = 0, double zang = 0);
	void set_pmt_angle(double ang);
	void set_pmt_plane_zs(double imin, double imax);
	void set_large_mirror_zs(double imin, double imax);
	void set_mirror_plane_offsets(double off_y, double off_z);

        bool convert_particle_kinematics(double &particle_x,\
                                        double &particle_y,\
                                        double &particle_theta,\
                                        double &particle_phi,\
                                        double &particle_bar,\
                                        double particle_x_hall,\
                                        double particle_y_hall,\
                                        double particle_theta_hall,\
                                        double particle_phi_hall);


	
	void set_store_optical_angles(bool ibool);
	std::vector<double> get_focus_photon_angles();
	std::vector<double> get_side_photon_angles();
	std::vector<double> get_large_flat_photon_angles();
	DircThreeSegBoxSim(\
		int rand_seed = 4357,\
		double ifoc_r = -1200, \
		double ifoc_mirror_size = 288, \
		double ifoc_rot = 74.11, \
		double isens_size = 600, \
		double isens_rot = 47.87,\
		const char* igeometry_infile = "",\
                double ibar_length=4900,\
		double ibar_width=35,\
                double ibar_depth=17.25,\
                double iupper_wedge_top = 178.6);

	void print_model();
	void set_geometry_outfile(const char* filename);
};
#endif
