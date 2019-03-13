#include "dirc_point.h"
#include "dirc_digitizer.h"
#include <vector>
#include <iostream>

//Root is infecting more files :(
#include <TRandom3.h>

#ifndef DIRC_RECT_DIGITIZER
#define DIRC_RECT_DIGITIZER
class DircRectDigitizer : public DircDigitizer
{
private:
	double minx,maxx,miny,maxy;
	double resx, resy;
	double t_unc;
	double t_bin_size;
	
	TRandom3 *dig_rand;
        int CHANNEL_PER_PMT = 64;
        int PMT_ROWS = 18;
        int PMT_COLUMNS = 6;

	
public:
	DircRectDigitizer(\
		double iminx,\
		double imaxx,\
		double iresx,\
		double iminy,\
		double imaxy,\
		double iresy,\
		double it_unc,\
		double it_bin_size);

	void digitize_point(dirc_point &pt);
	void digitize_points(std::vector<dirc_point> &points);
};
#endif
