#include "dirc_point.h"
#include "dirc_rect_digitizer.h"
#include <vector>
DircRectDigitizer::DircRectDigitizer(\
		double iminx,\
		double imaxx,\
		double iresx,\
		double iminy,\
		double imaxy,\
		double iresy,\
		double it_unc,\
		double it_bin_size) : DircDigitizer()
{
//Digitizes to a rectagular grid
	minx = iminx;
	maxx = imaxx;
	resx = iresx;
	miny = iminy;
	maxy = imaxy;
	resy = iresy;
	t_unc = it_unc;
	t_bin_size = it_bin_size;
	
	dig_rand = new TRandom3();

}
void DircRectDigitizer::digitize_point(dirc_point &pt)
{
		//overflow/underflow?
	double x = pt.x;
	double y = pt.y;
	int xdig = (x - minx)/resx;
	int ydig = (y - miny)/resy;
	
	double xout = resx*xdig + minx + resx/2;
	double yout = resy*ydig + miny + resy/2;
	
	if (x < minx)
	{
		xout = minx - resx/2;
	}
	else if (x > maxx)
	{
		xout = maxx + resx/2;
	}
	
	if (y < miny)
	{
		yout = miny - resy/2;
	}
	else if (y > maxy)
	{
		yout = maxy + resy/2;
	}
	
	pt.x = xout;
	pt.y = yout;
	pt.t += dig_rand->Gaus(0,t_unc);
	if (fabs(t_bin_size) > t_unc/100)
	{
		//Don't bother with binning if it's small
		//perhaps slow - a hard cut or a constant variable could be better here
		int tmp_t_bin = pt.t/(t_bin_size);
		pt.t = tmp_t_bin*t_bin_size + t_bin_size/2;
	}

	//Converting to as-built config
        int xhit = find_pixel_row(x);
        int yhit = find_pixel_col(y);
        //int yhit = (y - (-55.))/(PIXEL_SIZE + PIXEL_GAP);

	if (xhit == -1337 || yhit == -1337)
	{
		pt.ch = -999;
		pt.pixel_row = -999;
		pt.pixel_col = -999;
	}
	else
	{
		pt.pixel_row = xhit;
		pt.pixel_col = yhit;
		int pmt_id = (xhit/8) + (yhit/8) * PMT_ROWS;
		pt.ch = (pmt_id*CHANNEL_PER_PMT) + (yhit%8) * 8 + xhit%8;  
		//pt.pixel_col = PMT_COLUMNS*8 - yhit;
		//int pmt_id = (xhit/8) + (PMT_COLUMNS - 1 - (yhit/8)) * PMT_ROWS;
		//pt.ch = (pmt_id*CHANNEL_PER_PMT) + ( 7 - (yhit%8)) * 8 + xhit%8;  
	}

}
int DircRectDigitizer::find_pixel_row(double x)
{
	int pmt_row  = (24.904 - x)/PMT_SIZE_GAP;
	double pixel_x = 24.904 - x - pmt_row*PMT_SIZE_GAP;
	int pmt_pixel_row = (pixel_x + PIXEL_GAP)/(PIXEL_SIZE+PIXEL_GAP);
	if (pixel_x > PMT_SIZE)
		return -1337;
	else
		return pmt_row*8 + pmt_pixel_row;
}
int DircRectDigitizer::find_pixel_col(double y)
{
        int pmt_col  = y / PMT_SIZE_GAP;
        double pixel_y = y - pmt_col*PMT_SIZE_GAP;
        int pmt_pixel_col = (pixel_y + PIXEL_GAP)/(PIXEL_SIZE+PIXEL_GAP);
        if (pixel_y > PMT_SIZE)
                return -1337;
        else
                return pmt_col*8 + pmt_pixel_col;
}


void DircRectDigitizer::digitize_points(std::vector<dirc_point> &points)
{
	for (unsigned int i = 0; i < points.size(); i++)
	{
		//Wasting time Make it faster later TODO
		digitize_point(points[i]);
	}
	
}

	
