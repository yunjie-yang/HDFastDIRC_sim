#include "dirc_geometry_gluex.h"

dirc_geometry_gluex::dirc_geometry_gluex()
{
	CHANNEL_PER_PMT = 64;
	PMT_ROWS = 18;
	PMT_COLUMNS = 6;
}

int
dirc_geometry_gluex::GetPmtID( int channel ) const 
{
	int MAX_BOX_CHANNEL = CHANNEL_PER_PMT * PMT_ROWS * PMT_COLUMNS;
	if(channel < MAX_BOX_CHANNEL) 
		return channel/CHANNEL_PER_PMT;
	else 
		return (channel-MAX_BOX_CHANNEL)/CHANNEL_PER_PMT;
}

int
dirc_geometry_gluex::GetPixelID( int channel ) const 
{
	return channel%CHANNEL_PER_PMT;
}

int dirc_geometry_gluex::GetPmtColumn( int channel ) const 
{
	int pmt = GetPmtID(channel);
	return pmt/PMT_ROWS; //  0 - 5
}

int dirc_geometry_gluex::GetPmtRow( int channel ) const 
{
	int pmt = GetPmtID(channel);
	return pmt%PMT_ROWS; //  0 - 17
}

int dirc_geometry_gluex::GetPmtPixelColumn( int channel ) const 
{
	int pix = GetPixelID(channel);
	return pix/8; //  0 - 7
}

int dirc_geometry_gluex::GetPmtPixelRow( int channel ) const 
{
	int pix = GetPixelID(channel);
	return pix%8; //  0 - 7
}

int dirc_geometry_gluex::GetPixelRow( int channel ) const 
{
	int pmt_row = GetPmtRow(channel);
	int pixel_row = GetPmtPixelRow(channel);
	return 8*pmt_row + pixel_row; //  0 - 143
}

int dirc_geometry_gluex::GetPixelColumn( int channel ) const 
{
	int pmt_column = GetPmtColumn(channel);
	int pixel_column = GetPmtPixelColumn(channel);
	return 8*pmt_column + pixel_column; // 0 - 47
}
			
int dirc_geometry_gluex::GetPixelX( int channel ) const 
{
	return abs(GetPixelRow(channel) - 143); //  0 - 143
}

int dirc_geometry_gluex::GetPixelY( int channel ) const 
{
	return 47 - GetPixelColumn(channel); // 0 - 47
}
