#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
#include <Windows.h>
#pragma pack(1)
#define THR 3


struct f_info
{
	unsigned char signature[2];
	unsigned int sizefile;
	unsigned int reserved;
	unsigned int addr_offset;
};

struct pic_info
{
	unsigned int Size;
	unsigned int Width;
	unsigned int Height;
	unsigned short int  Planes;
	unsigned short int  BitCount;
	unsigned int Compression;
	unsigned int SizeImage;
	unsigned int XPelsPerMeter;
	unsigned int YPelsPerMeter;
	unsigned int ClrUsed;
	unsigned int ClrImportant;
};

void Marra_Xilderta(float* First, float* Second, int height, int width) {
	int window[9] = { 1,1,1,
		1,2,1,
		1,1,1 };

	int window1[9] = { 1,1,1,
		1,-8,1,
		1,1,1 };

	int Mx = 3, My = 3, r = 9, g = 9, winf, Pm, Qm, winf1;

	winf = 0;
	winf1 = 0;
	while (r--)winf += window[r];
	if (winf == 0)winf = 1;

	while (g--)winf1 += window1[r];
	if (winf1 == 0)winf1 = 1;

	Pm = Mx / 2;
	Qm = My / 2;
	for (unsigned int i = 1; i < height - 1; i++) {
		for (unsigned int j = 1; j < width - 1; j++) {
			r = 0;
			float SumI = 0;
			for (int l = -Qm; l <= Qm; l++) {
				for (int k = -Pm; k <= Pm; k++) {
					float I;
					I = First[(i + l)*width + (j + k)];
					SumI += I*window[r++];
				}
			}
			SumI = SumI / (winf);
			Second[i*width + j] = SumI;
		}
	}
	for (unsigned int i = 1; i < height - 1; i++) {
		for (unsigned int j = 1; j < width - 1; j++) {
			g = 0;
			float SumI = 0;
			for (int l = -Qm; l <= Qm; l++) {
				for (int k = -Pm; k <= Pm; k++) {
					float I;
					I = Second[(i + l)*width + (j + k)];
					SumI += I*window1[g++];
				}
			}

			First[i*width + j] = SumI;
		}
	}
}

void null_lavel(float* First, float* Second, int height, int width, int param) {
	int Mx = 3, My = 3, Pm, Qm, k;
	Pm = Mx / 2;
	Qm = My / 2;
	for (unsigned int i = 1; i<height - 1; i++) {
		for (unsigned int j = 1; j<width - 1; j++) {
			k = 0;
			if ((First[(i - 1)*width + (j - 1)] * First[(i + 1)*width + (j + 1)] < 0) && (fabs(First[(i - 1)*width + (j - 1)] - First[(i + 1)*width + (j + 1)]) > param))
				k++;
			if ((First[(i - 1)*width + (j)] * First[(i + 1)*width + (j)] < 0) && (fabs(First[(i - 1)*width + (j)] - First[(i + 1)*width + (j)]) > param))
				k++;
			if ((First[(i - 1)*width + (j + 1)] * First[(i + 1)*width + (j - 1)] < 0) && (fabs(First[(i - 1)*width + (j + 1)] - First[(i + 1)*width + (j - 1)]) > param))
				k++;
			if ((First[(i)*width + (j - 1)] * First[(i)*width + (j + 1)] < 0) && (fabs(First[(i)*width + (j - 1)] - First[(i)*width + (j + 1)]) > param))
				k++;
			if (k > 1) {
				Second[i*width + j] = 200;
			}
			else {
				Second[i*width + j] = 0;
			}
		}
	}


}

int main(int argc, char *argv[])
{
	int start_time, end_time;
	char infile[100] = "kartinka.bmp";
	char outfile[100] = "new1.bmp";
	if (argc>1) strcpy_s(infile, argv[1]);
	if (argc>2) strcpy_s(outfile, argv[2]);

	std::ifstream in(infile, std::ios::in | std::ios::binary);
	if (!in) { std::cout << "File not found..."; exit(1); }
	struct f_info f_i;
	struct pic_info pic_i;
	in.read((char*)&f_i, sizeof(f_info));
	in.read((char*)&pic_i, sizeof(pic_info));
	int ln_str = (f_i.sizefile - 54) / pic_i.Height;
	float *image = new float[pic_i.Height*pic_i.Width];
	std::cout << "Width - " << pic_i.Width << ", Height - " << pic_i.Height << std::endl;

	for (unsigned int i = 0; i<pic_i.Height; i++)
	{
		for (unsigned int j = 0; j<pic_i.Width; j++)
		{
			unsigned char R, G, B;
			char C;
			in.get(C);
			B = C;
			in.get(C);
			G = C;
			in.get(C);
			R = C;
			image[i*pic_i.Width + j] = B*0.114 + G*0.587 + R*0.299;
		}
		for (unsigned int k = 0; k<(ln_str - pic_i.Width * 3); k++)
		{
			char d;
			in.get(d);
		}
	}
	in.close();

	float *image2 = new float[pic_i.Height*pic_i.Width];
	float *image3 = new float[pic_i.Height*pic_i.Width];

	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	LARGE_INTEGER time1, time2;
	QueryPerformanceCounter(&time1);

	Marra_Xilderta(image, image2, pic_i.Height, pic_i.Width);
	null_lavel(image, image3, pic_i.Height, pic_i.Width, 20);

	QueryPerformanceCounter(&time2);
	time2.QuadPart -= time1.QuadPart;
	double span = (double)time2.QuadPart / freq.QuadPart;

	std::ofstream out;
	out.open(outfile, std::ios::out | std::ios::binary);
	out.write((char*)&f_i, sizeof(f_info));
	out.write((char*)&pic_i, sizeof(pic_info));
	for (unsigned int i = 0; i<pic_i.Height; i++)
	{
		for (unsigned int j = 0; j<pic_i.Width; j++)
		{
			out.put(image3[i*pic_i.Width + j]);
			out.put(image3[i*pic_i.Width + j]);
			out.put(image3[i*pic_i.Width + j]);
		}
		for (unsigned int k = 0; k<(ln_str - pic_i.Width * 3); k++) out.put(0);
	}

	std::cout << "Press any key...\n";
	std::cout << "Width - " << pic_i.Width << ", Height - " << pic_i.Height << std::endl;
	std::cout << "Time of working programm: " << span << " sec. \n";
	delete image2;
	delete image3;
	system("pause");
	return 0;
}
