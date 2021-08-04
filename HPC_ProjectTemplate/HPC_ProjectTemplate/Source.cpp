#include <iostream>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header
#pragma once
#define HISTOGRAM 256
#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*******Read Image and save it to local arrayss***	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int *Red = new int[BM.Height * BM.Width];
	int *Green = new int[BM.Height * BM.Width];
	int *Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height*BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i*BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	return input;
}


void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);

	
	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			
			
			if (image[i*width + j] < 0)
			{
				image[i*width + j] = 0;
				
			}
			if (image[i*width + j] > 255)
			{
				image[i*width + j] = 255;
				
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}


int main()
{
	
	int ImageWidth = 4, ImageHeight = 4;
	
	int start_s, stop_s, TotalTime = 0;
	int start_s_parallel, stop_s_parallel, TotalTime_parallel = 0;

	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test.png";
	imagePath = marshal_as<System::String^>(img);
	int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
	int ImageSize = ImageWidth * ImageHeight;
	
	MPI_Init(NULL, NULL);
	start_s_parallel = clock();
	int size, rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int count = ImageSize / size;
	
	
	if (rank == 0)
	{
		//save the original image as index 0
		createImage(imageData, ImageWidth, ImageHeight, 0);
		
	}

	//<------------------intialize------------------>//


	int *recevesd_img = new int[count];
	int global_histo[HISTOGRAM] = { 0 };
	int local_histo[HISTOGRAM] = { 0 };

	MPI_Scatter(imageData, count, MPI_INT, recevesd_img, count, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < count; i++)
	{
		local_histo[recevesd_img[i]] += 1;
	}

	MPI_Reduce(&local_histo, &global_histo, HISTOGRAM, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	

	//<--------------------- Calculate probabilty--------------------->\\-
	int count_prop = HISTOGRAM / size;
	int no_of_pixels = ImageSize;
	double *local_probailty_histo = new double[count_prop];
	double global_probabilty_histo[HISTOGRAM] = {0};
	int *receved_prop = new int[count_prop];
	MPI_Scatter(&global_histo, count_prop, MPI_INT, receved_prop, count_prop, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < count_prop; i++)
	{
		local_probailty_histo[i] = receved_prop[i] / double(no_of_pixels);
	}

	MPI_Gather(local_probailty_histo, count_prop, MPI_DOUBLE, &global_probabilty_histo, count_prop, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//<--------------------- Calculate Comulative probabilty--------------------->\\-
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double cumulative_prob[HISTOGRAM];
	if (rank == 0)
	{


		for (int i = 0; i < HISTOGRAM; i++) {
			if (i == 0) {
				cumulative_prob[i] = global_probabilty_histo[i];
			}
			else {
				cumulative_prob[i] = cumulative_prob[i - 1] + global_probabilty_histo[i];
			}
		}


	}

	//<--------------------- Calculate floor |_ Comulative probabilty * scale _|--------------------->\\-
	int  cumulative_prob_mulitply[HISTOGRAM];
	int *local_cumulative_prob_mulitply= new int[count_prop];
	double *receved_cum = new double[count_prop];
	MPI_Scatter(&cumulative_prob, count_prop, MPI_DOUBLE, receved_cum, count_prop, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = 0; i < count_prop; i++)
	{
		local_cumulative_prob_mulitply[i] = floor(receved_cum[i] * 255.0);
	}

	MPI_Gather(local_cumulative_prob_mulitply, count_prop, MPI_INT, &cumulative_prob_mulitply, count_prop, MPI_INT, 0, MPI_COMM_WORLD);
	
	//<--------------------- map to the new image values --------------------->\\-
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		for (int i = 0; i < ImageSize; i++) {
			imageData[i] = cumulative_prob_mulitply[imageData[i]];
		}
		

		start_s = clock();
		createImage(imageData, ImageWidth, ImageHeight, size);
		stop_s = clock();
		stop_s_parallel= clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
		TotalTime_parallel += (stop_s_parallel - start_s_parallel) / double(CLOCKS_PER_SEC) * 1000;
		cout << "create image time: " << TotalTime << endl;
		cout << "parallel code time: " << TotalTime_parallel << endl;

		free(imageData); 
	}
	//<------------------end------------------>//

	MPI_Finalize();
	return 0;


}