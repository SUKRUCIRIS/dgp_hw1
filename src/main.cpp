// main.cpp
#include "../third_party/Easy3D/easy3d/util/initializer.h"
#include "hw2.h"

int main(int argc, char **argv)
{
	easy3d::initialize();
	Hw2Viewer viewer("Sukru HW2 - Kernel & Clustering");
	return viewer.run();
}