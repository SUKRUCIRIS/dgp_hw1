#include "../third_party/Easy3D/easy3d/viewer/viewer.h"
#include "../third_party/Easy3D/easy3d/util/initializer.h"
#include "../third_party/Easy3D/easy3d/core/surface_mesh.h"
#include "../third_party/Easy3D/easy3d/fileio/surface_mesh_io.h"

int main(int argc, char **argv)
{
	std::cout << "Initializing Easy3D..." << std::endl;
	easy3d::initialize();

	easy3d::Viewer viewer("Sukru HW1");

	std::cout << "Attempting to load horse.obj..." << std::endl;

	easy3d::SurfaceMesh *mesh = easy3d::SurfaceMeshIO::load("./meshes1/1) use for geodesic/timing/dragon.obj");

	if (!mesh)
	{
		std::cerr << "[!] ERROR: Failed to load! The path is wrong or the file is missing." << std::endl;
		return -1;
	}

	std::cout << "Success! Loaded " << mesh->n_vertices() << " vertices and " << mesh->n_faces() << " faces." << std::endl;

	viewer.add_model(mesh);

	std::cout << "Launching viewer..." << std::endl;
	return viewer.run();
}