// hw1.h
#ifndef HW1_H
#define HW1_H

#include "../third_party/Easy3D/easy3d/viewer/viewer.h"
#include "../third_party/Easy3D/easy3d/core/surface_mesh.h"
#include "../third_party/Easy3D/easy3d/renderer/drawable_lines.h"

#include <vector>
#include <atomic>
#include <mutex>
#include <functional>
#include <string>

class Hw1Viewer : public easy3d::Viewer
{
public:
	Hw1Viewer(const std::string &title);
	~Hw1Viewer();

	void init() override;
	void post_draw() override;

private:
	char filepath[512] = "";
	char filepath_orig[512] = "";

	std::vector<easy3d::SurfaceMesh *> current_meshes;
	std::vector<easy3d::SurfaceMesh *> temp_meshes;

	easy3d::LinesDrawable *geodesic_drawable;

	int start_v = 0;
	int end_v = 1;
	int subdiv_type = 0;
	int subdiv_iters = 1;

	std::atomic<bool> is_processing{false};
	std::atomic<bool> task_finished{false};
	std::string processing_msg;
	std::mutex task_mutex;
	std::function<void()> on_task_finished;

	std::vector<float> temp_dist;
	std::vector<easy3d::SurfaceMesh::Vertex> temp_parent;

	template <typename Func, typename OnFinished>
	void run_async(std::string msg, Func background_task, OnFinished main_thread_update);

	void process_pending_tasks();
	void visualize_path(easy3d::SurfaceMesh::Vertex start, easy3d::SurfaceMesh::Vertex end);
};

#endif