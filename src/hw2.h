// hw2.h
#ifndef HW2_H
#define HW2_H

#include "../third_party/Easy3D/easy3d/viewer/viewer.h"
#include "../third_party/Easy3D/easy3d/core/surface_mesh.h"
#include "../third_party/Easy3D/easy3d/renderer/drawable_points.h"

#include <vector>
#include <atomic>
#include <mutex>
#include <functional>
#include <string>
#include <memory>

struct GridSample
{
	easy3d::vec3 p;
	int ix, iy, iz;
	int parity;
	int cluster;
};

class Hw2Viewer : public easy3d::Viewer
{
public:
	Hw2Viewer(const std::string &title);
	~Hw2Viewer();

	void init() override;
	void post_draw() override;

private:
	char filepath[512] = "";
	std::shared_ptr<easy3d::SurfaceMesh> current_mesh = nullptr;
	easy3d::PointsDrawable *points_drawable = nullptr;

	int grid_res = 20;
	int selected_cluster = -1;
	bool show_highest_parity = false;
	float discard_dist = 1.0f; // Stored as percentage of Bounding Box diagonal now

	std::vector<GridSample> valid_samples;
	std::vector<easy3d::vec3> cluster_colors;
	int max_parity = 0;
	int num_clusters = 0;

	std::atomic<bool> is_processing{false};
	std::atomic<bool> task_finished{false};
	std::string processing_msg;
	std::mutex task_mutex;
	std::function<void()> on_task_finished;

	template <typename Func, typename OnFinished>
	void run_async(std::string msg, Func background_task, OnFinished main_thread_update);

	void process_pending_tasks();
	void run_hw2_tasks();
	void update_visualization();
};

#endif