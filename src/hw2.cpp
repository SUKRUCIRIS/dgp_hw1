#include "hw2.h"
#include "../third_party/Easy3D/easy3d/util/initializer.h"
#include "../third_party/Easy3D/easy3d/fileio/surface_mesh_io.h"
#include "../third_party/Easy3D/3rd_party/imgui/imgui.h"
#include "../third_party/Easy3D/3rd_party/imgui/backends/imgui_impl_glfw.h"
#include "../third_party/Easy3D/3rd_party/imgui/backends/imgui_impl_opengl3.h"
#include "../third_party/Easy3D/3rd_party/glfw/include/GLFW/glfw3.h"
#include "../third_party/Easy3D/3rd_party/portable_file_dialogs/portable_file_dialogs.h"
#include <iostream>
#include <queue>
#include <limits>
#include <chrono>
#include <thread>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <random>

using namespace easy3d;

easy3d::vec3 closest_point_on_triangle(const easy3d::vec3 &p, const easy3d::vec3 &a, const easy3d::vec3 &b, const easy3d::vec3 &c)
{
	easy3d::vec3 ab = b - a;
	easy3d::vec3 ac = c - a;
	easy3d::vec3 ap = p - a;
	float d1 = easy3d::dot(ab, ap);
	float d2 = easy3d::dot(ac, ap);
	if (d1 <= 0.0f && d2 <= 0.0f)
		return a;

	easy3d::vec3 bp = p - b;
	float d3 = easy3d::dot(ab, bp);
	float d4 = easy3d::dot(ac, bp);
	if (d3 >= 0.0f && d4 <= d3)
		return b;

	float vc = d1 * d4 - d3 * d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
	{
		float v = d1 / (d1 - d3);
		return a + v * ab;
	}

	easy3d::vec3 cp = p - c;
	float d5 = easy3d::dot(ab, cp);
	float d6 = easy3d::dot(ac, cp);
	if (d6 >= 0.0f && d5 <= d6)
		return c;

	float vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
	{
		float w = d2 / (d2 - d6);
		return a + w * ac;
	}

	float va = d3 * d6 - d5 * d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
	{
		float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		return b + w * (c - b);
	}

	float denom = va + vb + vc;
	if (denom <= 1e-8f)
	{
		return (a + b + c) / 3.0f;
	}
	float v = vb / denom;
	float w = vc / denom;
	return a + ab * v + ac * w;
}

Hw2Viewer::Hw2Viewer(const std::string &title) : Viewer(title)
{
	points_drawable = nullptr;
}

Hw2Viewer::~Hw2Viewer()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}

void Hw2Viewer::init()
{
	Viewer::init();

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO &io = ImGui::GetIO();
	(void)io;
	ImGui::StyleColorsDark();

	ImGui_ImplGlfw_InitForOpenGL(window_, true);
	ImGui_ImplOpenGL3_Init("#version 330");
}

template <typename Func, typename OnFinished>
void Hw2Viewer::run_async(std::string msg, Func background_task, OnFinished main_thread_update)
{
	is_processing = true;
	processing_msg = msg;

	std::cout << "\n================================================================================\n";
	std::cout << "[PROCESS START] Initiating Routine: " << msg << "\n";
	std::cout << "================================================================================\n";

	std::thread([this, bg = background_task, fg = main_thread_update]()
				{
        auto start_time = std::chrono::high_resolution_clock::now();
        
        bg();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_time - start_time;

        std::lock_guard<std::mutex> lock(this->task_mutex);
        this->on_task_finished = fg;
        this->task_finished = true;

        std::cout << "================================================================================\n";
        std::cout << "[PROCESS END]   Routine '" << this->processing_msg << "' successfully completed in " << diff.count() << " seconds\n";
        std::cout << "================================================================================\n\n"; })
		.detach();
}

void Hw2Viewer::process_pending_tasks()
{
	if (task_finished.load())
	{
		std::lock_guard<std::mutex> lock(this->task_mutex);
		if (on_task_finished)
		{
			auto ui_start = std::chrono::high_resolution_clock::now();
			on_task_finished();
			auto ui_end = std::chrono::high_resolution_clock::now();
			std::cout << "[GUI THREAD] Main thread UI/Visualization update took "
					  << std::chrono::duration<double>(ui_end - ui_start).count() << " seconds.\n";

			on_task_finished = nullptr;
		}
		task_finished = false;
		is_processing = false;
	}
}

void Hw2Viewer::run_hw2_tasks()
{
	if (!current_mesh || current_mesh->n_faces() == 0)
		return;

	run_async("Computing Visibility Parity & Clustering", [this]()
			  {
				  auto t_start = std::chrono::high_resolution_clock::now();

				  easy3d::vec3 min_pt(1e9f, 1e9f, 1e9f), max_pt(-1e9f, -1e9f, -1e9f);
				  for (auto v : current_mesh->vertices())
				  {
					  easy3d::vec3 p = current_mesh->position(v);
					  min_pt.x = std::min(min_pt.x, p.x);
					  min_pt.y = std::min(min_pt.y, p.y);
					  min_pt.z = std::min(min_pt.z, p.z);
					  max_pt.x = std::max(max_pt.x, p.x);
					  max_pt.y = std::max(max_pt.y, p.y);
					  max_pt.z = std::max(max_pt.z, p.z);
				  }

				  float diag = easy3d::length(max_pt - min_pt);
				  float actual_discard = (discard_dist / 100.0f) * diag;

				  easy3d::vec3 step = (max_pt - min_pt) / std::max(1.0f, (float)(grid_res - 1));

				  std::vector<easy3d::vec3> face_centroids(current_mesh->n_faces());
				  std::vector<easy3d::vec3> face_normals(current_mesh->n_faces());
				  std::vector<std::vector<easy3d::vec3>> face_pts(current_mesh->n_faces(), std::vector<easy3d::vec3>(3));

				  int f_idx = 0;
                  for (auto f : current_mesh->faces())
                  {
                      easy3d::vec3 c(0, 0, 0);
                      int k = 0;
                      for (auto v : current_mesh->vertices(f))
                      {
                          c += current_mesh->position(v);
                          if (k < 3) {
                              face_pts[f_idx][k] = current_mesh->position(v);
                          }
                          k++;
                      }
                      
                      c /= (float)k; 
                      face_centroids[f_idx] = c;

                      easy3d::vec3 p0 = face_pts[f_idx][0];
                      easy3d::vec3 p1 = face_pts[f_idx][1];
                      easy3d::vec3 p2 = face_pts[f_idx][2];
                      easy3d::vec3 cr = cross(p1 - p0, p2 - p0);
                      float len = length(cr);
                      face_normals[f_idx] = len > 1e-8f ? cr / len : easy3d::vec3(0, 1, 0);
                      f_idx++;
                  }

				  auto t_bbox_normals = std::chrono::high_resolution_clock::now();
				  std::cout << "  -> [SUB-TASK 1] Computed bounding box and pre-calculated face normals/centroids in "
                            << std::chrono::duration<double>(t_bbox_normals - t_start).count() << " seconds.\n";

				  std::vector<GridSample> temp_samples;
				  int total_evaluated = 0;
				  for (int i = 0; i < grid_res; i++)
				  {
					  for (int j = 0; j < grid_res; j++)
					  {
						  for (int k = 0; k < grid_res; k++)
						  {
							  total_evaluated++;
							  easy3d::vec3 g = min_pt + easy3d::vec3(i * step.x, j * step.y, k * step.z);

							  float min_d = 1e9f;
							  int closest_f = -1;
							  easy3d::vec3 closest_pt(0, 0, 0);

							  for (int fi = 0; fi < current_mesh->n_faces(); fi++)
							  {
								  easy3d::vec3 cp = closest_point_on_triangle(g, face_pts[fi][0], face_pts[fi][1], face_pts[fi][2]);
								  float d = easy3d::length(g - cp);
								  if (d < min_d)
								  {
									  min_d = d;
									  closest_f = fi;
									  closest_pt = cp;
								  }
							  }

							  if (min_d < actual_discard)
								  continue;

							  easy3d::vec3 n = face_normals[closest_f];

							  if (dot(g - closest_pt, n) > 0)
								  continue;

							  int parity = 0;
							  for (int fi = 0; fi < current_mesh->n_faces(); fi++)
							  {
								  if (dot(g - face_centroids[fi], face_normals[fi]) < 0)
								  {
									  parity++;
								  }
							  }

							  temp_samples.push_back({g, i, j, k, parity, -1});
						  }
					  }
				  }

				  auto t_parity = std::chrono::high_resolution_clock::now();
				  std::cout << "  -> [SUB-TASK 2] Evaluated " << total_evaluated << " grid points.\n";
				  std::cout << "  -> [SUB-TASK 2] Completed inside/outside testing and parity scoring in "
                            << std::chrono::duration<double>(t_parity - t_bbox_normals).count() << " seconds.\n";
                  std::cout << "  -> [SUB-TASK 2] Filtered down to " << temp_samples.size() << " valid interior samples.\n";

				  std::vector<GridSample *> grid_ptrs(grid_res * grid_res * grid_res, nullptr);
				  for (auto &s : temp_samples)
				  {
					  grid_ptrs[s.ix + s.iy * grid_res + s.iz * grid_res * grid_res] = &s;
				  }

				  std::vector<GridSample *> sorted_ptrs;
				  for (auto &s : temp_samples)
					  sorted_ptrs.push_back(&s);

				  std::string heuristic_name = "";
				  if (bfs_heuristic == 0) {
                      heuristic_name = "Highest Parity First";
					  std::sort(sorted_ptrs.begin(), sorted_ptrs.end(), [](GridSample *A, GridSample *B)
								{ return A->parity > B->parity; });
				  } else if (bfs_heuristic == 1) {
                      heuristic_name = "Lowest Parity First";
					  std::sort(sorted_ptrs.begin(), sorted_ptrs.end(), [](GridSample *A, GridSample *B)
								{ return A->parity < B->parity; });
				  } else if (bfs_heuristic == 2) {
                      heuristic_name = "Random Seed";
					  auto rng = std::default_random_engine { std::random_device{}() };
					  std::shuffle(sorted_ptrs.begin(), sorted_ptrs.end(), rng);
				  } else if (bfs_heuristic == 3) {
                      heuristic_name = "Spatial (Lexicographical X, Y, Z)";
					  std::sort(sorted_ptrs.begin(), sorted_ptrs.end(), [](GridSample *A, GridSample *B) {
						  if (A->ix != B->ix) return A->ix < B->ix;
						  if (A->iy != B->iy) return A->iy < B->iy;
						  return A->iz < B->iz;
					  });
				  }

                  std::cout << "  -> [SUB-TASK 3] Starting BFS traversal with heuristic: " << heuristic_name << ".\n";

				  int c_id = 0;
				  for (auto s_ptr : sorted_ptrs)
				  {
					  if (s_ptr->cluster != -1)
						  continue;

					  std::queue<GridSample *> q;
					  q.push(s_ptr);
					  s_ptr->cluster = c_id;

					  while (!q.empty())
					  {
						  GridSample *curr = q.front();
						  q.pop();

						  for (int dx = -1; dx <= 1; dx++)
						  {
							  for (int dy = -1; dy <= 1; dy++)
							  {
								  for (int dz = -1; dz <= 1; dz++)
								  {
									  if (dx == 0 && dy == 0 && dz == 0)
										  continue;

									  int nx = curr->ix + dx;
									  int ny = curr->iy + dy;
									  int nz = curr->iz + dz;

									  if (nx >= 0 && nx < grid_res && ny >= 0 && ny < grid_res && nz >= 0 && nz < grid_res)
									  {
										  GridSample *n_ptr = grid_ptrs[nx + ny * grid_res + nz * grid_res * grid_res];
										  if (n_ptr && n_ptr->cluster == -1 && n_ptr->parity == curr->parity)
										  {
											  n_ptr->cluster = c_id;
											  q.push(n_ptr);
										  }
									  }
								  }
							  }
						  }
					  }
					  c_id++;
				  }

				  auto t_cluster = std::chrono::high_resolution_clock::now();
				  std::cout << "  -> [SUB-TASK 3] BFS grouping isolated " << c_id << " distinct clusters in " 
                            << std::chrono::duration<double>(t_cluster - t_parity).count() << " seconds.\n";

				  std::lock_guard<std::mutex> lock(this->task_mutex);
				  this->valid_samples = temp_samples;
				  this->num_clusters = c_id;

				  int current_max_parity = 0;
				  for (auto s_ptr : sorted_ptrs) {
                      if (s_ptr->parity > current_max_parity) {
                          current_max_parity = s_ptr->parity;
                      }
				  }
				  this->max_parity = current_max_parity;

				  this->cluster_colors.resize(c_id);
				  for (int i = 0; i < c_id; i++)
				  {
					  float hue = (i * 0.618033988749895f);
					  hue -= std::floor(hue);
					  float s = 0.8f, v = 0.9f;
					  int hi = (int)(hue * 6.0f);
					  float f = hue * 6.0f - hi;
					  float p = v * (1.0f - s);
					  float q = v * (1.0f - f * s);
					  float t = v * (1.0f - (1.0f - f) * s);

					  if (hi == 0)
						  this->cluster_colors[i] = easy3d::vec3(v, t, p);
					  else if (hi == 1)
						  this->cluster_colors[i] = easy3d::vec3(q, v, p);
					  else if (hi == 2)
						  this->cluster_colors[i] = easy3d::vec3(p, v, t);
					  else if (hi == 3)
						  this->cluster_colors[i] = easy3d::vec3(p, q, v);
					  else if (hi == 4)
						  this->cluster_colors[i] = easy3d::vec3(t, p, v);
					  else
						  this->cluster_colors[i] = easy3d::vec3(v, p, q);
				  } }, [this]()
			  { this->update_visualization(); });
}

void Hw2Viewer::update_visualization()
{
	if (!current_mesh)
		return;

	if (this->points_drawable)
	{
		this->delete_drawable(this->points_drawable);
		this->points_drawable = nullptr;
	}

	this->clear_scene();
	this->add_model(current_mesh);

	auto v_color = current_mesh->get_vertex_property<easy3d::vec3>("v:color");
	if (v_color)
		current_mesh->remove_vertex_property(v_color);

	auto f_color = current_mesh->get_face_property<easy3d::vec3>("f:color");
	if (!f_color)
		f_color = current_mesh->add_face_property<easy3d::vec3>("f:color");

	for (auto f : current_mesh->faces())
	{
		easy3d::vec3 centroid(0, 0, 0);
		int v_count = 0;
		for (auto v : current_mesh->vertices(f))
		{
			centroid += current_mesh->position(v);
			v_count++;
		}
		centroid /= (float)v_count;

		auto h = current_mesh->halfedge(f);
		easy3d::vec3 p0 = current_mesh->position(current_mesh->source(h));
		easy3d::vec3 p1 = current_mesh->position(current_mesh->target(h));
		easy3d::vec3 p2 = current_mesh->position(current_mesh->target(current_mesh->next(h)));
		easy3d::vec3 cr = easy3d::cross(p1 - p0, p2 - p0);
		float len = easy3d::length(cr);
		easy3d::vec3 normal = len > 1e-8f ? cr / len : easy3d::vec3(0, 1, 0);

		int best_c = -1;
		float best_d = 1e9f;

		for (auto &s : valid_samples)
		{
			if (easy3d::dot(s.p - centroid, normal) > 0.0f)
				continue;

			float d = easy3d::length(s.p - centroid);
			if (d < best_d)
			{
				best_d = d;
				best_c = s.cluster;
			}
		}

		if (best_c == -1)
		{
			for (auto &s : valid_samples)
			{
				float d = easy3d::length(s.p - centroid);
				if (d < best_d)
				{
					best_d = d;
					best_c = s.cluster;
				}
			}
		}

		if (selected_cluster != -1 && best_c != selected_cluster)
		{
			f_color[f] = easy3d::vec3(0.2f, 0.2f, 0.2f);
		}
		else
		{
			if (best_c != -1)
			{
				f_color[f] = cluster_colors[best_c];
			}
			else
			{
				f_color[f] = easy3d::vec3(0.5f, 0.5f, 0.5f);
			}
		}
	}

	std::vector<easy3d::vec3> pts;
	std::vector<easy3d::vec3> cols;

	for (auto &s : valid_samples)
	{
		if (show_highest_parity && s.parity != max_parity)
			continue;
		if (selected_cluster != -1 && s.cluster != selected_cluster)
			continue;

		pts.push_back(s.p);
		cols.push_back(cluster_colors[s.cluster]);
	}

	if (!pts.empty())
	{
		points_drawable = new easy3d::PointsDrawable("grid_points", current_mesh.get());
		points_drawable->update_vertex_buffer(pts);
		points_drawable->update_color_buffer(cols);
		points_drawable->set_point_size(5.0f);
		this->add_drawable(points_drawable);
	}

	this->update();
}

void Hw2Viewer::post_draw()
{
	Viewer::post_draw();
	process_pending_tasks();

	bool disabled = is_processing.load();
	bool mesh_unloaded = (current_mesh == nullptr);

	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	ImGui::Begin("Homework 2 Toolbox");

	if (disabled)
	{
		ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Working: %s...", processing_msg.c_str());
		ImGui::Separator();
	}

	ImGui::BeginDisabled(disabled);

	ImGui::Text("Load Mesh");
	ImGui::PushItemWidth(200);
	ImGui::InputText("##filepath", filepath, 512);
	ImGui::PopItemWidth();
	ImGui::SameLine();
	if (ImGui::Button("Browse..."))
	{
		auto sel = pfd::open_file("Select Mesh", ".", {"Mesh Files", "*.obj *.off *.ply *.stl *.sm", "All Files", "*"}).result();
		if (!sel.empty())
		{
			strncpy(filepath, sel[0].c_str(), 511);
			filepath[511] = '\0';
		}
	}
	ImGui::SameLine();
	if (ImGui::Button("Load"))
	{
		std::string fp(filepath);
		run_async("Loading Mesh Data", [this, fp]()
				  {
            auto m = SurfaceMeshIO::load(fp);
            std::lock_guard<std::mutex> lock(this->task_mutex);
            this->current_mesh = std::shared_ptr<easy3d::SurfaceMesh>(m); }, [this]()
				  {
            if (this->current_mesh) {
                this->valid_samples.clear();
                this->cluster_colors.clear();
                this->num_clusters = 0;
                this->selected_cluster = -1;
                this->show_highest_parity = false;

                if (this->points_drawable) {
                    this->delete_drawable(this->points_drawable);
                    this->points_drawable = nullptr;
                }

                this->clear_scene();
                this->add_model(this->current_mesh);
                this->fit_screen(); 
                this->update();
            } });
	}

	ImGui::EndDisabled();

	ImGui::BeginDisabled(disabled || mesh_unloaded);

	ImGui::Separator();
	ImGui::InputInt("Grid Resolution", &grid_res);
	if (grid_res < 5)
		grid_res = 5;
	if (grid_res > 100)
		grid_res = 100;

	ImGui::InputFloat("Discard Dist (%)", &discard_dist, 0.1f, 1.0f, "%.1f");
	if (discard_dist < 0.0f)
		discard_dist = 0.0f;
	if (discard_dist > 20.0f)
		discard_dist = 20.0f;

	ImGui::Combo("BFS Seed Heuristic", &bfs_heuristic, "Highest Parity\0Lowest Parity\0Random\0Spatial (X,Y,Z)\0");

	if (ImGui::Button("Compute Parity & Clusters"))
	{
		selected_cluster = -1;
		run_hw2_tasks();
	}

	ImGui::EndDisabled();

	ImGui::BeginDisabled(disabled || mesh_unloaded || valid_samples.empty());

	ImGui::Separator();
	ImGui::Text("Visualization Options");

	if (ImGui::Checkbox("Show Highest Parity Only", &show_highest_parity))
	{
		update_visualization();
	}

	ImGui::Text("Filter Cluster (-1 for all):");
	if (ImGui::SliderInt("##cluster_slider", &selected_cluster, -1, num_clusters - 1))
	{
		update_visualization();
	}

	ImGui::EndDisabled();

	ImGui::End();

	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

	if (disabled)
	{
		this->update();
	}
}