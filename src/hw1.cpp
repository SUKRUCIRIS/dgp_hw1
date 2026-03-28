// hw1.cpp
#include "hw1.h"

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

using namespace easy3d;

float point_triangle_distance_hw1(const vec3 &p, const vec3 &a, const vec3 &b, const vec3 &c)
{
	vec3 ab = b - a;
	vec3 ac = c - a;
	vec3 ap = p - a;
	float d1 = dot(ab, ap);
	float d2 = dot(ac, ap);
	if (d1 <= 0.0f && d2 <= 0.0f)
		return length(ap);

	vec3 bp = p - b;
	float d3 = dot(ab, bp);
	float d4 = dot(ac, bp);
	if (d3 >= 0.0f && d4 <= d3)
		return length(bp);

	float vc = d1 * d4 - d3 * d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
	{
		float denom = d1 - d3;
		float v = (denom == 0.0f) ? 0.0f : d1 / denom;
		return length(a + v * ab - p);
	}

	vec3 cp = p - c;
	float d5 = dot(ab, cp);
	float d6 = dot(ac, cp);
	if (d6 >= 0.0f && d5 <= d6)
		return length(cp);

	float vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
	{
		float denom = d2 - d6;
		float w = (denom == 0.0f) ? 0.0f : d2 / denom;
		return length(a + w * ac - p);
	}

	float va = d3 * d6 - d5 * d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
	{
		float denom = (d4 - d3) + (d5 - d6);
		float w = (denom == 0.0f) ? 0.0f : (d4 - d3) / denom;
		return length(b + w * (c - b) - p);
	}

	float denom = va + vb + vc;
	if (denom == 0.0f)
	{
		return std::min({length(ap), length(bp), length(cp)});
	}
	float v = vb / denom;
	float w = vc / denom;
	return length(a + ab * v + ac * w - p);
}

void dijkstra_array(SurfaceMesh *mesh, SurfaceMesh::Vertex start, std::vector<float> &dist, std::vector<SurfaceMesh::Vertex> &parent, const SurfaceMesh::EdgeProperty<float> &edge_lengths)
{
	int n = mesh->n_vertices();
	dist.assign(n, std::numeric_limits<float>::infinity());
	parent.assign(n, SurfaceMesh::Vertex());
	std::vector<bool> visited(n, false);

	dist[start.idx()] = 0.0f;

	for (int i = 0; i < n; ++i)
	{
		float min_d = std::numeric_limits<float>::infinity();
		int u_idx = -1;
		for (int j = 0; j < n; ++j)
		{
			if (!visited[j] && dist[j] < min_d)
			{
				min_d = dist[j];
				u_idx = j;
			}
		}

		if (u_idx == -1)
			break;

		visited[u_idx] = true;
		SurfaceMesh::Vertex u(u_idx);

		for (auto h : mesh->halfedges(u))
		{
			SurfaceMesh::Vertex v = mesh->target(h);
			float weight = edge_lengths[mesh->edge(h)];

			if (!visited[v.idx()] && dist[u.idx()] + weight < dist[v.idx()])
			{
				dist[v.idx()] = dist[u.idx()] + weight;
				parent[v.idx()] = u;
			}
		}
	}
}

void dijkstra_min_heap(SurfaceMesh *mesh, SurfaceMesh::Vertex start, std::vector<float> &dist, std::vector<SurfaceMesh::Vertex> &parent, const SurfaceMesh::EdgeProperty<float> &edge_lengths)
{
	int n = mesh->n_vertices();
	dist.assign(n, std::numeric_limits<float>::infinity());
	parent.assign(n, SurfaceMesh::Vertex());

	using Node = std::pair<float, int>;
	std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

	dist[start.idx()] = 0.0f;
	pq.push({0.0f, start.idx()});

	while (!pq.empty())
	{
		auto [d, u_idx] = pq.top();
		pq.pop();

		if (d > dist[u_idx])
			continue;

		SurfaceMesh::Vertex u(u_idx);

		for (auto h : mesh->halfedges(u))
		{
			SurfaceMesh::Vertex v = mesh->target(h);
			float weight = edge_lengths[mesh->edge(h)];

			if (dist[u.idx()] + weight < dist[v.idx()])
			{
				dist[v.idx()] = dist[u.idx()] + weight;
				parent[v.idx()] = u;
				pq.push({dist[v.idx()], v.idx()});
			}
		}
	}
}

void compute_and_dump_all_pairs(SurfaceMesh *mesh, const SurfaceMesh::EdgeProperty<float> &edge_lengths)
{
	int n = mesh->n_vertices();
	std::vector<std::vector<float>> dist_matrix(n, std::vector<float>(n, 0.0f));
	int num_threads = std::max(1u, std::thread::hardware_concurrency());
	std::vector<std::thread> threads;
	std::atomic<int> current_idx(0);

	for (int t = 0; t < num_threads; ++t)
	{
		threads.emplace_back([&]()
							 {
            std::vector<SurfaceMesh::Vertex> local_parent(n);
            while (true) {
                int i = current_idx.fetch_add(1);
                if (i >= n) break;
                SurfaceMesh::Vertex v(i);
                dijkstra_min_heap(mesh, v, dist_matrix[i], local_parent, edge_lengths);
            } });
	}

	for (auto &t : threads)
	{
		if (t.joinable())
		{
			t.join();
		}
	}

	FILE *outfile = fopen("geodesic_matrix.txt", "w");
	if (outfile)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				fprintf(outfile, "%.6f ", dist_matrix[i][j]);
			}
			fprintf(outfile, "\n");
		}
		fclose(outfile);
		std::cout << "Successfully wrote geodesic_matrix.txt to disk.\n";
	}
	else
	{
		std::cerr << "Error: Could not open geodesic_matrix.txt for writing.\n";
	}
}

SurfaceMesh *subdivide_4to1(SurfaceMesh *mesh)
{
	SurfaceMesh *out = new SurfaceMesh();
	auto v_map = mesh->add_vertex_property<SurfaceMesh::Vertex>("v:new_v");
	auto e_map = mesh->add_edge_property<SurfaceMesh::Vertex>("e:new_v");

	for (auto v : mesh->vertices())
	{
		if (mesh->is_border(v))
		{
			vec3 pos = mesh->position(v) * 0.75f;
			for (auto h : mesh->halfedges(v))
			{
				if (mesh->is_border(mesh->edge(h)))
				{
					pos += mesh->position(mesh->target(h)) * 0.125f;
				}
			}
			v_map[v] = out->add_vertex(pos);
		}
		else
		{
			int n = 0;
			vec3 neighbor_sum(0, 0, 0);
			for (auto h : mesh->halfedges(v))
			{
				neighbor_sum += mesh->position(mesh->target(h));
				n++;
			}

			float beta = 1.0f / n * (5.0f / 8.0f - std::pow(3.0f / 8.0f + 1.0f / 4.0f * std::cos(2.0f * M_PI / n), 2.0f));
			vec3 new_pos = mesh->position(v) * (1.0f - n * beta) + neighbor_sum * beta;
			v_map[v] = out->add_vertex(new_pos);
		}
	}

	for (auto e : mesh->edges())
	{
		SurfaceMesh::Halfedge h0 = mesh->halfedge(e, 0);
		vec3 p0 = mesh->position(mesh->source(h0));
		vec3 p1 = mesh->position(mesh->target(h0));

		if (mesh->is_border(e))
		{
			e_map[e] = out->add_vertex((p0 + p1) * 0.5f);
		}
		else
		{
			SurfaceMesh::Halfedge h1 = mesh->halfedge(e, 1);
			vec3 p2 = mesh->position(mesh->target(mesh->next(h0)));
			vec3 p3 = mesh->position(mesh->target(mesh->next(h1)));

			vec3 new_pos = 0.375f * (p0 + p1) + 0.125f * (p2 + p3);
			e_map[e] = out->add_vertex(new_pos);
		}
	}

	for (auto f : mesh->faces())
	{
		SurfaceMesh::Halfedge h0 = mesh->halfedge(f);
		SurfaceMesh::Halfedge h1 = mesh->next(h0);
		SurfaceMesh::Halfedge h2 = mesh->next(h1);

		SurfaceMesh::Vertex v0 = v_map[mesh->source(h0)];
		SurfaceMesh::Vertex v1 = v_map[mesh->source(h1)];
		SurfaceMesh::Vertex v2 = v_map[mesh->source(h2)];

		SurfaceMesh::Vertex m0 = e_map[mesh->edge(h0)];
		SurfaceMesh::Vertex m1 = e_map[mesh->edge(h1)];
		SurfaceMesh::Vertex m2 = e_map[mesh->edge(h2)];

		out->add_triangle(v0, m0, m2);
		out->add_triangle(m0, v1, m1);
		out->add_triangle(m2, m1, v2);
		out->add_triangle(m0, m1, m2);
	}

	mesh->remove_vertex_property(v_map);
	mesh->remove_edge_property(e_map);
	return out;
}

SurfaceMesh *subdivide_phong(SurfaceMesh *mesh)
{
	auto v_normal = mesh->add_vertex_property<vec3>("v:phong_normal", vec3(0, 0, 0));
	for (auto v : mesh->vertices())
		v_normal[v] = vec3(0, 0, 0);

	for (auto f : mesh->faces())
	{
		SurfaceMesh::Halfedge h0 = mesh->halfedge(f);
		vec3 p0 = mesh->position(mesh->source(h0));
		vec3 p1 = mesh->position(mesh->target(h0));
		vec3 p2 = mesh->position(mesh->target(mesh->next(h0)));

		vec3 cr = cross(p1 - p0, p2 - p0);
		float len = length(cr);
		vec3 n = len > 1e-8f ? cr / len : vec3(0, 0, 0);

		for (auto v : mesh->vertices(f))
			v_normal[v] += n;
	}

	for (auto v : mesh->vertices())
	{
		float len = length(v_normal[v]);
		v_normal[v] = len > 1e-8f ? v_normal[v] / len : vec3(0, 1, 0);
	}

	SurfaceMesh *out = new SurfaceMesh();
	auto v_map = mesh->add_vertex_property<SurfaceMesh::Vertex>("v:new_v");
	auto e_map = mesh->add_edge_property<SurfaceMesh::Vertex>("e:new_v");

	for (auto v : mesh->vertices())
	{
		v_map[v] = out->add_vertex(mesh->position(v));
	}

	for (auto e : mesh->edges())
	{
		SurfaceMesh::Halfedge he0 = mesh->halfedge(e, 0);
		SurfaceMesh::Vertex v0 = mesh->source(he0);
		SurfaceMesh::Vertex v1 = mesh->target(he0);

		vec3 p0 = mesh->position(v0);
		vec3 p1 = mesh->position(v1);
		vec3 n0 = v_normal[v0];
		vec3 n1 = v_normal[v1];

		vec3 mid = (p0 + p1) * 0.5f;
		vec3 proj0 = mid - dot(mid - p0, n0) * n0;
		vec3 proj1 = mid - dot(mid - p1, n1) * n1;
		vec3 phong_mid = (proj0 + proj1) * 0.5f;

		e_map[e] = out->add_vertex(phong_mid);
	}

	for (auto f : mesh->faces())
	{
		SurfaceMesh::Halfedge h0 = mesh->halfedge(f);
		SurfaceMesh::Halfedge h1 = mesh->next(h0);
		SurfaceMesh::Halfedge h2 = mesh->next(h1);

		SurfaceMesh::Vertex v0 = v_map[mesh->source(h0)];
		SurfaceMesh::Vertex v1 = v_map[mesh->source(h1)];
		SurfaceMesh::Vertex v2 = v_map[mesh->source(h2)];

		SurfaceMesh::Vertex m0 = e_map[mesh->edge(h0)];
		SurfaceMesh::Vertex m1 = e_map[mesh->edge(h1)];
		SurfaceMesh::Vertex m2 = e_map[mesh->edge(h2)];

		out->add_triangle(v0, m0, m2);
		out->add_triangle(m0, v1, m1);
		out->add_triangle(m2, m1, v2);
		out->add_triangle(m0, m1, m2);
	}

	mesh->remove_vertex_property(v_normal);
	mesh->remove_vertex_property(v_map);
	mesh->remove_edge_property(e_map);
	return out;
}

SurfaceMesh *subdivide_sqrt3(SurfaceMesh *mesh)
{
	SurfaceMesh *out = new SurfaceMesh();
	auto v_map = mesh->add_vertex_property<SurfaceMesh::Vertex>("v:new_v");
	auto f_map = mesh->add_face_property<SurfaceMesh::Vertex>("f:new_v");

	for (auto v : mesh->vertices())
	{
		if (mesh->is_border(v))
		{
			v_map[v] = out->add_vertex(mesh->position(v));
		}
		else
		{
			int n = 0;
			vec3 neighbor_sum(0, 0, 0);
			for (auto h : mesh->halfedges(v))
			{
				neighbor_sum += mesh->position(mesh->target(h));
				n++;
			}
			float alpha = (4.0f - 2.0f * std::cos(2.0f * M_PI / n)) / 9.0f;
			vec3 new_pos = (1.0f - alpha) * mesh->position(v) + (alpha / n) * neighbor_sum;
			v_map[v] = out->add_vertex(new_pos);
		}
	}

	for (auto f : mesh->faces())
	{
		vec3 centroid(0, 0, 0);
		int n = 0;
		for (auto v : mesh->vertices(f))
		{
			centroid += mesh->position(v);
			n++;
		}
		centroid /= n;
		f_map[f] = out->add_vertex(centroid);
	}

	for (auto e : mesh->edges())
	{
		if (mesh->is_border(e))
			continue;

		SurfaceMesh::Halfedge h0 = mesh->halfedge(e, 0);
		SurfaceMesh::Halfedge h1 = mesh->halfedge(e, 1);

		SurfaceMesh::Vertex c0 = f_map[mesh->face(h0)];
		SurfaceMesh::Vertex c1 = f_map[mesh->face(h1)];

		SurfaceMesh::Vertex v_end = v_map[mesh->target(h0)];
		SurfaceMesh::Vertex v_start = v_map[mesh->target(h1)];

		out->add_triangle(v_start, c1, c0);
		out->add_triangle(v_end, c0, c1);
	}

	for (auto h : mesh->halfedges())
	{
		if (mesh->is_border(h))
		{
			SurfaceMesh::Halfedge opp = mesh->opposite(h);
			SurfaceMesh::Face f = mesh->face(opp);
			SurfaceMesh::Vertex c = f_map[f];
			SurfaceMesh::Vertex v0 = v_map[mesh->source(h)];
			SurfaceMesh::Vertex v1 = v_map[mesh->target(h)];
			out->add_triangle(v1, v0, c);
		}
	}

	mesh->remove_vertex_property(v_map);
	mesh->remove_face_property(f_map);
	return out;
}

SurfaceMesh *run_subdiv(SurfaceMesh *m, int type, int iters)
{
	if (iters == 0)
		return nullptr;
	SurfaceMesh *curr = m;
	for (int i = 0; i < iters; ++i)
	{
		SurfaceMesh *next = nullptr;
		if (type == 0)
			next = subdivide_sqrt3(curr);
		else if (type == 1)
			next = subdivide_4to1(curr);
		else if (type == 2)
			next = subdivide_phong(curr);

		if (curr != m)
			delete curr;
		curr = next;
	}
	return curr;
}

void translate_mesh(SurfaceMesh *mesh, const vec3 &offset)
{
	if (!mesh)
		return;
	for (auto v : mesh->vertices())
	{
		mesh->position(v) += offset;
	}
}

void print_mesh_stats(SurfaceMesh *mesh, const std::string &name)
{
	if (!mesh)
		return;
	int faces = mesh->n_faces();
	float area = 0.0f;
	for (auto f : mesh->faces())
	{
		auto h = mesh->halfedge(f);
		vec3 p0 = mesh->position(mesh->source(h));

		auto h_curr = mesh->next(h);
		while (mesh->target(h_curr) != mesh->source(h))
		{
			vec3 p1 = mesh->position(mesh->source(h_curr));
			vec3 p2 = mesh->position(mesh->target(h_curr));
			area += 0.5f * length(cross(p1 - p0, p2 - p0));
			h_curr = mesh->next(h_curr);
		}
	}
	std::cout << name << " - Triangles/Polygons: " << faces << " | Total Area: " << area << "\n";
}

Hw1Viewer::Hw1Viewer(const std::string &title) : Viewer(title)
{
	geodesic_drawable = nullptr;
}

Hw1Viewer::~Hw1Viewer()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}

void Hw1Viewer::init()
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
void Hw1Viewer::run_async(std::string msg, Func background_task, OnFinished main_thread_update)
{
	is_processing = true;
	processing_msg = msg;

	std::cout << "\n========================================\n";
	std::cout << "[START] " << msg << "...\n";
	std::cout << "========================================\n";

	std::thread([this, bg = background_task, fg = main_thread_update]()
				{
        auto start_time = std::chrono::high_resolution_clock::now();
        
        bg();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_time - start_time;

        std::lock_guard<std::mutex> lock(this->task_mutex);
        this->on_task_finished = fg;
        this->task_finished = true;

        std::cout << "========================================\n";
        std::cout << "[DONE]  " << this->processing_msg << " finished in " << diff.count() << "s\n";
        std::cout << "========================================\n\n"; })
		.detach();
}

void Hw1Viewer::process_pending_tasks()
{
	if (task_finished.load())
	{
		std::lock_guard<std::mutex> lock(this->task_mutex);
		if (on_task_finished)
		{
			on_task_finished();
			on_task_finished = nullptr;
		}
		task_finished = false;
		is_processing = false;
	}
}

void Hw1Viewer::visualize_path(SurfaceMesh::Vertex start, SurfaceMesh::Vertex end)
{
	if (geodesic_drawable)
	{
		this->delete_drawable(geodesic_drawable);
		geodesic_drawable = nullptr;
	}

	if (this->current_meshes.empty())
		return;
	SurfaceMesh *render_mesh = this->current_meshes[0];

	std::vector<vec3> path_points;
	std::vector<unsigned int> path_indices;

	SurfaceMesh::Vertex curr = end;
	while (curr.is_valid())
	{
		path_points.push_back(render_mesh->position(curr));
		if (curr == start)
			break;
		curr = temp_parent[curr.idx()];
	}

	for (size_t i = 0; i + 1 < path_points.size(); ++i)
	{
		path_indices.push_back(i);
		path_indices.push_back(i + 1);
	}

	if (!path_indices.empty())
	{
		auto line_drawable = new LinesDrawable("geodesic_path", render_mesh);
		line_drawable->update_vertex_buffer(path_points);
		line_drawable->update_element_buffer(path_indices);
		line_drawable->set_impostor_type(LinesDrawable::CYLINDER);
		line_drawable->set_uniform_coloring(vec4(1.0f, 0.0f, 0.0f, 1.0f));
		line_drawable->set_line_width(3.0f);

		geodesic_drawable = line_drawable;
		this->add_drawable(line_drawable);
		this->update();
	}
	else
	{
		std::cout << "Warning: Could not extract valid path between vertices!\n";
	}
}

void Hw1Viewer::post_draw()
{
	Viewer::post_draw();

	process_pending_tasks();

	bool disabled = is_processing.load();
	bool mesh_unloaded = current_meshes.empty();

	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	ImGui::Begin("Homework 1 Toolbox");

	if (disabled)
	{
		ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Working: %s...", processing_msg.c_str());
		ImGui::Separator();
	}

	ImGui::BeginDisabled(disabled);

	ImGui::Text("Input Decimated Mesh (~500 v)");
	ImGui::PushItemWidth(200);
	ImGui::InputText("##filepath", filepath, 512);
	ImGui::PopItemWidth();
	ImGui::SameLine();
	if (ImGui::Button("Browse...##1"))
	{
		auto sel = pfd::open_file("Select Decimated Mesh", ".", {"Mesh Files", "*.obj *.off *.ply *.stl *.sm", "All Files", "*"}).result();
		if (!sel.empty())
		{
			strncpy(filepath, sel[0].c_str(), 511);
			filepath[511] = '\0';
		}
	}
	ImGui::SameLine();
	if (ImGui::Button("Load Mesh"))
	{
		std::string fp(filepath);
		run_async("Loading Decimated Mesh", [this, fp]()
				  { 
            SurfaceMesh* m = SurfaceMeshIO::load(fp); 
            if(m) this->temp_meshes.push_back(m); }, [this]()
				  {
            if (!this->temp_meshes.empty()) {
                this->clear_scene(); 
                geodesic_drawable = nullptr; 

                this->current_meshes = this->temp_meshes;
                this->temp_meshes.clear();

                for (auto m : this->current_meshes) this->add_model(m);
                
                this->fit_screen();
                this->update();
            } else {
                std::cerr << "Failed to load mesh: " << this->filepath << "\n";
            } });
	}

	ImGui::Separator();
	ImGui::Text("Original High-Res Mesh (For Comparison)");
	ImGui::PushItemWidth(200);
	ImGui::InputText("##filepath_orig", filepath_orig, 512);
	ImGui::PopItemWidth();
	ImGui::SameLine();
	if (ImGui::Button("Browse...##2"))
	{
		auto sel = pfd::open_file("Select Original Mesh", ".", {"Mesh Files", "*.obj *.off *.ply *.stl *.sm", "All Files", "*"}).result();
		if (!sel.empty())
		{
			strncpy(filepath_orig, sel[0].c_str(), 511);
			filepath_orig[511] = '\0';
		}
	}

	ImGui::EndDisabled();

	ImGui::BeginDisabled(disabled || mesh_unloaded);

	ImGui::Separator();
	ImGui::Text("Geodesic Distances");
	ImGui::InputInt("Start Vertex", &start_v);
	ImGui::InputInt("End Vertex", &end_v);

	if (ImGui::Button("Run Geodesic") && !current_meshes.empty())
	{
		SurfaceMesh::Vertex sv(start_v);
		SurfaceMesh::Vertex ev(end_v);

		if (sv.is_valid() && ev.is_valid())
		{
			run_async("Computing Geodesic Path & Timings", [this, fp = std::string(filepath), sv, ev]()
					  {
                SurfaceMesh* base_mesh = SurfaceMeshIO::load(fp);
                if (!base_mesh) return;

                if (start_v >= base_mesh->n_vertices() || end_v >= base_mesh->n_vertices()) {
                    delete base_mesh;
                    return;
                }

                auto edge_lengths = base_mesh->add_edge_property<float>("e:length", 0.0f);
                for (auto e : base_mesh->edges()) {
                    SurfaceMesh::Halfedge he = base_mesh->halfedge(e, 0);
                    SurfaceMesh::Vertex v0 = base_mesh->source(he);
                    SurfaceMesh::Vertex v1 = base_mesh->target(he);
                    edge_lengths[e] = length(base_mesh->position(v0) - base_mesh->position(v1));
                }

                auto start_array = std::chrono::high_resolution_clock::now();
                std::vector<float> dist_arr;
                std::vector<SurfaceMesh::Vertex> parent_arr;
                dijkstra_array(base_mesh, sv, dist_arr, parent_arr, edge_lengths);
                auto end_array = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> diff_array = end_array - start_array;

                auto start_heap = std::chrono::high_resolution_clock::now();
                dijkstra_min_heap(base_mesh, sv, this->temp_dist, this->temp_parent, edge_lengths);
                auto end_heap = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> diff_heap = end_heap - start_heap;

                std::cout << "\n----------------------------------------\n";
                std::cout << "Geodesic Algorithm Performance Timings:\n";
                std::cout << "Array priority queue:    " << diff_array.count() << " seconds\n";
                std::cout << "Min-Heap priority queue: " << diff_heap.count() << " seconds\n";
                std::cout << "----------------------------------------\n\n";

                this->temp_meshes.push_back(base_mesh); }, [this, sv, ev]()
					  {
                if (!this->temp_meshes.empty()) {
                    this->clear_scene();
                    geodesic_drawable = nullptr;

                    this->current_meshes = this->temp_meshes;
                    this->temp_meshes.clear();

                    for(auto m : this->current_meshes) this->add_model(m);

                    this->visualize_path(sv, ev);
                } });
		}
	}

	if (ImGui::Button("Dump Matrix") && !current_meshes.empty())
	{
		run_async("Dumping Geodesic Matrix (Multithreaded)", [fp = std::string(filepath)]()
				  {
            SurfaceMesh* base_mesh = SurfaceMeshIO::load(fp);
            if (!base_mesh) return;
            auto edge_lengths = base_mesh->add_edge_property<float>("e:length", 0.0f);
            for (auto e : base_mesh->edges()) {
                SurfaceMesh::Halfedge he = base_mesh->halfedge(e, 0);
                SurfaceMesh::Vertex v0 = base_mesh->source(he);
                SurfaceMesh::Vertex v1 = base_mesh->target(he);
                edge_lengths[e] = length(base_mesh->position(v0) - base_mesh->position(v1));
            }
            compute_and_dump_all_pairs(base_mesh, edge_lengths);
            delete base_mesh; }, []() {});
	}

	ImGui::Separator();
	ImGui::Text("Subdivision");
	ImGui::RadioButton("Sqrt(3)", &subdiv_type, 0);
	ImGui::SameLine();
	ImGui::RadioButton("4-to-1", &subdiv_type, 1);
	ImGui::SameLine();
	ImGui::RadioButton("Phong", &subdiv_type, 2);
	ImGui::SameLine();
	ImGui::RadioButton("All Three", &subdiv_type, 3);
	ImGui::InputInt("Iterations", &subdiv_iters);

	if (ImGui::Button("Run Subdivision") && !current_meshes.empty())
	{
		run_async("Running Subdivision", [this, fp = std::string(filepath), fp_orig = std::string(filepath_orig), type = subdiv_type, iters = subdiv_iters]()
				  {
            SurfaceMesh *base_mesh = SurfaceMeshIO::load(fp);
            if (!base_mesh)
                return;

            SurfaceMesh *orig_mesh = nullptr;
            if (!fp_orig.empty()) {
                orig_mesh = SurfaceMeshIO::load(fp_orig);
            }

            if (!orig_mesh) {
                orig_mesh = SurfaceMeshIO::load(fp);
            }

            std::cout << "\n----------------------------------------\n";
            std::cout << "Subdivision Area & Count Statistics:\n";
            print_mesh_stats(orig_mesh, "Original High-Res Ref.");
            print_mesh_stats(base_mesh, "Decimated Input");

            std::vector<vec3> reference_tris;
            for (auto f : orig_mesh->faces()) {
                auto h = orig_mesh->halfedge(f);
                vec3 p0 = orig_mesh->position(orig_mesh->source(h));
                auto h_curr = orig_mesh->next(h);
                while (orig_mesh->target(h_curr) != orig_mesh->source(h)) {
                    vec3 p1 = orig_mesh->position(orig_mesh->source(h_curr));
                    vec3 p2 = orig_mesh->position(orig_mesh->target(h_curr));
                    reference_tris.push_back(p0);
                    reference_tris.push_back(p1);
                    reference_tris.push_back(p2);
                    h_curr = orig_mesh->next(h_curr);
                }
            }

            if (iters == 0) {
                this->temp_meshes.push_back(base_mesh);
                if (orig_mesh && orig_mesh != base_mesh)
                    delete orig_mesh;
                std::cout << "----------------------------------------\n\n";
                return;
            }

            auto process_subdiv = [&](int t, const std::string &name) -> SurfaceMesh * {
                SurfaceMesh *m = SurfaceMeshIO::load(fp);
                SurfaceMesh *res = run_subdiv(m, t, iters);
                print_mesh_stats(res, name + " Result");

                int n_sub = res->n_vertices();
                std::vector<float> dists(n_sub, 0.0f);
                int n_tris = reference_tris.size() / 3;

                if (n_tris == 0) {
                    std::fill(dists.begin(), dists.end(), 0.0f);
                } else {
                    int num_threads = std::max(1u, std::thread::hardware_concurrency());
                    std::vector<std::thread> threads;
                    std::atomic<int> current_idx(0);

                    for (int th = 0; th < num_threads; ++th) {
                        threads.emplace_back([&]() {
                        while (true) {
                            int i = current_idx.fetch_add(1);
                            if (i >= n_sub) break;
                            
                            vec3 p = res->position(SurfaceMesh::Vertex(i));
                            float min_d = std::numeric_limits<float>::infinity();
                            
                            for (int j = 0; j < n_tris; ++j) {
                                float d = point_triangle_distance_hw1(p, reference_tris[j*3], reference_tris[j*3+1], reference_tris[j*3+2]);
                                if (d < min_d) min_d = d;
                            }
                            dists[i] = min_d;
                        } 
                        });
                    }

                    for (auto &th : threads) {
                        if (th.joinable())
                            th.join();
                    }
                }

                float max_d = 1e-8f;
                for (float d : dists) {
                    if (!std::isinf(d) && !std::isnan(d) && d > max_d)
                        max_d = d;
                }

                std::cout << "   -> Max distance to original ref: " << max_d << "\n";

                auto v_color = res->add_vertex_property<vec3>("v:color", vec3(1, 1, 1));
                for (int i = 0; i < n_sub; ++i) {
                    float d = dists[i];
                    if (std::isinf(d) || std::isnan(d))
                        d = max_d;
                    float normalized_d = std::max(0.0f, std::min(1.0f, d / max_d));
                    v_color[SurfaceMesh::Vertex(i)] = vec3(normalized_d, 0.0f, 1.0f - normalized_d);
                }
                return res;
            };

            if (type == 3) {                                 
                SurfaceMesh *m0 = base_mesh; 
                SurfaceMesh *m1 = process_subdiv(0, "Sqrt(3)");
                SurfaceMesh *m2 = process_subdiv(1, "4-to-1");
                SurfaceMesh *m3 = process_subdiv(2, "Phong");

                float min_x = std::numeric_limits<float>::infinity();
                float max_x = -std::numeric_limits<float>::infinity();
                for (auto v : m0->vertices()) {
                    min_x = std::min(min_x, m0->position(v).x);
                    max_x = std::max(max_x, m0->position(v).x);
                }
                float offset = (max_x - min_x) * 1.25f;

                translate_mesh(m1, vec3(offset, 0, 0));
                translate_mesh(m2, vec3(offset * 2, 0, 0));
                translate_mesh(m3, vec3(offset * 3, 0, 0));

                this->temp_meshes.push_back(m0);
                this->temp_meshes.push_back(m1);
                this->temp_meshes.push_back(m2);
                this->temp_meshes.push_back(m3);
            } else {
                SurfaceMesh *res = process_subdiv(type, "Subdivided");
                this->temp_meshes.push_back(res);
                delete base_mesh;
            }

            std::cout << "----------------------------------------\n\n";
            if (orig_mesh && orig_mesh != base_mesh)
                delete orig_mesh; }, [this]()
				  {
            if (!this->temp_meshes.empty()) {
                this->clear_scene();
                geodesic_drawable = nullptr; 

                this->current_meshes = this->temp_meshes;
                this->temp_meshes.clear();

                for(auto m : this->current_meshes) this->add_model(m);
                
                this->fit_screen();
                this->update();
            } });
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