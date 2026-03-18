#include "../third_party/Easy3D/easy3d/viewer/viewer.h"
#include "../third_party/Easy3D/easy3d/util/initializer.h"
#include "../third_party/Easy3D/easy3d/core/surface_mesh.h"
#include "../third_party/Easy3D/easy3d/fileio/surface_mesh_io.h"
#include "../third_party/Easy3D/easy3d/renderer/drawable_lines.h"

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <chrono>
#include <thread>
#include <atomic>
#include <cstdio>

using namespace easy3d;

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

void dijkstra_array(SurfaceMesh *mesh, SurfaceMesh::Vertex start, std::vector<float> &dist, std::vector<SurfaceMesh::Vertex> &parent, const SurfaceMesh::EdgeProperty<float> &edge_lengths)
{
	int n = mesh->n_vertices();
	dist.assign(n, std::numeric_limits<float>::infinity());
	parent.assign(n, SurfaceMesh::Vertex());
	std::vector<bool> visited(n, false);

	dist[start.idx()] = 0.0f;

	for (int i = 0; i < n; i++)
	{
		float min_dist = std::numeric_limits<float>::infinity();
		int u_idx = -1;

		for (int j = 0; j < n; j++)
		{
			if (!visited[j] && dist[j] < min_dist)
			{
				min_dist = dist[j];
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

void compute_and_dump_all_pairs(SurfaceMesh *mesh, bool use_heap, const SurfaceMesh::EdgeProperty<float> &edge_lengths)
{
	int n = mesh->n_vertices();
	std::cout << n << "\n";

	auto start_time = std::chrono::high_resolution_clock::now();

	std::vector<std::vector<float>> dist_matrix(n, std::vector<float>(n, 0.0f));

	int num_threads = std::thread::hardware_concurrency();
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
                if (use_heap) {
                    dijkstra_min_heap(mesh, v, dist_matrix[i], local_parent, edge_lengths);
                } else {
                    dijkstra_array(mesh, v, dist_matrix[i], local_parent, edge_lengths);
                }
            } });
	}

	for (auto &t : threads)
	{
		t.join();
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end_time - start_time;
	std::cout << diff.count() << "\n";

	FILE *outfile = fopen("geodesic_matrix.txt", "w");
	if (outfile)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n - 1; j++)
			{
				fprintf(outfile, "%f ", dist_matrix[i][j]);
			}
			if (n > 0)
			{
				fprintf(outfile, "%f\n", dist_matrix[i][n - 1]);
			}
		}
		fclose(outfile);
	}
}

void visualize_path(Viewer &viewer, SurfaceMesh *mesh, SurfaceMesh::Vertex start, SurfaceMesh::Vertex end, const std::vector<SurfaceMesh::Vertex> &parent)
{
	std::vector<vec3> path_points;
	std::vector<unsigned int> path_indices;

	SurfaceMesh::Vertex curr = end;
	while (curr.is_valid())
	{
		path_points.push_back(mesh->position(curr));

		if (curr == start)
			break;

		SurfaceMesh::Vertex next = parent[curr.idx()];
		if (next.is_valid())
		{
			int curr_idx = path_points.size() - 1;
			path_indices.push_back(curr_idx);
			path_indices.push_back(curr_idx + 1);
		}
		curr = next;
	}

	auto line_drawable = new LinesDrawable("geodesic_path", mesh);

	line_drawable->update_vertex_buffer(path_points);
	line_drawable->update_element_buffer(path_indices);

	line_drawable->set_impostor_type(LinesDrawable::CYLINDER);
	line_drawable->set_uniform_coloring(vec4(1.0f, 0.0f, 0.0f, 1.0f));
	line_drawable->set_line_width(3.0f);

	viewer.add_drawable(std::shared_ptr<LinesDrawable>(line_drawable));
}

int main(int argc, char **argv)
{
	easy3d::initialize();

	easy3d::Viewer viewer("Sukru HW1");

	easy3d::SurfaceMesh *mesh = easy3d::SurfaceMeshIO::load("./meshes1/1) use for geodesic/timing/dragon.obj");

	if (!mesh)
	{
		return -1;
	}

	auto edge_lengths = mesh->add_edge_property<float>("e:length", 0.0f);
	for (auto e : mesh->edges())
	{
		SurfaceMesh::Vertex v0 = mesh->vertex(e, 0);
		SurfaceMesh::Vertex v1 = mesh->vertex(e, 1);
		edge_lengths[e] = length(mesh->position(v0) - mesh->position(v1));
	}

	std::vector<float> dist;
	std::vector<SurfaceMesh::Vertex> parent;

	SurfaceMesh::Vertex start_v(0);
	SurfaceMesh::Vertex end_v(mesh->n_vertices() / 2);

	dijkstra_min_heap(mesh, start_v, dist, parent, edge_lengths);
	visualize_path(viewer, mesh, start_v, end_v, parent);

	viewer.add_model(mesh);

	return viewer.run();
}