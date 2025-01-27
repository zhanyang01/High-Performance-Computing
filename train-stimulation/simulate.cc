#include <mpi.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <string>
#include <algorithm>

#include "platform_load_time_gen.hpp"

using std::vector;
using std::unordered_map;
using std::priority_queue;
using std::string;
using std::sort;
using std::cout;
using std::endl;
using adjacency_matrix = std::vector<std::vector<size_t>>;

// Structs
struct Train {
	int id;
	char line;
	int direction;
	int position;
	int loading_time;
	int travel_time;
	bool on_platform;
	size_t arrival_time;
	int state; // "hold=0", "load=1", "travel=2" "dummy=3"
};

struct TrainComparator {
	bool operator()(const Train* t1, const Train* t2) const {
		if (t1->arrival_time == t2->arrival_time) {
			return t1->id > t2->id;
		}
		return t1->arrival_time > t2->arrival_time;
	}
};

struct Platform {
	Train* train_on_platform = nullptr;
	Train* train_on_link = nullptr;
	priority_queue<Train*, vector<Train*>, TrainComparator> holding_area;
	
	PlatformLoadTimeGen pltg;

	Platform() : pltg(0) {}
	Platform(const size_t popularity) : pltg(popularity) {}
};

struct Station {
	string name;
	size_t popularity;
	unordered_map<char, bool> is_terminal;
	unordered_map<string, Platform> platforms;  // Map each platform by its outgoing destination station name
};

// Helper function implementations
MPI_Datatype MPI_Train;

void create_train_mpi_type() {
	int block_lengths[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype types[9] = {
		MPI_INT, MPI_CHAR, MPI_INT, MPI_INT, MPI_INT,
		MPI_INT, MPI_C_BOOL, MPI_UNSIGNED_LONG, MPI_INT
	};
	MPI_Aint offsets[9];
	offsets[0] = offsetof(Train, id);
	offsets[1] = offsetof(Train, line);
	offsets[2] = offsetof(Train, direction);
	offsets[3] = offsetof(Train, position);
	offsets[4] = offsetof(Train, loading_time);
	offsets[5] = offsetof(Train, travel_time);
	offsets[6] = offsetof(Train, on_platform);
	offsets[7] = offsetof(Train, arrival_time);
	offsets[8] = offsetof(Train, state);

	MPI_Type_create_struct(9, block_lengths, offsets, types, &MPI_Train);
	MPI_Type_commit(&MPI_Train);
}

void send_to_process(size_t process_id, vector<Train*>& trains) {
	size_t num_trains = trains.size();
	size_t train_array_size = num_trains == 0 ? 1 : num_trains;
	Train* train_array = new Train[train_array_size];
	if (num_trains == 0) {
		Train* dummy_train = new Train;
		dummy_train->id = -1;
		dummy_train->line = 'g';
		dummy_train->direction = 1;
		dummy_train->position = 0;
		dummy_train->loading_time = 1;
		dummy_train->travel_time = 0;
		dummy_train->on_platform = false;
		dummy_train->arrival_time = 0;
		dummy_train->state = 3;
		train_array[0] = (*dummy_train);
	} else {
		for (size_t i = 0; i < num_trains; i++) {
			train_array[i] = (*trains[i]);
		}
	}

	MPI_Send(train_array, num_trains, MPI_Train, process_id, 1, MPI_COMM_WORLD);
}

void add_train_to_holding_area(Train* train, Station* station, string station_name, size_t tick) {
	station->platforms.at(station_name).holding_area.push(train);
}

void update_platform(Train* train, Station* station, string station_name) {
	station->platforms.at(station_name).train_on_platform = nullptr;
	station->platforms.at(station_name).train_on_link = train;
}

void move_train_to_link(Train* train, vector<Station*>& local_stations, size_t start_station_idx, const unordered_map<char, vector<string>>& station_lines,
						const adjacency_matrix &mat, unordered_map<string, int>& name_to_idx) {
	string from_station = station_lines.at(train->line)[train->position];
	string to_station = station_lines.at(train->line)[train->position + train->direction];

	train->travel_time = mat[name_to_idx.at(from_station)][name_to_idx.at(to_station)];

	update_platform(train, local_stations[name_to_idx[from_station] - start_station_idx], to_station);

	train->state = 2;
}

void manage_platform_queues(Station* station) {
	for (auto& [destination, platform] : station->platforms) {
		if (platform.train_on_platform == nullptr && !platform.holding_area.empty()) {
			Train* next_train = platform.holding_area.top();
			platform.holding_area.pop();
			platform.train_on_platform = next_train;
			next_train->loading_time = platform.pltg.next(next_train->id);
			next_train->on_platform = true;
			next_train->state = 1;
		}
	}
}

void spawn_trains(vector<Station*>& local_stations, vector<Train*>& trains, size_t tick, size_t& train_id, const unordered_map<char, vector<string>>& station_lines, 
				  const unordered_map<char, size_t>& num_trains, const unordered_map<size_t, size_t>& station_to_process_map,
				  unordered_map<string, int>& name_to_idx, unordered_map<char, size_t>& trains_spawned_per_line, size_t total_processes, size_t mpi_rank) {
	const vector<char> line_order = {'g', 'y', 'b'};
	unordered_map<size_t, vector<Train*>> process_send_trains_map;
	
	for (size_t i = 1; i < total_processes; i++) {
		vector<Train*> cur;
		process_send_trains_map[i] = cur;	
	}
	
	for (char line : line_order) {
		if (trains_spawned_per_line[line] >= num_trains.at(line)) {
			continue;
		}

	const auto& line_stations = station_lines.at(line);       
		for (size_t terminal = 0; terminal < 2; ++terminal) {
			if (trains_spawned_per_line[line] < num_trains.at(line)) {
				Train* new_train = new Train;
				new_train->id = train_id++;
				new_train->line = line;
				new_train->position = terminal == 0 ? 0 : line_stations.size() - 1;
				new_train->direction = terminal == 0 ? 1 : -1;
				new_train->loading_time = 1;
				new_train->travel_time = 0;
				new_train->on_platform = false;
				new_train->arrival_time = tick;
				new_train->state = 0;
				
				string terminal_station = line_stations[new_train->position];
				string to_station = line_stations[new_train->position + new_train->direction];
				size_t process_id = station_to_process_map.at(name_to_idx[terminal_station]);
				trains.push_back(new_train);
				
				if (process_id == 0) {
					add_train_to_holding_area(new_train, local_stations[name_to_idx[terminal_station]], to_station, tick);
				} else {
					process_send_trains_map[process_id].push_back(new_train);
				}
				trains_spawned_per_line[line]++;
			}
		}
	}
	
	for (size_t i = 1; i < total_processes; i++) {
		send_to_process(i, process_send_trains_map[i]);
	}
	
};

// Function to initialise the stations in their own processes
void initialise_stations(vector<Station*>& local_stations, const vector<string>& station_names, const vector<size_t>& popularities, 
						 size_t start_station_idx, size_t end_station_idx, unordered_map<string, int>& name_to_idx,
						 const unordered_map<char, vector<string>>& station_lines, const adjacency_matrix& mat) {
	for (size_t i = start_station_idx; i < end_station_idx; ++i) {
		Station* new_station = new Station;
		new_station->name = station_names[i];
		new_station->popularity = popularities[i];
		
		local_stations.push_back(new_station);
	}

	for (size_t i = start_station_idx; i < end_station_idx; ++i) {
		Station* station = local_stations[i - start_station_idx];
		for (size_t j = 0; j < station_names.size(); ++j) {
			if (mat[i][j] > 0) {
				string to_station = station_names[j];
				station->platforms.emplace(to_station, Platform(popularities[i]));
			}
		}
	}

	// Set terminal stations
	for (size_t i = start_station_idx; i < end_station_idx; ++i) {
		for (const auto& [line, line_stations] : station_lines) {
			if (station_names[i] == line_stations[0] || station_names[i] == line_stations[line_stations.size() - 1]) {
				local_stations[i - start_station_idx]->is_terminal[line] = true;
			} else {
				local_stations[i - start_station_idx]->is_terminal[line] = false;
			}
		}
	}
}

void receive_from_process(vector<Station*> local_stations, unordered_map<string, int> name_to_idx, const unordered_map<char, vector<string>>& station_lines, size_t start_station_idx,
							size_t tick) {
	MPI_Status status;
	MPI_Probe(0, 1, MPI_COMM_WORLD, &status);
	
	int count;
	MPI_Get_count(&status, MPI_Train, &count);
	
	Train* trains = new Train[count];
	MPI_Recv(trains, count, MPI_Train, 0, 1, MPI_COMM_WORLD, &status);
	
	if (trains[0].state == 3) return ;
	for (int i = 0; i < count; i++) {
		Train* train = &trains[i];
		Station* station = local_stations[name_to_idx[station_lines.at(train->line)[train->position]] - start_station_idx];
		string station_name = station_lines.at(train->line)[train->position + train->direction];
		add_train_to_holding_area(train, station, station_name, tick);
	}
}


// The main simulation function
void simulate(size_t num_stations, const vector<string> &station_names, const std::vector<size_t> &popularities,
			  const adjacency_matrix &mat, const unordered_map<char, vector<string>> &station_lines, size_t ticks,
			  const unordered_map<char, size_t> num_trains, size_t num_ticks_to_print, size_t mpi_rank,
			  size_t total_processes) {
	create_train_mpi_type();

	unordered_map<char, size_t> trains_spawned_per_line = {{'g', 0}, {'y', 0}, {'b', 0}};
	unordered_map<string, int> name_to_idx;

	for (size_t i = 0; i < num_stations; i++) {
		name_to_idx[station_names[i]] = i;
	}

	size_t train_id = 0;

	vector<Train*> trains; 
	
	vector<Station*> local_stations;

	size_t num_stations_per_process = num_stations / total_processes;
	size_t num_stations_plus_one = num_stations - num_stations_per_process * total_processes;
	
	size_t start_station_idx;
	size_t end_station_idx;

	unordered_map<size_t, size_t> station_to_process_map;

	for (size_t i = 0; i < total_processes; i++) {
		if (i < num_stations_plus_one) {
			start_station_idx = i * (num_stations_per_process + 1);
			end_station_idx = start_station_idx + num_stations_per_process + 1;
			for (size_t j = start_station_idx; j < end_station_idx; j++) {
				station_to_process_map[j] = i;
			}
		} else {
			start_station_idx = num_stations_plus_one * (num_stations_per_process + 1) + (i - num_stations_plus_one) * num_stations_per_process;
			end_station_idx = start_station_idx + num_stations_per_process;
			for (size_t j = start_station_idx; j < end_station_idx; j++) {
				station_to_process_map[j] = i;
			}
		}
	}

	if (mpi_rank < num_stations_plus_one) {
		start_station_idx = mpi_rank * (num_stations_per_process + 1);
		end_station_idx = start_station_idx + num_stations_per_process + 1;
	} else {
		start_station_idx = num_stations_plus_one * (num_stations_per_process + 1) + (mpi_rank - num_stations_plus_one) * num_stations_per_process;
		end_station_idx = start_station_idx + num_stations_per_process;
	}

	initialise_stations(local_stations, station_names, popularities, start_station_idx, end_station_idx, name_to_idx, station_lines, mat);

	MPI_Barrier(MPI_COMM_WORLD);
	
	for (size_t tick = 0; tick < ticks; ++tick) {
		// Phase 1: Spawn trains
		if (mpi_rank == 0) {
			spawn_trains(local_stations, trains, tick, train_id, station_lines, num_trains, station_to_process_map, name_to_idx, 
						trains_spawned_per_line, total_processes, mpi_rank);
		} else {
			// Receive incoming trains
			receive_from_process(local_stations, name_to_idx, station_lines, start_station_idx, tick);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// Phase 2: Update trains on links
		unordered_map<size_t, vector<Train*>> process_send_trains_map;

		for (size_t i = 0; i < total_processes; i++) {
			vector<Train*> cur;
			if (i != mpi_rank) process_send_trains_map[i] = cur;
		}

		for (size_t i = start_station_idx; i < end_station_idx; ++i) {
			for (auto& [destination, platform] : local_stations[i - start_station_idx]->platforms) {
				Train* train = platform.train_on_link;
				if (train == nullptr) continue;
				train->travel_time--;
				if (train->travel_time == 0) {
					string from_station = station_lines.at(train->line)[train->position];
					string to_station = station_lines.at(train->line)[train->position + train->direction];
					train->state = 0;
					train->position += train->direction;
					train->arrival_time = tick;
					train->on_platform = false;
					platform.train_on_link = nullptr;
					size_t process_id = station_to_process_map[name_to_idx[to_station]];
					if (process_id == mpi_rank) {
						Station* current_station = local_stations[name_to_idx[to_station] - start_station_idx];
						if (current_station->is_terminal.at(train->line)) {
							train->direction = -train->direction;
							add_train_to_holding_area(train, current_station, from_station, tick);
						} else {
							string to_to_station = station_lines.at(train->line)[train->position + train->direction];
							add_train_to_holding_area(train, current_station, to_to_station, tick);
						}
					} else {
						process_send_trains_map[process_id].push_back(train);
					}	
				}
			}
		}
		
		// send trains arriving at platforms handled by other processes
		for (size_t i = 0; i < total_processes; i++) {
			if (i == mpi_rank) continue;
			send_to_process(i, process_send_trains_map[i]);
		}

		// Phase 3: Update trains on platforms
		for (size_t i = start_station_idx; i < end_station_idx; ++i) {
			for (auto& [destination, platform] : local_stations[i - start_station_idx]->platforms) {
				Train* train = platform.train_on_platform;
				if (train == nullptr) continue;
				if (train->loading_time > 0) train->loading_time--;
				if (train->loading_time == 0) {
					if (platform.train_on_link == nullptr) {
						platform.train_on_platform = nullptr;
						move_train_to_link(train, local_stations, start_station_idx, station_lines, mat, name_to_idx);
						train->on_platform = false;
					}
				}
			}
		}

		// sending trains to other processes
		for (size_t i = 0; i < total_processes; i++) {
			if (i == mpi_rank) continue;
			MPI_Status status;
			MPI_Probe(i, 1, MPI_COMM_WORLD, &status);
			
			int count;
			MPI_Get_count(&status, MPI_Train, &count);
			
			Train* trains = new Train[count];
			MPI_Recv(trains, count, MPI_Train, i, 1, MPI_COMM_WORLD, &status);
			
			// dummy train means no trains
			if (trains[0].state == 3) continue;
			
			for (int i = 0; i < count; i++) {
				Train* train = &trains[i];
				Station* station = local_stations[name_to_idx[station_lines.at(train->line)[train->position]] - start_station_idx];
				if (station->is_terminal.at(train->line)) train->direction = -train->direction;
				string station_name = station_lines.at(train->line)[train->position + train->direction];
				add_train_to_holding_area(train, station, station_name, tick);
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);

		// Phase 4: Queue management - Move from holding area to platform if free(done by all processes)
		for (size_t i = start_station_idx; i < end_station_idx; ++i) {
			manage_platform_queues(local_stations[i - start_station_idx]);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// Collect all the trains in the simulation within each process
		if (tick >= ticks - num_ticks_to_print) {
			vector<Train*> process_trains;
			for (size_t i = start_station_idx; i < end_station_idx; i++) {
				Station* station = local_stations[i - start_station_idx];
				for (auto& [destination, platform] : station->platforms) {
					if (platform.train_on_platform != nullptr) {
						process_trains.push_back(platform.train_on_platform);
					}
					if (platform.train_on_link != nullptr) {
						process_trains.push_back(platform.train_on_link);
					}
					vector<Train*> temp_holding_area;
					while (!platform.holding_area.empty()) {
						Train* train = platform.holding_area.top();
						platform.holding_area.pop();
						process_trains.push_back(train);
						temp_holding_area.push_back(train);
					}
					for (Train* train : temp_holding_area) {
						add_train_to_holding_area(train, station, destination, tick);
					}
				}
			}

			/*for (size_t i = 0; i < total_processes; i++) {
				if (i == mpi_rank) {
					printf("The size of process_trains is %ld\n", process_trains.size());
				}
			}*/


			int local_num_trains = process_trains.size() == 0 ? 1 : process_trains.size();
			int local_num_trains_array[total_processes];

			MPI_Gather(&local_num_trains, 1, MPI_INT, local_num_trains_array, 1, MPI_INT, 0, MPI_COMM_WORLD);

			Train* local_trains = new Train[local_num_trains];
			if (process_trains.size() == 0) {
				Train* dummy_train = new Train;
				dummy_train->id = -1;
				dummy_train->line = 'g';
				dummy_train->direction = 1;
				dummy_train->position = 0;
				dummy_train->loading_time = 1;
				dummy_train->travel_time = 0;
				dummy_train->on_platform = false;
				dummy_train->arrival_time = 0;
				dummy_train->state = 3;
				local_trains[0] = (*dummy_train);
			} else {
				for (int i = 0; i < local_num_trains; i++) {
					local_trains[i] = (*process_trains[i]);
				}
			}


			MPI_Barrier(MPI_COMM_WORLD);
			int number_of_trains = 0;
			if (mpi_rank == 0) {
				for (size_t i = 0; i < total_processes; i++) {
					number_of_trains += local_num_trains_array[i];
				}
			}

			Train* all_trains = new Train[number_of_trains];
			int displacement[total_processes];
			int current_displacement = 0;
			if (mpi_rank == 0) {
				for (size_t i = 0; i < total_processes; i++) {
					displacement[i] = current_displacement;
					current_displacement += local_num_trains_array[i];
				}
			}
			MPI_Gatherv(local_trains, local_num_trains, MPI_Train, all_trains, local_num_trains_array, displacement, MPI_Train, 0, MPI_COMM_WORLD);

			if (mpi_rank == 0) {
				printf("%ld: ", tick);
				vector<string> station_names;
				for (int i = 0; i < number_of_trains; i++) {
					string currtrain;
					Train train = all_trains[i];
					string from_station = station_lines.at(train.line)[train.position];
					string to_station = station_lines.at(train.line)[train.position + train.direction];
					if (train.state == 2) {
						currtrain = train.line + std::to_string(train.id) + "-" + from_station + "->" + to_station;
						station_names.push_back(currtrain);
					} else if (train.state == 1) {
						currtrain = train.line + std::to_string(train.id) + "-" + from_station + "%";
						station_names.push_back(currtrain);
					} else if (train.state == 0) {
						currtrain = train.line + std::to_string(train.id) + "-" + from_station + "#";
						station_names.push_back(currtrain);
					} else {
						continue;
					}
				}

				// sort the station names
				std::sort(station_names.begin(), station_names.end());

				// print out the station names
				for (size_t i = 0; i < station_names.size(); i++) {
					printf("%s ", station_names[i].c_str());
				}

				printf("\n");
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);

	}
}
