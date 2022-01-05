// int printArray(std::ofstream& myfile, int *arr, int m, int n) {
//     int i, j;
//     for (i = 0; i < m; i++) {
//       for (j = 0; j < n; j++) {
//             myfile <<  *((arr+i*n) + j) << ' ';
//       }
//       myfile << std::endl;
//     }
// }

// /// TODO: ignore requests.

// int mpi_frontier(graph_t* graph, int start_vertex, int* result) {
//     int my_rank = -1;
//     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

//     std::ofstream myfile;
//     myfile.open (std::to_string(my_rank) + ".txt");

//     myfile <<  "rank: " << my_rank << std::endl;

//     int world_size = -1;
//     MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//     myfile <<  "world size: " << world_size << std::endl;

//     const int num_vertices = graph->num_vertices; 
//     myfile <<  "total vertex count " << num_vertices << std::endl;
//     fill_n(result, num_vertices, MAX_DIST);

//     auto start_time = Time::now();
    
//     int depth = 0;
//     result[start_vertex] = depth;

//     // int vertex_index_per_process[world_size] = {0};

//     int *frontier_in = new int[num_vertices];
//     int *frontier_out = new int[num_vertices];
//     frontier_in[0] = start_vertex;
//     int front_in_size = 1;
//     int front_out_size = 0;

//     int index[world_size - 1] = {0};

//     MPI_Request reqs[2*world_size - 2]; 
//     MPI_Status stats[2*world_size - 2];

//     while (front_in_size != 0) {
//         myfile << "-----------------------------" << std::endl;
//         myfile << "depth: " << depth << std::endl;
//         myfile << "front_in_size: " << front_in_size << std::endl;
//         myfile << "frontier_in contents: ";
//         for(int i = 0; i < front_in_size; ++i) {
//             myfile << frontier_in[i] << ' ';
//         }
//         myfile << std::endl;

//         front_out_size = 0;

//         const int vertex_count_per_process = (front_in_size % world_size == 0) ? (front_in_size / world_size) : (front_in_size / world_size + 1);
//         myfile <<  "vertex count per process: " << vertex_count_per_process << std::endl;

//         /// send and receive buffers
//         int send_buf[world_size - 1][vertex_count_per_process];
//         int recv_buf[world_size - 1][vertex_count_per_process];
//         std::fill((int*)send_buf, (int*)send_buf + sizeof(send_buf)/sizeof(int), -1);
//         std::fill((int*)recv_buf, (int*)recv_buf + sizeof(recv_buf)/sizeof(int), -1);

//         /// result arrays. what is the difference between the two?
//         int my_result[vertex_count_per_process];
//         std::fill_n(my_result, vertex_count_per_process, MAX_DIST);

//         int result_temp[world_size * vertex_count_per_process];
//         std::fill_n(result_temp, world_size * vertex_count_per_process, MAX_DIST);
        
//         for (int v = 0; v < vertex_count_per_process; v++) {
//             /// recv calls. 
//             for(int i = 0; i < world_size; ++i) {
//                 if(i != my_rank) {
//                     int temp = 0;
//                     if(i > my_rank) {
//                         temp = i-1;
//                     } else {
//                         temp = i;
//                     }

//                     myfile << "MPI_Irecv called -> requesting from rank: " << i << " with tag (depth): " << depth << std::endl; 
//                     MPI_Irecv(&recv_buf[temp], vertex_count_per_process, MPI_INT, i, depth, MPI_COMM_WORLD, &reqs[temp]);
//                 } else {
//                     // myfile << "rank " << i << " is not calling MPI_Irecv on itself!" << std::endl;
//                 }
//             }

//             /// core logic.
//             const int frontier_in_index = my_rank * vertex_count_per_process + v;
//             myfile << "frontier_in_index: " << frontier_in_index << std::endl;
//             if(frontier_in_index < front_in_size) {
//                 int vertex = frontier_in[frontier_in_index];
//                 myfile <<  "rank: " << my_rank << ", frontier_in_index: " << frontier_in_index << ", vertex: " << vertex << std::endl;

//                 for (int n = graph->v_adj_begin[vertex]; n < graph->v_adj_begin[vertex] + graph->v_adj_length[vertex]; n++) {
//                     int neighbor = graph->v_adj_list[n];
//                     int destination = neighbor / vertex_count_per_process; /// this is not modulo operator!  19/8=2.

//                     if(destination == my_rank) {
//                         if (my_result[neighbor % vertex_count_per_process] > depth+1) {
//                             my_result[neighbor % vertex_count_per_process] = depth+1;
//                             frontier_out[front_out_size] = neighbor;
//                             front_out_size++;
//                         }
//                     } else {
//                         int dest_rank = 0;
//                         if(destination > my_rank) {
//                             dest_rank = destination - 1;
//                         } else {
//                             dest_rank = destination;
//                         }

//                         int temp_index = index[dest_rank];
//                         send_buf[dest_rank][temp_index] = neighbor % vertex_count_per_process;
//                         index[dest_rank]++;
//                     }
//                 }
//             }

//             myfile << "process " << my_rank << " printing send_buf:" << std::endl;
//             printArray(myfile, (int*)send_buf, world_size - 1, vertex_count_per_process);

//             /// send calls.
//             for(int i = 0; i < world_size; ++i) {
//                 if(i != my_rank) {
//                     int temp = 0;
//                     if(i > my_rank) {
//                         temp = i - 1;
//                     } else {
//                         temp = i;
//                     }


//                     myfile << "send called -> sending to rank: " << i << " with tag (depth): " << depth << std::endl; 
//                     MPI_Isend(&send_buf[temp], vertex_count_per_process, MPI_INT, i, depth, MPI_COMM_WORLD, &reqs[world_size - 1 + temp]);
//                 } else {
//                     // myfile << "rank " << i << " is not calling MPI_Irecv on itself!" << std::endl;
//                 }
//             }

//             /// cleanup
//             std::fill( (int *)send_buf, (int *)send_buf + sizeof(send_buf)/sizeof(int), -1 );
//             std::fill( index, index + sizeof(index) / sizeof(int), 0);
      
//             /// sync point between recv and send calls.
//             myfile << "rank " << my_rank << " waits just before MPI_Waitall" << std::endl;
//             MPI_Waitall(2*world_size-2, reqs, stats);
//             myfile << "rank " << my_rank << " MPI_Waitall wait ended." << std::endl;

//             myfile << "process " << my_rank << " printing recv_buf:" << std::endl;
//             printArray(myfile, (int*)recv_buf, world_size - 1, vertex_count_per_process);

//             /// processing received values into my private array.
//             for(int i = 0; i < world_size; ++i) {
//                 if(i != my_rank) {
//                     int temp = 0;
//                     if(i > my_rank) {
//                         temp = i - 1;
//                     } else {
//                         temp = i;
//                     }

//                     int elt = 0;
//                     while(recv_buf[temp][elt] != -1) {
//                         int received_vertex = recv_buf[temp][elt];
//                         myfile << "received_vertex: " << received_vertex << std::endl;
                        
//                         if (my_result[received_vertex] > depth+1) {
//                             my_result[received_vertex] = depth+1;
//                             frontier_out[front_out_size] = received_vertex;
//                             front_out_size++;
//                         }

//                         ++elt;
//                     }
//                 }
//             }

//             myfile << "front_out_size: " << front_out_size << std::endl;
//         } /// end of for loop (ended iterating through frontier)

        

//         /// now, process handled all of its vertices, now we will look to see if we should go on or not. 
//         /// First, get the new frontier array.
//         front_in_size = front_out_size;
//         int* temp = frontier_in;
//         frontier_in = frontier_out;
//         frontier_out = temp;

// //         int global_front_in_size;
// //         MPI_Allreduce(&front_in_size, &global_front_in_size, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD); /// try MPI_BOR if this does not work.
// //         front_in_size = global_front_in_size;   
        
//         MPI_Allgather(my_result, vertex_count_per_process, MPI_INT,
//                     result_temp, vertex_count_per_process, MPI_INT,
//                     MPI_COMM_WORLD);
                     
//         MPI_Barrier(MPI_COMM_WORLD);
//         for (size_t i = 0; i < vertex_count_per_process; i++) {
//             /* code */
//             my_result[i] = result_temp[vertex_count_per_process*my_rank+i];
//         }

//         depth++;
//     }

//     return std::chrono::duration_cast<us>(Time::now()-start_time).count();
// }

std::string printVector(const std::vector<int>& vertices)
{
    std::string res;
    for(const auto vertex : vertices) {
        res += std::to_string(vertex);
        res += ' ';
    }
    return res;
}

void printMap(std::ofstream& myfile, const std::map<int, std::pair< std::vector<int>, std::vector<int> > >& rank_vertices)
{
    myfile << "rank -->   <received vertices, sent vertices>" << std::endl;
    for(const auto& rank_vertice : rank_vertices) {
        myfile << rank_vertice.first << " -->   <" << printVector(rank_vertice.second.first) << ">, <" <<  printVector(rank_vertice.second.second) << "> " << std::endl;
    }
}

int mpi_frontier(graph_t* graph, int start_vertex, int* result) {
    int world_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::ofstream myfile;
    myfile.open(std::to_string(world_rank) + ".txt");

    myfile <<  "world rank: " << world_rank << std::endl;

    int world_size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    myfile <<  "world size: " << world_size << std::endl;

    std::vector<int> ranks;
    myfile << "ranks are: ";
    for(auto i = 0; i < world_size; ++i) {
        myfile << i << ' ';
        ranks.push_back(i);
    }
    myfile << std::endl;

    int num_vertices = graph->num_vertices; 
    myfile << "total vertex count in the graph: " << num_vertices << std::endl;
    std::fill_n(result, num_vertices, MAX_DIST);

    auto start_time = Time::now();
    
    int depth = 0;
    result[start_vertex] = depth;

    int *frontier_in = new int[num_vertices];
    int *frontier_out = new int[num_vertices];
    frontier_in[0] = start_vertex;
    int front_in_size = 1;
    int front_out_size = 0;

    std::map<int, std::pair< std::vector<int>, std::vector<int> > > rank_vertices; /// key -> rank, value -> <received vertices, sent vertices>

    MPI_Request send_requests[world_size - 1];
    MPI_Request recv_requests[world_size - 1]; 
    int send_req_index = 0;
    int recv_req_index = 0;

    MPI_Status send_stats[world_size - 1];
    MPI_Status recv_stats[world_size - 1];

    while (front_in_size != 0) {
        myfile << "-----------------------------" << std::endl;
        myfile << "depth: " << depth << std::endl;
        myfile << "front_in_size: " << front_in_size << std::endl;
        myfile << "frontier_in contents: ";
        for(int i = 0; i < front_in_size; ++i) {
            myfile << frontier_in[i] << ' ';
        }
        myfile << std::endl;

        front_out_size = 0;

        const int vertex_per_process = (front_in_size % world_size == 0) ? front_in_size / world_size : front_in_size / world_size + 1;
        myfile << "vertex_per_process: " << vertex_per_process << std::endl;
        std::vector<int> my_results(vertex_per_process, MAX_DIST);

        for (int v = 0; v < vertex_per_process; v++) {
            /// nonblocking recv calls.
            
            recv_req_index = 0;
            for(const auto rank : ranks) { /// use lambdas here: https://stackoverflow.com/questions/11502523/how-to-use-lambda-in-for-each
                if(rank != world_rank) {
                    rank_vertices[rank].first.reserve(vertex_per_process);
                    MPI_Irecv(rank_vertices[rank].first.data(), vertex_per_process, MPI_INT, rank, depth, MPI_COMM_WORLD, &recv_requests[recv_req_index++]);
                    myfile << "MPI_Irecv called -> requesting from rank: " << rank << " with tag (depth): " << depth << std::endl; 
                }
            }

            /// core logic.
            const int frontier_in_index = world_rank * vertex_per_process + v;
            myfile << "frontier_in_index: " << frontier_in_index << std::endl;
            if(frontier_in_index < front_in_size) {
                int vertex = frontier_in[frontier_in_index];

                for (int n = graph->v_adj_begin[vertex]; n < graph->v_adj_begin[vertex] + graph->v_adj_length[vertex]; n++) {
                    int neighbor = graph->v_adj_list[n];
                    myfile << "neighbor: " << neighbor << std::endl;

                    int assigned_to_rank = neighbor / vertex_per_process;
                    myfile << "assigned_to_rank: " << assigned_to_rank << std::endl;

                    if(assigned_to_rank == world_rank) {
                        if(my_results[neighbor % vertex_per_process] > depth+1) {
                            my_results[neighbor % vertex_per_process] = depth+1;
                            frontier_out[front_out_size] = neighbor; /// did not use modulo here!
                            front_out_size++;
                        }
                    } else {
                        rank_vertices[assigned_to_rank].second.push_back(neighbor % vertex_per_process); /// prepare send buffer here.
                        myfile << "process " << assigned_to_rank << " will be sent " << neighbor % vertex_per_process << std::endl;
                    }   
                }
            }

            /// nonblocking send calls.
            send_req_index = 0;
            for(const auto rank : ranks) { /// use lambdas here: https://stackoverflow.com/questions/11502523/how-to-use-lambda-in-for-each
                if(rank != world_rank) {
                    rank_vertices[rank].second.reserve(vertex_per_process);
                    MPI_Isend(rank_vertices[rank].second.data(), vertex_per_process, MPI_INT, rank, depth, MPI_COMM_WORLD, &send_requests[send_req_index++]);
                    myfile << "MPI_Isend called -> sending to rank: " << rank << " with tag (depth): " << depth << std::endl; 
                }
            }

            myfile << "Just before send MPI_Waitall" << std::endl;
            MPI_Waitall(world_size - 1, send_requests, send_stats);
            myfile << "MPI_Waitall send  wait ended." << std::endl;

            myfile << "Just before recv MPI_Waitall" << std::endl;
            MPI_Waitall(world_size - 1, recv_requests, recv_stats);
            myfile << "MPI_Waitall recv  wait ended." << std::endl;

            myfile << "printing map..." << std::endl;
            printMap(myfile, rank_vertices);

            for(int i = 0; i < recv_req_index; ++i) {
                int number_amount = -1;
                MPI_Get_count(&recv_stats[i], MPI_INT, &number_amount);

                myfile << "process " << my_rank << " received " << number_amount << " number(s). Message source = " << recv_stats[i].MPI_SOURCE << ", tag = " << recv_stats[i].MPI_TAG << std::endl;
            }

            myfile << "rank_vertices.size(): " << rank_vertices.size() << std::endl;

            for(const auto& rank_vertice : rank_vertices) {
                auto received_vertices = rank_vertice.second.first;
                myfile << "received_vertices.size(): " << received_vertices.size() << std::endl;
                
                for(auto received_vertex : received_vertices) {
                    myfile << "received_vertex: " << received_vertex << std::endl;
                    if(my_results[received_vertex] > depth+1) {
                        my_results[received_vertex] = depth+1;
                    }
                }
            }
        }

        myfile << "Swapping frontier_in with frontier_out." << std::endl;
        
        front_in_size = front_out_size;
        int* temp = frontier_in;
        frontier_in = frontier_out;
        frontier_out = temp;

        myfile << "After swapping frontier_in: " << std::endl;
        myfile << "front_in_size: " << front_in_size << std::endl;
        myfile << "frontier_in contents: ";
        for(int i = 0; i < front_in_size; ++i) {
            myfile << frontier_in[i] << ' ';
        }
        myfile << std::endl;
        
        // MPI_Allgather(my_result, vertex_count_per_process, MPI_INT,
        //             result_temp, vertex_count_per_process, MPI_INT,
        //             MPI_COMM_WORLD);

        myfile << "Before barrier" << std::endl;             
        MPI_Barrier(MPI_COMM_WORLD);
        myfile << "After barrier" << std::endl;
        // for (size_t i = 0; i < vertex_count_per_process; i++) {
        //     /* code */
        //     my_result[i] = result_temp[vertex_count_per_process*my_rank+i];
        // }

        depth++;
    }

    myfile << "end of while frontier_in loop." << std::endl;

    // print_result(graph, result, depth);

    // std::memcpy(result, result_temp, sizeof(int)*num_vertices);
    return std::chrono::duration_cast<us>(Time::now()-start_time).count();
}
