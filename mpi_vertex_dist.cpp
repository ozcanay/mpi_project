// int mpi_vertex_dist(graph_t *graph, int start_vertex, int *result)
// {
//     int num_vertices = graph->num_vertices; 
//     fill_n(result, num_vertices, MAX_DIST);
    
//     auto start_time = Time::now();
    
//     int depth = 0;
//     result[start_vertex] = depth;

//     int keep_going = true;

//     while (keep_going)
//     {
//         keep_going = false;

//         for (int vertex = 0; vertex < num_vertices; vertex++)
//         {
//             if (result[vertex] == depth) {
//                 for (int n = graph->v_adj_begin[vertex]; 
//                     n < graph->v_adj_begin[vertex] + graph->v_adj_length[vertex]; 
//                     n++)
//                 {
//                     int neighbor = graph->v_adj_list[n];

//                     if (result[neighbor] > depth+1)
//                     {
//                         result[neighbor] = depth+1;
//                         keep_going = true;
//                     }
//                 }
//             }
//         }

//         depth++;
//     }

//     //print_result(graph, result, depth);
//     return std::chrono::duration_cast<us>(Time::now()-start_time).count();
// }

int mpi_vertex_dist(graph_t *graph, int start_vertex, int *result)
{
    int P, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // 2 * (P - 1) --> *2 because for recv and send. -1 because every process except me (current process).
    MPI_Request reqs[2*P -2]; 
    MPI_Status stats[2*P -2];
    //initBuffers(sendbuf);

    int num_vertices = graph->num_vertices; //total number of vertices
    int num_vertices_per_process = num_vertices%P == 0? num_vertices/P : num_vertices/P + 1;
    // printf("num vertices per process: %d\n",num_vertices_per_process);
    int dummy = -1;
    int sendbuf[P-1][num_vertices_per_process];
    int recvbuf[P-1][num_vertices_per_process];
    
    int index[P-1] = { 0 };
    //my vertices and results are length of N/P and private to each processes
    int my_result[num_vertices_per_process];
    int result_temp[P*num_vertices_per_process];
    fill_n(result_temp, P*num_vertices_per_process, MAX_DIST);
    fill_n(my_result, num_vertices_per_process, MAX_DIST);
    fill( (int *)sendbuf, (int *)sendbuf + sizeof(sendbuf)/sizeof(int), -1 );
    fill( (int *)recvbuf, (int *)recvbuf + sizeof(recvbuf)/sizeof(int), -1 );
    auto start_time = Time::now();

    int depth = 0;
    result_temp[start_vertex] = depth;
    if(start_vertex/num_vertices_per_process == my_rank) {
        my_result[start_vertex%num_vertices_per_process] = depth; //// ?????
    }
    
    int keep_going = true;

    while (keep_going) {
        keep_going = false;

        for (int vertex = 0; vertex < num_vertices_per_process; vertex++) {
            int abs_vertex = my_rank*num_vertices_per_process + vertex; /// offsetting.

            
            int temp = 0;
            for(int i = 0; i<P; i++) {
                if(i != my_rank) {
                    if(i>my_rank) {
                        temp = i-1;
                    } else {
                        temp = i;
                    }

                    MPI_Irecv(&recvbuf[temp], num_vertices_per_process, MPI_INT, i, depth, MPI_COMM_WORLD, &reqs[temp]);
                }
            }
            /* code */
            if(abs_vertex<num_vertices && my_result[vertex] == depth) {
                for (int n = graph->v_adj_begin[abs_vertex]; n < graph->v_adj_begin[abs_vertex] + graph->v_adj_length[abs_vertex]; n++) {
                    int neighbor = graph->v_adj_list[n];
                    int destination = neighbor/num_vertices_per_process; /// this is not modulo operator!  19/8=2.

                    if(destination == my_rank) {
                        if (my_result[neighbor%num_vertices_per_process] > depth + 1) {
                            my_result[neighbor%num_vertices_per_process] = depth+1;
                            keep_going = true;
                        }  
                    } else { 
                        int temp = 0;
                        if(destination>my_rank) {
                            temp = destination-1;
                        } else {
                            temp = destination;
                        }

                        int temp_index = index[temp];
                        sendbuf[temp][temp_index] = neighbor%num_vertices_per_process;
                        index[temp]++;      
                    }
                }   
            }

            for(int i = 0; i<P; i++) {
                if(i != my_rank) {
                    if(i>my_rank) {
                        temp = i-1;
                    } else {
                        temp = i;
                    }

                    MPI_Isend(&sendbuf[temp], num_vertices_per_process, MPI_INT, i, depth, MPI_COMM_WORLD, &reqs[P-1+temp]);
                    // printf("aydin depth: %d, from: %d, to: %d, index: %d ,items: %d, %d, %d \n", depth, my_rank, i, index[temp], sendbuf[temp][0], sendbuf[temp][1], sendbuf[temp][2]);
                }
            }
            
            fill( (int *)sendbuf, (int *)sendbuf + sizeof(sendbuf)/sizeof(int), -1 );
            fill( index, index + sizeof(index)/sizeof(int), 0);
            // printf("aydin depth: %d, from: %d, index: %d  \n", depth, my_rank, index[my_rank]);
      
            MPI_Waitall(2*P-2, reqs, stats);
            
            for(int i = 0; i<P; i++) {
                if(i != my_rank) {
                    if(i>my_rank) {
                        temp = i-1;
                    } else {
                        temp = i;
                    }
                    // printf("aydin depth: %d, from: %d, to: %d, index: %d ,items: %d, %d, %d \n", depth, i, my_rank, index[temp], recvbuf[temp][0], recvbuf[temp][1], recvbuf[temp][2]);
                    int elt = 0;
                    while(recvbuf[temp][elt]!=-1) {
                        int shared_vertice = recvbuf[temp][elt];
                        if (my_result[shared_vertice] > depth + 1) {
                            my_result[shared_vertice] = depth+1;
                            keep_going = true;
                        }  
                        elt++;
                    }
                }
            }
        }
        // printf("aydin depth = %d, rank = %d, depths = %d, %d, %d\n", depth, my_rank,  my_result[0],  my_result[1],  my_result[2] );
 
        bool global_keep_going;
        MPI_Allreduce(&keep_going, &global_keep_going, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        keep_going = global_keep_going;   
        
        MPI_Allgather(my_result, num_vertices_per_process, MPI_INT,
                    result_temp, num_vertices_per_process, MPI_INT,
                     MPI_COMM_WORLD);
                     
        MPI_Barrier(MPI_COMM_WORLD);
        for (size_t i = 0; i < num_vertices_per_process; i++) {
            /* code */
            my_result[i] = result_temp[num_vertices_per_process*my_rank+i];
        }
        
        depth++;
    }
    
    // print_result(graph, result, depth);
    memcpy(result, result_temp, sizeof(int)*num_vertices);
    return std::chrono::duration_cast<us>(Time::now()-start_time).count();
    
}