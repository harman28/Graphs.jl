# A-star algorithm for single-source shortest path
# A-star is an improvement over Dijkstra's Algorithm, that uses heuristics

###################################################################
#
#   Result
#
###################################################################

type AstarResult{V,D}
    seq::Array{V}      #result returns a sequence of nodes that make up the shortest path
    dist::D                 #final distance
end

immutable AstarHEntry{V,D}
    vertex::V
    g_score::D
    h_score::D
end

< (e1::AstarHEntry, e2::AstarHEntry) = e1.g_score+e1.h_score < e2.g_score+e2.h_score


###################################################################
#
#   Algorithm
#
###################################################################



function astar_shortest_paths!{V, D, Heap, H}(
    graph::AbstractGraph{V},                # the graph
    edge_dists::Vector{D},                  # distances associated with edges    
    source::V,             # the sources 
    target::V,              #goal
    heuristic::Function)    #estimator
    @graph_requires graph incidence_list                

    
    # Initialise
    # closedset - stores all nodes that have completed evaluation, initially empty
    # heap  -     Equivalent of priority queue of Dijkstra's algorithm
    # hmap -      index of elements in heap. Useful for updating
    # camefrom -  Stores parent nodes
    # openset -   Nodes ready to be evaulated, initially only source
    n = num_vertices(g)
    closedset = zeros(Int, n)
    hmap=zeros(Int,n)
    camefrom = Array(V, n)
    openset = zeros(Int, n)

    sindex::Int = vertex_index(source, graph)
    openset[sindex]=1
    camefrom[sindex]=typemin(D)

    heap = mutable_binary_minheap(AstarHEntry{V,D})
    hmap[sindex]=push!(heap,AstarHEntrysource(source,0,heuristic(source,target)))

    # main loop 
    
    while !isempty(heap)
        
        # pick next vertex to include
        entry = pop!(heap)      #lowest f_score = g_score+h_score, i.e. real dist plus est dist
        u::V = entry.vertex
        du::D = entry.g_score
        
        if(u==target)
            return AstarResult(reconstruct(camefrom,u),du) #target achieved

        u_index::Int=vertex_index(source,graph)
        openset[u_index]=0          #remove from openset
        closedset[u_index]=1        #put in closedset

        for e in out_edges(u, graph)        #eval neighours
            v::V = target(e, graph)
            v_index::Int = vertex_index(v, graph)

            if closedset[v_index]==1        #already evaluated neighbour. Ignore.
                continue
            end

            tentative_g_score::Int = du + edge_dists[edge_index(e, graph)]

            #check if neighbour is new node, or is previously discovered with higher g_score
            if openset[v_index]==0 || tentative_g_score<heap[hmap[v_index]].g_score
                camefrom[v_index]=u
                if openset[v_index]==0
                    #place in heap
                    hmap[v_index]=push!(heap,AstarHEntrysource(v,tentative_g_score,heuristic(v,target)))
                    openset[v_index]=1
                else
                    #update in heap
                    update!(heap, hmap[v_index], AstarHEntry(v,tentative_g_score,heuristic(v,target)))
                end
            end
        end

    return AstarResult(None,typemax(D))   #failure
end

function reconstruct{V}(g::AbstractGraph{V},camefrom::Array{V},current::V)
    
    result::Array{V}=[current]
    index::Int=vertex_index(current,g)

    #reverse path from target till source is reached
    while camefrom[index]!=typemin(D)
        result=cat(result,[camefrom[index]])
        index=vertex_index(camefrom[index],g)
    end

    return result
end
