/******************************************************************************
 * File: Trailblazer.cpp
 *
 * Implementation of the graph algorithms that comprise the Trailblazer
 * assignment.
 */

#include "Trailblazer.h"
#include "TrailblazerGraphics.h"
#include "TrailblazerTypes.h"
#include "TrailblazerPQueue.h"
#include "random.h"

using namespace std;

/* Function: shortestPath
 * 
 * Finds the shortest path between the locations given by start and end in the
 * specified world.	 The cost of moving from one edge to the next is specified
 * by the given cost function.	The resulting path is then returned as a
 * Vector<Loc> containing the locations to visit in the order in which they
 * would be visited.	If no path is found, this function should report an
 * error.
 *
 * In Part Two of this assignment, you will need to add an additional parameter
 * to this function that represents the heuristic to use while performing the
 * search.  Make sure to update both this implementation prototype and the
 * function prototype in Trailblazer.h.
 */
Vector<Loc>
shortestPath(Loc start, Loc end,
			 Grid<double>& world,
			 double costFunction(Loc one, Loc two, Grid<double>& world),
			 double heuristic(Loc start, Loc end, Grid<double>& world)){

	int x = world.numRows();
	int y = world.numCols();

	Vector<Loc> answer;

	Grid<Loc> path(x, y);
	Grid<double> distance(x, y);
	for(int i = 0; i < x; i++){
		for(int j = 0; j < y; j++){
			distance[i][j] = -1;
		}
	}	
	distance[start.row][start.col] = 0;


	TrailblazerPQueue<Loc> yellNodes;
	yellNodes.enqueue(start, 0);
	colorCell(world,start,YELLOW);
		
	while(!yellNodes.isEmpty()){
		Loc curr = yellNodes.dequeueMin();
		colorCell(world, curr, GREEN);

		if(curr == end){
			Stack<Loc> tmp;
			tmp.push(end);
			while(curr != start){
				tmp.push(path[curr.row][curr.col]);
				curr = path[curr.row][curr.col];
			}
			while(!tmp.isEmpty()){
				answer.add(tmp.pop());
			}
			break;
		}

		//Iterate over neighbours (+ itself) 
		for(int i = -1; i <=1; i++){
			for(int j = -1; j <=1; j++){
				int nextX = curr.row + i;
				int nextY = curr.col + j;
				if(nextX >= 0 && nextX < x && nextY >= 0 && nextY < y){
					Loc nextLoc = makeLoc(nextX, nextY);
					double newDistance = distance[curr.row][curr.col] + costFunction(curr, nextLoc, world);

					if(distance[nextX][nextY] == -1) {
						colorCell(world,nextLoc,YELLOW);
						path[nextX][nextY] = curr;
						distance[nextX][nextY] = newDistance;
						yellNodes.enqueue(nextLoc, newDistance + heuristic(nextLoc, end,world));
					}

					else if(distance[nextX][nextY] > newDistance){
						path[nextX][nextY] = curr;
						distance[nextX][nextY] = newDistance;
						yellNodes.decreaseKey(nextLoc, newDistance + heuristic(nextLoc, end,world));
					}
				}
			}
		}

	}
	


    return answer;
}


Set<Edge> createMaze(int numRows, int numCols) {
	
	// Direction arrays.
	int di[2] = {1, 0};
	int dj[2] = {0, 1};
	//

	Map<Loc, HashSet<Loc> *> clusters;
	TrailblazerPQueue<Edge> edges;

	//Create all clusters, all edges and randomise their weights.
	for(int x = 0; x < numRows; x++){
		for(int y = 0; y < numCols; y++){
			// Create new cluster, containing single node.
			Loc newLoc = makeLoc(x, y);
			HashSet<Loc> * newCluster = new HashSet<Loc>;
			newCluster->add(newLoc);
			clusters.put(newLoc, newCluster);
			// Iterate over current node's neighbours (not all of them tho).
			for(int k = 0; k < 2; k++){
				int nextX = x + di[k];
				int nextY = y + dj[k];
				//Create edges.
				if(nextX >= 0 && nextX < numRows && nextY >= 0 && nextY < numCols){
					Loc nextLoc = makeLoc(nextX, nextY);
					Edge newEdge = makeEdge(newLoc, nextLoc);
					//Enques created edge with random weight
					edges.enqueue(newEdge, randomInteger(0, 20));
				}
			}
		}
	}
	

	Set<Edge> answer;

	while(!edges.isEmpty()){
		Edge tempEdge = edges.dequeueMin();
		if(mergeClusters(clusters, tempEdge.start, tempEdge.end)){
			answer.add(tempEdge);
		}
	}
	// As all nodes are united, deleting one node's HashSet will do the job.
	delete clusters[makeLoc(0,0)];
	

    return answer;
}

bool mergeClusters(Map<Loc, HashSet<Loc> *> & clusters, Loc firstLoc, Loc secondLoc){
	if(clusters[firstLoc]->contains(secondLoc)) return false;
	

	HashSet<Loc> * united = new HashSet<Loc>;

	foreach(Loc lc in *clusters[secondLoc]){
		united->add(lc);
	}
	foreach(Loc lc in *clusters[firstLoc]){
		united->add(lc);
	}
	HashSet<Loc> * firstCluster = clusters[firstLoc];
	HashSet<Loc> * secondCluster = clusters[secondLoc];
	foreach(Loc lc in *united){
		clusters[lc] = united;
	}
	
	delete firstCluster;
	delete secondCluster;

	return true;
}
