/**
 * Implementation of a Graph, which is a collection of nodes and edges
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Random;

public class Graph {

	/**
	 * Graph attributes
	 */
	public ArrayList<String> nodes; //all node IDs
	public static Hashtable<String, Node> nodesInfoTable; //all node objects
	public ArrayList<String> activeNodes; //active node IDs
	public ArrayList<String> dynamicNodes;
	public Hashtable<String, Edge> edges; //edge objects
	ArrayList<Edge> edgeList;
	public int numActivities; //total number of activities 
	public double totalWeight; //total weight of edges in graph
	public ArrayList<Activity> activities; //all activities
	public Hashtable<String, Activity> activitiesTable; //activities table for easy lookup
	public Hashtable<String, ArrayList<Edge>> nodeAdjacenciesTable; //the nodes that are adjacent to each node
	public Hashtable<String, ArrayList<Edge>> nodeInEdgesTable; //the nodes that are adjacent to each node
	public double threshold;
	public RandomCollection<Edge> edgeCollection;
	public double minEdgeWeight;
	public double minActivityDistance;
	
	public Hashtable<String, Edge> originalEdges; //original edge objects
	public ArrayList<String> originalNodes; //original node ids
	
	/**
	 * Constructs a Graph.
	 * Initialize graph object based on local file
	 * @param edgesFilename
	 */
	public Graph(String nodesFile, String edgesFile, String activitiesFile, double threshold) {
		// nodes
		nodes = new ArrayList<String>(); // nodes are added into here in addEdge()
		nodesInfoTable = new Hashtable<String, Node>(); // nodes are added into here in constructor
		//String: node id
		activeNodes = new ArrayList<String>();
		dynamicNodes = new ArrayList<String>();
		
		// edges
		edges = new Hashtable<String, Edge>(); // edges are added into here in addEdge()
		// String: "edge.node1,edge.node2"
		edgeList = new ArrayList<Edge>(); // edges are added into here in addEdge()
		
		// activities
		numActivities = 0;
		totalWeight = 0;
		activities = new ArrayList<Activity>();
		activitiesTable = new Hashtable<String, Activity>();
		
		// relationship between nodes, edges, activities
		nodeAdjacenciesTable = new Hashtable<String, ArrayList<Edge>>(); // edges going out of a node
		nodeInEdgesTable = new Hashtable<String, ArrayList<Edge>>(); // edges going in a node
		
		this.threshold = threshold;
		edgeCollection = new RandomCollection<Edge>();
		minEdgeWeight = Double.MAX_VALUE;
		minActivityDistance = Double.MAX_VALUE;

		originalEdges = new Hashtable<String, Edge>();
		originalNodes = new ArrayList<String>();
		
		//read in nodes info
		try {	
			String line = new String();
			BufferedReader r = new BufferedReader(new FileReader(new File(nodesFile)));
			
			while((line = r.readLine()) != null) { // line in fomrat "id; x; y"
				
				String[] lineArray = line.split(";"); // input one string, return multiple strings
				if (lineArray.length < 3) {
					System.err.println("Invalid line in nodes file: " + line);
					System.exit(1);
				}
				
				String nodeID = lineArray[0].trim();
				double x = Double.parseDouble(lineArray[1].trim());
				double y = Double.parseDouble(lineArray[2].trim());
				//System.out.println("nodeX = "+x + ", nodeY = "+y);
				nodesInfoTable.put(nodeID, new Node(nodeID,x,y));
				if (!originalNodes.contains(nodeID))
					originalNodes.add(nodeID);
			}
			r.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
		
		//read in edges info
		try {
			
			String line = new String();
			BufferedReader r = new BufferedReader(new FileReader(new File(edgesFile)));
			
			while((line = r.readLine()) != null) {
				
				//TODO: Figure out why edge 2->3 is not present in the edges table by the time we get to dynamicSegmentation
				/*if (line.equals("2;3;1")) {
					System.out.println("found line 2;3;1");
				}*/
				
				String[] lineArray = line.split(";"); // edge in format "node1; node2; weight"
				if (lineArray.length < 3) {
					System.err.println("Invalid line in edges file: " + line);
					System.exit(1);
				}
				
				String node1 = lineArray[0].trim();
				String node2 = lineArray[1].trim();
				double weight = Double.parseDouble(lineArray[2].trim());

				totalWeight += weight;
				if (weight < minEdgeWeight)
					minEdgeWeight = weight;
				
				Edge edge = new Edge(node1,              //from node
									 node2,              //to node
									 weight);            //weight
						
				addEdge(edge);
				backupOriginalEdge(edge);
			}
			r.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
		
		//System.out.println("check point 1");
		
		
		//read in activities info
		try {	
			String line = new String();
			BufferedReader r = new BufferedReader(new FileReader(new File(activitiesFile)));
			
			while((line = r.readLine()) != null) {
				
				String[] lineArray = line.split(";");
				if (lineArray.length < 5) {
					System.err.println("Invalid line in activities file: " + line);
					System.exit(1);
				}
				
				String activityID  = lineArray[0].trim();
				double x = Double.parseDouble(lineArray[1].trim());
				double y = Double.parseDouble(lineArray[2].trim());
				String edgeNode1ID = lineArray[3].trim();
				String edgeNode2ID = lineArray[4].trim();
				Edge edge = edges.get(edgeNode1ID+","+edgeNode2ID);
				Edge backupEdge = originalEdges.get(edgeNode1ID+","+edgeNode2ID);
				
//				for (String key : edges.keySet()) {
//					System.out.println(key);
//				}
//				
				//System.out.println("edge size: "+edges.keySet().size());
				if (edge!=null) {
					Activity activity = new Activity(activityID,x,y,edge);
					activities.add(activity);
					activitiesTable.put(activity.id, activity);
					numActivities++;
					edge.activities.add(activity);
					edge.numActivities++;
					backupEdge.activities.add(activity);
					backupEdge.numActivities++;
					
					//store active nodes based on activities
					if (!activeNodes.contains(edge.node1))
						activeNodes.add(edge.node1);
					if (!activeNodes.contains(edge.node2))
						activeNodes.add(edge.node2);	
				}
				else{
					System.out.println("No edge of this activity: "+ activityID);
				}
			}
			
			//calculate the minimum activity distance
			minActivityDistance = getMinActivityDistance(activities);
			//minActivityDistance = (minEdgeWeight/4);
			//minActivityDistance = getMedActivityDistance(activities);
			r.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * Adds an edge to the graph
	 * @param edge
	 */
	public void addEdge(Edge edge) {

		edges.put(edge.node1+","+edge.node2, edge);
		edgeList.add(edge);
		//edgeCollection.add(edge.weight, edge);
		
		//add the nodes in a separate data structure
		if(!nodes.contains(edge.node1)) 
			nodes.add(edge.node1);
		if(!nodes.contains(edge.node2))
			nodes.add(edge.node2);		

		//add node information to the node adjacencies table
		ArrayList<Edge> nodeAdjacencies = nodeAdjacenciesTable.get(edge.node1);
		if (nodeAdjacencies == null || nodeAdjacencies.size() == 0)
			nodeAdjacencies = new ArrayList<Edge>();
		nodeAdjacencies.add(edge);		 // put edge into its node1's adjacenies	
		nodeAdjacenciesTable.put(edge.node1, nodeAdjacencies); // 
		
		//add node information to the node in edges table
		ArrayList<Edge> nodeInEdges = nodeInEdgesTable.get(edge.node2);
		if (nodeInEdges == null || nodeInEdges.size() == 0)
			nodeInEdges = new ArrayList<Edge>();
		nodeInEdges.add(edge);			
		nodeInEdgesTable.put(edge.node2, nodeInEdges);			
	}	

	/**
	 * Backup the original edge information in the graph
	 * @param edge
	 * @throws CloneNotSupportedException 
	 */
	public void backupOriginalEdge(Edge edge) throws CloneNotSupportedException {
		Edge backupEdge = edge.clone();		
		originalEdges.put(backupEdge.node1+","+backupEdge.node2, backupEdge);
		edgeCollection.add(backupEdge.weight, backupEdge);	
	}
	
	
	/**
	 * deletes the given edge from the graph
	 * @param edge
	 */
	public void deleteEdge(Edge edge) {
	
		try {
		if (edges.get(edge.node1+","+edge.node2)!=null)
			edges.remove(edge.node1+","+edge.node2);
		}
		catch(NullPointerException e) {
			e.printStackTrace();
		}
		
		ArrayList<Edge> nodeAdjacencies = nodeAdjacenciesTable.get(edge.node1);
		if (nodeAdjacencies != null && nodeAdjacencies.size() > 0) {
			nodeAdjacencies.remove(edge);
		}	
		nodeAdjacenciesTable.put(edge.node1, nodeAdjacencies);	
		
		ArrayList<Edge> nodeInEdges = nodeInEdgesTable.get(edge.node2);
		if (nodeInEdges != null && nodeInEdges.size() > 0) {
			nodeInEdges.remove(edge);
		}	
		nodeInEdgesTable.put(edge.node2, nodeInEdges);		
	}
	
	
	/**
	 * Calculates the active nodes in this graph
	 */
	public void calculateActiveNodes() {
		activeNodes = new ArrayList<String>();
		Enumeration<String> en = this.edges.keys();		
		while (en.hasMoreElements()) {	
			Edge currentEdge = this.edges.get(en.nextElement());
			if (currentEdge.numActivities > 0) {
				if (!activeNodes.contains(currentEdge.node1))
					activeNodes.add(currentEdge.node1);
				if (!activeNodes.contains(currentEdge.node2))
					activeNodes.add(currentEdge.node2);
			}
		}
	}	
	
	/**
	 * Calculate the likelihood ratio of the given path
	 * @param path
	 * @param dynamicSegmentation
	 * @return the likelihood ratio of the given path
	 */
	double getStatistic(Path path, boolean dynamicSegmentation) {
		double a = 0;
		
		if (dynamicSegmentation) {
			a = path.calculateNumActivitiesDynamic();
		}
		else {
			for (int i = 0; i < path.edges.size(); i++) {
				Edge edge = this.edges.get(path.edges.get(i).node1+","+path.edges.get(i).node2);
				a+=edge.numActivities;	
			}
		}
		
		if (a == this.numActivities || path.weight == this.totalWeight)
			return 0;
					
		path.numActivities = a;
		
		double w = path.weight;
		double div1 = a/w;
		double aComp = (this.numActivities/1) - a;
		double wComp = (this.totalWeight/2) - w; //account for undirected edges
		//double wComp = (this.totalWeight/1) - w; //directed edges only
		double div2 = aComp/wComp;		
		if (wComp == 0) return 0;
		if (div2 == 0) return 0;			
		return div1/div2;
	}	
	
	void shuffleGraphActivities () {		
		//initialize number of activities on all edges to zero
		Enumeration<String> en = this.edges.keys();		
		while (en.hasMoreElements()) {	
			String edgeID = en.nextElement();
			Edge currentEdge = this.edges.get(edgeID);
			currentEdge.numActivities = 0;
			edges.put(edgeID, currentEdge);
		}
		
		
		//assign activity to new edge, weighted by edge weight
		for (int i = 0; i < activities.size(); i++) {		
			Edge edge = edgeCollection.next();
			//Collections.shuffle(edgeList);
			//Edge edge = edgeList.get(0);	
			edge.numActivities += 1;
		}
		
		
		/*
		Collections.shuffle(edgeList);
		if (activities.size() < edgeList.size()) {
			for (int i = 0; i < activities.size(); i++) {		
				Edge edge = edgeList.get(i);
				edge.numActivities += 1;
			}	
		}
		else {
			//assign activity to new edge, weighted by edge weight
			for (int i = 0; i < activities.size(); i++) {		
				int low = 0;
				int high = edgeList.size();
				Random rand = new Random();
				int randomIndex = rand.nextInt(high-low) + low;			
				Edge edge = edgeList.get(randomIndex);
				//Edge edge = edgeList.get(0);	
				edge.numActivities += 1;
			}		
		}
		*/
		
		/*
		for (int i = 0; i < 100; i++) {		
			Edge edge = edgeCollection.next();
			edge.numActivities += 1;
		}		
		try {
			BufferedWriter w = new BufferedWriter(new FileWriter(new File("data/graphfile.txt")));
			Enumeration<String> en2 = this.edges.keys();		
			while (en2.hasMoreElements()) {	
				String edgeID = en2.nextElement();
				Edge currentEdge = this.edges.get(edgeID);
				//System.out.println(currentEdge.node1 + ";" + currentEdge.node2 + ";" + currentEdge.weight + ";" + currentEdge.numActivities);
				w.write(currentEdge.node1 + ";" + currentEdge.node2 + ";" + currentEdge.weight + ";" + currentEdge.numActivities + "\n");
			}
			w.close();
			System.out.println("new graph file written");
			System.exit(0);
			
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		*/
	}
	
	/**
	 * Reset the graph to its original state (original means before dynamic segmentation), no activities involved
	 */
	void resetGraph() {
		//reset the graph to its original state	
		nodes = new ArrayList<String>();
		dynamicNodes = new ArrayList<String>();
		edges = new Hashtable<String, Edge>();
		edgeList = new ArrayList<Edge>();
		nodeAdjacenciesTable = new Hashtable<String, ArrayList<Edge>>();	
		nodeInEdgesTable = new Hashtable<String, ArrayList<Edge>>();
		
		Enumeration<String> en = this.originalEdges.keys();		
		while (en.hasMoreElements()) {	
			String edgeID = en.nextElement();
			Edge currentEdge = this.originalEdges.get(edgeID);
			//clear out the activities on each edge
			currentEdge.numActivities = 0;
			currentEdge.activities = new ArrayList<Activity>();
			addEdge(currentEdge);
		}			
	}
	
	/**
	 * Shuffle activities on the graph
	 */
	void shuffleGraphActivitiesDynamic2 () {		
		//reset the graph information
		resetGraph();
		
		ArrayList<Activity> shuffledActivities = new ArrayList<Activity>();
		
		//add each activity to a random edge such that it is at least a minimum distance from every other activity
		for (int i = 0; i < activities.size(); i++) {
			Activity activity = activities.get(i);
			Edge edge = null;
			while ((edge = edgeCollection.next())!=null) {
				
				Node edgeNode1 = Graph.nodesInfoTable.get(edge.node1);
				double edgeNode1X = edgeNode1.x;
				double edgeNode1Y = edgeNode1.y;			
				Node edgeNode2 = Graph.nodesInfoTable.get(edge.node2);
				double edgeNode2X = edgeNode2.x;
				double edgeNode2Y = edgeNode2.y;
				
				Random rand2 = new Random();
				double t = rand2.nextDouble();
				double newActivityX = edgeNode1X + t * (edgeNode2X - edgeNode1X);
				double newActivityY = edgeNode1Y + t * (edgeNode2Y - edgeNode1Y);
				
				boolean distanceGreater = true;
				//ensure that the new location of the activity is at least a certain min distance from every other activity 
				for (int j = 0; j < shuffledActivities.size(); j++) {
					Activity shuffledActivity = shuffledActivities.get(j);
					double distBetweenActivities = Math.sqrt(Math.pow((newActivityX-shuffledActivity.x),2)+Math.pow((newActivityY-shuffledActivity.y),2));	
					if (distBetweenActivities < minActivityDistance) {
						distanceGreater = false;
						break;
					}	
				}

				if (distanceGreater) {
					shuffledActivities.add(activity);
					activity.updateActivityInfo(newActivityX, newActivityY, edge);
					edge.numActivities += 1;
					edge.activities.add(activity);	
					break;
				}
			}
		}
		double minDist = getMinActivityDistance(shuffledActivities);
		//System.out.println("min original dist: " + minActivityDistance + ", min shuffled distance: " + minDist);
		if (minActivityDistance > minDist) 
			System.err.println("error shuffling activities!");
	}
	
	/**
	 * Create a random new activity set for this graph
	 */
	void createNewActivitySetForGraph (int numActivities) {		
		//reset the graph information
		resetGraph();
		
		ArrayList<Activity> shuffledActivities = new ArrayList<Activity>();
		
		//add each activity to a random edge such that it is at least a minimum distance from every other activity
		for (int i = 0; i < numActivities; i++) {
			Activity activity = new Activity(Integer.toString(i+1),0,0,edgeList.get(0));
			Edge edge = null;
			while ((edge = edgeCollection.next())!=null) {
				
				Node edgeNode1 = Graph.nodesInfoTable.get(edge.node1);
				double edgeNode1X = edgeNode1.x;
				double edgeNode1Y = edgeNode1.y;			
				Node edgeNode2 = Graph.nodesInfoTable.get(edge.node2);
				double edgeNode2X = edgeNode2.x;
				double edgeNode2Y = edgeNode2.y;
				
				Random rand2 = new Random();
				double t = rand2.nextDouble();
				double newActivityX = edgeNode1X + t * (edgeNode2X - edgeNode1X);
				double newActivityY = edgeNode1Y + t * (edgeNode2Y - edgeNode1Y);
				
				boolean distanceGreater = true;
				//ensure that the new location of the activity is at least a certain min distance from every other activity 
				for (int j = 0; j < shuffledActivities.size(); j++) {
					Activity shuffledActivity = shuffledActivities.get(j);
					double distBetweenActivities = Math.sqrt(Math.pow((newActivityX-shuffledActivity.x),2)+Math.pow((newActivityY-shuffledActivity.y),2));	
					if (distBetweenActivities < minActivityDistance) {
						distanceGreater = false;
						break;
					}	
				}

				if (distanceGreater) {
					shuffledActivities.add(activity);
					activity.updateActivityInfo(newActivityX, newActivityY, edge);
					edge.numActivities += 1;
					edge.activities.add(activity);	
					break;
				}
			}
		}
		double minDist = getMinActivityDistance(shuffledActivities);
		//System.out.println("min original dist: " + minActivityDistance + ", min shuffled distance: " + minDist);
		if (minActivityDistance > minDist) 
			System.err.println("error shuffling activities!");
		
		System.out.println("id,x,y");
		for (int j = 0; j < shuffledActivities.size(); j++) {
			Activity a = shuffledActivities.get(j);
			System.out.println(a.id + "," + a.x + "," + a.y);
		}
	}
	
	/**
	 * @return the minimum pairwise distance between activities
	 * @param listOfActivities - a list of activities
	 */
	double getMinActivityDistance(ArrayList<Activity> listOfActivities) {
		double minActivityDistace = Double.MAX_VALUE;
		for (int i = 0; i < listOfActivities.size(); i++) {
			Activity a = listOfActivities.get(i);
			for (int j = 0; j < listOfActivities.size(); j++) {
				Activity b = listOfActivities.get(j);
				if (!a.id.equals(b.id)) {
					double dist = Math.sqrt(Math.pow((a.x-b.x),2)+Math.pow((a.y-b.y),2));	
					if (dist < minActivityDistace) minActivityDistace = dist;
				}
			}
		}
		return minActivityDistace;
	}
	
	double getAvgActivityDistance(ArrayList<Activity> listOfActivities) {
		
		double sum = 0;
		ArrayList<Double> elements = new ArrayList<Double>();
		for (int i = 0; i < listOfActivities.size(); i++) {
			Activity a = listOfActivities.get(i);
			for (int j = 0; j < listOfActivities.size(); j++) {
				Activity b = listOfActivities.get(j);
				if (!a.id.equals(b.id)) {
					double dist = Math.sqrt(Math.pow((a.x-b.x),2)+Math.pow((a.y-b.y),2));	
					sum += dist;
					elements.add(dist);
				}
			}
		}	
		
		return (double) sum/elements.size();
	}
	
	double getMedActivityDistance(ArrayList<Activity> listOfActivities) {
		
		ArrayList<Double> elements = new ArrayList<Double>();
		for (int i = 0; i < listOfActivities.size(); i++) {
			Activity a = listOfActivities.get(i);
			for (int j = 0; j < listOfActivities.size(); j++) {
				Activity b = listOfActivities.get(j);
				if (!a.id.equals(b.id)) {
					double dist = Math.sqrt(Math.pow((a.x-b.x),2)+Math.pow((a.y-b.y),2));	
					elements.add(dist);
				}
			}
		}	
		
		Collections.sort(elements);
		return elements.get(elements.size()/2);
	}
	
	/**
	 * Initialize the number of activities on all edges to zero
	 */
	void resetEdgeActivitiesToOriginal () {		
		Enumeration<String> en = this.edges.keys();		
		while (en.hasMoreElements()) {	
			String edgeID = en.nextElement();
			Edge currentEdge = this.edges.get(edgeID);
			currentEdge.numActivities = currentEdge.originalNumActivities;
		}
	}
	
	/**
	 * Returns a psuedo-random number between min and max, inclusive.
	 * The difference between min and max can be at most
	 * <code>Integer.MAX_VALUE - 1</code>.
	 *
	 * @param min Minimim value
	 * @param max Maximim value.  Must be greater than min.
	 * @return Integer between min and max, inclusive.
	 * @see java.util.Random#nextInt(int)
	 */
	public static int randInt(int min, int max) {
	    Random rand = new Random();
	    int randomNum = rand.nextInt((max - min) + 1) + min;
	    return randomNum;
	}
	
	/**
	 * Returns a list of edges going out from given node
	 * @param nodeID - the given node id
	 * @return a list of edges going out from given node
	 */
	public ArrayList<Edge> getOutEdges(String nodeID) {		
		return nodeAdjacenciesTable.get(nodeID);
	}
	
	/**
	 * Returns a list of edges going into given node
	 * @param nodeID - the given node id
	 * @return a list of edges going out from given node
	 */
	public ArrayList<Edge> getInEdges(String nodeID) {		
		return nodeInEdgesTable.get(nodeID);
	}
	
	
	@Override
	public String toString() {
		String s = new String();	
		Enumeration<String> en = edges.keys();
		//s += "Graph Edges: \n";
		while (en.hasMoreElements()) {	
			Edge currentEdge = edges.get(en.nextElement());
			s += currentEdge.node1 + " -> " + currentEdge.node2 + "(" + currentEdge.weight + "," + currentEdge.numActivities + ")";
			s += currentEdge.activities + "\n";
		}	
		/*s+="\n";
		Enumeration<String> en2 = originalEdges.keys();
		while (en2.hasMoreElements()) {	
			Edge currentEdge = originalEdges.get(en2.nextElement());
			s += currentEdge.node1 + " -> " + currentEdge.node2 + "(" + currentEdge.weight + "," + currentEdge.numActivities + ")";
			s += currentEdge.activities + "\n";
		}*/			
		return s;
	}
	
	/**
	 * Dynamically segment the network (add new nodes and edges) based on activity locations
	 */
	void dynamicSegmentation() {
		
		ArrayList<Edge> edgesToAdd = new ArrayList<Edge>();
		ArrayList<Edge> edgesToAddRev = new ArrayList<Edge>(); // reversed edge, because all the edges are bi-directioned
		ArrayList<Edge> edgesToDelete = new ArrayList<Edge>();
		
		//for each edge in the graph
		Enumeration<String> en = this.originalEdges.keys();		
		while (en.hasMoreElements()) {	
			String edgeID = en.nextElement();
			Edge edge = this.originalEdges.get(edgeID); // get one original edge object (original means from the static road network)
			
			//if the edge has 1 or more activities
			if (edge.numActivities > 0) {
				
				//sort activities on an edge based on their distance to node 1 of the edge
				ArrayList<Activity> edgeActivities = edge.activities;
				Collections.sort(edgeActivities, activityDistanceComparator);
				// a. is a prefix of dynamic edge
				//create a dynamic edge between edge node 1 and first activity in list
				Edge dynamicEdgeFirst = new Edge(edge.node1, "a." + edgeActivities.get(0).id, edgeActivities.get(0).distanceToEdgeNode1);			
				//addEdge(dynamicEdgeFirst);
				edgesToAdd.add(dynamicEdgeFirst);
				
				Edge dynamicEdgeFirstRev = new Edge("a."+edgeActivities.get(0).id, edge.node1, edgeActivities.get(0).distanceToEdgeNode1);			
				//addEdge(dynamicEdgeFirstRev);
				edgesToAddRev.add(dynamicEdgeFirstRev);
				
				//create edges between each activity based on its distance to node 1 of the edge
				for (int i = 0; i < edgeActivities.size()-1; i++) {
					
					//note this activity as a special node for later
					dynamicNodes.add("a."+edgeActivities.get(i).id);
					
					//create a dynamic edge between activities
					Edge dynamicEdge = new Edge("a."+edgeActivities.get(i).id,
												"a."+edgeActivities.get(i+1).id,
												edgeActivities.get(i+1).distanceToEdgeNode1 - edgeActivities.get(i).distanceToEdgeNode1);
					//addEdge(dynamicEdge);
					edgesToAdd.add(dynamicEdge);
					
					Edge dynamicEdgeRev = new Edge("a."+edgeActivities.get(i+1).id,
				                        		   "a."+edgeActivities.get(i).id,
							                       edgeActivities.get(i+1).distanceToEdgeNode1 - edgeActivities.get(i).distanceToEdgeNode1);
					//addEdge(dynamicEdgeRev);
					edgesToAddRev.add(dynamicEdgeRev);
				}
				
				//note this activity as a special node for later
				dynamicNodes.add("a."+edgeActivities.get(edgeActivities.size()-1).id); // the last activities on the edge
				
				//create a dynamic edge between last activity in list and edge node 2
				Edge dynamicEdgeLast = new Edge("a."+edgeActivities.get(edgeActivities.size()-1).id,
												edge.node2,
												edgeActivities.get(edgeActivities.size()-1).distanceToEdgeNode2);
				//addEdge(dynamicEdgeLast);
				edgesToAdd.add(dynamicEdgeLast);
				
				Edge dynamicEdgeLastRev = new Edge(edge.node2,
						                           "a."+edgeActivities.get(edgeActivities.size()-1).id,
						                           edgeActivities.get(edgeActivities.size()-1).distanceToEdgeNode2);
				//addEdge(dynamicEdgeLastRev);				
				edgesToAddRev.add(dynamicEdgeLastRev);
				
				//remove the original edge from the graph
				//deleteEdge(edge);
				edgesToDelete.add(edge);
				
				Edge edgeRev = edges.get(edge.node2+","+edge.node1);
				//deleteEdge(edgeRev);
				edgesToDelete.add(edgeRev);
			}		
				
			//add dynamic forward edges
			for (int i = 0; i < edgesToAdd.size(); i++) {	
				Edge edgeToAdd = edgesToAdd.get(i);
				addEdge(edgeToAdd);

				//if we found an activity
				if (edgeToAdd.node2.startsWith("a.")) {
					//get the activity
					Activity activity = null;
					String[] na = edgeToAdd.node2.split("\\.");
					if (na.length > 1) {
						activity = activitiesTable.get(na[1]); //lookup activity in the activities table
					}
					
					if (activity!=null) {						
						ArrayList<Edge> pathEdgesToNode1 = new ArrayList<Edge>(edgesToAddRev.subList(0, i+1));		
						ArrayList<Edge> pathEdgesToNode2 = new ArrayList<Edge>(edgesToAdd.subList(i+1, edgesToAdd.size()));	
						//set paths to node 1 and node 2
						if (pathEdgesToNode1.size()>0 && pathEdgesToNode2.size()>0) {
							activity.pathToNode1 = new Path(pathEdgesToNode1);
							activity.pathToNode2 = new Path(pathEdgesToNode2);
						}
					}
				}
			}
			
			//add dynamic reverse edges
			for (int i = 0; i < edgesToAddRev.size(); i++) {
				addEdge(edgesToAddRev.get(i));
			}
			
			//reset forward and reverse dynamic edge lists
			edgesToAdd = new ArrayList<Edge>();
			edgesToAddRev = new ArrayList<Edge>();
		}
		
		//delete all edges that have been segmented
		for (int i = 0; i < edgesToDelete.size(); i++) {
			deleteEdge(edgesToDelete.get(i));
		}
		
	}
	
    private final Comparator<Activity> activityDistanceComparator = new Comparator<Activity>() {
	    public int compare(Activity a1, Activity a2) {	        
	        if (a1.distanceToEdgeNode1 > a2.distanceToEdgeNode1) return 1;	        
	        else if (a1.distanceToEdgeNode1 < a2.distanceToEdgeNode1) return -1;
	        else return 0;                       
	    }
    };				
	
	public static void main(String[] args) {
		
		if (args.length < 1) {
			System.err.println("Input file location must be specified");
			System.exit(1);
		}
		
		Graph g = new Graph(args[0],args[1],args[2],0);
		//g.dynamicSegmentation();
		//g.createNewActivitySetForGraph(995);
		//g.shuffleGraphActivitiesDynamic2();
		//System.out.println(g);
		
		
		/*System.out.println("id,x,y");
		for (int i = 0; i < g.activities.size(); i++) {
			Activity a = g.activities.get(i);
			System.out.println(a.id + "," + a.x + "," + a.y);
		}*/
		
		System.out.println("all nodes: " + g.nodes.size() + ", edges: " + g.edges.size() + ", activities: " + g.activities.size() + ", weight: " + g.totalWeight);
		//System.out.println(g);
		//g.dynamicSegmentation();
		//g.shuffleGraphActivitiesDynamic();
		//System.out.println("-------------------------------------------------------------------");
		//System.out.println("-------------------------------------------------------------------");
		//System.out.println(g);
		
		/*g.dynamicSegmentation();
		System.out.println("-------------------------------------------------------------------");
		System.out.println("-------------------------------------------------------------------");
		System.out.println(g);*/
		System.out.println("done!");
	}
}
