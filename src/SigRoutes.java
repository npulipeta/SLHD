import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;


public class SigRoutes {

	public double mcIterations;
	public double pvalue;
	
	public Graph graph;
	public ShortestPath sp; 
	public Hashtable<String, Path> spTable;
	public ArrayList<Path> dominantPaths;
	public ArrayList<Double> statistics;
	public double threshold;
	
	static String nodesFile = new String();
	static String edgesFile = new String();
	static String activitiesFile = new String();
	static double lrThreshold = 0;
	
	/**
	 * Construct a SigRoutes object from the given graph
	 * @param graph
	 */
	public SigRoutes(Graph graph, double mcIterations, double pvalue) {		
		this.graph = graph;
		this.threshold = graph.threshold;
		spTable = new Hashtable<String, Path>();
		sp = new ShortestPath(graph);
		dominantPaths = new ArrayList<Path>();
		statistics = new ArrayList<Double>();
		this.mcIterations = mcIterations;
		this.pvalue = pvalue;
	}
	
	public ArrayList<Path> getSignificantRoutes(boolean likelihoodPruning,
												boolean monteCarloSpeedup,
												boolean dynamicSegmentation, //shortest paths calculated between all activities
												boolean dynamicSegmentationHeirarchical, //shortest paths calculated between all nodes - paths between activities derived
												boolean dynamicSegmentationHeirarchicalActive)       //shortest paths calculated between active nodes - paths between activities derived
														throws CloneNotSupportedException { 
		
		
		long ts = System.currentTimeMillis();
		sp.dynamicSegmentation = dynamicSegmentation;
		
		ArrayList<Path> candidatePaths = new ArrayList<Path>();	
		ArrayList<Path> significantPaths = new ArrayList<Path>();	

		Path maxPath = null;
		
		//heirarchical dynamic segmentation
		if (dynamicSegmentationHeirarchical) {	
			//ts = System.currentTimeMillis();
			maxPath = dynamicSegmentationHeirarchical (candidatePaths, likelihoodPruning, false);
			//System.out.println("time for dynamicSegmentationHeirarchical: " + (System.currentTimeMillis() - ts));
		}
		//heirarchical dynamic segmentation with active node pruning
		else if (dynamicSegmentationHeirarchicalActive) { 
			//ts = System.currentTimeMillis();
			maxPath = dynamicSegmentationHeirarchical (candidatePaths, likelihoodPruning, true);
			//System.out.println("time for dynamicSegmentationHeirarchicalActive: " + (System.currentTimeMillis() - ts));
		}
		//basic dynamic segmentation or no dynamic segmentation			
		else {	
			//ts = System.currentTimeMillis();
			ArrayList<String> nodes = new ArrayList<String>();	
			if (dynamicSegmentation) {			
				graph.dynamicSegmentation(); // dynamically segment the graph
				nodes = graph.dynamicNodes;
			}
			else {
				nodes = graph.nodes;
			}
			
			for (int i = 0; i < nodes.size(); i++) {			
				String currentNode = nodes.get(i);
				if (likelihoodPruning) {
					sp.processAllSPNodesLP(currentNode);
				}
				else {
					sp.processAllSPNodes(currentNode); //calculating shortest paths from current node to all other nodes			
				}
				for (int j = 0; j < nodes.size(); j++) {
					if (nodes.get(j) != currentNode) {
						
						if (dynamicSegmentation) {
							if (!graph.dynamicNodes.contains(nodes.get(j)))
								continue;
						}
						
						ArrayList<Edge> a = sp.getShortestPathEdgesNoProcessing(currentNode, nodes.get(j));	
						if (a.size() > 0) {
							Path p = new Path(a);					
							p.statistic = graph.getStatistic(p, dynamicSegmentation);
							spTable.put("node"+currentNode+",node"+nodes.get(j), p); //store shortest path for later					
							//test each shortest path to see if it meets the scan statistic
							if ((p.statistic >= threshold) && (p.weight >= 0.01)) {
								if (p.numActivities > 1) {
									candidatePaths.add(p.clone());	
									if (maxPath == null || maxPath.statistic < p.statistic)
										maxPath = p;
								}
								//System.out.println("path: "+p);
							}
						}
					}
				}		
			}	
			//System.out.println("time for basic or no dynamicSegmentation " + (System.currentTimeMillis() - ts));
		}

		//if there is no max path then terminate
		if (maxPath == null)
			return significantPaths;
		
		maxPath = maxPath.clone();
		System.out.println("original max path: " + maxPath);
		//System.out.println("candidatePaths: " + candidatePaths);
		System.out.println("candidatePaths size: " + candidatePaths.size());
		
		
		int maxmax = 0;
		double maxmaxLR = 0;
		for (int i = 0; i < candidatePaths.size(); i++){
			//System.out.println("No." + i + " path: "+candidatePaths.get(i).toString());
			if (candidatePaths.get(i).edges.size() > maxmax){
				maxmax = candidatePaths.get(i).edges.size();
			}
			if (candidatePaths.get(i).statistic > maxmaxLR) {
				maxmaxLR = candidatePaths.get(i).statistic;
			}
		}
		System.out.println("candidiate path edge size(): " + maxmax);
		System.out.println("max LR in candidate paths: " + maxmaxLR);
		
		
		System.out.println("time for generating candidate paths: " + (System.currentTimeMillis() - ts));
		
		//System.out.println("Calculating significance through Monte Carlo simulations...");
		//test significance of candidate paths
		//significantPaths = getSignificance(candidatePaths, maxPath, likelihoodPruning, monteCarloSpeedup,dynamicSegmentation);	
		long ts2 = System.currentTimeMillis();
		significantPaths = getSignificanceSameRoutes(candidatePaths, maxPath, likelihoodPruning, monteCarloSpeedup,dynamicSegmentation);	
		System.out.println("time for MC simulation: " + (System.currentTimeMillis() - ts2));
		//significantPaths = candidatePaths;
		
		long ts3 = System.currentTimeMillis();
		//get dominant paths
		for (int i = 0; i < significantPaths.size(); i++) {
			Path currentPath = significantPaths.get(i);			
			boolean isSubPath = false;
			for (int j = 0; j < significantPaths.size(); j++) {		
				Path otherPath = significantPaths.get(j);	
				if (i != j) {
					if (currentPath.isSubPathOf(otherPath)) {
						isSubPath = true;
						break;
					}
				}
			}
			if (!isSubPath)
				addDominantPath(currentPath);
		}
		System.out.println("time for dominant paths: " + (System.currentTimeMillis() - ts3));
		
		return dominantPaths;
	}	
	
	/**
	 * Calculate the significance of the given set of paths considering all shortest paths
	 * @param paths
	 * @return
	 */
	ArrayList<Path> getSignificanceSameRoutes (ArrayList<Path> paths, 
									 Path originalMaxPath,
									 boolean likelihoodPruning, 
									 boolean monteCarloSpeedup,
									 boolean dynamicSegmentation) {

		//long ts1 = System.currentTimeMillis();
		
		ArrayList<Path> significantPaths = new ArrayList<Path>();
		
		if (paths.size() == 0 || paths == null)
			return significantPaths;
		
		int numTimesPValueBeaten = 0;
		graph.resetGraph(); //reset the graph to its original state
		sp = new ShortestPath(graph);
		spTable = new Hashtable<String, Path>();
		for (int x = 0; x < mcIterations; x++) {
			
			//shuffle number of activities on each edge			
			ArrayList<String> nodes = new ArrayList<String>();	

			graph.shuffleGraphActivities();
			nodes = graph.nodes;
			/*
			if (dynamicSegmentation){
				graph.dynamicSegmentation();
			}
			*/
			//System.out.println("MC nodes size: " + nodes.size());
			
			Path maxPath = null;	
			boolean speedup = false;
			
			
			for (int i = 0; i < nodes.size(); i++) {			
				String currentNode = nodes.get(i);
			
				//if (x == 0) {
					if (likelihoodPruning)
						sp.processAllSPNodesLP(currentNode);
					else
						sp.processAllSPNodes(currentNode); //calculating shortest paths from current node to all other nodes			
				//}
				
				for (int j = 0; j < nodes.size(); j++) {
					if (nodes.get(j) != currentNode) {
						
						Path p = null;
						//if x = 0, store the path for future simulations
						//if (x==0) {	
							ArrayList<Edge> a = sp.getShortestPathEdgesNoProcessing(currentNode, nodes.get(j));	
							if (a.size() > 0) {
								p = new Path(a);	
								spTable.put("node"+currentNode+",node"+nodes.get(j), p); //store shortest path for later	
							}
						//}
						//else {
						//	p = spTable.get("node"+currentNode+",node"+nodes.get(j));
						//}
						
						if (p != null) {
							p.statistic = graph.getStatistic(p,false);
							
							//if (p.numActivities <= 1) continue;
							
							//test each shortest path to see if it meets the likelihood ratio threshold			
							if (maxPath == null || maxPath.statistic < p.statistic)
								maxPath = p;	
							
							if (monteCarloSpeedup && p.statistic > originalMaxPath.statistic) {
								System.out.println("MC speedup!!");
								addStatistic(p.statistic);
								speedup = true;
								numTimesPValueBeaten++;
								
								//no paths will be significant after consistently finding paths in each simulation
								if ((numTimesPValueBeaten/mcIterations) > pvalue) {
									System.out.println("MC speedup 2!!");
									return significantPaths; 
								}
								break;
							}
						}		
			
						//}
					}
				}
				if (speedup) break;
			}

			if (!speedup) {
				if (maxPath != null)
					addStatistic(maxPath.statistic);
				else
					addStatistic(0.0);
			}
			System.out.println("Simulation " + (x+1) + " complete. Max Path is: " + maxPath); 			
		}

		//System.out.println("time for MC1: " + (System.currentTimeMillis() - ts1));
		
		//long ts2 = System.currentTimeMillis();
		//calculate p-value
		for (int i = 0; i < paths.size(); i++) {
			Path path = paths.get(i);
			double pathStat = path.statistic;		
			
			/*int pNum = 1;
			for (int j = 0; j < statistics.size()-1; j++) {
				if (pathStat < statistics.get(j)) {
					pNum = pNum + 1;
				}
			}*/
			
			
			
			int pNum = Collections.binarySearch(statistics, (pathStat+0.001));
			if (pNum < 0) pNum = -pNum-1;
			pNum = (statistics.size() - pNum) + 1;
			//System.out.println("pNum: " + pNum);
			
			
			double pathPVal = (double) (pNum/mcIterations);
			//System.out.println(path + ", pVal: " + pathPVal);
			if (pathPVal <= pvalue) {			
				path.pVal = pathPVal;
				significantPaths.add(path);		
			}
		}
		
		//System.out.println("time for MC2: " + (System.currentTimeMillis() - ts2));
		return significantPaths;
	}
	
	/**
	 * Add to the collection of statistics for Monte Carlo simulations
	 * Statistics will be sorted smallest to largest
	 * @param statistic
	 */
    public void addStatistic(Double statistic) {
    	int index = Collections.binarySearch(statistics, statistic);
    	
    	if (index < 0)
    		statistics.add(-index-1, statistic);	
    	else 
    		statistics.add(index, statistic);
    }

    /**
     * Add to the collection of dominant paths, sorted by their likelihood
     * @param path
     */
    public void addDominantPath(Path path) {
    	int index = Collections.binarySearch(dominantPaths, path, pathStatisticComparator);
    	if (index < 0)
    		dominantPaths.add(-index-1, path);	
    	else 
    		dominantPaths.add(index, path);
    }
	
    /**
     * Comparator object for paths - used to sort by likelihood
     */
    private final Comparator<Path> pathStatisticComparator = new Comparator<Path>() {
	    public int compare(Path path1, Path path2) {	        	
	        if (path1.statistic > path2.statistic) return 1;	        
	        else if (path1.statistic < path2.statistic) return -1;
	        else return 0;
	                       
	    }
    };	

    /**
     * Dynamically segment the graph using the heirarchical dynamic segmentation approach
     * @param candidatePaths
     * @param maxPath
     * @throws CloneNotSupportedException
     */
    Path dynamicSegmentationHeirarchical(ArrayList<Path> candidatePaths, 
    									 boolean likelihoodPruning, 
    									 boolean active) 
    											 throws CloneNotSupportedException {
    	
    	Path maxPath = null;
    	//if activities exceed nodes
    	if (graph.numActivities > graph.originalNodes.size()) {
    		
    		//algorithmic refinement for only considering active nodes
    		ArrayList<String> nodes = new ArrayList<String>();
    		if (active)
    			nodes = graph.activeNodes;
    		else
    			nodes = graph.originalNodes;
    		
    		
			sp = new ShortestPath(graph);
			
			for (int i = 0; i < nodes.size(); i++) {		
				String currentNode = nodes.get(i);
				sp.processAllSPNodes(currentNode);
				for (int j = 0; j < nodes.size(); j++) {		
					if (nodes.get(j) != currentNode) {
						ArrayList<Edge> a = sp.getShortestPathEdgesNoProcessing(currentNode, nodes.get(j));				
						if (a.size() > 0) {
							Path p = new Path(a);					
							spTable.put("node"+currentNode+",node"+nodes.get(j), p); //store shortest path for later					
						}				
					}	
				}
			}		
			
			//dynamically segment the graph
			graph.dynamicSegmentation();
			
			//for each activity, evaluate the shortest path between each pair of activities
			for (int i = 0; i < graph.activities.size(); i++) {
				Activity sourceActivity = graph.activities.get(i);
				for (int j = 0; j < graph.activities.size(); j++) {
					Activity destActivity = graph.activities.get(j);
					if (i!=j) {						
						//get 4 paths - a1n1 -> a2n1, a1n1 -> a2n2, a1n2 -> a2n1, a1n2 -> a2n2
						//String midp1Key = ;
						Path midp1 = spTable.get("node"+sourceActivity.edge.node1+",node"+destActivity.edge.node1);
						Path midp2 = spTable.get("node"+sourceActivity.edge.node1+",node"+destActivity.edge.node2);
						Path midp3 = spTable.get("node"+sourceActivity.edge.node2+",node"+destActivity.edge.node1);
						Path midp4 = spTable.get("node"+sourceActivity.edge.node2+",node"+destActivity.edge.node2);
						
						//initialize shortest path between activity nodes
						Path shortestMidPath = null;	
						if (midp1 != null) shortestMidPath = midp1;
						else if (midp2 != null) shortestMidPath = midp2; 
						else if (midp3 != null) shortestMidPath = midp3; 
						else if (midp4 != null) shortestMidPath = midp4; 
						
						if (midp2 != null) 		
							if (midp2.weight < shortestMidPath.weight) shortestMidPath = midp2;
						
						if (midp3 != null) 		
							if (midp3.weight < shortestMidPath.weight) shortestMidPath = midp3;
						
						if (midp4 != null) 		
							if (midp4.weight < shortestMidPath.weight) shortestMidPath = midp4;
						
						if (shortestMidPath != null) {
							String shortestMidPathFirstNode = shortestMidPath.edges.get(0).node1;
							String shortestMidPathLastNode = shortestMidPath.edges.get(shortestMidPath.edges.size()-1).node2;
							
							//get shortest path to between activities based on shortestMidPath
							ArrayList<Edge> pathEdges = new ArrayList<Edge>();
							if (shortestMidPathFirstNode.equals(sourceActivity.edge.node1)) {
								pathEdges.addAll(sourceActivity.pathToNode1.edges);
							}		
							else if (shortestMidPathFirstNode.equals(sourceActivity.edge.node2)) {
								pathEdges.addAll(sourceActivity.pathToNode2.edges);
							}
													
							if (shortestMidPathLastNode.equals(destActivity.edge.node1)) {
								pathEdges.addAll(destActivity.pathToNode1.edges);
							}		
							else if (shortestMidPathLastNode.equals(destActivity.edge.node2)) {
								pathEdges.addAll(destActivity.pathToNode2.edges);
							}
							
							//path between activities
							Path p = new Path(pathEdges);
							
							p.statistic = graph.getStatistic(p,true);				
							//test each shortest path to see if it meets the scan statistic
							if (p.statistic >= threshold) {
								if (p.numActivities > 1) {
									candidatePaths.add(p.clone());	
									if (maxPath == null || maxPath.statistic < p.statistic)
										maxPath = p;
								}
							}
						}
					}				
				}		
			}
    	}
    	//nodes outnumber activities
    	else {
			graph.dynamicSegmentation(); // dynamically segment the graph	
			ArrayList<String> nodes = graph.dynamicNodes;
	
			for (int i = 0; i < nodes.size(); i++) {			
				String currentNode = nodes.get(i);
				if (likelihoodPruning)
					sp.processAllSPNodesLP(currentNode);
				else
					sp.processAllSPNodes(currentNode); //calculating shortest paths from current node to all other nodes			
				for (int j = 0; j < nodes.size(); j++) {
					if (nodes.get(j) != currentNode) {
						
						if (!graph.dynamicNodes.contains(nodes.get(j)))
							continue;
									
						ArrayList<Edge> a = sp.getShortestPathEdgesNoProcessing(currentNode, nodes.get(j));	
						if (a.size() > 0) {
							Path p = new Path(a);					
							p.statistic = graph.getStatistic(p,true);
							spTable.put("node"+currentNode+",node"+nodes.get(j), p); //store shortest path for later					
							//test each shortest path to see if it meets the scan statistic
							if (p.statistic >= threshold) {
								if (p.numActivities > 1) {
									candidatePaths.add(p.clone());	
									if (maxPath == null || maxPath.statistic < p.statistic)
										maxPath = p;
								}
							}
						}
					}
				}		
			}
    	}
    	return maxPath;
    }
    
 
    
    /**
     * Main driver
     * java SigRoutes graphFile numSimulations pvalue lrThreshold likelihoodPruning monteCarloSpeedup
     * @param args
     * @throws CloneNotSupportedException
     * @throws IOException
     */
	public static void main (String[] args) throws CloneNotSupportedException, IOException {		
		
		long ts1 = System.currentTimeMillis();
		
		nodesFile = "/Users/xuntang/Desktop/orlando_nodes.txt";
		edgesFile = "/Users/xuntang/Desktop/orlando_edges.txt";
		activitiesFile = "/Users/xuntang/Desktop/orlando_activities.txt";

		double numSimulations = 100;
		double pvalue = 0.05;
		lrThreshold = 20;		
		boolean likelihoodPruning = false;
		boolean monteCarloSpeedup = false;
		boolean dynamicSegmentation = false;
		boolean dynamicSegmentationHeirarchical = false;
		boolean dynamicSegmentationHeirarchicalActive = false;
		
		Graph g = new Graph(nodesFile, edgesFile, activitiesFile, lrThreshold);
			
		System.out.println("------------------------------------------------------------------------------------------------------------------------------");
		System.out.println("Nodes File: " + nodesFile + ", Edges File: " + edgesFile + ", Activities File: " + activitiesFile);
		System.out.println("Nodes: " + g.nodes.size() + ", edges: " + g.edges.size() + ", activities: " + g.activities.size() + ", weight: " + g.totalWeight + ", anr: " + (double) g.activeNodes.size()/g.nodes.size());			
		System.out.println("Simulations: " + numSimulations + ", p-value: " + pvalue + ", LR: " + lrThreshold + ", LP: " + likelihoodPruning + 
							", MS: " + monteCarloSpeedup + ", DS: " + dynamicSegmentation + ", DSH: " + dynamicSegmentationHeirarchical + ", DSA: " + dynamicSegmentationHeirarchicalActive);
		System.out.println("------------------------------------------------------------------------------------------------------------------------------");
		
		/*
		for (int i = 0; i < g.edges.size(); i++){
			System.out.println("i: "+ i + ": " + g.edges.get(g.edgeList.get(i).node1+","+g.edgeList.get(i).node2).numActivities);
			System.out.println("i: "+ i + ": " + g.originalEdges.get(g.edgeList.get(i).node1+","+g.edgeList.get(i).node2).numActivities+" original");
		}
		*/
		
		SigRoutes routeScan = new SigRoutes(g, numSimulations, pvalue);		
		ArrayList<Path> paths2 = routeScan.getSignificantRoutes (likelihoodPruning, 
																 monteCarloSpeedup, 
																 dynamicSegmentation,
																 dynamicSegmentationHeirarchical,
																 dynamicSegmentationHeirarchicalActive);
	
		System.out.println();
		System.out.println("Significant Routes output: ");
		for (int i = 0; i < paths2.size(); i++) {
			Path path = paths2.get(i);
			System.out.println("Route " + (i+1) + ": " + path);

			/*for (int j = 0; j < path.edges.size(); j++) {
				System.out.println("(\"from_node\" = " + path.edges.get(j).node1 + " AND \"to_node\" = " + path.edges.get(j).node2 + ") or");
				System.out.println("(\"from_node\" = " + path.edges.get(j).node2 + " AND \"to_node\" = " + path.edges.get(j).node1 + ") or"); 
			}*/
		}

		//System.out.println("done");
	    System.out.println("total time: " + (System.currentTimeMillis() - ts1));
	}
}
