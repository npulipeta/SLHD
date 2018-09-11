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


public class SigRoutes2 {

	public double mcIterations;
	public double pvalue;
	
	public Graph graph;
	public ShortestPath sp; 
	public Hashtable<String, Path> spTable;
	public ArrayList<Path> dominantPaths;
	public ArrayList<Double> statistics;
	public double threshold;
	
	/**
	 * Construct a SigRoutes object from the given graph
	 * @param graph
	 */
	public SigRoutes2(Graph graph, double mcIterations, double pvalue) {		
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
												boolean dynamicSegmentation) throws CloneNotSupportedException {
		
		sp.dynamicSegmentation = dynamicSegmentation;
		
		//System.out.println("Generating candidate paths...");
		ArrayList<Path> candidatePaths = new ArrayList<Path>();	
		ArrayList<Path> significantPaths = new ArrayList<Path>();	

		Path maxPath = null;

		long ts1 = System.currentTimeMillis();
		//if dynamic segmentation and more activities than nodes
		if (dynamicSegmentation && graph.numActivities > graph.originalNodes.size()) {
			sp = new ShortestPath(graph);
			System.out.println("Activities outnumber nodes...");
			ArrayList<String> nodes = graph.activeNodes;
			for (int i = 0; i < nodes.size(); i++) {		
				String currentNode = nodes.get(i);
				sp.processAllSPNodes(currentNode);
				for (int j = 0; j < graph.nodes.size(); j++) {		
					if (graph.nodes.get(j) != currentNode) {
						ArrayList<Edge> a = sp.getShortestPathEdgesNoProcessing(currentNode, graph.nodes.get(j));				
						if (a.size() > 0) {
							Path p = new Path(a);					
							spTable.put("node"+currentNode+",node"+graph.nodes.get(j), p); //store shortest path for later					
						}						
					}	
				}
			}		
			
			graph.dynamicSegmentation();
			
			//for each activity, evaluate the shortest path between each pair of activities
			for (int i = 0; i < graph.activities.size(); i++) {
				Activity sourceActivity = graph.activities.get(i);
				for (int j = 0; j < graph.activities.size(); j++) {
					Activity destActivity = graph.activities.get(j);
					if (i!=j) {						
						//get 4 paths - a1n1 -> a2n1, a1n1 -> a2n2, a1n2 -> a2n1, a1n2 -> a2n2
						String midp1Key = "node"+sourceActivity.edge.node1+",node"+destActivity.edge.node1;
						Path midp1 = spTable.get(midp1Key);
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
							
							p.statistic = graph.getStatistic(p,dynamicSegmentation);				
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
		else {	
			System.out.println("Nodes outnumber activities");		
			ArrayList<String> nodes = new ArrayList<String>();
			//nodes = graph.nodes;			
			if (dynamicSegmentation) {
				graph.dynamicSegmentation();
				nodes = graph.dynamicNodes;
			}
			else {
				nodes = graph.nodes;
			}
			
			for (int i = 0; i < nodes.size(); i++) {			
				String currentNode = nodes.get(i);
				if (likelihoodPruning)
					sp.processAllSPNodesLP(currentNode);
				else
					sp.processAllSPNodes(currentNode); //calculating shortest paths from current node to all other nodes			
				for (int j = 0; j < nodes.size(); j++) {
					if (nodes.get(j) != currentNode) {
						
						if (dynamicSegmentation) {
							if (!graph.dynamicNodes.contains(nodes.get(j)))
								continue;
						}
						
						ArrayList<Edge> a = sp.getShortestPathEdgesNoProcessing(currentNode, nodes.get(j));	
						if (a.size() > 0) {
							Path p = new Path(a);					
							p.statistic = graph.getStatistic(p,dynamicSegmentation);
							spTable.put("node"+currentNode+",node"+nodes.get(j), p); //store shortest path for later					
							//test each shortest path to see if it meets the scan statistic
							if (p.statistic >= threshold) {
								//if (p.numActivities > 1) {
									candidatePaths.add(p.clone());	
									if (maxPath == null || maxPath.statistic < p.statistic)
										maxPath = p;
								//}
							}
						}
					}
				}		
			}
		
		}

		//if there is no max path then terminate
		if (maxPath == null)
			return significantPaths;
		
		maxPath = maxPath.clone();
		System.out.println("original max path: " + maxPath);
		//System.out.println("candidatePaths: " + candidatePaths);
		System.out.println("candidatePaths size: " + candidatePaths.size());
		System.out.println("Time to calculate candidate paths: " + (System.currentTimeMillis() - ts1));
		
		
		//System.out.println("Calculating significance through Monte Carlo simulations...");
		//test significance of candidate paths
		significantPaths = getSignificance(candidatePaths, maxPath, likelihoodPruning, monteCarloSpeedup,dynamicSegmentation);	
		//significantPaths = getSignificanceSameRoutes(candidatePaths, maxPath, likelihoodPruning, monteCarloSpeedup,dynamicSegmentation);	
		//significantPaths = candidatePaths;
		
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
		
		return dominantPaths;
	}	
	
	/**
	 * Calculate the significance of the given set of paths considering all shortest paths
	 * @param paths
	 * @return
	 */
	//TODO: Optimize for cases when number of activities exceeds number of nodes
	ArrayList<Path> getSignificance (ArrayList<Path> paths, 
									 Path originalMaxPath,
									 boolean likelihoodPruning, 
									 boolean monteCarloSpeedup,
									 boolean dynamicSegmentation) {

		ArrayList<Path> significantPaths = new ArrayList<Path>();
		
		if (paths.size() == 0 || paths==null)
			return significantPaths;
		
		int numTimesPValueBeaten = 0;
		//graph.resetGraph(); //reset the graph to its original state
		for (int x = 0; x < mcIterations; x++) {
			//shuffle number of activities on each edge
			
			ArrayList<String> nodes = new ArrayList<String>();	
			
			graph.shuffleGraphActivitiesDynamic2();
			if (dynamicSegmentation) {
				graph.dynamicSegmentation();
				nodes = graph.dynamicNodes;
			}
			else {
				nodes = graph.nodes;	
			}
			//graph.shuffleGraphActivities();
			//nodes = graph.nodes;
			
			Path maxPath = null;	
			boolean speedup = false;
			
			
			for (int i = 0; i < nodes.size(); i++) {			
				String currentNode = nodes.get(i);
				if (likelihoodPruning)
					sp.processAllSPNodesLP(currentNode);
				else
					sp.processAllSPNodes(currentNode); //calculating shortest paths from current node to all other nodes			
				for (int j = 0; j < nodes.size(); j++) {
					if (nodes.get(j) != currentNode) {
						
						if (dynamicSegmentation) {
							if (!graph.dynamicNodes.contains(nodes.get(j)))
								continue;
						}
						
						ArrayList<Edge> a = sp.getShortestPathEdgesNoProcessing(currentNode, nodes.get(j));	
						if (a.size() > 0) {
							Path p = new Path(a);					
							p.statistic = graph.getStatistic(p,dynamicSegmentation);
							//p.statistic = graph.getStatistic(p,false);
							
							//if (p.numActivities <= 1) continue;
							
							//test each shortest path to see if it meets the likelihood ratio threshold			
							if (maxPath == null || maxPath.statistic < p.statistic)
								maxPath = p;	
							
							if (monteCarloSpeedup && p.statistic > originalMaxPath.statistic) {
								//System.out.println("MC speedup!!");
								addStatistic(p.statistic);
								speedup = true;
								numTimesPValueBeaten++;
								
								//no paths will be significant after consistently finding paths in each simulation
								if ((numTimesPValueBeaten/mcIterations) > pvalue) {
									//System.out.println("MC speedup 2!!");
									return significantPaths; 
								}
								break;
							}
															
						}
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

		//calculate p-value
		for (int i = 0; i < paths.size(); i++) {
			Path path = paths.get(i);
			double pathStat = path.statistic;		
			int pNum = 1;
			for (int j = 0; j < statistics.size()-1; j++) {
				if (pathStat < statistics.get(j)) {
					pNum = pNum + 1;
				}
			}	
			double pathPVal = (double) (pNum/mcIterations);
			//System.out.println(path + ", pVal: " + pathPVal);
			if (pathPVal <= pvalue) {			
				path.pVal = pathPVal;
				significantPaths.add(path);		
			}
		}
		return significantPaths;
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

		ArrayList<Path> significantPaths = new ArrayList<Path>();
		
		if (paths.size() == 0 || paths==null)
			return significantPaths;
		
		int numTimesPValueBeaten = 0;
		graph.resetGraph(); //reset the graph to its original state
		for (int x = 0; x < mcIterations; x++) {
			//shuffle number of activities on each edge
			
			ArrayList<String> nodes = new ArrayList<String>();	
			
			/*graph.shuffleGraphActivitiesDynamic2();
			if (dynamicSegmentation) {
				graph.dynamicSegmentation();
				nodes = graph.dynamicNodes;
			}
			else {
				nodes = graph.nodes;	
			}*/
			graph.shuffleGraphActivities();
			nodes = graph.nodes;
			
			Path maxPath = null;	
			boolean speedup = false;
			
			
			for (int i = 0; i < nodes.size(); i++) {			
				String currentNode = nodes.get(i);
				if (likelihoodPruning)
					sp.processAllSPNodesLP(currentNode);
				else
					sp.processAllSPNodes(currentNode); //calculating shortest paths from current node to all other nodes			
				for (int j = 0; j < nodes.size(); j++) {
					if (nodes.get(j) != currentNode) {
						
						/*if (dynamicSegmentation) {
							if (!graph.dynamicNodes.contains(nodes.get(j)))
								continue;
						}*/
						
						ArrayList<Edge> a = sp.getShortestPathEdgesNoProcessing(currentNode, nodes.get(j));	
						if (a.size() > 0) {
							Path p = new Path(a);					
							//p.statistic = graph.getStatistic(p,dynamicSegmentation);
							p.statistic = graph.getStatistic(p,false);
							
							if (p.numActivities <= 1) continue;
							
							if (monteCarloSpeedup && p.statistic > originalMaxPath.statistic) {
								//System.out.println("MC speedup!!");
								addStatistic(p.statistic);
								speedup = true;
								numTimesPValueBeaten++;
								
								//no paths will be significant after consistently finding paths in each simulation
								if ((numTimesPValueBeaten/mcIterations) > pvalue) {
									//System.out.println("MC speedup 2!!");
									return significantPaths; 
								}
								break;
							}
												
							//test each shortest path to see if it meets the likelihood ratio threshold			
							if (maxPath == null || maxPath.statistic < p.statistic)
								maxPath = p;				
						}
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

		//calculate p-value
		for (int i = 0; i < paths.size(); i++) {
			Path path = paths.get(i);
			double pathStat = path.statistic;		
			int pNum = 1;
			for (int j = 0; j < statistics.size()-1; j++) {
				if (pathStat < statistics.get(j)) {
					pNum = pNum + 1;
				}
			}	
			double pathPVal = (double) (pNum/mcIterations);
			//System.out.println(path + ", pVal: " + pathPVal);
			if (pathPVal <= pvalue) {			
				path.pVal = pathPVal;
				significantPaths.add(path);		
			}
		}
		return significantPaths;
	}
	
	/**
	 * Add to the collection of statistics for Monte Carlo simulations
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
     * Main driver
     * java SigRoutes graphFile numSimulations pvalue lrThreshold likelihoodPruning monteCarloSpeedup
     * @param args
     * @throws CloneNotSupportedException
     * @throws IOException
     */
	public static void main (String[] args) throws CloneNotSupportedException, IOException {		
		
		long ts1 = System.currentTimeMillis();
		
		if (args.length != 9) {
			System.err.println("Insufficient arguments provided. Example arguments: java SigRoutes graphFile numSimulations pvalue lrThreshold likelihoodPruning monteCarloSpeedup");
			System.exit(1);
		}
			
		String nodesFile = args[0];
		String edgesFile = args[1];
		String activitiesFile = args[2];
		double numSimulations = Double.parseDouble(args[3]);
		double pvalue = Double.parseDouble(args[4]);;
		double lrThreshold = Double.parseDouble(args[5]);;		
		boolean likelihoodPruning = Boolean.parseBoolean(args[6]);
		boolean monteCarloSpeedup = Boolean.parseBoolean(args[7]);
		boolean dynamicSegmentation = Boolean.parseBoolean(args[8]);
		
		
		Graph g = new Graph(nodesFile, edgesFile, activitiesFile, lrThreshold);
			
		System.out.println("------------------------------------------------------------------------------------------------------------------------------");
		System.out.println("Nodes File: " + nodesFile + ", Edges File: " + edgesFile + ", Activities File: " + activitiesFile);
		System.out.println("Nodes: " + g.nodes.size() + ", edges: " + g.edges.size() + ", activities: " + g.activities.size() + ", weight: " + g.totalWeight + ", anr: " + (double) g.activeNodes.size()/g.nodes.size());			
		System.out.println("Simulations: " + numSimulations + ", p-value: " + pvalue + ", LR: " + lrThreshold + ", LP: " + likelihoodPruning + ", MS: " + monteCarloSpeedup + ", DS: " + dynamicSegmentation);
		System.out.println("------------------------------------------------------------------------------------------------------------------------------");
		
		
		SigRoutes2 routeScan = new SigRoutes2 (g,numSimulations,pvalue);		
		ArrayList<Path> paths2 = routeScan.getSignificantRoutes (likelihoodPruning, monteCarloSpeedup, dynamicSegmentation);

		
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
	    System.out.println("time: " + (System.currentTimeMillis() - ts1));
	}
}
