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


public class SigRoutes_backup {

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
	public SigRoutes_backup(Graph graph, double mcIterations, double pvalue) {		
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
		
		ArrayList<String> nodes = new ArrayList<String>();
		
		if (dynamicSegmentation) {
			graph.dynamicSegmentation();
			nodes = graph.dynamicNodes;
		}
		else {
			//nodes = graph.activeNodes;
			nodes = graph.nodes;
		}
		
		
		//TODO: Optimize for dynamic segmentation
		//if there are more activities than nodes, first calculate all pair shortest paths between nodes
		//for each pair of activities, the shortest path is the shortest path between each pair of nodes 
		
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

		//if there is no max path then terminate
		if (maxPath == null)
			return significantPaths;
		
		maxPath = maxPath.clone();
		System.out.println("original max path: " + maxPath);
		//System.out.println("candidatePaths: " + candidatePaths);
		System.out.println("candidatePaths size: " + candidatePaths.size());
		
		//System.out.println("Calculating significance through Monte Carlo simulations...");
		//test significance of candidate paths
		significantPaths = getSignificance(candidatePaths, maxPath, likelihoodPruning, monteCarloSpeedup,dynamicSegmentation);	
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
		
		
		SigRoutes_backup routeScan = new SigRoutes_backup(g,numSimulations,pvalue);		
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
