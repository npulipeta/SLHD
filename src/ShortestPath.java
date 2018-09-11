import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.PriorityQueue;

/**
 * 
 * Implementation of Dijkstra's Shortest Path Algorithm.<br><br>
 * 
 * @see 
 * Dijkstra, E. W. (1959). "A note on two problems in connexion with graphs". 
 * Numerische Mathematik 1: 269�271. http://www-m3.ma.tum.de/twiki/pub/MN0506/WebHome/dijkstra.pdf.
 */
public class ShortestPath {

	/**
	 * Attributes
	 */
	public Graph graph;
	public PriorityQueue<String> unvisitedNodes;
	public ArrayList<String> visitedNodes;
	public Hashtable<String, Double> lowestWeights; 
	public Hashtable<String, Double> activitiesForLowestWeights;           //keep activities for each shortest path found
	public Hashtable<String, Double> statisticUpperboundForLowestWeights;  //keep statistic upper bound for each shortest path found 
	public double shortestTreeActivities;
	public boolean dynamicSegmentation;
	
	public Hashtable<String, String> previousNodes;
	public Hashtable<String, Edge> previousEdges;
    private final Comparator<String> lowestWeightComparator = new Comparator<String>() {
	    public int compare(String node1ID, String node2ID) {	        
	        double lowestWeightNode1 = getLowestWeight(node1ID);
	        double lowestWeightNode2 = getLowestWeight(node2ID);	
	        if (lowestWeightNode1 > lowestWeightNode2) return 1;	        
	        else if (lowestWeightNode1  < lowestWeightNode2) return -1;
	        else return 0;	                       
	    }
    };			 
	
	/**
	 * Constructs a ShortestPath object
	 * @param graph
	 */
	public ShortestPath(Graph graph) {
		this.graph = graph;
		unvisitedNodes = new PriorityQueue<String>(10000, lowestWeightComparator);
		visitedNodes = new ArrayList<String>(); 
		lowestWeights = new Hashtable<String, Double>();
		activitiesForLowestWeights  = new Hashtable<String, Double>();
		statisticUpperboundForLowestWeights = new Hashtable<String, Double>();
		shortestTreeActivities = 0;
		dynamicSegmentation = false;
		
		previousNodes = new Hashtable<String, String>();
		previousEdges = new Hashtable<String, Edge>(); 
	}
	
	/**
	 * Processes shortest path from source to destination
	 * @param source - source node
	 * @param destination - destination node
	 */
	public void processShortestPathNodes(String sourceID, String destinationID) {		
		clearOld();
		setLowestWeight(sourceID, 0);
		unvisitedNodes.add(sourceID);	
		if (dynamicSegmentation) shortestTreeActivities = 1;	
        String currentNodeID;
        while ((currentNodeID = unvisitedNodes.poll()) != null) {            
            if (currentNodeID == destinationID) break;
            visitedNodes.add(currentNodeID);
            
            Edge edge = previousEdges.get(currentNodeID);
            if (edge != null)
            	updateShortestTreeActivities(edge);
            	//shortestTreeActivities += edge.numActivities;
            
            processNeighborNodes(currentNodeID);
        }        
	}	
	
	/**
	 * Processes shortest path from source to all other nodes
	 * @param source - source node
	 * @param destination - destination node
	 */
	public void processAllSPNodes(String sourceID) {		
		clearOld();		
		setLowestWeight(sourceID, 0);		
		unvisitedNodes.add(sourceID);		
		if (dynamicSegmentation) shortestTreeActivities = 1;	
        String currentNodeID;        
        while ((currentNodeID = unvisitedNodes.poll()) != null) {            
            visitedNodes.add(currentNodeID);    
            
            Edge edge = previousEdges.get(currentNodeID);
            if (edge != null)
            	updateShortestTreeActivities(edge);
            	//shortestTreeActivities += edge.numActivities;
            
            double currentStat = getStatisticUpperBound(currentNodeID);
            processNeighborNodes(currentNodeID);                                    
        }        
	}		
	
	/**
	 * Processes shortest path from source to all other nodes with likelihood pruning
	 * @param source - source node
	 * @param destination - destination node
	 */
	public void processAllSPNodesLP(String sourceID) {		
		clearOld();		
		setLowestWeight(sourceID, 0);
		unvisitedNodes.add(sourceID);	
		if (dynamicSegmentation) shortestTreeActivities = 1;	
        String currentNodeID;        
        while ((currentNodeID = unvisitedNodes.poll()) != null) {            
            visitedNodes.add(currentNodeID);    
           
            
            Edge edge = previousEdges.get(currentNodeID);
            if (edge != null)
            	updateShortestTreeActivities(edge);
            	//shortestTreeActivities += edge.numActivities;
 
            double currentStat = getStatisticUpperBound(currentNodeID);
            if (currentStat >= graph.threshold)
            	processNeighborNodes(currentNodeID);      
        }        
	}
	
	
    /**
	 * Calculates the new lowest weight for neighboring nodes. Weights
	 * might need to be updated if a lower weight is found
	 * 
	 * @param nodeID
	 */
    private void processNeighborNodes(String nodeID) {  	
    	ArrayList<Edge> outEdges = graph.getOutEdges(nodeID);

    	if (outEdges != null) { 
	    	for (int i = 0; i < outEdges.size(); i++){
	    		
	    		String neighbor = outEdges.get(i).node2;
	    		
	    		//if neighbor node has already been visited then skip it
	    		if (visitedNodes.contains(neighbor)) continue;
	    		
	    		double lowestWeight = getLowestWeight(nodeID) + outEdges.get(i).weight;		
	    		
	    		if (lowestWeight < getLowestWeight(neighbor)) {
	            	//assign new lowest weight and mark unvisited
	                setLowestWeight(neighbor, lowestWeight);
	                                
	                //assign previous node in shortest path
	                setPreviousNode(neighbor, nodeID);
	                
	                //assign previous edge in shortest path
	                setPreviousEdge(neighbor, outEdges.get(i));
		    		
		    		setStatisticUpperBound(neighbor, calculateStatisticUpperBound(nodeID,neighbor,outEdges.get(i),lowestWeight));	   
	    		}    		
	    	}
    	}
    }   
	
    
    double calculateStatisticUpperBound(String nodeID, String neighbor, Edge neighborEdge, double lowestWeight) {

    	 //assign number of activities in shortest path
    	double pathActivities = getActivities(nodeID);
    	if (dynamicSegmentation && neighborEdge.node2.startsWith("a.")) {
    		pathActivities += 1;
    	}
    	else {
    		pathActivities += neighborEdge.numActivities;      
    	}
      	
    	setActivities(neighbor, pathActivities); 
        
        //assign statistic upper bound in shortest path       
        double a = graph.numActivities - 1 - shortestTreeActivities;
        
        double w = lowestWeight;
        double minNeigborWeight = Double.MAX_VALUE;
        ArrayList<Edge> neighborNeighbors = graph.nodeAdjacenciesTable.get(neighbor);
        if (neighborNeighbors != null) {
            for (int j = 0; j < neighborNeighbors.size(); j++) {
            	if (neighborNeighbors.get(j).weight < minNeigborWeight)
            		minNeigborWeight = neighborNeighbors.get(j).weight;
            }
		}
        if (minNeigborWeight != Double.MAX_VALUE)
        	w += minNeigborWeight;
                 
		double div1 = a/w;
		double aComp = graph.numActivities - a;
		double wComp = (graph.totalWeight/2) - w; //account for undirected edges
		double div2 = aComp/wComp;		
		double stat = 0;
		if (wComp > 0 && div2 > 0) 
			stat = div1/div2;
		
		return stat;
    }
 
	/**
	 * Returns shortest path from source to destination as a list of paths<br>
	 * NOTE: call processAllSPEdges before calling this method
	 * @return shortest path from source to destination as a list of paths
	 * @param source
	 * @param destination
	 */
	public ArrayList<Edge> getShortestPathEdgesNoProcessing(String sourceID, String destinationID) {	
		ArrayList<Edge> shortestPathList = new ArrayList<Edge>();
		
		if (getPreviousNode(destinationID) == null)
			return shortestPathList;
		
		try {
			for (String nodeID = destinationID; nodeID != null; nodeID = getPreviousNode(nodeID)) {
				Edge edge = getPreviousEdge(nodeID);
				if (edge != null) {
					shortestPathList.add(edge);
				}
			}			
		}
		catch(Exception e) { }		
		Collections.reverse(shortestPathList); 		
		return shortestPathList;	
	}	
	
	/**
	 * Update number of activities for shortest path tree based on whether the graph is dynamically segmented
	 * @param edge
	 */
	void updateShortestTreeActivities(Edge edge) {
		if (dynamicSegmentation)
			shortestTreeActivities += edge.calculateNumActivitiesDynamicForSPTree();
		else
			shortestTreeActivities += edge.numActivities;	
	}
	
	
	/**
	 * Clears old value from relevant globals
	 */
	public void clearOld(){
		visitedNodes.clear();
		unvisitedNodes.clear();
        lowestWeights.clear();
        activitiesForLowestWeights.clear();
        statisticUpperboundForLowestWeights.clear();
        shortestTreeActivities = 0;
        previousNodes.clear();
        previousEdges.clear();
	}
	
	/**
	 * Sets the lowest weight for the given node. Note that we remove and re-add
	 * the node to the unvisited list so that it gets reordered in the PriorityQueue 
	 * properly in the event that the node is already in the tree and now has 
	 * a lower weight.
	 * @param node
	 * @param weight
	 */
	public void setLowestWeight(String nodeID, double weight){
		unvisitedNodes.remove(nodeID);
		lowestWeights.put(nodeID, weight);
        unvisitedNodes.add(nodeID);       		
	}
	
    /**
     * Returns the lowest weight from the source node to the given node. If there is no weight
     * recorded then return infinite weight which is represented by Integer.MAX_VALUE
     * @param node - the node id for which we want to get the lowest weight
     * @return the lowest weight from the source node to the given node
     */    
    public double getLowestWeight(String nodeID) {
        Double cost = lowestWeights.get(nodeID);
        
        if (cost == null) {
        	return Integer.MAX_VALUE;
        }
        else {
        	return cost;
        }        
    }	
    
    /**
     * Returns the number of activities from the source node to the given node. 
     * @param nodeID
     * @return
     */
    public double getActivities(String nodeID) { 	
    	Double activities = activitiesForLowestWeights.get(nodeID);
    	if (activities == null) {
    		return 0;
    	}
    	else {
    		return activities;
    	}    	
    }
    
    /**
     * Set the number of activities in the shortest path
     * @param node1ID
     * @param activities
     */
    public void setActivities(String node1ID, double activities) {
    	activitiesForLowestWeights.put(node1ID, activities);
    } 
    
    /**
     * Returns the upper bound likelihood ratio in the shortest path
     * @param nodeID
     * @return
     */
    public double getStatisticUpperBound(String nodeID) { 	
    	Double statisticUpperBound = statisticUpperboundForLowestWeights.get(nodeID);
    	if (statisticUpperBound == null) {
    		return Double.MAX_VALUE;
    	}
    	else {
    		return statisticUpperBound;
    	}    	
    }
    
    /**
     * Set the upper bound likelihood ratio in the shortest path
     * @param node1ID
     * @param statisticUpperBound
     */
    public void setStatisticUpperBound(String node1ID, double statisticUpperBound) {
    	statisticUpperboundForLowestWeights.put(node1ID, statisticUpperBound);
    } 
    
    /**
     * Returns previous node in shortest path
     * @param node
     * @return previous node in shortest path
     */
    public String getPreviousNode(String nodeID) {
        return previousNodes.get(nodeID);
    }
    
    /**
     * Sets previous node in shortest path
     * 
     * @param node1
     * @param node2
     */
    public void setPreviousNode(String node1ID, String node2ID) {
    	previousNodes.put(node1ID, node2ID);
    }    
    
    /**
     * Returns previous edge in shortest path
     * @param node
     * @return previous edge in shortest path
     */
    public Edge getPreviousEdge(String nodeID) {
        return previousEdges.get(nodeID);
    }
    
    /**
     * Sets previous edge in shortest path
     * @param node1
     * @param node2
     */
    public void setPreviousEdge(String nodeID, Edge edge) {
    	previousEdges.put(nodeID, edge);
    }   
    
    
}
